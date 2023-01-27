"""modvege_lib.py

Main ModVege functions.
"""

import numpy as np

np.seterr("raise")


def leaf_area_index(
    ts_vals: dict[str, float], params: dict[str, float]
) -> float:
    """
    Calculate the leaf area index (LAI)

    Equation (12) in Jouven et al. (2006a)

    Parameters
    ----------
    params : A dictionary containing these model parameters:
        - pct_lam: Percentage of laminae in the green vegetative (GV) biomass
            compartment; default is 0.68 (%LAM) [dimensionless]
        - sla: Specific leaf area; default is 0.033 (SLA) [m² g⁻¹]
    ts_vals : A dictionary with intermediate time series values for:
        - bm_gv: Standing biomass of the green vegetative (GV) compartment
            (BM_GV) [kg DM ha⁻¹]

    Returns
    -------
    - Leaf area index (LAI) [dimensionless]
    """

    return params["sla"] * ts_vals["bm_gv"] / 10.0 * params["pct_lam"]


def actual_evapotranspiration(pet: float, ts_vals: dict[str, float]) -> float:
    """
    Calculate the actual evapotranspiration (AET)

    AET is equivalent to potential evapotranspiration (PET) when the cover
    intercepts approximately 0.95 of the incident photosynthetically active
    radiation (PAR), i.e., when the leaf area index (LAI) > 3,
    based on Johnson and Parsons (1985).
    AET is proportional to LAI when the proportion of intercepted radiation
    is lower than 0.95, i.e. LAI < 3.

    See Equation (14) in Jouven et al. (2006a)

    Parameters
    ----------
    pet : Potential evapotranspiration (PET) [mm]
    ts_vals : A dictionary with intermediate time series values for:
        - lai: Leaf area index (LAI) [dimensionless]

    Returns
    -------
    - Actual evapotranspiration (AET) [mm]
    """

    return min(pet, pet * ts_vals["lai"] / 3.0)


def potential_growth(
    par_i: float, ts_vals: dict[str, float], params: dict[str, float]
) -> float:
    """
    Calculate potential growth (PGRO)

    See Equation (12) in Jouven et al. (2006a)

    Based on Schapendonk et al. (1998).

    The model extinction coefficient is set to a constant value of 0.6
    according to Schapendonk et al. (1998) and Bonesmo and Bélanger (2002).

    The maximum radiation use efficiency is 3 g DM MJ⁻¹ based on
    Schapendonk et al. (1998).

    Parameters
    ----------
    par_i : Incident photosynthetically active radiation (PAR_*i*) [MJ m⁻²]
    params : A dictionary containing model parameters:
        - rue_max: Maximum radiation use efficiency (RUE_max); default is 3
            [g DM MJ⁻¹]
    ts_vals : A dictionary with intermediate time series values for:
        - lai: Leaf area index (LAI) [dimensionless]

    Returns
    -------
    - Potential growth (PGRO) [kg DM ha⁻¹]
    """

    return (
        par_i * params["rue_max"] *
        (1.0 - np.exp(-0.6 * ts_vals["lai"])) * 10.0
    )


def par_function(par_i: float) -> float:
    """
    Incident photosynthetically active radiation (PAR_*i*) function
    (*f*(PAR_*i*)) needed to calculate the environmental limitation of growth
    (ENV).

    The definition has been derived from Schapendonk et al. (1998).
    This function accounts for the decrease in radiation use efficiency (RUE)
    at light intensities higher than 5 MJ m⁻².

    See Figure 2(a), Equation (13), and the section on "Growth functions" in
    Jouven et al. (2006a).

    Parameters
    ----------
    par_i : Incident photosynthetically active radiation (PAR_*i*) [MJ m⁻²]

    Returns
    -------
    - PAR_i function (*f*(PAR_*i*)) [dimensionless]
    """

    if par_i < 5.0:
        val = 1.0
    else:
        # linear gradient
        gradient = 1.0 / (5.0 - (25.0 + 30.0) / 2.0)
        intercept = 1.0 - gradient * 5.0
        val = max(gradient * par_i + intercept, 0.0)
    return val


def sum_of_temperatures(
    params: dict[str, float], ts_vals: dict[str, float],
    t_ts: list[float], day: int
) -> float:
    """
    Return the sum of temperatures for each day of the year above the minimum
    temperature for growth (ST)

    Parameters
    ----------
    t_ts : Temperature (*T*) field of the input time series data (temperature
        should be in °C)
    day : Day number
    params : A dictionary containing model parameters:
        - t_0: Minimum temperature for growth (*T*₀); default is 4 [°C]
    ts_vals : A dictionary with intermediate time series values for:
        - st: Sum of temperatures value for the previous data row (ST)
            [°C d]

    Returns
    -------
    - Sum of temperatures above *T*₀ corresponding to each day of the year
        (ST) [°C d]

    Notes
    -----
    - Degree days are measures of how cold or warm a location is
    - A *degree day* compares the mean (the average of the high and low)
      outdoor temperatures recorded for a location to a
      *standard temperature*
    - Also known as heat units or thermal units
    - All species of plants have a cutoff temperature below which no
      development occurs (developmental threshold)
    - Degree days are accumulated whenever the temperature exceeds the
      predetermined developmental threshold
    - Calculate degree days by subtracting the developmental threshold from
      the average daily temperature
    - If the average degree day value for a given day is less than zero, just
      record zero, not a negative number

    References
    ----------
    - https://hort.extension.wisc.edu/articles/degree-day-calculation/
    - https://www.eia.gov/energyexplained/units-and-calculators/degree-days.php
    """

    for i in range(day):
        if t_ts[i] > params["t_0"]:
            val = ts_vals["st"] + t_ts[i] - params["t_0"]
    return val


def ten_day_moving_avg_temperature(day: int, t_ts: list[float]) -> float:
    """
    Calculate the 10-d moving average temperature.

    See sec. "Growth functions", par. above Equation (13) in Jouven et al.
    (2006a).

    Parameters
    ----------
    t_ts : Temperature (*T*) field of the input time series data (temperature
        should be in °C)
    day : Day number

    Returns
    -------
    - 10-d moving average temperature [°C]
    """

    if (day - 1) < (10 - 1):
        # ** USING THE TEMP, NOT 10-d MOVING AVG!
        val = t_ts[day - 1]
    else:
        val = np.mean([
            t_ts[(day - 1) - j] for j in range(10 - 1, 0 - 1, -1)
        ])
    return val


def temperature_function(
    ts_vals: dict[str, float], params: dict[str, float]
) -> float:
    """
    Temperature function, *f*(*T*)

    See Figure 2(b) of Jouven et al. (2006a) and the accompanying text for more
    info; *f*(*T*) has been derived based on Schapendonk et al. (1998)

    Assume no growth takes place after a maximum temperature

    Parameters
    ----------
    ts_vals : A dictionary with intermediate time series values for:
        - t_m10: 10-d moving average temperature [°C]
    params : A dictionary containing these model parameters:
        - t_0: Minimum temperature for growth (*T*₀); default is 4 [°C]
        - t_1: Minimum temperature for optimal growth; default is 10 [°C]
        - t_2: Maximum temperature for optimal growth; default is 20 [°C]
        - t_max: Maximum temperature for growth; default is 40 [°C]

    Returns
    -------
    - Temperature function (*f*(*T*)) [dimensionless]
    """

    if (
        ts_vals["t_m10"] <=
        params["t_0"] or ts_vals["t_m10"] >=
        params["t_max"]
    ):
        val = 0.0
    elif params["t_0"] < ts_vals["t_m10"] < params["t_1"]:
        # linear relationship
        gradient = 1.0 / (params["t_1"] - params["t_0"])
        intercept = 1.0 - gradient * params["t_1"]
        val = gradient * ts_vals["t_m10"] + intercept
    elif params["t_1"] <= ts_vals["t_m10"] <= params["t_2"]:
        val = 1.0
    elif params["t_2"] < ts_vals["t_m10"] < params["t_max"]:
        # linear relationship
        gradient = 1.0 / (params["t_2"] - params["t_max"])
        intercept = 1.0 - gradient * params["t_2"]
        val = gradient * ts_vals["t_m10"] + intercept
    return val


def seasonal_effect(
    ts_vals: dict[str, float], params: dict[str, float]
) -> float:
    """
    Calculate seasonal effect (SEA) on growth, driven by the sum of
    temperatures

    SEA > 1 indicates above-ground stimulation by mobilisation of reserves;
    SEA < 1 indicates growth limitation by storage of reserves

    SEA = minSEA when ST < 200°C d, then increases and reaches maxSEA when
    (ST₁ - 200) < ST < (ST₁ - 100) (ST = ST₁ at the beginning of the
    reproductive period); during summer, SEA decreases, returning to minSEA at
    ST₂ (ST = ST₂ at the end of the reproductive period)

    "T-sum 200 date" had its origins in the Netherlands and reflects the
    amount of heat absorbed and hence the amount of energy available to
    promote grass growth, and is used by farmers as an indication of
    conditions suitable for nitrogen application to grass swards
    (Collins and Cummins, 1998).

    See Figure 3 of Jouven et al. (2006a) and the accompanying paragraphs for
    more info

    minSEA and maxSEA are functional traits arranged symmetrically around 1:
    (minSEA + maxSEA) / 2 = 1

    Parameters
    ----------
    params : A dictionary containing these model parameters:
        - max_sea: Maximum seasonal effect (maxSEA); default is 1.2
            [dimensionless]
        - min_sea: Minimum seasonal effect (minSEA); default is 0.8
            [dimensionless]
        - st_1: Sum of temperatures at the beginning of the reproductive
            period (ST₁); default is 600 [°C d]
        - st_2: Sum of temperatures at the end of the reproductive period
            (ST₂); default is 1200 [°C d]
    ts_vals : A dictionary with intermediate time series values for:
        - st: Sum of temperatures (ST) [°C d]

    Returns
    -------
    - Seasonal effect [dimensionless]
    """

    if params["st_1"] <= 200.0:
        # use a constant value if the sum of temperatures at the
        # beginning of the reproductive period is lower than the sum of
        # temperatures at the onset of reproductive growth
        val = np.mean([params["max_sea"], params["min_sea"]])
    elif ts_vals["st"] <= 200.0 or ts_vals["st"] >= params["t_2"]:
        val = params["min_sea"]
    elif (
        (params["st_1"] - 200.0) <=
        ts_vals["st"] <=
        (params["st_1"] - 100.0)
    ):
        val = params["max_sea"]
    elif 200.0 < ts_vals["st"] < (params["st_1"] - 200.0):
        # assume SEA increases linearly from minSEA at the onset of
        # growth to maxSEA
        gradient = (
            (params["max_sea"] - params["min_sea"]) /
            ((params["st_1"] - 200.0) - 200.0)
        )
        intercept = params["min_sea"] - gradient * 200.0
        val = max(gradient * ts_vals["st"] + intercept, params["min_sea"])
    elif (params["st_1"] - 100.0) < ts_vals["st"] < params["st_2"]:
        # SEA decreases linearly from maxSEA to minSEA at ST_2
        gradient = (
            (params["max_sea"] - params["min_sea"]) /
            ((params["st_1"] - 100.0) - params["st_2"])
        )
        intercept = params["min_sea"] - gradient * params["st_2"]
        val = max(gradient * ts_vals["st"] + intercept, params["min_sea"])
    return val


def water_reserves(
    ts_vals: dict[str, float], params: dict[str, float], precipitation: float
) -> float:
    """
    Calculate the water reserves (WR).

    WR vary between zero and the soil water-holding capacity (WHC).
    Precipitation (PP) fill the WHC, increasing WR, while actual
    evapotranspiration (AET) empties it.

    See Equation (14) in Jouven et al. (2006a).

    Parameters
    ----------
    precipitation : Precipitation (PP) [mm]
    ts_vals : A dictionary with intermediate time series values for:
        - wr: Water reserves (WR) [mm]
        - aet: Actual evapotranspiration (AET) [mm]
    params : A dictionary containing model parameters:
        - whc: Soil water-holding capacity (WHC) [mm]

    Returns
    -------
    - Water reserves (WR) [mm]
    """

    return min(
        max(0.0, ts_vals["wr"] + precipitation - ts_vals["aet"]),
        params["whc"]
    )


def water_stress(
    ts_vals: dict[str, float], params: dict[str, float]
) -> float:
    """
    Calculate the water stress (*W*).

    See Equation (14) in Jouven et al. (2006a)

    Parameters
    ----------
    ts_vals : A dictionary with intermediate time series values for:
        - wr: Water reserves (WR) [mm]
    params : A dictionary containing model parameters:
        - whc: Soil water-holding capacity (WHC) [mm]

    Returns
    -------
    - Water stress (*W*) [dimensionless]
    """

    return min(ts_vals["wr"] / params["whc"], 1.0)


def water_stress_function(ts_vals: dict[str, float], pet) -> float:
    """
    Water stress function (*f*(*W*)).

    See Figure 2(c) and Equation (14) of Jouven et al. (2006a).

    Based on McCall and Bishop-Hurley (2003).

    Parameters
    ----------
    ts_vals : A dictionary with intermediate time series values for:
        - w: Water stress (*W*) [dimensionless]
    pet : Potential evapotranspiration (PET) [mm]

    Returns
    -------
    - Water stress function (*f*(*W*)) [dimensionless]
    """

    if pet < 3.8:
        # linear gradients
        if ts_vals["w"] < 0.2:
            gradient = 0.8 / 0.2
            val = gradient * ts_vals["w"]
        elif ts_vals["w"] < 0.4:
            gradient = (0.95 - 0.8) / (0.4 - 0.2)
            intercept = 0.8 - gradient * 0.2
            val = gradient * ts_vals["w"] + intercept
        elif ts_vals["w"] < 0.6:
            gradient = (1.0 - 0.95) / (0.6 - 0.4)
            intercept = 1.0 - gradient * 0.6
            val = gradient * ts_vals["w"] + intercept
        else:
            val = 1.0
    elif pet <= 6.5:
        if ts_vals["w"] < 0.2:
            gradient = 0.4 / 0.2
            val = gradient * ts_vals["w"]
        elif ts_vals["w"] < 0.4:
            gradient = (0.7 - 0.4) / (0.4 - 0.2)
            intercept = 0.4 - gradient * 0.2
            val = gradient * ts_vals["w"] + intercept
        elif ts_vals["w"] < 0.6:
            gradient = (0.9 - 0.7) / (0.6 - 0.4)
            intercept = 0.9 - gradient * 0.6
            val = gradient * ts_vals["w"] + intercept
        elif ts_vals["w"] < 0.8:
            gradient = (1.0 - 0.9) / (0.8 - 0.6)
            intercept = 1.0 - gradient * 0.8
            val = 0.5 * ts_vals["w"] + 0.6
        else:
            val = 1.0
    else:
        val = ts_vals["w"]
    return val


def reproductive_function(
    params: dict[str, float], ts_vals: dict[str, float]
) -> float:
    """
    Reproductive function (REP).
    REP is zero when there is a cut due to grazing or harvesting.
    REP is also zero before and after the period of reproductive growth.

    See Equation (15) in Jouven et al. (2006a)

    Parameters
    ----------
    ts_vals : A dictionary with intermediate time series values for:
        - st: Sum of temperatures [°C d]
    params : A dictionary containing model parameters:
        - ni: Nitrogen nutritional index (NI) [dimensionless]
        - st_1: Sum of temperatures at the beginning of the reproductive
            period [°C d]
        - st_2: Sum of temperatures at the end of the reproductive period
            [°C d]
        - sr: Stocking rate [LU ha⁻¹]
        - h_grass: Minimum residual grass height; default is 0.05 [m]

    Returns
    -------
    - Reproductive function [dimensionless]
    """

    if (
        ts_vals["st"] < params["st_1"] or ts_vals["st"] > params["st_2"]
    ):
        val = 0.0
    elif (
        params["sr"] > 0.0 and
        params["st_g1"] <= ts_vals["st"] <= params["st_2"]
    ):
        val = 0.0
    elif (
        params["h_grass"] > 0.0 and
        params["st_h1"] <= ts_vals["st"] <= params["st_2"]
    ):
        val = 0.0
    else:
        val = 0.25 + ((1.0 - 0.25) * (params["ni"] - 0.35)) / (1.0 - 0.35)
    return val


def environmental_limitation(
    ts_vals: dict[str, float], params: dict[str, float], par_i: float
) -> float:
    """
    Environmental limitation of growth (ENV).

    See Equation (13) of Jouven et al. (2006a).

    Parameters
    ----------
    ts_vals : A dictionary with intermediate time series values for:
        - f_t: temperature function (*f*(*T*)) [dimensionless]
        - f_w: Water stress function (*f*(*W*)) [dimensionless]
    params : A dictionary containing model parameters:
        - ni: Nutritional index of pixel (NI) [dimensionless]
    par_i : Incident photosynthetically active radiation (PAR_i) [MJ m⁻²]

    Returns
    -------
    - Environmental limitation of growth (ENV) [dimensionless]
    """

    return (
        ts_vals["f_t"] * params["ni"] *
        par_function(par_i=par_i) * ts_vals["f_w"]
    )


def total_growth(ts_vals: dict[str, float]) -> float:
    """
    Calculate the total biomass growth (GRO)

    See Equation (11) in Jouven et al. (2006a)

    Parameters
    ----------
    ts_vals : A dictionary with intermediate time series values for:
        - pgro: Potential growth (PGRO) [kg DM ha⁻¹]
        - env: Environmental limitation of growth (ENV) [dimensionless]
        - sea: Seasonal effect (SEA) [dimensionless]

    Returns
    -------
    - Total biomass growth (GRO) [kg DM ha⁻¹]
    """

    return ts_vals["pgro"] * ts_vals["env"] * ts_vals["sea"]


def abscission(
    ts_vals: dict[str, float], params: dict[str, float], temperature: float
):
    """
    Compute abscission biomass for the dead vegetative (DV) and dead
    reproductive (DR) compartments.
    See Equation (18) and Figure 4(c) and (d) in in Jouven et al. (2006a).

    Note that abscission only occurs when T > 0.

    Parameters
    ----------
    params : A dictionary containing model parameters:
        - lls: Leaf lifespan (LLS) [500 °C d]
        - kl_dv: Basic abscission rate for the dead vegetative compartment;
            default is 0.001 (Kl_DV) [dimensionless]
        - kl_dr: Basic abscission rate for the dead reproductive compartment;
            default is 0.0005 (Kl_DR) [dimensionless]
        - st_1: Sum of temperatures at the beginning of the reproductive
            period; default is 600 (ST₁) [°C d]
        - st_2: Sum of temperatures at the end of the reproductive period;
            default is 1200 (ST₂) [°C d]
    temperature : Mean daily temperature (T) [°C]
    ts_vals : A dictionary with intermediate time series values for:
        - bm_dv: DV biomass (BM_DV) [kg DM ha⁻¹]
        - age_dv: Age of the DV compartment (AGE_DV) [°C d]
        - bm_dr: DR biomass (BM_DR) [kg DM ha⁻¹]
        - age_dr: Age of the DR compartment (AGE_DR) [°C d]

    Returns
    -------
    - An updated `ts_vals` dictionary with:
        - abs_dv: Abscission of the DV biomass (ABS_DV) [kg DM ha⁻¹]
        - abs_dr: Abscission of the DR biomass (ABS_DR) [kg DM ha⁻¹]
    """

    # abscission of the DV biomass
    if ts_vals["age_dv"] / params["lls"] < 1.0 / 3.0:
        f_age = 1.0
    elif ts_vals["age_dv"] / params["lls"] < 2.0 / 3.0:
        f_age = 2.0
    else:
        f_age = 3.0
    if temperature > 0.0:
        ts_vals["abs_dv"] = (
            params["kl_dv"] * ts_vals["bm_dv"] * temperature * f_age
        )
    else:
        ts_vals["abs_dv"] = 0.0

    # abscission of the DR biomass
    if ts_vals["age_dr"] / (params["st_2"] - params["st_1"]) < 1.0 / 3.0:
        f_age = 1.0
    elif ts_vals["age_dr"] / (params["st_2"] - params["st_1"]) < 2.0 / 3.0:
        f_age = 2.0
    else:
        f_age = 3.0
    if temperature > 0.0:
        ts_vals["abs_dr"] = (
            params["kl_dr"] * ts_vals["bm_dr"] * temperature * f_age
        )
    else:
        ts_vals["abs_dr"] = 0.0


def senescence(
    ts_vals: dict[str, float], params: dict[str, float], temperature: float
):
    """
    Senescing biomass for the GV and GR compartments.
    See Equations (16) and (17) and Figure 4(a) and (b) in Jouven et al.
    (2006a).

    No senescence occurs when *T* is between zero and *T*₀.
    When T drops below zero, senescence is driven by freezing effects and is
    proportional to |*T*|.

    Parameters
    ----------
    params : A dictionary containing model parameters:
        - k_gv: Basic senescence rate for the green vegetative compartment;
            default is 0.002 (K_GV) [dimensionless]
        - k_gr: Basic senescence rate for the green reproductive compartment;
            default is 0.001 (K_GR) [dimensionless]
        - t_0: Minimum temperature for growth; default is 4 (*T*₀) [°C]
        - lls: Leaf lifespan; default is 500 (LLS) [°C d]
        - st_1 : Sum of temperatures at the beginning of the reproductive
            period; default is 600 (ST₁) [°C d]
        - st_2 : Sum of temperatures at the end of the reproductive period;
            default is 1200 (ST₂) [°C d]
    ts_vals : A dictionary with intermediate time series values for:
        - bm_gv: GV biomass (BM_GV) [kg DM ha⁻¹]
        - age_gv: Age of the GV compartment (AGE_GV) [°C d]
        - bm_gr: Biomass available for GR (BM_GR) [kg DM ha⁻¹]
        - age_gr: Age of the GR compartment (AGE_GR) [°C d]
    temperature : Mean daily temperature (*T*) [°C]

    Returns
    -------
    - An updated `ts_vals` dictionary with:
        - Senescing GV biomass (SEN_GV) [kg DM ha⁻¹]
        - Senescing GR biomass (SEN_GR) [kg DM ha⁻¹]
    """

    # senescing GV biomass
    if ts_vals["age_gv"] / params["lls"] < 1.0 / 3.0:
        f_age = 1.0
    elif ts_vals["age_gv"] / params["lls"] < 1.0:
        # linear gradient
        gradient = (3.0 - 1.0) / (1.0 - 1.0 / 3.0)
        intercept = 3.0 - gradient * 1.0
        f_age = gradient * ts_vals["age_gv"] / params["lls"] + intercept
    else:
        f_age = 3.0
    if temperature > params["t_0"]:
        ts_vals["sen_gv"] = (
            params["k_gv"] * ts_vals["bm_gv"] * temperature * f_age
        )
    elif temperature < 0.0:
        ts_vals["sen_gv"] = (
            params["k_gv"] * ts_vals["bm_gv"] * abs(temperature)
        )
    else:
        ts_vals["sen_gv"] = 0.0

    # senescing GR biomass
    if ts_vals["age_gr"] / (params["st_2"] - params["st_1"]) < 1.0 / 3.0:
        f_age = 1.0
    elif ts_vals["age_gr"] / (params["st_2"] - params["st_1"]) < 1.0:
        # linear gradient
        gradient = (3.0 - 1.0) / (1.0 - 1.0 / 3.0)
        intercept = 3.0 - gradient * 1.0
        f_age = (
            gradient * ts_vals["age_gr"] /
            (params["st_2"] - params["st_1"]) + intercept
        )
    else:
        f_age = 3.0
    if temperature > params["t_0"]:
        ts_vals["sen_gr"] = (
            params["k_gr"] * ts_vals["bm_gr"] * temperature * f_age
        )
    elif temperature < 0.0:
        ts_vals["sen_gr"] = (
            params["k_gr"] * ts_vals["bm_gr"] * abs(temperature)
        )
    else:
        ts_vals["sen_gr"] = 0.0


def biomass_growth(ts_vals: dict[str, float]) -> float:
    """
    Calculate the growth of the GV and GR biomass compartments.
    See Equations (1) and (2) in Jouven et al. (2006a)

    Parameters
    ----------
    ts_vals : A dictionary with intermediate time series values for:
        - gro: Total growth (GRO) [kg DM ha⁻¹]
        - rep: Reproductive function (REP) [dimensionless]

    Returns
    -------
    - Growth of the GV compartment (GRO_GV) [kg DM ha⁻¹]
    - Growth of the GR compartment (GRO_GR) [kg DM ha⁻¹]
    """

    return (
        ts_vals["gro"] * (1.0 - ts_vals["rep"]),
        ts_vals["gro"] * ts_vals["rep"]
    )


def standing_biomass(ts_vals: dict[str, float], params: dict[str, float]):
    """
    Update the standing biomass for each compartment.
    See Equations (1), (2), (3), and (4) in Jouven et al. (2006a).

    Parameters
    ----------
    params : A dictionary containing model parameters:
        - sigma_gv: Rate of biomass loss with respiration for GV
            [dimensionless]
        - sigma_gr: Rate of biomass loss with respiration for GR
            [dimensionless]
    ts_vals : A dictionary with intermediate time series values for:
        - bm_gv: GV biomass (BM_GV) [kg DM ha⁻¹]
        - sen_gv: Senescence of GV compartment (SEN_GV) [kg DM ha⁻¹]
        - bm_gr: GR biomass (BM_GR) [kg DM ha⁻¹]
        - sen_gr: Senescence of GR compartment (SEN_GR) [kg DM ha⁻¹]
        - bm_dv: DV biomass (BM_DV) [kg DM ha⁻¹]
        - abs_dv: Abscission of the DV compartment (ABS_DV) [kg DM ha⁻¹]
        - bm_dr: DR biomass (BM_DR) [kg DM ha⁻¹]
        - abs_dr: Abscission of the DR compartment (ABS_DR) [kg DM ha⁻¹]
        - gro: Total growth (GRO) [kg DM ha⁻¹]
        - rep: Reproductive function (REP) [dimensionless]

    Returns
    -------
    - An updated `ts_vals` dictionary with:
        - bm_gv: GV biomass (BM_GV) [kg DM ha⁻¹]
        - bm_gr: GR biomass (BM_GR) [kg DM ha⁻¹]
        - bm_dv: DV biomass (BM_DV) [kg DM ha⁻¹]
        - bm_dr: DR biomass (BM_DR) [kg DM ha⁻¹]
    """

    # GV compartment
    ts_vals["bm_gv"] = (
        ts_vals["bm_gv"] +
        biomass_growth(ts_vals=ts_vals)[0] - ts_vals["sen_gv"]
    )

    # GR compartment
    ts_vals["bm_gr"] = (
        ts_vals["bm_gr"] +
        biomass_growth(ts_vals=ts_vals)[1] - ts_vals["sen_gr"]
    )

    # DV compartment
    ts_vals["bm_dv"] = (
        ts_vals["bm_dv"] +
        (1.0 - params["sigma_gv"]) * ts_vals["sen_gv"] - ts_vals["abs_dv"]
    )

    # DR compartment
    ts_vals["bm_dr"] = (
        ts_vals["bm_dr"] +
        (1.0 - params["sigma_gr"]) * ts_vals["sen_gr"] - ts_vals["abs_dr"]
    )


def biomass_age(
    temperature: float, params: dict[str, float], ts_vals: dict[str, float]
):
    """
    Update the age of each biomass compartment.
    See Equations (5), (6), (7), and (8) in Jouven et al. (2006a).

    The age of the residual biomass is increased daily by the mean daily
    temperature, if this temperature is positive.

    Parameters
    ----------
    params : A dictionary containing model parameters:
        - sigma_gv: Rate of biomass loss with respiration for GV
            [dimensionless]
        - sigma_gr: Rate of biomass loss with respiration for GR
            [dimensionless]
    ts_vals : A dictionary with intermediate time series values for:
        - bm_gv: GV biomass (BM_GV) [kg DM ha⁻¹]
        - age_gv: Age of the GV compartment (AGE_GV) [°C d]
        - sen_gv: Senescence of GV compartment (SEN_GV) [kg DM ha⁻¹]
        - bm_gr: GR biomass (BM_GR) [kg DM ha⁻¹]
        - age_gr: Age of the GR compartment (AGE_GR) [°C d]
        - sen_gr: Senescence of GR compartment (SEN_GR) [kg DM ha⁻¹]
        - bm_dv: DV biomass (BM_DV) [kg DM ha⁻¹]
        - age_dv: Age of the DV compartment (AGE_DV) [°C d]
        - abs_dv: Abscission of the DV compartment (ABS_DV) [kg DM ha⁻¹]
        - bm_dr: DR biomass (BM_DR) [kg DM ha⁻¹]
        - age_dr: Age of the DR compartment (AGE_DR) [°C d]
        - abs_dr: Abscission of the DR compartment (ABS_DR) [kg DM ha⁻¹]
        - gro: Total growth (GRO) [kg DM ha⁻¹]
        - rep: Reproductive function (REP) [dimensionless]
    temperature : Mean daily temperature (*T*) [°C]

    Returns
    -------
    - An updated `ts_vals` dictionary with:
        - age_gv: Age of the GV compartment (AGE_GV) [°C d]
        - age_gr: Age of the GR compartment (AGE_GR) [°C d]
        - age_dv: Age of the DV compartment (AGE_DV) [°C d]
        - age_dr: Age of the DR compartment (AGE_DR) [°C d]
    """

    # GV compartment
    if temperature > 0.0 and ts_vals["bm_gv"] > 0.0:
        ts_vals["age_gv"] = (
            (ts_vals["bm_gv"] - ts_vals["sen_gv"]) /
            (
                ts_vals["bm_gv"] - ts_vals["sen_gv"] +
                biomass_growth(ts_vals=ts_vals)[0]
            ) *
            (ts_vals["age_gv"] + temperature)
        )

    # GR compartment
    if temperature > 0.0 and ts_vals["bm_gr"] > 0.0:
        ts_vals["age_gr"] = (
            (ts_vals["bm_gr"] - ts_vals["sen_gr"]) /
            (
                ts_vals["bm_gr"] - ts_vals["sen_gr"] +
                biomass_growth(ts_vals=ts_vals)[1]
            ) *
            (ts_vals["age_gr"] + temperature)
        )

    # DV compartment
    if temperature > 0.0 and ts_vals["bm_dv"] > 0.0:
        ts_vals["age_dv"] = (
            (ts_vals["bm_dv"] - ts_vals["abs_dv"]) /
            (
                ts_vals["bm_dv"] - ts_vals["abs_dv"] +
                (1.0 - params["sigma_gv"]) * ts_vals["sen_gv"]
            ) *
            (ts_vals["age_dv"] + temperature)
        )

    # DR compartment
    if temperature > 0.0 and ts_vals["bm_dr"] > 0.0:
        ts_vals["age_dr"] = (
            (ts_vals["bm_dr"] - ts_vals["abs_dr"]) /
            (
                ts_vals["bm_dr"] - ts_vals["abs_dr"] +
                (1.0 - params["sigma_gr"]) * ts_vals["sen_gr"]
            ) *
            (ts_vals["age_dr"] + temperature)
        )
