"""modvege_lib.py

Main ModVege functions.
"""

from dataclasses import dataclass
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


@dataclass
class PARFunction:

    par_i: float

    def __call__(self) -> float:
        if self.par_i < 5.0:
            val = 1.0
        else:
            # linear gradient
            gradient = 1.0 / (5.0 - (25.0 + 30.0) / 2.0)
            intercept = 1.0 - gradient * 5.0
            val = max(gradient * self.par_i + intercept, 0.0)
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
        - t_sum: Sum of temperatures value for the previous data row (ST)
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
            val = ts_vals["t_sum"] + t_ts[i] - params["t_0"]
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


@dataclass
class SeasonalEffect:
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
    max_sea : Maximum seasonal effect (maxSEA); default is 1.2 [dimensionless]
    min_sea : Minimum seasonal effect (minSEA); default is 0.8 [dimensionless]
    t_sum : Sum of temperatures (ST) [°C d]
    st_1 : Sum of temperatures at the beginning of the reproductive period
        (ST₁); default is 600 [°C d]
    st_2 : Sum of temperatures at the end of the reproductive period
        (ST₂); default is 1200 [°C d]

    Returns
    -------
    - Seasonal effect [dimensionless]
    """

    t_sum: float
    min_sea: float = 0.8
    max_sea: float = 1.2
    st_1: float = 600.0
    st_2: float = 1200.0

    def __call__(self) -> float:
        if self.st_1 <= 200.0:
            # use a constant value if the sum of temperatures at the
            # beginning of the reproductive period is lower than the sum of
            # temperatures at the onset of reproductive growth
            val = np.mean([self.max_sea, self.min_sea])
        elif self.t_sum <= 200.0 or self.t_sum >= self.st_2:
            val = self.min_sea
        elif (
            (self.st_1 - 200.0) <= self.t_sum <= (self.st_1 - 100.0)
        ):
            val = self.max_sea
        elif 200.0 < self.t_sum < (self.st_1 - 200.0):
            # assume SEA increases linearly from minSEA at the onset of
            # growth to maxSEA
            gradient = (
                (self.max_sea - self.min_sea) / ((self.st_1 - 200.0) - 200.0)
            )
            intercept = self.min_sea - gradient * 200.0
            val = max(gradient * self.t_sum + intercept, self.min_sea)
        elif (self.st_1 - 100.0) < self.t_sum < self.st_2:
            # SEA decreases linearly from maxSEA to minSEA at ST_2
            gradient = (
                (self.max_sea - self.min_sea) /
                ((self.st_1 - 100.0) - self.st_2)
            )
            intercept = self.min_sea - gradient * self.st_2
            val = max(gradient * self.t_sum + intercept, self.min_sea)
        return val


def water_reserves(
    ts_vals: dict[str, float], params: dict[str, float],
    precipitation: float
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


@dataclass
class ReproductiveFunction:
    """
    Reproductive function (REP).
    REP is zero when there is a cut due to grazing or harvesting.
    REP is also zero before and after the period of reproductive growth.

    See Equation (15) in Jouven et al. (2006a)

    Parameters
    ----------
    n_index : Nitrogen nutritional index (NI) [dimensionless]
    t_sum : Sum of temperatures [°C d]
    st_1 : Sum of temperatures at the beginning of the reproductive period
        [°C d]
    st_2 : Sum of temperatures at the end of the reproductive period [°C d]
    stocking_rate : Stocking rate [LU ha⁻¹]
    cut_height : Grass cut height [m]

    Returns
    -------
    - Reproductive function [dimensionless]
    """

    n_index: float
    t_sum: float
    stocking_rate: float
    st_1: float = 600.0
    st_2: float = 1200.0
    cut_height: float = 0.05

    def __call__(self) -> float:
        if (
            self.stocking_rate > 0.0 and
            self.st_1 + 50.0 <= self.t_sum <= self.st_2
        ):
            val = 0.0
        elif (
            self.cut_height > 0.0 and
            self.st_2 >= self.t_sum >= self.st_2 - 50.0
        ):
            val = 0.0
        elif self.t_sum < self.st_1 or self.t_sum > self.st_2:
            val = 0.0
        else:
            val = 0.25 + ((1.0 - 0.25) * (self.n_index - 0.35)) / (1.0 - 0.35)
        return val


# def environmental_limitation(
#     ts_vals: dict[str, float], params: dict[str, float]
# ) -> float:
#     """
#     Environmental limitation of growth (ENV).

#     See Equation (13) of Jouven et al. (2006a).

#     Parameters
#     ----------
#     t_fn : temperature function (*f*(*T*)) [dimensionless]
#     n_index : Nutritional index of pixel (NI) [dimensionless]
#     par_i : Incident photosynthetically active radiation (PAR_i) [MJ m⁻²]
#     w_fn : Water stress function (*f*(*W*)) [dimensionless]

#     Returns
#     -------
#     - Environmental limitation of growth (ENV) [dimensionless]
#     """

#     return (
#         ts_vals["t_fn"] * params["ni"] *
#         PARFunction(par_i=self.par_i)() * ts_vals["w_fn"]
#     )


@dataclass
class EnvironmentalLimitation:
    """
    Environmental limitation of growth (ENV).

    See Equation (13) of Jouven et al. (2006a).

    Parameters
    ----------
    t_fn : temperature function (*f*(*T*)) [dimensionless]
    n_index : Nutritional index of pixel (NI) [dimensionless]
    par_i : Incident photosynthetically active radiation (PAR_i) [MJ m⁻²]
    w_fn : Water stress function (*f*(*W*)) [dimensionless]

    Returns
    -------
    - Environmental limitation of growth (ENV) [dimensionless]
    """

    t_fn: float
    n_index: float
    par_i: float
    w_fn: float

    def __call__(self) -> float:
        return (
            self.t_fn * self.n_index *
            PARFunction(par_i=self.par_i)() * self.w_fn
        )


@dataclass
class TotalGrowth:
    """
    Calculate the total biomass growth (GRO)

    See Equation (11) in Jouven et al. (2006a)

    Parameters
    ----------
    - pgro : Potential growth (PGRO) [kg DM ha⁻¹]
    - env : environmental limitation of growth (ENV) [dimensionless]
    - sea : seasonal effect (SEA) [dimensionless]

    Returns
    -------
    - Total biomass growth (GRO) [kg DM ha⁻¹]
    """

    pgro: float
    env: float
    sea: float

    def __call__(self) -> float:
        return self.pgro * self.env * self.sea


@dataclass
class AbscissionDV:
    """
    Compute abscission biomass for the dead vegetative (DV) compartment.
    See Equation (18) and Figure 4(c) in in Jouven et al. (2006a).

    Note that abscission only occurs when T > 0.

    Parameters
    ----------
    lls : Leaf lifespan (LLS) [500 °C d]
    kl_dv : Basic abscission rate for the dead vegetative compartment; default
        is 0.001 (Kl_DV) [dimensionless]
    temperature : Mean daily temperature (T) [°C]
    bm_dv : DV biomass (BM_DV) [kg DM ha⁻¹]
    age_dv : Average DV age (AGE_DV) [°C d]

    Returns
    -------
    - Abscission biomass (ABS_DV) [kg DM ha⁻¹]
    """

    temperature: float
    bm_dv: float
    age_dv: float
    lls: float = 500.0
    kl_dv: float = 0.001

    def __call__(self) -> float:
        if self.age_dv / self.lls < 1.0 / 3.0:
            age_fn = 1.0
        elif self.age_dv / self.lls < 2.0 / 3.0:
            age_fn = 2.0
        else:
            age_fn = 3.0
        if self.temperature > 0.0:
            val = self.kl_dv * self.bm_dv * self.temperature * age_fn
        else:
            val = 0.0
        return val


@dataclass
class AbscissionDR:
    """
    Compute abscission biomass for the dead reproductive (DR) compartment.
    See Equation (18) and Figure 4(d) in Jouven et al. (2006a).

    Note that abscission only occurs when T > 0.

    Parameters
    ----------
    kl_dr : Basic abscission rate for the dead reproductive compartment;
        default is 0.0005 (Kl_DR) [dimensionless]
    bm_dr : DR biomass (BM_DR) [kg DM ha⁻¹]
    temperature : Mean daily temperature (*T*) [°C]
    age_dr : Average DR age (AGE_DR) [°C d]
    st_1 : Sum of temperatures at the beginning of the reproductive period;
        default is 600 (ST₁) [°C d]
    st_2 : Sum of temperatures at the end of the reproductive period; default
        is 1200 (ST₂) [°C d]

    Returns
    -------
    - Abscission biomass for DR (ABS_DR) [kg DM ha⁻¹]
    """

    bm_dr: float
    temperature: float
    age_dr: float
    kl_dr: float = 0.0005
    st_1: float = 600.0
    st_2: float = 1200.0

    def __call__(self) -> float:
        if self.age_dr / (self.st_2 - self.st_1) < 1.0 / 3.0:
            age_fn = 1.0
        elif self.age_dr / (self.st_2 - self.st_1) < 2.0 / 3.0:
            age_fn = 2.0
        else:
            age_fn = 3.0
        if self.temperature > 0.0:
            val = self.kl_dr * self.bm_dr * self.temperature * age_fn
        else:
            val = 0.0
        return val


@dataclass
class SenescenceGV:
    """
    Senescent biomass for GV compartment.
    See Equations (16) and (17) and Figure 4(a) in Jouven et al. (2006a).

    No senescence occurs when *T* is between zero and *T*₀.
    When T drops below zero, senescence is driven by freezing effects and is
    proportional to |*T*|.

    Parameters
    ----------
    k_gv : Basic senescence rate for the green vegetative compartment; default
        is 0.002 (K_GV) [dimensionless]
    bm_gv : GV biomass (BM_GV) [kg DM ha⁻¹]
    temperature : Mean daily temperature (*T*) [°C]
    t_0 : Minimum temperature for growth; default is 4 (*T*₀) [°C]
    lls : Leaf lifespan; default is 500 (LLS) [°C d]
    age_gv : Average GV age (AGE_GV) [°C d]

    Returns
    -------
    - Senescent biomass for GV (SEN_GV) [kg DM ha⁻¹]
    """

    temperature: float
    age_gv: float
    bm_gv: float
    k_gv: float = 0.002
    t_0: float = 4.0
    lls: float = 500.0

    def __call__(self) -> float:
        if self.age_gv / self.lls < 1.0 / 3.0:
            age_fn = 1.0
        elif self.age_gv / self.lls < 1.0:
            # linear gradient
            gradient = (3.0 - 1.0) / (1.0 - 1.0 / 3.0)
            intercept = 3.0 - gradient * 1.0
            age_fn = gradient * self.age_gv / self.lls + intercept
        else:
            age_fn = 3.0
        if self.temperature > self.t_0:
            val = self.k_gv * self.bm_gv * self.temperature * age_fn
        elif self.temperature < 0.0:
            val = self.k_gv * self.bm_gv * abs(self.temperature)
        else:
            val = 0.0
        return val


@dataclass
class SenescenceGR:
    """
    See Equations (16) and (17) and Figure 4(b) in Jouven et al. (2006a).

    Parameters
    ----------
    k_gr : Basic senescence rate for the green reproductive compartment;
        default is 0.001 (K_GR) [dimensionless]
    bm_gr : Biomass available for GR (BM_GR) [kg DM ha⁻¹]
    temperature : Mean daily temperature (*T*) [°C]
    t_0 : Minimum temperature for growth; default is 4 (*T*₀) [°C]
    age_gr : Average GR age (AGE_GR) [°C d]
    st_1 : Sum of temperatures at the beginning of the reproductive period;
        default is 600 (ST₁) [°C d]
    st_2 : Sum of temperatures at the end of the reproductive period; default
        is 1200 (ST₂) [°C d]

    Returns
    -------
    - Senescent biomass (SEN_GR) [kg DM ha⁻¹]
    """

    age_gr: float
    temperature: float
    bm_gr: float
    k_gr: float = 0.001
    t_0: float = 4.0
    st_1: float = 600.0
    st_2: float = 1200.0

    def __call__(self) -> float:
        if self.age_gr / (self.st_2 - self.st_1) < 1.0 / 3.0:
            age_fn = 1.0
        elif self.age_gr / (self.st_2 - self.st_1) < 1.0:
            # linear gradient
            gradient = (3.0 - 1.0) / (1.0 - 1.0 / 3.0)
            intercept = 3.0 - gradient * 1.0
            age_fn = (
                gradient * self.age_gr / (self.st_2 - self.st_1) + intercept
            )
        else:
            age_fn = 3.0
        if self.temperature > self.t_0:
            val = self.k_gr * self.bm_gr * self.temperature * age_fn
        elif self.temperature < 0.0:
            val = self.k_gr * self.bm_gr * abs(self.temperature)
        else:
            val = 0.0
        return val


@dataclass
class BiomassDV:
    """
    Update dead vegetative compartment.
    See Equations (3) and (7) in Jouven et al. (2006a).

    The age of the residual biomass is increased daily by the mean daily
    temperature, if this temperature is positive

    Parameters
    ----------
    sigma_gv : Respiratory C loss during senescence of GV [dimensionless]
    sen_gv : Senescence of GV compartment (SEN_GV)
    abs_dv : Abscission of the DV compartment (ABS_DV)
    temperature : Mean daily temperature (*T*) [°C]
    bm_dv : DV biomass (BM_DV) [kg DM ha⁻¹]
    age_dv : Average DV age (AGE_DV) [°C d]

    Returns
    -------
    - Dead vegetative biomass
    - Average DV age [°C d]
    """

    bm_dv: float
    abs_dv: float
    sen_gv: float
    age_dv: float
    temperature: float
    sigma_gv: float = 0.4

    def __call__(self) -> float:
        bm_dv = self.bm_dv + (1.0 - self.sigma_gv) * self.sen_gv - self.abs_dv
        if self.temperature > 0.0 and self.bm_dv > 0.0:
            age_dv = (
                (self.bm_dv - self.abs_dv) /
                (
                    self.bm_dv - self.abs_dv +
                    (1.0 - self.sigma_gv) * self.sen_gv
                ) *
                (self.age_dv + self.temperature)
            )
        else:
            age_dv = self.age_dv
        return bm_dv, age_dv


@dataclass
class BiomassDR:
    """
    Update dead reproductive compartment.
    See Equations (4) and (8) in Jouven et al. (2006a).

    The age of the residual biomass is increased daily by the mean daily
    temperature, if this temperature is positive

    Parameters
    ----------
    sigma_gr : Respiratory C loss during senescence for GR
    sen_gr : Senescence of GR compartment
    abs_dr : Abscission of the DR compartment
    temperature : Mean daily temperature (*T*) [°C]
    bm_dr : DR biomass (BM_DR) [kg DM ha⁻¹]
    age_dr : Average DR age (AGE_DR) [°C d]

    Returns
    -------
    - Dead reproductive biomass
    - Average DR age
    """

    bm_dr: float
    abs_dr: float
    sen_gr: float
    age_dr: float
    temperature: float
    sigma_gr: float = 0.2

    def __call__(self) -> float:
        bm_dr = self.bm_dr + (1.0 - self.sigma_gr) * self.sen_gr - self.abs_dr
        if self.temperature > 0.0 and self.bm_dr > 0.0:
            age_dr = (
                (self.bm_dr - self.abs_dr) /
                (
                    self.bm_dr - self.abs_dr +
                    (1.0 - self.sigma_gr) * self.sen_gr
                ) *
                (self.age_dr + self.temperature)
            )
        else:
            age_dr = self.age_dr
        return bm_dr, age_dr


@dataclass
class GrowthGV:
    """
    Calculate the growth proportion of the GV compartment.
    See Equation (1) in Jouven et al. (2006a)

    Parameters
    ----------
    gro : Growth
    rep : Reproductive function
    """

    gro: float
    rep: float

    def __call__(self) -> float:
        return self.gro * (1.0 - self.rep)


@dataclass
class BiomassGV:
    """
    Update green vegetative compartment.
    See Equations (1) and (5) in Jouven et al. (2006a).

    The age of the residual biomass is increased daily by the mean daily
    temperature, if this temperature is positive

    Parameters
    ----------
    gro_gv : Growth GV
    sen_gv : Senescence GV
    bm_gv : GV biomass
    age_gv : GV age
    temperature : Mean daily temperature (*T*) [°C]

    Returns
    -------
    - Green vegetative biomass
    - Average GV age
    """

    gro_gv: float
    sen_gv: float
    bm_gv: float
    age_gv: float
    temperature: float

    def __call__(self) -> float:
        bm_gv = self.bm_gv + self.gro_gv - self.sen_gv
        if self.temperature > 0.0 and self.bm_gv > 0.0:
            age_gv = (
                (self.bm_gv - self.sen_gv) /
                (self.bm_gv - self.sen_gv + self.gro_gv) *
                (self.age_gv + self.temperature)
            )
        else:
            age_gv = self.age_gv
        return bm_gv, age_gv


@dataclass
class GrowthGR:
    """
    See Equation (2) in Jouven et al. (2006a)

    Parameters
    ----------
    gro : Growth
    rep : Reproductive function
    """

    gro: float
    rep: float

    def __call__(self) -> float:
        return self.gro * self.rep


@dataclass
class BiomassGR:
    """
    Update green reproductive compartment.
    See Equations (2) and (6) in Jouven et al. (2006a).

    The age of the residual biomass is increased daily by the mean daily
    temperature, if this temperature is positive

    Parameters
    ----------
    gro_gr : Growth GV
    sen_gr : Senescence GV
    bm_gr : GV biomass
    age_gr : GV age
    temperature : Mean daily temperature (*T*) [°C]

    Returns
    -------
    - Green reproductive biomass
    - Average GR age
    """

    gro_gr: float
    sen_gr: float
    bm_gr: float
    age_gr: float
    temperature: float

    def __call__(self) -> float:
        bm_gr = self.bm_gr + self.gro_gr - self.sen_gr
        if self.temperature > 0.0 and self.bm_gr > 0.0:
            age_gr = (
                (self.bm_gr - self.sen_gr) /
                (self.bm_gr - self.sen_gr + self.gro_gr) *
                (self.age_gr + self.temperature)
            )
        else:
            age_gr = self.age_gr
        return bm_gr, age_gr
