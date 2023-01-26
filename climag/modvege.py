"""modvege.py

This is a Python implementation of the ModVege pasture model, a modified
version of the original Java to Python translation by Chemin (2022).
The Java model was provided by Raphael Martin, INRAE UREP Clermont-Ferrand
for the original Python implementation.
The original ModVege pasture model was developed by Jouven et al. (2006a).

Chemin, Y. (2022). 'modvege', Python. [Online]. Available at
https://github.com/YannChemin/modvege (Accessed 6 September 2022).

References
----------
- Jouven, M., Carrère, P., and Baumont, R. (2006a). 'Model predicting dynamics
  of biomass, structure and digestibility of herbage in managed permanent
  pastures. 1. Model description', Grass and Forage Science, vol. 61, no. 2,
  pp. 112-124. DOI: 10.1111/j.1365-2494.2006.00515.x.
- Jouven, M., Carrère, P., and Baumont, R. (2006b). 'Model predicting dynamics
  of biomass, structure and digestibility of herbage in managed permanent
  pastures. 2. Model evaluation', Grass and Forage Science, vol. 61, no. 2,
  pp. 125-133. DOI: 10.1111/j.1365-2494.2006.00517.x.

Full list of references are available here:
https://www.zotero.org/groups/4706660/climag/collections/ZNH6N9Y8

Model definition
----------------
Compartments representing the structural components of herbage:
- Green vegetative (GV) - green leaves and sheath
- Dead vegetative (DV) - dead leaves and sheath
- Green reproductive (GR) - green stems and flowers
- Dead reproductive (DR) - dead stems and flowers

Compartment description:
- Standing biomass (BM) [kg DM ha⁻¹]
- Age of biomass (AGE) [°C d] (degree day - units of thermal time)
- Organic matter digestibility (OMD)

List of time series variables
----------------------------
- T: Mean daily temperature (*T*) [°C]
- PAR_i: Incident photosynthetically active radiation (PAR_i) [MJ m⁻²]
- PP: Precipitation (PP) [mm]
- PET: Potential evapotranspiration (PET) [mm]

List of parameters (constants)
------------------------------
See Tables 2 and 3 in Jouven et al. (2006a).
The functional groups and traits were parameterised for the Auvergne region in
France, which has a temperate climate.
Functional traits for Group A is used (species found in fertile sites, adapted
to frequent defoliation (e.g. perennial ryegrass; *Lolium perenne*).
Characteristic traits include high specific leaf area, high digestibility,
short leaf lifespan, and early reproductive growth and flowering).
- Specific leaf area (SLA) [0.033 m² g⁻¹]
- Percentage of laminae in the green vegetative compartment (%LAM) [0.68]
- Sum of temperatures at the beginning of the reproductive period (ST₁)
  [600 °C d]
- Sum of temperatures at the end of the reproductive period (ST₂)
  [1200 °C d]
- Maximum seasonal effect (maxSEA) [1.20]
- Minimum seasonal effect (minSEA) [0.80]
- Leaf lifespan (LLS) [500 °C d]
- Maximum organic matter digestibility of the green vegetative compartment
  (maxOMD_GV) [0.90]
- Minimum organic matter digestibility of the green vegetative compartment
  (minOMD_GV) [0.75]
- Maximum organic matter digestibility of the green reproductive compartment
  (maxOMD_GR) [0.90]
- Minimum organic matter digestibility of the green reproductive compartment
  (minOMD_GR) [0.65]
- Bulk density of the green vegetative compartment (BD_GV) [850 g DM m⁻³]
- Bulk density of the dead vegetative compartment (BD_DV) [500 g DM m⁻³]
- Bulk density of the green reproductive compartment (BD_GR) [300 g DM m⁻³]
- Bulk density of the dead reproductive compartment (BD_DR) [150 g DM m⁻³]
- Rate of biomass loss with respiration for the green vegetative compartment
  (σ_GV) [0.4]
- Rate of biomass loss with respiration for the green reproductive compartment
  (σ_GR) [0.2]
- Minimum temperature for growth (*T*₀) [4 °C]
- Minimum temperature for optimal growth (*T*₁) [10 °C]
- Maximum temperature for optimal growth (*T*₂) [20 °C]
- Basic senescence rate for the green vegetative compartment (K_GV) [0.002]
- Basic senescence rate for the green reproductive compartment (K_GR) [0.001]
- Basic abscission rate for the dead vegetative compartment (Kl_DV) [0.001]
- Basic abscission rate for the dead reproductive compartment (Kl_DR) [0.0005]
- Organic matter digestibility for the dead vegetative compartment (OMD_DV)
  [0.45]
- Organic matter digestibility for the dead reproductive compartment (OMD_DR)
  [0.40]
- Maximum radiation use efficiency (RUE_max) [3 g DM MJ⁻¹]
"""

import climag.modvege_lib as lm
import climag.modvege_consumption as cm


def modvege(params, tseries, endday=365) -> dict[str, float]:
    """
    **ModVege** model as a function

    Jouven, M., Carrère, P., and Baumont, R. (2006). 'Model predicting dynamics
    of biomass, structure and digestibility of herbage in managed permanent
    pastures. 1. Model description', Grass and Forage Science, vol. 61, no. 2,
    pp. 112-124. DOI: 10.1111/j.1365-2494.2006.00515.x.

    Parameters
    ----------
    params : Parameters (constants)
    tseries : Time series meteorological data
    endday : Number of days of the year (default is 365)

    Returns
    -------
    - Green vegetative biomass [kg DM ha⁻¹]
    - Dead vegetative biomass [kg DM ha⁻¹]
    - Green reproductive biomass [kg DM ha⁻¹]
    - Dead reproductive biomass [kg DM ha⁻¹]
    - Total standing biomass [kg DM ha⁻¹]
    - Potential growth [kg DM ha⁻¹]
    - Total growth [kg DM ha⁻¹]
    - Ingested biomass [kg DM ha⁻¹]
    - Harvested biomass [kg DM ha⁻¹]
    """

    params["sr"] = cm.stocking_rate(params=params)

    # dictionary of outputs
    outputs_dict = {
        "time": [],
        "bm_gv": [],
        "bm_gr": [],
        "bm_dv": [],
        "bm_dr": [],
        # "age_gv": [],
        # "age_gr": [],
        # "age_dv": [],
        # "age_dr": [],
        "bm": [],
        "pgro": [],
        "gro": [],
        "i_bm": [],
        "h_bm": [],
        # "st": [],
        # "sea": [],
        "f_t": [],
        "env": [],
        "rep": [],
        "lai": [],
        "aet": [],
        "wr": []
    }

    # dictionary to store intermediate time series values
    ts_vals = {}

    # daily loop
    for i in range(endday):
        # initialise starting parameters
        # standing biomass, biomass age, and water reserves
        # assume that the initial value of WR = WHC
        if i == 0:
            (
                ts_vals["bm_gv"], ts_vals["bm_gr"],
                ts_vals["bm_dv"], ts_vals["bm_dr"],
                ts_vals["age_gv"], ts_vals["age_gr"],
                ts_vals["age_dv"], ts_vals["age_dr"],
                ts_vals["wr"]
            ) = (
                params["bm_gv"], params["bm_gr"],
                params["bm_dv"], params["bm_dr"],
                params["age_gv"], params["age_gr"],
                params["age_dv"], params["age_dr"],
                params["whc"]
            )

        # initialise ingested/harvested biomass and temperature sum
        if tseries["time"][i].dayofyear == 1:
            # reset to zero on the first day of the year
            ts_vals["i_bm"], ts_vals["h_bm"], ts_vals["st"] = 0.0, 0.0, 0.0

        # 10-d moving average temperature (T_m10)
        ts_vals["t_m10"] = lm.ten_day_moving_avg_temperature(
            day=(i + 1), t_ts=tseries["T"]
        )

        # sum of temperatures (ST)
        ts_vals["st"] = lm.sum_of_temperatures(
            params=params, ts_vals=ts_vals,
            t_ts=tseries["T"], day=(i + 1)
        )

        # temperature function (f(T))
        ts_vals["f_t"] = lm.temperature_function(
            ts_vals=ts_vals, params=params
        )

        # seasonal effect (SEA)
        ts_vals["sea"] = lm.seasonal_effect(ts_vals=ts_vals, params=params)

        # leaf area index (LAI)
        ts_vals["lai"] = lm.leaf_area_index(ts_vals=ts_vals, params=params)

        # actual evapotranspiration (AET)
        ts_vals["aet"] = lm.actual_evapotranspiration(
            pet=tseries["PET"][i], ts_vals=ts_vals
        )

        # water reserves (WR)
        ts_vals["wr"] = lm.water_reserves(
            ts_vals=ts_vals, params=params, precipitation=tseries["PP"][i]
        )

        # water stress (W)
        ts_vals["w"] = lm.water_stress(ts_vals=ts_vals, params=params)

        # water stress function (f(W))
        ts_vals["f_w"] = lm.water_stress_function(
            ts_vals=ts_vals, pet=tseries["PET"][i]
        )

        # environmental limitation of growth (ENV)
        ts_vals["env"] = lm.environmental_limitation(
            ts_vals=ts_vals, params=params, par_i=tseries["PAR_i"][i]
        )

        # potential growth (PGRO)
        ts_vals["pgro"] = lm.potential_growth(
            par_i=tseries["PAR_i"][i], ts_vals=ts_vals, params=params
        )

        # total growth (GRO)
        ts_vals["gro"] = lm.total_growth(ts_vals=ts_vals)

        # reproductive function (REP)
        # grazing always takes place during the grazing season if the
        # stocking rate is > 0
        ts_vals["rep"] = lm.reproductive_function(
            params=params, ts_vals=ts_vals
        )

        # allocation to reproductive
        # TO-DO: when to change NI, and by how much?
        # NI        A2R         NI = [0.35 - 1.2] A2R = [0.3 - 1.23]
        # 0.4       0.30769
        # 0.5       0.42307
        # 0.6       0.53846
        # 0.7       0.65384
        # 0.8       0.76923
        # 0.9       0.88461
        # 1.0       1
        # 1.1       1.11538
        # 1.2       1.23076

        # senescence (SEN)
        lm.senescence(
            ts_vals=ts_vals, params=params, temperature=tseries["T"][i]
        )

        # abscission (ABS)
        lm.abscission(
            ts_vals=ts_vals, params=params, temperature=tseries["T"][i]
        )

        # standing biomass (BM) and biomass age (AGE)
        lm.standing_biomass(ts_vals=ts_vals, params=params)
        lm.biomass_age(
            temperature=tseries["T"][i], params=params, ts_vals=ts_vals
        )

        # organic matter digestibility (OMD)
        cm.organic_matter_digestibility(ts_vals=ts_vals, params=params)

        # ingested biomass
        cm.biomass_ingestion(ts_vals=ts_vals, params=params)

        # harvested biomass
        cm.biomass_harvest(ts_vals=ts_vals, params=params)

        # total standing biomass
        ts_vals["bm"] = (
            ts_vals["bm_gv"] + ts_vals["bm_gr"] +
            ts_vals["bm_dv"] + ts_vals["bm_dr"]
        )

        # recover output streams
        outputs_dict["time"].append(tseries["time"][i])
        outputs_dict["bm_gv"].append(ts_vals["bm_gv"])
        outputs_dict["bm_dv"].append(ts_vals["bm_dv"])
        outputs_dict["bm_gr"].append(ts_vals["bm_gr"])
        outputs_dict["bm_dr"].append(ts_vals["bm_dr"])
        outputs_dict["h_bm"].append(ts_vals["h_bm"])
        outputs_dict["i_bm"].append(ts_vals["i_bm"])
        outputs_dict["gro"].append(ts_vals["gro"])
        outputs_dict["bm"].append(ts_vals["bm"])
        # outputs_dict["st"].append(ts_vals["st"])
        # outputs_dict["age_gv"].append(ts_vals["age_gv"])
        # outputs_dict["age_gr"].append(ts_vals["age_gr"])
        # outputs_dict["age_dv"].append(ts_vals["age_dv"])
        # outputs_dict["age_dr"].append(ts_vals["age_dr"])
        outputs_dict["f_t"].append(ts_vals["f_t"])
        # outputs_dict["sea"].append(ts_vals["sea"])
        outputs_dict["lai"].append(ts_vals["lai"])
        outputs_dict["wr"].append(ts_vals["wr"])
        outputs_dict["aet"].append(ts_vals["aet"])
        outputs_dict["env"].append(ts_vals["env"])
        outputs_dict["pgro"].append(ts_vals["pgro"])
        outputs_dict["rep"].append(ts_vals["rep"])

    return outputs_dict
