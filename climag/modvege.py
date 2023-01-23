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
- tas: Mean daily temperature (*T*) [°C]
- par: Incident photosynthetically active radiation (PAR_i) [MJ m⁻²]
- pr: Precipitation (PP) [mm]
- evspsblpot: Potential evapotranspiration (PET) [mm]

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


def modvege(params, tseries, endday=365):
    """**ModVege** model as a function

    Jouven, M., Carrère, P. and Baumont, R. (2006). 'Model predicting dynamics
    of biomass, structure and digestibility of herbage in managed permanent
    pastures. 1. Model description', Grass and Forage Science, vol. 61, no. 2,
    pp. 112-124. DOI: 10.1111/j.1365-2494.2006.00515.x.

    Parameters
    ----------
    params : Parameters (constants)
    tseries : Time series data (weather, grass cut, grazing)
    endday : Number of days of the year (default 365)

    Returns
    -------
    - Green vegetative biomass [kg DM ha⁻¹]
    - Dead vegetative biomass [kg DM ha⁻¹]
    - Green reproductive biomass [kg DM ha⁻¹]
    - Dead reproductive biomass [kg DM ha⁻¹]
    - Harvested biomass [kg DM ha⁻¹]
    - Ingested biomass [kg DM ha⁻¹]
    - GRO biomass [kg DM ha⁻¹]
    - Available biomass for [kg DM ha⁻¹]
    """

    params["sr"] = cm.stocking_rate(params=params)

    # # nitrogen nutritional index (NI)
    # params["ni"] = lm.nitrogen_index(params=params)

    # dictionary of outputs
    outputs_dict = {
        "bm_gv": [],
        "bm_dv": [],
        "bm_gr": [],
        "bm_dr": [],
        "biomass_growth": [],
        "biomass_harvested": [],
        "biomass_ingested": [],
        "biomass_available": [],
        "temperature_sum": [],
        "age_gv": [],
        "age_gr": [],
        "age_dv": [],
        "age_dr": [],
        "seasonality": [],
        "temperature_fn": [],
        "env": [],
        "biomass_growth_pot": [],
        "reproductive_fn": [],
        "leaf_area_index": [],
        "actual_evapotranspiration": [],
        "water_reserves": []
    }

    # dictionary to store intermediate time series values
    ts_vals = {}

    # daily loop
    for i in range(endday):
        # initialise starting parameters
        # standing biomass, biomass age, and water reserves
        if i == 0:
            (
                ts_vals["bm_gv"], ts_vals["bm_gr"],
                ts_vals["bm_dv"], ts_vals["bm_dr"]
            ) = (
                params["bm_gv"], params["bm_gr"],
                params["bm_dv"], params["bm_dr"]
            )

            (
                ts_vals["age_gv"], ts_vals["age_gr"],
                ts_vals["age_dv"], ts_vals["age_dr"]
            ) = (
                params["age_gv"], params["age_gr"],
                params["age_dv"], params["age_dr"]
            )

            # assume that the initial value of WR = WHC
            ts_vals["wr"] = params["whc"]
        else:
            (
                ts_vals["bm_gv"], ts_vals["bm_gr"],
                ts_vals["bm_dv"], ts_vals["bm_dr"]
            ) = (
                outputs_dict["bm_gv"][i - 1],
                outputs_dict["bm_gr"][i - 1],
                outputs_dict["bm_dv"][i - 1],
                outputs_dict["bm_dr"][i - 1]
            )

            (
                ts_vals["age_gv"], ts_vals["age_gr"],
                ts_vals["age_dv"], ts_vals["age_dr"]
            ) = (
                outputs_dict["age_gv"][i - 1],
                outputs_dict["age_gr"][i - 1],
                outputs_dict["age_dv"][i - 1],
                outputs_dict["age_dr"][i - 1]
            )

            ts_vals["wr"] = outputs_dict["water_reserves"][i - 1]

        # initialise ingested/harvested biomass and temperature sum
        if tseries["time"][i].dayofyear == 1:
            # reset to zero on the first day of the year
            ts_vals["i_bm"] = 0.0
            ts_vals["h_bm"] = 0.0
            ts_vals["t_sum"] = 0.0
        else:
            ts_vals["i_bm"] = outputs_dict["biomass_ingested"][i - 1]
            ts_vals["h_bm"] = outputs_dict["biomass_harvested"][i - 1]
            ts_vals["t_sum"] = outputs_dict["temperature_sum"][i - 1]

        # total standing biomass at the beginning of the day
        biomass_available = (
            ts_vals["bm_gv"] + ts_vals["bm_gr"] +
            ts_vals["bm_dv"] + ts_vals["bm_dr"]
        )

        # temperature (T)
        temperature = tseries["T"][i]

        # 10-d moving average temperature (T_m10)
        temperature_m10 = lm.TenDayMovingAverageTemperature(
            t_ts=tseries["T"], day=(i + 1)
        )()

        # sum of temperatures (ST)
        ts_vals["temperature_sum"] = lm.SumOfTemperatures(
            t_ts=tseries["T"], day=(i + 1), t_sum=ts_vals["t_sum"]
        )()

        # temperature function (f(T))
        temperature_fn = lm.TemperatureFunction(t_m10=temperature_m10)()
        outputs_dict["temperature_fn"].append(temperature_fn)

        # seasonal effect (SEA)
        seasonality = lm.SeasonalEffect(
            t_sum=ts_vals["temperature_sum"],
            st_1=params["st_1"], st_2=params["st_2"],
            min_sea=params["min_sea"], max_sea=params["max_sea"]
        )()
        outputs_dict["seasonality"].append(seasonality)

        # leaf area index (LAI)
        # lai = lm.LeafAreaIndex(bm_gv=ts_vals["bm_gv"])()
        # outputs_dict["leaf_area_index"].append(lai)
        ts_vals["lai"] = lm.leaf_area_index(ts_vals=ts_vals, params=params)
        outputs_dict["leaf_area_index"].append(ts_vals["lai"])

        # potential evapotranspiration (PET)
        # pet = tseries["PET"][i]

        # actual evapotranspiration (AET)
        # eta = lm.ActualEvapotranspiration(pet=pet, lai=ts_vals["lai"])()
        # outputs_dict["actual_evapotranspiration"].append(eta)
        ts_vals["eta"] = lm.actual_evapotranspiration(
            pet=tseries["PET"][i], ts_vals=ts_vals
        )
        outputs_dict["actual_evapotranspiration"].append(ts_vals["eta"])

        # precipitation (PP)
        # precipitation = tseries["PP"][i]

        # water reserves (WR)
        ts_vals["wr"] = lm.WaterReserves(
            precipitation=tseries["PP"][i], w_reserves=ts_vals["wr"],
            aet=ts_vals["eta"], whc=params["whc"]
        )()

        # water stress (W)
        waterstress = lm.WaterStress(
            w_reserves=ts_vals["wr"], whc=params["whc"]
        )()

        # water stress function (f(W))
        waterstress_fn = lm.WaterStressFunction(
            w_stress=waterstress, pet=tseries["PET"][i]
        )()

        # incident photosynthetically active radiation (PAR_i)
        pari = tseries["PAR"][i]

        # environmental limitation of growth (ENV)
        env = lm.EnvironmentalLimitation(
            t_fn=temperature_fn,
            par_i=pari,
            n_index=params["ni"],
            w_fn=waterstress_fn
        )()

        # potential growth (PGRO)
        biomass_growth_pot = lm.PotentialGrowth(
            par_i=pari, lai=ts_vals["lai"]
        )()

        # total biomass growth (GRO)
        gro = lm.TotalGrowth(
            pgro=biomass_growth_pot, env=env, sea=seasonality
        )()

        # reproductive function (REP)
        # grazing always takes place during the grazing season if the
        # stocking rate is > 0
        rep_f = lm.ReproductiveFunction(
            n_index=params["ni"], t_sum=ts_vals["temperature_sum"],
            st_1=params["st_1"], st_2=params["st_2"],
            stocking_rate=params["sr"],
            cut_height=params["h_grass"]
        )()

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

        # senescence (SEN) and abscission (ABS)
        gv_senescent_biomass = lm.SenescenceGV(
            temperature=temperature, age_gv=ts_vals["age_gv"],
            bm_gv=ts_vals["bm_gv"]
        )()
        gr_senescent_biomass = lm.SenescenceGR(
            temperature=temperature, age_gr=ts_vals["age_gr"],
            bm_gr=ts_vals["bm_gr"],
            st_1=params["st_1"], st_2=params["st_2"]
        )()
        dv_abscission_biomass = lm.AbscissionDV(
            temperature=temperature, bm_dv=ts_vals["bm_dv"],
            age_dv=ts_vals["age_dv"]
        )()
        dr_abscission_biomass = lm.AbscissionDR(
            temperature=temperature, bm_dr=ts_vals["bm_dr"],
            age_dr=ts_vals["age_dr"],
            st_1=params["st_1"], st_2=params["st_2"]
        )()

        # standing biomass (BM) and biomass age (AGE)
        ts_vals["bm_gv"], ts_vals["age_gv"] = lm.BiomassGV(
            gro_gv=lm.GrowthGV(gro=gro, rep=rep_f)(),
            sen_gv=gv_senescent_biomass,
            bm_gv=ts_vals["bm_gv"],
            age_gv=ts_vals["age_gv"],
            temperature=temperature
        )()

        ts_vals["bm_gr"], ts_vals["age_gr"] = lm.BiomassGR(
            gro_gr=lm.GrowthGR(gro=gro, rep=rep_f)(),
            sen_gr=gr_senescent_biomass,
            bm_gr=ts_vals["bm_gr"],
            age_gr=ts_vals["age_gr"],
            temperature=temperature
        )()

        ts_vals["bm_dv"], ts_vals["age_dv"] = lm.BiomassDV(
            bm_dv=ts_vals["bm_dv"],
            abs_dv=dv_abscission_biomass,
            sen_gv=gv_senescent_biomass,
            age_dv=ts_vals["age_dv"],
            temperature=temperature
        )()

        ts_vals["bm_dr"], ts_vals["age_dr"] = lm.BiomassDR(
            bm_dr=ts_vals["bm_dr"],
            abs_dr=dr_abscission_biomass,
            sen_gr=gr_senescent_biomass,
            age_dr=ts_vals["age_dr"],
            temperature=temperature
        )()

        # organic matter digestibility (OMD)
        cm.organic_matter_digestibility(ts_vals=ts_vals, params=params)

        # ingested biomass
        cm.biomass_ingestion(ts_vals=ts_vals, params=params)

        # harvested biomass
        cm.biomass_harvest(ts_vals=ts_vals, params=params)

        # recover output streams
        outputs_dict["bm_gv"].append(ts_vals["bm_gv"])
        outputs_dict["bm_dv"].append(ts_vals["bm_dv"])
        outputs_dict["bm_gr"].append(ts_vals["bm_gr"])
        outputs_dict["bm_dr"].append(ts_vals["bm_dr"])
        outputs_dict["biomass_harvested"].append(ts_vals["h_bm"])
        outputs_dict["biomass_ingested"].append(ts_vals["i_bm"])
        outputs_dict["biomass_growth"].append(gro)
        outputs_dict["biomass_available"].append(biomass_available)
        outputs_dict["temperature_sum"].append(ts_vals["temperature_sum"])
        outputs_dict["age_gv"].append(ts_vals["age_gv"])
        outputs_dict["age_gr"].append(ts_vals["age_gr"])
        outputs_dict["age_dv"].append(ts_vals["age_dv"])
        outputs_dict["age_dr"].append(ts_vals["age_dr"])
        outputs_dict["water_reserves"].append(ts_vals["wr"])
        outputs_dict["env"].append(env)
        outputs_dict["biomass_growth_pot"].append(biomass_growth_pot)
        outputs_dict["reproductive_fn"].append(rep_f)

    return (
        outputs_dict["bm_gv"],
        outputs_dict["bm_dv"],
        outputs_dict["bm_gr"],
        outputs_dict["bm_dr"],
        outputs_dict["biomass_harvested"],
        outputs_dict["biomass_ingested"],
        outputs_dict["biomass_growth"],
        outputs_dict["biomass_available"],
        outputs_dict["temperature_sum"],
        outputs_dict["age_gv"],
        outputs_dict["age_dv"],
        outputs_dict["age_gr"],
        outputs_dict["age_dr"],
        outputs_dict["seasonality"],
        outputs_dict["temperature_fn"],
        outputs_dict["env"],
        outputs_dict["biomass_growth_pot"],
        outputs_dict["reproductive_fn"],
        outputs_dict["leaf_area_index"],
        outputs_dict["actual_evapotranspiration"],
        outputs_dict["water_reserves"]
    )
