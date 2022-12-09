"""modvege.py

This is a Python implementation of the ModVege pasture model, a modified
version of the original Java to Python translation by Chemin (2022).
The Java model was provided by Raphael Martin, INRAE UREP Clermont-Ferrand
for the original Python implementation.
The original ModVege pasture model was developed by Jouven et al. (2006).

Chemin, Y. (2022). 'modvege', Python. [Online]. Available at
https://github.com/YannChemin/modvege (Accessed 6 September 2022).

References
----------
- Jouven, M., Carrère, P. and Baumont, R. (2006). 'Model predicting dynamics
  of biomass, structure and digestibility of herbage in managed permanent
  pastures. 1. Model description', Grass and Forage Science, vol. 61, no. 2,
  pp. 112-124. DOI: 10.1111/j.1365-2494.2006.00515.x.
- Jouven, M., Carrère, P. and Baumont, R. (2006). 'Model predicting dynamics
  of biomass, structure and digestibility of herbage in managed permanent
  pastures. 2. Model evaluation', Grass and Forage Science, vol. 61, no. 2,
  pp. 125-133. DOI: 10.1111/j.1365-2494.2006.00517.x.
- Lardy, R., Bellocchi, G. G., Bachelet, B. and Hill, D. (2011). 'Climate
  Change Vulnerability Assessment with Constrained Design of Experiments,
  Using a Model-Driven Engineering Approach', Guimaraes, Portugal, pp.
  354-362. [Online]. Available at
  https://hal.archives-ouvertes.fr/hal-02802688 (Accessed 4 November 2022).

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

Notes
-----
- SEA > 1 indicates above-ground growth stimulation by mobilisation of
  reserves; SEA is equal to minSEA when ST < 200 °C d, then increases and
  reaches maxSEA when (ST₁ - 200) < ST < (ST₁ - 100) (ST = ST₁ at the
  beginning of the reproductive period); during summer, SEA decreases,
  returning to minSEA at ST₂ (ST = ST₂ at the end of the reproductive period);
  minSEA and maxSEA are functional traits, arranged symmetrically around 1:
  (minSEA + maxSEA)/2 = 1.

List of parameters (constants)
------------------------------
See Tables 2 and 3 in Jouven et al. (2006).
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

from dataclasses import dataclass
import climag.modvege_lib as lm
import climag.modvege_consumption as cm


def modvege(params, tseries, endday=365):
    """**ModVege** model as a function

    ! This model cannot regenerate reproductive growth after a cut !

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

    # # senescent biomass for compartments
    # gv_senescent_biomass = 0
    # gr_senescent_biomass = 0

    cut_height = params["cutHeight"]

    is_harvested = bool(cut_height != 0)
    is_grazed = bool(
        params["livestock_units"] != 0 and params["grazing_area"] != 0
    )

    # nitrogen nutritional index (NI)
    # if NI is below 0.35, force it to 0.35 (Bélanger et al. 1994)
    params["NI"] = max(params["NI"], 0.35)

    # outputs
    outputs_dict = {
        "biomass_gv": [],
        "biomass_dv": [],
        "biomass_gr": [],
        "biomass_dr": [],
        "biomass_growth": [],
        "biomass_harvested": [],
        "biomass_ingested": [],
        "biomass_cut_avail": [],
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
        "actual_evapotranspiration": []
    }

    @dataclass
    class StandingBiomass:
        """Standing biomass"""
        gv: float
        gr: float
        dv: float
        dr: float

    @dataclass
    class BiomassAge:
        """Biomass age"""
        gv: float
        gr: float
        dv: float
        dr: float

    # daily loop
    for i in range(endday):
        # initialise standing biomass and biomass age
        if i == 0:
            biomass = StandingBiomass(
                gv=params["W_GV"], gr=params["W_GR"],
                dv=params["W_DV"], dr=params["W_DR"]
            )
            age = BiomassAge(
                gv=params["init_AGE_GV"], gr=params["init_AGE_GR"],
                dv=params["init_AGE_DV"], dr=params["init_AGE_DR"]
            )
        else:
            biomass = StandingBiomass(
                gv=outputs_dict["biomass_gv"][i - 1],
                gr=outputs_dict["biomass_gr"][i - 1],
                dv=outputs_dict["biomass_dv"][i - 1],
                dr=outputs_dict["biomass_dr"][i - 1]
            )
            age = BiomassAge(
                gv=outputs_dict["age_gv"][i - 1],
                gr=outputs_dict["age_gr"][i - 1],
                dv=outputs_dict["age_dv"][i - 1],
                dr=outputs_dict["age_dr"][i - 1]
            )

        # temperature (T)
        temperature = tseries["T"][i]

        # 10-d moving average temperature (T_m10)
        temperature_m10 = lm.TenDayMovingAverageTemperature(
            t_ts=tseries["T"], day=(i + 1)
        )()

        # sum of temperatures (ST)
        if tseries["time"][i].dayofyear == 1:
            # reset sum to zero on the first day of the year
            temperature_sum = lm.SumOfTemperatures(
                t_ts=tseries["T"], day=(i + 1), t_sum=0.0
            )()
        else:
            temperature_sum = lm.SumOfTemperatures(
                t_ts=tseries["T"], day=(i + 1),
                t_sum=outputs_dict["temperature_sum"][i - 1]
            )()

        # temperature function (f(T))
        temperature_fn = lm.TemperatureFunction(t_m10=temperature_m10)()
        outputs_dict["temperature_fn"].append(temperature_fn)

        # seasonal effect (SEA)
        seasonality = lm.SeasonalEffect(
            t_sum=temperature_sum, st_1=params["ST1"], st_2=params["ST2"]
        )()
        outputs_dict["seasonality"].append(seasonality)

        # leaf area index (LAI)
        lai = lm.LeafAreaIndex(bm_gv=biomass.gv)()
        outputs_dict["leaf_area_index"].append(lai)

        # potential evapotranspiration (PET)
        pet = tseries["PET"][i]

        # actual evapotranspiration (AET)
        eta = lm.ActualEvapotranspiration(pet=pet, lai=lai)()
        outputs_dict["actual_evapotranspiration"].append(eta)

        # precipitation (PP)
        precipitation = tseries["PP"][i]

        # water reserves (WR)
        params["WR"] = lm.WaterReserves(
            precipitation=precipitation, wreserves=params["WR"],
            aet=eta, whc=params["WHC"]
        )()

        # water stress (W)
        waterstress = lm.WaterStress(
            wreserves=params["WR"], whc=params["WHC"]
        )()

        # water stress function (f(W))
        waterstress_fn = lm.WaterStressFunction(wstress=waterstress, pet=pet)()

        # incident photosynthetically active radiation (PAR_i)
        pari = tseries["PAR"][i]

        # environmental limitation of growth (ENV)
        env = lm.EnvironmentalLimitation(
            t_fn=temperature_fn,
            par_i=pari,
            n_index=params["NI"],
            w_fn=waterstress_fn
        )()
        outputs_dict["env"].append(env)

        # potential growth (PGRO)
        biomass_growth_pot = lm.PotentialGrowth(par_i=pari, lai=lai)()
        outputs_dict["biomass_growth_pot"].append(biomass_growth_pot)

        # total biomass growth (GRO)
        gro = lm.TotalGrowth(
            pgro=biomass_growth_pot, env=env, sea=seasonality
        )()

        # reproductive function (REP)
        # grazing always takes place during the grazing season if the
        # stocking rate is > 0
        # if params["ST2"] > temperature_sum > params["ST1"] and is_grazed:
        #     rep_f = 0
        # else:
        #     rep_f = lm.ReproductiveFunction(n_index=params["NI"])()

        rep_f = lm.ReproductiveFunction(
            n_index=params["NI"], t_sum=temperature_sum,
            st_1=params["ST1"], st_2=params["ST2"]
        )()
        outputs_dict["reproductive_fn"].append(rep_f)

        # senescence (SEN) and abscission (ABS)
        gv_senescent_biomass = lm.SenescenceGV(
            temperature=temperature, age_gv=age.gv, bm_gv=biomass.gv
        )()
        gr_senescent_biomass = lm.SenescenceGR(
            temperature=temperature, age_gr=age.gr, bm_gr=biomass.gr,
            st_1=params["ST1"], st_2=params["ST2"]
        )()
        dv_abscission_biomass = lm.AbscissionDV(
            temperature=temperature, bm_dv=biomass.dv, age_dv=age.dv
        )()
        dr_abscission_biomass = lm.AbscissionDR(
            temperature=temperature, bm_dr=biomass.dr, age_dr=age.dr,
            st_1=params["ST1"], st_2=params["ST2"]
        )()

        # standing biomass (BM) and biomass age (AGE)
        biomass.gv, age.gv = lm.BiomassGV(
            gro_gv=lm.GrowthGV(gro=gro, rep=rep_f)(),
            sen_gv=gv_senescent_biomass,
            bm_gv=biomass.gv,
            age_gv=age.gv,
            temperature=temperature
        )()

        biomass.gr, age.gr = lm.BiomassGR(
            gro_gr=lm.GrowthGR(gro=gro, rep=rep_f)(),
            sen_gr=gr_senescent_biomass,
            bm_gr=biomass.gr,
            age_gr=age.gr,
            temperature=temperature
        )()

        biomass.dv, age.dv = lm.BiomassDV(
            bm_dv=biomass.dv,
            abs_dv=dv_abscission_biomass,
            sen_gv=gv_senescent_biomass,
            age_dv=age.dv,
            temperature=temperature
        )()

        biomass.dr, age.dr = lm.BiomassDR(
            bm_dr=biomass.dr,
            abs_dr=dr_abscission_biomass,
            sen_gr=gr_senescent_biomass,
            age_dr=age.dr,
            temperature=temperature
        )()

        # compute grazing and harvesting
        harvested_biomass_part = 0
        ingested_biomass_part = 0

        if is_grazed:
            # organic matter digestibility
            omd_gv = cm.OrganicMatterDigestibilityGV(age_gv=age.gv)()
            omd_gr = cm.OrganicMatterDigestibilityGR(age_gr=age.gr)()

            # available biomass per compartment
            bm_gv_max = cm.MaximumAvailableBiomass(
                bulk_density=params["rho_GV"], standing_biomass=biomass.gv
            )()
            bm_gr_max = cm.MaximumAvailableBiomass(
                bulk_density=params["rho_GR"], standing_biomass=biomass.gr
            )()
            bm_dv_max = cm.MaximumAvailableBiomass(
                bulk_density=params["rho_DV"], standing_biomass=biomass.dv
            )()
            bm_dr_max = cm.MaximumAvailableBiomass(
                bulk_density=params["rho_DR"], standing_biomass=biomass.dr
            )()

            # stocking rate and max ingestion
            stocking_rate = cm.StockingRate(
                livestock_units=params["livestock_units"],
                grazing_area=params["grazing_area"]
            )()
            bm_ing_max = cm.MaximumIngestedBiomass(
                stocking_rate=stocking_rate
            )()

            # actual ingestion
            ingestion = cm.Ingestion(
                bm_gv_av=bm_gv_max[0], bm_gr_av=bm_gr_max[0],
                bm_dv_av=bm_dv_max[0], bm_dr_av=bm_dr_max[0],
                max_ingested_biomass=bm_ing_max,
                omd_gv=omd_gv, omd_gr=omd_gr
            )()

            # total ingestion
            ingested_biomass_part += (
                ingestion["gv"] + ingestion["gr"] +
                ingestion["dv"] + ingestion["dr"]
            )

            # update biomass compartments
            biomass.gv -= ingestion["gv"]
            biomass.gr -= ingestion["gr"]
            biomass.dv -= ingestion["dv"]
            biomass.dr -= ingestion["dr"]

        # is_harvested = bool(cut_height != 0)
        # is_grazed = bool(
        #     params["livestock_units"] != 0 and params["grazing_area"] != 0
        # )

        #     # look for flags to indicate mechanical cut
        #     if is_harvested:
        #         harvested_biomass_part = [0, 0, 0, 0]
        #         harvested_biomass_part[0], biomass.gv = lm.harvest_biomass(
        #             cutHeight=cut_height, bulkDensity=params["rho_GV"],
        #             biomass=biomass.gv
        #         )
        #         harvested_biomass_part[1], biomass.dv = lm.harvest_biomass(
        #             cutHeight=cut_height, bulkDensity=params["rho_DV"],
        #             biomass=biomass.dv
        #         )
        #         harvested_biomass_part[2], biomass.gr = lm.harvest_biomass(
        #             cutHeight=cut_height, bulkDensity=params["rho_GR"],
        #             biomass=biomass.gr
        #         )
        #         harvested_biomass_part[3], biomass.dr = lm.harvest_biomass(
        #             cutHeight=cut_height, bulkDensity=params["rho_DR"],
        #             biomass=biomass.dr
        #         )
        #         harvested_biomass_part = sum(harvested_biomass_part)
        #         # (
        #         #     harvested_biomass_part,
        #         #     biomass.gv, biomass.dv, biomass.gr, biomass.dr
        #         # ) = (
        #         #     lm.cut(
        #         #         cutHeight=cut_height, rhogv=params["rho_GV"],
        #         #         rhodv=params["rho_DV"], rhogr=params["rho_GR"],
        #         #         rhodr=params["rho_DR"], gvb=biomass.gv,
        #         #         dvb=biomass.dv, grb=biomass.gr, drb=biomass.dr
        #         #     )
        #         # )

        #     # look for flags to indicate livestock ingestion
        #     if is_grazed:
        #         # ingested_biomass_part = lm.defoliation(
        #         #     biomass.gv=biomass.gv, biomass.dv=biomass.dv,
        #         #     biomass.gr=biomass.gr, biomass.dr=biomass.dr,
        #         #     cutHeight=cut_height, rhogv=params["rho_GV"],
        #         #     rhodv=params["rho_DV"], rhogr=params["rho_GR"],
        #         #     rhodr=params["rho_DR"]
        #         # )  # ** MODIFIED -- NEED TO CHECK!
        #         # ingested biomass based on stocking rate
        #         # for bd in ["rho_GV", "rho_GR", "rho_DV", "rho_DR"]:
        #         #     ingested_biomass_part += lm.ingested_biomass(
        #         #         livestock_units=params["livestock_units"],
        #         #         grazing_area=params["grazing_area"],
        #         #         bulk_density=params[bd]
        #         #     )
        #         # ingested_biomass_part = [0, 0, 0, 0]
        #         ingested_biomass_part = lm.ingested_biomass(
        #             livestock_units=params["livestock_units"],
        #             grazing_area=params["grazing_area"],
        #             bulk_density=params["rho_GV"],
        #             min_cut_height=params["cutHeight"]
        #         )
        #         # ingested_biomass_part = sum(ingested_biomass_part)
        #     # allocation to reproductive
        #     a2r = lm.ReproductiveFunction(n_index=params["NI"])()
        #     # TO-DO: when to change NI, and by how much?
        #     # NI        A2R         NI = [0.35 - 1.2] A2R = [0.3 - 1.23]
        #     # 0.4       0.30769
        #     # 0.5       0.42307
        #     # 0.6       0.53846
        #     # 0.7       0.65384
        #     # 0.8       0.76923
        #     # 0.9       0.88461
        #     # 1.0       1
        #     # 1.1       1.11538
        #     # 1.2       1.23076
        # else:
        #     # if (temperature_sum < st1 or st2 < temperature_sum)
        #     a2r = 0
        #     is_harvested = False
        #     is_grazed = False

        # # isCut = bool(is_grazed is True or is_harvested is True)

        # if bool(is_grazed is True or is_harvested is True):
        #     # permanently stop reproduction
        #     a2r = 0

        # gro_gr = lm.GrowthGR(gro=gro, rep=rep_f)
        # biomass.gv, age.gv = lm.BiomassGV(
        #     gro_gv=lm.GrowthGV(gro=gro, rep=rep_f),
        #     sen_gv=lm.SenescenceGV(temperature, age_gv, bm_gv), bm_gv,
        #     age_gv, temperature
        # )

        # # update the state of the vegetative parts
        # biomass.gv, age.gv, gv_senescent_biomass = lm.gv_update(
        #     gro=gro, a2r=rep_f, temperature=temperature,
        #     t0=params["T0"], biomass.gv=biomass.gv,
        #     age.gv=age.gv
        # )
        # biomass.dv, age.dv = lm.dv_update(
        #     gv_gamma=params["sigmaGV"],
        #     gv_senescent_biomass=gv_senescent_biomass,
        #     temperature=temperature,
        #     biomass.dv=biomass.dv, age.dv=age.dv
        # )

        # # start the reproductive phase of the vegetation
        # biomass.gr, age.gr, gr_senescent_biomass = lm.gr_update(
        #     temperature=temperature, a2r=a2r, gro=gro, t0=params["T0"],
        #     biomass.gr=biomass.gr, age.gr=age.gr
        # )  # ** UNUSED ARGUMENTS REMOVED!
        # biomass.dr, age.dr = lm.dr_update(
        #     gr_gamma=params["sigmaGR"],
        #     gr_senescent_biomass=gr_senescent_biomass,
        #     temperature=temperature,
        #     biomass.dr=biomass.dr, age.dr=age.dr
        # )

        # # compute available biomass for cut (output comparison requirement)
        # biomass_cut_avail = lm.getAvailableBiomassForCut(
        #     biomass.gv=biomass.gv, biomass.dv=biomass.dv,
        #     biomass.gr=biomass.gr, biomass.dr=biomass.dr,
        #     cutHeight=cut_height,
        #     rhogv=params["rho_GV"], rhodv=params["rho_DV"],
        #     rhogr=params["rho_GR"], rhodr=params["rho_DR"]
        # )

        # # accumulate harvestedBiomass
        # if is_harvested:
        #     # harvestedBiomass += outputs_dict["biomass_cut_avail"][-1]
        #     harvestedBiomass += harvested_biomass_part
        # # Accumulate ingestedBiomass
        # if is_grazed:
        #     ingestedBiomass += ingested_biomass_part

        biomass_cut_avail = 0

        # Recover output streams
        outputs_dict["biomass_gv"].append(biomass.gv)
        outputs_dict["biomass_dv"].append(biomass.dv)
        outputs_dict["biomass_gr"].append(biomass.gr)
        outputs_dict["biomass_dr"].append(biomass.dr)
        outputs_dict["biomass_harvested"].append(harvested_biomass_part)
        outputs_dict["biomass_ingested"].append(ingested_biomass_part)
        outputs_dict["biomass_growth"].append(gro)
        outputs_dict["biomass_cut_avail"].append(biomass_cut_avail)
        outputs_dict["temperature_sum"].append(temperature_sum)
        outputs_dict["age_gv"].append(age.gv)
        outputs_dict["age_gr"].append(age.gr)
        outputs_dict["age_dv"].append(age.dv)
        outputs_dict["age_dr"].append(age.dr)

    return (
        outputs_dict["biomass_gv"], outputs_dict["biomass_dv"],
        outputs_dict["biomass_gr"], outputs_dict["biomass_dr"],
        outputs_dict["biomass_harvested"], outputs_dict["biomass_ingested"],
        outputs_dict["biomass_growth"], outputs_dict["biomass_cut_avail"],
        outputs_dict["temperature_sum"], outputs_dict["age_gv"],
        outputs_dict["age_dv"], outputs_dict["age_gr"], outputs_dict["age_dr"],
        outputs_dict["seasonality"], outputs_dict["temperature_fn"],
        outputs_dict["env"], outputs_dict["biomass_growth_pot"],
        outputs_dict["reproductive_fn"],
        outputs_dict["leaf_area_index"],
        outputs_dict["actual_evapotranspiration"]
    )
