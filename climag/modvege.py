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
- Calanca, P., Deléglise, C., Martin, R., Carrère, P. and Mosimann, E. (2016).
  'Testing the ability of a simple grassland model to simulate the seasonal
  effects of drought on herbage growth', Field Crops Research, vol. 187, pp.
  12-23. DOI: 10.1016/j.fcr.2015.12.008.

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

import numpy as np
import climag.modvege_lib as lm


def modvege(params, tseries, enddoy=365):
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
    enddoy : End day of the year (default 365)

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

    # Pixel area [m2]
    # cellSurfaceMeter = 10000 * params["cellSurface"]

    # distance minimum between center of the cell and border of the cell
    # r = np.sqrt(6 * np.sqrt(3) * cellSurfaceMeter) / 6

    # Initialise state parameters
    # This is a status flag changed in lib_cell.updateCell()
    # isCut = False
    # This is an actionable flag modified by gcut_height presence
    # isHarvested = False
    # This is an actionable flag modified by grazing* presences
    # isGrazed = False
    # permanently stop reproduction after the first cut (isCut is True)
    # p116, Jouven et al. (2006)
    # a2rFlag = False
    # biomass for compartments
    gv_biomass = params["W_GV"]
    gr_biomass = params["W_GR"]
    dv_biomass = params["W_DV"]
    dr_biomass = params["W_DR"]
    # an==gro: biomass growth
    # an = 0.0
    # allocate to reproductive
    a2r = 0.0
    # senescent biomass for compartments
    gv_senescent_biomass = 0.0
    gr_senescent_biomass = 0.0
    # average age of grass
    gv_avg_age = params["init_AGE_GV"]
    dv_avg_age = params["init_AGE_DV"]
    gr_avg_age = params["init_AGE_GR"]
    dr_avg_age = params["init_AGE_DR"]

    # outputs
    # GV biomass
    gvb = []
    # DV biomass
    dvb = []
    # GR biomass
    grb = []
    # DR biomass
    drb = []
    # GRO biomass
    g = []
    # harvested biomass
    hb = []
    # ingested biomass
    ib = []
    # available biomass for cut
    abc = []
    # sum of temperature
    stp = []
    # ages
    gva = []
    gra = []
    dva = []
    dra = []
    # fSEA
    sea = []
    # fTemperature
    ftm = []
    # mk_env
    env = []
    # pGRO
    pgr = []
    # a2r
    atr = []

    # daily loop
    for i in range(enddoy):
        #######################################################
        # Load additional input arrays into variables
        #######################################################
        temperature = tseries["tas"][i]
        # mean ten days temperature
        if i < (10 - 1):
            meanTenDaysT = temperature  # ** USING THE TEMP, NOT 10-d AVG!
        else:
            meanTenDaysT = np.mean(
                [tseries["tas"][i - j] for j in range(10 - 1, 0 - 1, -1)]
            )

        pari = tseries["par"][i]
        pmm = tseries["pr"][i]
        pet = tseries["evspsblpot"][i]
        eta = tseries["eta"][i]
        lai = tseries["lai"][i]
        cutHeight = tseries["gcut_height"][i]
        # use default cut height if time series value is zero
        # see sec. "Harvested biomass" in Jouven et al. (2006)
        if cutHeight == 0.0:
            cutHeight = params["cutHeight"]
        # grazing_animal_count = tseries["grazing_animal_count"][i]
        # grazing_avg_animal_weight = tseries["grazing_avg_animal_weight"][i]

        #######################################################
        # Prepare additional variables
        #######################################################
        # mk sumTemperature uses T0=0 and not T0
        sumT = lm.getSumTemperature(timeseries=tseries, doy=(i + 1), t0=0)
        # fSEA array for graphs
        sea.append(
            lm.fsea(
                maxsea=params["maxSEA"], minsea=params["minSEA"], sumT=sumT,
                st2=params["ST2"], st1=params["ST1"]
            )
        )
        # fTemperature the array for graphs (** UNUSED ARGUMENT REMOVED!)
        ftm.append(
            lm.fTemperature(
                meanTenDaysT=meanTenDaysT, t0=params["T0"], t1=params["T1"],
                t2=params["T2"]
            )
        )

        # grass cut flag modification if time series file has grass cut for
        # that day
        # isHarvested = bool(cutHeight != 0.0)
        # # grazing flag modification if time series file has BOTH animal
        # # related values
        # isGrazed = bool(
        #     params["livestock_units"] != 0 and params["grazing_area"] != 0
        # )
        # reset the flag isCut
        # isCut = bool(isGrazed is True or isHarvested is True)
        # if isGrazed is False and isHarvested is False:
        #     isCut = False
        # else:
        #     isCut = True
        # TO-DO: This variable is not found in the manual/sourcecode, yet is
        # used widely
        correctiveFactorForAn = 1
        # if the Nitrogen Nutrition Index (NI) is below 0.35, force it to 0.35
        # (Belanger et al. 1994)
        params["NI"] = max(params["NI"], 0.35)

        #####################################################################
        # The model starts here really
        #####################################################################
        # if ETA from remote sensing not available, then compute it
        if int(eta) == 0:
            # if LAI from remote sensing not available, then compute it
            if int(lai) == 0:
                lai = lm.fclai(
                    pctlam=params["pctLAM"], sla=params["SLA"],
                    gv_biomass=gv_biomass
                )
            eta = lm.aet(
                pet=pet, pctlam=params["pctLAM"], sla=params["SLA"],
                gv_biomass=gv_biomass, waterReserve=params["WR"],
                waterHoldingCapacity=params["WHC"], lai=lai
            )
        # compute WR
        params["WR"] = min(max(0, params["WR"] + pmm - eta), params["WHC"])

        # compute CUT
        harvestedBiomassPart = 0
        ingestedBiomassPart = 0

        # are we in vegetative growth period?
        if params["ST2"] > sumT > params["ST1"]:
            isHarvested = bool(cutHeight != 0.0)
            isGrazed = bool(
                params["livestock_units"] != 0 and params["grazing_area"] != 0
            )

            # look for flags to indicate mechanical cut
            if isHarvested:
                # The Holy Grail: The Holy Hand Grenade:
                # "Thou Shalst Make the CUT!"
                (
                    harvestedBiomassPart,
                    gv_biomass, dv_biomass, gr_biomass, dr_biomass
                ) = (
                    lm.cut(
                        cutHeight=cutHeight, rhogv=params["rho_GV"],
                        rhodv=params["rho_DV"], rhogr=params["rho_GR"],
                        rhodr=params["rho_DR"], gvb=gv_biomass,
                        dvb=dv_biomass, grb=gr_biomass, drb=dr_biomass
                    )
                )

            # look for flags to indicate livestock ingestion
            if isGrazed:
                # The Holy Grail: The Holy Hand Grenade: "Thou Shalst be wary
                # of this henceforth wicked rabbit!"
                # ingestedBiomassPart = lm.defoliation(
                #     gv_biomass=gv_biomass, dv_biomass=dv_biomass,
                #     gr_biomass=gr_biomass, dr_biomass=dr_biomass,
                #     cutHeight=cutHeight, rhogv=params["rho_GV"],
                #     rhodv=params["rho_DV"], rhogr=params["rho_GR"],
                #     rhodr=params["rho_DR"]
                # )  # ** MODIFIED -- NEED TO CHECK!
                # ingested biomass based on stocking rate
                # for bd in ["rho_GV", "rho_GR", "rho_DV", "rho_DR"]:
                #     ingestedBiomassPart += lm.ingested_biomass(
                #         livestock_units=params["livestock_units"],
                #         grazing_area=params["grazing_area"],
                #         bulk_density=params[bd]
                #     )
                ingestedBiomassPart = lm.ingested_biomass(
                    livestock_units=params["livestock_units"],
                    grazing_area=params["grazing_area"],
                    bulk_density=params["rho_GV"]
                )
            # allocation to reproductive
            a2r = lm.rep(ni=params["NI"])
            # TO-DO: When to change NI, and by how much?
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
        else:
            # if (sumT < st1 or st2 < sumT)
            a2r = 0
            isHarvested = False
            isGrazed = False

        # isCut = bool(isGrazed is True or isHarvested is True)

        if bool(isGrazed is True or isHarvested is True):
            # permanently stop reproduction
            a2r = 0

        atr.append(a2r)

        # compute biomass growth
        env.append(
            lm.mk_env(
                meanTenDaysT=meanTenDaysT, t0=params["T0"], t1=params["T1"],
                t2=params["T2"], ni=params["NI"], pari=pari,
                alphapar=params["alpha_PAR"], pet=pet,
                waterReserve=params["WR"], waterHoldingCapacity=params["WHC"]
            )  # ** UNUSED ARGUMENT REMOVED!
        )
        pgr.append(
            lm.pgro(
                pari=pari, ruemax=params["RUEmax"], pctlam=params["pctLAM"],
                sla=params["SLA"], gv_biomass=gv_biomass, lai=lai
            )
        )
        gro = (
            lm.mk_env(
                meanTenDaysT=meanTenDaysT, t0=params["T0"], t1=params["T1"],
                t2=params["T2"], ni=params["NI"], pari=pari,
                alphapar=params["alpha_PAR"], pet=pet,
                waterReserve=params["WR"], waterHoldingCapacity=params["WHC"]
            )  # ** UNUSED ARGUMENT REMOVED!
            * lm.pgro(
                pari=pari, ruemax=params["RUEmax"], pctlam=params["pctLAM"],
                sla=params["SLA"], gv_biomass=gv_biomass, lai=lai
            )
            * lm.fsea(
                maxsea=params["maxSEA"], minsea=params["minSEA"], sumT=sumT,
                st2=params["ST2"], st1=params["ST1"]
            )
            * correctiveFactorForAn
        )

        # egro = lm.mk_env(
        #     meanTenDaysT, params["T0"], params["T1"], params["T2"], sumT,
        #     params["NI"], pari, params["alpha_PAR"], pet, params["WR"],
        #     params["WHC"]
        # )
        # ggro = lm.pgro(
        #     pari, params["RUEmax"], params["pctLAM"], params["SLA"],
        #     gv_biomass, lai
        # )
        # sgro = lm.fsea(
        #     params["maxSEA"], params["minSEA"], sumT, params["ST2"],
        #     params["ST1"]
        # )
        # cgro = correctiveFactorForAn
        # gro = egro * ggro * sgro * cgro

        # Update the state of the vegetative parts
        # Used T0 = 0 instead of T0 to match output data!
        gv_biomass, gv_avg_age, gv_senescent_biomass = lm.gv_update(
            gro=gro, a2r=a2r, lls=params["LLS"], temperature=temperature,
            kgv=params["K_GV"], t0=0, gv_biomass=gv_biomass,
            gv_avg_age=gv_avg_age
        )
        dv_biomass, dv_avg_age = lm.dv_update(
            gv_gamma=params["sigmaGV"],
            gv_senescent_biomass=gv_senescent_biomass, lls=params["LLS"],
            kldv=params["Kl_DV"], temperature=temperature,
            dv_biomass=dv_biomass, dv_avg_age=dv_avg_age
        )
        # Start the reproductive phase of the vegetation
        gr_biomass, gr_avg_age, gr_senescent_biomass = lm.gr_update(
            temperature=temperature, a2r=a2r, gro=gro, st1=params["ST1"],
            st2=params["ST2"], kgr=params["K_GR"], t0=params["T0"],
            gr_biomass=gr_biomass, gr_avg_age=gr_avg_age
        )  # ** UNUSED ARGUMENTS REMOVED!
        dr_biomass, dr_avg_age = lm.dr_update(
            gr_gamma=params["sigmaGR"],
            gr_senescent_biomass=gr_senescent_biomass, st1=params["ST1"],
            st2=params["ST2"], temperature=temperature, kldr=params["Kl_DR"],
            dr_biomass=dr_biomass, dr_avg_age=dr_avg_age
        )

        # Compute available biomass for cut (output comparison requirement)
        avBiom4cut = lm.getAvailableBiomassForCut(
            gv_biomass=gv_biomass, dv_biomass=dv_biomass,
            gr_biomass=gr_biomass, dr_biomass=dr_biomass, cutHeight=cutHeight,
            rhogv=params["rho_GV"], rhodv=params["rho_DV"],
            rhogr=params["rho_GR"], rhodr=params["rho_DR"]
        )

        #####################################################################
        # The model stops here really
        #####################################################################
        # # Accumulate harvestedBiomass
        # if isHarvested:
        #     # harvestedBiomass += abc[-1]
        #     harvestedBiomass += harvestedBiomassPart
        # # Accumulate ingestedBiomass
        # if isGrazed:
        #     ingestedBiomass += ingestedBiomassPart
        # Recover output streams
        gvb.append(gv_biomass)
        dvb.append(dv_biomass)
        grb.append(gr_biomass)
        drb.append(dr_biomass)
        hb.append(harvestedBiomassPart)
        ib.append(ingestedBiomassPart)
        g.append(gro)
        abc.append(avBiom4cut)
        stp.append(sumT)
        gva.append(gv_avg_age)
        gra.append(gr_avg_age)
        dva.append(dv_avg_age)
        dra.append(dr_avg_age)

    return (
        gvb, dvb, grb, drb, hb, ib, g, abc, stp, gva, dva, gra, dra, sea, ftm,
        env, pgr, atr
    )
