"""modvege.py

This is the Python implementation of the ModVege pasture model, translated
from Java to Python by Chemin (2022).
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

List of variables
-----------------
- Compartments representing the structural components of herbage
    - Green vegetative (GV) - green leaves and sheath
    - Dead vegetative (DV) - dead leaves and sheath
    - Green reproductive (GR) - green stems and flowers
    - Dead reproductive (DR) - dead stems and flowers
- Compartment description
    - Standing biomass (BM) [kg DM ha⁻¹]
    - Age of biomass (AGE) [°C d] (degree day - units of thermal time)
    - Organic matter digestibility (OMD)
- Quantities
    - Mean daily temperature (*T*) [°C]
    - Sum of temperatures (ST) [°C d]
    - Seasonal effect (SEA)
    - Specific leaf area (SLA)
    - Percentage of laminae in the green vegetative compartment (%LAM)
    - Leaf lifespan (LLS) [°C d]

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
"""

import numpy as np
# import ModVege libraries
import climag.modvege_lib as lm


def modvege(params, weather, startdoy, enddoy, default_cut_height=0.05):
    """**ModVege** model as a function

    ! This model cannot regenerate reproductive growth after a cut !

    Jouven, M., Carrère, P. and Baumont, R. (2006). 'Model predicting dynamics
    of biomass, structure and digestibility of herbage in managed permanent
    pastures. 1. Model description', Grass and Forage Science, vol. 61, no. 2,
    pp. 112-124. DOI: 10.1111/j.1365-2494.2006.00515.x.

    Parameters
    ----------
    params : parameters of this run
    weather : weather data (and grass cut and grazing)
    startdoy : day of year when the simulation starts
    enddoy : day of year when simulation stops
    default_cut_height : defined as 0.05

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
    #######################################################
    # Load input parameters into variables
    #######################################################
    # Onset of reproductive growth (degree day)
    # st1 = params["ST1"]
    # End of reproductive growth (degree day)
    # st2 = params["ST2"]
    # Initial nutritional index of cell
    # ni = params["NI"]
    # Soil water-holding capacity (mm)
    # waterHoldingCapacity = params["WHC"]
    # Soil water reserve (mm)
    # waterReserve = params["WR"]
    # Growth increase in winter
    # minsea = params["minSEA"]
    # Growth increase in summer
    # maxsea = params["maxSEA"]
    # Biomass of GV (kg DM ha-1)
    # wgv = params["W_GV"]
    # Light Use Interception
    # alphapar = params["alpha_PAR"]
    # Temperature threshold: photosynthesis activation (degC)
    # t0 = params["T0"]
    # Temp threshold: stable growth (degC)
    # t1 = params["T1"]
    # Temp threshold: growth decline (degC)
    # t2 = params["T2"]
    # beta_T
    # betaT = params["beta_T"]
    # b_IN
    # b_IN = params["b_IN"]
    # Specific leaf area (m2 g-1)
    # sla = params["SLA"]
    # Leaf lifespan (degree day)
    # lls = params["LLS"]
    # Volume GV (g m-3)
    # rhogv = params["rho_GV"]
    # % leaf of laminae in GV
    # pctlam = params["pctLAM"]
    # Biomass of GR (kg ha-1)
    # wgr = params["W_GR"]
    # Value of ALLOC at NI=0
    # allocNI = params["a_IN"]
    # max of fNI
    # maxFNI = params["max_fIN"]
    # Volume GR (g m-3)
    # rhogr = params["rho_GR"]
    # Biomass of DV (kg ha-1)
    # wdv = params["W_DV"]
    # Senescence coefficient GV
    # kgv = params["K_GV"]
    # Abscission coefficient DV
    # kldv = params["Kl_DV"]
    # Volume DV (g m-3)
    # rhodv = params["rho_DV"]
    # Biomass of DR (kg ha-1)
    # wdr = params["W_DR"]
    # Senescence coefficient GR
    # kgr = params["K_GR"]
    # Abscission coefficient DR
    # kldr = params["Kl_DR"]
    # Volume DR (g m-3)
    # rhodr = params["rho_DR"]
    # Initial value of age of compartment GV
    # gv_init_age = params["init_AGE_GV"]
    # Initial value of age of compartment GR
    # gr_init_age = params["init_AGE_GR"]
    # Initial value of age of compartment DV
    # dv_init_age = params["init_AGE_DV"]
    # Initial value of age of compartment DR
    # dr_init_age = params["init_AGE_DR"]
    # Max of R.U.E.
    # ruemax = params["RUEmax"]
    # rates of biomass loss with respiration for GV
    # sigma_gv = params["sigmaGV"]
    # rates of biomass loss with respiration for GR
    # sigma_gr = params["sigmaGR"]
    # maximum OMD green veg
    # maxOMDgv = params["maxOMDgv"]
    # minimum OMD green veg
    # minOMDgv = params["minOMDgv"]
    # maximum OMD green rep
    # maxOMDgr = params["maxOMDgr"]
    # minimum OMD green rep
    # minOMDgr = params["minOMDgr"]
    # mean OMD dead veg
    # meanOMDdv = params["meanOMDdv"]
    # mean OMD dead rep
    # meanOMDdr = params["meanOMDdr"]
    # Pixel area [Ha]
    # cellSurface = params["cellSurface"]

    # Pixel area [m2]
    # cellSurfaceMeter = 10000 * params["cellSurface"]

    # distance minimum between center of the cell and border of the cell
    # r = np.sqrt(6 * np.sqrt(3) * cellSurfaceMeter) / 6

    # Initialise state parameters
    # This is a status flag changed in lib_cell.updateCell()
    isCut = False
    # This is an actionable flag modified by weather.gcut_height presence
    isHarvested = False
    # This is an actionable flag modified by weather.grazing* presences
    isGrazed = False
    # permanently stop reproduction after the first cut (isCut is True)
    # p116, Jouven et al. (2006)
    a2rFlag = False
    # Harvested biomass
    harvestedBiomass = 0
    # Harvested biomass
    ingestedBiomass = 0
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
    for i in range(startdoy, enddoy, 1):
        #######################################################
        # Load additional input arrays into variables
        #######################################################
        # arr[0][0] = DOY[0] = 1
        # arr[0][1] = Temperature[0] = -0.84125
        temperature = weather[i - 1][1]
        # mean Ten Days Temperature
        if i < 10:
            listA = [weather[i - j][1] for j in range(1, 10, 1)]
        else:
            listA = [weather[i - j][1] for j in range(10, 1, -1)]
        meanTenDaysT = np.mean(listA)

        # arr[0][2] = PARi[0] = 2.22092475
        pari = weather[i - 1][2]
        # arr[0][3] = PP[0] = 0.119
        pmm = weather[i - 1][3]
        # arr[0][4] = PET[0] = 0.602689848
        pet = weather[i - 1][4]
        # arr[0][5] = ETA[0] = 0.4 [RS data, optional]
        eta = weather[i - 1][5]
        # arr[0][6] = LAI[0] = 0.02 [RS data, optional]
        lai = weather[i - 1][6]
        # arr[0][7] = gcut_height[0] = 0.0 [default is 0.05 if cut]
        cutHeight = weather[i - 1][7]
        # arr[0][8] = grazing_animal_count[0] = 0 [default is 1 for test]
        grazing_animal_count = weather[i - 1][8]
        # arr[0][9] = grazing_avg_animal_weight[0] = 0 [default is 400 for cow]
        grazing_avg_animal_weight = weather[i - 1][9]
        #######################################################
        # Prepare additional variables
        #######################################################
        # mk sumTemperature Uses T0=0 and not T0
        sumT = lm.getSumTemperature(weather, i, 0.55)
        # fSEA array for graphs
        sea.append(
            lm.fsea(
                params["maxSEA"], params["minSEA"], sumT, params["ST2"],
                params["ST1"]
            )
        )
        # fTemperature the array for graphs (** UNUSED ARGUMENT REMOVED!)
        ftm.append(
            lm.fTemperature(
                meanTenDaysT, params["T0"], params["T1"], params["T2"]
            )
        )

        # Grass cut flag modification if weather file has grass cut for that
        # day
        isHarvested = bool(cutHeight != 0.0)
        # grazing flag modification if weather file has BOTH animal related
        # values
        isGrazed = bool(
            grazing_animal_count != 0 and grazing_avg_animal_weight != 0
        )
        # Reset the flag isCut
        if isGrazed is False and isHarvested is False:
            isCut = False
        # TO-DO: This variable is not found in the manual/sourcecode, yet is
        # used widely
        correctiveFactorForAn = 1
        # If the Nitrogen Nutrition Index (NI) is below 0.35, force it to 0.35
        # (Belanger et al. 1994)
        params["NI"] = max(params["NI"], 0.35)

        #####################################################################
        # The model starts here really
        #####################################################################
        # If ETA from remote sensing not available, then compute it
        if int(eta) == 0:
            # If LAI from remote sensing not available, then compute it
            if int(lai) == 0:
                lai = lm.fclai(params["pctLAM"], params["SLA"], gv_biomass)
            eta = lm.aet(
                pet, params["pctLAM"], params["SLA"], gv_biomass, params["WR"],
                params["WHC"], lai
            )
        # Compute WR
        params["WR"] = min(
            max(0, params["WR"] + pmm - eta), params["WHC"]
        )

        # Compute CUT
        # harvestedBiomassPart = 0
        ingestedBiomassPart = 0
        # Are we in vegetative growth period?
        if params["ST2"] > sumT > params["ST1"]:
            # Look for flags to indicate mechanical cut
            if isHarvested:
                # Change status flag
                isCut = True
                # The Holy Grail: The Holy Hand Grenade:
                # "Thou Shalst Make the CUT!"
                (
                    isHarvested, harvestedBiomassPart,
                    gv_biomass, dv_biomass, gr_biomass, dr_biomass
                ) = lm.cut(
                    cutHeight, params["rho_GV"], params["rho_DV"],
                    params["rho_GR"], params["rho_DR"], gv_biomass,
                    dv_biomass, gr_biomass, dr_biomass, params["cellSurface"],
                    isHarvested
                )

            # Look for flags to indicate livestock ingestion
            if isGrazed:
                # change status flag
                isCut = True
                # The Holy Grail: The Holy Hand Grenade: "Thou Shalst be wary
                # of this henceforth wicked rabbit!"
                ingestedBiomassPart = lm.defoliation(
                    gv_biomass, dv_biomass, gr_biomass, dr_biomass, cutHeight,
                    params["rho_GV"], params["rho_DV"], params["rho_GR"],
                    params["rho_DR"]
                )  # ** MODIFIED -- NEED TO CHECK!
            # allocation to reproductive
            a2r = lm.rep(params["NI"])
            # TO-DO: When to change NI, and by how much?
            # NI        A2R         NI = [0.35 - 1.2] A2R = [0.3 - 1.23]
            # 0.4       0.30769
            # 0.5       0.42307
            # 0.6       0.53846
            # 0.7       0.65384
            # 0.8       0.76923
            # 0.9       0.88461
            # 1         1
            # 1.1       1.11538
            # 1.2       1.23076
        else:
            # if (sumT < st1 or st2 < sumT)
            a2r = 0

        if isCut:
            # permanently stop reproduction
            a2rFlag = True

        if a2rFlag:
            a2r = 0

        atr.append(a2r)
        # compute biomass growth
        env.append(
            lm.mk_env(
                meanTenDaysT, params["T0"], params["T1"], params["T2"],
                params["NI"], pari, params["alpha_PAR"], pet, params["WR"],
                params["WHC"]
            )  # ** UNUSED ARGUMENT REMOVED!
        )
        pgr.append(
            lm.pgro(
                pari, params["RUEmax"], params["pctLAM"], params["SLA"],
                gv_biomass, lai
            )
        )
        gro = (
            lm.mk_env(
                meanTenDaysT, params["T0"], params["T1"], params["T2"],
                params["NI"], pari, params["alpha_PAR"], pet, params["WR"],
                params["WHC"]
            )  # ** UNUSED ARGUMENT REMOVED!
            * lm.pgro(
                pari, params["RUEmax"], params["pctLAM"], params["SLA"],
                gv_biomass, lai
            )
            * lm.fsea(
                params["maxSEA"], params["minSEA"], sumT, params["ST2"],
                params["ST1"]
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
            gro, a2r, params["LLS"], temperature, params["K_GV"], 0,
            gv_biomass, gv_avg_age
        )
        dv_biomass, dv_avg_age = lm.dv_update(
            params["sigmaGV"], gv_senescent_biomass, params["LLS"],
            params["Kl_DV"], temperature, dv_biomass, dv_avg_age
        )
        # Start the reproductive phase of the vegetation
        gr_biomass, gr_avg_age, gr_senescent_biomass = lm.gr_update(
            temperature, a2r, gro, params["ST1"], params["ST2"],
            params["K_GR"], params["T0"], gr_biomass, gr_avg_age
        )  # ** UNUSED ARGUMENTS REMOVED!
        dr_biomass, dr_avg_age = lm.dr_update(
            params["sigmaGR"], gr_senescent_biomass, params["ST1"],
            params["ST2"], temperature, params["Kl_DR"], dr_biomass,
            dr_avg_age
        )
        # If we do not cut the grass, ensure default estimation is created
        if not isCut:
            cutHeight = default_cut_height
        # Compute available biomass for cut (output comparison requirement)
        avBiom4cut = lm.getAvailableBiomassForCut(
            gv_biomass, dv_biomass, gr_biomass, dr_biomass, cutHeight,
            params["rho_GV"], params["rho_DV"], params["rho_GR"],
            params["rho_DR"]
        )

        #####################################################################
        # The model stops here really
        #####################################################################
        # Accumulate harvestedBiomass
        if isCut:
            harvestedBiomass += abc[-1]
        # Accumulate ingestedBiomass
        ingestedBiomass += ingestedBiomassPart
        # Recover output streams
        gvb.append(gv_biomass)
        dvb.append(dv_biomass)
        grb.append(gr_biomass)
        drb.append(dr_biomass)
        hb.append(harvestedBiomass)
        ib.append(ingestedBiomass)
        g.append(gro)
        abc.append(avBiom4cut)
        stp.append(sumT)
        gva.append(gv_avg_age)
        gra.append(gr_avg_age)
        dva.append(dv_avg_age)
        dra.append(dr_avg_age)

    return (
        gvb, dvb, grb, drb, hb, ib, g, abc, stp, gva,
        gra, dva, dra, sea, ftm, env, pgr, atr
    )
