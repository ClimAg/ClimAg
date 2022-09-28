"""modvege_lib.py

https://github.com/YannChemin/modvege
"""

import numpy as np


# def getAverageHeight(biomass, bulkDensity):
#     """
#     Return the average height based on biomass and bulk density.
#     ** THIS FUNCTION IS UNUSED!

#     Parameters
#     ----------
#     biomass : Biomass (BM) [kg DM ha⁻¹]
#     bulkDensity : Bulk density (BD) [g DM m⁻³]

#     Returns
#     -------
#     - Average height [m]
#     """
#     avg_height = biomass / (bulkDensity * 10)
#     return avg_height


def avDefoliationBiomass(biomass, cutHeight, bulkDensity):
    """
    Estimate biomass for ingestion.

    Parameters
    ----------
    biomass : Biomass [kg DM ha⁻¹]
    cutHeight : Default cut height in cut event [m]
    bulkDensity : Bulk density [g DM m⁻³]

    Returns
    -------
    - Biomass for ingestion [kg DM ha⁻¹]
    """
    biomassAfterCut = cutHeight * bulkDensity * 10
    return max(0, biomass - biomassAfterCut)


# def exeDefoliation(biomass, cut_biomass, area):
#     """
#     Defoliation method.
#     ** THIS FUNCTION IS UNUSED!

#     Parameters
#     ----------
#     biomass : Biomass [kg DM ha⁻¹] ?
#     cut_biomass : Biomass removed by cut [kg DM] ?
#     area : Area studied (i.e. pixel area) [?]

#     Returns
#     -------
#     - Updated biomass after defoliation [kg DM ha⁻¹] ?
#     """
#     biomass = biomass - cut_biomass / area
#     if biomass < 0 | np.isnan(biomass):
#         biomass = 0
#     return biomass


def exeCut(cutHeight, bulkDensity, biomass):
    """
    Realise a cut so that the average height is under cutHeight.
    This height is calculated by using the bulkDensity given in parameter.

    Parameters
    ----------
    cutHeight : Average height after the cut [m]
    bulkDensity : Bulk density [g DM m⁻³]
        (biomass after cut = height * 10 * bulk density)
    biomass : Biomass [kg DM ha⁻¹]

    Returns
    -------
    - Biomass taken [kg DM ha⁻¹]
    """
    biomassAfterCut = cutHeight * bulkDensity * 10
    if biomassAfterCut < biomass:
        takenBiomass = biomass - biomassAfterCut
    else:
        takenBiomass = 0
        biomassAfterCut = biomass
    return (takenBiomass, biomassAfterCut)


def exeDefoliationByBiomass(biomass, biomassToIngest):
    """
    Realise a defoliation by providing the biomass to remove.

    Parameters
    ----------
    biomass : Biomass [kg DM ha⁻¹]
    biomassToIngest : Biomass to remove [kg DM ha⁻¹]

    Returns
    -------
    - Biomass left [kg DM ha⁻¹]
    """
    biomass -= biomassToIngest
    return biomassToIngest


# DEAD VEGETATIVE FUNCTIONS
def mk_dv_abscission(kldv, dv_biomass, temperature, dv_avg_age, lls):
    """
    Compute abscission biomass for the DV compartment.
    See Equation (18) in Jouven et al. (2006).

    Parameters
    ----------
    lls : Leaf lifespan (LLS) [°C d]
    kldv : Abscission coefficient for DV (Kl_DV) (Ducrocq 1996)
    temperature : Mean daily temperature (T) [°C]
    dv_biomass : DV biomass (BM_DV) [kg DM ha⁻¹]
    dv_avg_age : Average DV age (AGE_DV) [°C d]

    Returns
    -------
    - Abscission biomass (ABS_DV) [kg DM ha⁻¹]
    """
    # method to compute the age of DV for computing abcission of DV
    # f(AGE_DV) in Equation (18)
    if dv_avg_age / lls < 1.0 / 3.0:
        age = 1
    elif dv_avg_age / lls < 2.0 / 3.0:
        age = 2
    else:
        age = 3
    # compute the abscission for dead vegetative part
    # abscission only occurs when T > 0
    if temperature > 0:
        abscission_biomass = kldv * dv_biomass * temperature * age
    else:
        abscission_biomass = 0
    return abscission_biomass


def dv_update(
    gv_gamma, gv_senescent_biomass, lls, kldv, temperature,
    dv_biomass, dv_avg_age
):
    """
    Update DV compartment.
    See Equation (3) in Jouven et al. (2006).

    Parameters
    ----------
    gv_gamma : Respiratory C loss during senescence (DV)
        (1 - gv_gamma = dv_gamma)
    gv_senescent_biomass : Senescence of GV compartment
    lls : Leaf lifespan (LLS) [°C d]
    kldv : Abscission coefficient for DV (Kl_DV)
    temperature : Temperature (T) [°C]
    dv_biomass : DV biomass (BM_DV) [kg DM ha⁻¹]
    dv_avg_age : Average DV age (AGE_DV) [°C d]

    Returns
    -------
    - Dead vegetation biomass
    - Average DV age
    """
    abscissionBiomass = mk_dv_abscission(
        kldv, dv_biomass, temperature, dv_avg_age, lls
    )
    dv_biomass -= abscissionBiomass
    # at this point the biomass include cut, ingestion, and abscission, not
    # growth
    growthBiomass = (1.0 - gv_gamma) * gv_senescent_biomass
    if dv_biomass + growthBiomass > 0:
        dv_avg_age = (max(0, temperature) + dv_avg_age) * (
            dv_biomass / (dv_biomass + growthBiomass)
        )
    else:
        dv_avg_age = 0
    dv_biomass += growthBiomass
    return (dv_biomass, dv_avg_age)


# DEAD REPRODUCTIVE FUNCTIONS
def mk_dr_abscission(kldr, dr_biomass, temperature, dr_avg_age, st1, st2):
    """
    Compute abscission biomass for the DR compartment.
    See Equation (18) in Jouven et al. (2006).

    Parameters
    ----------
    kldr : Abscission coefficient for DR (Kl_DR) (Ducrocq 1996)
    dr_biomass : DR biomass (BM_DR) [kg DM ha⁻¹]
    temperature : Mean daily temperature (T) [°C]
    dr_avg_age : Average DR age (AGE_DR) [°C d]
    st1 : Sum of temperature at the beginning (ST_1) [°C d]
    st2 : Sum of temperature in the end (ST_2) [°C d]

    Returns
    -------
    - Abscission biomass for DR (ABS_DR) [kg DM ha⁻¹]
    """
    # method to compute the age of DR for computing abscission of DR
    # f(AGE_DV) in Equation (18)
    if dr_avg_age / (st2 - st1) < 1.0 / 3.0:
        age = 1
    elif dr_avg_age / (st2 - st1) < 2.0 / 3.0:
        age = 2
    else:
        age = 3
    # compute abscission for dead reproductive
    # abscission only occurs when T > 0
    if temperature > 0:
        abscission_biomass = kldr * dr_biomass * temperature * age
    else:
        abscission_biomass = 0
    return abscission_biomass


def dr_update(
    gr_gamma, gr_senescent_biomass, st1, st2,
    temperature, kldr, dr_biomass, dr_avg_age
):
    """
    Update DR compartment.

    Parameters
    ----------
    gr_gamma : Parameter to compute growth of DR (1 - gr_gamma = dr_gamma)
    gr_senescent_biomass : Senescence of GR, computed in GR
    st1 : Sum of temperature at the beginning (ST_1) [°C d]
    st2 : Sum of temperature in the end (ST_2) [°C d]
    temperature : Temperature [°C]
    kldr : Basic rates of abscission in DR
    dr_biomass : DR biomass
    dr_avg_age : Average DR age [°C d]

    Returns
    -------
    - Updated biomass for DR
    - Average DR age
    """
    abscissionBiomass = mk_dr_abscission(
        kldr, dr_biomass, temperature, dr_avg_age, st1, st2
    )
    dr_biomass -= abscissionBiomass
    # at this point the biomass include cut, ingestion and abscission, not
    # growth
    growthBiomass = (1 - gr_gamma) * gr_senescent_biomass
    if dr_biomass + growthBiomass > 0:
        dr_avg_age = max(0, temperature) + dr_avg_age * dr_biomass / (
            dr_biomass + growthBiomass
        )
    else:
        dr_avg_age = 0
    dr_biomass += growthBiomass
    return (dr_biomass, dr_avg_age)


# GREEN VEGETATIVE FUNCTIONS
def mk_gv_senescence(kgv, gv_biomass, temperature, t0, lls, gv_avg_age):
    """
    Extract about 2-6% (kDV=0.002, T=10C, gv_fAge=[1-3]) of gv_biomass as
    senescent.
    See Equations (16) and (17) in Jouven et al. (2006).

    Parameters
    ----------
    kgv : Senescence coefficient (K_GV) (Ducrocq 1996)
    gv_biomass : GV biomass (BM_GV) [kg DM ha⁻¹]
    temperature : Temperature (T) [°C]
    t0 : Minimum temperature for growth (T_0) [°C]
    lls : Leaf lifespan (LLS) [°C d]
    gv_avg_age : Average GV age [°C d]

    Returns
    -------
    - Senescent biomass for GV (SEN_GV) [kg DM ha⁻¹]
    """
    # method to compute the age of GV for computing senescence of GV
    # f(AGE_GV) in Equation (16)
    if gv_avg_age / lls < 1.0 / 3.0:
        age = 1
    elif gv_avg_age / lls < 1:
        age = 3 * gv_avg_age / lls
    else:
        age = 3
    # compute senescence of GV
    if temperature > t0:
        # Equation (16)
        senescence_biomass = kgv * gv_biomass * temperature * age
    elif temperature < 0:
        # Equation (17)
        senescence_biomass = kgv * gv_biomass * abs(temperature)
    else:
        senescence_biomass = 0
    return senescence_biomass


def gv_update(gro, a2r, lls, temperature, kdv, t0, gv_biomass, gv_avg_age):
    """
    Update GV compartment.

    Parameters
    ----------
    gro : in Jouven et al. (2006), total growth (GRO)
    a2r : Allocate to reproductive
        (REP in Jouven et al. (2006), reproductive function)
    lls : Leaf lifespan (LLS) [°C d]
    temperature : Temperature [°C]
    kdv : Senescence coefficient DV [°C d]
    t0 : Minimum temperature for growth (T_0) [°C]
    gv_biomass : Updated biomass
    gv_avg_age : Average GV age [°C day]

    Returns
    -------
    - Senescent biomass
    """
    senescentBiomass = mk_gv_senescence(
        kdv, gv_biomass, temperature, t0, lls, gv_avg_age
    )
    gv_biomass -= senescentBiomass
    # at this point the biomass include cut, ingestion and senescence, not
    # growth
    if temperature > t0:
        growthBiomass = gro * (1 - a2r)
    else:
        growthBiomass = 0
    if gv_biomass + growthBiomass > 0:
        gv_avg_age = (max(0, temperature) + gv_avg_age) * (
            gv_biomass / (gv_biomass + growthBiomass)
        )
    else:
        gv_avg_age = 0
    gv_biomass += growthBiomass
    return (gv_biomass, gv_avg_age, senescentBiomass)


# GREEN REPRODUCTIVE FUNCTIONS
def mk_gr_senescence(
    kgr, gr_biomass, temperature, t0, gr_avg_age, st1, st2
):
    """
    See Equations (16) and (17) in Jouven et al. (2006).

    Parameters
    ----------
    kgr : Senescence coefficient DV [°C d] ** CHECK PARAM NAME!
    gr_biomass : Biomass available for GR (BM_GR) [kg DM ha⁻¹]
    temperature : Temperature (T) [°C]
    t0 : Minimum temperature for growth (T_0) [°C]
    gr_avg_age : Average GR age (AGE_GR) [°C d]
    st1 : Onset of reproductive growth (ST_1) [°C d]
    st2 : End of reproductive growth (ST_2) [°C d]
    lls : Leaf lifespan (LLS) [°C d] (** UNUSED ARGUMENT!)

    Returns
    -------
    - Senescent biomass [kg DM ha⁻¹]
    """
    # method to compute the age of GR for computing senescence of GR
    if gr_avg_age / (st2 - st1) < 1.0 / 3.0:
        age = 1
    elif gr_avg_age / (st2 - st1) < 1.0:
        age = 3 * gr_avg_age / (st2 - st1)
    else:
        age = 3
    # Compute senescence of GV
    if temperature > t0:
        # T=10C kgr = 0.001 gr_fAge=[1-3] => 1-3% of gr_biomass
        senescence_biomass = kgr * gr_biomass * temperature * age
    elif temperature < 0:
        senescence_biomass = kgr * gr_biomass * abs(temperature)
    else:
        senescence_biomass = 0
    return senescence_biomass


def gr_update(
    temperature, a2r, gro, st1, st2, kgr, t0, gr_biomass, gr_avg_age
):
    """
    Update GR compartment.

    Parameters
    ----------
    temperature : Temperature (T) [°C]
    a2r : Allocate to reproductive
        (REP in Jouven et al. (2006), reproductive function)
    gro : Total growth (GRO) [kg DM ha⁻¹]
    st1 : Onset of reproductive growth (ST₁) [°C d]
    st2 : End of reproductive growth (ST₂) [°C d]
    kgr : Basic rates of senescence in compartment GR (K_GR)
    t0 : Minimum temperature for growth (T₀) [°C]
    gr_biomass : GR biomass (BM_GR) [kg DM ha⁻¹]
    gr_avg_age : Average GR age (AGE_GR) [°C d]
    lls : Leaf lifespan (LLS) [°C d] (** UNUSED ARGUMENT!)
    rhogr : Volume GR [g m⁻³] (** UNUSED ARGUMENT!)

    Returns
    -------
    - Updated GR biomass
    - Average GR age
    """
    senescentBiomass = mk_gr_senescence(
        kgr, gr_biomass, temperature, t0, gr_avg_age, st1, st2
    )  # ** UNUSED ARGUMENT REMOVED!
    gr_biomass -= senescentBiomass
    # at this point the biomass include cut, ingestion and senescence, not
    # growth
    if temperature > t0:
        growthBiomass = gro * a2r
    else:
        growthBiomass = 0
    if gr_biomass + growthBiomass > 0:
        gr_avg_age = (max(0, temperature) + gr_avg_age) * (
            gr_biomass / (gr_biomass + growthBiomass)
        )
    else:
        gr_avg_age = 0
    gr_biomass += growthBiomass
    return (gr_biomass, gr_avg_age, senescentBiomass)


#############################################################################
# Nutrition Index                   # double ni
# water reserves WR                 # double waterReserve
# Compartment green vegetative      # CompartmentGreenVegetative cGV
# Compartment green reproductive    # CompartmentGreenReproductive cGR
# Compartment dry vegetative        # CompartmentDryVegetative cDV
# Compartment dry reproductive      # CompartmentDryReproductive cDR
# Indicate if the cell was previously cut during a certain period
#     # boolean isCut
# Cut tag shows if there is a cut event
#     # boolean isHarvested
# Cut tag shows if there is a grazed event
#     # boolean isGrazed
# Sum of GRO                        # double sumGRO
# GRO                               # double gro


# def getHeight(gv_avg_h, gr_avg_h, dv_avg_h, dr_avg_h):
#     """
#     Return the height of this cell
#     ** THIS FUNCTION IS UNUSED!

#     Parameters
#     ----------
#     gv_avg_h : ?
#     gr_avg_h : ?
#     dv_avg_h : ?
#     dr_avg_h : ?

#     Returns
#     -------
#     - the maximum height of the 4 cs
#     """
#     return max(max(gv_avg_h, gr_avg_h), np.max(dv_avg_h, dr_avg_h))


def cut(
    cutHeight, rhogv, rhodv, rhogr, rhodr, gvb, dvb,
    grb, drb, cellSurface, isHarvested
):
    """
    Realise the harvest on each c. If the amount of cut biomass is not null,
    then the flag isHarvested is set to True

    Parameters
    ----------
    cutHeight : the height of the cut [m]
    rhogv : rho green vegetation
    rhodv : rho dry vegetation
    rhogr : rho green reproduction
    rhodr : rho dry reproduction
    gvb : the biomass of green vegetation
    dvb : the biomass of dry vegetation
    grb : the biomass of green reproductive
    drb : the biomass of dry reproductive
    cellSurface : Surface of the pixel [ha]
    isHarvested : Status flag indicating harvest happened

    Returns
    -------
    - Status flag indicating harvest happened
    - the total amount of biomass cut [kg DM]
    - the amount of GV biomass cut [kg DM]
    - the amount of DV biomass cut [kg DM]
    - the amount of GR biomass cut [kg DM]
    - the amount of DR biomass cut [kg DM]
    """
    # exeCut returns harvested biomass in [kg DM m⁻²]
    gv_h, gv_b = exeCut(rhogv, cutHeight, gvb)
    dv_h, dv_b = exeCut(rhodv, cutHeight, dvb)
    gr_h, gr_b = exeCut(rhogr, cutHeight, grb)
    dr_h, dr_b = exeCut(rhodr, cutHeight, drb)
    # sum of harvested biomass [kg DM m⁻²]
    sumBiomassHarvested = gv_h + dv_h + gr_h + dr_h
    if sumBiomassHarvested > 0:
        isHarvested = True
    return (
        isHarvested, sumBiomassHarvested * cellSurface, gv_b, dv_b, gr_b, dr_b
    )


def mk_env(
    meanTenDaysT, t0, t1, t2, ni, pari, alphapar,
    pet, waterReserve, waterHoldingCapacity
):
    """
    Environmental stress.

    Parameters
    ----------
    meanTenDaysT : the mean of the ten days of temperature
    t0 : minimum temperature for growth (T_0) [°C]
    t1 : sum of temperature at the beginning (growth activation threshold)
    t2 : sum of temperature in the end (growth decline threshold)
    ni : Nutritional index of pixel (NI)
    pari : Incident photosynthetically active radiation (PAR_i) [MJ m⁻²]
    alphapar : Light use interception
    pet : Potential evapotranspiration (PET) [mm]
    waterReserve : Water reserves (WR) [mm]
    waterHoldingCapacity : Soil water-holding capacity (WHC) [mm]
    sumT : sum of temperatures (** UNUSED ARGUMENT!)

    Returns
    -------
    - Environmental stress
    """
    return (
        fTemperature(meanTenDaysT, t0, t1, t2)  # ** UNUSED ARGUMENT REMOVED!
        * ni
        * fPARi(pari, alphapar)
        * fWaterStress(waterReserve, waterHoldingCapacity, pet)
    )


# def getTotalBiomass(gv_biomass, dv_biomass, gr_biomass, dr_biomass):
#     """
#     Return the total biomass of the cell (by adding the biomass of the 4 cs)
#     ** THIS FUNCTION IS UNUSED!

#     Parameters
#     ----------
#     gv_biomass : Green vegetation biomass
#     dv_biomass : Dry vegetation biomass
#     gr_biomass : Green reproduction biomass
#     dr_biomass : Dry reproduction biomass

#     Returns
#     -------
#     - Total biomass
#     """
#     return gv_biomass + dv_biomass + gr_biomass + dr_biomass


def fTemperature(meanTenDaysT, t0, t1, t2):
    """
    f of temperature to compute ENV

    Parameters
    ----------
    meanTenDaysT : Mean of the ten days of temperature
    t0 : Minimum temperature for growth
    t1 : Sum of temperature at the beginning (growth activation threshold)
    t2 : Sum of temperature in the end (growth decline threshold)
    sumT : Sum of temperatures (** UNUSED ARGUMENT!)

    Returns
    -------
    - the value given by the temperature f
    """
    if meanTenDaysT < t0 or meanTenDaysT >= 40:
        f_temp = 0
    elif t1 > meanTenDaysT >= t0:
        f_temp = (meanTenDaysT - t0) / (t1 - t0)
    elif t2 > meanTenDaysT >= t1:
        f_temp = 1
    else:
        f_temp = (40 - meanTenDaysT) / (40 - t2)
    return f_temp


def fsea(maxsea, minsea, sumT, st2, st1):
    """
    Function for seasonality (SEA) to compute the potential growth

    Parameters
    ----------
    maxsea : Growth increase in summer (minSEA)
    minsea : Growth increase in winter (maxSEA)
    sumT : Sum of temperature [°C d]
    st1 : Sum of temperature at the beginning of growth (ST_1) [°C d]
    st2 : Sum of temperature in the end of growth (ST_2) [°C d]

    Returns
    -------
    - Value given by the sea f
    """
    if sumT < 200 or sumT >= st2:
        f_sea = minsea
    elif sumT < st1 - 200:
        f_sea = minsea + (maxsea - minsea) * (sumT - 200) / (st1 - 400)
    elif sumT < st1 - 100:
        f_sea = maxsea
    else:
        f_sea = (
            maxsea + (minsea - maxsea) * (sumT - st1 + 100) / (st2 - st1 + 100)
        )
    return f_sea


def fPARi(pari, alphapar):
    """
    Function of PAR interception (PAR_i) to compute ENV

    Parameters
    ----------
    pari : Photosynthetic radiation incident (PAR_i) [MJ m⁻²]
    alphapar : the light use interception

    Returns
    -------
    - the value given by the PAR_i [0-1]
    """
    if pari < 5:
        f_pari = 1
    else:
        f_pari = max(1 - alphapar * (pari - 5), 0)
    return f_pari


def fWaterStress(waterReserve, waterHoldingCapacity, pet):
    """
    f of water stress to compute ENV

    Parameters
    ----------
    waterReserve : reserve of water in the soil
    waterHoldingCapacity : capacity of the soil to hold a certain volume
        of water
    pet : potential evapotranspiration

    Returns
    -------
    - the value given by the waterstress f
    """
    waterStress = min(waterReserve / waterHoldingCapacity, 1)
    if pet <= 3.8:
        if waterStress <= 0.2:
            f_waterstress = 4 * waterStress
        elif waterStress <= 0.4:
            f_waterstress = 0.75 * waterStress + 0.65
        elif waterStress <= 0.6:
            f_waterstress = 0.25 * waterStress + 0.85
        else:
            f_waterstress = 1
    elif pet <= 6.5:
        if waterStress <= 0.2:
            f_waterstress = 2 * waterStress
        elif waterStress <= 0.4:
            f_waterstress = 1.5 * waterStress + 0.1
        elif waterStress <= 0.6:
            f_waterstress = waterStress + 0.3
        elif waterStress <= 0.8:
            f_waterstress = 0.5 * waterStress + 0.6
        else:
            f_waterstress = 1
    else:
        f_waterstress = waterStress
    return f_waterstress


def rep(ni):
    """
    Replace the value of NI (Nutrition Index (Belanger et al. 1992))

    Parameters
    ----------
    ni : The Nutrition Index (Belanger et al. 1992)

    Returns
    -------
    - the value of the rep f
    """
    return 0.25 + ((0.75 * (ni - 0.35)) / 0.65)


def pgro(pari, ruemax, pctlam, sla, gv_biomass, lai):
    """
    Compute and return potential growth

    Parameters
    ----------
    pari : Incident PAR (PAR_i) [MJ m⁻²]
    ruemax : Maximum radiation use efficiency (RUE_max) [3 g DM MJ⁻¹]
    pctlam : Percentage of laminae in GV (%LAM)
    sla : Specific leaf area (SLA) [m² g⁻¹]
    gv_biomass : Green vegetation biomass (BM_GV) [kg DM ha⁻¹]
    lai : the LAI from remote sensing (if available)

    Returns
    -------
    - the calculated pGRO [kg DM ha⁻¹]
    """
    if int(lai) == 0:
        try:
            lai = sla * pctlam * (gv_biomass / 10)
        except (IOError, ValueError):
            # in case of input malfunction
            lai = sla * pctlam * 1.0
    lightInterceptionByPlant = 1 - np.exp(-0.6 * lai)
    p_gro = pari * ruemax * lightInterceptionByPlant * 10
    return p_gro


def fclai(pctlam, sla, gv_biomass):
    """
    Compute and return the leaf area index

    Parameters
    ----------
    pctlam : Percentage of laminae in GV (%LAM)
    sla : Specific leaf area (SLA) [m² g⁻¹]
    gv_biomass : GV biomass (BM_GV) [kg DM ha⁻¹]

    Returns
    -------
    - the calculated LAI
    """
    return sla * (gv_biomass / 10) * pctlam


def aet(pet, pctlam, sla, gv_biomass, waterReserve, waterHoldingCapacity, lai):
    """
    Return the actual evapotranspiration (AET)

    pet : Potential evapotranspiration (PET) [mm]
    pctlam : Percentage of laminae in GV (%LAM)
    sla : Specific leaf area (SLA) [m² g⁻¹]
    gv_biomass : GV biomass
    waterReserve : Water reserves (WR) [mm]
    waterHoldingCapacity : Soil water-holding capacity (WHC) [mm]
    lai : Leaf area index (LAI)

    Returns
    -------
    - Actual evapotranspiration (AET) [mm]
    """
    if int(lai) == 0:
        lai = sla * pctlam * (gv_biomass / 10)
    lightInterceptionByPlant = 1 - np.exp(-0.6 * lai)
    pt = pet * lightInterceptionByPlant
    pe = pet - pt
    ta = pt * fWaterStress(waterReserve, waterHoldingCapacity, pet)
    ea = pe * min(waterReserve / waterHoldingCapacity, 1)
    return ta + ea


# def updateSumTemperature(temperature, t0, sumT, tbase):
#     """
#     Add the daily temperature to the sum temperature if the daily one is
#     positive
#     ** THIS FUNCTION IS UNUSED!

#     Parameters
#     ----------
#     temperature : ?
#     t0 : minimum temperature for growth
#     sumT : actual sum of temperature
#     tbase : base temperature (subtracted each day for the calculation of
#         the ST)
#     currentDay : the current DOY (** UNUSED ARGUMENT!)

#     Returns
#     -------
#     - Sum of temperatures
#     """
#     # TO-DO: return DOY in doc, but return sumT in code O_O?????
#     if temperature >= t0:
#         sumT += np.max(temperature - tbase, 0)
#     return sumT


def getAvailableBiomassForCut(
    gv_biomass, dv_biomass, gr_biomass, dr_biomass,
    cutHeight, rhogv, rhodv, rhogr, rhodr
):
    """
    Return the amount of biomass av for cut

    Parameters
    ----------
    gv_biomass : Biomass of green vegetation
    dv_biomass : Biomass of dry vegetation
    gr_biomass : Biomass of green reproduction
    dr_biomass : Biomass of dry reproduction
    cutHeight : Height of the cut
    rhogv : Volume VV [g m⁻³]
    rhodv : Volume DV [g m⁻³]
    rhogr : Volume GR [g m⁻³]
    rhodr : Volume DR [g m⁻³]

    Returns
    -------
    - Amount of biomass av for cut
    """
    avDefoliationBiomassGV = avDefoliationBiomass(gv_biomass, cutHeight, rhogv)
    avDefoliationBiomassDV = avDefoliationBiomass(dv_biomass, cutHeight, rhodv)
    avDefoliationBiomassGR = avDefoliationBiomass(gr_biomass, cutHeight, rhogr)
    avDefoliationBiomassDR = avDefoliationBiomass(dr_biomass, cutHeight, rhodr)
    return (
        avDefoliationBiomassGV
        + avDefoliationBiomassDV
        + avDefoliationBiomassGR
        + avDefoliationBiomassDR
    )


def defoliation(
    gv_biomass, dv_biomass, gr_biomass, dr_biomass, cutHeight,
    rhogv, rhodv, rhogr, rhodr, maxAmountToIngest=9999
):
    """
    Defoliation method

    Parameters
    ----------
    gv_biomass : biomass of green vegetation
    dv_biomass : biomass of dry vegetation
    gr_biomass : biomass of green reproduction
    dr_biomass : biomass of dry reproduction
    cutHeight : height of the cut
    rhogv : Volume VV [g m⁻³]
    rhodv : Volume DV [g m⁻³]
    rhogr : Volume GR [g m⁻³]
    rhodr : Volume DR [g m⁻³]
    maxAmountToIngest : The maximum amount of biomass to ingest
        (** USING A DEFAULT OF 9999 - NEED TO CHECK!)

    Returns
    -------
    - the sum of ingested biomass
    """
    avDefoliationBiomassGV = avDefoliationBiomass(gv_biomass, cutHeight, rhogv)
    avDefoliationBiomassDV = avDefoliationBiomass(dv_biomass, cutHeight, rhodv)
    avDefoliationBiomassGR = avDefoliationBiomass(gr_biomass, cutHeight, rhogr)
    avDefoliationBiomassDR = avDefoliationBiomass(dr_biomass, cutHeight, rhodr)
    sumAvailable = (
        avDefoliationBiomassGV
        + avDefoliationBiomassDV
        + avDefoliationBiomassGR
        + avDefoliationBiomassDR
    )
    sumBiomassIngested = 0
    if sumAvailable > 0:
        if sumAvailable <= maxAmountToIngest:
            sumBiomassIngested = (
                exeDefoliationByBiomass(
                    avDefoliationBiomassGV, maxAmountToIngest
                ) + exeDefoliationByBiomass(
                    avDefoliationBiomassDV, maxAmountToIngest
                ) + exeDefoliationByBiomass(
                    avDefoliationBiomassGR, maxAmountToIngest
                ) + exeDefoliationByBiomass(
                    avDefoliationBiomassDR, maxAmountToIngest
                )
            )
        else:
            sumBiomassIngested = (
                exeDefoliationByBiomass(
                    maxAmountToIngest * avDefoliationBiomassGV / sumAvailable,
                    maxAmountToIngest
                ) + exeDefoliationByBiomass(
                    maxAmountToIngest * avDefoliationBiomassDV / sumAvailable,
                    maxAmountToIngest
                ) + exeDefoliationByBiomass(
                    maxAmountToIngest * avDefoliationBiomassGR / sumAvailable,
                    maxAmountToIngest
                ) + exeDefoliationByBiomass(
                    maxAmountToIngest * avDefoliationBiomassDR / sumAvailable,
                    maxAmountToIngest
                )
            )
    return sumBiomassIngested


# def greenLimbsMass(gv_gamma, gv_biomass):
#     """
#     Calculation of mass of green limbs for the cell
#     (old name: mlv)
#     ** THIS FUNCTION IS UNUSED!

#     Parameters
#     ----------
#     gv_biomass : biomass of green vegetation
#     gv_gamma : ?

#     Returns
#     -------
#     - mass of green limbs
#     """
#     return gv_gamma * gv_biomass


# def getOMDgv(gv_min_omd, gv_max_omd, gv_avg_age, lls):
#     """
#     Compute the green vegetation actual organic matter digestibility
#     ** THIS FUNCTION IS UNUSED!

#     Parameters
#     ----------
#     gv_min_omd : The minimum green vegetation organic matter digestibility
#     gv_max_omd : The maximum green vegetation organic matter digestibility
#     gv_avg_age : The average age of the green vegetation
#     lls : the leaf lifespan [°C d]

#     Returns
#     -------
#     - the green vegetation organic matter digestibility
#     """
#     return max(
#         gv_min_omd, gv_max_omd - gv_avg_age * (gv_max_omd - gv_min_omd) / lls
#     )


# def getOMDgr(gr_min_omd, gr_max_omd, gr_avg_age, st1, st2):
#     """
#     Compute the green reproduction actual organic matter digestibility
#     ** THIS FUNCTION IS UNUSED!

#     Parameters
#     ----------
#     gr_min_omd : The minimum green reproduction organic matter digestibility
#     gr_max_omd : The maximum green reproduction organic matter digestibility
#     gr_avg_age : The average age of the green reproduction
#     st1 : sum of temperature to begin vegetative activity
#     st2 : sum of temperature to end vegetative activity

#     Returns
#     -------
#     - the green vegetation organic matter digestibility
#     """
#     return max(
#         gr_min_omd,
#         gr_max_omd - gr_avg_age * (gr_max_omd - gr_min_omd) / (st2 - st1)
#     )


def getSumTemperature(weather, doy, t0):
    """
    Return the sum temperature corresponding to the DOY

    Parameters
    ----------
    weather : Weather array
    doy : Day of the year [1-366]
    t0 : Minimum temperature for growth [°C]

    Returns
    -------
    - Sum of temperatures above t0 corresponding to the DOY
    """
    sumTemperature = 0
    for i in range(doy):
        if weather[i][1] > t0:
            sumTemperature += weather[i][1] - t0
    return sumTemperature


# TO-DO: This set of functions are either not used or not useful
# def addNI(ni, amountToIncrease):
#     return max(0, min(amountToIncrease + ni, 1.2))
