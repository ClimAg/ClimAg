"""lib_modvege.py

https://github.com/YannChemin/modvege
"""

import numpy as np


def getAverageHeight(biomass, bulkDensity):
    """
    Return the height based on biomass and bulk density

    Parameters
    ----------
    biomass : biomass
    bulkDensity : bulk density

    Returns
    -------
    - the average height
    """
    return biomass/(bulkDensity)


def avDefoliationBiomass(biomass, cutHeight, bulkDensity):
    """
    Estimate av biomass for ingestion

    Parameters
    ----------
    biomass : av biomass
    cutHeight : default cut height in cutEvent
    bulkDensity : bulk density

    Returns
    -------
    - av biomass for ingestion
    """
    biomassAfterCut = cutHeight * bulkDensity * 10
    return max(0, biomass - biomassAfterCut)


def exeDefoliation(biomass, cut_biomass, area):
    """
    Defoliation method

    Parameters
    ----------
    biomass : av biomass
    cut_biomass : biomass removed by cut
    area : area studied (i.e. pixel area)

    Returns
    -------
    - updated biomass after defoliation
    """
    biomass = biomass - cut_biomass/area
    if biomass < 0 | np.isnan(biomass):
        biomass = 0
    return biomass


def exeCut(cutHeight, bulkDensity, biomass):
    """
    Realize a cut in order that the average height is under cutHeight.
    This height is calculated by using the bulkDensity given in parameter

    Parameters
    ----------
    cutHeight : the average height after the cut in m
    bulkDensity : the bulk density
        (biomass after cut = height * 10 * bulk density)
    biomass : biomass

    Returns
    -------
    - the biomass taken in kg DM / m²
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
    Realize a defoliation by providing the biomass to remove

    Parameters
    ----------
    biomass : av biomass
    biomassToIngest : the biomass to remove

    Returns
    -------
    - the biomass left in kg DM / m²
    """
    biomass -= biomassToIngest
    return biomassToIngest


# Dry Vegetative Functions
def dv_update(
    gv_gamma,
    gv_senescent_biomass,
    lls,
    kldv,
    temperature,
    dv_biomass,
    dv_avg_age
):
    """
    Update DV

    Parameters
    ----------
    gv_gamma : Respiratory C loss during senescence (DV)
        (1 - gv_gamma = dv_gamma)
    gv_senescent_biomass : senescence of compartment GV
    lls : Leaf lifespan (degree day)
    kldv : Abscission coefficient DV (degree day)
    temperature : temperature
    dv_biomass : av DV biomass
    dv_avg_age : the average DV age

    Returns
    -------
    - the Dry Vegetation biomass
    - the average DV age
    """
    abscissionBiomass = mk_dv_abscission(
        kldv, dv_biomass, temperature, dv_avg_age, lls
    )
    dv_biomass -= abscissionBiomass
    # at this point the biomass include cut, ingestion, and abscission, not
    # growth
    growthBiomass = (1.0 - gv_gamma) * gv_senescent_biomass
    if dv_biomass + growthBiomass > 0:
        dv_avg_age = (
            (max(0, temperature) + dv_avg_age)
            * (dv_biomass/(dv_biomass + growthBiomass))
        )
    else:
        dv_avg_age = 0
    dv_biomass += growthBiomass
    return (dv_biomass, dv_avg_age)


def mk_dv_abscission(kldv, dv_biomass, temperature, dv_avg_age, lls):
    """
    Compute abscission biomass

    Parameters
    ----------
    lls : Leaf lifespan (degreeday)
    kldv : Abscission coefficient DV (degree day) (Ducroq, 1996)
    temperature : temperature
    dv_biomass : av biomass
    dv_avg_age : the average DV age

    Returns
    -------
    - the abscission biomass
    """
    # method to compute the age of DV for computing abcission of DV
    if dv_avg_age/lls < 1.0/3.0:
        age = 1
    elif dv_avg_age/lls < 2.0/3.0:
        age = 2
    else:
        age = 3
    # Compute the abscission for Dead Vegetative part
    if temperature > 0:
        abscission_biomass = kldv * dv_biomass * temperature * age
    else:
        abscission_biomass = 0
    return abscission_biomass


# Dead Reproductive Functions
def dr_update(
    gr_gamma, gr_senescent_biomass, st1, st2,
    temperature, kldr, dr_biomass, dr_avg_age
):
    """
    Update Dead Reproductive compartment

    Parameters
    ----------
    gr_gamma : a parameter to compute growth of DR (1-gr_gamma)=dr_gamma
    gr_senescent_biomass : Senescence of GR, computed in GR
    st1 : sum of temperature at the beginning
    st2 : sum of temperature in the end
    temperature : temperature
    kldr : basic rates of abscission in DR
    dr_biomass : av DR biomass
    dr_avg_age : the average DR age

    Returns
    -------
    - updated biomass for DR
    - the average DR age
    """
    abscissionBiomass = mk_dr_abscission(
        kldr, dr_biomass, temperature, dr_avg_age, st1, st2
    )
    dr_biomass -= abscissionBiomass
    # at this point the biomass include cut, ingestion and abscission, not
    # growth
    growthBiomass = (1 - gr_gamma) * gr_senescent_biomass
    if dr_biomass+growthBiomass > 0:
        # print(temperature, dr_avg_age, dr_biomass, growthBiomass)
        dr_avg_age = (
            max(0, temperature) +
            dr_avg_age * dr_biomass/(dr_biomass + growthBiomass)
        )
    else:
        dr_avg_age = 0
    dr_biomass += growthBiomass
    return (dr_biomass, dr_avg_age)


def mk_dr_abscission(kldr, dr_biomass, temperature, dr_avg_age, st1, st2):
    """
    Compute abscission biomass

    Parameters
    ----------
    kldr : basic rates of abscission in DR (Ducroq, 1996)
    dr_biomass : av biomass
    temperature : temperature
    dr_avg_age : the average DR age
    st1 : sum of temperature at the beginning
    st2 : sum of temperature in the end

    Returns
    -------
    - the abscission biomass
    """
    # method to compute the age of DR for computing abscission of DR
    if dr_avg_age/(st2 - st1) < 1.0/3.0:
        age = 1
    elif dr_avg_age/(st2 - st1) < 2.0/3.0:
        age = 2
    else:
        age = 3
    # Compute abscission for Dead Reproductive
    if temperature > 0:
        # print(kldr, dr_biomass, temperature,age)
        abscission_biomass = kldr * dr_biomass * temperature * age
    else:
        abscission_biomass = 0
    return abscission_biomass


# Green Vegetative Functions
def gv_update(gro, a2r, lls, temperature, kdv, t0, gv_biomass, gv_avg_age):
    """
    Update Green Vegetation

    Parameters
    ----------
    gro : in Jouven et al. (2006), total growth
    a2r : Allocate to reproductive
        (REP in Jouven et al. (2006), reproductive function)
    lls : Leaf lifespan (degree day)
    temperature : temperature
    kdv : Senescence coefficient DV (degreeday)
    t0 : minimum temperature for growth
    gv_biomass : Updated biomass
    gv_avg_age : the average GV age

    Returns
    -------
    - senescent biomass
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
        gv_avg_age = (
            (max(0, temperature) + gv_avg_age) *
            (gv_biomass/(gv_biomass+growthBiomass))
        )
    else:
        gv_avg_age = 0
    gv_biomass += growthBiomass
    return (gv_biomass, gv_avg_age, senescentBiomass)


def mk_gv_senescence(kgv, gv_biomass, temperature, t0, lls, gv_avg_age):
    """
    Extract about 2-6% (kDV=0.002, T=10C, gv_fAge=[1-3]) of gv_biomass as
    senescent

    Parameters
    ----------
    kgv : Senescence coefficient (degree day) (Ducroq, 1996)
    gv_biomass : the Green Vegetation biomass
    temperature : Temperature
    t0 : minimum temperature for growth
    lls : Leaf lifespan (degree day)
    gv_avg_age : the average GV age

    Returns
    -------
    - the biomass that is senescent
    """
    # method to compute the age of GV for computing senescence of GV
    if gv_avg_age/lls < 1.0/3.0:
        age = 1
    elif gv_avg_age/lls < 1:
        age = 3 * gv_avg_age / lls
    else:
        age = 3
    # Compute senescence of GV
    if temperature > t0:
        senescence_biomass = kgv * gv_biomass * temperature * age
    elif temperature < 0:
        senescence_biomass = kgv * gv_biomass * abs(temperature)
    else:
        senescence_biomass = 0
    return senescence_biomass


# Green Reproductive Functions
def gr_update(
    temperature, a2r, gro, st1, st2, kdr, lls,
    rhogr, t0, gr_biomass, gr_avg_age
):
    """
    Update Green Reproductive

    Parameters
    ----------
    temperature : temperature
    a2r : allocate to reproductive
        (REP in Jouven et al. (2006), reproductive function)
    gro : in Jouven et al. (2006), total growth
    st1 : Onset of reproductive growth (degree day)
    st2 : End of reproductive growth (degree day)
    kdr : basic rates of  in compartment GR
    lls : Leaf lifespan (degreeday)
    rhogr : Volume GR (g m-3)
    t0 : minimum temperature for growth
    gr_biomass : the av GR biomass
    gr_avg_age : the average GR age

    Returns
    -------
    - Updated GR biomass
    - the average GR age
    """
    senescentBiomass = mk_gr_senescence(
        kdr, gr_biomass, temperature, t0, lls, gr_avg_age, st1, st2
    )
    gr_biomass -= senescentBiomass
    # print("senescentBiomass: = %.2f" % (senescentBiomass))
    # at this point the biomass include cut, ingestion and senescence, not
    # growth
    if temperature > t0:
        growthBiomass = gro*(a2r)
        # print("growthBiomass: t>t0 = %.2f" % (growthBiomass))
    else:
        growthBiomass = 0
        # print("growthBiomass: t<t0 = %.2f" % (growthBiomass))
    if gr_biomass+growthBiomass > 0:
        gr_avg_age = (
            (max(0, temperature) + gr_avg_age) *
            (gr_biomass/(gr_biomass + growthBiomass))
        )
    else:
        gr_avg_age = 0
    gr_biomass += growthBiomass
    return (gr_biomass, gr_avg_age, senescentBiomass)


def mk_gr_senescence(
    kdr, gr_biomass, temperature, t0, lls, gr_avg_age, st1, st2
):
    """
    Parameters
    ----------
    kdr : Senescence coefficient DV (degree day) **CHECK PARAM NAME!
    gr_biomass : the biomass available for GR
    temperature : Temperature
    t0 : minimum temperature for growth
    lls : Leaf lifespan (degreeday)
    gr_avg_age : the average GR age
    st1 : Onset of reproductive growth (degreeday)
    st2 : End of reproductive growth (degreeday)

    Returns
    -------
    - the senescent biomass
    """
    # method to compute the age of GR for computing senescence of GR
    if gr_avg_age/(st2 - st1) < 1.0/3.0:
        age = 1
    elif gr_avg_age/(st2 - st1) < 1.0:
        age = 3 * gr_avg_age/(st2 - st1)
    else:
        age = 3
    # Compute senescence of GV
    if temperature > t0:
        # T=10C kdr = 0.001 gr_fAge=[1-3] => 1-3% of gr_biomass
        senescence_biomass = kdr * gr_biomass * temperature * age
    elif temperature < 0:
        senescence_biomass = kdr * gr_biomass * abs(temperature)
    else:
        senescence_biomass = 0
    return senescence_biomass

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


def getHeight(gv_avg_h, gr_avg_h, dv_avg_h, dr_avg_h):
    """
    Return the height of this cell

    Parameters
    ----------
    gv_avg_h : ?
    gr_avg_h : ?
    dv_avg_h : ?
    dr_avg_h : ?

    Returns
    -------
    - the maximum height of the 4 cs
    """
    return max(max(gv_avg_h, gr_avg_h), np.max(dv_avg_h, dr_avg_h))


def cut(
    cutHeight, rhogv, rhodv, rhogr, rhodr, gvb, dvb,
    grb, drb, cellSurface, isHarvested
):
    """
    Realize the harvest on each c. If the amount of cut biomass is not null,
    then the flag isHarvested is set to True

    Parameters
    ----------
    cutHeight : the height of the cut (m)
    rhogv : rho green vegetation
    rhodv : rho dry vegetation
    rhogr : rho green reproduction
    rhodr : rho dry reproduction
    gvb : the biomass of Green Vegetation
    dvb : the biomass of Dry Vegetation
    grb : the biomass of Green Reproductive
    drb : the biomass of Dry Reproductive
    cellSurface : Surface of the pixel (Ha)
    isHarvested : Status flag indicating harvest happened

    Returns
    -------
    - Status flag indicating harvest happened
    - the total amount of biomass cut (in kg DM)
    - the amount of GV biomass cut (in kg DM)
    - the amount of DV biomass cut (in kg DM)
    - the amount of GR biomass cut (in kg DM)
    - the amount of DR biomass cut (in kg DM)
    """
    # exeCut returns harvested biomass in [kg DM m-2]
    gv_h, gv_b = exeCut(rhogv, cutHeight, gvb)
    dv_h, dv_b = exeCut(rhodv, cutHeight, dvb)
    gr_h, gr_b = exeCut(rhogr, cutHeight, grb)
    dr_h, dr_b = exeCut(rhodr, cutHeight, drb)
    # sum of harvested biomass [kg DM m-2]
    sumBiomassHarvested = gv_h + dv_h + gr_h + dr_h
    if sumBiomassHarvested > 0:
        isHarvested = True
    return (
        isHarvested, sumBiomassHarvested * cellSurface, gv_b, dv_b, gr_b, dr_b
    )


def mk_env(
    meanTenDaysT, t0, t1, t2, sumT, ni, pari, alphapar,
    pet, waterReserve, waterHoldingCapacity
):
    """
    Environmental stress

    Parameters
    ----------
    meanTenDaysT : the mean of the ten days of temperature
    t0 : minimum temperature for growth
    t1 : sum of temperature at the beginning (growth activation threshold)
    t2 : sum of temperature in the end (growth decline threshold)
    sumT : sum of temperatures
    ni : Nutritional index of pixel
    pari : incident photosynthetic active radiation (PARi)
    alphapar : the Light Use Interception
    pet : potential evapotranspiration
    waterReserve : reserve of water in the soil
    waterHoldingCapacity : capacity of the soil to hold a certain volume
        of water

    Returns
    -------
    - the environmental stress
    """
    return (
        fTemperature(meanTenDaysT, t0, t1, t2, sumT) *
        ni * fPARi(pari, alphapar) *
        fWaterStress(waterReserve, waterHoldingCapacity, pet)
    )


def getTotalBiomass(gv_biomass, dv_biomass, gr_biomass, dr_biomass):
    """
    Return the total biomass of the cell (by adding the biomass of the 4 cs)

    Parameters
    ----------
    gv_biomass : the Green Vegetation biomass
    dv_biomass : the Dry Vegetation biomass
    gr_biomass : the Green Reproduction biomass
    dr_biomass : the Dry Reproduction biomass

    Returns
    -------
    - total biomass
    """
    return gv_biomass + dv_biomass + gr_biomass + dr_biomass


def fTemperature(meanTenDaysT, t0, t1, t2, sumT):
    """
    f of temperature to compute ENV

    Parameters
    ----------
    meanTenDaysT : the mean of the ten days of temperature
    t0 : minimum temperature for growth
    t1 : sum of temperature at the beginning (growth activation threshold)
    t2 : sum of temperature in the end (growth decline threshold)
    sumT : sum of temperatures

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
    Function for seasonality (SEA) to compute the Potential Growth

    Parameters
    ----------
    maxsea : growth increase in summer
    minsea : growth increase in winter
    sumT : sum of temperature
    st1 : sum of temperature at the beginning of growth
    st2 : sum of temperature in the end of growth

    Returns
    -------
    - the value given by the sea f
    """
    if sumT < 200 or sumT >= st2:
        f_sea = minsea
    elif sumT < st1 - 200:
        f_sea = minsea + (maxsea - minsea) * (sumT - 200)/(st1 - 400)
    elif sumT < st1 - 100:
        f_sea = maxsea
    else:
        f_sea = (
            maxsea + (minsea - maxsea) * (sumT - st1 + 100)/(st2 - st1 + 100)
        )
    return f_sea


def fPARi(pari, alphapar):
    """
    Function of PAR interception (PARi) to compute ENV

    Parameters
    ----------
    pari : Photosynthetic radiation incident (PARi)
    alphapar : the Light Use Interception

    Returns
    -------
    - the value given by the PARi [0-1]
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
    waterStress = min(waterReserve/waterHoldingCapacity, 1)
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
    Replace the value of NI (Nutrition Index (Belanger et al., 1992))

    Parameters
    ----------
    ni : The Nutrition Index (Belanger et al., 1992)

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
    pari : the incident PAR
    ruemax : the maximum Radiation Use Efficiency
    pctlam : % leaf of laminae in Green Vegetation
    sla : the specific leaf area (m2 g-1)
    gv_biomass : the Green Vegetation biomass
    lai : the LAI from remote sensing (if available)

    Returns
    -------
    - the calculated pGRO (kg DM ha-1)
    """
    # print("pgro**************************")
    # print("pgro: pariIn = %.2f" % (pari))
    # print("pgro: ruemaxIn = %.2f" % (ruemax))
    # print("pgro: pctlamIn = %.2f" % (pctlam))
    # print("pgro: SLAIn = %.2f" % (sla))
    # print("pgro: gv_biomassIn = %.2f" % (gv_biomass))
    if int(lai) == 0:
        try:
            lai = sla * pctlam * (gv_biomass/10)
        except:
            # In case of input malfunction
            lai = sla * pctlam * 1.0
    lightInterceptionByPlant = (1 - np.exp(-0.6 * lai))
    # print("pgro: LightInter.byPlant = %.2f" % (lightInterceptionByPlant))
    p_gro = (pari * ruemax * lightInterceptionByPlant * 10)
    # print("pgro: pgro = %.2f" % (pgro))
    return p_gro


def fclai(pctlam, sla, gv_biomass):
    """
    Compute and return the leaf area index

    Parameters
    ----------
    pctlam : % leaf of laminae in Green Vegetation
    sla : the specific leaf area (m2 g-1)
    gv_biomass : the Green Vegetation biomass

    Returns
    -------
    - the calculated LAI
    """
    return sla * (gv_biomass/10) * pctlam


def aet(pet, pctlam, sla, gv_biomass, waterReserve, waterHoldingCapacity, lai):
    """
    Return the actual evapotranspiration (AET)

    pet : the daily potential evapotranspiration (PET)
    pctlam : % leaf of laminae in Green Vegetation
    sla : the specific leaf area (m2 g-1)
    gv_biomass : the Green Vegetation biomass
    waterReserve : reserve of water in the soil
    waterHoldingCapacity : capacity of the soil to hold a certain volume
        of water
    lai : ?

    Returns
    -------
    - the actual evapotranspiration (AET)
    """
    if int(lai) == 0:
        lai = sla * pctlam * (gv_biomass/10)
        # print("aet mk LAI: LAI = %.2f" %(lai))
    lightInterceptionByPlant = (1 - np.exp(-0.6 * lai))
    pt = pet * lightInterceptionByPlant
    pe = pet - pt
    ta = pt * fWaterStress(waterReserve, waterHoldingCapacity, pet)
    ea = pe * min(waterReserve/waterHoldingCapacity, 1)
    return ta + ea


def updateSumTemperature(temperature, t0, sumT, tbase):
    """
    Add the daily temperature to the sum temperature if the daily one is
        positive

    Parameters
    ----------
    temperature : ?
    t0 : minimum temperature for growth
    sumT : actual sum of temperature
    tbase : base temperature (subtracted each day for the calculation of
        the ST)
    currentDay : the current DOY **CHECK WHY THIS ISN'T INCLUDED!

    Returns
    -------
    - sum T
    """
    # TO-DO return DOY in doc, but return sumT in code O_O ?????
    if temperature >= t0:
        sumT += np.max(temperature - tbase, 0)
    return sumT


def getAvailableBiomassForCut(
    gv_biomass, dv_biomass, gr_biomass, dr_biomass,
    cutHeight, rhogv, rhodv, rhogr, rhodr
):
    """
    Return the amount of biomass av for cut

    Parameters
    ----------
    gv_biomass : biomass of Green Vegetation
    dv_biomass : biomass of Dry Vegetation
    gr_biomass : biomass of Green Reproduction
    dr_biomass : biomass of Dry Reproduction
    cutHeight : height of the cut
    rhogv : Volume VV (g m-3)
    rhodv : Volume DV (g m-3)
    rhogr : Volume GR (g m-3)
    rhodr : Volume DR (g m-3)

    Returns
    -------
    - the amount of biomass av for cut
    """
    avDefoliationBiomassGV = avDefoliationBiomass(gv_biomass, cutHeight, rhogv)
    avDefoliationBiomassDV = avDefoliationBiomass(dv_biomass, cutHeight, rhodv)
    avDefoliationBiomassGR = avDefoliationBiomass(gr_biomass, cutHeight, rhogr)
    avDefoliationBiomassDR = avDefoliationBiomass(dr_biomass, cutHeight, rhodr)
    return (
        avDefoliationBiomassGV + avDefoliationBiomassDV +
        avDefoliationBiomassGR + avDefoliationBiomassDR
    )


def defoliation(
    gv_biomass, dv_biomass, gr_biomass, dr_biomass, cutHeight,
    rhogv, rhodv, rhogr, rhodr, maxAmountToIngest
):
    """
    Defoliation method

    Parameters
    ----------
    gv_biomass : biomass of Green Vegetation
    dv_biomass : biomass of Dry Vegetation
    gr_biomass : biomass of Green Reproduction
    dr_biomass : biomass of Dry Reproduction
    cutHeight : height of the cut
    rhogv : Volume VV (g m-3)
    rhodv : Volume DV (g m-3)
    rhogr : Volume GR (g m-3)
    rhodr : Volume DR (g m-3)
    maxAmountToIngest : The maximum amount of ingest

    Returns
    -------
    - the sum of ingested biomass
    - isGrazed = True
    """
    avDefoliationBiomassGV = avDefoliationBiomass(gv_biomass, cutHeight, rhogv)
    avDefoliationBiomassDV = avDefoliationBiomass(dv_biomass, cutHeight, rhodv)
    avDefoliationBiomassGR = avDefoliationBiomass(gr_biomass, cutHeight, rhogr)
    avDefoliationBiomassDR = avDefoliationBiomass(dr_biomass, cutHeight, rhodr)
    sumAvailable = (
        avDefoliationBiomassGV + avDefoliationBiomassDV +
        avDefoliationBiomassGR + avDefoliationBiomassDR
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
                    maxAmountToIngest * avDefoliationBiomassGV/sumAvailable,
                    maxAmountToIngest
                ) + exeDefoliationByBiomass(
                    maxAmountToIngest * avDefoliationBiomassDV/sumAvailable,
                    maxAmountToIngest
                ) + exeDefoliationByBiomass(
                    maxAmountToIngest * avDefoliationBiomassGR/sumAvailable,
                    maxAmountToIngest
                ) + exeDefoliationByBiomass(
                    maxAmountToIngest * avDefoliationBiomassDR/sumAvailable,
                    maxAmountToIngest
                )
            )
    return sumBiomassIngested


def greenLimbsMass(gv_gamma, gv_biomass):
    """
    Calculation of mass of green limbs for the cell
    (old name: mlv)

    Parameters
    ----------
    gv_biomass : biomass of Green Vegetation
    gv_gamma : ?

    Returns
    -------
    - mass of green limbs
    """
    return gv_gamma * gv_biomass


def getOMDgv(gv_min_omd, gv_max_omd, gv_avg_age, lls):
    """
    Compute the Green Vegetation actual Organic Matter Digestibility

    Parameters
    ----------
    gv_min_omd : The minimum Green Vegetation Org. Mat. Digestibility
    gv_max_omd : The maximum Green Vegetation Org. Mat. Digestibility
    gv_avg_age : The average age of the Green Vegetation
    lls : the leaf lifespan (degree Celsius per day-1)

    Returns
    -------
    - the Green Vegetation Organic Matter Digestibility
    """
    return max(
        gv_min_omd, gv_max_omd - gv_avg_age *
        (gv_max_omd - gv_min_omd)/lls
    )


def getOMDgr(gr_min_omd, gr_max_omd, gr_avg_age, st1, st2):
    """
    Compute the Green Reproduction actual Organic Matter Digestibility

    Parameters
    ----------
    @param gr_min_omd The minimum Green Reproduction Org. Mat. Digestibility
    @param gr_max_omd The maximum Green Reproduction Org. Mat. Digestibility
    @param gr_avg_age The average age of the Green Reproduction
    @param st1 sum of temperature to begin vegeative activity
    @param st2 sum of temperature to end vegetative activity

    Returns
    -------
    - the Green Vegetation Organic Matter Digestibility
    """
    return max(
        gr_min_omd,
        gr_max_omd - gr_avg_age * (gr_max_omd - gr_min_omd)/(st2 - st1)
    )


def getSumTemperature(weather, doy, t0):
    """
    Return the sum temperature corresponding to the DOY

    Parameters
    ----------
    weather : the weather array
    doy : the day of year wanted [1-366]
    t0 : minimum temperature for growth

    Returns
    -------
    - the sum temperature above t0 corresponding to the DOY
    """
    sumTemperature = 0
    for i in range(doy):
        if weather[i][1] > t0:
            sumTemperature += (weather[i][1] - t0)
    return sumTemperature


# TO-DO This set of functions are either not used or not useful
def addNI(ni, amountToIncrease):
    return max(0, min(amountToIncrease + ni, 1.2))
