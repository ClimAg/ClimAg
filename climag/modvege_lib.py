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
    Estimate biomass available for ingestion by livestock.

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
    See Jouven et al. (2006), sec. "Harvested biomass", equation (19).
    Assumption: during harvest, 10% of the harvestable biomass in each
    structural component is lost.

    Parameters
    ----------
    cutHeight : Average height after the cut [m]
    bulkDensity : Bulk density [g DM m⁻³]
    biomass : Biomass [kg DM ha⁻¹]

    Returns
    -------
    - Harvested biomass [kg DM ha⁻¹]
    - Residual biomass [kg DM ha⁻¹]
    """

    residual_biomass = cutHeight * bulkDensity * 10
    if residual_biomass < biomass:
        harvested_biomass = (biomass - residual_biomass) * .9
    else:
        harvested_biomass = 0
        residual_biomass = biomass
    return (harvested_biomass, residual_biomass)


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
    return biomass


# DEAD VEGETATIVE FUNCTIONS
def mk_dv_abscission(kldv, dv_biomass, temperature, dv_avg_age, lls):
    """
    Compute abscission biomass for the DV compartment.
    See Equation (18) in Jouven et al. (2006).

    Parameters
    ----------
    lls : Leaf lifespan (LLS) [500 °C d]
    kldv : Basic abscission rate for the dead vegetative compartment
        (Kl_DV) [0.001]
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
    gv_gamma, gv_senescent_biomass, lls, kldv, temperature, dv_biomass,
    dv_avg_age
):
    """
    Update dead vegetative compartment.
    See Equation (3) in Jouven et al. (2006).

    Parameters
    ----------
    gv_gamma : Respiratory C loss during senescence (DV)
        (1 - gv_gamma = dv_gamma)
    gv_senescent_biomass : Senescence of GV compartment
    lls : Leaf lifespan (LLS) [500 °C d]
    kldv : Basic abscission rate for the dead vegetative compartment
        (Kl_DV) [0.001]
    temperature : Mean daily temperature (*T*) [°C]
    dv_biomass : DV biomass (BM_DV) [kg DM ha⁻¹]
    dv_avg_age : Average DV age (AGE_DV) [°C d]

    Returns
    -------
    - Dead vegetative biomass
    - Average DV age
    """

    abscissionBiomass = mk_dv_abscission(
        kldv=kldv, dv_biomass=dv_biomass, temperature=temperature,
        dv_avg_age=dv_avg_age, lls=lls
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
    kldr : Basic abscission rate for the dead reproductive compartment
        (Kl_DR) [0.0005]
    dr_biomass : DR biomass (BM_DR) [kg DM ha⁻¹]
    temperature : Mean daily temperature (*T*) [°C]
    dr_avg_age : Average DR age (AGE_DR) [°C d]
    st1 : Sum of temperatures at the beginning of the reproductive period
        (ST₁) [600 °C d]
    st2 : Sum of temperatures at the end of the reproductive period
        (ST₂) [1200 °C d]

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
    gr_gamma, gr_senescent_biomass, st1, st2, temperature, kldr, dr_biomass,
    dr_avg_age
):
    """
    Update dead reproductive compartment.

    Parameters
    ----------
    gr_gamma : Parameter to compute growth of DR (1 - gr_gamma = dr_gamma)
    gr_senescent_biomass : Senescence of GR, computed in GR
    st1 : Sum of temperatures at the beginning of the reproductive period
        (ST₁) [600 °C d]
    st2 : Sum of temperatures at the end of the reproductive period
        (ST₂) [1200 °C d]
    temperature : Mean daily temperature (*T*) [°C]
    kldr : Basic abscission rate for the dead reproductive compartment
        (Kl_DR) [0.0005]
    dr_biomass : DR biomass
    dr_avg_age : Average DR age [°C d]

    Returns
    -------
    - Updated biomass for DR
    - Average DR age
    """

    abscissionBiomass = mk_dr_abscission(
        kldr=kldr, dr_biomass=dr_biomass, temperature=temperature,
        dr_avg_age=dr_avg_age, st1=st1, st2=st2
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
    kgv : Basic senescence rate for the green vegetative compartment
        (K_GV) [0.002]
    gv_biomass : GV biomass (BM_GV) [kg DM ha⁻¹]
    temperature : Mean daily temperature (*T*) [°C]
    t0 : Minimum temperature for growth (*T*₀) [°C]
    lls : Leaf lifespan (LLS) [500 °C d]
    gv_avg_age : Average GV age (AGE_GV) [°C d]

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


def gv_update(gro, a2r, lls, temperature, kgv, t0, gv_biomass, gv_avg_age):
    """
    Update green vegetative compartment.

    Parameters
    ----------
    gro : in Jouven et al. (2006), total growth (GRO)
    a2r : Allocate to reproductive
        (REP in Jouven et al. (2006), reproductive function)
    lls : Leaf lifespan (LLS) [500 °C d]
    temperature : Mean daily temperature (*T*) [°C]
    kgv : Senescence coefficient GV [°C d]
    t0 : Minimum temperature for growth (*T*₀) [4 °C]
    gv_biomass : Updated biomass
    gv_avg_age : Average GV age [°C day]

    Returns
    -------
    - GV biomass
    - GV average age
    - Senescent biomass
    """

    senescentBiomass = mk_gv_senescence(
        kgv=kgv, gv_biomass=gv_biomass, temperature=temperature, t0=t0,
        lls=lls, gv_avg_age=gv_avg_age
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
def mk_gr_senescence(kgr, gr_biomass, temperature, t0, gr_avg_age, st1, st2):
    """
    See Equations (16) and (17) in Jouven et al. (2006).

    Parameters
    ----------
    kgr : Basic senescence rate for the green reproductive compartment
        (K_GR) [0.001]
    gr_biomass : Biomass available for GR (BM_GR) [kg DM ha⁻¹]
    temperature : Mean daily temperature (*T*) [°C]
    t0 : Minimum temperature for growth (*T*₀) [4 °C]
    gr_avg_age : Average GR age (AGE_GR) [°C d]
    st1 : Sum of temperatures at the beginning of the reproductive period
        (ST₁) [600 °C d]
    st2 : Sum of temperatures at the end of the reproductive period
        (ST₂) [1200 °C d]
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
    Update green reproductive compartment.

    Parameters
    ----------
    temperature : Mean daily temperature (*T*) [°C]
    a2r : Allocate to reproductive
        (REP in Jouven et al. (2006), reproductive function)
    gro : Total growth (GRO) [kg DM ha⁻¹]
    st1 : Onset of reproductive growth (ST₁) [°C d]
    st2 : End of reproductive growth (ST₂) [°C d]
    kgr : Basic rates of senescence in compartment GR (K_GR)
    t0 : Minimum temperature for growth (*T*₀) [4 °C]
    gr_biomass : GR biomass (BM_GR) [kg DM ha⁻¹]
    gr_avg_age : Average GR age (AGE_GR) [°C d]
    lls : Leaf lifespan (LLS) [°C d] (** UNUSED ARGUMENT!)
    rhogr : Bulk density of GR [g m⁻³] (** UNUSED ARGUMENT!)

    Returns
    -------
    - Updated GR biomass
    - Average GR age
    - Senescent biomass
    """

    senescentBiomass = mk_gr_senescence(
        kgr=kgr, gr_biomass=gr_biomass, temperature=temperature, t0=t0,
        gr_avg_age=gr_avg_age, st1=st1, st2=st2
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
# Compartment dead vegetative       # CompartmentDeadVegetative cDV
# Compartment dead reproductive     # CompartmentDeadReproductive cDR
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
    cutHeight, rhogv, rhodv, rhogr, rhodr, gvb, dvb, grb, drb
):
    """
    Realise the harvest on each compartment.

    Parameters
    ----------
    cutHeight : the height of the cut [m]
    rhogv : bulk density of green vegetative [g m⁻³]
    rhodv : bulk density of dead vegetative [g m⁻³]
    rhogr : bulk density of green reproductive [g m⁻³]
    rhodr : bulk density of dead reproductive [g m⁻³]
    gvb : standing biomass of green vegetative [kg DM ha⁻¹]
    dvb : standing biomass of dead vegetative [kg DM ha⁻¹]
    grb : standing biomass of green reproductive [kg DM ha⁻¹]
    drb : standing biomass of dead reproductive [kg DM ha⁻¹]

    Returns
    -------
    - the total amount of biomass cut [kg DM ha⁻¹]
    - the amount of GV biomass cut [kg DM ha⁻¹]
    - the amount of DV biomass cut [kg DM ha⁻¹]
    - the amount of GR biomass cut [kg DM ha⁻¹]
    - the amount of DR biomass cut [kg DM ha⁻¹]
    """

    gv_h, gv_b = exeCut(bulkDensity=rhogv, cutHeight=cutHeight, biomass=gvb)
    dv_h, dv_b = exeCut(bulkDensity=rhodv, cutHeight=cutHeight, biomass=dvb)
    gr_h, gr_b = exeCut(bulkDensity=rhogr, cutHeight=cutHeight, biomass=grb)
    dr_h, dr_b = exeCut(bulkDensity=rhodr, cutHeight=cutHeight, biomass=drb)
    sumBiomassHarvested = gv_h + dv_h + gr_h + dr_h
    # if sumBiomassHarvested > 0:
    #     isHarvested = True
    return (sumBiomassHarvested, gv_b, dv_b, gr_b, dr_b)


def mk_env(
    meanTenDaysT, t0, t1, t2, ni, pari, alphapar, pet, waterReserve,
    waterHoldingCapacity
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
        fTemperature(
            meanTenDaysT=meanTenDaysT, t0=t0, t1=t1, t2=t2
        )  # ** UNUSED ARGUMENT REMOVED!
        * ni
        * fPARi(pari=pari, alphapar=alphapar)
        * fWaterStress(
            waterReserve=waterReserve,
            waterHoldingCapacity=waterHoldingCapacity,
            pet=pet
        )
    )


# def getTotalBiomass(gv_biomass, dv_biomass, gr_biomass, dr_biomass):
#     """
#     Return the total biomass of the cell (by adding the biomass of the 4 cs)
#     ** THIS FUNCTION IS UNUSED!

#     Parameters
#     ----------
#     gv_biomass : Green vegetative biomass
#     dv_biomass : Dead vegetative biomass
#     gr_biomass : Green reproductive biomass
#     dr_biomass : Dead reproductive biomass

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
    Function for seasonality (SEA) to compute the potential growth.

    Parameters
    ----------
    maxsea : Maximum seasonal effect (maxSEA) [1.20]
    minsea : Minimum seasonal effect (minSEA) [0.80]
    sumT : Sum of temperatures (ST) [°C d]
    st1 : Sum of temperatures at the beginning of the reproductive period
        (ST₁) [600 °C d]
    st2 : Sum of temperatures at the end of the reproductive period
        (ST₂) [1200 °C d]

    Returns
    -------
    - Value given by *f*(SEA)
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
    waterReserve : Reserve of water in the soil
    waterHoldingCapacity : Capacity of the soil to hold a certain volume
        of water
    pet : Potential evapotranspiration

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
    gv_biomass : Green vegetative biomass (BM_GV) [kg DM ha⁻¹]
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
    ta = pt * fWaterStress(
        waterReserve=waterReserve,
        waterHoldingCapacity=waterHoldingCapacity,
        pet=pet
    )
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
    gv_biomass, dv_biomass, gr_biomass, dr_biomass, cutHeight, rhogv, rhodv,
    rhogr, rhodr
):
    """
    Return the amount of biomass av for cut

    Parameters
    ----------
    gv_biomass : Biomass of green vegetative
    dv_biomass : Biomass of dead vegetative
    gr_biomass : Biomass of green reproductive
    dr_biomass : Biomass of dead reproductive
    cutHeight : Height of the cut [m]
    rhogv : Volume VV [g m⁻³]
    rhodv : Volume DV [g m⁻³]
    rhogr : Volume GR [g m⁻³]
    rhodr : Volume DR [g m⁻³]

    Returns
    -------
    - Amount of biomass available for cut
    """

    avDefoliationBiomassGV = avDefoliationBiomass(
        biomass=gv_biomass, cutHeight=cutHeight, bulkDensity=rhogv
    )
    avDefoliationBiomassDV = avDefoliationBiomass(
        biomass=dv_biomass, cutHeight=cutHeight, bulkDensity=rhodv
    )
    avDefoliationBiomassGR = avDefoliationBiomass(
        biomass=gr_biomass, cutHeight=cutHeight, bulkDensity=rhogr
    )
    avDefoliationBiomassDR = avDefoliationBiomass(
        biomass=dr_biomass, cutHeight=cutHeight, bulkDensity=rhodr
    )
    return (
        avDefoliationBiomassGV
        + avDefoliationBiomassDV
        + avDefoliationBiomassGR
        + avDefoliationBiomassDR
    )


def defoliation(
    gv_biomass, dv_biomass, gr_biomass, dr_biomass, cutHeight, rhogv, rhodv,
    rhogr, rhodr, maxAmountToIngest=9999999
):
    """
    Defoliation method

    Parameters
    ----------
    gv_biomass : biomass of green vegetative
    dv_biomass : biomass of dead vegetative
    gr_biomass : biomass of green reproductive
    dr_biomass : biomass of dead reproductive
    cutHeight : height of the cut [m]
    rhogv : Volume VV [g m⁻³]
    rhodv : Volume DV [g m⁻³]
    rhogr : Volume GR [g m⁻³]
    rhodr : Volume DR [g m⁻³]
    maxAmountToIngest : The maximum amount of biomass to ingest
        (** USING A DEFAULT OF 9999999 - NEED TO CHECK!)

    Returns
    -------
    - the sum of ingested biomass
    """

    avDefoliationBiomassGV = avDefoliationBiomass(
        biomass=gv_biomass, cutHeight=cutHeight, bulkDensity=rhogv
    )
    avDefoliationBiomassDV = avDefoliationBiomass(
        biomass=dv_biomass, cutHeight=cutHeight, bulkDensity=rhodv
    )
    avDefoliationBiomassGR = avDefoliationBiomass(
        biomass=gr_biomass, cutHeight=cutHeight, bulkDensity=rhogr
    )
    avDefoliationBiomassDR = avDefoliationBiomass(
        biomass=dr_biomass, cutHeight=cutHeight, bulkDensity=rhodr
    )
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
                    biomass=avDefoliationBiomassGV,
                    biomassToIngest=maxAmountToIngest
                )
                + exeDefoliationByBiomass(
                    biomass=avDefoliationBiomassDV,
                    biomassToIngest=maxAmountToIngest
                )
                + exeDefoliationByBiomass(
                    biomass=avDefoliationBiomassGR,
                    biomassToIngest=maxAmountToIngest
                )
                + exeDefoliationByBiomass(
                    biomass=avDefoliationBiomassDR,
                    biomassToIngest=maxAmountToIngest
                )
            )
        else:
            sumBiomassIngested = (
                exeDefoliationByBiomass(
                    biomass=(
                        maxAmountToIngest
                        * avDefoliationBiomassGV / sumAvailable
                    ),
                    biomassToIngest=maxAmountToIngest
                )
                + exeDefoliationByBiomass(
                    biomass=(
                        maxAmountToIngest
                        * avDefoliationBiomassDV / sumAvailable
                    ),
                    biomassToIngest=maxAmountToIngest
                )
                + exeDefoliationByBiomass(
                    biomass=(
                        maxAmountToIngest
                        * avDefoliationBiomassGR / sumAvailable
                    ),
                    biomassToIngest=maxAmountToIngest
                )
                + exeDefoliationByBiomass(
                    biomass=(
                        maxAmountToIngest
                        * avDefoliationBiomassDR / sumAvailable
                    ),
                    biomassToIngest=maxAmountToIngest
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
#     gv_biomass : Biomass of green vegetative
#     gv_gamma : ?

#     Returns
#     -------
#     - Mass of green limbs
#     """

#     return gv_gamma * gv_biomass


# def getOMDgv(gv_min_omd, gv_max_omd, gv_avg_age, lls):
#     """
#     Compute the green vegetative actual organic matter digestibility
#     ** THIS FUNCTION IS UNUSED!

#     Parameters
#     ----------
#     gv_min_omd : The minimum green vegetative organic matter digestibility
#     gv_max_omd : The maximum green vegetative organic matter digestibility
#     gv_avg_age : The average age of the green vegetative
#     lls : the leaf lifespan [°C d]

#     Returns
#     -------
#     - the green vegetative organic matter digestibility
#     """

#     return max(
#         gv_min_omd, gv_max_omd - gv_avg_age * (gv_max_omd - gv_min_omd) / lls
#     )


# def getOMDgr(gr_min_omd, gr_max_omd, gr_avg_age, st1, st2):
#     """
#     Compute the green reproductive actual organic matter digestibility
#     ** THIS FUNCTION IS UNUSED!

#     Parameters
#     ----------
#     gr_min_omd : The minimum green reproductive organic matter digestibility
#     gr_max_omd : The maximum green reproductive organic matter digestibility
#     gr_avg_age : The average age of the green reproductive
#     st1 : Sum of temperature to begin vegetative activity
#     st2 : Sum of temperature to end vegetative activity

#     Returns
#     -------
#     - Green vegetative organic matter digestibility
#     """

#     return max(
#         gr_min_omd,
#         gr_max_omd - gr_avg_age * (gr_max_omd - gr_min_omd) / (st2 - st1)
#     )


def getSumTemperature(timeseries, doy, t0):
    """
    Return the sum temperature corresponding to the DOY

    Parameters
    ----------
    timeseries : Input time series data
    doy : Day of the year [1-366]
    t0 : Minimum temperature for growth [°C]

    Returns
    -------
    - Sum of temperatures above t0 corresponding to the DOY
    """

    sum_temperature = 0
    for i in range(doy):
        if timeseries["tas"][i] > t0:
            sum_temperature += timeseries["tas"][i] - t0
    return sum_temperature


# TO-DO: This set of functions are either not used or not useful
# def addNI(ni, amountToIncrease):
#     return max(0, min(amountToIncrease + ni, 1.2))


def stocking_rate(livestock_units, grazing_area):
    """
    Calculate the stocking rate

    Parameters
    ----------
    livestock_units : total number of livestock units [LU]
    grazing_area : total grazing area (i.e. grassland available for grazing)
        [ha]

    Returns
    -------
    stocking rate [LU ha⁻¹]
    """

    try:
        stocking_rate_ha = livestock_units / grazing_area
    except ZeroDivisionError:
        stocking_rate_ha = 0
    return stocking_rate_ha


def ingested_biomass(
    livestock_units, grazing_area, ingestion_per_livestock_unit=13
):
    """
    Return the amount of biomass ingested by the livestock units based on the
    stocking rate.

    Parameters
    ----------
    livestock_units : total number of livestock units [LU]
    grazing_area : total grazing area (i.e. grassland available for grazing)
        [ha]
    ingestion_per_livestock_unit : average ingestion of grass by a dairy cow
        [kg DM LU⁻¹]; default is 13 based on Teagasc data

    Returns
    -------
    total grass ingestion by livestock per hectare [kg DM ha⁻¹]

    Notes
    -----
    - Grass10: https://www.teagasc.ie/crops/grassland/grass10/
        - average ingestion of grass for dairy cows is 13.41 kg DM LU⁻¹
        - average ingestion of supplements (meal, concentrate, silage) is
        4.25 kg DM LU⁻¹
    - Teagasc Dairy Manual:
      https://www.teagasc.ie/publications/2016/teagasc-dairy-manual.php
        - 8-13 kg DM grass per cow in the spring
        - increase of 0.75-1.0 kg DM until peak intake is reached
        - peak intake of 16-18 kg DM
        - average intake is 13 kg DM ((8 + 18) / 2)
    - one dairy cow is equivalent to one livestock unit; see
      https://cap-calculators.apps.rhos.agriculture.gov.ie/stocking-rate
    """

    total_ingestion = (
        stocking_rate(
            livestock_units=livestock_units, grazing_area=grazing_area
        ) * ingestion_per_livestock_unit
    )
    return total_ingestion
