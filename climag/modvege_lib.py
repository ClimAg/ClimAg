"""modvege_lib.py

https://github.com/YannChemin/modvege

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


def harvest_biomass(bulkDensity, biomass, cutHeight=0.05):
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
    See Equation (18) and Figure 4(c) in in Jouven et al. (2006).

    Note that abscission only occurs when T > 0.

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
    if dv_avg_age / lls < 1 / 3:
        age = 1
    elif dv_avg_age / lls < 2 / 3:
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
    See Equation (18) and Figure 4(d) in Jouven et al. (2006).

    Note that abscission only occurs when T > 0.

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
    if dr_avg_age / (st2 - st1) < 1 / 3:
        age = 1
    elif dr_avg_age / (st2 - st1) < 2 / 3:
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
    Senescent biomass for GV compartment.
    See Equations (16) and (17) and Figure 4(a) in Jouven et al. (2006).

    No senescence occurs when *T* is between zero and T0.
    When T drops below zero, senescence is driven by freezing effects and is
    proportional to |*T*|.

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
    if gv_avg_age / lls < 1 / 3:
        age = 1
    elif gv_avg_age / lls < 1:
        # linear gradient
        gradient = (3 - 1) / (1 - 1 / 3)
        intercept = 3 - gradient * 1
        age = gradient * gv_avg_age / lls + intercept
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
    gro : total biomass growth (GRO) [kg DM ha⁻¹]
    a2r : Allocate to reproductive
        (REP in Jouven et al. (2006), reproductive function)
    lls : Leaf lifespan (LLS) [500 °C d]
    temperature : Mean daily temperature (*T*) [°C]
    kgv : Senescence coefficient GV [°C d]
    t0 : Minimum temperature for growth (*T*₀) [4 °C]
    gv_biomass : Updated biomass [kg DM ha⁻¹]
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
    See Equations (16) and (17) and Figure 4(b) in Jouven et al. (2006).

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

    if gr_avg_age / (st2 - st1) < 1 / 3:
        age = 1
    elif gr_avg_age / (st2 - st1) < 1:
        # linear gradient
        gradient = (3 - 1) / (1 - 1 / 3)
        intercept = 3 - gradient * 1
        age = gradient * gr_avg_age / (st2 - st1) + intercept
    else:
        age = 3

    if temperature > t0:
        # Equation (16)
        senescence_biomass = kgr * gr_biomass * temperature * age
    elif temperature < 0:
        # Equation (17)
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


# def cut(
#     cutHeight, rhogv, rhodv, rhogr, rhodr, gvb, dvb, grb, drb
# ):
#     """
#     Realise the harvest on each compartment.

#     Parameters
#     ----------
#     cutHeight : the height of the cut [m]
#     rhogv : bulk density of green vegetative [g m⁻³]
#     rhodv : bulk density of dead vegetative [g m⁻³]
#     rhogr : bulk density of green reproductive [g m⁻³]
#     rhodr : bulk density of dead reproductive [g m⁻³]
#     gvb : standing biomass of green vegetative [kg DM ha⁻¹]
#     dvb : standing biomass of dead vegetative [kg DM ha⁻¹]
#     grb : standing biomass of green reproductive [kg DM ha⁻¹]
#     drb : standing biomass of dead reproductive [kg DM ha⁻¹]

#     Returns
#     -------
#     - the total amount of biomass cut [kg DM ha⁻¹]
#     - the amount of GV biomass cut [kg DM ha⁻¹]
#     - the amount of DV biomass cut [kg DM ha⁻¹]
#     - the amount of GR biomass cut [kg DM ha⁻¹]
#     - the amount of DR biomass cut [kg DM ha⁻¹]
#     """

#     gv_h, gv_b = harvest_biomass(bulkDensity=rhogv, cutHeight=cutHeight, biomass=gvb)
#     dv_h, dv_b = harvest_biomass(bulkDensity=rhodv, cutHeight=cutHeight, biomass=dvb)
#     gr_h, gr_b = harvest_biomass(bulkDensity=rhogr, cutHeight=cutHeight, biomass=grb)
#     dr_h, dr_b = harvest_biomass(bulkDensity=rhodr, cutHeight=cutHeight, biomass=drb)
#     sumBiomassHarvested = gv_h + dv_h + gr_h + dr_h
#     return (sumBiomassHarvested, gv_b, dv_b, gr_b, dr_b)


def environmental_limitation(
    temperature_fn, ni, pari, pet, waterReserve,
    waterHoldingCapacity
):
    """
    Environmental limitation of growth (ENV).

    See Equation (13) of Jouven et al. (2006).

    Parameters
    ----------
    temperature_fn : temperature function (f(T)) [dimensionless]
    ni : Nutritional index of pixel (NI)
    pari : Incident photosynthetically active radiation (PAR_i) [MJ m⁻²]
    pet : Potential evapotranspiration (PET) [mm]
    waterReserve : Water reserves (WR) [mm]
    waterHoldingCapacity : Soil water-holding capacity (WHC) [mm]

    Returns
    -------
    - Environmental limitation of growth (ENV) [dimensionless]
    """

    return (
        temperature_fn * ni * par_function(pari=pari)
        * water_stress_function(
            waterReserve=waterReserve,
            waterHoldingCapacity=waterHoldingCapacity,
            pet=pet
        )
    )


def temperature_function(meanTenDaysT, t0=4, t1=10, t2=20, tmax=40):
    """
    Temperature function, *f*(*T*)

    See Figure 2(b) of Jouven et al. (2006) and the accompanying text for more
    info; *f*(*T*) has been derived based on Schapendonk et al. (1998)

    Assume no growth takes place after a maximum temperature

    Schapendonk, A. H. C. M., Stol, W., van Kraalingen, D. W. G. and Bouman,
    B. A. M. (1998). 'LINGRA, a sink/source model to simulate grassland
    productivity in Europe', European Journal of Agronomy, vol. 9, no. 2,
    pp. 87-100. DOI: 10.1016/S1161-0301(98)00027-6.

    Parameters
    ----------
    meanTenDaysT : Mean of the ten days of temperature [°C]
    t0 : Minimum temperature for growth; default is 4 [°C]
    t1 : Minimum temperature for optimal growth; default is 10 [°C]
    t2 : Maximum temperature for optimal growth; default is 20 [°C]
    tmax : Maximum temperature for growth; default is 40 [°C]

    Returns
    -------
    - Temperature function [dimensionless]
    """

    if meanTenDaysT < t0 or meanTenDaysT >= tmax:
        f_temp = 0
    elif t0 <= meanTenDaysT < t1:
        # linear relationship
        gradient = 1 / (t1 - t0)
        intercept = 1 - gradient * t1
        f_temp = gradient * meanTenDaysT + intercept
    elif t1 <= meanTenDaysT <= t2:
        f_temp = 1
    elif t2 < meanTenDaysT < tmax:
        # linear relationship
        gradient = 1 / (t2 - tmax)
        intercept = 1 - gradient * t2
        f_temp = gradient * meanTenDaysT + intercept
    return f_temp


def seasonal_effect(sumT, maxsea=1.2, minsea=0.8, st2=1200, st1=600):
    """
    Calculate seasonal effect (SEA) on growth, driven by the sum of
    temperatures

    SEA > 1 indicates above-ground stimulation by mobilisation of reserves;
    SEA < 1 indicates growth limitation by storage of reserves

    See Figure 3 of Jouven et al. (2006) and the accompanying paragraphs for
    more info

    minSEA and maxSEA are functional traits arranged symmetrically around 1:
    (minSEA + maxSEA) / 2 = 1

    Parameters
    ----------
    maxsea : Maximum seasonal effect (maxSEA); default is 1.2 [dimensionless]
    minsea : Minimum seasonal effect (minSEA); default is 0.8 [dimensionless]
    sumT : Sum of temperatures (ST) [°C d]
    st1 : Sum of temperatures at the beginning of the reproductive period
        (ST₁); default is 600 [°C d]
    st2 : Sum of temperatures at the end of the reproductive period
        (ST₂); default is 1200 [°C d]

    Returns
    -------
    - Seasonal effect [dimensionless]
    """

    if sumT < 200 or sumT > st2:
        f_sea = minsea
    elif (st1 - 200) <= sumT <= (st1 - 100):
        f_sea = maxsea
    elif 200 <= sumT < (st1 - 200):
        # assume SEA increases linearly from minSEA at 200 °C d to maxSEA
        gradient = (maxsea - minsea) / ((st1 - 200) - 200)
        intercept = minsea - gradient * 200
        f_sea = max(gradient * sumT + intercept, minsea)
    elif (st1 - 100) < sumT <= st2:
        # SEA decreases linearly from maxSEA to minSEA at ST_2
        gradient = (maxsea - minsea) / ((st1 - 100) - st2)
        intercept = minsea - gradient * st2
        f_sea = max(gradient * sumT + intercept, minsea)
    return f_sea


def par_function(pari):
    """
    Incident photosynthetically active radiation (PARi) function (fPARi)
    needed to calculate the environmental limitation of growth (ENV).

    The definition has been derived from Schapendonk et al. (1998).
    This function accounts for the decrease in radiation use efficiency (RUE)
    at light intensities higher than 5 MJ m⁻².

    See Figure 2(a), Equation (13), and the section on "Growth functions" in
    Jouven et al. (2006).

    Parameters
    ----------
    pari : Photosynthetic radiation incident (PAR_i) [MJ m⁻²]

    Returns
    -------
    - PARi function [dimensionless]
    """

    if pari < 5:
        f_pari = 1
    else:
        # linear gradient
        gradient = 1 / (5 - 27.5 / 2)
        intercept = 1 - gradient * 5
        f_pari = max(gradient * pari + intercept, 0)
    return f_pari


def water_reserves(precipitation, water_reserve, actual_et, soil_whc):
    """
    Calculate the water reserves (WR).

    WR vary between zero and the soil water-holding capacity (WHC).
    Precipitation (PP) fill the WHC, increasing WR, while actual
    evapotranspiration (AET) empties it.

    See Equation (14) in Jouven et al. (2006).

    Parameters
    ----------
    precipitation : Precipitation (PP) [mm]
    water_reserve : Water reserve (WR) [mm]
    actual_et : Actual evapotranspiration (AET) [mm]
    soil_whc : Soil water-holding capacity (WHC) [mm]

    Returns
    -------
    - Water reserves [mm]
    """

    water_reserve = min(
        max(0, water_reserve + precipitation - actual_et), soil_whc
    )
    return water_reserve


def water_stress_function(waterReserve, waterHoldingCapacity, pet):
    """
    Water stress function.

    See Figure 2(c) and Equation (14) of Jouven et al. (2006).

    Based on McCall and Bishop-Hurley (2003).

    Parameters
    ----------
    waterReserve : Water reserves (WR) [mm]
    waterHoldingCapacity : Soil water-holding capacity (WHC) [mm]
    pet : Potential evapotranspiration (PET) [mm]

    Returns
    -------
    - water stress function [dimensionless]
    """

    waterStress = min(waterReserve / waterHoldingCapacity, 1)

    if pet < 3.8:
        # linear gradients
        if waterStress < 0.2:
            gradient = 0.8 / 0.2
            f_waterstress = gradient * waterStress
        elif waterStress < 0.4:
            gradient = (0.95 - 0.8) / (0.4 - 0.2)
            intercept = 0.8 - gradient * 0.2
            f_waterstress = gradient * waterStress + intercept
        elif waterStress < 0.6:
            gradient = (1 - 0.95) / (0.6 - 0.4)
            intercept = 1 - gradient * 0.6
            f_waterstress = gradient * waterStress + intercept
        else:
            f_waterstress = 1
    elif pet <= 6.5:
        if waterStress < 0.2:
            gradient = 0.4 / 0.2
            f_waterstress = gradient * waterStress
        elif waterStress < 0.4:
            gradient = (0.7 - 0.4) / (0.4 - 0.2)
            intercept = 0.4 - gradient * 0.2
            f_waterstress = gradient * waterStress + intercept
        elif waterStress < 0.6:
            gradient = (0.9 - 0.7) / (0.6 - 0.4)
            intercept = 0.9 - gradient * 0.6
            f_waterstress = gradient * waterStress + intercept
        elif waterStress < 0.8:
            gradient = (1 - 0.9) / (0.8 - 0.6)
            intercept = 1 - gradient * 0.8
            f_waterstress = 0.5 * waterStress + 0.6
        else:
            f_waterstress = 1
    else:
        f_waterstress = waterStress
    return f_waterstress


def reproductive_function(n_index):
    """
    Reproductive function.

    See Equation (15) in Jouven et al. (2006)

    Parameters
    ----------
    n_index : Nitrogen nutritional index (NI) [dimensionless]

    Returns
    -------
    - Reproductive function [dimensionless]
    """

    rep_fn = (0.25 + ((1 - 0.25) * (n_index - 0.35)) / (1 - 0.35))
    return rep_fn


def potential_growth(pari, lai, ruemax=3):
    """
    Calculate potential growth (PGRO)

    See Equation (12) in Jouven et al. (2006)

    Based on Schapendonk et al. (1998).

    The model extinction coefficient is set to a constant value of 0.6
    according to Schapendonk et al. (1998) and Bonesmo and Bélanger (2002).

    The maximum radiation use efficiency is 3 g DM MJ⁻¹ based on
    Schapendonk et al. (1998).

    Parameters
    ----------
    pari : Incident PAR (PAR_i) [MJ m⁻²]
    ruemax : Maximum radiation use efficiency (RUE_max); default is 3
        [g DM MJ⁻¹]
    lai : Leaf area index (LAI) [dimensionless]

    Returns
    -------
    - potential growth (PGRO) [kg DM ha⁻¹]
    """

    p_gro = pari * ruemax * (1 - np.exp(-0.6 * lai)) * 10
    return p_gro


def total_growth(biomass_growth_pot, env, seasonality):
    """
    Calculate the total biomass growth (GRO)

    See Equation (11) in Jouven et al. (2006)

    Parameters
    ----------
    - biomass_growth_pot : Potential growth (PGRO) [kg DM ha⁻¹]
    - env : environmental limitation of growth (ENV) [dimensionless]
    - seasonality : seasonal effect (SEA) [dimensionless]

    Returns
    -------
    - total biomass growth [kg DM ha⁻¹]
    """

    biomass_growth = biomass_growth_pot * env * seasonality
    return biomass_growth


def leaf_area_index(gv_biomass, pctlam=0.68, sla=0.033):
    """
    Calculate the leaf area index

    Equation (12) in Jouven et al. (2006)

    Parameters
    ----------
    pctlam : Percentage of laminae in GV (%LAM); default is 0.68
        [dimensionless]
    sla : Specific leaf area (SLA); default is 0.033 [m² g⁻¹]
    gv_biomass : GV biomass (BM_GV) [kg DM ha⁻¹]

    Returns
    -------
    - Leaf area index (LAI) [dimensionless]
    """

    return sla * (gv_biomass / 10) * pctlam


def actual_evapotranspiration(pet, lai):
    """
    Calculate the actual evapotranspiration (AET)

    AET is equivalent to potential evapotranspiration (PET) when the cover
    intercepts approximately 0.95 of the incident photosynthetically active
    radiation (PAR), i.e., when the leaf area index (LAI) > 3,
    based on Johnson and Parsons (1985).
    AET is proportional to LAI when the proportion of intercepted radiation
    is lower than 0.95, i.e. LAI < 3.

    See Equation (14) in Jouven et al. (2006)

    pet : Potential evapotranspiration (PET) [mm]]
    lai : Leaf area index (LAI) [dimensionless]

    Returns
    -------
    - Actual evapotranspiration (AET) [mm]
    """

    return min(pet, pet * (lai / 3))


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
    Return the amount of biomass available for cut

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
#         gv_min_omd,
#         gv_max_omd - gv_avg_age * (gv_max_omd - gv_min_omd) / lls
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


def sum_of_temperatures(timeseries, doy, t0=4):
    """
    Return the sum of temperatures for the day of the year

    Parameters
    ----------
    timeseries : Input time series data
    doy : Day of the year [1-366]
    t0 : Minimum temperature for growth; default is 4 [°C]

    Returns
    -------
    - Sum of temperatures above t0 corresponding to the DOY [°C d]

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
    livestock_units, grazing_area, bulk_density, min_cut_height=0.05,
    ingestion_per_livestock_unit=13
):
    """
    Return the amount of biomass ingested by the livestock units based on the
    stocking rate.

    Parameters
    ----------
    livestock_units : total number of livestock units [LU]
    grazing_area : total grazing area (i.e. grassland available for grazing)
        [ha]
    ingestion_per_livestock_unit : average ingestion of grass by a livestock
        unit (e.g. dairy cow); default is 13 based on Teagasc data
        [kg DM LU⁻¹]
    min_cut_height : minimum residual grass height to be maintained; default
        is 0.05 [m]
    bulk_density : bulk density of the biomass compartment [g DM m⁻³]

    Returns
    -------
    - total grass ingestion by livestock per hectare [kg DM ha⁻¹]

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
    # if the total ingestion exceeds the minimum cut height of 5 cm, reduce it
    if total_ingestion > min_cut_height * bulk_density * 10:
        total_ingestion = total_ingestion - min_cut_height * bulk_density * 10
    return total_ingestion
