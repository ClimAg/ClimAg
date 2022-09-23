"""modvege_read_files.py

https://github.com/YannChemin/modvege

ModVege is using four parts:
Green Vegetative      (GV)
Green Reproductive    (GR)
Dead Vegetative       (DV)
Dead Reproductive     (DR)
"""

import numpy as np


def read_params(filename):
    """
    Read the input parameters CSV file ("params.csv")

    Definition of input parameters
    ------------------------------
    ST1            Onset of reproductive growth [degree day]
    ST2            End of reproductive growth [degree day]
    INcell         Initial Nutritional index of cell - NNI
    WHC            Soil water-holding capacity [mm]
    WR             Soil water reserve [mm]
    minSEA         Growth increase in winter
    maxSEA         Growth increase in summer
    W_GV           Initial Biomass of GV [kg ha-1]
    alpha_PAR      Light Extinction Coefficient
    T0             Temperature threshold: photosynthesis activation [°C]
    T1             Temp threshold: stable growth [°C]
    T2             Temp threshold: growth decline [°C]
    beta_T         Decrease in LUE after T2
    b_IN           Impact of IN on LUE at IN=0
    SLA            Specific leaf area [m2 g-1]
    LLS            Leaf lifespan [degree day]
    rho_GV         Volume GV [g m-3]
    percentLAM     Fraction of leaf of laminae in GV
    W_GR           Biomass of GR [kg ha-1]
    a_IN           Value of ALLOC at IN=0
    max_fIN        Max of fNI
    rho_GR         Volume GR [g m-3]
    W_DV           Biomass of DV [kg ha-1]
    K_DV           Senescence coefficient DV [degree day]
    Kl_DV          Abscission coefficient DV [degree day]
    rho_DV         Volume DV [g m-3]
    W_DR           Biomass of DR [kg ha-1]
    K_DR           Senescence coefficient DR [degree day]
    Kl_DR          Abscission coefficient DR [degree day]
    rho_DR         Volume DR [g m-3]
    init_AGE_GV    Initial value of age GV
    init_AGE_GR    Initial value of age VR
    init_AGE_DV    Initial value of age DV
    init_AGE_DR    Initial value of age DR
    RUEmax         Maximum radiation use efficiency [g MJ-1]
    gammaGV        Respiratory C loss during senescence (DV)
    gammaGR        Respiratory C loss during senescence (DR)
    maxOMDgv       maximum OMD green veg
    minOMDgv       minimum OMD green veg
    maxOMDgr       maximum OMD green rep
    minOMDgr       minimum OMD green rep
    meanOMDdv      mean OMD dry veg (digestibility of dead part is constant)
    meanOMDdr      mean OMD dry rep (digestibility of dead part constant)
    cellSurface    Pixel area [ha]

    Parameters
    ----------
    file : the input file named param.csv

    Returns
    -------
    arr : the returning array of [DOY, Temperature, PARi, PP, PET]
    """
    arr = np.genfromtxt(filename, delimiter=",", dtype=float, usecols=(-1))
    return arr


def read_weather(filename):
    """
    Read the weather CSV file ("weather.csv")

    Definition of input parameters
    ------------------------------
    DOY            day of year
    Temperature    temperature [°C]
    PARi           photosynthetic radiation incident [MJ m-2]
    PP             Precipitation [mm]
    PET            Potential ET [mm/day]
    ETA            Actual ET from remote sensing [mm/day]
    LAI            Leaf Area Index from remote sensing
    gcut           Grass cut event cutHeight [m]
    grazing        Grazing animal count
    grazingw       Grazing average animal weight [kg]

    Parameters
    ----------
    file : the input file named weather.csv

    Returns
    -------
    arr : the returning array of [DOY, Temperature, PARi, PP, PET]
    """
    arr = np.genfromtxt(filename, delimiter=",", skip_header=0, names=True)
    return arr


def read_out(filename):
    """Read the "out_cut.csv" file

    This is used only for development, to remove for operational mode

    Definition of columns in out_cut.csv
    ------------------------------------
    day
    Mean biomass                      [kg DM ha-1]
    Mean green vegetative biomass     [kg DM ha-1]
    Mean green reproductive biomass   [kg DM ha-1]
    Mean dry vegetative biomass       [kg DM ha-1]
    Mean dry reproductive biomass     [kg DM ha-1]
    Harvested Biomass                 [kg DM ha-1]
    Ingested Biomass                  [kg DM ha-1]
    Mean GRO biomass                  [kg DM ha-1]
    Mean available biomass for cut    [kg DM ha-1]
    """
    arr = np.genfromtxt(filename, delimiter=",", skip_header=0, names=True)
    return arr
