"""modvege_read_files.py

https://github.com/YannChemin/modvege

ModVege has four compartments:
Green Vegetative      (GV)
Green Reproductive    (GR)
Dead Vegetative       (DV)
Dead Reproductive     (DR)
"""

import numpy as np


def read_params(filename):
    """Read the input parameters CSV file

    See Tables 2 and 3 in Jouven (2006). Temperate grasses have been
    classified into four groups based on their functional traits. The four
    groups have been parameterised and the estimates are given in these
    tables.

    Functional group A is relevant to Ireland: Species found in fertile sites,
    adapted to frequent defoliation **perennial ryegrass (*Lolium perenne*)**;
    has high specific leaf area, high digestibility, short leaf lifespan,
    early reproductive growth and flowering.

    The nutritional index (NI) is site-specific.

    Definition of input parameters
    ------------------------------
    - ST1         : Onset of reproductive growth (ST₁) [°C d]
    - ST2         : End of reproductive growth (ST₂) [°C d]
    - INcell      : Initial Nutritional index of cell - NNI
    - WHC         : Soil water-holding capacity (WHC) [mm]
    - WR          : Water reserve (WR) [mm]
    - minSEA      : Growth increase in winter
    - maxSEA      : Growth increase in summer
    - W_GV        : Initial biomass of GV [kg ha⁻¹]
    - alpha_PAR   : Light extinction coefficient
    - T0          : Temperature threshold: photosynthesis activation (T₀) [°C]
    - T1          : Temp threshold: stable growth [°C]
    - T2          : Temp threshold: growth decline [°C]
    - beta_T      : Decrease in LUE after T2
    - b_IN        : Impact of IN on LUE at IN=0
    - SLA         : Specific leaf area [m² g⁻¹]
    - LLS         : Leaf lifespan [°C d]
    - rho_GV      : Volume GV [g m⁻³]
    - percentLAM  : Fraction of leaf of laminae in GV
    - W_GR        : Biomass of GR [kg ha⁻¹]
    - a_IN        : Value of ALLOC at IN=0
    - max_fIN     : Max of fNI
    - rho_GR      : Volume GR [g m⁻³]
    - W_DV        : Biomass of DV [kg ha⁻¹]
    - K_DV        : Senescence coefficient DV [°C d]
    - Kl_DV       : Abscission coefficient DV [°C d]
    - rho_DV      : Volume DV [g m⁻³]
    - W_DR        : Biomass of DR [kg DM ha⁻¹]
    - K_DR        : Senescence coefficient DR [°C d]
    - Kl_DR       : Abscission coefficient DR [°C d]
    - rho_DR      : Volume DR [g m⁻³]
    - init_AGE_GV : Initial GV age [°C d]
    - init_AGE_GR : Initial GR age [°C d]
    - init_AGE_DV : Initial DV age [°C d]
    - init_AGE_DR : Initial DR age [°C d]
    - RUEmax      : Maximum radiation use efficiency (RUE_max) [g MJ⁻¹]
    - sigmaGV     : Respiratory C loss during senescence (GV) (σ_GV)
    - sigmaGR     : Respiratory C loss during senescence (GR) (σ_GR)
    - maxOMDgv    : Maximum OMD green veg
    - minOMDgv    : Minimum OMD green veg
    - maxOMDgr    : Maximum OMD green rep
    - minOMDgr    : Minimum OMD green rep
    - meanOMDdv   : Mean OMD dry veg (digestibility of dead part is constant)
    - meanOMDdr   : Mean OMD dry rep (digestibility of dead part constant)
    - cellSurface : Pixel area [ha]

    Parameters
    ----------
    filename : Path to the parameter input file

    Returns
    -------
    - An array of the input parameters
    """
    arr = np.genfromtxt(filename, delimiter=",", dtype=float, usecols=(-1))
    return arr


def read_weather(filename):
    """Read the weather data CSV file

    Definition of input parameters
    ------------------------------
    - DOY         : Day of year
    - Temperature : Temperature [°C]
    - PARi        : Incident photosynthetically active radiation
        (PAR_i) [MJ m⁻²]
    - PP          : Precipitation (PP) [mm]
    - PET         : Potential evapotranspiration (PET) [mm/day]
    - ETA         : Actual evapotranspiration (AET); from remote sensing
        [mm/day]
    - LAI         : Leaf area index (LAI); from remote sensing
    - gcut        : Grass cut event cutHeight [m]
    - grazing     : Grazing animal count
    - grazingw    : Grazing average animal weight [kg]

    Parameters
    ----------
    filename : Path to the input weather data file

    Returns
    -------
    - An array of the input weather data
    """
    arr = np.genfromtxt(filename, delimiter=",", skip_header=0, names=True)
    return arr


def read_out(filename):
    """Read the "out_cut.csv" file

    This is used only for development, to remove for operational mode

    Definition of columns
    ---------------------
    - Day
    - Mean biomass                      [kg DM ha⁻¹]
    - Mean green vegetative biomass     [kg DM ha⁻¹]
    - Mean green reproductive biomass   [kg DM ha⁻¹]
    - Mean dry vegetative biomass       [kg DM ha⁻¹]
    - Mean dry reproductive biomass     [kg DM ha⁻¹]
    - Harvested biomass                 [kg DM ha⁻¹]
    - Ingested biomass                  [kg DM ha⁻¹]
    - Mean GRO biomass                  [kg DM ha⁻¹]
    - Mean available biomass for cut    [kg DM ha⁻¹]
    """
    arr = np.genfromtxt(filename, delimiter=",", skip_header=0, names=True)
    return arr
