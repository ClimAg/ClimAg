"""modvege_read_files.py

https://github.com/YannChemin/modvege

ModVege has four compartments:
Green vegetative      (GV)
Green reproductive    (GR)
Dead vegetative       (DV)
Dead reproductive     (DR)
"""

import numpy as np


def read_params(filename):
    """Read the input parameters CSV file

    See Tables 2 and 3 in Jouven et al. (2006) for estimates of these
    parameters. Temperate grasses have been classified into four groups based
    on their functional traits. The four groups have been parameterised for
    the Auvergne region in France, which has a temperate climate.

    Functional group A is relevant to Ireland: Species found in fertile sites,
    adapted to frequent defoliation **perennial ryegrass (*Lolium perenne*)**;
    has high specific leaf area, high digestibility, short leaf lifespan,
    early reproductive growth and flowering.

    The nutritional index (NI) is site-specific.

    Definition of input parameters
    ------------------------------
    - SLA         : Specific leaf area (SLA) [0.033 m² g⁻¹]
    - pctLAM      : Percentage of laminae in the green vegetative compartment
                    (%LAM) [0.68]
    - ST1         : Sum of temperatures at the beginning of the reproductive
                    period (ST₁) [600 °C d]
    - ST2         : Sum of temperatures at the end of the reproductive period
                    (ST₂) [1200 °C d]
    - maxSEA      : Maximum seasonal effect (maxSEA) [1.20]
    - minSEA      : Minimum seasonal effect (minSEA) [0.80]
    - LLS         : Leaf lifespan (LLS) [500 °C d]
    - maxOMDgv    : Maximum organic matter digestibility of the green
                    vegetative compartment (maxOMD_GV) [0.90]
    - minOMDgv    : Minimum organic matter digestibility of the green
                    vegetative compartment (minOMD_GV) [0.75]
    - maxOMDgr    : Maximum organic matter digestibility of the green
                    reproductive compartment (maxOMD_GR) [0.90]
    - minOMDgr    : Minimum organic matter digestibility of the green
                    reproductive compartment (minOMD_GR) [0.65]
    - rho_GV      : Bulk density of the green vegetative compartment (BD_GV)
                    [850 g DM m⁻³]
    - rho_DV      : Bulk density of the dead vegetative compartment (BD_DV)
                    [500 g DM m⁻³]
    - rho_GR      : Bulk density of the green reproductive compartment (BD_GR)
                    [300 g DM m⁻³]
    - rho_DR      : Bulk density of the dead reproductive compartment (BD_DR)
                    [150 g DM m⁻³]
    - sigmaGV     : Rate of biomass loss with respiration for the green
                    vegetative compartment (σ_GV) [0.4]
    - sigmaGR     : Rate of biomass loss with respiration for the green
                    reproductive compartment (σ_GR) [0.2]
    - T0          : Minimum temperature for growth (T₀) [4 °C]
    - T1          : Minimum temperature for optimal growth (T₁) [10 °C]
    - T2          : Maximum temperature for optimal growth (T₂) [20 °C]
    - K_GV        : Basic senescence rate for the green vegetative compartment
                    (K_GV) [0.002]
    - K_GR        : Basic senescence rate for the green reproductive
                    compartment (K_GR) [0.001]
    - Kl_DV       : Basic abscission rate for the dead vegetative compartment
                    (Kl_DV) [0.001]
    - Kl_DR       : Basic abscission rate for the dead reproductive
                    compartment (Kl_DR) [0.0005]
    - meanOMDdv   : Organic matter digestibility for the dead vegetative
                    compartment (OMD_DV) [0.45]
    - meanOMDdr   : Organic matter digestibility for the dead reproductive
                    compartment (OMD_DR) [0.40]

    - INcell      : Initial nutritional index of cell - NNI
    - WHC         : Soil water-holding capacity (WHC) [mm]
    - WR          : Water reserve (WR) [mm]
    - W_GV        : Initial biomass of GV [kg ha⁻¹]
    - alpha_PAR   : Light extinction coefficient
    - beta_T      : Decrease in LUE after T2
    - b_IN        : Impact of IN on LUE at IN=0
    - W_GR        : Biomass of GR [kg ha⁻¹]
    - a_IN        : Value of ALLOC at IN=0
    - max_fIN     : Max of fNI
    - W_DV        : Biomass of DV [kg ha⁻¹]
    - W_DR        : Biomass of DR [kg DM ha⁻¹]
    - init_AGE_GV : Initial GV age [°C d]
    - init_AGE_GR : Initial GR age [°C d]
    - init_AGE_DV : Initial DV age [°C d]
    - init_AGE_DR : Initial DR age [°C d]
    - RUEmax      : Maximum radiation use efficiency (RUE_max) [g MJ⁻¹]
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
    - DOY         : Day of the year
    - Temperature : Temperature (T) [°C]
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
    - Mean dead vegetative biomass      [kg DM ha⁻¹]
    - Mean dead reproductive biomass    [kg DM ha⁻¹]
    - Harvested biomass                 [kg DM ha⁻¹]
    - Ingested biomass                  [kg DM ha⁻¹]
    - Mean GRO biomass                  [kg DM ha⁻¹]
    - Mean available biomass for cut    [kg DM ha⁻¹]
    """
    arr = np.genfromtxt(filename, delimiter=",", skip_header=0, names=True)
    return arr
