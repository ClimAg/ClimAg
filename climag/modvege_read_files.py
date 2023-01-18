"""modvege_read_files.py

https://github.com/YannChemin/modvege

ModVege has four compartments:
Green vegetative          (GV)
Green reproductive        (GR)
Dead vegetative           (DV)
Dead reproductive         (DR)

References
----------
- Jouven, M., Carrère, P. and Baumont, R. (2006). 'Model predicting dynamics
  of biomass, structure and digestibility of herbage in managed permanent
  pastures. 1. Model description', Grass and Forage Science, vol. 61, no. 2,
  pp. 112-124. DOI: 10.1111/j.1365-2494.2006.00515.x.
"""

# import numpy as np
import pandas as pd


def read_params(filename):
    """
    Read the input parameters (constants) file

    See Tables 2 and 3 in Jouven et al. (2006) for estimates of these
    parameters. Temperate grasses have been classified into four groups based
    on their functional traits. The four groups have been parameterised for
    the Auvergne region in France, which has a temperate climate.

    Functional group A is relevant to Ireland: Species found in fertile sites,
    adapted to frequent defoliation **perennial ryegrass (*Lolium perenne*)**;
    has high specific leaf area, high digestibility, short leaf lifespan,
    early reproductive growth and flowering.

    The nutritional index (NI) is site-specific.

    Definition of inputs
    --------------------
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
    - cut_height  : Grass cut height; see sec. "Harvested biomass" in Jouven
                    et al. (2005)  [0.05 m]
    - RUEmax      : Maximum radiation use efficiency (RUE_max) [3 g DM MJ⁻¹]
    - maxOMDgv    : Maximum organic matter digestibility of the green
                    vegetative compartment (maxOMD_GV) [0.90]
    - minOMDgv    : Minimum organic matter digestibility of the green
                    vegetative compartment (minOMD_GV) [0.75]
    - maxOMDgr    : Maximum organic matter digestibility of the green
                    reproductive compartment (maxOMD_GR) [0.90]
    - minOMDgr    : Minimum organic matter digestibility of the green
                    reproductive compartment (minOMD_GR) [0.65]
    - meanOMDdv   : Organic matter digestibility for the dead vegetative
                    compartment (OMD_DV) [0.45]
    - meanOMDdr   : Organic matter digestibility for the dead reproductive
                    compartment (OMD_DR) [0.40]
    - n_index     : Nitrogen nutritional index
    - whc         : Soil water-holding capacity (WHC) [mm]
    - wr_init     : Initial water reserves (WR) [mm]
    - bm_gv_init  : Initial biomass of GV [kg DM ha⁻¹]
    - bm_gr_init  : Initial biomass of GR [kg DM ha⁻¹]
    - bm_dv_init  : Initial biomass of DV [kg DM ha⁻¹]
    - bm_dr_init  : Initial biomass of DR [kg DM ha⁻¹]
    - age_gv_init : Initial GV age [°C d]
    - age_gr_init : Initial GR age [°C d]
    - age_dv_init : Initial DV age [°C d]
    - age_dr_init : Initial DR age [°C d]
    - livestock_units
    - grazing_area

    Parameters
    ----------
    filename : Path to the parameter input file

    Returns
    -------
    - A dictionary of the input parameters
    """

    params = pd.read_csv(
        filename, header=None, index_col=0
    ).squeeze().to_dict()
    return params


def read_timeseries(filename):
    """
    Read the time series input data

    Definition of inputs
    --------------------
    - day        : Day number
    - tas        : Temperature (*T*) [°C]
    - par        : Incident photosynthetically active radiation
                   (PAR_i) [MJ m⁻²]
    - pr         : Precipitation (PP) [mm]
    - evspsblpot : (Potential or reference) evapotranspiration (ET) [mm]
    - gcut       : Grass cut height [m] (if cut, the default is 0.05)
    - grazing    : Grazing animal count
    - grazingw   : Grazing animal average weight [kg]

    Parameters
    ----------
    filename : Path to the input time series data file

    Returns
    -------
    - A dataframe of the input time series data
    - Length of the data (number of days)
    """

    timeseries = pd.read_csv(filename, parse_dates=["time"])
    timeseries.sort_values(by=["time"], inplace=True)
    timeseries.reset_index(inplace=True)
    endday = len(timeseries)
    return timeseries, endday
