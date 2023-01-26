"""modvege_read_files.py

https://github.com/YannChemin/modvege

ModVege has four compartments:
Green vegetative          (GV)
Green reproductive        (GR)
Dead vegetative           (DV)
Dead reproductive         (DR)
"""

import pandas as pd


def read_params(filename):
    """
    Read the input parameters (constants) file

    See Tables 2 and 3 in Jouven et al. (2006a) for estimates of these
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
    - sla         : Specific leaf area (SLA) [0.033 m² g⁻¹]
    - pct_lam     : Percentage of laminae in the green vegetative compartment
                    (%LAM) [0.68]
    - st_1        : Sum of temperatures at the beginning of the reproductive
                    period (ST₁) [600 °C d]
    - st_2        : Sum of temperatures at the end of the reproductive period
                    (ST₂) [1200 °C d]
    - max_sea     : Maximum seasonal effect (maxSEA) [1.20]
    - min_sea     : Minimum seasonal effect (minSEA) [0.80]
    - lls         : Leaf lifespan (LLS) [500 °C d]
    - bd_gv       : Bulk density of the green vegetative compartment (BD_GV)
                    [850 g DM m⁻³]
    - bd_dv       : Bulk density of the dead vegetative compartment (BD_DV)
                    [500 g DM m⁻³]
    - bd_gr       : Bulk density of the green reproductive compartment (BD_GR)
                    [300 g DM m⁻³]
    - bd_dr       : Bulk density of the dead reproductive compartment (BD_DR)
                    [150 g DM m⁻³]
    - sigma_gv    : Rate of biomass loss with respiration for the green
                    vegetative compartment (σ_GV) [0.4]
    - sigma_gr    : Rate of biomass loss with respiration for the green
                    reproductive compartment (σ_GR) [0.2]
    - t_0         : Minimum temperature for growth (T₀) [4 °C]
    - t_1         : Minimum temperature for optimal growth (T₁) [10 °C]
    - t_2         : Maximum temperature for optimal growth (T₂) [20 °C]
    - t_max       : Maximum temperature for growth (T_max) [40 °C]
    - k_gv        : Basic senescence rate for the green vegetative compartment
                    (K_GV) [0.002]
    - k_gr        : Basic senescence rate for the green reproductive
                    compartment (K_GR) [0.001]
    - kl_dv       : Basic abscission rate for the dead vegetative compartment
                    (Kl_DV) [0.001]
    - kl_dr       : Basic abscission rate for the dead reproductive
                    compartment (Kl_DR) [0.0005]
    - h_grass     : Minimum residual grass height after cutting; see sec.
                    "Harvested biomass" in Jouven et al. (2005)  [0.05 m]
    - rue_max     : Maximum radiation use efficiency (RUE_max) [3 g DM MJ⁻¹]
    - max_omd_gv  : Maximum organic matter digestibility of the green
                    vegetative compartment (maxOMD_GV) [0.90]
    - min_omd_gv  : Minimum organic matter digestibility of the green
                    vegetative compartment (minOMD_GV) [0.75]
    - max_omd_gr  : Maximum organic matter digestibility of the green
                    reproductive compartment (maxOMD_GR) [0.90]
    - min_omd_gr  : Minimum organic matter digestibility of the green
                    reproductive compartment (minOMD_GR) [0.65]
    - omd_dv      : Organic matter digestibility for the dead vegetative
                    compartment (OMD_DV) [0.45]
    - omd_dr      : Organic matter digestibility for the dead reproductive
                    compartment (OMD_DR) [0.40]
    - ni          : Nitrogen nutritional index
    - whc         : Soil water-holding capacity (WHC) [mm]
    - bm_gv       : Initial biomass of GV [kg DM ha⁻¹]
    - bm_gr       : Initial biomass of GR [kg DM ha⁻¹]
    - bm_dv       : Initial biomass of DV [kg DM ha⁻¹]
    - bm_dr       : Initial biomass of DR [kg DM ha⁻¹]
    - age_gv      : Initial GV age [°C d]
    - age_gr      : Initial GR age [°C d]
    - age_dv      : Initial DV age [°C d]
    - age_dr      : Initial DR age [°C d]
    - lu          : Number of livestock units in a grazing area [LU]
    - area        : Total grazing area (i.e. grassland available for grazing)
                    [ha]
    - i_bm_lu     : Maximum biomass ingestion per livestock unit
                    [13 kg DM LU⁻¹]

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
    - day   : Day number
    - T     : Temperature (*T*) [°C]
    - PAR_i : Incident photosynthetically active radiation (PAR_i) [MJ m⁻²]
    - PP    : Precipitation (PP) [mm]
    - PET   : (Potential or reference) evapotranspiration (ET) [mm]

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
