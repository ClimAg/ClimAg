"""modvege_run.py

https://github.com/YannChemin/modvege

This code runs a single geographical "cell" (as the Java code was trying
to do on a grid).
"""

import matplotlib.pyplot as plt
import pandas as pd
# import the model function
from climag.modvege import modvege
# file reading
from climag.modvege_read_files import read_params, read_timeseries


def run_modvege(input_params_file, input_timeseries_file, out_file):
    """
    Preprocess the inputs to run ModVege as a function and save the results
    as a CSV file

    Definition of columns in the output
    -----------------------------------
    - day                                                   doy
    - Mean green vegetative biomass         [kg DM ha⁻¹]    gv_b
    - Mean green reproductive biomass       [kg DM ha⁻¹]    gr_b
    - Mean dead vegetative biomass          [kg DM ha⁻¹]    dv_b
    - Mean dead reproductive biomass        [kg DM ha⁻¹]    dr_b
    - Harvested biomass                     [kg DM ha⁻¹]    h_b
    - Ingested biomass                      [kg DM ha⁻¹]    i_b
    - Mean GRO biomass                      [kg DM ha⁻¹]    gro
    - Mean available biomass for cut        [kg DM ha⁻¹]    abc
      (gv_b + gr_b + dv_b + dr_b)
    - Sum of temperatures                   [°C d]          sumT
    - GV biomass age                        [°C d]          gva
    - GR biomass age                        [°C d]          gra
    - DV biomass age                        [°C d]          dva
    - DR biomass age                        [°C d]          dra
    - Seasonal effect                                       sea
    - Temperature function*                                 ftm
    - Environmental limitation of growth                    env
    - Potential growth                      [kg DM ha⁻¹]    pgr
    - Reproductive function*                                atr

    Parameters
    ----------
    input_params_file : File path for the input parameters
    input_timeseries_file : File path for the input timeseries
    out_file : File path for the output
    """

    # read parameter files into array
    params = read_params(input_params_file)

    tseries = read_timeseries(input_timeseries_file)

    # initialise the run and return arrays
    mod_out = modvege(params, tseries)

    # convert output to dataframe and save as CSV
    data = tuple([list(range(1, len(mod_out[0]) + 1))]) + mod_out

    colnames = [
        "doy", "gv_b", "dv_b", "gr_b", "dr_b", "h_b", "i_b", "gro", "abc",
        "sumT", "gva", "gra", "dva", "dra", "sea", "ftm", "env", "pgr", "atr"
    ]

    output_df = pd.DataFrame(zip(*data), columns=colnames)

    output_df.to_csv(out_file, index=False)

    # PLOT
    # plot all columns
    output_df.set_index("doy", inplace=True)
    plot_title = [
        "GV biomass [kg DM ha⁻¹]",
        "DV biomass [kg DM ha⁻¹]",
        "GR biomass [kg DM ha⁻¹]",
        "DR biomass [kg DM ha⁻¹]",
        "Harvested biomass [kg DM ha⁻¹]",
        "Ingested biomass [kg DM ha⁻¹]",
        "Biomass growth [kg DM ha⁻¹]",
        "Available biomass [kg DM ha⁻¹]",
        "Sum of temperatures [°C d]",
        "GV biomass age [°C d]",
        "GR biomass age [°C d]",
        "DV biomass age [°C d]",
        "DR biomass age [°C d]",
        "Seasonal effect",
        "Temperature function*",
        "Environmental limitation of growth*",
        "Potential growth [kg DM ha⁻¹]",
        "Reproductive function*"
    ]
    output_df.plot(
        subplots=True, layout=(6, 3), figsize=(15, 14),
        xlabel="Day of the year", title=plot_title, legend=False
    )
    plt.tight_layout()

    plt.show()
