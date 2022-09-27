"""modvege_run.py

https://github.com/YannChemin/modvege

ModVege main code.
This code runs a single geographical "cell" (as the Java code was trying
to do on a grid).
This function *should* be self sustaining, nothing else needed.
"""

import matplotlib.pyplot as plt
import pandas as pd
# import the model function
from climag.modvege import modvege
# file reading
from climag.modvege_read_files import read_params, read_weather


def run_modvege(
    input_params_csv, input_weather_csv, out_csv, startdoy=1, enddoy=365
):
    """
    Preprocess the inputs to run ModVege as a function and save the results
    as a CSV file

    Parameters
    ----------
    input_params_csv : File name for the input parameters CSV
    input_weather_csv : File name for the input weather CSV
    out_csv : File name for the output CSV
    startdoy : start day of the year (default 1)
    enddoy : end day of the year (default 365)
    """

    # read parameter files into array
    params = read_params(input_params_csv)

    # read weather file into array
    # arr[0][0] = DOY[0] = 1
    # arr[0][1] = Temperature[0] = -0.84125
    # arr[0][2] = PARi[0] = 2.22092475
    # arr[0][3] = PP[0] = 0.119
    # arr[0][4] = PET[0] = 0.602689848
    # arr[0][5] = ETA[0] = 0.301344  # RS simulated
    # arr[0][6] = LAI[0] = 0.864162  # RS simulated
    # arr[0][7] = gcut_height[0] = 0.0  # [default is 0.05 if cut]
    # arr[0][8] = grazing_animal_count[0] = 0  # [default is 1 for test]
    # arr[0][9] = grazing_avg_animal_weight[0] = 0  # [default is 400 for cow]

    weather = read_weather(input_weather_csv)

    # initialise the run and return arrays
    mod_out = modvege(params, weather, startdoy, enddoy)

    # convert output to dataframe and save as CSV
    data = tuple([list(range(1, len(mod_out[0]) + 1))]) + mod_out

    colnames = [
        "doy", "gv_b", "dv_b", "gr_b", "dr_b", "h_b", "i_b", "gro", "abc",
        "sumT", "gva", "gra", "dva", "dra", "sea", "ftm", "env", "pgr", "atr"
    ]

    output_df = pd.DataFrame(zip(*data), columns=colnames)

    output_df.to_csv(out_csv, index=False)

    # ###############################################  ###################
    # Definition of columns in out_cut.csv             Eq. from output run
    # ###############################################  ###################
    # 0 day
    # 1 Mean biomass                     [kg DM ha⁻¹]  gv_b+gr_b+dv_b+dr_b
    # 2 Mean green vegetative biomass    [kg DM ha⁻¹]  gv_b
    # 3 Mean green reproductive biomass  [kg DM ha⁻¹]  gr_b
    # 4 Mean dry vegetative biomass      [kg DM ha⁻¹]  dv_b
    # 5 Mean dry reproductive biomass    [kg DM ha⁻¹]  dr_b
    # 6 Harvested Biomass                [kg DM ha⁻¹]  h_b
    # 7 Ingested Biomass                 [kg DM ha⁻¹]  i_b
    # 8 Mean GRO biomass                 [kg DM ha⁻¹]  gro
    # 9 Mean available biomass for cut   [kg DM ha⁻¹]  abc

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
        "Available biomass for cut [kg DM ha⁻¹]",
        "Sum of temperatures [°C d]",
        "GV biomass age [°C d]",
        "GR biomass age [°C d]",
        "DV biomass age [°C d]",
        "DR biomass age [°C d]",
        "Seasonal effect*",
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

    # # biomass compartments
    # output_df.plot(y=["gv_b", "dv_b", "gr_b", "dr_b"], figsize=(15, 4))
    # plt.title("Biomass compartments")
    # plt.xlabel("Day of the year")
    # plt.ylabel("[kg DM ha⁻¹]")

    # # growth
    # output_df.plot(y=["pgr", "gro"], figsize=(15, 4))
    # plt.title("Potential growth and growth of biomass")
    # plt.xlabel("Day of the year")
    # plt.ylabel("[kg DM ha⁻¹]")

    # # age and sum of temperatures
    # output_df.plot(y=["sumT", "gva", "dva", "gra", "dra"], figsize=(15, 4))
    # plt.title("Sum of temperatures and age of biomass compartments")
    # plt.xlabel("Day of the year")
    # plt.ylabel("Degree day [°C d]")

    plt.show()
