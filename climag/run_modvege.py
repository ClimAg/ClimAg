"""run_modvege.py

https://github.com/YannChemin/modvege

ModVege main code, had to rewrite most of the functions as the Java code
was a complete mess
This code runs a single geographical "Cell" (as the Java code was trying
to do on a grid)
This function *should* be self sustaining, nothing else needed.
"""

import pandas as pd
# import the model function
from climag.modvege import modvege
# file reading (read_out is ONLY FOR DEV)
from climag.modvege_read_files import read_params, read_weather, read_out


def run_modvege(input_params_csv, input_weather_csv, out_csv, out_dev=None):
    """
    Preprocess the inputs to run ModVege model as a function

    Parameters
    ----------
    input_params_csv : File name for the input parameters CSV
    input_weather_csv : File name for the input weather CSV
    out_csv : File name for the output CSV
    out_dev : Test output data (FOR DEV ONLY)
    """

    # Read parameter files into array
    params = read_params(input_params_csv)

    # Read weather file into array
    # arr[0][0] = DOY[0] = 1
    # arr[0][1] = Temperature[0] = -0.84125
    # arr[0][2] = PARi[0] = 2.22092475
    # arr[0][3] = PP[0] = 0.119
    # arr[0][4] = PET[0] = 0.602689848
    # arr[0][5] = ETA[0] = 0.301344 # RS simulated
    # arr[0][6] = LAI[0] = 0.864162 # RS simulated
    # arr[0][7] = gcut_height[0] = 0.0 [default is 0.05 if cut]
    # arr[0][8] = grazing_animal_count[0] = 0 [default is 1 for test]
    # arr[0][9] = grazing_avg_animal_weight[0] = 0 [default is 400 for cow]

    weather = read_weather(input_weather_csv)

    startdoy = 1
    enddoy = 365

    # Initialize the run and return arrays
    (
        gv_b, dv_b, gr_b, dr_b, h_b, i_b, gro, abc, sumT,
        gva, gra, dva, dra, sea, ftm, env, pgr, atr
    ) = modvege(params, weather, startdoy, enddoy)

    # convert output to dataframe and save as CSV
    data = zip(
        list(range(1, len(gv_b) + 1)), gv_b, dv_b, gr_b, dr_b, h_b, i_b,
        gro, abc, sumT, gva, gra, dva, dra, sea, ftm, env, pgr, atr
    )

    colnames = [
        "doy", "gv_b", "dv_b", "gr_b", "dr_b", "h_b", "i_b", "gro", "abc",
        "sumT", "gva", "gra", "dva", "dra", "sea", "ftm", "env", "pgr", "atr"
    ]

    output_df = pd.DataFrame(data, columns=colnames)

    output_df.to_csv(out_csv, index=False)

    # Print the output
    # print(output)

    # ###############################################  ###################
    # Definition of columns in out_cut.csv             Eq. from output run
    # ###############################################  ###################
    # 0 day
    # 1 Mean biomass                     [kg DM ha-1]  gv_b+gr_b+dv_b+dr_b
    # 2 Mean green vegetative biomass    [kg DM ha-1]  gv_b
    # 3 Mean green reproductive biomass  [kg DM ha-1]  gr_b
    # 4 Mean dry vegetative biomass      [kg DM ha-1]  dv_b
    # 5 Mean dry reproductive biomass    [kg DM ha-1]  dr_b
    # 6 Harvested Biomass                [kg DM ha-1]  h_b
    # 7 Ingested Biomass                 [kg DM ha-1]  i_b
    # 8 Mean GRO biomass                 [kg DM ha-1]  gro
    # 9 Mean available biomass for cut   [kg DM ha-1]  abc

    # ONLY FOR DEV
    if out_dev is not None:
        import matplotlib.pyplot as plt
        import climag.plot_configs

        out = read_out(out_dev)

        # PLOT
        out_doy = [out[i][0] for i in range(len(out) - 1)]
        out_gvb = [out[i][2] for i in range(len(out) - 1)]
        out_grb = [out[i][3] for i in range(len(out) - 1)]
        out_dvb = [out[i][4] for i in range(len(out) - 1)]
        out_drb = [out[i][5] for i in range(len(out) - 1)]
        out_hb = [out[i][6] for i in range(len(out) - 1)]
        out_ib = [out[i][7] for i in range(len(out) - 1)]
        out_gro = [out[i][8] for i in range(len(out) - 1)]
        out_abc = [out[i][9] for i in range(len(out) - 1)]

        plt.figure(figsize=(15, 7))

        plt.subplot(331)
        plt.plot(out_doy, gv_b, label="gv_b")
        plt.plot(out_doy, out_gvb, label="out_gvb")
        plt.title("Green vegetative biomass (kg DM/ha)")
        plt.legend()

        plt.subplot(332)
        plt.plot(out_doy, gr_b, label="gr_b")
        plt.plot(out_doy, out_grb, label="out_grb")
        plt.title("Green reproductive biomass (kg DM/ha)")
        plt.legend()

        plt.subplot(333)
        plt.plot(out_doy, sumT, label="sumT")
        plt.plot(out_doy, gva, label="gv_age")
        plt.plot(out_doy, gra, label="gr_age")
        plt.plot(out_doy, dva, label="dv_age")
        plt.plot(out_doy, dra, label="dr_age")
        plt.title("Sum of temperatures (°C)")
        plt.legend()

        plt.subplot(334)
        plt.plot(out_doy, dv_b, label="dv_b")
        plt.plot(out_doy, out_dvb, label="out_dvb")
        plt.title("Dead vegetative biomass (kg DM/ha)")
        plt.legend()

        plt.subplot(335)
        plt.plot(out_doy, dr_b, label="dr_b")
        plt.plot(out_doy, out_drb, label="out_drb")
        plt.title("Dead reproductive biomass (kg DM/ha)")
        plt.legend()

        plt.subplot(336)
        plt.plot(out_doy, pgr, label="pot. growth")
        plt.plot(out_doy, gro, label="gro")
        plt.plot(out_doy, out_gro, label="out_gro")
        plt.title("GRO biomass (kg DM/ha)")
        plt.legend()

        plt.subplot(337)
        plt.plot(out_doy, abc, label="abc")
        plt.plot(out_doy, out_abc, label="out_abc")
        plt.title("Mean available biomass for cut (kg DM/ha)")
        plt.legend()

        # Harvested Biomass Plot
        plt.subplot(338)
        plt.plot(out_doy, h_b, label="h_b")
        plt.plot(out_doy, out_hb, label="out_hb")
        plt.title("Harvested biomass (kg DM/ha)")
        plt.legend()

        plt.subplot(339)
        plt.plot(out_doy, atr, label="a2r")
        plt.plot(out_doy, sea, label="Season")
        plt.plot(out_doy, ftm, label="Temperature")
        plt.plot(out_doy, env, label="Environmental")
        plt.title("ENV and other factors")
        plt.legend()

        plt.tight_layout()
        plt.show()
