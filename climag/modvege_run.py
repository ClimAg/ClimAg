"""modvege_run.py

https://github.com/YannChemin/modvege

This code runs a single geographical "cell" (as the Java code was trying
to do on a grid).
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from climag.modvege import modvege
from climag.modvege_read_files import read_params  # read_timeseries


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

    # read parameter files into a dataframe
    params = read_params(filename=input_params_file)

    # if input_timeseries_file.endswith(".csv"):
    #     tseries, enddoy = read_timeseries(filename=input_timeseries_file)
    # else:
    #     open the climate model dataset

    tseries = xr.open_dataset(
        input_timeseries_file,
        # chunks="auto",  # the following operations do not work with Dask
        #                   yet, so chunking is disabled...
        decode_coords="all"
    )

    # assign new variables for the outputs
    tseries["bm_gv"] = xr.full_like(tseries["pr"], fill_value=np.nan)
    tseries["bm_gv"].attrs = {
        "standard_name": "green_vegetative_biomass",
        "long_name": "Green vegetative biomass",
        "units": "kg DM ha⁻¹"
    }

    # create a dictionary to store the timeseries output dataframes
    data_df = {}
    for rlon, rlat in [(20, 20), (21, 21)]:
        tseries_loc = tseries.isel(rlon=rlon, rlat=rlat)

        # ignore NaN cells
        if not tseries_loc["evspsblpot"].isnull().all():
            for year in [2050]:
                tseries_y = tseries_loc.sel(
                    time=slice(f"{year}-01-01", f"{year}-12-31")
                )

                # extract the end day of the year
                enddoy = tseries_y["time"].dt.dayofyear.values.max()

                data_df[f"{rlon}_{rlat}_{year}"] = pd.DataFrame(
                    {"time": tseries_y["time"]}
                )  # create a dataframe using the time array

                # assign the variables to columns
                for v in tseries_y.data_vars:
                    data_df[f"{rlon}_{rlat}_{year}"][v] = tseries_y[v]

                # assign other variables
                data_df[f"{rlon}_{rlat}_{year}"]["pari"] = 2.0
                data_df[f"{rlon}_{rlat}_{year}"]["eta"] = 0.0
                data_df[f"{rlon}_{rlat}_{year}"]["lai"] = 0.0
                data_df[f"{rlon}_{rlat}_{year}"]["gcut_height"] = 0.0
                data_df[f"{rlon}_{rlat}_{year}"]["grazing_animal_count"] = 0.0
                data_df[f"{rlon}_{rlat}_{year}"][
                    "grazing_avg_animal_weight"
                ] = 0.0

                # initialise the run
                data_df[f"{rlon}_{rlat}_{year}"] = modvege(
                    params=params,
                    tseries=data_df[f"{rlon}_{rlat}_{year}"],
                    enddoy=enddoy
                )

                # convert output to dataframe
                data_df[f"{rlon}_{rlat}_{year}"] = tuple([list(
                    range(1, len(data_df[f"{rlon}_{rlat}_{year}"][0]) + 1)
                )]) + data_df[f"{rlon}_{rlat}_{year}"]

                colnames = [
                    "doy", "gv_b", "dv_b", "gr_b", "dr_b", "h_b", "i_b",
                    "gro", "abc", "sumT", "gva", "dva", "gra", "dra",
                    "sea", "ftm", "env", "pgr", "atr"
                ]

                data_df[f"{rlon}_{rlat}_{year}"] = pd.DataFrame(
                    zip(*data_df[f"{rlon}_{rlat}_{year}"]), columns=colnames
                )

                # PLOT
                # plot all columns
                data_df[f"{rlon}_{rlat}_{year}"].set_index(
                    "doy", inplace=True
                )

                plot_title = [
                    "Green vegetative biomass [kg DM ha⁻¹]",
                    "Dead vegetative biomass [kg DM ha⁻¹]",
                    "Green reproductive biomass [kg DM ha⁻¹]",
                    "Dead reproductive biomass [kg DM ha⁻¹]",
                    "Harvested biomass [kg DM ha⁻¹]",
                    "Ingested biomass [kg DM ha⁻¹]",
                    "Biomass growth [kg DM ha⁻¹]",
                    "Available biomass [kg DM ha⁻¹]",
                    "Sum of temperatures [°C d]",
                    "Green vegetative biomass age [°C d]",
                    "Dead vegetative biomass age [°C d]",
                    "Green reproductive biomass age [°C d]",
                    "Dead reproductive biomass age [°C d]",
                    "Seasonal effect",
                    "Temperature function *",
                    "Environmental limitation of growth *",
                    "Potential growth [kg DM ha⁻¹]",
                    "Reproductive function *"
                ]

                data_df[f"{rlon}_{rlat}_{year}"].plot(
                    subplots=True, layout=(6, 3), figsize=(15, 14),
                    xlabel="Day of the year", title=plot_title, legend=False
                )

                plt.tight_layout()

                plt.show()

                # assign the outputs to the main xarray dataset
                tseries["bm_gv"].loc[dict(
                    rlon=tseries_y.coords["rlon"],
                    rlat=tseries_y.coords["rlat"],
                    time=tseries_y.coords["time"]
                )] = np.array(data_df[f"{rlon}_{rlat}_{year}"]["gv_b"])

    tseries.to_netcdf(out_file)
    # print(tseries)
