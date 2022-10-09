"""modvege_run.py

https://github.com/YannChemin/modvege

This code runs a single geographical "cell" (as the Java code was trying
to do on a grid).
"""

import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from climag.modvege import modvege
from climag.modvege_read_files import read_params, read_timeseries


def run_modvege(input_params_file, input_timeseries_file, out_file):
    """
    Preprocess the inputs to run ModVege as a function and save the results
    as a CSV file

    Definition of columns in the output
    -----------------------------------
    - Day of the year                                       doy
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

    # read parameter file into a dataframe
    params = read_params(filename=input_params_file)

    outputs = {
        "gv_b": ["Green vegetative biomass", "kg DM ha⁻¹"],
        "dv_b": ["Dead vegetative biomass", "kg DM ha⁻¹"],
        "gr_b": ["Green reproductive biomass", "kg DM ha⁻¹"],
        "dr_b": ["Dead reproductive biomass", "kg DM ha⁻¹"],
        "h_b": ["Harvested biomass", "kg DM ha⁻¹"],
        "i_b": ["Ingested biomass", "kg DM ha⁻¹"],
        "gro": ["Biomass growth", "kg DM ha⁻¹"],
        "abc": ["Available biomass", "kg DM ha⁻¹"],
        "sumT": ["Sum of temperatures", "°C d"],
        "gva": ["Green vegetative biomass age", "°C d"],
        "dva": ["Dead vegetative biomass age", "°C d"],
        "gra": ["Green reproductive biomass age", "°C d"],
        "dra": ["Dead reproductive biomass age", "°C d"],
        "sea": ["Seasonal effect"],
        "ftm": ["Temperature function"],
        "env": ["Environmental limitation of growth"],
        "pgr": ["Potential growth", "kg DM ha⁻¹"],
        "atr": ["Reproductive function"]
    }

    if input_timeseries_file.endswith(".csv"):
        tseries, enddoy = read_timeseries(filename=input_timeseries_file)

        # initialise the run
        data_df = modvege(params=params, tseries=tseries, enddoy=enddoy)

        # convert output to dataframe and save as CSV
        data_df = tuple([list(range(1, len(data_df[0]) + 1))]) + data_df

        data_df = pd.DataFrame(
            zip(*data_df), columns=["doy"] + list(outputs.keys())
        )

        data_df.to_csv(out_file, index=False)

        # PLOT
        # plot all columns
        data_df.set_index("doy", inplace=True)

        plot_title = []
        for val in outputs.values():
            if len(val) > 1:
                val = " [".join(val) + "]"
            else:
                val = val[0]
            plot_title.append(val)

        data_df.plot(
            subplots=True, layout=(6, 3), figsize=(15, 14),
            xlabel="Day of the year", title=plot_title, legend=False
        )

        plt.tight_layout()

        plt.show()
    else:  # open the climate model dataset
        tseries = xr.open_dataset(
            input_timeseries_file,
            # chunks="auto",  # the following operations do not work with Dask
            #                   yet, so chunking is disabled...
            decode_coords="all"
        )

        # assign new variables for the outputs
        for key, val in outputs.items():
            tseries[key] = xr.full_like(tseries["pr"], fill_value=np.nan)
            if len(val) > 1:
                tseries[key].attrs = {
                    "standard_name": val[0].lower().replace(" ", "_"),
                    "long_name": val[0],
                    "units": val[1]
                }
            else:
                tseries[key].attrs = {
                    "standard_name": val[0].lower().replace(" ", "_"),
                    "long_name": val[0],
                    "units": "dimensionless"
                }

        # create a dictionary to store the timeseries output dataframes
        data_df = {}

        # loop through each grid cell
        # for rlon, rlat in [(20, 20), (21, 21)]:
        for rlon, rlat in itertools.product(
            range(len(tseries.coords["rlon"])),
            range(len(tseries.coords["rlat"]))
        ):
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
                    for var in tseries_y.data_vars:
                        data_df[f"{rlon}_{rlat}_{year}"][var] = tseries_y[var]

                    # assign other variables
                    data_df[f"{rlon}_{rlat}_{year}"]["pari"] = 2.0
                    data_df[f"{rlon}_{rlat}_{year}"]["eta"] = 0.0
                    data_df[f"{rlon}_{rlat}_{year}"]["lai"] = 0.0
                    data_df[f"{rlon}_{rlat}_{year}"]["gcut_height"] = 0.0
                    data_df[f"{rlon}_{rlat}_{year}"][
                        "grazing_animal_count"
                    ] = 0.0
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

                    data_df[f"{rlon}_{rlat}_{year}"] = pd.DataFrame(
                        zip(*data_df[f"{rlon}_{rlat}_{year}"]),
                        columns=["doy"] + list(outputs.keys())
                    )

                    # assign the outputs to the main xarray dataset
                    for key in outputs:
                        tseries[key].loc[dict(
                            rlon=tseries_y.coords["rlon"],
                            rlat=tseries_y.coords["rlat"],
                            time=tseries_y.coords["time"]
                        )] = np.array(data_df[f"{rlon}_{rlat}_{year}"][key])

        tseries = tseries.drop_vars(["evspsblpot", "pr", "tas"])
        tseries.to_netcdf(out_file)
