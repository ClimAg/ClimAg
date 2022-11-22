"""modvege_run.py

https://github.com/YannChemin/modvege

This code runs a single geographical "cell" (as the Java code was trying
to do on a grid).

References
----------
- Jouven, M., Carrère, P. and Baumont, R. (2006). 'Model predicting dynamics
  of biomass, structure and digestibility of herbage in managed permanent
  pastures. 1. Model description', Grass and Forage Science, vol. 61, no. 2,
  pp. 112-124. DOI: 10.1111/j.1365-2494.2006.00515.x.
"""

import itertools
import os
from datetime import datetime, timezone
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from climag.modvege import modvege
from climag.modvege_read_files import read_params, read_timeseries
from climag.plot_configs import ie_cordex_modvege_ncfile_name

output_vars = {
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
    "sea": ["Seasonal effect", "dimensionless"],
    "ftm": ["Temperature function", "dimensionless"],
    "env": ["Environmental limitation of growth", "dimensionless"],
    "pgr": ["Potential growth", "kg DM ha⁻¹"],
    "atr": ["Reproductive function", "dimensionless"]
}


def run_modvege_csv(input_timeseries_file, input_params_file, out_dir):
    """
    Input time series: CSV
    Also creates time series plots
    """

    # read parameter file into a dataframe
    params = read_params(filename=input_params_file)

    tseries, enddoy = read_timeseries(filename=input_timeseries_file)

    # initialise the run
    data_df = modvege(params=params, tseries=tseries, enddoy=enddoy)

    # convert output to dataframe and save as CSV
    data_df = tuple([list(range(1, len(data_df[0]) + 1))]) + data_df

    data_df = pd.DataFrame(
        zip(*data_df), columns=(["doy"] + list(output_vars.keys()))
    )

    data_df.to_csv(os.path.join(out_dir, "output.csv"), index=False)

    # plot all columns
    data_df.set_index("doy", inplace=True)

    plot_title = []
    for val in output_vars.values():
        val = " [".join(val) + "]"
        plot_title.append(val)

    data_df.plot(
        subplots=True, layout=(6, 3), figsize=(15, 14),
        xlabel="Day of the year", title=plot_title, legend=False
    )

    plt.tight_layout()

    plt.show()


def run_modvege_nc(input_timeseries_file, input_params_file, out_dir):
    """
    Input time series: NetCDF (climate data)
    """

    # read parameter file into a dataframe
    params = read_params(filename=input_params_file)

    tseries = xr.open_dataset(
        input_timeseries_file,
        decode_coords="all",
        # chunks="auto"  # the following operations do not work with Dask
        #                # yet, so chunking must be disabled...
    )

    # get the CRS
    data_crs = tseries.rio.crs

    # list of input variables
    input_vars = sorted(list(tseries.data_vars))

    # use rsds as par for now for the HiResIreland data
    if "par" not in tseries.data_vars:
        tseries = tseries.rename({"rsds": "par"})

    # loop through each year
    # for year in set(tseries_loc["time"].dt.year.values):
    for year in [2055]:
        tseries_y = tseries.sel(
            time=slice(f"{year}-01-01", f"{year}-12-31")
        )

        # extract the end day of the year
        # enddoy = tseries_y["time"].dt.dayofyear.values.max()

        # assign the outputs as new variables
        for key, val in output_vars.items():
            tseries_y[key] = xr.full_like(
                tseries_y["pr"], fill_value=np.nan
            )
            tseries_y[key].attrs = {
                # "standard_name": val[0].lower().replace(
                #     " ", "_"
                # ),
                "long_name": val[0],
                "units": val[1]
            }

        # create a dictionary to store the time series output
        # dataframes
        data_df = {}

        # loop through each grid cell
        # for rlon, rlat in [(20, 20), (21, 21)]:
        for rlon, rlat in itertools.product(
            range(len(tseries.coords["rlon"])),
            range(len(tseries.coords["rlat"]))
        ):
            tseries_l = tseries_y.isel(rlon=rlon, rlat=rlat)

            # ignore null cells
            if not tseries_l["pr"].isnull().all():
                data_df[f"{rlon}_{rlat}_{year}"] = pd.DataFrame(
                    {"time": tseries_l["time"]}
                )  # create a dataframe using the time array

                # assign the variables to columns
                for var in tseries_l.data_vars:
                    data_df[f"{rlon}_{rlat}_{year}"][var] = tseries_l[var]

                # assign other variables
                data_df[f"{rlon}_{rlat}_{year}"]["eta"] = 0.0
                data_df[f"{rlon}_{rlat}_{year}"]["lai"] = 0.0
                data_df[f"{rlon}_{rlat}_{year}"]["gcut_height"] = 0.0
                # data_df[f"{rlon}_{rlat}_{year}"][
                #     "grazing_animal_count"
                # ] = 0.0
                # data_df[f"{rlon}_{rlat}_{year}"][
                #     "grazing_avg_animal_weight"
                # ] = 0.0

                # initialise the run
                data_df[f"{rlon}_{rlat}_{year}"] = modvege(
                    params=params,
                    tseries=data_df[f"{rlon}_{rlat}_{year}"],
                    enddoy=tseries_y["time"].dt.dayofyear.values.max()
                )

                # convert output to dataframe
                data_df[f"{rlon}_{rlat}_{year}"] = tuple([list(
                    range(1, len(data_df[f"{rlon}_{rlat}_{year}"][0]) + 1)
                )]) + data_df[f"{rlon}_{rlat}_{year}"]

                data_df[f"{rlon}_{rlat}_{year}"] = pd.DataFrame(
                    zip(*data_df[f"{rlon}_{rlat}_{year}"]),
                    columns=(["doy"] + list(output_vars.keys()))
                )

                # assign the output variables to the main xarray dataset
                for key in output_vars:
                    tseries_y[key].loc[dict(
                        rlon=tseries_l.coords["rlon"],
                        rlat=tseries_l.coords["rlat"],
                        time=tseries_l.coords["time"]
                    )] = np.array(data_df[f"{rlon}_{rlat}_{year}"][key])

        # delete input variables
        tseries_y = tseries_y.drop_vars(list(tseries.data_vars))

        # assign attributes for the data
        tseries_y.attrs = {
            "creation_date": str(datetime.now(tz=timezone.utc)),
            "contact": "nstreethran@ucc.ie",
            "frequency": "day",
            "references": "https://github.com/ClimAg",
            "input_data": str(tseries.attrs),
            "input_variables": input_vars
        }

        # save as a NetCDF file
        os.makedirs(out_dir, exist_ok=True)

        # reassign CRS
        tseries_y.rio.write_crs(data_crs, inplace=True)

        if tseries.attrs["contact"] == "rossby.cordex@smhi.se":
            file_name = os.path.join(
                out_dir,
                ie_cordex_modvege_ncfile_name(
                    cordex_data=tseries, output_data=tseries_y
                )
            )
        else:
            file_name = os.path.join(
                out_dir,
                f"modvege_{tseries.attrs['title']}_{year}.nc"
            )

        tseries_y.to_netcdf(file_name)


def run_modvege(input_params_file, input_timeseries_file, out_dir):
    """
    Preprocess the inputs to run ModVege as a function and save the results
    as a CSV file

    Outputs
    -------
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
    input_timeseries_file : File path for the input time series
    out_dir : Directory to store output file(s)
    """

    if input_timeseries_file.endswith(".csv"):
        run_modvege_csv(input_timeseries_file, input_params_file, out_dir)
    else:  # open the climate model dataset
        run_modvege_nc(input_timeseries_file, input_params_file, out_dir)
