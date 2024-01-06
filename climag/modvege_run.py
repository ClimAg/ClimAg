"""modvege_run.py

https://code.europa.eu/agri4cast/modvege
"""

import itertools
import os
from datetime import datetime, timezone

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from climag.modvege import modvege
from climag.modvege_read_files import read_params, read_timeseries

np.seterr("raise")

output_vars = {
    "bm_gv": ["Green vegetative biomass", "kg DM ha⁻¹"],
    "bm_gr": ["Green reproductive biomass", "kg DM ha⁻¹"],
    "bm_dv": ["Dead vegetative biomass", "kg DM ha⁻¹"],
    "bm_dr": ["Dead reproductive biomass", "kg DM ha⁻¹"],
    "age_gv": ["Green vegetative biomass age", "kg DM ha⁻¹"],
    "age_gr": ["Green reproductive biomass age", "kg DM ha⁻¹"],
    "age_dv": ["Dead vegetative biomass age", "kg DM ha⁻¹"],
    "age_dr": ["Dead reproductive biomass age", "kg DM ha⁻¹"],
    "bm": ["Standing biomass", "kg DM ha⁻¹"],
    "pgro": ["Potential growth", "kg DM ha⁻¹ day⁻¹"],
    "gro": ["Total growth", "kg DM ha⁻¹ day⁻¹"],
    "i_bm": ["Ingested biomass", "kg DM ha⁻¹"],
    "h_bm": ["Harvested biomass", "kg DM ha⁻¹"],
    "c_bm": ["Daily ingested biomass", "kg DM ha⁻¹ day⁻¹"],
    "env": ["Environmental limitation of growth", "dimensionless"],
    "lai": ["Leaf area index", "dimensionless"],
    "aet": ["Actual evapotranspiration", "mm day⁻¹"],
    "wr": ["Water reserves", "mm day⁻¹"],
    "sen_gv": ["Senescence of green vegetative biomass", "kg DM ha⁻¹"],
    "sen_gr": ["Senescence of green reproductive biomass", "kg DM ha⁻¹"],
    "abs_dv": ["Abscission of dead vegetative biomass", "kg DM ha⁻¹"],
    "abs_dr": ["Abscission of dead reproductive biomass", "kg DM ha⁻¹"],
    "omd_gv": [
        "Green vegetative biomass organic matter digestibility",
        "kg DM ha⁻¹",
    ],
    "omd_gr": [
        "Green reproductive biomass organic matter digestibility",
        "kg DM ha⁻¹",
    ],
}


def run_modvege_csv(input_timeseries_file, input_params_file, out_dir):
    """Input time series: CSV

    Also creates time series plots
    """
    # read parameter file into a dataframe
    params = read_params(filename=input_params_file)

    tseries, endday = read_timeseries(filename=input_timeseries_file)

    # assume the initial water reserves is equivalent to the water holding
    # capacity
    params["wr"] = params["whc"]

    # initialise the run
    data_df = modvege(params=params, tseries=tseries, endday=endday)

    # convert output to dataframe
    data_df = pd.DataFrame(data_df)

    # # drop biomass age columns
    # data_df.drop(
    #     columns=["age_gv", "age_gr", "age_dv", "age_dr"], inplace=True
    # )

    # save as CSV
    data_df.to_csv(os.path.join(out_dir, "output.csv"), index=False)

    # plot all columns
    data_df.set_index("time", inplace=True)

    # configure plot title
    plot_title = []
    for val in output_vars.values():
        val = " [".join(val) + "]"
        plot_title.append(val)

    # plot data
    data_df.plot(
        subplots=True,
        layout=(9, 3),
        figsize=(18, 18),
        xlabel="",
        title=plot_title,
        legend=False,
    )
    plt.tight_layout()
    plt.show()


def site_specific_params_file(input_params_vector, tseries, params):
    """
    Load the site-specific characteristics layers that vary spatially
    """
    # site-specific characteristics that vary spatially
    if input_params_vector is not None:
        if "EURO-CORDEX" in tseries.attrs["dataset"]:
            params["gpkg"] = gpd.read_file(
                input_params_vector, layer="eurocordex"
            )
        elif "HiResIreland" in tseries.attrs["dataset"]:
            params["gpkg"] = gpd.read_file(
                input_params_vector, layer="hiresireland"
            )
        else:
            params["gpkg"] = gpd.read_file(input_params_vector, layer="mera")


def run_modvege_nc(
    input_timeseries_file, input_params_file, out_dir, input_params_vector=None
):
    """
    Input time series: NetCDF (climate data)
    """
    print(
        f"Running simulations for input file '{input_timeseries_file}'...",
        datetime.now(tz=timezone.utc),
    )

    # dictionary to store parameters
    params = {}

    # dictionary to store intermediate time series values
    model_vals = {}

    # read parameter file into a dataframe
    params["csv"] = read_params(filename=input_params_file)

    # use x and y coordinates depending on the dataset
    if "MERA" in input_timeseries_file:
        xcoord, ycoord = "x", "y"
    else:
        xcoord, ycoord = "rlon", "rlat"

    tseries = xr.open_dataset(
        input_timeseries_file,
        decode_coords="all",
        # chunks="auto"  # the following operations do not work with Dask
        #                # yet, so chunking must be disabled...
    )

    # site-specific characteristics that vary spatially
    site_specific_params_file(
        input_params_vector=input_params_vector, tseries=tseries, params=params
    )

    # get the CRS
    model_vals["data_crs"] = tseries.rio.crs

    # loop through each year
    model_vals["year_list"] = list(sorted(set(tseries["time"].dt.year.values)))

    for year in model_vals["year_list"]:
        tseries_y = tseries.sel(time=slice(f"{year}-01-01", f"{year}-12-31"))

        # assign the outputs as new variables
        for key, val in output_vars.items():
            tseries_y[key] = xr.full_like(tseries_y["PP"], fill_value=np.nan)
            tseries_y[key].attrs = {"long_name": val[0], "units": val[1]}

        # create a dictionary to store the time series output
        # dataframes
        data_df = {}

        # loop through each grid cell
        for rlon, rlat in itertools.product(
            range(len(tseries_y.coords[xcoord])),
            range(len(tseries_y.coords[ycoord])),
        ):
            if "MERA" in input_timeseries_file:
                tseries_l = tseries_y.isel(x=rlon, y=rlat)
            else:
                tseries_l = tseries_y.isel(rlon=rlon, rlat=rlat)

            # ignore null cells
            if not tseries_l["PP"].isnull().all():
                # create a dataframe using the time array and variables
                data_df[f"{rlon}_{rlat}_{year}"] = pd.DataFrame(
                    {
                        "time": tseries_l["time"],
                        "PP": tseries_l["PP"],
                        "PAR": tseries_l["PAR"],
                        "PET": tseries_l["PET"],
                        "T": tseries_l["T"],
                    }
                )

                # site-specific characteristics
                if input_params_vector is not None:
                    for key in ["sr", "ni", "whc"]:
                        params["csv"][key] = float(
                            params["gpkg"][
                                (
                                    params["gpkg"][xcoord]
                                    == float(tseries_l[xcoord].values)
                                )
                                & (
                                    params["gpkg"][ycoord]
                                    == float(tseries_l[ycoord].values)
                                )
                            ][key]
                        )

                # temperatures for calculating ten-day moving averages in the
                # following year when day < 10
                model_vals[f"{rlon}_{rlat}_{year}"] = {}
                model_vals[f"{rlon}_{rlat}_{year}"]["t_init"] = data_df[
                    f"{rlon}_{rlat}_{year}"
                ]["T"].iloc[-10:-1]

                # starting values
                if year > model_vals["year_list"][0]:
                    (
                        params["csv"]["bm_gv"],
                        params["csv"]["bm_gr"],
                        params["csv"]["bm_dv"],
                        params["csv"]["bm_dr"],
                        params["csv"]["age_gv"],
                        params["csv"]["age_gr"],
                        params["csv"]["age_dv"],
                        params["csv"]["age_dr"],
                        params["csv"]["wr"],
                        model_vals["t_init"],
                    ) = (
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["bm_gv"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["bm_gr"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["bm_dv"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["bm_dr"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["age_gv"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["age_gr"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["age_dv"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["age_dr"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["wr"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["t_init"],
                    )
                else:
                    params["csv"]["wr"] = params["csv"]["whc"]
                    model_vals["t_init"] = None

                # initialise the run
                data_df[f"{rlon}_{rlat}_{year}"] = modvege(
                    params=params["csv"],
                    tseries=data_df[f"{rlon}_{rlat}_{year}"],
                    endday=len(tseries_l["time"]),
                    t_init=model_vals["t_init"],
                )

                # convert output to dataframe
                data_df[f"{rlon}_{rlat}_{year}"] = pd.DataFrame(
                    data_df[f"{rlon}_{rlat}_{year}"]
                )

                # save starting values for the next simulation year
                # round to 7 decimal places to prevent FloatingPointError
                (
                    model_vals[f"{rlon}_{rlat}_{year}"]["bm_gv"],
                    model_vals[f"{rlon}_{rlat}_{year}"]["bm_gr"],
                    model_vals[f"{rlon}_{rlat}_{year}"]["bm_dv"],
                    model_vals[f"{rlon}_{rlat}_{year}"]["bm_dr"],
                    model_vals[f"{rlon}_{rlat}_{year}"]["age_gv"],
                    model_vals[f"{rlon}_{rlat}_{year}"]["age_gr"],
                    model_vals[f"{rlon}_{rlat}_{year}"]["age_dv"],
                    model_vals[f"{rlon}_{rlat}_{year}"]["age_dr"],
                    model_vals[f"{rlon}_{rlat}_{year}"]["wr"],
                ) = (
                    round(
                        data_df[f"{rlon}_{rlat}_{year}"]["bm_gv"].iat[-1], 7
                    ),
                    round(
                        data_df[f"{rlon}_{rlat}_{year}"]["bm_gr"].iat[-1], 7
                    ),
                    round(
                        data_df[f"{rlon}_{rlat}_{year}"]["bm_dv"].iat[-1], 7
                    ),
                    round(
                        data_df[f"{rlon}_{rlat}_{year}"]["bm_dr"].iat[-1], 7
                    ),
                    round(
                        data_df[f"{rlon}_{rlat}_{year}"]["age_gv"].iat[-1], 7
                    ),
                    round(
                        data_df[f"{rlon}_{rlat}_{year}"]["age_gr"].iat[-1], 7
                    ),
                    round(
                        data_df[f"{rlon}_{rlat}_{year}"]["age_dv"].iat[-1], 7
                    ),
                    round(
                        data_df[f"{rlon}_{rlat}_{year}"]["age_dr"].iat[-1], 7
                    ),
                    round(data_df[f"{rlon}_{rlat}_{year}"]["wr"].iat[-1], 7),
                )

                # # drop biomass age columns
                # data_df[f"{rlon}_{rlat}_{year}"].drop(
                #     columns=["age_gv", "age_gr", "age_dv", "age_dr"],
                #     inplace=True
                # )

                # assign the output variables to the main Xarray dataset
                for key in output_vars:
                    tseries_y[key].loc[
                        {
                            xcoord: tseries_l.coords[xcoord],
                            ycoord: tseries_l.coords[ycoord],
                            "time": tseries_l.coords["time"],
                        }
                    ] = np.array(data_df[f"{rlon}_{rlat}_{year}"][key])

        # delete input variables
        tseries_y = tseries_y.drop_vars(list(tseries.data_vars))

        # assign attributes for the data
        tseries_y.attrs = {
            "creation_date": str(datetime.now(tz=timezone.utc)),
            "contact": "nstreethran@ucc.ie",
            "frequency": "day",
            "references": "https://github.com/ClimAg",
            "input_dataset": tseries.attrs["dataset"],
        }

        # reassign CRS
        tseries_y.rio.write_crs(model_vals["data_crs"], inplace=True)

        # save as a NetCDF file
        if "MERA" in input_timeseries_file:
            model_vals["out_dir"] = os.path.join(out_dir, "MERA")
        else:
            model_vals["out_dir"] = os.path.join(
                out_dir,
                tseries.attrs["dataset"].split("_")[1],
                tseries.attrs["dataset"].split("_")[4],
                tseries.attrs["dataset"].split("_")[3],
            )
        os.makedirs(model_vals["out_dir"], exist_ok=True)
        tseries_y.to_netcdf(
            os.path.join(
                model_vals["out_dir"],
                f"modvege_{tseries.attrs['dataset']}_{year}.nc",
            )
        )

        print(f"{year} complete...", datetime.now(tz=timezone.utc))


def run_modvege(
    input_params_file, input_timeseries_file, out_dir, input_params_vector=None
):
    """
    Preprocess the inputs to run ModVege as a function and save the results
    as a CSV file

    Parameters
    ----------
    input_params_file : File path for the input parameters
    input_timeseries_file : File path for the input time series
    out_dir : Directory to store output file(s)
    """
    if input_timeseries_file.endswith(".csv"):
        run_modvege_csv(input_timeseries_file, input_params_file, out_dir)
    else:  # open the climate model dataset
        run_modvege_nc(
            input_timeseries_file,
            input_params_file,
            out_dir,
            input_params_vector,
        )
