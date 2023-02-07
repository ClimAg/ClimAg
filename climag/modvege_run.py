"""modvege_run.py

https://github.com/YannChemin/modvege

Run in a Python interpreter
exec(open("climag/modvege_run.py").read())
"""

import itertools
import os
from datetime import datetime, timezone
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from climag.modvege import modvege
from climag.modvege_read_files import read_params, read_timeseries

output_vars = {
    "bm_gv": ["Green vegetative biomass", "kg DM ha⁻¹"],
    "bm_gr": ["Green reproductive biomass", "kg DM ha⁻¹"],
    "bm_dv": ["Dead vegetative biomass", "kg DM ha⁻¹"],
    "bm_dr": ["Dead reproductive biomass", "kg DM ha⁻¹"],
    "bm": ["Total standing biomass", "kg DM ha⁻¹"],
    "pgro": ["Potential growth", "kg DM ha⁻¹"],
    "gro": ["Total growth", "kg DM ha⁻¹"],
    "i_bm": ["Ingested biomass", "kg DM ha⁻¹"],
    "h_bm": ["Harvested biomass", "kg DM ha⁻¹"],
    "env": ["Environmental limitation of growth", "dimensionless"],
    "rep": ["Reproductive function", "dimensionless"],
    "lai": ["Leaf area index", "dimensionless"],
    "aet": ["Actual evapotranspiration", "mm"],
    "wr": ["Water reserves", "mm"]
}


def run_modvege_csv(input_timeseries_file, input_params_file, out_dir):
    """
    Input time series: CSV
    Also creates time series plots
    """

    # read parameter file into a dataframe
    params = read_params(filename=input_params_file)

    tseries, endday = read_timeseries(filename=input_timeseries_file)

    # initial water reserves
    params["wr"] = params["whc"]

    # initialise the run
    data_df = modvege(params=params, tseries=tseries, endday=endday)

    # convert output to dataframe
    data_df = pd.DataFrame(data_df)

    # drop biomass age columns
    data_df.drop(
        columns=["age_gv", "age_gr", "age_dv", "age_dr"],
        inplace=True
    )

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
        subplots=True, layout=(7, 3), figsize=(15, 14),
        xlabel="", title=plot_title, legend=False
    )
    plt.tight_layout()
    plt.show()


def run_modvege_nc(
    input_timeseries_file, input_params_file, out_dir,
    input_params_vector=None
):
    """
    Input time series: NetCDF (climate data)
    """

    params = {}
    model_vals = {}

    # read parameter file into a dataframe
    params["csv"] = read_params(filename=input_params_file)

    tseries = xr.open_dataset(
        input_timeseries_file,
        decode_coords="all",
        # chunks="auto"  # the following operations do not work with Dask
        #                # yet, so chunking must be disabled...
    )

    if input_params_vector is not None:
        if tseries.attrs["contact"] == "rossby.cordex@smhi.se":
            params["gpkg"] = gpd.read_file(
                os.path.join("data", "ModVege", "params.gpkg"),
                layer="eurocordex"
            )
        else:
            params["gpkg"] = gpd.read_file(
                os.path.join("data", "ModVege", "params.gpkg"),
                layer="hiresireland"
            )

    # get the CRS
    model_vals["data_crs"] = tseries.rio.crs

    # loop through each year
    model_vals["year_list"] = [2054, 2055, 2056]
    # model_vals["year_list"] = list(
    #     sorted(set(tseries["time"].dt.year.values))
    # )
    for year in model_vals["year_list"]:
        tseries_y = tseries.sel(
            time=slice(
                f"{year}-01-01",
                f"{year}-12-{int(tseries['time'][-1].dt.day)}"
            )
        )

        # assign the outputs as new variables
        for key, val in output_vars.items():
            tseries_y[key] = xr.full_like(
                tseries_y["PP"], fill_value=np.nan
            )
            tseries_y[key].attrs = {
                "long_name": val[0],
                "units": val[1]
            }

        # create a dictionary to store the time series output
        # dataframes
        data_df = {}

        # loop through each grid cell
        # for rlon, rlat in [(20, 20), (21, 21)]:
        # for rlon, rlat in itertools.product(range(10), range(10)):
        for rlon, rlat in itertools.product(
            range(len(tseries_y.coords["rlon"])),
            range(len(tseries_y.coords["rlat"]))
        ):
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
                        "T": tseries_l["T"]
                    }
                )

                # site-specific characteristics
                if input_params_vector is not None:
                    for key in ["sr", "ni", "whc"]:
                        params["csv"][key] = float(
                            params["gpkg"][
                                (
                                    params["gpkg"]["rlon"] == float(
                                        tseries_l["rlon"].values
                                    )
                                ) &
                                (
                                    params["gpkg"]["rlat"] == float(
                                        tseries_l["rlat"].values
                                    )
                                )
                            ][key]
                        )

                # starting values
                if year > model_vals["year_list"][0]:
                    (
                        params["csv"]["bm_gv"], params["csv"]["bm_gr"],
                        params["csv"]["bm_dv"], params["csv"]["bm_dr"],
                        params["csv"]["age_gv"], params["csv"]["age_gr"],
                        params["csv"]["age_dv"], params["csv"]["age_dr"],
                        params["csv"]["wr"]
                    ) = (
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["bm_gv"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["bm_gr"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["bm_dv"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["bm_dr"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["age_gv"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["age_gr"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["age_dv"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["age_dr"],
                        model_vals[f"{rlon}_{rlat}_{year - 1}"]["wr"]
                    )
                else:
                    params["csv"]["wr"] = params["csv"]["whc"]

                # initialise the run
                data_df[f"{rlon}_{rlat}_{year}"] = modvege(
                    params=params["csv"],
                    tseries=data_df[f"{rlon}_{rlat}_{year}"],
                    endday=len(tseries_l["time"])
                )

                # convert output to dataframe
                data_df[f"{rlon}_{rlat}_{year}"] = pd.DataFrame(
                    data_df[f"{rlon}_{rlat}_{year}"]
                )

                # save starting values for the next simulation year
                model_vals[f"{rlon}_{rlat}_{year}"] = {}
                (
                    model_vals[f"{rlon}_{rlat}_{year}"]["bm_gv"],
                    model_vals[f"{rlon}_{rlat}_{year}"]["bm_gr"],
                    model_vals[f"{rlon}_{rlat}_{year}"]["bm_dv"],
                    model_vals[f"{rlon}_{rlat}_{year}"]["bm_dr"],
                    model_vals[f"{rlon}_{rlat}_{year}"]["age_gv"],
                    model_vals[f"{rlon}_{rlat}_{year}"]["age_gr"],
                    model_vals[f"{rlon}_{rlat}_{year}"]["age_dv"],
                    model_vals[f"{rlon}_{rlat}_{year}"]["age_dr"],
                    model_vals[f"{rlon}_{rlat}_{year}"]["wr"]
                ) = (
                    data_df[f"{rlon}_{rlat}_{year}"]["bm_gv"].iat[-1],
                    data_df[f"{rlon}_{rlat}_{year}"]["bm_gr"].iat[-1],
                    data_df[f"{rlon}_{rlat}_{year}"]["bm_dv"].iat[-1],
                    data_df[f"{rlon}_{rlat}_{year}"]["bm_dr"].iat[-1],
                    data_df[f"{rlon}_{rlat}_{year}"]["age_gv"].iat[-1],
                    data_df[f"{rlon}_{rlat}_{year}"]["age_gr"].iat[-1],
                    data_df[f"{rlon}_{rlat}_{year}"]["age_dv"].iat[-1],
                    data_df[f"{rlon}_{rlat}_{year}"]["age_dr"].iat[-1],
                    data_df[f"{rlon}_{rlat}_{year}"]["wr"].iat[-1]
                )

                # drop biomass age columns
                data_df[f"{rlon}_{rlat}_{year}"].drop(
                    columns=["age_gv", "age_gr", "age_dv", "age_dr"],
                    inplace=True
                )

                # assign the output variables to the main Xarray dataset
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
            "references": "https://github.com/ClimAg"
        }

        # reassign CRS
        tseries_y.rio.write_crs(model_vals["data_crs"], inplace=True)

        # save as a NetCDF file
        if tseries.attrs["contact"] == "rossby.cordex@smhi.se":
            model_vals["out_dir"] = os.path.join(
                out_dir, "EURO-CORDEX",
                tseries.attrs["experiment_id"],
                tseries.attrs["driving_model_id"]
            )
            os.makedirs(model_vals["out_dir"], exist_ok=True)
            tseries_y.to_netcdf(os.path.join(
                model_vals["out_dir"],
                f"modvege_IE_EURO-CORDEX_{tseries.attrs['model_id']}_"
                f"{tseries.attrs['driving_model_id']}_"
                f"{tseries.attrs['experiment_id']}_"
                f"{tseries.attrs['CORDEX_domain']}_{year}.nc"
            ))
        else:
            model_vals["out_dir"] = os.path.join(
                out_dir, "HiResIreland",
                tseries.attrs["title"].split("_")[2],
                tseries.attrs["title"].split("_")[1]
            )
            os.makedirs(model_vals["out_dir"], exist_ok=True)
            tseries_y.to_netcdf(os.path.join(
                model_vals["out_dir"],
                f"modvege_IE_HiResIreland_{tseries.attrs['title']}_{year}.nc"
            ))


def run_modvege(
    input_params_file, input_timeseries_file, out_dir,
    input_params_vector=None
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
            input_timeseries_file, input_params_file, out_dir,
            input_params_vector
        )
