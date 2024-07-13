"""Utility functions

"""

import glob
import os
import warnings
from datetime import datetime
from itertools import product

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio as rio
import seaborn as sns
import xarray as xr
from matplotlib import patheffects

# from dateutil.parser import parse

warnings.filterwarnings(
    action="ignore", category=RuntimeWarning, module="dask"
)

# Irish Transverse Mercator
ITM_EPSG = 2157

# set plot projection to the projection of the HiResIreland dataset
projection_hiresireland = ccrs.RotatedPole(
    pole_longitude=172.100006103516, pole_latitude=36.5999984741211
)

# define Lambert Conformal Conic projection for plots and transformations
# using MÉRA GRIB metadata
projection_lambert_conformal = ccrs.LambertConformal(
    false_easting=1481641.6769636814,
    false_northing=537326.0638850163,
    standard_parallels=[53.5],
    central_longitude=5.0,
    central_latitude=53.5,
)

projection_eurocordex = ccrs.RotatedPole(
    pole_longitude=-162.0, pole_latitude=39.25
)

# Ireland boundary
ie_bbox = gpd.read_file(
    os.path.join("data", "boundaries", "boundaries_all.gpkg"),
    layer="ne_10m_land_2157_IE_BBOX_DIFF",
)
ie = gpd.read_file(
    os.path.join("data", "boundaries", "boundaries_all.gpkg"),
    layer="NUTS_RG_01M_2021_2157_IE",
)

# mask out non-pasture areas
mask = gpd.read_file(
    os.path.join("data", "boundaries", "boundaries_all.gpkg"),
    layer="CLC_2018_MASK_PASTURE_2157_IE",
)

season_list = ["DJF", "MAM", "JJA", "SON"]
exp_list = ["historical", "rcp45", "rcp85"]
model_list = ["CNRM-CM5", "EC-EARTH", "HadGEM2-ES", "MPI-ESM-LR"]


def rotated_pole_point(data, lon, lat):
    """Convert lat/lon of a specific point to rotated pole coordinates

    Parameters
    ----------
    data : xarray.Dataset
        input climate data which uses rotated pole coordinates
    lon : float
        longitude of the point
    lat : float
        latitude of the point

    Returns
    -------
    tuple[float, float]
        transformed longitude and latitude in rotated pole coordinates
    """
    if data.rio.crs is None:
        pole_longitude = data["rotated_pole"].attrs[
            "grid_north_pole_longitude"
        ]
        pole_latitude = data["rotated_pole"].attrs["grid_north_pole_latitude"]
    else:
        pole_longitude = data.rio.crs.to_dict(projjson=True)["conversion"][
            "parameters"
        ][1]["value"]
        pole_latitude = data.rio.crs.to_dict(projjson=True)["conversion"][
            "parameters"
        ][0]["value"]
    rp_cds = ccrs.RotatedGeodetic(
        pole_longitude=pole_longitude,
        pole_latitude=pole_latitude,
    ).transform_point(x=lon, y=lat, src_crs=ccrs.Geodetic())
    return rp_cds[0], rp_cds[1]


def rotated_pole_transform(data):
    """Rotated pole transform for plotting CORDEX data.

    Parameters
    ----------
    data : xarray.Dataset
        Input CORDEX data

    Returns
    -------
    cartopy.crs.RotatedPole
        rotated pole transform
    """
    if data.rio.crs is None:
        pole_longitude = data["rotated_pole"].attrs[
            "grid_north_pole_longitude"
        ]
        pole_latitude = data["rotated_pole"].attrs["grid_north_pole_latitude"]
    else:
        pole_longitude = data.rio.crs.to_dict(projjson=True)["conversion"][
            "parameters"
        ][1]["value"]
        pole_latitude = data.rio.crs.to_dict(projjson=True)["conversion"][
            "parameters"
        ][0]["value"]
    transform = ccrs.RotatedPole(
        pole_longitude=pole_longitude, pole_latitude=pole_latitude
    )
    return transform


def weighted_average(data, time_groupby: str, months=None):
    """
    Calculate the weighted average

    - https://docs.xarray.dev/en/stable/user-guide/computation.html
    - https://ncar.github.io/esds/posts/2021/yearly-averages-xarray/
    - https://docs.xarray.dev/en/stable/examples/monthly-means.html
    - https://ncar.github.io/esds/posts/2020/Time/
    """
    if months:
        data = data.sel(time=data["time"].dt.month.isin(months))
    # calculate the weights by grouping month length by season or month
    weights = (
        data.time.dt.days_in_month.groupby(f"time.{time_groupby}")
        / data.time.dt.days_in_month.groupby(f"time.{time_groupby}").sum()
    )
    # test that the sum of weights for each year/season/month is one
    np.testing.assert_allclose(
        weights.groupby(f"time.{time_groupby}").sum().values,
        np.ones(len(set(weights[time_groupby].values))),
    )
    # calculate the weighted average
    data = (data * weights).groupby(f"time.{time_groupby}").sum(dim="time")
    return data


def stats(stat, data, time_groupby: str, months=None):
    if months:
        data = data.sel(time=data["time"].dt.month.isin(months))
    if stat == "std":
        data = data.groupby(f"time.{time_groupby}").std(dim="time", ddof=1)
    elif stat == "median":
        data = data.groupby(f"time.{time_groupby}").median(dim="time")
    elif stat == "0.9q":
        data = data.groupby(f"time.{time_groupby}").quantile(0.9, dim="time")
    elif stat == "0.1q":
        data = data.groupby(f"time.{time_groupby}").quantile(0.1, dim="time")
    return data


def keep_minimal_vars(data):
    """
    Drop variables that are not needed
    """
    data = data.drop_vars(
        [
            "bm_gv",
            "bm_gr",
            "bm_dv",
            "bm_dr",
            "age_gv",
            "age_gr",
            "age_dv",
            "age_dr",
            "omd_gv",
            "omd_gr",
            "lai",
            "env",
            "wr",
            "aet",
            "sen_gv",
            "sen_gr",
            "abs_dv",
            "abs_dr",
            "c_bm",
            # "i_bm",
            # "h_bm",
            # "pgro",
            # "bm",
            # "gro",
        ]
    )
    return data


def combine_datasets(dataset_dict, dataset_crs):
    """combine a dict of datasets into a single dataset"""
    dataset = xr.combine_by_coords(
        dataset_dict.values(), combine_attrs="override"
    )
    # dataset = xr.combine_nested(dataset_dict.values(), concat_dim=["model", "exp", "rlat", "rlon", "time"])
    dataset.rio.write_crs(dataset_crs, inplace=True)
    return dataset


def datasets_dict(clim_dataset):
    ds = {}
    # ds_res = {}
    # ds_mam = {}
    # ds_jja = {}
    # ds_son = {}

    for exp, model in product(exp_list, model_list):
        # auto-rechunking may cause NotImplementedError with object dtype
        # where it will not be able to estimate the size in bytes of object
        # data
        if model == "HadGEM2-ES":
            CHUNKS = 300
        else:
            CHUNKS = "auto"

        ds[f"{model}_{exp}"] = xr.open_mfdataset(
            glob.glob(
                os.path.join(
                    "data",
                    "ModVege",
                    clim_dataset,
                    exp,
                    model,
                    f"*{clim_dataset}*{model}*{exp}*.nc",
                )
            ),
            chunks=CHUNKS,
            decode_coords="all",
            # decode_cf=False,
        )

        # copy CRS
        crs_ds = ds[f"{model}_{exp}"].rio.crs

        # remove spin-up year
        if exp == "historical":
            ds[f"{model}_{exp}"] = ds[f"{model}_{exp}"].sel(
                time=slice("1976", "2005")
            )
        else:
            ds[f"{model}_{exp}"] = ds[f"{model}_{exp}"].sel(
                time=slice("2041", "2070")
            )

        # convert HadGEM2-ES data back to 360-day calendar
        # this ensures that the correct weighting is applied when
        # calculating the weighted average
        # if model == "HadGEM2-ES":
        #     ds[f"{model}_{exp}"] = ds[f"{model}_{exp}"].convert_calendar(
        #         "360_day", align_on="year"
        #     )

        # assign new coordinates and dimensions
        ds[f"{model}_{exp}"] = ds[f"{model}_{exp}"].assign_coords(exp=exp)
        ds[f"{model}_{exp}"] = ds[f"{model}_{exp}"].expand_dims(dim="exp")
        ds[f"{model}_{exp}"] = ds[f"{model}_{exp}"].assign_coords(model=model)
        ds[f"{model}_{exp}"] = ds[f"{model}_{exp}"].expand_dims(dim="model")

        # calculate cumulative biomass
        ds[f"{model}_{exp}"] = ds[f"{model}_{exp}"].assign(
            bm_t=(
                ds[f"{model}_{exp}"]["bm"]
                + ds[f"{model}_{exp}"]["i_bm"]
                + ds[f"{model}_{exp}"]["h_bm"]
            )
        )

        # drop unnecessary variables
        ds[f"{model}_{exp}"] = keep_minimal_vars(data=ds[f"{model}_{exp}"])
        # drop time_bnds
        ds[f"{model}_{exp}"] = ds[f"{model}_{exp}"].drop_vars(["time_bnds"])
        # # weighted mean - yearly, MAM
        # ds_mam[f"{model}_{exp}"] = weighted_average(ds[f"{model}_{exp}"], "year", [3, 4, 5])
        # # weighted mean - yearly, JJA
        # ds_jja[f"{model}_{exp}"] = weighted_average(ds[f"{model}_{exp}"], "year", [6, 7, 8])
        # # weighted mean - yearly, SON
        # ds_son[f"{model}_{exp}"] = weighted_average(ds[f"{model}_{exp}"], "year", [9, 10, 11])
        # # weighted annual average
        # ds[f"{model}_{exp}"] = weighted_average(ds[f"{model}_{exp}"], "year")

    # combine data
    ds = xr.merge(ds.values())
    # reassign CRS
    ds.rio.write_crs(crs_ds, inplace=True)
    # ds = combine_datasets(ds, crs_ds)
    # ds_mam = combine_datasets(ds_mam, crs_ds)
    # ds_jja = combine_datasets(ds_jja, crs_ds)
    # ds_son = combine_datasets(ds_son, crs_ds)

    # ensemble mean
    # ds_son = ds_son.mean(dim="model", skipna=True)
    # ds_mam = ds_mam.mean(dim="model", skipna=True)
    # ds_jja = ds_jja.mean(dim="model", skipna=True)

    # long-term average
    # ds_son_lta = ds_son.mean(dim="year", skipna=True)

    # # shift SON year by one
    # ds_son["year"] = ds_son["year"] + 1

    return ds, crs_ds  #, ds_mam, ds_jja, ds_son


def calculate_stats(ds, crs_ds, clim_dataset, stat):
    ds_stats = {}
    land_mask = gpd.read_file(os.path.join("data", "ModVege", "params.gpkg"), layer=clim_dataset.replace("-", "").lower())
    ds_stats["year"] = {}
    # weighted mean
    ds_stats[season]["mean"][f"{model}_{exp}"] = weighted_average(ds[f"{model}_{exp}"], "season")


    for season in ["year", "MAM", "JJA", "SON"]:
        ds_stats[season] = {}



        for stat in ["mean", "std", "median"]:
            ds_stats[season][stat] = {}
        for exp, model in product(exp_list, model_list):
            if season == "year":
                time_groupby = "year"
            else:
                # weighted mean
                ds_stats[season]["mean"][f"{model}_{exp}"] = weighted_average(ds[f"{model}_{exp}"], "season")
                # standard deviation and quantiles
                for stat in ["std", "median"]:
                    ds_stats[season][stat][f"{model}_{exp}"] = stats(stat, ds[f"{model}_{exp}"], "season")
            # weighted mean
            ds_stats[season]["mean"][f"{model}_{exp}"] = weighted_average(ds[f"{model}_{exp}"], time_groupby)
            # standard deviation and quantiles
            for stat in ["std", "median"]:
                ds_stats[season][stat][f"{model}_{exp}"] = stats(stat, ds[f"{model}_{exp}"], "season")
        # for season, months in zip(["year", "MAM", "JJA", "SON"], [None, [3, 4, 5], [6, 7, 8], [9, 10, 11]]):
        #     # weighted mean
        #     ds_stats[season]["mean"][f"{model}_{exp}"] = weighted_average(ds[f"{model}_{exp}"], "year", months)
        #     # standard deviation and quantiles
        #     for stat in ["std", "median"]:
        #         ds_stats[season][stat][f"{model}_{exp}"] = stats(stat, ds[f"{model}_{exp}"], "year", months)
    # combine data
    for season, stat in product(["year", "MAM", "JJA", "SON"], ["mean", "std", "median"]):
        ds_stats[season][stat] = combine_datasets(ds_stats[season][stat], crs_ds)
        # crop offshore cells
        ds_stats[season][stat] = ds_stats[season][stat].rio.clip(land_mask.to_crs(crs_ds).dissolve()["geometry"], all_touched=True)
    return ds_stats
