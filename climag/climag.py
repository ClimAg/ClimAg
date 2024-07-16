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
        data = data.groupby(f"time.{time_groupby}").std(dim="time", ddof=1, skipna=True)
    elif stat == "median":
        data = data.groupby(f"time.{time_groupby}").median(dim="time", skipna=True)
    elif stat == "0.9q":
        data = data.groupby(f"time.{time_groupby}").quantile(0.9, dim="time", skipna=True)
    elif stat == "0.1q":
        data = data.groupby(f"time.{time_groupby}").quantile(0.1, dim="time", skipna=True)
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


def load_modvege_dataset(clim_dataset, exp):
    ds = {}

    for model in model_list:
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

    return ds, crs_ds


def results_mean(clim_dataset):
    ds = {}
    for exp in exp_list:
        ds[exp], crs_ds = load_modvege_dataset(clim_dataset=clim_dataset, exp=exp)
    ds = xr.combine_by_coords(ds.values(), combine_attrs="override")
    ds.rio.write_crs(crs_ds, inplace=True)
    # seasonal and annual weighted means
    ds_season = weighted_average(ds, "season")
    ds_annual = weighted_average(ds, "year")
    # land mask geometry for cropping
    land_mask = gpd.read_file(os.path.join("data", "ModVege", "params.gpkg"), layer=clim_dataset.replace("-", "").lower())
    land_mask = land_mask.to_crs(crs_ds).dissolve()["geometry"]
    # calculate ensemble means
    # crop offshore cells
    for data in [ds_season, ds_annual]:
        data_e = data.mean(dim="model", skipna=True)
        data_e = data_e.assign_coords(model="Ensemble")
        data_e = data_e.expand_dims(dim="model")
        data = xr.combine_by_coords([data, data_e], combine_attrs="override")
        data.rio.write_crs(crs_ds, inplace=True)
        data = data.rio.clip(land_mask, all_touched=True)
    # sort seasons in the right order
    ds_season = ds_season.reindex(season=season_list)
    return ds_season, ds_year, crs_ds


def results_normalised(data_season, data_year):
    mean_hist = data_year.sel(exp="historical").drop_vars("exp").mean(dim="year", skipna=True)
    std_hist = data_year.sel(exp="historical").drop_vars("exp").std(dim="year", ddof=1, skipna=True)
    ds_season_norm = (ds_season - mean_hist) / std_hist
    ds_year_norm = (ds_year - mean_hist) / std_hist
    return ds_season_norm, ds_year_norm


def calculate_stats(ds, crs_ds, clim_dataset, stat):
    ds_stats = {}
    if stat == "mean":
        ds_stats["season"] = weighted_average(ds, "season")
        ds_stats["year"] = weighted_average(ds, "year")
        ds_stats["lt"] = weighted_average(ds, "year").mean(dim=["year", "model"], skipna=True)
        ds_stats["season_ensemble"] = weighted_average(ds, "season").mean(dim="model", skipna=True)
        ds_stats["year_ensemble"] = weighted_average(ds, "year").mean(dim="model", skipna=True)
        ds_stats["lt_ensemble"] = weighted_average(ds, "year").mean(dim=["year", "model"], skipna=True)
    elif stat == "std":
        ds_stats["season"] = ds.groupby("time.season").std(dim="time", ddof=1, skipna=True)
        ds_stats["year"] = ds.groupby("time.year").std(dim="time", ddof=1, skipna=True)
        ds_stats["lt"] = ds.std(dim="time", ddof=1, skipna=True)
        ds_stats["season_ensemble"] = ds.groupby("time.season").std(dim=["time", "model"], ddof=1, skipna=True)
        ds_stats["year_ensemble"] = ds.groupby("time.year").std(dim=["time", "model"], ddof=1, skipna=True)
        ds_stats["lt_ensemble"] = ds.std(dim=["time", "model"], ddof=1, skipna=True)
    # sort seasons in the right order
    for key in ds_stats:
        if "season" in key:
            ds_stats[key] = ds_stats[key].reindex(season=season_list)
    # land mask geometry for cropping
    land_mask = gpd.read_file(os.path.join("data", "ModVege", "params.gpkg"), layer=clim_dataset.replace("-", "").lower())
    land_mask = land_mask.to_crs(crs_ds).dissolve()["geometry"]
    # crop offshore cells
    for key in ds_stats:
        ds_stats[key].rio.write_crs(crs_ds, inplace=True)
        ds_stats[key] = ds_stats[key].rio.clip(land_mask, all_touched=True)
    return ds_stats
