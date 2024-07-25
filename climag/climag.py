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
ie = gpd.read_file(
    os.path.join("data", "boundaries", "boundaries_all.gpkg"),
    layer="NUTS_RG_01M_2021_2157_IE",
)
# mask for offshore areas
ie_bbox = gpd.read_file(
    os.path.join("data", "boundaries", "boundaries_all.gpkg"),
    layer="ne_10m_land_2157_IE_BBOX_DIFF",
)
# mask for non-pasture areas
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
    """Rotated pole transform for plotting EURO-CORDEX data.

    Parameters
    ----------
    data : xarray.Dataset
        Input EURO-CORDEX data

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


def weighted_average(data, time_groupby: str):
    """
    Calculate the weighted average

    - https://docs.xarray.dev/en/stable/user-guide/computation.html
    - https://ncar.github.io/esds/posts/2021/yearly-averages-xarray/
    - https://docs.xarray.dev/en/stable/examples/monthly-means.html
    - https://ncar.github.io/esds/posts/2020/Time/
    """
    # if months:
    #     data = data.sel(time=data["time"].dt.month.isin(months))
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


def weighted_average_year(data):
    ds_yearmon = {}
    for year in set(data.time.dt.year.values):
        ds_yearmon[year] = ds.sel(time=slice(str(year), str(year)))
        ds_yearmon[year] = weighted_average(
            data=ds_yearmon[year], time_groupby="month"
        )
        ds_mon = {}
        for month in ds_yearmon[year].month.values:
            ds_mon[month] = (
                ds_yearmon[year]
                .sel(month=month)
                .assign_coords(time=datetime.datetime(year, month, 15))
                .expand_dims(dim="time")
            )
        ds_yearmon[year] = xr.combine_by_coords(ds_mon.values())
    ds_yearmon = xr.combine_by_coords(ds_yearmon.values())
    return ds_yearmon


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
            "i_bm",
            "h_bm",
            "pgro",
            # "bm",
            # "gro",
        ]
    )
    return data


def load_all_data(clim_dataset, hist_only=False):
    """Load all results datasets

    Parameters
    ----------

    clim_dataset : str
        EURO-CORDEX or HiResIreland
    """
    ds = {}

    if hist_only:
        model_exp_list = product(model_list, ["historical"])
    else:
        model_exp_list = product(model_list, exp_list)

    for model, exp in model_exp_list:
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
        if model == "HadGEM2-ES":
            ds[f"{model}_{exp}"] = ds[f"{model}_{exp}"].convert_calendar(
                "360_day", align_on="year"
            )

        # assign new coordinates and dimensions
        ds[f"{model}_{exp}"] = (
            ds[f"{model}_{exp}"].assign_coords(exp=exp).expand_dims(dim="exp")
        )
        ds[f"{model}_{exp}"] = (
            ds[f"{model}_{exp}"]
            .assign_coords(model=model)
            .expand_dims(dim="model")
        )

        # calculate ingested + harvested biomass
        ds[f"{model}_{exp}"] = ds[f"{model}_{exp}"].assign(
            bm_c=(ds[f"{model}_{exp}"]["i_bm"] + ds[f"{model}_{exp}"]["h_bm"])
        )
        # # total biomass
        # ds[f"{model}_{exp}"] = ds[f"{model}_{exp}"].assign(
        #     bm_t=(ds[f"{model}_{exp}"]["bm"] + ds[f"{model}_{exp}"]["bm_c"])
        # )

        # drop unnecessary variables
        ds[f"{model}_{exp}"] = keep_minimal_vars(data=ds[f"{model}_{exp}"])
        # # drop time_bnds
        # ds[f"{model}_{exp}"] = ds[f"{model}_{exp}"].drop_vars(["time_bnds"])
        # if historical data only, keep overlapping years with obs data
        if hist_only:
            ds[f"{model}_{exp}"] = ds[f"{model}_{exp}"].sel(
                time=slice("1981", "2005")
            )
        # reassign CRS
        ds[f"{model}_{exp}"] = ds[f"{model}_{exp}"].rio.write_crs(
            crs_ds, inplace=True
        )

    return ds


def load_obs_data():
    mera = xr.open_mfdataset(
        glob.glob(
            os.path.join("data", "ModVege", "MERA", "*MERA_FC3hr_3_day*.nc")
        ),
        chunks="auto",
        decode_coords="all",
    )
    # remove spin-up year and years outside historical reference period
    mera = mera.sel(time=slice("1981", "2005"))
    # calculate ingested + harvested biomass
    mera = mera.assign(bm_c=(mera["i_bm"] + mera["h_bm"]))
    # drop unnecessary variables
    mera = keep_minimal_vars(data=mera)
    # reassign CRS
    mera = mera.rio.write_crs(projection_lambert_conformal, inplace=True)
    return mera


def calc_annual_mean(data_dict, seasonal, skipna, hist_only=False):
    ds_mean_ann = {}
    if hist_only:
        model_exp_list = product(model_list, ["historical"])
    else:
        model_exp_list = product(model_list, exp_list)
    for model, exp in model_exp_list:
        if seasonal:
            ds_seas_year = {}
            for year in set(data_dict[f"{model}_{exp}"].time.dt.year.values):
                ds_seas_year[year] = data_dict[f"{model}_{exp}"].sel(
                    time=slice(str(year), str(year))
                )
                ds_seas_year[year] = (
                    ds_seas_year[year]
                    .groupby("time.season")
                    .mean(dim="time", skipna=skipna)
                )
                ds_seas_year[year] = (
                    ds_seas_year[year]
                    .assign_coords(year=year)
                    .expand_dims(dim="year")
                )
            ds_mean_ann[f"{model}_{exp}"] = xr.combine_by_coords(
                ds_seas_year.values(), combine_attrs="override"
            )
        else:
            ds_mean_ann[f"{model}_{exp}"] = (
                data_dict[f"{model}_{exp}"]
                .groupby("time.year")
                .mean(dim="time", skipna=skipna)
            )
    # combine data
    ds_mean_ann = xr.combine_by_coords(
        ds_mean_ann.values(), combine_attrs="override"
    )
    # sort seasons in the right order
    if seasonal:
        ds_mean_ann = ds_mean_ann.reindex(season=season_list)
    return ds_mean_ann


def calc_obs_annual_mean(obs_data, seasonal, skipna):
    if seasonal:
        ds_seas_year = {}
        for year in set(obs_data.time.dt.year.values):
            ds_seas_year[year] = obs_data.sel(time=slice(str(year), str(year)))
            ds_seas_year[year] = (
                ds_seas_year[year]
                .groupby("time.season")
                .mean(dim="time", skipna=skipna)
            )
            ds_seas_year[year] = (
                ds_seas_year[year]
                .assign_coords(year=year)
                .expand_dims(dim="year")
            )
        obs_data_mean = xr.combine_by_coords(
            ds_seas_year.values(), combine_attrs="override"
        )
    else:
        obs_data_mean = obs_data.groupby("time.year").mean(
            dim="time", skipna=skipna
        )
    return obs_data_mean


def regrid_climate_model_data(
    obs_data, data_dict, seasonal=False, skipna=None
):
    # annual mean for climate data
    ds_calc = calc_annual_mean(
        data_dict=data_dict, seasonal=seasonal, skipna=skipna, hist_only=True
    )
    # annual mean for observational data
    obs_calc = calc_obs_annual_mean(
        obs_data=obs_data, seasonal=seasonal, skipna=skipna
    )
    clim_data = {}
    for model in model_list:
        # remove model and exp dimensions
        clim_data[model] = ds_calc.sel(model=model, exp="historical")
        clim_data[model] = clim_data[model].drop(["lat", "lon"])
        clim_data[model] = clim_data[model].rename({"rlon": "x", "rlat": "y"})
        clim_data[model] = clim_data[model].rio.reproject_match(
            obs_calc, resampling=rio.enums.Resampling.bilinear
        )
        clim_data[model] = clim_data[model].assign_coords(
            {"x": obs_calc["x"], "y": obs_calc["y"]}
        )
        # reassign coordinates and dimensions
        clim_data[model] = (
            clim_data[model]
            .assign_coords(model=model)
            .expand_dims(dim="model")
        )
        # reassign projection
        clim_data[model].rio.write_crs(
            projection_lambert_conformal, inplace=True
        )
    clim_data = xr.combine_by_coords(
        clim_data.values(), combine_attrs="override"
    )
    return clim_data, obs_calc


def calc_bias(obs_data, clim_model_data):
    bias_abs = clim_model_data - obs_data
    bias_rel = bias_abs / obs_data * 100
    return bias_abs, bias_rel


def calc_normalised_std(data_dict, seasonal=False, skipna=None):
    ds_calc = calc_annual_mean(
        data_dict=data_dict, seasonal=seasonal, skipna=skipna
    )
    # historical mean
    hist_mean = (
        ds_calc.sel(exp="historical")
        .drop_vars("exp")
        .mean(dim="year", skipna=skipna)
    )
    # historical standard deviation
    hist_std = (
        ds_calc.sel(exp="historical")
        .drop_vars("exp")
        .std(dim="year", skipna=skipna, ddof=1)
    )
    # normalise
    ds_norm = (ds_calc - hist_mean) / hist_std
    return ds_norm


def calc_normalised_relative(data_dict, seasonal=False, skipna=None):
    ds_calc = calc_annual_mean(
        data_dict=data_dict, seasonal=seasonal, skipna=skipna
    )
    # historical mean
    hist_mean = (
        ds_calc.sel(exp="historical")
        .drop_vars("exp")
        .mean(dim="year", skipna=skipna)
    )
    # normalise
    ds_norm = (ds_calc - hist_mean) / hist_mean * 100
    return ds_norm


def calc_anomaly_absolute(data_dict, seasonal=False, skipna=None):
    ds_calc = calc_annual_mean(
        data_dict=data_dict, seasonal=seasonal, skipna=skipna
    )
    # historical mean
    hist_mean = (
        ds_calc.sel(exp="historical")
        .drop_vars("exp")
        .mean(dim="year", skipna=skipna)
    )
    # calculate anomaly
    ds_anom = ds_calc - hist_mean
    return ds_anom


def calc_event_frequency_intensity(data_dict, seasonal=False, skipna=None):
    ds_calc = calc_annual_mean(
        data_dict=data_dict, seasonal=seasonal, skipna=skipna
    )
    # historical 10th percentile
    hist_p10 = (
        ds_calc.sel(exp="historical")
        .drop_vars("exp")
        .chunk(dict(year=-1))
        .quantile(dim="year", skipna=skipna, q=0.1)
    )
    # historical standard deviation
    hist_std = (
        ds_calc.sel(exp="historical")
        .drop_vars("exp")
        .std(dim="year", skipna=skipna, ddof=1)
    )
    # calculate difference
    ds_anom = ds_calc - hist_p10
    # set negative values to 1, positive to 0
    ds_freq = xr.where(ds_anom < 0, 1, 0)
    # sum for all years
    ds_freq = ds_freq.sum(dim="year", skipna=skipna)
    # historical frequency
    hist_freq = ds_freq.sel(exp="historical").drop_vars("exp")
    # number of times more frequent in future
    ds_freq = ds_freq / hist_freq
    # intensity - keep only negative values
    ds_int = xr.where(ds_anom < 0, ds_anom, 0)
    ds_int = -ds_int / hist_std
    return ds_anom, ds_freq, ds_int


def plot_stats(dataset, transform, levels=14, seasonal=False, cmap="BrBG"):
    if seasonal:
        row = "season"
        figsize = (9, 16.25)
        pad = 0.015
    else:
        row = None
        figsize = (9, 4.75)
        pad = 0.075
    for v in list(dataset.data_vars):
        fig = dataset[v].plot.contourf(
            x="rlon",
            y="rlat",
            col="exp",
            row=row,
            subplot_kws={"projection": projection_hiresireland},
            transform=transform,
            xlim=(-1.775, 1.6),
            ylim=(-2.1, 2.1),
            cmap=cmap,
            robust=True,
            cbar_kwargs={"location": "bottom", "aspect": 30, "pad": pad},
            figsize=figsize,
            levels=levels,
        )
        for axis in fig.axs.flat:
            mask.to_crs(projection_hiresireland).plot(
                ax=axis, color="white", linewidth=0
            )
            ie_bbox.to_crs(projection_hiresireland).plot(
                ax=axis,
                edgecolor="darkslategrey",
                color="white",
                linewidth=0.5,
            )
        fig.set_titles("{value}", weight="semibold", fontsize=14)
        plt.show()
