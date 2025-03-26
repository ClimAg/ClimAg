"""Utility functions"""

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
from rioxarray.exceptions import MissingSpatialDimensionError

# from dateutil.parser import parse

warnings.filterwarnings("once")
warnings.simplefilter("once")
# warnings.filterwarnings(
#     action="once", category=RuntimeWarning, module="dask"
# )
# warnings.filterwarnings(
#     action="once", category=RuntimeWarning, module="numpy"
# )
# warnings.filterwarnings(
#     action="once", category=UserWarning, module="gribapi"
# )

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

# # Ireland boundary
# ie = gpd.read_file(
#     os.path.join("data", "boundaries", "boundaries_all.gpkg"),
#     layer="NUTS_RG_01M_2021_2157_IE",
# )
# # mask for offshore areas
# ie_bbox = gpd.read_file(
#     os.path.join("data", "boundaries", "boundaries_all.gpkg"),
#     layer="ne_10m_land_2157_IE_BBOX_DIFF",
# )
# # mask for non-pasture areas
# mask = gpd.read_file(
#     os.path.join("data", "boundaries", "boundaries_all.gpkg"),
#     layer="CLC_2018_MASK_PASTURE_2157_IE",
# )

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


def calc_annual_mean(data_dict, seasonal, skipna, var_avg, hist_only=False):
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
                ds_seas_year[year] = getattr(ds_seas_year[year].groupby("time.season"), var_avg)(dim="time", skipna=skipna)
                ds_seas_year[year] = (
                    ds_seas_year[year]
                    .assign_coords(year=year)
                    .expand_dims(dim="year")
                )
            ds_mean_ann[f"{model}_{exp}"] = xr.combine_by_coords(
                ds_seas_year.values(), combine_attrs="override"
            )
        else:
            ds_mean_ann[f"{model}_{exp}"] = getattr(data_dict[f"{model}_{exp}"].groupby("time.year"), var_avg)(dim="time", skipna=skipna)
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
    obs_data, data_dict, seasonal=False, skipna=None, var_avg="mean"
):
    # annual mean for climate data
    ds_calc = calc_annual_mean(
        data_dict=data_dict, seasonal=seasonal, skipna=skipna, var_avg=var_avg, hist_only=True
    )
    # annual mean for observational data
    obs_calc = calc_obs_annual_mean(
        obs_data=obs_data, seasonal=seasonal, skipna=skipna
    )
    clim_data = {}
    for model in model_list:
        # remove model and exp dimensions
        clim_data[model] = ds_calc.sel(model=model, exp="historical")
        clim_data[model] = clim_data[model].drop_vars(["lat", "lon"])
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


def calc_normalised_std(data_dict, seasonal=False, skipna=None, var_avg="mean"):
    ds_calc = calc_annual_mean(
        data_dict=data_dict, seasonal=seasonal, skipna=skipna, var_avg=var_avg
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


def calc_normalised_relative(data_dict, seasonal=False, skipna=None, var_avg="mean"):
    ds_calc = calc_annual_mean(
        data_dict=data_dict, seasonal=seasonal, skipna=skipna, var_avg=var_avg
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


def calc_anomaly_absolute(data_dict, seasonal=False, skipna=None, var_avg="mean"):
    ds_calc = calc_annual_mean(
        data_dict=data_dict, seasonal=seasonal, skipna=skipna, var_avg=var_avg
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


def calc_event_frequency_intensity(data_dict, variable, seasonal=False, skipna=None, var_avg="mean"):
    # https://docs.xarray.dev/en/stable/examples/visualization_gallery.html#Control-the-plot's-colorbar
    ds_calc = calc_annual_mean(
        data_dict=data_dict, seasonal=seasonal, skipna=skipna, var_avg=var_avg
    )[variable]
    # historical 10th percentile
    hist_p10 = (
        ds_calc.sel(exp="historical")
        .drop_vars("exp")
        .chunk(dict(year=-1))
        .quantile(dim="year", skipna=skipna, q=0.1)
    )
    hist_p10.rio.write_crs(ds_calc.rio.crs, inplace=True)
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
    ds_freq_norm = ds_freq / hist_freq
    # intensity - keep only negative values
    ds_int = xr.where(ds_anom < 0, -ds_anom, 0)
    ds_int.rio.write_crs(ds_calc.rio.crs, inplace=True)
    # intensity compared to historical std
    ds_int_std = ds_int / hist_std
    # intensity compared to historical p10
    ds_int_p10 = ds_int / hist_p10 * 100
    return ds_anom, ds_freq, ds_freq_norm, ds_int, ds_int_std, ds_int_p10


def describe_dataset(dataset, pastures, model=False, exp=False, season=False, xrdataset=True):
    try:
        dataset_df = dataset.rio.clip(pastures.to_crs(dataset.rio.crs), all_touched=True).to_dataframe()
    except MissingSpatialDimensionError:
        dataset_df = dataset.rename({"rlon": "x", "rlat": "y"}).rio.clip(pastures.to_crs(dataset.rio.crs), all_touched=True).to_dataframe()
    cols = []
    if exp:
        cols.append("exp")
    if model:
        cols.append("model")
    if season:
        cols.append("season")
    if xrdataset:
        vl = list(dataset.data_vars)
    else:
        vl = [dataset.name]
    cols += vl
    dataset_df = dataset_df.reset_index()[cols]
    dataset_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    dataset_df.dropna(subset=vl, how="all", inplace=True)
    if model:
        return dataset_df.groupby("model").describe().T
    elif exp:
        return dataset_df.groupby("exp").describe().T
    elif season:
        return dataset_df.groupby("season").describe().T
    else:
        return dataset_df.describe()


def plot_stats(dataset, transform, mask, ie_bbox, label, row=None, col="exp", levels=14, cmap="BrBG", extend="both"):
    if row == "season":
        figsize = (9, 16.25)
        pad = 0.015
    elif col == "season":
        figsize = (10, 8.5)
        pad = 0.045
    else:
        figsize = (9, 4.75)
        pad = 0.075
    # format exp names: https://github.com/pydata/xarray/discussions/9097
    dataset["exp"] = dataset["exp"].where(dataset["exp"] != "historical", "Historical")
    dataset["exp"] = dataset["exp"].where(dataset["exp"] != "rcp45", "RCP4.5")
    dataset["exp"] = dataset["exp"].where(dataset["exp"] != "rcp85", "RCP8.5")
    fig = dataset.plot.contourf(
        x="rlon",
        y="rlat",
        col=col,
        row=row,
        subplot_kws={"projection": projection_hiresireland},
        transform=transform,
        xlim=(-1.775, 1.6),
        ylim=(-2.1, 2.1),
        cmap=cmap,
        robust=True,
        extend=extend,
        cbar_kwargs={"location": "bottom", "aspect": 30, "pad": pad, "label": label},
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


def event_duration(data_dict, variable, skipna=None, var_avg="mean", rolling_time=7):
    data_model = {}
    data_count = {}
    data_val = {}
    for model, exp in product(model_list, exp_list):
        if model == "HadGEM2-ES":
            data_model[f"{model}_{exp}"] = data_dict[f"{model}_{exp}"][variable].convert_calendar("standard", align_on="year")
        else:
            data_model[f"{model}_{exp}"] = data_dict[f"{model}_{exp}"][variable]
        data_model[f"{model}_{exp}"] = data_model[f"{model}_{exp}"].rolling(time=rolling_time).mean(skipna=skipna).groupby(data_model[f"{model}_{exp}"]["time"].dt.isocalendar().week).mean(skipna=skipna)
        d = data_model[f"{model}_{exp}"] - data_model[f"{model}_historical"].sel(exp="historical").drop_vars("exp")
        data_count[f"{model}_{exp}"] = xr.where(d < 0, 1, 0).sum(dim="week", skipna=skipna)
        data_val[f"{model}_{exp}"] = getattr(xr.where(d < 0, -d, 0), var_avg)(dim="week", skipna=skipna)
    # combine data
    data_count = xr.combine_by_coords(
        data_count.values(), combine_attrs="override"
    )
    data_val = xr.combine_by_coords(
        data_val.values(), combine_attrs="override"
    )
    data_val.rio.write_crs(data_dict[f"{model}_historical"].rio.crs, inplace=True)
    return data_count, data_val


def calc_event_duration(data_dict, variable, seasonal=False, skipna=None, var_avg="mean"):
    # merge data of same driving model
    data_model = {}
    data_count = {}
    data_val = {}
    for model, exp in product(model_list, exp_list):
        # d = hist.sel(time=slice(year, year))
        # d = d.groupby(d["time"].dt.isocalendar().week).mean(skipna=True).assign_coords(year=int(year)).expand_dims(dim="year")
        data_model[f"{model}_{exp}"] = data_dict[f"{model}_{exp}"].rolling(time=30).mean(skipna=skipna)[variable]#.sel(exp="historical").drop_vars("exp")
        yr_dict = {}
        for year in data_model[f"{model}_{exp}"].time.dt.year.values:
            yr_dict[year] = data_model[f"{model}_{exp}"].sel(time=slice(str(year), str(year)))
            yr_dict[year] = yr_dict[year].groupby(yr_dict[year]["time"].dt.isocalendar().week).mean(skipna=True)
            if exp == "historical":
                yrval = year + 65
            else:
                yrval = year
            yr_dict[year] = yr_dict[year].assign_coords(year=int(yrval)).expand_dims(dim="year")
        data_model[f"{model}_{exp}"] = xr.combine_by_coords(yr_dict.values(), combine_attrs="override")
        d = data_model[f"{model}_{exp}"] - data_model[f"{model}_historical"].sel(exp="historical").drop_vars("exp")
        # hist = data_model[f"{model}_historical"].sel(exp="historical").drop_vars("exp")

        # d = data_dict[f"{model}_{exp}"].rolling(time=30).mean(skipna=skipna)[variable] - hist
        # rcp45 = data_dict[f"{model}_rcp45"].rolling(time=30).mean(skipna=skipna)[variable]
        # rcp85 = data_dict[f"{model}_rcp85"].rolling(time=30).mean(skipna=skipna)[variable]
        # diff45 = rcp45 - hist
        # diff85 = rcp85 - hist
        # keep only negative values and aggregate as yearly data
        data_count[f"{model}_{exp}"] = xr.where(d < 0, 1, 0).sum(dim="week", skipna=skipna)
        data_val[f"{model}_{exp}"] = getattr(xr.where(d < 0, -d, 0), var_avg)(dim="week", skipna=skipna)
        # data_model[model] = xr.combine_by_coords([data_dict[f"{model}_historical"], data_dict[f"{model}_rcp45"], data_dict[f"{model}_rcp85"]], combine_attrs="override")
    # ds_calc = calc_annual_mean(
    #     data_dict=data_dict, seasonal=seasonal, skipna=skipna, var_avg=var_avg
    # )[variable]
    # # historical 10th percentile
    # # hist_p10 = (
    # #     ds_calc.sel(exp="historical")
    # #     .drop_vars("exp")
    # #     .chunk(dict(year=-1))
    # #     .quantile(dim="year", skipna=skipna, q=0.1)
    # # )
    # # historical mean
    # hist_mean = (
    #     ds_calc.sel(exp="historical")
    #     .drop_vars("exp")
    #     .mean(dim="year", skipna=skipna)
    # )
    # 30-day rolling mean of daily data

    # # for key in data_model:
    # #     d = data_model[key].rolling(time=30).mean(skipna=skipna)[variable]
    # #     # difference between daily rolling average and historical
    # #     # 10th percentile
    # #     # d = d - hist_p10
    # #     if seasonal:
    # #         ddict = {}
    # #         for s, m in zip(season_list, [[12, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]):
    # #             # ddict[s] = d.sel(time=(d.time.dt.month.isin(m))) - hist_mean.sel(season=s).drop_vars("season")
    # #             ddict[s] = d.sel(time=(d.time.dt.month.isin(m))) - d.sel(time=(d.time.dt.month.isin(m))).sel(exp="historical").drop_vars("exp")
    # #             ddict[s] = ddict[s].assign_coords(season=s).expand_dims(dim="season")
    # #         # d = xr.merge(ddict.values())
    # #         d = xr.combine_by_coords(ddict.values(), combine_attrs="override")
    # #     else:
    # #         d = d - d.sel(exp="historical").drop_vars("exp")
    # #     # keep only negative values and aggregate as yearly data
    # #     data_count[key] = xr.where(d < 0, 1, 0).groupby("time.year").sum(dim="time", skipna=skipna)
    # #     data_val[key] = getattr(xr.where(d < 0, -d, 0).groupby("time.year"), var_avg)(dim="time", skipna=skipna)

    # combine data
    data_count = xr.combine_by_coords(
        data_count.values(), combine_attrs="override"
    )
    data_val = xr.combine_by_coords(
        data_val.values(), combine_attrs="override"
    )
    data_val.rio.write_crs(data_dict[f"{model}_historical"].rio.crs, inplace=True)
    return data_count, data_val


def duration_surplus_biomass(data_dict, variable="bm_c", seasonal=False, skipna=None):
    data_count = {}
    for key in data_dict:
        d = data_dict[key][variable]
        data_count[key] = xr.where(d > 0, 1, 0).groupby("time.year").sum(dim="time", skipna=skipna)
    # combine data
    data_count = xr.combine_by_coords(
        data_count.values(), combine_attrs="override"
    )
    return data_count
