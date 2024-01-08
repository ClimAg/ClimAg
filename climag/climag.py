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

import climag.climag as cplt

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


def weighted_average(data, averages: str):
    """
    Calculate the weighted average

    - https://docs.xarray.dev/en/stable/user-guide/computation.html
    - https://ncar.github.io/esds/posts/2021/yearly-averages-xarray/
    - https://docs.xarray.dev/en/stable/examples/monthly-means.html
    - https://ncar.github.io/esds/posts/2020/Time/
    """

    # calculate the weights by grouping month length by season or month
    weights = (
        data.time.dt.days_in_month.groupby(f"time.{averages}")
        / data.time.dt.days_in_month.groupby(f"time.{averages}").sum()
    )

    # test that the sum of weights for each year/season/month is one
    np.testing.assert_allclose(
        weights.groupby(f"time.{averages}").sum().values,
        np.ones(len(set(weights[averages].values))),
    )

    # calculate the weighted average
    data_weighted = (
        (data * weights).groupby(f"time.{averages}").sum(dim="time")
    )

    return data_weighted
