"""Helper functions to plot datasets

"""

import climag.climag as cplt
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import xarray as xr
import cartopy.crs as ccrs
import numpy as np
import glob
import os
import warnings
from itertools import product
import geopandas as gpd
import rasterio as rio
from matplotlib import patheffects
# from dateutil.parser import parse

warnings.filterwarnings(
    action="ignore", category=RuntimeWarning, module="dask"
)


def colormap_configs(var):
    """
    Configure colourmap for each variable
    """

    if var in ("PP", "TOT_PREC", "pr", "tp"):
        cmap = "mako_r"  # precipitation
    elif var in ("wr", "r", "u", "v"):
        cmap = "GnBu"  # wind speed, humidity, water reserves
    elif var in (
        "T",
        "PAR",
        "ASWDIR_S",
        "ASWDIFD_S",
        "ASWDIFU_S",
        "ASOB_S",
        "T_2M",
        "rsds",
        "tas",
        "t",
        "grad",
        "tmax",
        "tmin",
        "nswrs",
        "nlwrs",
    ):
        cmap = "Spectral_r"  # temperature and radiation
    elif var in ("PET", "aet", "ET", "evspsblpot"):
        cmap = "BrBG"  # evapotranspiration
    elif var in ("ALB_RAD", "pres", "env"):
        cmap = "flare"  # albedo, pressure, environmental limitation
    else:
        cmap = "PRGn"
    return cmap


def plot_single_map(data, var, boundary_data=None, cbar_levels=None, contour=False):
    """
    Create an individual plot of a climate data variable covering the Island
    of Ireland.

    Parameters
    ----------
    data : climate model dataset (loaded using Xarray)
    var : The variable to plot
    cbar_levels : Number of colour bar levels
    title : Plot title; if "default", use the default plot title
    contour : Create a filled contour plot
    """

    plot_transform = cplt.rotated_pole_transform(data)

    cbar_label = f"{data[var].attrs['long_name']} [{data[var].attrs['units']}]"

    cmap = colormap_configs(var)

    plt.figure(figsize=(7, 7))

    axs = plt.axes(projection=cplt.projection_hiresireland)

    # plot data for the variable
    if contour:
        data[var].plot.contourf(
            ax=axs,
            cmap=cmap,
            x="rlon",
            y="rlat",
            robust=True,
            cbar_kwargs={"label": cbar_label},
            transform=plot_transform,
            levels=cbar_levels,
        )
    else:
        data[var].plot(
            ax=axs,
            cmap=cmap,
            x="rlon",
            y="rlat",
            robust=True,
            cbar_kwargs={"label": cbar_label},
            transform=plot_transform,
            levels=cbar_levels,
        )

    # add boundaries
    if boundary_data is None:
        axs.coastlines(resolution="10m", color="darkslategrey", linewidth=0.75)
    else:
        boundary_data.to_crs(projection_hiresireland).plot(
            ax=axs, edgecolor="darkslategrey", color="white", linewidth=0.75
        )

    axs.set_title(None)

    plt.axis("equal")
    plt.tight_layout()
    plt.xlim(-1.5, 1.33)
    plt.ylim(-2.05, 2.05)

    # specify gridline spacing and labels
    axs.gridlines(
        draw_labels={"bottom": "x", "left": "y"},
        xlocs=range(-180, 180, 2),
        ylocs=range(-90, 90, 1),
        color="lightslategrey",
        linewidth=0.5,
        x_inline=False,
        y_inline=False,
    )

    plt.show()


def plot_averages(
    data, var: str, averages: str, boundary_data=None, cbar_levels=None
):
    """Monthly, yearly, or seasonal averages plots

    - https://docs.xarray.dev/en/stable/examples/monthly-means.html
    - https://ncar.github.io/esds/posts/2021/yearly-averages-xarray/

    Parameters
    ----------
    data : Xarray dataset
    var : The variable to plot (e.g. "T")
    boundary_data : GeoPandas boundary vector data
    averages : Seasonal ("season"), annual ("year"), or monthly ("month")
        averages
    cbar_levels : Number of discrete colour bar levels; if None, use a
        continuous colour bar
    """

    # calculate the weighted average
    data_weighted = cplt.weighted_average(data=data, averages=averages)

    plot_transform = cplt.rotated_pole_transform(data)

    cmap = colormap_configs(var)

    # configure number of columns of the plot and aspect of the colour bar
    if averages == "month":
        columns_cbar_aspect = 4, 25
    elif averages == "year":
        columns_cbar_aspect = 6, 35
    else:
        columns_cbar_aspect = 2, 20
        # sort seasons
        data = data.reindex(season=["DJF", "MAM", "JJA", "SON"])

    fig = (
        data_weighted[var]
        .where(pd.notnull(data[var][0]))
        .plot(
            x="rlon",
            y="rlat",
            col=averages,
            col_wrap=columns_cbar_aspect[0],
            cmap=cmap,
            robust=True,
            cbar_kwargs={
                "aspect": columns_cbar_aspect[1],
                "label": (
                    f"{data[var].attrs['long_name']} [{data[var].attrs['units']}]"
                ),
            },
            transform=plot_transform,
            subplot_kws={"projection": cplt.projection_hiresireland},
            levels=cbar_levels,
            xlim=(-1.9, 1.6),
            ylim=(-2.1, 2.1),
            aspect=0.9,
        )
    )

    for i, axs in enumerate(fig.axs.flat):
        # boundary_data.to_crs(projection_hiresireland).boundary.plot(
        #     ax=axs, color="darkslategrey", linewidth=.5
        # )
        if boundary_data is None:
            axs.coastlines(
                resolution="10m", color="darkslategrey", linewidth=0.5
            )
        else:
            boundary_data.to_crs(cplt.projection_hiresireland).plot(
                ax=axs, color="white", edgecolor="darkslategrey", linewidth=0.5
            )

        # if averages == "month":
        #     axs.set_title(
        #         datetime.strptime(
        #             str(data_weighted.isel({averages: i})[averages].values),
        #             "%m",
        #         )
        #         .strftime("%-b")
        #         .upper()
        #     )
        # else:
        #     axs.set_title(data[averages][i].values)
        # elif averages == "season":
        #     seasons = {
        #         "DJF": "Winter",
        #         "MAM": "Spring",
        #         "JJA": "Summer",
        #         "SON": "Autumn"
        #     }
        #     season = str(data_weighted.isel({averages: i})[averages].values)
        #     axs.set_title(f"{seasons[season]} ({season})")
        # else:
        #     axs.set_title(
        #         str(data_weighted.isel({averages: i})[averages].values)
        #     )

    plt.show()


def boxplot_data(datasets, varlist, lonlat):
    """
    Process time series data for a given location to create box plots

    Parameters
    ----------
    data : Dictionary of Xarray climate model datasets
    var : The variable to be plotted from the dataset
    lonlat : A tuple of the location's coordinates (longitude, latitude)

    Returns
    -------
    - A dictionary of dataframes (one for each variable) of the boxplot data
      to be plotted
    """

    data_all = {}

    for var in varlist:
        data_all[var] = {}

    for key in datasets.keys():
        cds = rotated_pole_point(
            data=datasets[key], lon=lonlat[0], lat=lonlat[1]
        )
        data_all[key] = datasets[key].sel(
            {"rlon": cds[0], "rlat": cds[1]}, method="nearest"
        )

        for var in varlist:
            data_all[var][key] = pd.DataFrame(
                {var: data_all[key][var]}
            ).assign(
                dataset=f"{key.split('_')[0]}\n{key.split('_')[1]}",
                exp=key.split("_")[2],
                legend=(
                    f"{key.split('_')[0]}\n{key.split('_')[1]}\n"
                    f"{key.split('_')[2]}"
                ),
            )

    for var in varlist:
        data_all[var] = pd.concat([v for k, v in data_all[var].items()])

    return data_all


def boxplot_all(data, var, title, showfliers=False, figsize=(12, 5)):
    """
    Generate box plots for Xarray datasets stored in a dictionary

    Parameters
    ----------
    data : A dictionary of Xarray datasets
    var : The variable from the Xarray dataset to plot
    title : Plot title
    showfliers : Show outliers (default is False)
    figsize : Size of the plot figure
    """

    plt.figure(figsize=figsize)
    sns.boxplot(
        data,
        x="dataset",
        y=var,
        hue="exp",
        showfliers=showfliers,
        showmeans=True,
        palette="Pastel1",
        meanprops={
            "markeredgecolor": "darkslategrey",
            "marker": "d",
            "markerfacecolor": "white",
            "markersize": 7.5,
        },
        boxprops={"edgecolor": "darkslategrey"},
        medianprops={"color": "darkslategrey"},
        whiskerprops={"color": "darkslategrey"},
        capprops={"color": "darkslategrey"},
    )
    plt.xlabel("")
    plt.ylabel("")
    plt.title(title)
    plt.legend(title=None, loc="upper right")
    plt.tight_layout()
    plt.show()


def colorbar_levels(maximum, levels=22):
    """
    Create a list of diverging colourbar levels based on the maximum value and
    number of levels
    """

    # levels = [
    #     i for i in [
    #         -maximum + num * n for n in range(int(maximum / num * 2 + 1))
    #     ] if i != 0
    # ]
    return [-maximum + maximum / ((levels - 1) / 2) * n for n in range(levels)]


def colorbar_ticks(maximum):
    """
    Colourbar tick labels
    """

    return [-maximum, 0, maximum]
