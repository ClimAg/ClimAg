"""plot_configs.py

Functions to plot climate model datasets, e.g. CORDEX
"""

from datetime import datetime
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from dateutil.parser import parse

# set plot projection to the projection of the HiResIreland dataset
plot_projection = ccrs.RotatedPole(
    pole_longitude=172.100006103516, pole_latitude=36.5999984741211
)

# define Lambert Conformal Conic projection for plots and transformations
# using MERA GRIB metadata
lambert_conformal = ccrs.LambertConformal(
    false_easting=1481641.6769636814,
    false_northing=537326.0638850163,
    standard_parallels=[53.5],
    central_longitude=5.0,
    central_latitude=53.5
)

# seaborn colourmaps
cmap_mako_r = sns.color_palette("mako_r", as_cmap=True)
# cmap_crest = sns.color_palette("crest", as_cmap=True)
cmap_flare = sns.color_palette("flare", as_cmap=True)


def colormap_configs(var):
    """
    Configure colourmap for each variable
    """

    if var in ("PP", "TOT_PREC", "pr"):
        cmap = cmap_mako_r
    elif var in ("wr", "env"):
        cmap = "GnBu"
    elif var in (
        "T", "PAR",
        "ASWDIR_S", "ASWDIFD_S", "ASWDIFU_S", "ASOB_S", "T_2M",
        "rsds", "tas"
    ):
        cmap = "Spectral_r"
    elif var in ("PET", "aet", "ALB_RAD", "w", "ET", "evspsblpot"):
        cmap = cmap_flare
    else:
        cmap = "YlGn"
    return cmap


def longitude_tick_format(x, pos):
    """
    Return the longitude in degrees west.
    The two arguments are the value and tick position.
    https://matplotlib.org/stable/gallery/ticks/tick-formatters.html
    """

    return "{:,.0f}°W".format(x * -1)


def latitude_tick_format(x, pos):
    """
    Return the latitude in degrees north.
    The two arguments are the value and tick position.
    https://matplotlib.org/stable/gallery/ticks/tick-formatters.html
    """

    return "{:.0f}°N".format(x)


def rotated_pole_point(data, lon, lat):
    """
    Convert the longitude and latitude of a specific point to rotated pole
    coordinates used in the input data.

    Parameters
    ----------
    data : input climate data which uses rotated pole coordinates
    lon : longitude of the point
    lat : latitude of the point

    Returns
    -------
    - transformed longitude and latitude in rotated pole coordinates
    """

    if data.rio.crs is None:
        pole_longitude = (
            data["rotated_pole"].attrs["grid_north_pole_longitude"]
        )
        pole_latitude = data["rotated_pole"].attrs["grid_north_pole_latitude"]
    else:
        pole_longitude = (
            data.rio.crs.to_dict(
                projjson=True
            )["conversion"]["parameters"][1]["value"]
        )
        pole_latitude = (
            data.rio.crs.to_dict(
                projjson=True
            )["conversion"]["parameters"][0]["value"]
        )
    rp_cds = ccrs.RotatedGeodetic(
        pole_longitude=pole_longitude, pole_latitude=pole_latitude,
    ).transform_point(x=lon, y=lat, src_crs=ccrs.Geodetic())
    return rp_cds[0], rp_cds[1]


def rotated_pole_transform(data):
    """
    Rotated pole transform for plotting CORDEX data.

    Parameters
    ----------
    data : input CORDEX data

    Returns
    -------
    - rotated pole transform
    """

    if data.rio.crs is None:
        pole_longitude = (
            data["rotated_pole"].attrs["grid_north_pole_longitude"]
        )
        pole_latitude = data["rotated_pole"].attrs["grid_north_pole_latitude"]
    else:
        pole_longitude = (
            data.rio.crs.to_dict(
                projjson=True
            )["conversion"]["parameters"][1]["value"]
        )
        pole_latitude = (
            data.rio.crs.to_dict(
                projjson=True
            )["conversion"]["parameters"][0]["value"]
        )
    transform = ccrs.RotatedPole(
        pole_longitude=pole_longitude, pole_latitude=pole_latitude
    )
    return transform


def hiresireland_date_format(data):
    """
    Format date
    """

    date = datetime.strftime(parse(str(data["time"].values)), "%-d %b %Y")
    return date


# def cordex_date_format(data):
#     """
#     Format date
#     """

#     # if data.attrs["frequency"] == "mon":
#     #     date_format = "%b %Y"
#     # elif data.attrs["frequency"] == "day":
#     #     date_format = "%-d %b %Y"
#     # else:
#     #     date_format = "%Y-%m-%d %H:%M:%S"
#     date = datetime.strftime(parse(str(data["time"].values)), "%-d %b %Y")
#     return date


def plot_facet_map(data, var, boundary_data, cbar_levels=None, ticks=False):
    """
    Create a facet plot of a variable from an xarray dataset covering the
    Island of Ireland.

    Parameters
    ----------
    data : dataset (loaded using xarray)
    boundary_data : Ireland boundary data (vector, loaded as a GeoPandas
        dataframe), e.g. from Ordnance Survey Ireland, NUTS (Eurostat)
    var : The variable to plot
    cbar_levels : Number of colour bar levels
    ticks : Whether to display latitude and longitude axes tick values
    """

    plot_transform = rotated_pole_transform(data)

    cbar_label = (
        f"{data[var].attrs['long_name']} [{data[var].attrs['units']}]"
    )  # colorbar label

    cmap = colormap_configs(var)

    if len(data["time"]) == 12:
        col_wrap = 4
        y_ticks = [0, 4, 8]  # index of subplots with y tick labels
        x_ticks = [8, 9, 10, 11]  # index of subplots with x tick labels
    elif len(data["time"]) == 30:
        col_wrap = 5
        y_ticks = [0, 5, 10, 15, 20, 25]
        x_ticks = [25, 26, 27, 28, 29]
    else:
        col_wrap = None
        y_ticks = []
        x_ticks = []

    fig = data[var].plot(
        x="rlon", y="rlat", col="time",
        col_wrap=col_wrap,
        cmap=cmap,
        robust=True,
        cbar_kwargs={"aspect": 40, "label": cbar_label},
        transform=plot_transform,
        subplot_kws={"projection": plot_projection},
        levels=cbar_levels
    )

    fig.set_xlabels("")
    fig.set_ylabels("")

    for i, axs in enumerate(fig.axs.flat):
        boundary_data.to_crs(plot_projection).boundary.plot(
            ax=axs, color="darkslategrey", linewidth=.5
        )
        # axs.set_title(cordex_date_format(data.isel(time=i)))
        axs.set_xlim(-1.9, 1.6)
        axs.set_ylim(-2.1, 2.1)
        # use gridlines to add tick labels (lon/lat)
        if ticks and i in y_ticks:
            axs.gridlines(
                draw_labels=["y", "left"],
                ylocs=range(-90, 90, 1),
                color="None",
                linewidth=.5,
                x_inline=False,
                y_inline=False
            )
        if ticks and i in x_ticks:
            axs.gridlines(
                draw_labels=["x", "bottom"],
                xlocs=range(-180, 180, 2),
                color="None",
                linewidth=.5,
                x_inline=False,
                y_inline=False
            )

    plt.show()


def plot_map(data, var, cbar_levels=None, title="default"):
    """
    Create an individual plot of a climate data variable covering the Island
    of Ireland.

    Parameters
    ----------
    data : climate model dataset (loaded using xarray)
    var : The variable to plot
    cbar_levels : Number of colour bar levels
    title : Plot title; if "default", use the default plot title
    """

    plot_transform = rotated_pole_transform(data)

    cbar_label = (
        f"{data[var].attrs['long_name']} [{data[var].attrs['units']}]"
    )  # colorbar label

    cmap = colormap_configs(var)

    plt.figure(figsize=(7, 7))

    axs = plt.axes(projection=plot_projection)

    # plot data for the variable
    data[var].plot(
        ax=axs,
        cmap=cmap,
        x="rlon",
        y="rlat",
        robust=True,
        cbar_kwargs={"label": cbar_label},
        transform=plot_transform,
        levels=cbar_levels
    )

    # add boundaries
    axs.coastlines(resolution="10m", color="darkslategrey", linewidth=.75)

    if title != "default":
        axs.set_title(title)

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
        linewidth=.5,
        x_inline=False,
        y_inline=False
    )

    plt.show()


def plot_averages(
    data, var: str, averages: str, boundary_data, cbar_levels=None
):
    """
    Monthly or seasonal averages plots

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
    title : Main plot title
    """

    # calculate the weights by grouping month length by year/season/month
    weights = (
        data.time.dt.days_in_month.groupby(f"time.{averages}") /
        data.time.dt.days_in_month.groupby(f"time.{averages}").sum()
    )

    # test that the sum of weights for each year/season/month is one
    np.testing.assert_allclose(
        weights.groupby(f"time.{averages}").sum().values,
        np.ones(len(set(weights[averages].values)))
    )

    # calculate the weighted average
    data_weighted = (
        (data * weights).groupby(f"time.{averages}").sum(dim="time")
    )

    plot_transform = rotated_pole_transform(data)

    cmap = colormap_configs(var)

    # configure number of columns of the plot and aspect of the colour bar
    if averages == "month":
        columns_cbar_aspect = 4, 25
    elif averages == "year":
        columns_cbar_aspect = 6, 35
    else:
        columns_cbar_aspect = 2, 20

    fig = data_weighted[var].where(pd.notnull(data[var][0])).plot(
        x="rlon", y="rlat", col=averages,
        col_wrap=columns_cbar_aspect[0],
        cmap=cmap,
        robust=True,
        cbar_kwargs={
            "aspect": columns_cbar_aspect[1],
            "label": (
                f"{data[var].attrs['long_name']} [{data[var].attrs['units']}]"
            )
        },
        transform=plot_transform,
        subplot_kws={"projection": plot_projection},
        levels=cbar_levels,
        xlim=(-1.9, 1.6),
        ylim=(-2.1, 2.1),
        aspect=.9
    )

    for i, axs in enumerate(fig.axs.flat):
        # boundary_data.to_crs(plot_projection).boundary.plot(
        #     ax=axs, color="darkslategrey", linewidth=.5
        # )
        boundary_data.to_crs(plot_projection).plot(
            ax=axs, color="white", edgecolor="darkslategrey", linewidth=.5
        )
        if averages == "month":
            axs.set_title(
                datetime.strptime(
                    str(data_weighted.isel({averages: i})[averages].values),
                    "%m"
                ).strftime("%-b").upper()
            )
        elif averages == "season":
            seasons = {
                "DJF": "Winter",
                "MAM": "Spring",
                "JJA": "Summer",
                "SON": "Autumn"
            }
            season = str(data_weighted.isel({averages: i})[averages].values)
            axs.set_title(
                f"{seasons[season]} ({season})"
            )
        else:
            axs.set_title(
                str(data_weighted.isel({averages: i})[averages].values)
            )

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
                {var: data_all[key][var]}).assign(
                    dataset=f"{key.split('_')[0]}\n{key.split('_')[1]}",
                    exp=key.split("_")[2]
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
        data, x="dataset", y=var, hue="exp", showfliers=showfliers,
        showmeans=True, notch=True, palette="viridis",
        meanprops={
            "markeredgecolor": "darkslategrey",
            "marker": "d",
            "markerfacecolor": "white",
            "markersize": 7.5
        },
        boxprops={"edgecolor": "white"},
        medianprops={"color": (1, 1, 0, 0)},  # transparent
        # whiskerprops={"color": "darkslategrey"},
        # capprops={"color": "darkslategrey"},
    )
    plt.xlabel("")
    plt.ylabel("")
    plt.title(title)
    plt.legend(title=None)
    plt.tight_layout()
    plt.show()
