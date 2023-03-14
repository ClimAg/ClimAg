"""plot_facet_map.py

Helper functions to plot facet maps
"""

from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import climag.plot_configs as cplt


def plot_facet_map(
    data, var: str, boundary_data=None, cbar_levels=None,  # ticks=False
):
    """
    Create a facet plot of a variable from an Xarray dataset covering the
    Island of Ireland.

    Parameters
    ----------
    data : dataset (loaded using Xarray)
    boundary_data : Ireland boundary data (vector, loaded as a GeoPandas
        dataframe), e.g. from Ordnance Survey Ireland, NUTS (Eurostat)
    var : The variable to plot
    cbar_levels : Number of colour bar levels
    ticks : Whether to display latitude and longitude axes tick values
    """

    plot_transform = cplt.rotated_pole_transform(data)

    cbar_label = (
        f"{data[var].attrs['long_name']} [{data[var].attrs['units']}]"
    )  # colorbar label

    cmap = cplt.colormap_configs(var)

    if len(data["time"]) == 12:
        col_wrap = 4
        # y_ticks = [0, 4, 8]  # index of subplots with y tick labels
        # x_ticks = [8, 9, 10, 11]  # index of subplots with x tick labels
    elif len(data["time"]) == 30:
        col_wrap = 5
        # y_ticks = [0, 5, 10, 15, 20, 25]
        # x_ticks = [25, 26, 27, 28, 29]
    else:
        col_wrap = None
        # y_ticks = []
        # x_ticks = []

    fig = data[var].plot(
        x="rlon", y="rlat", col="time",
        col_wrap=col_wrap,
        cmap=cmap,
        robust=True,
        cbar_kwargs={"aspect": 40, "label": cbar_label},
        transform=plot_transform,
        subplot_kws={"projection": cplt.plot_projection},
        levels=cbar_levels
    )

    fig.set_xlabels("")
    fig.set_ylabels("")

    # for i, axs in enumerate(fig.axs.flat):
    for axs in fig.axs.flat:
        if boundary_data is None:
            axs.coastlines(
                resolution="10m", color="darkslategrey", linewidth=.5
            )
        else:
            boundary_data.to_crs(cplt.plot_projection).plot(
                ax=axs, edgecolor="darkslategrey", color="white",
                linewidth=.5
            )

        axs.set_xlim(-1.9, 1.6)
        axs.set_ylim(-2.1, 2.1)

        # axs.set_title(cordex_date_format(data.isel(time=i)))
        # # use gridlines to add tick labels (lon/lat)
        # if ticks and i in y_ticks:
        #     axs.gridlines(
        #         draw_labels=["y", "left"],
        #         ylocs=range(-90, 90, 1),
        #         color="None",
        #         linewidth=.5,
        #         x_inline=False,
        #         y_inline=False
        #     )
        # if ticks and i in x_ticks:
        #     axs.gridlines(
        #         draw_labels=["x", "bottom"],
        #         xlocs=range(-180, 180, 2),
        #         color="None",
        #         linewidth=.5,
        #         x_inline=False,
        #         y_inline=False
        #     )

    plt.show()


def plot_averages(
    data, var: str, averages: str, boundary_data=None, cbar_levels=None
):
    """
    Monthly, yearly, or seasonal averages plots

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

    cmap = cplt.colormap_configs(var)

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
        subplot_kws={"projection": cplt.plot_projection},
        levels=cbar_levels,
        xlim=(-1.9, 1.6),
        ylim=(-2.1, 2.1),
        aspect=.9
    )

    for i, axs in enumerate(fig.axs.flat):
        # boundary_data.to_crs(plot_projection).boundary.plot(
        #     ax=axs, color="darkslategrey", linewidth=.5
        # )
        if boundary_data is None:
            axs.coastlines(
                resolution="10m", color="darkslategrey", linewidth=.5
            )
        else:
            boundary_data.to_crs(cplt.plot_projection).plot(
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
            axs.set_title(f"{seasons[season]} ({season})")
        else:
            axs.set_title(
                str(data_weighted.isel({averages: i})[averages].values)
            )

    plt.show()


def plot_seasonal(
    data, var: str, stat="mean", cbar_levels=None, contour=False,
    boundary_data=None
):
    """
    Seasonal plots

    - https://docs.xarray.dev/en/stable/examples/monthly-means.html
    - https://ncar.github.io/esds/posts/2021/yearly-averages-xarray/

    Parameters
    ----------
    data : Xarray dataset
    var : The variable to plot (e.g. "T")
    stat : Statistic to plot; default is "mean" (weighted mean); other options
        are "std" (standard deviation), "0.9q" (90th percentile), "0.1q"
        (10th percentile), "min" (minimum), "max" (maximum), and "median"
        (median). Note that only the mean is weighted, i.e. the number of days
        in each month have been taken into account. The standard deviation is
        unbiased (using Delta Degrees of Freedom of 1). See the Xarray API
        docs for "DatasetGroupBy" for more info.
    cbar_levels : Number of discrete colour bar levels; if None, use a
        continuous colour bar
    boundary_data : Boundary data as a GeoPandas geodataframe; if None,
        Cartopy's coastlines are used
    """

    if stat == "mean":
        data_stat = cplt.weighted_average(data=data, averages="season")
    elif stat == "std":
        data_stat = data.groupby("time.season").std(dim="time", ddof=1)
    elif stat == "0.9q":
        data_stat = data.groupby("time.season").quantile(.9, dim="time")
    elif stat == "0.1q":
        data_stat = data.groupby("time.season").quantile(.1, dim="time")
    elif stat == "min":
        data_stat = data.groupby("time.season").min(dim="time")
    elif stat == "max":
        data_stat = data.groupby("time.season").max(dim="time")
    elif stat == "median":
        data_stat = data.groupby("time.season").median(dim="time")
    # elif stat == "sum":
    #     data_stat = data.groupby("time.season").sum(dim="time")

    plot_transform = cplt.rotated_pole_transform(data)

    cmap = cplt.colormap_configs(var)

    if contour:
        fig = data_stat[var].where(pd.notnull(data[var][0])).plot.contourf(
            x="rlon", y="rlat", col="season",
            col_wrap=2,
            cmap=cmap,
            robust=True,
            cbar_kwargs={
                "aspect": 20,
                "label": (
                    f"{data[var].attrs['long_name']} "
                    f"[{data[var].attrs['units']}]"
                )
            },
            transform=plot_transform,
            subplot_kws={"projection": cplt.plot_projection},
            levels=cbar_levels,
            xlim=(-1.9, 1.6),
            ylim=(-2.1, 2.1),
            aspect=.9
        )
    else:
        fig = data_stat[var].where(pd.notnull(data[var][0])).plot(
            x="rlon", y="rlat", col="season",
            col_wrap=2,
            cmap=cmap,
            robust=True,
            cbar_kwargs={
                "aspect": 20,
                "label": (
                    f"{data[var].attrs['long_name']} "
                    f"[{data[var].attrs['units']}]"
                )
            },
            transform=plot_transform,
            subplot_kws={"projection": cplt.plot_projection},
            levels=cbar_levels,
            xlim=(-1.9, 1.6),
            ylim=(-2.1, 2.1),
            aspect=.9
        )

    for i, axs in enumerate(fig.axs.flat):
        # add boundary
        if boundary_data is None:
            axs.coastlines(
                resolution="10m", color="darkslategrey", linewidth=.5
            )
        else:
            boundary_data.to_crs(cplt.plot_projection).plot(
                ax=axs, edgecolor="darkslategrey", color="white", linewidth=.5
            )

        # specify subplot titles
        axs.set_title(str(data_stat.isel({"season": i})["season"].values))


def plot_season_diff(data, var, boundary_data=None, stat="mean"):
    """
    Plot differences between weighted/unweighted mean and unbiased/biased
    standard deviation

    Parameters
    ----------
    data : Xarray dataset
    var : Variable to plot
    boundary_data : Boundary data as a GeoPandas geodataframe; if None,
        Cartopy's coastlines are used
    stat : Statistic; either "mean" (default) or "std" for standard deviation

    - https://docs.xarray.dev/en/stable/user-guide/computation.html
    - https://ncar.github.io/esds/posts/2021/yearly-averages-xarray/
    - https://docs.xarray.dev/en/stable/examples/monthly-means.html
    - https://ncar.github.io/esds/posts/2020/Time/
    """

    notnull = pd.notnull(data[var][0])

    if stat == "mean":
        # weighted average
        data_1 = cplt.weighted_average(data=data, averages="season")
        # unweighted average
        data_2 = data.groupby("time.season").mean(dim="time")
        titles = "Weighted mean", "Unweighted mean"
    elif stat == "std":
        # unbiased SD
        data_1 = data.groupby("time.season").std(dim="time", ddof=1)
        # biased SD
        data_2 = data.groupby("time.season").std(dim="time")
        titles = "Unbiased SD", "Biased SD"
    data_diff = data_1 - data_2  # difference between the stats

    cmap = cplt.colormap_configs(var)

    fig, axs = plt.subplots(
        nrows=4, ncols=3, figsize=(10, 10),
        subplot_kw={"projection": cplt.plot_projection}
    )

    for i, season in enumerate(("DJF", "MAM", "JJA", "SON")):
        data_1[var].sel(season=season).where(notnull).plot(
            ax=axs[i, 0],
            cmap=cmap,
            # extend="both",
            robust=True,
            transform=cplt.rotated_pole_transform(data),
            cbar_kwargs={"label": None},
            xlim=(-1.9, 1.6),
            ylim=(-2.1, 2.1)
        )

        data_2[var].sel(season=season).where(notnull).plot(
            ax=axs[i, 1],
            cmap=cmap,
            # extend="both",
            robust=True,
            transform=cplt.rotated_pole_transform(data),
            cbar_kwargs={"label": None},
            xlim=(-1.9, 1.6),
            ylim=(-2.1, 2.1)
        )

        data_diff[var].sel(season=season).where(notnull).plot(
            ax=axs[i, 2],
            cmap="RdBu_r",
            # extend="both",
            robust=True,
            transform=cplt.rotated_pole_transform(data),
            cbar_kwargs={"label": None},
            xlim=(-1.9, 1.6),
            ylim=(-2.1, 2.1)
        )

        axs[i, 0].set_ylabel(season)
        axs[i, 0].set_yticks([])

    for axis in axs.flat:
        if boundary_data is None:
            axis.coastlines(
                resolution="10m", color="darkslategrey", linewidth=.5
            )
        else:
            boundary_data.to_crs(cplt.plot_projection).plot(
                ax=axis, edgecolor="darkslategrey", color="white",
                linewidth=.5
            )

        axis.set_title(None)

    axs[0, 0].set_title(titles[0])
    axs[0, 1].set_title(titles[1])
    axs[0, 2].set_title("Difference")

    plt.tight_layout()

    fig.suptitle(
        f"{data[var].attrs['long_name']} [{data[var].attrs['units']}]",
        fontsize=16, y=1.02
    )
    plt.show()


def plot_season_diff_hist_rcp(data, var, boundary_data=None, stat="mean"):
    """
    Plot differences between historical and rcp45/rcp8.5 data

    Parameters
    ----------
    data : a tuple of two Xarray datasets to compare
    var : Variable to plot
    boundary_data : Boundary data as a GeoPandas geodataframe; if None,
        Cartopy's coastlines are used
    stat : Statistic to plot; default is "mean" (weighted mean); other options
        are "std" (standard deviation), "0.9q" (90th percentile), "0.1q"
        (10th percentile), "min" (minimum), "max" (maximum), and "median"
        (median). Note that only the mean is weighted, i.e. the number of days
        in each month have been taken into account. The standard deviation is
        unbiased (using Delta Degrees of Freedom of 1). See the Xarray API
        docs for "DatasetGroupBy" for more info.

    - https://docs.xarray.dev/en/stable/user-guide/computation.html
    - https://ncar.github.io/esds/posts/2021/yearly-averages-xarray/
    - https://docs.xarray.dev/en/stable/examples/monthly-means.html
    - https://ncar.github.io/esds/posts/2020/Time/
    """

    notnull = pd.notnull(data[0][var][0])

    if stat == "mean":
        data_h = cplt.weighted_average(data=data[0], averages="season")
        data_f = cplt.weighted_average(data=data[1], averages="season")
    elif stat == "std":
        data_h = data[0].groupby("time.season").std(dim="time", ddof=1)
        data_f = data[1].groupby("time.season").std(dim="time", ddof=1)
    # elif stat == "0.9q":
    #     data_h = data[0].groupby("time.season").quantile(.9, dim="time")
    #     data_f = data[1].groupby("time.season").quantile(.9, dim="time")
    # elif stat == "0.1q":
    #     data_h = data[0].groupby("time.season").quantile(.1, dim="time")
    #     data_f = data[1].groupby("time.season").quantile(.1, dim="time")
    elif stat == "min":
        data_h = data[0].groupby("time.season").min(dim="time")
        data_f = data[1].groupby("time.season").min(dim="time")
    elif stat == "max":
        data_h = data[0].groupby("time.season").max(dim="time")
        data_f = data[1].groupby("time.season").max(dim="time")
    elif stat == "median":
        data_h = data[0].groupby("time.season").median(dim="time")
        data_f = data[1].groupby("time.season").median(dim="time")
    data_diff = data_f - data_h

    cmap = cplt.colormap_configs(var)

    fig, axs = plt.subplots(
        nrows=4, ncols=3, figsize=(10, 10),
        subplot_kw={"projection": cplt.plot_projection}
    )

    for i, season in enumerate(("DJF", "MAM", "JJA", "SON")):
        data_h[var].sel(season=season).where(notnull).plot.contourf(
            ax=axs[i, 0],
            cmap=cmap,
            # extend="both",
            robust=True,
            transform=cplt.rotated_pole_transform(data[0]),
            cbar_kwargs={"label": None},
            xlim=(-1.9, 1.6),
            ylim=(-2.1, 2.1)
        )

        data_f[var].sel(season=season).where(notnull).plot.contourf(
            ax=axs[i, 1],
            cmap=cmap,
            robust=True,
            transform=cplt.rotated_pole_transform(data[0]),
            cbar_kwargs={"label": None},
            xlim=(-1.9, 1.6),
            ylim=(-2.1, 2.1)
        )

        data_diff[var].sel(season=season).where(notnull).plot.contourf(
            ax=axs[i, 2],
            cmap=cmap,
            robust=True,
            transform=cplt.rotated_pole_transform(data[0]),
            cbar_kwargs={"label": None},
            xlim=(-1.9, 1.6),
            ylim=(-2.1, 2.1)
        )

        axs[i, 0].set_ylabel(season, fontweight="semibold")
        axs[i, 0].set_yticks([])

    for axis in axs.flat:
        if boundary_data is None:
            axis.coastlines(
                resolution="10m", color="darkslategrey", linewidth=.5
            )
        else:
            boundary_data.to_crs(cplt.plot_projection).plot(
                ax=axis, edgecolor="darkslategrey", color="white",
                linewidth=.5
            )

        axis.set_title(None)

    try:
        axs[0, 0].set_title(data[0].attrs["dataset"].split("_")[4])
        axs[0, 1].set_title(data[1].attrs["dataset"].split("_")[4])
    except KeyError:
        axs[0, 0].set_title(data[0].attrs["input_dataset"].split("_")[4])
        axs[0, 1].set_title(data[1].attrs["input_dataset"].split("_")[4])
    axs[0, 2].set_title("difference")

    plt.tight_layout()

    # try:
    #     suptitle = (
    #         f"{data[0][var].attrs['long_name']} "
    #         f"[{data[0][var].attrs['units']}] ("
    #         + ", ".join(data[0].attrs["dataset"].split("_")[1:4]) + ")"
    #     )
    # except KeyError:
    #     suptitle = (
    #         f"{data[0][var].attrs['long_name']} "
    #         f"[{data[0][var].attrs['units']}] ("
    #         + ", ".join(data[0].attrs["input_dataset"].split("_")[1:4]) + ")"
    #     )

    suptitle = (
        f"{data[0][var].attrs['long_name']} [{data[0][var].attrs['units']}]"
    )

    fig.suptitle(suptitle, fontsize=16, y=1.02)
    plt.show()
