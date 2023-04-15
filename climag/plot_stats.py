"""plot_stats.py

Helper functions to plot statistics
"""

import glob
import os
from itertools import product
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio as rio
import xarray as xr
from matplotlib import patheffects
import climag.plot_configs as cplt

# Ireland boundary
ie_bbox = gpd.read_file(
    os.path.join("data", "boundaries", "boundaries.gpkg"),
    layer="ne_10m_land_2157_IE_BBOX_DIFF"
)
ie = gpd.read_file(
    os.path.join("data", "boundaries", "boundaries.gpkg"),
    layer="NUTS_RG_01M_2021_2157_IE"
)

season_list = ["DJF", "MAM", "JJA", "SON"]
model_list = ["CNRM-CM5", "EC-EARTH", "HadGEM2-ES", "MPI-ESM-LR"]


##########################################################################
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


##########################################################################
def hist_rcp_diff(data):
    """
    Calculate difference between historical and rcp45/rcp85
    """

    data_out = xr.combine_by_coords([
        (
            data.sel(exp="rcp45") - data.sel(exp="historical")
        ).assign_coords(exp="rcp45 - historical").expand_dims(dim="exp"),
        (
            data.sel(exp="rcp85") - data.sel(exp="historical")
        ).assign_coords(exp="rcp85 - historical").expand_dims(dim="exp")
    ])
    # reassign attributes and CRS
    for var in data_out.data_vars:
        data_out[var].attrs = data[var].attrs
    data_out.rio.write_crs(data.rio.crs)

    return data_out


def hist_rcp_stats_data(dataset, stat, diff=True):
    """
    Prepare simulation results from climate model datasets for plotting
    differences between historical and rcp45/rcp85
    """

    data = {}

    for x in ["season", "cumulative"]:
        data[f"{dataset}_{x[0]}"] = xr.open_mfdataset(
            glob.glob(
                os.path.join(
                    "data", "ModVege", "stats", f"IE_{dataset}*{stat}_{x}.nc"
                )
            ),
            decode_coords="all", chunks="auto"
        )
        if diff:
            # difference between hist and rcp45/rcp85
            data[f"{dataset}_{x[0]}"] = hist_rcp_diff(
                data=data[f"{dataset}_{x[0]}"]
            )

    return data


def plot_all(data, var, season, levels=None, ticks=None):
    """
    Helper function to plot facet maps
    """

    # dummy data to crop out null grid cells
    if round(data.rio.resolution()[0], 2) == 0.11:
        notnull = xr.open_dataset(
            os.path.join(
                "data", "ModVege", "EURO-CORDEX", "historical", "CNRM-CM5",
                "modvege_IE_EURO-CORDEX_RCA4_CNRM-CM5_historical_1976.nc"
            ),
            decode_coords="all", chunks="auto"
        )
    else:
        notnull = xr.open_dataset(
            os.path.join(
                "data", "ModVege", "HiResIreland", "historical", "CNRM-CM5",
                "modvege_IE_HiResIreland_COSMO5_CNRM-CM5_historical_1976.nc"
            ),
            decode_coords="all", chunks="auto"
        )
    notnull = pd.notnull(notnull["gro"].isel(time=0))

    plot_transform = cplt.rotated_pole_transform(data)

    cbar_kwargs = {
        "label": f"Difference [{data[var].attrs['units']}]",
        "aspect": 30,
        "location": "bottom",
        "fraction": 0.085,
        "shrink": 0.85,
        "pad": 0.05,
        "extendfrac": "auto"
    }

    if ticks is not None:
        cbar_kwargs["ticks"] = ticks

    if len(data["exp"]) == 3:
        cmap = cplt.colormap_configs(var)
        figsize = (12, 13.75)
        extend = "max"
        robust = False
    else:
        cmap = "BrBG"
        figsize = (12, 9.25)
        extend = "both"
        robust = True

    if season is not None:
        data = data.sel(season=season)

    fig = data[var].where(notnull).plot.contourf(
        x="rlon", y="rlat", col="model", row="exp",
        cmap=cmap,
        extend=extend,
        robust=robust,
        cbar_kwargs=cbar_kwargs,
        transform=plot_transform,
        subplot_kws={"projection": cplt.plot_projection},
        levels=levels,
        xlim=(-1.775, 1.6),
        ylim=(-2.1, 2.1),
        figsize=figsize
    )

    fig.set_titles("{value}", weight="semibold", fontsize=14)

    title_print = (
        "_" * 80 +
        f"\nDifference in average {data[var].attrs['long_name'].lower()}"
    )
    if season is not None:
        print(title_print + f" in {season}\n" + "_" * 80)
    else:
        print(title_print + "\n" + "_" * 80)
    # plt.suptitle(title_print, size=16, y=0.8)

    # add boundary
    for (col, model), (row, exp) in product(
        enumerate(model_list), enumerate(data["exp"].values)
    ):
        try:
            ie_bbox.to_crs(cplt.plot_projection).plot(
                ax=fig.axs[row][col], edgecolor="darkslategrey", color="white",
                linewidth=.5
            )
        except NameError:
            fig.axs[row][col].coastlines(
                resolution="10m", color="darkslategrey", linewidth=.5
            )

        # add max/min values as annotations
        da = data.sel(model=model, exp=exp)[var]
        # clip to remove cells over the coastline
        da = da.rio.clip(ie.buffer(1).to_crs(data.rio.crs))
        da_max = da.isel(da.compute().argmax(dim=["rlon", "rlat"]))
        da_min = da.isel(da.compute().argmin(dim=["rlon", "rlat"]))
        da = gpd.GeoDataFrame(
            {
                "name": ["max", "min"],
                "value": [da_max.values, da_min.values],
                "geometry": gpd.GeoSeries.from_wkt([
                    f"POINT ({da_max['rlon'].values} {da_max['rlat'].values})",
                    f"POINT ({da_min['rlon'].values} {da_min['rlat'].values})"
                ])
            },
            crs=data.rio.crs
        )

        da = da.to_crs(cplt.plot_projection)
        da_max = da.loc[[0]]
        da_min = da.loc[[1]]
        da_max.plot(
            ax=fig.axs[row][col],
            color="#01665e", marker="o", edgecolor="white"
        )
        da_min.plot(
            ax=fig.axs[row][col],
            color="#8c510a", marker="o", edgecolor="white"
        )
        da_max["value"] = da_max["value"].astype(float).round(1)
        da_min["value"] = da_min["value"].astype(float).round(1)

        for xy, lab in zip(
            zip(da_max["geometry"].x - 0.3, da_max["geometry"].y),
            da_max["value"]
        ):
            fig.axs[row][col].annotate(
                text=lab, xy=xy, ha="center", va="center", size=12.5,
                weight="semibold", color="#01665e",
                path_effects=[
                    patheffects.withStroke(linewidth=3, foreground="white")
                ]
            )
        for xy, lab in zip(
            zip(da_min["geometry"].x + 0.3, da_min["geometry"].y),
            da_min["value"]
        ):
            fig.axs[row][col].annotate(
                text=lab, xy=xy, ha="center", va="center", size=12.5,
                weight="semibold", color="#8c510a",
                path_effects=[
                    patheffects.withStroke(linewidth=3, foreground="white")
                ]
            )

    plt.show()


##########################################################################
def hist_obs_diff(stat, dataset):
    """
    Prepare data for plotting difference between simulation results for MERA
    (observations) and climate model datasets for the historical period
    (1981-2005)
    """

    data = {}

    for x in ["season", "cumulative"]:
        data[f"MERA_{x[0]}"] = xr.open_mfdataset(
            glob.glob(
                os.path.join(
                    "data", "ModVege", "stats", f"*MERA*{stat}_{x}.nc"
                )
            ),
            decode_coords="all", chunks="auto"
        )

        # reassign projection
        data[f"MERA_{x[0]}"].rio.write_crs(
            cplt.lambert_conformal, inplace=True
        )

        data[f"{dataset}_{x[0]}"] = xr.open_mfdataset(
            glob.glob(
                os.path.join(
                    "data", "ModVege", "stats",
                    f"*{dataset}*{stat}_{x}_MERA.nc"
                )
            ),
            decode_coords="all", chunks="auto"
        )

        data[f"{dataset}_{x[0]}"] = data[f"{dataset}_{x[0]}"].isel(exp=0)

        # regrid climate model data
        data[f"{dataset}_{x[0]}"] = data[f"{dataset}_{x[0]}"].drop([
            "lat", "lon", "exp"
        ])
        data[f"{dataset}_{x[0]}"] = data[f"{dataset}_{x[0]}"].rename({
            "rlon": "x", "rlat": "y"
        })
        if x == "season":
            # split by season first
            for season in season_list:
                data[season] = data[f"{dataset}_{x[0]}"].sel(season=season)
                data[season] = data[season].rio.reproject_match(
                    data[f"MERA_{x[0]}"],
                    resampling=rio.enums.Resampling.bilinear
                )
                data[season] = data[season].assign_coords({
                    "x": data[f"MERA_{x[0]}"]["x"],
                    "y": data[f"MERA_{x[0]}"]["y"]
                })

            # combine seasons
            data[f"{dataset}_{x[0]}"] = xr.combine_by_coords([
                data["DJF"].expand_dims(dim="season"),
                data["MAM"].expand_dims(dim="season"),
                data["JJA"].expand_dims(dim="season"),
                data["SON"].expand_dims(dim="season")
            ])
        else:
            data[f"{dataset}_{x[0]}"] = (
                data[f"{dataset}_{x[0]}"].rio.reproject_match(
                    data[f"MERA_{x[0]}"],
                    resampling=rio.enums.Resampling.bilinear
                )
            )
            data[f"{dataset}_{x[0]}"] = (
                data[f"{dataset}_{x[0]}"].assign_coords({
                    "x": data[f"MERA_{x[0]}"]["x"],
                    "y": data[f"MERA_{x[0]}"]["y"]
                })
            )

        # clip to Ireland's boundary
        data[f"{dataset}_{x[0]}"] = data[f"{dataset}_{x[0]}"].rio.clip(
            ie.buffer(1).to_crs(cplt.lambert_conformal), all_touched=True
        )

        # calculate difference
        data[f"MERA_{x[0]}_diff"] = np.subtract(
            data[f"{dataset}_{x[0]}"], data[f"MERA_{x[0]}"]
        )

        # percentage difference (issues w/ values close to zero)
        data[f"MERA_{x[0]}_diff_pct"] = np.divide(
            data[f"MERA_{x[0]}_diff"], data[f"MERA_{x[0]}"]
        )

        # reassign attributes
        for var in data[f"MERA_{x[0]}_diff"].data_vars:
            data[f"MERA_{x[0]}_diff"][var].attrs = (
                data[f"MERA_{x[0]}"][var].attrs
            )
            data[f"MERA_{x[0]}_diff_pct"][var].attrs = (
                data[f"MERA_{x[0]}"][var].attrs
            )
            data[f"MERA_{x[0]}_diff_pct"][var].attrs["units"] = "%"
        data[f"MERA_{x[0]}_diff"].attrs["dataset"] = dataset
        data[f"MERA_{x[0]}_diff_pct"].attrs["dataset"] = dataset

        if x == "season":
            # sort seasons in the correct order
            data[f"MERA_{x[0]}_diff"] = data[f"MERA_{x[0]}_diff"].reindex(
                season=season_list
            )
            data[f"MERA_{x[0]}_diff_pct"] = (
                data[f"MERA_{x[0]}_diff_pct"].reindex(
                    season=season_list
                )
            )

    return data


def plot_obs_diff_all(data, var, season, levels=None, ticks=None):
    """
    Helper function to plot facet maps
    """

    cbar_kwargs = {
        "label": f"Difference [{data[var].attrs['units']}]",
        "aspect": 30,
        "location": "bottom",
        "fraction": 0.085,
        "shrink": 0.85,
        "pad": 0.05,
        "extendfrac": "auto"
    }

    if ticks is not None:
        cbar_kwargs["ticks"] = ticks

    if season is not None:
        data = data.sel(season=season)

    fig = data[var].plot.contourf(
        x="x", y="y", col="model",
        cmap="BrBG",
        extend="both",
        robust=True,
        cbar_kwargs=cbar_kwargs,
        transform=cplt.lambert_conformal,
        subplot_kws={"projection": cplt.plot_projection},
        levels=levels,
        xlim=(-1.775, 1.6),
        ylim=(-2.1, 2.1),
        figsize=(12, 6.25)
    )

    fig.set_titles("{value}", weight="semibold", fontsize=14)

    title_print = (
        "_" * 80 +
        f"\nDifference in average {data[var].attrs['long_name'].lower()} "
        f"from MÉRA-driven simulations for {data.attrs['dataset']}"
    )
    if season is not None:
        print(title_print + f"\n{season}, 1981-2005\n" + "_" * 80)
    else:
        print(title_print + "\n1981-2005\n" + "_" * 80)
    # plt.suptitle(title_print, size=16, y=0.8)

    # add boundary
    for axis, model in zip(fig.axs.flat, model_list):
        try:
            ie_bbox.to_crs(cplt.plot_projection).plot(
                ax=axis, edgecolor="darkslategrey", color="white",
                linewidth=.5
            )
        except NameError:
            axis.coastlines(
                resolution="10m", color="darkslategrey", linewidth=.5
            )

        # add max/min values as annotations
        da = data.sel(model=model)[var]
        # clip to remove cells over the coastline
        da = da.rio.clip(ie.buffer(1).to_crs(data.rio.crs))
        da_max = da.isel(da.compute().argmax(dim=["x", "y"]))
        da_min = da.isel(da.compute().argmin(dim=["x", "y"]))
        da = gpd.GeoDataFrame(
            {
                "name": ["max", "min"],
                "value": [da_max.values, da_min.values],
                "geometry": gpd.GeoSeries.from_wkt([
                    f"POINT ({da_max['x'].values} {da_max['y'].values})",
                    f"POINT ({da_min['x'].values} {da_min['y'].values})"
                ])
            },
            crs=cplt.lambert_conformal
        )

        da = da.to_crs(cplt.plot_projection)
        da_max = da.loc[[0]]
        da_min = da.loc[[1]]
        da_max.plot(ax=axis, color="#01665e", marker="o", edgecolor="white")
        da_min.plot(ax=axis, color="#8c510a", marker="o", edgecolor="white")
        da_max["value"] = da_max["value"].astype(float).round(1)
        da_min["value"] = da_min["value"].astype(float).round(1)

        for xy, lab in zip(
            zip(da_max["geometry"].x - 0.3, da_max["geometry"].y),
            da_max["value"]
        ):
            axis.annotate(
                text=lab, xy=xy, ha="center", va="center", size=12.5,
                weight="semibold", color="#01665e",
                path_effects=[
                    patheffects.withStroke(linewidth=3, foreground="white")
                ]
            )
        for xy, lab in zip(
            zip(da_min["geometry"].x + 0.3, da_min["geometry"].y),
            da_min["value"]
        ):
            axis.annotate(
                text=lab, xy=xy, ha="center", va="center", size=12.5,
                weight="semibold", color="#8c510a",
                path_effects=[
                    patheffects.withStroke(linewidth=3, foreground="white")
                ]
            )

    plt.show()
