"""plot_stats.py

Helper functions to plot statistics
"""

import glob
import os
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import rasterio as rio
import xarray as xr
import climag.plot_configs as cplt

# Ireland boundary
ie_bbox = gpd.read_file(
    os.path.join("data", "boundaries", "boundaries.gpkg"),
    layer="ne_10m_land_2157_IE_BBOX_DIFF"
)
ie = gpd.read_file(
    os.path.join("data", "boundaries", "boundaries.gpkg"),
    layer="ne_10m_land_2157_IE"
)

season_list = ["DJF", "MAM", "JJA", "SON"]
model_list = ["CNRM-CM5", "EC-EARTH", "HadGEM2-ES", "MPI-ESM-LR"]


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


def hist_obs_diff(stat, dataset):
    """
    Difference between historical and observational data
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

        for model in model_list:
            # calculate difference
            data[model] = data[f"{dataset}_{x[0]}"].sel(model=model)
            data[model] = data[model] - data[f"MERA_{x[0]}"]

        # combine
        data[f"MERA_{x[0]}_diff"] = xr.combine_by_coords([
            data["CNRM-CM5"].expand_dims(dim="model"),
            data["EC-EARTH"].expand_dims(dim="model"),
            data["HadGEM2-ES"].expand_dims(dim="model"),
            data["MPI-ESM-LR"].expand_dims(dim="model")
        ])

        # # calculate difference
        # data[f"MERA_{x[0]}_diff"] = (
        #     data[f"{dataset}_{x[0]}"] - data[f"MERA_{x[0]}"]
        # )

        # reassign attributes
        for var in data[f"MERA_{x[0]}_diff"].data_vars:
            data[f"MERA_{x[0]}_diff"][var].attrs = (
                data[f"MERA_{x[0]}"][var].attrs
            )
        data[f"MERA_{x[0]}_diff"].attrs["dataset"] = dataset

        if x == "season":
            # sort seasons in the correct order
            data[f"MERA_{x[0]}_diff"] = data[f"MERA_{x[0]}_diff"].reindex(
                season=season_list
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
        "label": (
            f"{data[var].attrs['long_name']} [{data[var].attrs['units']}]"
        )
    }

    if ticks is not None:
        cbar_kwargs["ticks"] = ticks

    if len(data["exp"]) == 3:
        cmap = cplt.colormap_configs(var)
        cbar_kwargs["aspect"] = 30
        figsize = (12.45, 9.25)
        extend = "max"
        robust = False
    else:
        cmap = "BrBG"
        cbar_kwargs["aspect"] = 19
        figsize = (12.35, 6.25)
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

    # add boundary
    for axis in fig.axs.flat:
        try:
            ie_bbox.to_crs(cplt.plot_projection).plot(
                ax=axis, edgecolor="darkslategrey", color="white",
                linewidth=.5
            )
        except NameError:
            axis.coastlines(
                resolution="10m", color="darkslategrey", linewidth=.5
            )

    plt.show()


def plot_obs_diff_all(data, var, season, levels=None, ticks=None):
    """
    Helper function to plot facet maps
    """

    plot_transform = cplt.lambert_conformal

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
        transform=plot_transform,
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
    for axis in fig.axs.flat:
        try:
            ie_bbox.to_crs(cplt.plot_projection).plot(
                ax=axis, edgecolor="darkslategrey", color="white",
                linewidth=.5
            )
        except NameError:
            axis.coastlines(
                resolution="10m", color="darkslategrey", linewidth=.5
            )

    plt.show()


def colorbar_levels(maximum, num=26):
    levels = [-maximum + maximum / ((num - 1) / 2) * n for n in range(num)]
    return levels


def colorbar_ticks(maximum):
    ticks = [-maximum + maximum / 3 * n for n in range(7)]
    return ticks


# def plot_obs_diff_all_new(data, var, season, levels=None, ticks=None):

# def mera_climate_stats_data():
#     """
#     Prepare data for plotting difference between simulation results for MERA
#     (observations) and climate model datasets for the historical period
#     (1981-2005)
#     """
