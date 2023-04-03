"""plot_stats.py

Helper functions to plot statistics
"""

import glob
import os
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
# import rasterio as rio
import xarray as xr
import climag.plot_configs as cplt

# Ireland boundary
ie_bbox = gpd.read_file(
    os.path.join("data", "boundaries", "boundaries.gpkg"),
    layer="ne_10m_land_2157_IE_BBOX_DIFF"
)


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


# def hist_obs_diff(stat, diff=True):
#     """
#     Difference between historical and observational data
#     """

#     data = {}
#     temp = {}

#     for x in ["season", "cumulative"]:
#         data[f"MERA_{x[0]}"] = xr.open_mfdataset(
#             glob.glob(
#                 os.path.join(
#                     "data", "ModVege", "stats", f"IE_MERA*{stat}_{x}.nc"
#                 )
#             ),
#             decode_coords="all", chunks="auto"
#         )

#         # reassign projection
#         data[f"MERA_{x[0]}"].rio.write_crs(
#             cplt.lambert_conformal, inplace=True
#         )

#         if diff:
#             for dataset in ["EURO-CORDEX", "HiResIreland"]:
#                 temp[f"{dataset}_{x[0]}"] = xr.open_mfdataset(
#                     glob.glob(
#                         os.path.join(
#                             "data", "ModVege", "stats",
#                             f"IE_{dataset}*historical_{stat}_{x}_MERA.nc"
#                         )
#                     ),
#                     decode_coords="all", chunks="auto"
#                 )
#                 temp[f"{dataset}_{x[0]}"] = (
#                     temp[f"{dataset}_{x[0]}"].isel(exp=0)
#                 )

#                 if x == "season":
#                     # regrid climate model data to match observations
#                     for season in ["DJF", "MAM", "JJA", "SON"]:
#                         temp[season] = temp[f"{dataset}_{x[0]}"].drop([
#                             "lat", "lon", "exp"
#                         ])
#                         temp[season] = temp[season].rename({
#                             "rlon": "x", "rlat": "y"
#                         })
#                         temp[season] = temp[season].sel(season=season)
#                         temp[season] = temp[season].rio.reproject_match(
#                             data[f"MERA_{x[0]}"],
#                             resampling=rio.enums.Resampling.bilinear
#                         )
#                         temp[season] = temp[season].assign_coords({
#                             "x": data[f"MERA_{x[0]}"]["x"],
#                             "y": data[f"MERA_{x[0]}"]["y"]
#                         })

#                     temp[f"{dataset}_{x[0]}"] = xr.combine_by_coords([
#                         temp["DJF"].expand_dims(dim="season"),
#                         temp["MAM"].expand_dims(dim="season"),
#                         temp["JJA"].expand_dims(dim="season"),
#                         temp["SON"].expand_dims(dim="season")
#                     ])

#             # combine
#             temp[f"MERA_{x[0]}"] = xr.combine_by_coords(
#                 [
#                     (
#                         temp[f"EURO-CORDEX_{x[0]}"] - data[f"MERA_{x[0]}"]
#                     ).assign_coords(exp="EURO-CORDEX - MERA").expand_dims(
#                         dim="exp"
#                     ),
#                     (
#                         temp[f"HiResIreland_{x[0]}"] - data[f"MERA_{x[0]}"]
#                     ).assign_coords(exp="HiResIreland - MERA").expand_dims(
#                         dim="exp"
#                     )
#                 ],
#                 combine_attrs="drop_conflicts"
#             )

#             # reassign variable attributes
#             for var in temp[f"MERA_{x[0]}"].data_vars:
#                 temp[f"MERA_{x[0]}"][var].attrs = (
#                     data[f"MERA_{x[0]}"][var].attrs
#                 )

#             # sort seasons in the correct order
#             data[f"MERA_{x[0]}"] = temp[f"MERA_{x[0]}"].reindex(
#                 season=["DJF", "MAM", "JJA", "SON"]
#             )

#     return data


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
        cmap = "RdBu"
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


# def mera_climate_stats_data():
#     """
#     Prepare data for plotting difference between simulation results for MERA
#     (observations) and climate model datasets for the historical period
#     (1981-2005)
#     """
