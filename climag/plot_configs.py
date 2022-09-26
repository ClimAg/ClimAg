"""plot_configs.py

Functions to plot climate model datasets, e.g. CORDEX
"""

from datetime import datetime
import cartopy.crs as ccrs
from dateutil.parser import parse


# convert lat/lon to rotated pole coordinates
def rotated_pole_point(data, lon, lat):
    """
    Convert the latitude and longitude of a specific point to rotated pole
    coordinates used in the input data.

    Parameters
    ----------
    data : input climate data which uses rotated pole coordinates
    lon : longitude of the point
    lat : latitude of the point
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


def cordex_plot_title(data, lon=None, lat=None):
    """
    Define the map plot title for CORDEX data.

    Parameters
    ----------
    data : input CORDEX data
    """
    if lon is None and lat is None:
        if data.attrs["frequency"] == "mon":
            date_format = "%b %Y"
        elif data.attrs["frequency"] == "day":
            date_format = "%-d %b %Y"
        else:
            date_format = "%Y-%m-%d %H:%M:%S"
        end_str = datetime.strftime(
            parse(str(data["time"].values)), date_format
        )
    else:
        end_str = "(" + str(lon) + ", " + str(lat) + ")"
    plot_title = (
        data.attrs["project_id"] + ", " +
        data.attrs["CORDEX_domain"] + ", " +
        data.attrs["driving_model_id"] + ", " +
        data.attrs["driving_model_ensemble_member"] + ", " +
        data.attrs["driving_experiment_name"] + ", " +
        data.attrs["model_id"] + ", " +
        data.attrs["rcm_version_id"] + ", " +
        data.attrs["frequency"] + ", " +
        end_str
    )
    return plot_title


def rotated_pole_transform(data):
    """
    Rotated pole transform for plotting CORDEX data.

    Parameters
    ----------
    data : input CORDEX data
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


# def data_plot(
#     data,
#     cmap="terrain",
#     vmin=None,
#     vmax=None,
#     grid_color="lightslategrey",
#     border_color="darkslategrey",
#     border_width=.5,
#     border_res="50m",
#     cbar_label=None,
#     transform=None,
#     grid_xlocs=range(-180, 180, 10),
#     grid_ylocs=range(-90, 90, 5),
#     plot_title=None,
#     plot_figsize=(20, 10)
# ):
#     """
#     Custom plot function
#     """
#     plt.figure(figsize=plot_figsize)
#     ax = plt.axes(projection=transform)
#     ax.gridlines(
#         draw_labels=True,
#         linewidth=.5,
#         color=grid_color,
#         xlocs=grid_xlocs,
#         ylocs=grid_ylocs
#     )
#     data.plot(
#         ax=ax,
#         cmap=cmap,
#         transform=transform,
#         vmin=vmin,
#         vmax=vmax,
#         x="rlon",
#         y="rlat",
#         cbar_kwargs=dict(label=cbar_label)
#     )
#     ax.coastlines(
#         resolution=border_res, color=border_color, linewidth=border_width
#     )
#     if plot_title is not None:
#         ax.set_title(plot_title)
