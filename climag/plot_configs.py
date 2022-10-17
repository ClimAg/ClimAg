"""plot_configs.py

Functions to plot climate model datasets, e.g. CORDEX
"""

from datetime import datetime
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from dateutil.parser import parse


# convert lat/lon to rotated pole coordinates
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


def rotated_pole_transform_base(pole_longitude, pole_latitude):
    """
    Rotated pole transform for a given pole longitude and latitude.
    """
    transform = ccrs.RotatedPole(
        pole_longitude=pole_longitude, pole_latitude=pole_latitude
    )
    return transform


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
    transform = rotated_pole_transform_base(pole_longitude, pole_latitude)
    return transform


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


def cordex_plot_title_main(data):
    """
    Define the map plot title for CORDEX data.

    Parameters
    ----------
    data : input CORDEX data

    Returns
    -------
    - plot title
    """
    plot_title = (
        data.attrs["project_id"] + ", " +
        data.attrs["CORDEX_domain"] + ", " +
        data.attrs["driving_model_id"] + ", " +
        data.attrs["driving_model_ensemble_member"] + ", " +
        data.attrs["driving_experiment_name"] + ", " +
        data.attrs["model_id"] + ", " +
        data.attrs["rcm_version_id"] + ", " +
        data.attrs["frequency"]
    )
    return plot_title


def hiresireland_date_format(data):
    """
    Format date
    """
    date = datetime.strftime(parse(str(data["time"].values)), "%-d %b %Y")
    return date


def cordex_date_format(data):
    """
    Format date
    """
    if data.attrs["frequency"] == "mon":
        date_format = "%b %Y"
    elif data.attrs["frequency"] == "day":
        date_format = "%-d %b %Y"
    else:
        date_format = "%Y-%m-%d %H:%M:%S"
    date = datetime.strftime(parse(str(data["time"].values)), date_format)
    return date


def cordex_plot_title(data, lon=None, lat=None):
    """
    Define the map plot title for CORDEX data with information about the time
    or point subset.

    Parameters
    ----------
    data : input CORDEX data
    lon : longitude of the point subset (optional)
    lat : latitude of the point subset (optional)

    Returns
    -------
    - plot title
    """
    if lon is None and lat is None:
        end_str = cordex_date_format(data)
    else:
        end_str = "(" + str(lon) + ", " + str(lat) + ")"
    plot_title = cordex_plot_title_main(data) + ", " + end_str
    return plot_title


def ie_cordex_ncfile_name(data):
    """
    Define the NetCDF file name for the CORDEX data that has been subset for
    Ireland.

    Parameters
    ----------
    data : input CORDEX data

    Returns
    -------
    - file name
    """
    filename = (
        "_".join(list(data.data_vars)) + "_" +
        data.attrs["CORDEX_domain"] + "_" +
        data.attrs["driving_model_id"] + "_" +
        data.attrs["driving_experiment_name"] + "_" +
        data.attrs["driving_model_ensemble_member"] + "_" +
        data.attrs["model_id"] + "_" +
        data.attrs["rcm_version_id"] + "_" +
        data.attrs["frequency"] + "_" +
        datetime.strftime(parse(str(data["time"][0].values)), "%Y%m%d") + "-" +
        datetime.strftime(parse(str(data["time"][-1].values)), "%Y%m%d")
        + "_IE.nc"
    )
    return filename


def ie_cordex_modvege_ncfile_name(cordex_data, output_data):
    """
    Define the NetCDF file name for the ModVege grass growth model outputs.

    Parameters
    ----------
    data : input CORDEX data

    Returns
    -------
    - file name
    """
    filename = (
        "modvege_" +
        cordex_data.attrs["CORDEX_domain"] + "_" +
        cordex_data.attrs["driving_model_id"] + "_" +
        cordex_data.attrs["driving_experiment_name"] + "_" +
        cordex_data.attrs["driving_model_ensemble_member"] + "_" +
        cordex_data.attrs["model_id"] + "_" +
        cordex_data.attrs["rcm_version_id"] + "_" +
        cordex_data.attrs["frequency"] + "_" +
        datetime.strftime(parse(str(output_data["time"][0].values)), "%Y%m%d")
        + "-" +
        datetime.strftime(parse(str(output_data["time"][-1].values)), "%Y%m%d")
        + "_IE.nc"
    )
    return filename


def plot_facet_map_variables(data, boundary_data):
    """
    Create a facet plot of variables from an xarray dataset covering the
    Island of Ireland.

    Parameters
    ----------
    data : dataset (loaded using xarray)
    boundary_data : Ireland boundary data (vector, loaded as a GeoPandas
        dataframe), e.g. from Ordnance Survey Ireland, NUTS (Eurostat)
    """
    for v in data.data_vars:
        cbar_label = (
            data[v].attrs["long_name"] + " [" + data[v].attrs["units"] + "]"
        )  # colorbar label

        if v == "pr":
            cmap = "mako_r"
        elif v == "evspsblpot":
            cmap = "BrBG_r"
        elif v in ("tas", "rsun"):
            cmap = "Spectral_r"
        else:
            cmap = "YlGn"

        fig = data[v].plot(
            x="lon", y="lat", col="time", col_wrap=5, cmap=cmap, levels=15,
            robust=True, cbar_kwargs=dict(aspect=40, label=cbar_label)
        )

        fig.set_xlabels("")
        fig.set_ylabels("")

        for i, ax in enumerate(fig.axes.flat):
            boundary_data.to_crs(4326).boundary.plot(
                ax=ax, color="darkslategrey", linewidth=.5
            )
            ax.set_title(cordex_date_format(data.isel(time=i)))
            ax.xaxis.set_major_formatter(longitude_tick_format)
            ax.yaxis.set_major_formatter(latitude_tick_format)

        plt.show()


def plot_map_variables(data):
    """
    Create individual plots of the climate data variables covering the Island
    of Ireland.

    Parameters
    ----------
    data : climate model dataset (loaded using xarray)
    """
    for v in data.data_vars:
        cbar_label = (
            data[v].attrs["long_name"] + " [" + data[v].attrs["units"] + "]"
        )  # colorbar label
        if v == "pr":
            cmap = "GnBu"
        elif v == "evspsblpot":
            cmap = "BrBG_r"
        else:
            cmap = "Spectral_r"

        plt.figure(figsize=(7.5, 7))

        plot_transform = rotated_pole_transform(data)
        ax = plt.axes(projection=plot_transform)

        # specify gridline spacing and labels
        ax.gridlines(
            draw_labels=True,
            xlocs=range(-180, 180, 2),
            ylocs=range(-90, 90, 1),
            color="lightslategrey",
            linewidth=.5
        )

        # plot data for the variable
        data[v].plot(
            ax=ax,
            cmap=cmap,
            x="rlon",
            y="rlat",
            levels=15,
            cbar_kwargs=dict(label=cbar_label),
            robust=True
        )

        # add boundaries
        ax.coastlines(resolution="10m", color="darkslategrey", linewidth=.75)

        # ax.set_title(cplt.cordex_plot_title(data_ie))  # set plot title
        ax.set_title(None)

        plt.axis("equal")
        plt.tight_layout()
        plt.show()


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
