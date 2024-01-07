"""Utility functions

"""

import cartopy.crs as ccrs

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
