"""hiresireland.py

Subset HiResIreland data for the Island of Ireland

Run the following in a Python interpreter in the project's directory and Conda
environment:

import os
exec(
    open(
        os.path.join("scripts", "data", "hiresireland.py"),
        encoding="utf-8"
    ).read()
)
"""

# import libraries
import glob
import itertools
import os
from datetime import datetime, timezone
import geopandas as gpd
import xarray as xr

DATA_DIR_BASE = os.path.join("data", "HiResIreland")

# directory to store outputs
DATA_DIR = os.path.join(DATA_DIR_BASE, "IE")
os.makedirs(DATA_DIR, exist_ok=True)

# Ireland boundary
GPKG_BOUNDARY = os.path.join(
    "data", "boundaries", "NUTS2021", "NUTS_2021.gpkg"
)
ie = gpd.read_file(GPKG_BOUNDARY, layer="NUTS_RG_01M_2021_2157_IE")

for exp, model in itertools.product(
    ["historical", "rcp45", "rcp85"],
    ["CNRM-CM5", "EC-EARTH", "HadGEM2-ES", "MPI-ESM-LR"]
):
    data = xr.open_mfdataset(
        list(itertools.chain(*list(
            glob.glob(os.path.join(
                DATA_DIR_BASE, "COSMO5-CLM", exp, model, e
            ))
            for e in [
                "*mean_T_2M*.nc", "*ASOB_S*.nc", "*ET*.nc", "*TOT_PREC*.nc"
            ]
        ))),
        # chunks="auto",
        decode_coords="all"
    )
    # disable auto-rechunking; may cause NotImplementedError with object dtype
    # where it will not be able to estimate the size in bytes of object data

    # copy time_bnds
    data_time_bnds = data.coords["time_bnds"]

    # copy CRS
    data_crs = data.rio.crs

    # ### Ireland subset
    # clip to Ireland's boundary
    data = data.rio.clip(ie.buffer(1).to_crs(data_crs))

    # reassign time_bnds
    data.coords["time_bnds"] = data_time_bnds

    # ### Calculate photosynthetically active radiation
    # Papaioannou et al. (1993) - irradiance ratio
    data = data.assign(par=(data["ASOB_S"] * 0.473))

    # ### Convert units and rename variables
    for v in data.data_vars:
        var_attrs = data[v].attrs  # extract attributes
        if v == "T_2M":
            var_attrs["units"] = "°C"  # convert K to deg C
            data[v] = data[v] - 273.15
            var_attrs["note"] = (
                f"Original name is '{v}'; converted from K to °C by "
                "subtracting 273.15"
            )
            var_attrs["long_name"] = "Near-Surface Air Temperature"
        elif v in ("ASOB_S", "par"):
            var_attrs["units"] = "MJ m⁻² day⁻¹"
            # convert W m-2 to MJ m-2 day-1
            # Allen (1998) - FAO Irrigation and Drainage Paper No. 56 (p. 45)
            # (per second to per day; then convert to mega)
            data[v] = data[v] * (60 * 60 * 24 / 1e6)
            if v == "par":
                var_attrs["long_name"] = (
                    "Surface Photosynthetically Active Radiation"
                )
                var_attrs["note"] = (
                    "Calculated by multiplying 'rsds' with an irradiance "
                    "ratio of 0.473 based on Papaioannou et al. (1993); "
                    "converted from W m⁻² to MJ m⁻² day⁻¹ by multiplying "
                    "0.0864 based on the FAO Irrigation and Drainage Paper "
                    "No. 56 (Allen et al., 1998, p. 45)"
                )
            else:
                var_attrs["long_name"] = (
                    "Surface Downwelling Shortwave Radiation"
                )
                var_attrs["note"] = (
                    f"Original name is '{v}'; converted from W m⁻² to "
                    "MJ m⁻² day⁻¹ by multiplying 0.0864 based on the FAO "
                    "Irrigation and Drainage Paper No. 56 "
                    "(Allen et al., 1998, p. 45)"
                )
        elif v in ("TOT_PREC", "w"):
            var_attrs["units"] = "mm day⁻¹"  # kg m-2 is the same as mm day-1
            var_attrs["note"] = (
                f"Original name is '{v}'; kg m⁻² is equivalent to mm day⁻¹, "
                "assuming a water density of 1,000 kg m⁻³"
            )
            if v == "w":
                var_attrs["long_name"] = "Potential Evapotranspiration"
            else:
                var_attrs["long_name"] = "Precipitation"
        data[v].attrs = var_attrs  # reassign attributes

    # rename
    data = data.rename({
        "T_2M": "T", "ASOB_S": "RG", "TOT_PREC": "PP",
        "w": "PET", "par": "PAR"
    })

    # assign attributes for the data
    data.attrs["comment"] = (
        "This dataset has been clipped with the Island of Ireland's boundary "
        "and units have been converted. "
        "Last updated: " + str(datetime.now(tz=timezone.utc)) +
        " by nstreethran@ucc.ie."
    )

    # remove dataset history
    del data.attrs["history"]

    # ### Export data
    # reassign CRS
    data.rio.write_crs(data_crs, inplace=True)

    # export to NetCDF
    FILE_NAME = "IE_" + data.attrs["title"] + ".nc"
    data.to_netcdf(os.path.join(DATA_DIR, FILE_NAME))
