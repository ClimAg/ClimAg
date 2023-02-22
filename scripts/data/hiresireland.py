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
import pandas as pd

DATA_DIR_BASE = os.path.join("data", "HiResIreland")

# directory to store outputs
DATA_DIR = os.path.join(DATA_DIR_BASE, "IE")
os.makedirs(DATA_DIR, exist_ok=True)

# Ireland boundary
GPKG_BOUNDARY = os.path.join("data", "boundaries", "boundaries.gpkg")
ie = gpd.read_file(GPKG_BOUNDARY, layer="NUTS_RG_01M_2021_2157_IE")

for exp, model in itertools.product(
    ["historical", "rcp45", "rcp85"],
    ["CNRM-CM5", "EC-EARTH", "HadGEM2-ES", "MPI-ESM-LR"]
):
    # auto-rechunking may cause NotImplementedError with object dtype
    # where it will not be able to estimate the size in bytes of object data
    if model == "HadGEM2-ES":
        chunks = 300
    else:
        chunks = "auto"

    data = xr.open_mfdataset(
        list(itertools.chain(*list(
            glob.glob(os.path.join(
                DATA_DIR_BASE, "COSMO5-CLM", exp, model, e
            ))
            for e in [
                "*mean_T_2M*.nc", "*ASOB_S*.nc", "*ET*.nc", "*TOT_PREC*.nc"
            ]
        ))),
        chunks=chunks,
        decode_coords="all"
    )

    # copy CRS
    data_crs = data.rio.crs

    # adjustments for 360-day calendar
    if model == "HadGEM2-ES":
        data = data.convert_calendar("standard", align_on="year")

    # copy time_bnds
    # data_time_bnds = data.coords["time_bnds"]

    # ### Ireland subset
    # clip to Ireland's boundary
    data = data.rio.clip(ie.buffer(1).to_crs(data_crs))

    # reassign time_bnds
    # data.coords["time_bnds"] = data_time_bnds

    # ### Calculate downward shortwave radiation
    # assume an albedo of 0.23 for grass (Allen et al., 1998)
    data = data.assign(RS=data["ASOB_S"] / (1 - 0.23))

    # ### Calculate photosynthetically active radiation
    # Papaioannou et al. (1993) - irradiance ratio
    data = data.assign(PAR=data["RS"] * 0.473)

    # keep only required variables
    data = data.drop_vars(["RS", "ASOB_S"])

    # ### Convert units and rename variables
    for v in data.data_vars:
        var_attrs = data[v].attrs  # extract attributes
        if v == "T_2M":
            var_attrs["units"] = "°C"
            data[v] = data[v] - 273.15
            var_attrs["note"] = "Converted from K to °C by subtracting 273.15"
            var_attrs["long_name"] = "Near-Surface Air Temperature"
        elif v == "PAR":
            var_attrs["units"] = "MJ m⁻² day⁻¹"
            data[v] = data[v] * (60 * 60 * 24 / 1e6)
            var_attrs["long_name"] = (
                "Surface Photosynthetically Active Radiation"
            )
            var_attrs["note"] = (
                "Calculated by dividing the surface net downward "
                "shortwave radiation with (1 - 0.23), where 0.23 is"
                " assumed to be the albedo of grass (Allen et al., 1998),"
                " then multiplying it with an irradiance ratio of 0.473"
                " based on Papaioannou et al. (1993); converted from "
                "W m⁻² to MJ m⁻² day⁻¹ by multiplying 0.0864 as "
                "documented in the FAO Irrigation and Drainage Paper No. "
                "56 (Allen et al., 1998, p. 45)"
            )
        elif v in ("TOT_PREC", "w"):
            var_attrs["units"] = "mm day⁻¹"
            var_attrs["note"] = (
                "kg m⁻² is equivalent to mm day⁻¹, assuming a water density "
                "of 1,000 kg m⁻³"
            )
            if v == "w":
                var_attrs["long_name"] = "Potential Evapotranspiration"
            else:
                var_attrs["long_name"] = "Precipitation"
        data[v].attrs = var_attrs  # reassign attributes

    # rename
    data = data.rename({"T_2M": "T", "TOT_PREC": "PP", "w": "PET"})

    # remove dataset history
    del data.attrs["history"]

    # assign dataset name
    data.attrs["dataset"] = f"IE_HiResIreland_{data.attrs['title'][:-4]}"

    # ### Extend data to a spin-up year
    data_interp = data.interp(
        time=pd.date_range(
            f"{int(data['time'][0].dt.year) - 1}-01-01T10:30:00",
            f"{int(data['time'][0].dt.year) - 1}-12-31T10:30:00",
            freq="D"
        ),
        kwargs={"fill_value": None}
    )

    data_interp.rio.write_crs(data_crs, inplace=True)

    # merge spin-up year with first two years of the main data
    data_interp = xr.combine_by_coords([
        data_interp,
        data.sel(
            time=slice(
                str(int(data["time"][0].dt.year)),
                str(int(data["time"][0].dt.year) + 1)
            )
        )
    ])

    # shift first year of the main data to the spin-up year
    data_interp = data_interp.shift(
        time=-data_interp.sel(
            time=str(int(data_interp["time"][0].dt.year))
        ).dims["time"]
    )

    # keep only spin-up year
    data_interp = data_interp.sel(
        time=str(int(data_interp["time"][0].dt.year))
    )

    # merge with main dataset
    data = xr.combine_by_coords([data, data_interp])

    # assign attributes to the data
    data.attrs["comment"] = (
        "This dataset has been clipped with the Island of Ireland's boundary"
        " and units have been converted. "
        "Last updated: " + str(datetime.now(tz=timezone.utc)) +
        " by nstreethran@ucc.ie."
    )

    # ### Export data
    # reassign CRS
    data.rio.write_crs(data_crs, inplace=True)

    # export to netCDF
    data.to_netcdf(os.path.join(DATA_DIR, f"{data.attrs['dataset']}.nc"))
