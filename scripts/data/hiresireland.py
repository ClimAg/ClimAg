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
import sys
from datetime import datetime, timezone

import geopandas as gpd
import pandas as pd
import xarray as xr

# from dask.distributed import Client

# client = Client(n_workers=2, threads_per_worker=4, memory_limit="3GB")

DATA_DIR_BASE = os.path.join("data", "HiResIreland")

# directory to store outputs
DATA_DIR = os.path.join(DATA_DIR_BASE, "IE")
os.makedirs(DATA_DIR, exist_ok=True)

# Ireland boundary
GPKG_BOUNDARY = os.path.join("data", "boundaries", "boundaries.gpkg")
ie = gpd.read_file(GPKG_BOUNDARY, layer="NUTS_RG_01M_2021_2157_IE")

for exp, model in itertools.product(
    ["historical", "rcp45", "rcp85"],
    ["CNRM-CM5", "EC-EARTH", "HadGEM2-ES", "MPI-ESM-LR"],
):
    # auto-rechunking may cause NotImplementedError with object dtype
    # where it will not be able to estimate the size in bytes of object data
    if model == "HadGEM2-ES":
        CHUNKS = 300
    else:
        CHUNKS = "auto"

    # rename evapotranspiration field
    et_file = os.path.join(
        DATA_DIR_BASE,
        "COSMO5-CLM",
        exp,
        model,
        f"day_sum_ET_COSMO5_{model}_{exp}_4km",
    )

    os.system(f"cdo -s chname,w,ET {et_file}.nc {et_file}_ref.nc")

    data = xr.open_mfdataset(
        list(
            itertools.chain(
                *list(
                    glob.glob(
                        os.path.join(
                            DATA_DIR_BASE, "COSMO5-CLM", exp, model, e
                        )
                    )
                    for e in [
                        f"*mean_T_2M*{model}*{exp}*.nc",
                        f"S1_daymean*{model}*{exp}*.nc",
                        f"*ET*{model}*{exp}*4km_*.nc",
                        f"*TOT_PREC*{model}*{exp}*.nc",
                    ]
                )
            )
        ),
        chunks=CHUNKS,
        decode_coords="all",
    )

    # copy CRS
    data_crs = data.rio.crs

    # adjustments for 360-day calendar
    if model == "HadGEM2-ES":
        data = data.convert_calendar("standard", align_on="year")

    # copy time_bnds
    data_time_bnds = data.coords["time_bnds"]

    # ### Ireland subset
    # clip to Ireland's boundary
    data = data.rio.clip(ie.buffer(1).to_crs(data_crs), all_touched=True)

    # ### Calculate photosynthetically active radiation
    # Papaioannou et al. (1993) - irradiance ratio
    data = data.assign(PAR=(data["ASWDIR_S"] + data["ASWDIFD_S"]) * 0.473)

    # keep only required variables
    data = data.drop_vars(["ASWDIR_S", "ASWDIFD_S", "ASWDIFU_S", "ALB_RAD"])

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
            var_attrs[
                "long_name"
            ] = "Surface Photosynthetically Active Radiation"
            var_attrs["note"] = (
                "Calculated by multiplying the surface downwelling shortwave "
                "radiation (calculated by summing the direct and diffuse "
                "downward shortwave radiation components) with an irradiance "
                "ratio of 0.473 based on Papaioannou et al. (1993); "
                "converted from W m⁻² to MJ m⁻² day⁻¹ by multiplying 0.0864 "
                "as documented in the FAO Irrigation and Drainage Paper No. "
                "56 (Allen et al., 1998, p. 45)"
            )
        elif v in ("TOT_PREC", "ET"):
            var_attrs["units"] = "mm day⁻¹"
            if v == "ET":
                var_attrs["long_name"] = "Potential Evapotranspiration"
            else:
                var_attrs["long_name"] = "Precipitation"
                var_attrs["note"] = (
                    "kg m⁻² is equivalent to mm day⁻¹, assuming a water "
                    "density of 1,000 kg m⁻³"
                )
        data[v].attrs = var_attrs  # reassign attributes

    # rename variables
    data = data.rename({"T_2M": "T", "TOT_PREC": "PP", "ET": "PET"})

    # remove dataset history
    del data.attrs["history"]

    # assign dataset name
    data.attrs["dataset"] = f"IE_HiResIreland_{data.attrs['title'][:-4]}"

    # create spin-up year for PET for the historical data
    if exp == "historical":
        data_et = data.drop_vars(["T", "PP", "PAR"]).sel(
            time=slice("1976", "2005")
        )
        data = data.drop_vars(["PET"])

        # ### Extend data to a spin-up year
        data_interp = data_et.interp(
            time=pd.date_range(
                "1975-01-01T10:30:00", "1975-12-31T10:30:00", freq="D"
            ),
            kwargs={"fill_value": None},
        )

        # adjustments for 360-day calendar
        if model == "HadGEM2-ES":
            data_interp = data_interp.convert_calendar(
                "360_day", align_on="year"
            )
            data_interp = data_interp.convert_calendar(
                "standard", align_on="year"
            )

        # merge spin-up year with first two years of PET of the main data
        data_interp = xr.combine_by_coords(
            [data_interp, data_et.sel(time=slice("1976", "1977"))]
        )

        # shift first year of the main data to the spin-up year
        data_interp = data_interp.shift(
            time=-data_interp.sel(time="1975").dims["time"]
        )

        # keep only spin-up year
        data_interp = data_interp.sel(time="1975")

        # merge with main dataset
        data_et = xr.combine_by_coords([data_et, data_interp])
        data = xr.combine_by_coords([data, data_et])

    # reassign time_bnds
    data.coords["time_bnds"] = data_time_bnds

    # assign attributes to the data
    data.attrs["comment"] = (
        "This dataset has been clipped with the Island of Ireland's boundary "
        "and units have been converted. Last updated: "
        + str(datetime.now(tz=timezone.utc))
        + " by nstreethran@ucc.ie."
    )

    # ### Export data
    # reassign CRS
    data.rio.write_crs(data_crs, inplace=True)

    # export to netCDF
    data.to_netcdf(os.path.join(DATA_DIR, f"{data.attrs['dataset']}.nc"))

    print(f"{data.attrs['dataset']} done!")

sys.exit()
