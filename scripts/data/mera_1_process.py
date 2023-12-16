"""mera_process.py

Process Met Éireann Reanalysis data

This is script 1 out of 4 to be run. This converts the GRIB files to netCDF
while keeping only the third time step of the FC3hr data. It also clips the
data to the boundary of the Island of Ireland.

https://confluence.ecmwf.int/pages/viewpage.action?pageId=197702790
"""

# import libraries
import bz2
import glob
import os
import shutil
import sys
from datetime import datetime, timezone
from itertools import chain
import geopandas as gpd
import xarray as xr
import climag.plot_configs as cplt

print("Begin MÉRA data processing...", datetime.now(tz=timezone.utc))

# Ireland boundary (derived from NUTS 2021)
GPKG_BOUNDARY = os.path.join("data", "boundaries", "boundaries.gpkg")
ie = gpd.read_file(GPKG_BOUNDARY, layer="NUTS_RG_01M_2021_2157_IE")

# directory of MÉRA GRIB files
DATA_DIR = os.path.join("/run/media/nms/Elements", "MERA", "grib")
NC_DIR = os.path.join("/run/media/nms/MyPassport", "MERA", "netcdf")

# list of folders containing variables
var_dirs = [
    "1_105_0_0",    # surface pressure
    "11_105_2_0",   # 2 m temperature
    "15_105_2_2",   # max temperature
    "16_105_2_2",   # min temperature
    "33_105_10_0",  # u-component of 10 m wind
    "34_105_10_0",  # v-component of 10 m wind
    "52_105_2_0",   # 2 m relative humidity
    "61_105_0_4",   # total precipitation
    "111_105_0_4",  # net shortwave irradiance
    "112_105_0_4",  # net longwave irradiance
    "117_105_0_4"   # global irradiance
]

for var in var_dirs:
    # directory to store temporary files
    TEMP_DIR = os.path.join(DATA_DIR, f"{var}_FC3hr", "temp")
    os.makedirs(TEMP_DIR, exist_ok=True)
    # directory to store netCDFs
    os.makedirs(os.path.join(NC_DIR, f"{var}_FC3hr"), exist_ok=True)

    # list of files; 1981 to 2020
    file_list = list(chain(*list(
        glob.glob(os.path.join(DATA_DIR, f"{var}_FC3hr", e))
        for e in [f"*{i}*{var}_FC3hr*" for i in range(1981, 2020)]
    )))

    for f in file_list:
        # if the GRIB file is archived, decompress it
        if f.endswith(".bz2"):
            # GRIB file path
            file_name = os.path.split(f)[1][:-4]
            gf = os.path.join(TEMP_DIR, file_name)
            # extract file
            with bz2.BZ2File(f, "rb") as in_file:
                with open(gf, "wb") as out_file:
                    shutil.copyfileobj(in_file, out_file)
        else:
            # if not compressed, copy GRIB file to temp dir
            os.system(f"cp {f} {TEMP_DIR}")
            file_name = os.path.split(f)[1]
            gf = os.path.join(TEMP_DIR, file_name)

        # open the GRIB file to get data length and variable attributes
        data = xr.open_dataset(
            gf, decode_coords="all", chunks="auto", engine="cfgrib"
        )
        data_varname = list(data.data_vars)[0]
        data_attrs = data[data_varname].attrs

        # keep only the third forecast step, shift the time by -3 hours to get
        # the time at the start of the forecast, and convert to netCDF
        os.system(
            "cdo -s -f nc4c -shifttime,-3hour -copy "
            f"-seltimestep,3/{len(data['time']) * 3}/3 {gf} {gf}.nc"
        )

        # assign variable attributes and clip to Ireland's boundary
        data = xr.open_dataset(f"{gf}.nc", decode_coords="all", chunks="auto")
        data = data.rename({list(data.data_vars)[0]: data_varname})
        data[data_varname].attrs = data_attrs
        data = data.rio.clip(
            ie.buffer(1).to_crs(cplt.projection_lambert_conformal), all_touched=True
        )
        data.to_netcdf(
            os.path.join(NC_DIR, f"{var}_FC3hr", file_name + ".nc")
        )

        # delete intermediate files
        os.system(f"rm -r {os.path.join(TEMP_DIR, '*')}")
        print(f, "done!")

    print(f"{var} done!", datetime.now(tz=timezone.utc))

sys.exit()
