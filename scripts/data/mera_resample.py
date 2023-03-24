"""mera_process.py

Resample Met Éireann Reanalysis data to a daily frequency

https://confluence.ecmwf.int/pages/viewpage.action?pageId=197702790
"""

# import libraries
import glob
import os
import sys
from datetime import datetime, timezone
import geopandas as gpd
import xarray as xr
import climag.plot_configs as cplt

print("Started MÉRA data processing...", datetime.now(tz=timezone.utc))

# directory of MERA netCDF files
DATA_DIR = os.path.join("/run/media/nms/MyPassport", "MERA", "netcdf")

# directory to store outputs
OUT_DIR = os.path.join("/run/media/nms/MyPassport", "MERA", "netcdf_day")
os.makedirs(OUT_DIR, exist_ok=True)

# Ireland boundary
GPKG_BOUNDARY = os.path.join("data", "boundaries", "boundaries.gpkg")
ie = gpd.read_file(GPKG_BOUNDARY, layer="NUTS_RG_01M_2021_2157_IE")

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
    data = xr.open_mfdataset(
        glob.glob(os.path.join(DATA_DIR, f"{var}_FC3hr", f"*{var}_FC3hr.nc")),
        chunks="auto", decode_coords="all"
    )

    varname = list(data.data_vars)[0]

    data_attrs = data[varname].attrs  # copy var attributes
    data_crs = data.rio.crs  # copy CRS
    time_attrs = data["time"].attrs  # copy time attributes

    # convert radiation to W m-2
    if var in ("111_105_0_4", "112_105_0_4", "117_105_0_4"):
        data[varname] = data[varname] / (3 * 3600)

    # resample using sum
    if var == "61_105_0_4":
        data = data.resample(time="D").sum()
    # resample using max
    elif var == "15_105_2_2":
        data = data.resample(time="D").max()
    # resample using min
    elif var == "16_105_2_2":
        data = data.resample(time="D").min()
    # resample using mean
    else:
        data = data.resample(time="D").mean()

    # clip to Ireland's boundary to remove NaNs after summing
    if var == "61_105_0_4":
        data.rio.write_crs(data_crs, inplace=True)  # reassign CRS
        data = data.rio.clip(
            ie.buffer(1).to_crs(cplt.lambert_conformal), all_touched=True
        )

    # convert units
    if var in ("111_105_0_4", "112_105_0_4", "117_105_0_4"):
        data[varname] = data[varname] * (60 * 60 * 24 / 1e6)  # MJ m-2 day-1
    elif var in ("11_105_2_0", "15_105_2_2", "16_105_2_2"):
        data[varname] = data[varname] - 273.15  # deg C
    elif var == "1_105_0_0":
        data[varname] = data[varname] / 1000  # kPa

    data[varname].attrs = data_attrs  # reassign attributes
    data["time"].attrs = time_attrs
    data.rio.write_crs(data_crs, inplace=True)  # reassign CRS

    # update units
    if var in ("111_105_0_4", "112_105_0_4", "117_105_0_4"):
        data[varname].attrs["units"] = "MJ m⁻² day⁻¹"
    elif var in ("11_105_2_0", "15_105_2_2", "16_105_2_2"):
        data[varname].attrs["units"] = "°C"
    elif var == "61_105_0_4":
        data[varname].attrs["units"] = "mm day⁻¹"
    elif var in ("33_105_10_0", "34_105_10_0"):
        data[varname].attrs["units"] = "m s⁻¹"
    elif var == "1_105_0_0":
        data[varname].attrs["units"] = "kPa"

    # save data
    data.to_netcdf(os.path.join(OUT_DIR, f"MERA_{var}_day.nc"))

    print(f"{var} done!", datetime.now(tz=timezone.utc))

sys.exit()
