"""eurocordex_ie.py

Subset EURO-CORDEX data for the Island of Ireland

Run the following in a Python interpreter in the project's directory and Conda
environment:

import os
exec(
    open(
        os.path.join("scripts", "data", "eurocordex_ie.py"),
        encoding="utf-8"
    ).read()
)
"""

# import libraries
import itertools
import os
from datetime import datetime, timezone
import geopandas as gpd
import intake
import xarray as xr

DATA_DIR_BASE = os.path.join("data", "EURO-CORDEX")

# directory to store outputs
DATA_DIR = os.path.join(DATA_DIR_BASE, "IE")
os.makedirs(DATA_DIR, exist_ok=True)

# Ireland boundary
GPKG_BOUNDARY = os.path.join("data", "boundaries", "boundaries.gpkg")
ie = gpd.read_file(GPKG_BOUNDARY, layer="NUTS_RG_01M_2021_2157_IE")

# reading the local catalogue
JSON_FILE_PATH = os.path.join(
    DATA_DIR_BASE, "eurocordex_eur11_local_disk.json"
)

cordex_eur11_cat = intake.open_esm_datastore(JSON_FILE_PATH)

# subset data for each experiment and driving model
driving_model = list(cordex_eur11_cat.df["driving_model"].unique())

experiment_id = list(cordex_eur11_cat.df["experiment_id"].unique())

for exp, model in itertools.product(experiment_id, driving_model):
    cordex_eur11 = cordex_eur11_cat.search(
        experiment_id=exp,
        driving_model=model
    )

    data = xr.open_mfdataset(
        list(cordex_eur11.df["uri"]),
        chunks="auto",
        decode_coords="all"
    )
    # auto-rechunking may cause NotImplementedError with object dtype
    # where it will not be able to estimate the size in bytes of object data

    # copy time_bnds coordinates
    data_time_bnds = data.coords["time_bnds"]

    # copy CRS
    data_crs = data.rio.crs

    # subset for reference period and spin-up year
    if exp == "historical":
        data = data.sel(time=slice("1975", "2005"))
    else:
        data = data.sel(time=slice("2040", "2070"))

    # clip to Ireland's boundary
    data = data.rio.clip(ie.buffer(500).to_crs(data.rio.crs))

    # reassign time_bnds
    data.coords["time_bnds"] = data_time_bnds

    # calculate photosynthetically active radiation (PAR)
    # Papaioannou et al. (1993) - irradiance ratio
    data = data.assign(PAR=data["rsds"] * 0.473)

    # convert variable units and assign attributes
    for v in data.data_vars:
        var_attrs = data[v].attrs  # extract attributes
        if v == "tas":
            var_attrs["units"] = "°C"
            data[v] = data[v] - 273.15
            var_attrs["note"] = "Converted from K to °C by subtracting 273.15"
        elif v in ("PAR", "rsds"):
            var_attrs["units"] = "MJ m⁻² day⁻¹"
            data[v] = data[v] * (60 * 60 * 24 / 1e6)
            if v == "PAR":
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
        elif v in ("pr", "evspsblpot"):
            var_attrs["units"] = "mm day⁻¹"
            data[v] = data[v] * 60 * 60 * 24
            var_attrs["note"] = (
                "Converted from kg m⁻² s⁻¹ to mm day⁻¹ by multiplying 86,400,"
                " assuming a water density of 1,000 kg m⁻³"
            )
        data[v].attrs = var_attrs  # reassign attributes

    # rename variables
    data = data.rename({"tas": "T", "pr": "PP", "evspsblpot": "PET"})

    # assign dataset name
    data.attrs["dataset"] = f"IE_EURO-CORDEX_RCA4_{model}_{exp}"

    # keep only relevant variables
    data = data.drop_vars(["rsds"])

    # assign attributes to the data
    data.attrs["comment"] = (
        "This dataset has been clipped with the Island of Ireland's boundary "
        "and units have been converted. "
        "Last updated: " + str(datetime.now(tz=timezone.utc)) +
        " by nstreethran@ucc.ie."
    )

    data.rio.write_crs(data_crs, inplace=True)

    # export to netCDF
    data.to_netcdf(os.path.join(DATA_DIR, f"{data.attrs['dataset']}.nc"))
