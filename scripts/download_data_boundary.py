# import libraries
import os
import geopandas as gpd
from src import download_data as dd

# base data download directory
DATA_DIR = os.path.join("data", "boundary")

GPKG_BOUNDARY = os.path.join("data", "boundary", "boundaries.gpkg")

######################################################################
# Administrative Areas - OSi National Statutory Boundaries - 2019

URL = (
    "https://data-osi.opendata.arcgis.com/datasets/" +
    "d81188d16e804bde81548e982e80c53e_0.geojson"
)

payload = {
    "outSR": {
        "latestWkid": "2157",
        "wkid": "2157"
    }
}

SUB_DIR = os.path.join(DATA_DIR, "admin-osi", "raw")

dd.download_data(server=URL, ddir=SUB_DIR, params=payload)

DATA_FILE = os.path.join(SUB_DIR, "data.geojson")

osi = gpd.read_file(DATA_FILE)

osi.to_file(GPKG_BOUNDARY, layer="Admin_Areas_OSi")

######################################################################
# OSNI Open Data - Largescale Boundaries - County Boundaries

URL = (
    "https://osni-spatialni.opendata.arcgis.com/datasets/spatialni::" +
    "osni-open-data-largescale-boundaries-county-boundaries-.geojson"
)

payload = {
    "outSR": {
        "latestWkid": "29902",
        "wkid": "29900"
    }
}

SUB_DIR = os.path.join(DATA_DIR, "admin-osni", "raw")

dd.download_data(server=URL, ddir=SUB_DIR, params=payload)

DATA_FILE = os.path.join(SUB_DIR, "data.geojson")

osni = gpd.read_file(DATA_FILE)

osni.to_file(GPKG_BOUNDARY, layer="Admin_Areas_OSNI")
