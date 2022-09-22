"""download_data_boundary.py

Download Ireland boundary data and save as layers in a GeoPackage file
"""

# import libraries
import os
from zipfile import BadZipFile, ZipFile
import geopandas as gpd
from climag.download_data import download_data

# base data download directory
DATA_DIR = os.path.join("data", "boundary")

GPKG_BOUNDARY = os.path.join("data", "boundary", "boundaries.gpkg")

######################################################################
# Administrative Areas - OSi National Statutory Boundaries - 2019

SUB_DIR = os.path.join(DATA_DIR, "admin-osi", "raw")

# download data if necessary
URL = (
    "https://data-osi.opendata.arcgis.com/datasets/"
    "osi::counties-osi-national-statutory-boundaries-2019.zip"
)

payload = {
    "outSR": {
        "latestWkid": "2157",
        "wkid": "2157"
    }
}

download_data(server=URL, dl_dir=SUB_DIR, params=payload)

ZIP_FILE = os.path.join(
    SUB_DIR, "Counties_-_OSi_National_Statutory_Boundaries_-_2019.zip"
)

osi = gpd.read_file(
    f"zip://{ZIP_FILE}!Counties___OSi_National_Statutory_Boundaries_.shp"
)

osi.to_file(GPKG_BOUNDARY, layer="OSi_Counties")

######################################################################
# OSNI Open Data - Largescale Boundaries - County Boundaries

SUB_DIR = os.path.join(DATA_DIR, "admin-osni", "raw")

# download data if necessary
URL = (
    "https://osni-spatialni.opendata.arcgis.com/datasets/spatialni::"
    "osni-open-data-largescale-boundaries-county-boundaries-.zip"
)

payload = {
    "outSR": {
        "latestWkid": "29902",
        "wkid": "29900"
    }
}

download_data(server=URL, dl_dir=SUB_DIR, params=payload)

ZIP_FILE = os.path.join(
    SUB_DIR, "OSNI_Open_Data_-_Largescale_Boundaries_-_County_Boundaries_.zip"
)

osni = gpd.read_file(
    f"zip://{ZIP_FILE}!OSNI_Open_Data_-_Largescale_Boundaries_-_"
    "County_Boundaries_.shp"
)

osni.to_file(GPKG_BOUNDARY, layer="OSNI_Counties")

######################################################################
# OS County boundaries - Island of Ireland

osi_counties = osi[["CONTAE", "COUNTY", "PROVINCE", "geometry"]]

osni_counties = osni.rename(columns={"CountyName": "COUNTY"})

osni_counties = osni_counties[["geometry", "COUNTY"]]

# https://en.wikipedia.org/wiki/Counties_of_Ireland
contae = {
    "ANTRIM": "Aontroim",
    "ARMAGH": "Ard Mhacha",
    "DOWN": "An Dún",
    "FERMANAGH": "Fear Manach",
    "LONDONDERRY": "Doire",
    "TYRONE": "Tír Eoghain"
}

osni_counties["CONTAE"] = osni_counties["COUNTY"].map(contae)

osni_counties["PROVINCE"] = "Ulster"

# reproject to Irish Transverse Mercator
osi_counties = osi_counties.to_crs(2157)

osni_counties = osni_counties.to_crs(2157)

# remove overlapping areas in OSi layer
osi_counties = osi_counties.overlay(osni_counties, how="difference")

# merge county layers
ie_counties = osi_counties.merge(osni_counties, how="outer")

ie_counties.to_file(GPKG_BOUNDARY, layer="OS_IE_Counties_ITM")

# reproject to EPSG:4326
ie_counties = ie_counties.to_crs(4326)

ie_counties.to_file(GPKG_BOUNDARY, layer="OS_IE_Counties")

######################################################################
# NUTS (Nomenclature of territorial units for statistics)

SUB_DIR = os.path.join(DATA_DIR, "nuts-2021", "raw")

# download data if necessary
URL = (
    "https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/"
    "ref-nuts-2021-01m.shp.zip"
)

download_data(server=URL, dl_dir=SUB_DIR)

DATA_FILE = os.path.join(SUB_DIR, "ref-nuts-2021-01m.shp.zip")

# extract the archive
try:
    z = ZipFile(DATA_FILE)
    z.extractall(SUB_DIR)
except BadZipFile:
    print("There were issues with the file", DATA_FILE)

# NUTS1

DATA_FILE = os.path.join(SUB_DIR, "NUTS_RG_01M_2021_4326_LEVL_1.shp.zip")

nuts1 = gpd.read_file(f"zip://{DATA_FILE}!NUTS_RG_01M_2021_4326_LEVL_1.shp")

# filter for Ireland and Northern Ireland
nuts1 = nuts1[nuts1["NUTS_ID"].isin(["IE0", "UKN"])]

nuts1 = nuts1.drop(columns="FID")

nuts1.to_file(GPKG_BOUNDARY, layer="NUTS1")

# NUTS2

DATA_FILE = os.path.join(SUB_DIR, "NUTS_RG_01M_2021_4326_LEVL_2.shp.zip")

nuts2 = gpd.read_file(f"zip://{DATA_FILE}!NUTS_RG_01M_2021_4326_LEVL_2.shp")

nuts2 = nuts2[nuts2["NUTS_ID"].str.contains("IE|UKN")]

nuts2 = nuts2.drop(columns="FID")

nuts2.to_file(GPKG_BOUNDARY, layer="NUTS2")

# NUTS 3

DATA_FILE = os.path.join(SUB_DIR, "NUTS_RG_01M_2021_4326_LEVL_3.shp.zip")

nuts3 = gpd.read_file(f"zip://{DATA_FILE}!NUTS_RG_01M_2021_4326_LEVL_3.shp")

nuts3 = nuts3[nuts3["NUTS_ID"].str.contains("IE|UKN")]

nuts3 = nuts3.drop(columns="FID")

nuts3.to_file(GPKG_BOUNDARY, layer="NUTS3")

# Island of Ireland boundary

ie = nuts1.copy()

ie = ie.dissolve(by="LEVL_CODE", as_index=False)

ie = ie[["geometry"]]

ie = ie.assign(NAME="Ireland")

description = (
    "Boundary for the Island of Ireland generated using NUTS1 boundaries"
)

ie = ie.assign(DESCRIPTION=description)

ie.to_file(GPKG_BOUNDARY, layer="NUTS_Ireland")

# Island of Ireland in Irish transverse mercator

ie.to_crs(2157, inplace=True)

ie.to_file(GPKG_BOUNDARY, layer="NUTS_Ireland_ITM")
