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

# Boundary

osi_roi = osi[["geometry"]].copy()

osi_roi["NAME"] = "Republic of Ireland"

osi_roi = osi_roi.dissolve(by="NAME")

osi_roi.to_file(GPKG_BOUNDARY, layer="Boundary_ROI")

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

# Boundary

osni_ni = osni[["geometry"]].copy()

osni_ni["NAME"] = "Northern Ireland"

osni_ni = osni_ni.dissolve(by="NAME")

osni_ni.to_file(GPKG_BOUNDARY, layer="Boundary_NI")

######################################################################
# All-Ireland counties

osi_counties = osi.dissolve(by="COUNTY")

osi_counties = osi_counties[["CONTAE", "PROVINCE", "geometry"]]

osi_counties.reset_index(inplace=True)

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

ie_counties = osi_counties.merge(osni_counties, how="outer")

ie_counties.to_file(GPKG_BOUNDARY, layer="Counties_IE")

# All-Ireland boundary

ie = ie_counties[["geometry"]].copy()

ie["NAME"] = "Ireland"

ie = ie.dissolve(by="NAME")

ie.to_file(GPKG_BOUNDARY, layer="Boundary_IE")

######################################################################
# NUTS (Nomenclature of territorial units for statistics)

URL = (
    "https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/" +
    "ref-nuts-2021-01m.geojson.zip"
)
SUB_DIR = os.path.join(DATA_DIR, "nuts-2021", "raw")

dd.download_data(server=URL, ddir=SUB_DIR)

DATA_FILE = os.path.join(SUB_DIR, "data.zip")

# NUTS2

nuts = gpd.read_file(
    "zip://" + DATA_FILE + "!" + "NUTS_RG_01M_2021_4326_LEVL_2.geojson"
)

nuts = nuts[nuts["CNTR_CODE"].isin(["IE", "UK"])]

nuts_ie = nuts[nuts["CNTR_CODE"].isin(["IE"])]

nuts_ni = nuts[nuts["NUTS_NAME"].isin(["Northern Ireland"])]

nuts2 = nuts_ie.merge(nuts_ni, how="outer")

nuts2.drop(columns="FID", inplace=True)

nuts2.to_file(GPKG_BOUNDARY, layer="NUTS2_IE")

# NUTS 3

nuts = gpd.read_file(
    "zip://" + DATA_FILE + "!" + "NUTS_RG_01M_2021_4326_LEVL_3.geojson"
)

nuts = nuts[nuts["CNTR_CODE"].isin(["IE", "UK"])]

nuts_ie = nuts[nuts["CNTR_CODE"].isin(["IE"])]

nuts_ni = nuts[nuts["NUTS_ID"].str.contains("UKN0")]

nuts3 = nuts_ie.merge(nuts_ni, how="outer")

nuts3.drop(columns="FID", inplace=True)

nuts3.to_file(GPKG_BOUNDARY, layer="NUTS3_IE")

# Boundaries

ie = nuts3.dissolve(by="CNTR_CODE")

ie = ie[["geometry"]]

ie.to_file(GPKG_BOUNDARY, layer="NUTS_ROI_NI")

ie["NAME"] = "Ireland"

ie = ie.dissolve(by="NAME")

ie.to_file(GPKG_BOUNDARY, layer="NUTS_IE")
