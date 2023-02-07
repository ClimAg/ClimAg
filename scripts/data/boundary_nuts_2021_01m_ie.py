"""boundary_nuts_2021_01m_ie.py

NUTS (Nomenclature of territorial units for statistics) 2021 boundaries
for the Island of Ireland

Run the following in a Python interpreter in the project's directory and Conda
environment:

import os
exec(
    open(
        os.path.join("scripts", "data", "boundary_nuts_2021_01m_ie.py"),
        encoding="utf-8"
    ).read()
)
"""

# import libraries
import os
from datetime import datetime, timezone
from zipfile import BadZipFile, ZipFile
import geopandas as gpd
import pooch

# base data download directory
SUB_DIR = os.path.join("data", "boundaries", "NUTS2021")
os.makedirs(SUB_DIR, exist_ok=True)

URL = (
    "https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/"
    "ref-nuts-2021-01m.shp.zip"
)
KNOWN_HASH = None
FILE_NAME = "ref-nuts-2021-01m.shp.zip"

# file name for the GeoPackage where the boundary vector layers will be saved
GPKG_BOUNDARY = os.path.join("data", "boundaries", "boundaries.gpkg")

DATA_DIR_TEMP = os.path.join(SUB_DIR, "temp")

os.makedirs(DATA_DIR_TEMP, exist_ok=True)

# download data if necessary
if not os.path.isfile(os.path.join(SUB_DIR, FILE_NAME)):
    pooch.retrieve(
        url=URL,
        known_hash=KNOWN_HASH,
        fname=FILE_NAME,
        path=SUB_DIR
    )

    with open(
        os.path.join(SUB_DIR, f"{FILE_NAME[:-8]}.txt"), "w", encoding="utf-8"
    ) as outfile:
        outfile.write(
            f"Data downloaded on: {datetime.now(tz=timezone.utc)}\n"
            f"Download URL: {URL}"
        )

DATA_FILE = os.path.join(SUB_DIR, "ref-nuts-2021-01m.shp.zip")

# extract the archive
try:
    z = ZipFile(DATA_FILE)
    z.extractall(DATA_DIR_TEMP)
except BadZipFile:
    print("There were issues with the file", DATA_FILE)

# NUTS1
DATA_FILE = os.path.join(
    DATA_DIR_TEMP, "NUTS_RG_01M_2021_4326_LEVL_1.shp.zip"
)

ie = gpd.read_file(f"zip://{DATA_FILE}!NUTS_RG_01M_2021_4326_LEVL_1.shp")

# filter for Ireland and UK
ie = ie[ie["CNTR_CODE"].isin(["IE", "UK"])]

# filter for Ireland and Northern Ireland
ie = ie[ie["NUTS_ID"].isin(["IE0", "UKN"])]

# Island of Ireland boundary
ie = ie.dissolve(by="LEVL_CODE", as_index=False)

ie = ie[["geometry"]]

ie = ie.assign(NAME="Ireland")

DESCRIPTION = (
    "Boundary for the Island of Ireland generated using NUTS 2021 Level 1 "
    "boundaries"
)

ie = ie.assign(DESCRIPTION=DESCRIPTION)

ie.to_file(GPKG_BOUNDARY, layer="NUTS_RG_01M_2021_4326_IE")

# Island of Ireland in Irish Transverse Mercator, EPSG:2157
ie.to_crs(2157, inplace=True)

ie.to_file(GPKG_BOUNDARY, layer="NUTS_RG_01M_2021_2157_IE")
