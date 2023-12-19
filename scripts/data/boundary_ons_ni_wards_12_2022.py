"""boundary_ons_ni_wards_12_2022.py

Northern Ireland electoral wards data from ONS Geography

Run the following in a Python interpreter in the project's directory and Conda
environment:

import os
exec(
    open(
        os.path.join("scripts", "data", "boundary_ons_ni_wards_12_2022.py"),
        encoding="utf-8"
    ).read()
)
"""

import os
from datetime import datetime, timezone

import geopandas as gpd
import pooch

FILE_NAME = "wards-uk-12-2022.zip"
URL = (
    "https://opendata.arcgis.com/api/v3/datasets/"
    "a2c204fedefe4120ac93f062c647bdcb_0/downloads/data?"
    "format=shp&spatialRefId=27700&where=1%3D1"
)
KNOWN_HASH = None
SUB_DIR = os.path.join("data", "boundaries", "ONS")
DATA_FILE = os.path.join(SUB_DIR, FILE_NAME)
os.makedirs(SUB_DIR, exist_ok=True)

# download data if necessary
if not os.path.isfile(os.path.join(SUB_DIR, FILE_NAME)):
    pooch.retrieve(
        url=URL, known_hash=KNOWN_HASH, fname=FILE_NAME, path=SUB_DIR
    )

    with open(
        os.path.join(SUB_DIR, f"{FILE_NAME[:-4]}.txt"), "w", encoding="utf-8"
    ) as outfile:
        outfile.write(
            f"Data downloaded on: {datetime.now(tz=timezone.utc)}\n"
            f"Download URL: {URL}"
        )

data = gpd.read_file(f"zip://{DATA_FILE}!WD_DEC_2022_UK_BFC.shp")

# filter NI data
data = data[data["WD22CD"].str.contains("N")]

data.to_file(
    os.path.join("data", "boundaries", "boundaries.gpkg"),
    layer="ONS_NI_wards_12_2022_27700",
)

# reproject to Irish Transverse Mercator
data.to_crs(2157, inplace=True)

data.to_file(
    os.path.join("data", "boundaries", "boundaries.gpkg"),
    layer="ONS_NI_wards_12_2022_2157",
)
