"""boundary_osi_ie_electoral_divisions_2019.py

OSi National Statutory Boundaries - electoral divisions

import os
exec(
    open(
        os.path.join(
            "scripts", "data", "boundary_osi_ie_electoral_divisions_2019.py"
        ),
        encoding="utf-8"
    ).read()
)
"""

import os
from datetime import datetime, timezone

import geopandas as gpd
import pooch

URL = (
    "https://opendata.arcgis.com/api/v3/datasets/"
    "429c839036934413bb740bea190f2596_0/downloads/data?"
    "format=shp&spatialRefId=2157&where=1%3D1"
)

KNOWN_HASH = None
FILE_NAME = "electoral-divisions-2019.zip"
SUB_DIR = os.path.join("data", "boundaries", "OSi")
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

data = gpd.read_file(
    f"zip://{DATA_FILE}!"
    "Electoral_Divisions_-_OSi_National_Statutory_Boundaries_-_2019.shp"
)

data.to_file(
    os.path.join(SUB_DIR, "osi_national_statutory_boundaries.gpkg"),
    layer="OSi_IE_electoral_divisions_2019",
)
