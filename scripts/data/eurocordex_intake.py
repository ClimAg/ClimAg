"""eurocordex_intake.py

Create a data catalogue for EURO-CORDEX at EUR-11 using intake-esm and the
DKRZ data store

<https://intake.readthedocs.io/>
<https://intake-esm.readthedocs.io/>
<https://gitlab.dkrz.de/data-infrastructure-services/intake-esm/>
"""

# import libraries
import json
import os
from datetime import datetime, timezone

import intake
import pooch

DATA_DIR_BASE = os.path.join("data", "EURO-CORDEX")
os.makedirs(DATA_DIR_BASE, exist_ok=True)

timerange = [
    "19710101-19751231",
    "19760101-19801231",
    "19810101-19851231",
    "19860101-19901231",
    "19910101-19951231",
    "19960101-20001231",
    "20010101-20051231",
    "20360101-20401231",
    "20410101-20451231",
    "20460101-20501231",
    "20510101-20551231",
    "20560101-20601231",
    "20610101-20651231",
    "20660101-20701231",
]

# add additional time ranges for the MOHC datasets, which use a 360-day year
timerange = timerange + [t.replace("1231", "1230") for t in timerange]

variables = ["evspsblpot", "pr", "rsds", "tas"]

driving_model_id = [
    "CNRM-CERFACS-CNRM-CM5",
    "ICHEC-EC-EARTH",
    "MPI-M-MPI-ESM-LR",
    "MOHC-HadGEM2-ES",
]

# create local catalogue
dkrz_cat = intake.open_catalog(["https://dkrz.de/s/intake"])

server = dkrz_cat._entries["dkrz_cordex_disk"]._open_args["esmcol_obj"]

dkrz_cordex = intake.open_esm_datastore(
    server,
    read_csv_kwargs={"dtype": {"time_min": "string", "time_max": "string"}},
)

# download data if necessary
FILE_NAME = "dkrz_cordex_disk.json"
KNOWN_HASH = None
if not os.path.isfile(os.path.join(DATA_DIR_BASE, FILE_NAME)):
    pooch.retrieve(
        url=server, known_hash=KNOWN_HASH, fname=FILE_NAME, path=DATA_DIR_BASE
    )

    with open(
        os.path.join(DATA_DIR_BASE, f"{FILE_NAME[:-5]}.txt"),
        "w",
        encoding="utf-8",
    ) as outfile:
        outfile.write(
            f"Data downloaded on: {datetime.now(tz=timezone.utc)}\n"
            f"Download URL: {server}"
        )

# keep data for the relevant variables and time ranges
query = {
    "CORDEX_domain": "EUR-11",
    "experiment_id": ["historical", "rcp45", "rcp85"],
    "frequency": "day",
    "variable_id": variables,
    "time_range": timerange,
    "model_id": "SMHI-RCA4",
    "driving_model_id": driving_model_id,
    "member": ["r1i1p1", "r12i1p1"],
}

cordex_eur11 = dkrz_cordex.search(**query)

# drop v1 of MPI-M-MPI-ESM-LR outputs
cordex_eur11_df = cordex_eur11.df.drop(
    cordex_eur11.df[
        (cordex_eur11.df["driving_model_id"] == "MPI-M-MPI-ESM-LR")
        & (cordex_eur11.df["rcm_version_id"] == "v1")
    ].index
)

# keep only r12i1p1 outputs of ICHEC-EC-EARTH
cordex_eur11_df = cordex_eur11_df.drop(
    cordex_eur11_df[
        (cordex_eur11_df["driving_model_id"] == "ICHEC-EC-EARTH")
        & (cordex_eur11_df["member"] == "r1i1p1")
    ].index
)

# extract driving model name, without the institution
cordex_eur11_df["driving_model"] = cordex_eur11_df["driving_model_id"].replace(
    to_replace={
        "CNRM-CERFACS-CNRM-CM5": "CNRM-CM5",
        "ICHEC-EC-EARTH": "EC-EARTH",
        "MOHC-HadGEM2-ES": "HadGEM2-ES",
        "MPI-M-MPI-ESM-LR": "MPI-ESM-LR",
    }
)

# replace URI to path to downloaded data
cordex_eur11_df["uri"] = (
    DATA_DIR_BASE
    + os.sep
    + "RCA4"
    + os.sep
    + cordex_eur11_df["experiment_id"]
    + os.sep
    + cordex_eur11_df["driving_model"]
    + os.sep
    + cordex_eur11_df["uri"].str.split("/").str[-1]
)

cordex_eur11_df["path"] = cordex_eur11_df["uri"]

# export catalogue
CSV_FILE_PATH = os.path.join(DATA_DIR_BASE, "eurocordex_eur11_catalogue.csv")

cordex_eur11_df.to_csv(CSV_FILE_PATH, index=False)

JSON_FILE_PATH = os.path.join(DATA_DIR_BASE, "dkrz_cordex_disk.json")

# modify the JSON catalogue
with open(JSON_FILE_PATH, encoding="utf-8") as json_file:
    cordex_eur11_cat = json.load(json_file)
    json_file.close()

GITHUB_CSV_LINK = (
    "https://raw.githubusercontent.com/ClimAg/data/main/eurocordex/"
    "eurocordex_eur11_catalogue.csv"
)

cordex_eur11_cat["catalog_file"] = GITHUB_CSV_LINK

cordex_eur11_cat["id"] = "eurocordex_eur11"

cordex_eur11_cat["description"] = (
    "This is an ESM collection for EURO-CORDEX data used in the ClimAg "
    "project. Data has been generated using the DKRZ intake-esm stores. "
    "Data is filtered for the EUR-11 CORDEX domain at the daily frequency, "
    "the historical (1976-2005) and future (rcp45, rcp85) (2041-2070) "
    "experiments, and the following variables: "
    + ", ".join(variables)
    + ". The SMHI-RCA4 RCM and four GCMs ("
    + ", ".join(driving_model_id)
    + ") are the models used to generate these data. Last updated: "
    + str(datetime.now(tz=timezone.utc))
    + "."
)

# save the modified JSON file
JSON_FILE_PATH = os.path.join(DATA_DIR_BASE, "eurocordex_eur11_local.json")

with open(JSON_FILE_PATH, "w", encoding="utf-8") as json_file:
    json.dump(cordex_eur11_cat, json_file, ensure_ascii=False, indent=4)

# create a copy that reads the CSV file from disk
cordex_eur11_cat["catalog_file"] = CSV_FILE_PATH
JSON_FILE_PATH = os.path.join(
    DATA_DIR_BASE, "eurocordex_eur11_local_disk.json"
)
with open(JSON_FILE_PATH, "w", encoding="utf-8") as json_file:
    json.dump(cordex_eur11_cat, json_file, ensure_ascii=False, indent=4)
