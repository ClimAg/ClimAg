"""modvege_mera.py

Simulate grass growth using the MÉRA dataset

Run the following in a Python interpreter in the project's directory and Conda
environment:

import os
exec(
    open(
        os.path.join("scripts", "modvege_mera.py"),
        encoding="utf-8"
    ).read()
)
"""

import os
import sys
from climag.modvege_run import run_modvege

DATA_DIR = os.path.join("data", "ModVege")

# define the name of the input params file
PARAMS_FILE = os.path.join(DATA_DIR, "params.csv")
PARAMS_GPKG_FILE = os.path.join(DATA_DIR, "params.gpkg")

# run the main function
run_modvege(
    input_params_file=PARAMS_FILE,
    input_timeseries_file=os.path.join(
        "data", "MERA", "IE_MERA_FC3hr_3_day.nc"
    ),
    out_dir=DATA_DIR,
    input_params_vector=PARAMS_GPKG_FILE
)

sys.exit()
