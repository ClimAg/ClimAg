"""modvege_eurocordex.py

Simulate grass growth using the EURO-CORDEX climate model dataset

Run the following in a Python interpreter in the project's directory and Conda
environment:

import os
exec(
    open(
        os.path.join("scripts", "modvege_eurocordex.py",
        encoding="utf-8"
    ).read()
)
"""

import os
import glob
from climag.modvege_run import run_modvege

DATA_DIR = os.path.join("data", "ModVege")

# define the name of the input params file
PARAMS_FILE = os.path.join(DATA_DIR, "params.csv")
PARAMS_GPKG_FILE = os.path.join(DATA_DIR, "params.gpkg")

# list of input time series files
ts_files = glob.glob(os.path.join("data", "EURO-CORDEX", "IE", "*.nc"))

# run the main function
for ts in ts_files:
    run_modvege(
        input_params_file=PARAMS_FILE,
        input_timeseries_file=ts,
        out_dir=DATA_DIR,
        input_params_vector=PARAMS_GPKG_FILE
    )
