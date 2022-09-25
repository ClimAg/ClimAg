"""modvege_example.py

Running ModVege using the example data provided in
https://github.com/YannChemin/modvege
"""

# %%
import os
from climag.modvege_run import run_modvege

# %%
DATA_PATH = os.path.join("data", "grass-growth", "modvege")

# define the name of the input params file
PARAMS_FILE = os.path.join(DATA_PATH, "params.csv")
# define the name of the input environment file
WEATHER_FILE = os.path.join(DATA_PATH, "weather.csv")
# outputs
OUT_FILE = os.path.join(DATA_PATH, "output.csv")

# %%
# run the main function
run_modvege(
    input_params_csv=PARAMS_FILE,
    input_weather_csv=WEATHER_FILE,
    out_csv=OUT_FILE
)
