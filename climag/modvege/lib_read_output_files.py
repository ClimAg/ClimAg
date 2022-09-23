"""lib_read_output_files.py

https://github.com/YannChemin/modvege

This is used only for development, to remove for operational mode
"""

import numpy as np


def read_out(filename):
    """Read the "out_cut.csv" file

    Definition of columns in out_cut.csv
    ------------------------------------
    day
    Mean biomass                      [kg DM ha-1]
    Mean green vegetative biomass     [kg DM ha-1]
    Mean green reproductive biomass   [kg DM ha-1]
    Mean dry vegetative biomass       [kg DM ha-1]
    Mean dry reproductive biomass     [kg DM ha-1]
    Harvested Biomass                 [kg DM ha-1]
    Ingested Biomass                  [kg DM ha-1]
    Mean GRO biomass                  [kg DM ha-1]
    Mean available biomass for cut    [kg DM ha-1]
    """
    arr = np.genfromtxt(filename, delimiter=",", skip_header=0, names=True)
    return arr
