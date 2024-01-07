"""Functions for reading input parameters and time series data

Some code has been adapted from:
https://code.europa.eu/agri4cast/modvege
"""

import pandas as pd


def read_params(filename):
    """Read the input parameters (constants) file.

    Parameters
    ----------
    filename : str
        Path to the parameter input file

    Returns
    -------
    dict[str, float]
        A dictionary of the input parameters
    """
    params = (
        pd.read_csv(filename, header=None, index_col=0).squeeze().to_dict()
    )
    return params


def read_timeseries(filename):
    """Read the time series input data

    Parameters
    ----------
    filename : str
        Path to the input time series data file

    Returns
    -------
    tuple[pandas.DataFrame, int]
        A dataframe of the input time series data;
        Length of the data (total number of days)
    """
    timeseries = pd.read_csv(filename, parse_dates=["time"])
    timeseries.sort_values(by=["time"], inplace=True)
    timeseries["doy"] = timeseries.set_index("time").index.dayofyear
    timeseries.reset_index(inplace=True)
    endday = len(timeseries)
    return timeseries, endday
