"""Functions for reading input parameters and time series data

Some code has been adapted from:
https://code.europa.eu/agri4cast/modvege
"""

import pandas as pd
import pooch
import os
from datetime import datetime, timezone


def download_data(url, data_dir, file_name, known_hash=None):
    """Download data and store it in the specified directory using Pooch.

    Parameters
    ----------
    url : str
        URL from which the data will be downloaded
    data_dir : str
        Directory to store the downloaded data
    file_name : str
        Name of the downloaded data file with its extension (not full path)
    known_hash : str
        SHA256 hash of downloaded file

    Notes
    -----
    This only downloads data if necessary, i.e. if the data file does not
    already exist in the directory.
    """
    os.makedirs(data_dir, exist_ok=True)
    data_file = os.path.join(data_dir, file_name)
    if not os.path.isfile(data_file):
        pooch.retrieve(
            url=url, known_hash=known_hash, fname=file_name, path=data_dir
        )
        print(f"Data downloaded on: {datetime.now(tz=timezone.utc)}")
        with open(f"{data_file}.txt", "w", encoding="utf-8") as outfile:
            outfile.write(
                f"Data downloaded on: {datetime.now(tz=timezone.utc)}\n"
                f"Download URL: {url}\n"
                f"SHA256 hash: {pooch.file_hash(data_file)}\n"
            )
    else:
        print(f"Data '{file_name}' already exists in '{data_dir}'.")
        with open(f"{data_file}.txt", encoding="utf-8") as f:
            print(f.read())


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
