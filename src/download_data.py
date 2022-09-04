"""download_data.py

Download data using specified URL and optional parameters into
the specified directory
"""

import os
import requests

BASE_DIR_DEF = "data"


# define function to download data
def download_data(
    server, ddir="data", params=None, chunk_size=1048676
):
    """
    Download data using specified URL and optional parameters into
    the specified directory

    Keyword arguments:
    server -- data download URL
    ddir -- directory where the downloaded files will be stored (default
            "data")
    params -- optional request parameters (default None)
    chunk_size -- number of bytes of data per downloaded chunk (default
                  1048676)
    """
    # get request
    r = requests.get(server, params=params, stream=True)
    # create directory to store files
    os.makedirs(ddir, exist_ok=True)
    # download data to directory
    if r.status_code == 200:
        if r.headers["content-type"] == "application/zip":
            file_name = "data.zip"
        elif r.headers["content-type"] == "application/json":
            file_name = "data.geojson"
        elif r.headers["content-type"] == "text/xml":
            file_name = "data.gml"
        else:
            file_name = "data"
        with open(os.path.join(ddir, file_name), "wb") as filedl:
            for chunk in r.iter_content(chunk_size=chunk_size):
                filedl.write(chunk)
        print("Data successfully downloaded to", ddir)
    else:
        print("Data not downloaded to", ddir, "\nStatus code:", r.status_code)
