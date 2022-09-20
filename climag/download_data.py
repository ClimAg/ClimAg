"""download_data.py

Download data using specified URL and optional parameters into the specified
directory
"""

import os
import re
from datetime import datetime, timezone
import requests


# define function to download data
def download_data(
    server, dl_dir="data", params=None, chunk_size=1048676
):
    """
    Download data using specified URL and optional parameters into
    the specified directory

    Parameters
    ----------
    server : data download URL
    dl_dir : directory where the downloaded files will be stored (default
        "data")
    params : optional request parameters (default None)
    chunk_size : number of bytes of data per downloaded chunk (default
        1048676)
    """
    # create directory to store files
    os.makedirs(dl_dir, exist_ok=True)
    # download data to directory
    # https://stackoverflow.com/a/53299682
    try:
        r = requests.get(server, params=params, stream=True, timeout=3000)
        r.raise_for_status()  # raise exceptions in case of HTTP errors
        if (
            "Content-Disposition" in r.headers.keys() and
            "filename" in r.headers["Content-Disposition"]
        ):
            file_name = re.findall(
                "filename=(.+)", r.headers["Content-Disposition"]
            )[0].replace('"', '')
        else:
            file_name = server.split("/")[-1]
        with open(os.path.join(dl_dir, file_name), "wb") as file_dl:
            for chunk in r.iter_content(chunk_size=chunk_size):
                file_dl.write(chunk)
        print(
            "Data successfully downloaded to", dl_dir,
            "\nLast downloaded:", datetime.now(tz=timezone.utc)
        )
    except requests.exceptions.RequestException as e:
        print("Data download unsuccessful!", e)
