"""stats.py

Calculate and save summary statistics as separate files

To-do: RuntimeWarning handling
- All-NaN slice encountered (max, min)
- Invalid value encountered in divide (std)
"""

import glob
import itertools
import os
import sys
import numpy as np
import xarray as xr

season_list = ["DJF", "MAM", "JJA", "SON"]
exp_list = ["historical", "rcp45", "rcp85"]
model_list = ["CNRM-CM5", "EC-EARTH", "HadGEM2-ES", "MPI-ESM-LR"]
dataset_list = ["EURO-CORDEX", "HiResIreland"]
stat_list = ["mean", "std", "max", "min"]

STATS_DIR = os.path.join("data", "ModVege", "stats")
os.makedirs(STATS_DIR, exist_ok=True)


def drop_unneeded_vars(data):
    """
    Drop variables that are not needed
    """

    data = data.drop_vars([
        "bm_gv", "bm_gr", "bm_dv", "bm_dr", "age_gv", "age_gr", "age_dv",
        "age_dr", "omd_gv", "omd_gr", "lai", "env", "wr", "aet"
    ])

    return data


def manage_season_vars(data):
    """
    Keep only relevant variables
    """

    # manage variables
    data["c_bm"].attrs["long_name"] = "Ingested biomass"
    data = data.assign(sen=data["sen_gv"] + data["sen_gr"])
    data["sen"].attrs["long_name"] = "Senescence"
    data["sen"].attrs["units"] = "kg DM ha⁻¹ day⁻¹"
    data = data.assign(abs=data["abs_dv"] + data["abs_dr"])
    data["abs"].attrs["long_name"] = "Abscission"
    data["abs"].attrs["units"] = "kg DM ha⁻¹ day⁻¹"
    data = data.drop_vars([
        "sen_gv", "sen_gr", "abs_dv", "abs_dr", "h_bm", "i_bm"
    ])

    return data


def manage_cumulative_vars(data):
    """
    Cumulative variables
    """

    data = data.drop_vars([
        "sen_gv", "sen_gr", "abs_dv", "abs_dr", "bm", "c_bm", "pgro", "gro"
    ])

    return data


def generate_season_stats(data, stat):
    """
    Generate stats with seasonal groupings
    """

    if stat == "mean":  # weighted mean
        # calculate the weights by grouping month length by season
        weights = (
            data["time"].dt.days_in_month.groupby("time.season") /
            data["time"].dt.days_in_month.groupby("time.season").sum()
        )

        # test that the sum of weights for each season is one
        np.testing.assert_allclose(
            weights.groupby("time.season").sum().values,
            np.ones(len(set(weights["season"].values)))
        )

        # calculate the weighted average
        data_st = (
            (data * weights).groupby("time.season").sum(dim="time")
        )

    elif stat == "std":  # unbiased standard deviation
        data_st = data.groupby("time.season").std(dim="time", ddof=1)

    elif stat == "min":  # minimum
        data_st = data.groupby("time.season").min(dim="time", skipna=True)

    elif stat == "max":  # maximum
        data_st = data.groupby("time.season").max(dim="time", skipna=True)

    # sort seasons in the correct order
    data_st = data_st.reindex(season=season_list)

    # reassign attributes
    data_st.rio.write_crs(data.rio.crs, inplace=True)
    for var in data_st.data_vars:
        data_st[var].attrs = data[var].attrs

    return data_st


def generate_cumulative_stats(data, stat):
    """
    Generate stats for cumulative variables
    """

    data_st = data.groupby("time.year").max(dim="time", skipna=True)

    if stat == "mean":
        data_st = data_st.mean(dim="year", skipna=True)
    elif stat == "std":  # unbiased standard deviation
        data_st = data_st.std(dim="year", skipna=True, ddof=1)
    elif stat == "min":  # minimum
        data_st = data_st.min(dim="year", skipna=True)
    elif stat == "max":  # maximum
        data_st = data_st.max(dim="year", skipna=True)

    # reassign attributes
    data_st.rio.write_crs(data.rio.crs, inplace=True)
    for var in data_st.data_vars:
        data_st[var].attrs = data[var].attrs

    return data_st


# climate model data
for exp, model, dataset in itertools.product(
    exp_list, model_list, dataset_list
):
    # auto-rechunking may cause NotImplementedError with object dtype where it
    # will not be able to estimate the size in bytes of object data
    if model == "HadGEM2-ES":
        CHUNKS = 300
    else:
        CHUNKS = "auto"

    ds = xr.open_mfdataset(
        glob.glob(
            os.path.join(
                "data", "ModVege", dataset, exp, model,
                f"*{dataset}*{model}*{exp}*.nc"
            )
        ),
        chunks=CHUNKS,
        decode_coords="all"
    )

    # convert HadGEM2-ES data back to 360-day calendar
    # this ensures that the correct weighting is applied when
    # calculating the weighted average
    if model == "HadGEM2-ES":
        ds = ds.convert_calendar("360_day", align_on="year")

    # remove spin-up year
    if exp == "historical":
        ds = ds.sel(time=slice("1976", "2005"))
    else:
        ds = ds.sel(time=slice("2041", "2070"))

    # assign new coordinates and dimensions
    ds = ds.assign_coords(exp=exp)
    ds = ds.expand_dims(dim="exp")
    ds = ds.assign_coords(model=model)
    ds = ds.expand_dims(dim="model")

    ds = drop_unneeded_vars(data=ds)

    # seasonal stats
    for s in stat_list:
        ds_sub = manage_season_vars(data=ds)
        ds_sub = generate_season_stats(data=ds_sub, stat=s)
        # save as a new file
        ds_sub.to_netcdf(
            os.path.join(
                STATS_DIR, f"{ds.attrs['input_dataset']}_{s}_season.nc"
            )
        )
        print(f"{ds.attrs['input_dataset']}_{s}_season done!")

        # cumulative stats
        ds_sub = manage_cumulative_vars(data=ds)
        ds_sub = generate_cumulative_stats(data=ds_sub, stat=s)
        # save as a new file
        ds_sub.to_netcdf(
            os.path.join(
                STATS_DIR, f"{ds.attrs['input_dataset']}_{s}_cumulative.nc"
            )
        )
        print(f"{ds.attrs['input_dataset']}_{s}_cumulative done!")

    # second set of outputs for comparison with MERA data
    if exp == "historical":
        for s in stat_list:
            # seasonal stats
            ds_sub = ds.sel(time=slice("1981", "2005"))
            ds_sub = manage_season_vars(data=ds_sub)
            ds_sub = generate_season_stats(data=ds_sub, stat=s)
            # save as a new file
            ds_sub.to_netcdf(
                os.path.join(
                    STATS_DIR,
                    f"{ds.attrs['input_dataset']}_{s}_season_MERA.nc"
                )
            )
            print(f"{ds.attrs['input_dataset']}_{s}_season_MERA done!")

            # cumulative stats
            ds_sub = ds.sel(time=slice("1981", "2005"))
            ds_sub = manage_cumulative_vars(data=ds_sub)
            ds_sub = generate_cumulative_stats(data=ds_sub, stat=s)
            # save as a new file
            ds_sub.to_netcdf(
                os.path.join(
                    STATS_DIR,
                    f"{ds.attrs['input_dataset']}_{s}_cumulative_MERA.nc"
                )
            )
            print(f"{ds.attrs['input_dataset']}_{s}_cumulative_MERA done!")

# MERA
ds = xr.open_mfdataset(
    glob.glob(
        os.path.join("data", "ModVege", "MERA", "*MERA*FC3hr*day*.nc")
    ),
    chunks="auto", decode_coords="all"
)

ds = drop_unneeded_vars(data=ds)

# remove spin-up year
ds = ds.sel(time=slice("1981", "2005"))

for s in stat_list:
    # seasonal stats
    ds_sub = manage_season_vars(data=ds)
    ds_sub = generate_season_stats(data=ds_sub, stat=s)
    # save as a new file
    ds_sub.to_netcdf(
        os.path.join(
            STATS_DIR, f"{ds.attrs['input_dataset']}_{s}_season.nc"
        )
    )
    print(f"{ds.attrs['input_dataset']}_{s}_season done!")

    # cumulative stats
    ds_sub = manage_cumulative_vars(data=ds)
    ds_sub = generate_cumulative_stats(data=ds_sub, stat=s)
    # save as a new file
    ds_sub.to_netcdf(
        os.path.join(
            STATS_DIR, f"{ds.attrs['input_dataset']}_{s}_cumulative.nc"
        )
    )
    print(f"{ds.attrs['input_dataset']}_{s}_cumulative done!")

sys.exit()
