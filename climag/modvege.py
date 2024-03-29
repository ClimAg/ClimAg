"""modvege.py

"""

import numpy as np

import climag.modvege_consumption as cm
import climag.modvege_lib as lm

np.seterr("raise")


def sum_of_temperature_thresholds(timeseries, params) -> dict[str, float]:
    """Calculate sum of temperatures

    Parameters
    ----------
    timeseries : pandas.DataFrame
        Input meteorological time series data
    params : dict
        A dictionary of input parameters

    Returns
    -------
    dict[str, float]
        Dictionary of the sum of temperature thresholds

    Notes
    -----
    Calculate sum of temperatures at:
        - the beginning of the reproductive period (ST₁) [°C d]
        - the end of the reproductive period (ST₂) [°C d]
        - the beginning of the grazing season (STg₁) [°C d]
        - the end of the grazing season (STg₂) [°C d]
        - the beginning of harvest (STh₁) [°C d]

    Nolan and Flanagan (2020) define the thermal growing season length as the
    number of days between the first occurrence of at least six consecutive
    days with a daily mean temperature of > 5°C and the first occurrence of at
    least six consecutive days with a daily mean temperature of < 5°C.

    If the temperatures are too low (i.e. no six consecutive days > 5°C) to
    calculate the start date of the growing season, assume it is the 15th
    March, which is the median date for the midlands and part of northern
    Ireland (Collins and Cummins, 1996; Connaughton, 1973). This is also
    the date when cows are out full time according to Teagasc recommendations
    (Kavanagh, 2016).

    Calculating the end date of the growing season is not straightforward, as
    the temperatues may be too high for there to be six consecutive days
    < 5°C, or these six consecutive days may happen very early in the year,
    and may even be before the growing season start date.

    Grazing season calculations based solely on temperature do not consider
    the delay before sufficient plant cover is available to support grazing
    animals or the ability of animals and machinery to pass over land (Nolan
    and Flanagan, 2020; Collins and Cummins, 1996; based on Keane, 1982).

    The beginning of the grazing season has a delay of 5-15 days after the
    start of the growing season based on Broad and Hough (1993).
    A delay of 10 days is used to allow sufficient reproductive growth.
    ~The end date of the grazing season is determined using the Smith formula
    for calculating the grazing season length (Collins and Cummins, 1996;
    based on Smith, 1976).~

    The end date of the grazing season cannot exceed 1st December. Livestock
    are assumed to be fully housed by 22nd November based on Teagasc
    recommendations (Kavanagh, 2016). The mean latest autumn closing date is
    3rd December, with a two-week interval of 26th November to 9th December
    (Looney et al., 2021).

    Grazing will only take place if there is sufficient biomass available; if
    the residual biomass has a height of less than 5 cm, no grazing or
    harvesting will take place. Therefore, the amount of ingested and
    harvested biomass are mainly influenced by the environmental factors that
    affect growth, such as temperature, radiation, and precipitation.

    The beginning of harvest is assumed to be one day before the grazing
    season ends. Grazing costs less than indoor feeding, so starting the
    harvest just a day before the end of the grazing season ensures grazing is
    maximised.
    The end of harvest is the same as the end of the grazing season.
    """
    st_thresholds = {}

    timeseries.sort_values(by=["time"], inplace=True)
    timeseries.set_index("time", inplace=True)

    # generate index
    try:
        # for multi-year time series, this should correspond to the day
        # of the year
        timeseries["idx"] = timeseries["doy"] - 1
    except KeyError:
        # some datasets do not use standard calendars, so day of the year
        # should not be used
        timeseries["idx"] = range(len(timeseries))

    # return only mean values above t_0, and subtract by t_0
    timeseries.loc[(timeseries["T"] >= params["t_0"]), "Tg"] = (
        timeseries["T"] - params["t_0"]
    )
    # fill NaN values
    timeseries[["Tg"]] = timeseries[["Tg"]].fillna(value=0)

    for year in timeseries.index.year.unique():
        st_thresholds[year] = {}

        # # grazing season length using the Smith formula
        # grazing_season = round(
        #     29.3 * np.mean(timeseries.loc[str(year)]["T"]) -
        #     0.1 * np.sum(timeseries.loc[str(year)]["PP"]) +
        #     19.5
        # )
        # # adjust the length if the dataset has 360 days/year
        # if len(timeseries.loc[str(year)]) == 360:
        #     grazing_season -= 5

        # sum of temperatures at the start
        try:
            start = list(
                timeseries.loc[str(year)]["T"]
                .rolling(6)
                .apply(lambda x: all(x > 5.0))
            ).index(1.0)
            # the lowest possible start index would be 5
            # so force the lowest value to zero
            if start == 5:
                start = 0
        except ValueError:
            # if the temperatures are too low for the start of the growing
            # season to be calculated, assume it is on 15th March
            start = int(timeseries.loc[f"{year}-03-15"]["idx"])

        # beginning of the reproductive period
        st_thresholds[year]["st_1"] = timeseries.loc[str(year)]["Tg"].cumsum()[
            start
        ]

        # beginning of the grazing season
        st_thresholds[year]["st_g1"] = timeseries.loc[str(year)][
            "Tg"
        ].cumsum()[start + 10]

        # end of the reproductive period
        st_thresholds[year]["st_2"] = timeseries.loc[str(year)]["Tg"].cumsum()[
            -1
        ]

        # end of the grazing and harvesting season
        # use the calculated grazing season length
        # if growing season continues in December, end grazing on 1st December
        # grazing_end = min(
        #     int(timeseries.loc[f"{year}-12-01"]["idx"]),
        #     start + 10 + grazing_season
        # )
        grazing_end = int(timeseries.loc[f"{year}-12-01"]["idx"])
        st_thresholds[year]["st_g2"] = timeseries.loc[str(year)][
            "Tg"
        ].cumsum()[grazing_end]

        # beginning of harvest
        st_thresholds[year]["st_h1"] = timeseries.loc[str(year)][
            "Tg"
        ].cumsum()[grazing_end - 1]

    timeseries.reset_index(inplace=True)
    timeseries.drop(columns=["Tg", "idx"], inplace=True)

    return st_thresholds


def modvege(params, tseries, endday=365, t_init=None) -> dict[str, float]:
    """ModVege model as a function

    Parameters
    ----------
    params : dict
        Parameters (constants)
    tseries : pandas.DataFrame
        Time series meteorological data
    endday : int
        Number of days of the year (default is 365)

    Returns
    -------
    dict[str, float]
        Dictionary of results

    Notes
    -----
    Results:
        - Green vegetative biomass [kg DM ha⁻¹]
        - Dead vegetative biomass [kg DM ha⁻¹]
        - Green reproductive biomass [kg DM ha⁻¹]
        - Dead reproductive biomass [kg DM ha⁻¹]
        - Total standing biomass [kg DM ha⁻¹]
        - Potential growth [kg DM ha⁻¹]
        - Total growth [kg DM ha⁻¹]
        - Ingested biomass [kg DM ha⁻¹]
        - Harvested biomass [kg DM ha⁻¹]
        - Leaf area index [dimensionless]
        - Water reserves [mm]
        - Actual evapotranspiration [mm]
        - Environmental limitation of growth [dimensionless]
        - Reproductive function [dimensionless]
    """
    st_thresholds = sum_of_temperature_thresholds(
        timeseries=tseries, params=params
    )

    # dictionary of outputs
    outputs_dict = {
        "time": [],
        "bm_gv": [],
        "bm_gr": [],
        "bm_dv": [],
        "bm_dr": [],
        "age_gv": [],
        "age_gr": [],
        "age_dv": [],
        "age_dr": [],
        "bm": [],
        "pgro": [],
        "gro": [],
        "i_bm": [],
        "h_bm": [],
        "c_bm": [],
        "env": [],
        "lai": [],
        "aet": [],
        "wr": [],
        "sen_gv": [],
        "sen_gr": [],
        "abs_dv": [],
        "abs_dr": [],
        "omd_gv": [],
        "omd_gr": [],
    }

    # dictionary to store intermediate time series values
    ts_vals = {}

    # daily loop
    for i in range(endday):
        # initialise starting parameters
        # standing biomass, biomass age, and water reserves
        # assume that the initial value of WR = WHC
        if i == 0:
            (
                ts_vals["bm_gv"],
                ts_vals["bm_gr"],
                ts_vals["bm_dv"],
                ts_vals["bm_dr"],
                ts_vals["age_gv"],
                ts_vals["age_gr"],
                ts_vals["age_dv"],
                ts_vals["age_dr"],
                ts_vals["wr"],
            ) = (
                params["bm_gv"],
                params["bm_gr"],
                params["bm_dv"],
                params["bm_dr"],
                params["age_gv"],
                params["age_gr"],
                params["age_dv"],
                params["age_dr"],
                params["wr"],
            )

        # initialise ingested/harvested biomass and temperature sum
        if tseries["time"][i].dayofyear == 1:
            # reset to zero on the first day of the year
            ts_vals["i_bm"], ts_vals["h_bm"], ts_vals["st"] = 0.0, 0.0, 0.0
            # sum of temperature thresholds
            for key in st_thresholds[tseries["time"][i].year]:
                params[key] = st_thresholds[tseries["time"][i].year][key]

        # 10-d moving average temperature (T_m10)
        ts_vals["t_m10"] = lm.ten_day_moving_avg_temperature(
            day=(i + 1), t_ts=tseries["T"], t_init=t_init
        )

        # sum of temperatures (ST)
        ts_vals["st"] = lm.sum_of_temperatures(
            params=params, ts_vals=ts_vals, t_ts=tseries["T"], day=(i + 1)
        )

        # temperature function (f(T))
        ts_vals["f_t"] = lm.temperature_function(
            ts_vals=ts_vals, params=params
        )

        # seasonal effect (SEA)
        ts_vals["sea"] = lm.seasonal_effect(params=params)

        # leaf area index (LAI)
        ts_vals["lai"] = lm.leaf_area_index(ts_vals=ts_vals, params=params)

        # actual evapotranspiration (AET)
        ts_vals["aet"] = lm.actual_evapotranspiration(
            pet=tseries["PET"][i], ts_vals=ts_vals
        )

        # water reserves (WR)
        ts_vals["wr"] = lm.water_reserves(
            ts_vals=ts_vals, params=params, precipitation=tseries["PP"][i]
        )

        # water stress (W)
        ts_vals["w"] = lm.water_stress(ts_vals=ts_vals, params=params)

        # water stress function (f(W))
        ts_vals["f_w"] = lm.water_stress_function(
            ts_vals=ts_vals, pet=tseries["PET"][i]
        )

        # environmental limitation of growth (ENV)
        ts_vals["env"] = lm.environmental_limitation(
            ts_vals=ts_vals, params=params, par_i=tseries["PAR"][i]
        )

        # potential growth (PGRO)
        ts_vals["pgro"] = lm.potential_growth(
            par_i=tseries["PAR"][i], ts_vals=ts_vals, params=params
        )

        # total growth (GRO)
        ts_vals["gro"] = lm.total_growth(ts_vals=ts_vals)

        # reproductive function (REP)
        # grazing always takes place during the grazing season if the
        # stocking rate is > 0
        ts_vals["rep"] = lm.reproductive_function(
            params=params, ts_vals=ts_vals
        )

        # senescence (SEN)
        lm.senescence(
            ts_vals=ts_vals, params=params, temperature=tseries["T"][i]
        )

        # abscission (ABS)
        lm.abscission(
            ts_vals=ts_vals, params=params, temperature=tseries["T"][i]
        )

        # standing biomass (BM) and biomass age (AGE)
        lm.standing_biomass(ts_vals=ts_vals, params=params)
        lm.biomass_age(
            temperature=tseries["T"][i], params=params, ts_vals=ts_vals
        )

        # organic matter digestibility (OMD)
        cm.organic_matter_digestibility(ts_vals=ts_vals, params=params)

        # ingested biomass
        cm.biomass_ingestion(ts_vals=ts_vals, params=params)

        # harvested biomass
        cm.biomass_harvest(ts_vals=ts_vals, params=params)

        # recover output streams
        outputs_dict["time"].append(tseries["time"][i])

        # net standing biomass
        outputs_dict["bm"].append(
            ts_vals["bm_gv"]
            + ts_vals["bm_gr"]
            + ts_vals["bm_dv"]
            + ts_vals["bm_dr"]
        )

        for out in [
            "bm_gv",
            "bm_gr",
            "bm_dv",
            "bm_dr",
            "age_gv",
            "age_gr",
            "age_dv",
            "age_dr",
            "pgro",
            "gro",
            "i_bm",
            "h_bm",
            "c_bm",
            "env",
            "lai",
            "aet",
            "wr",
            "sen_gv",
            "sen_gr",
            "abs_dv",
            "abs_dr",
            "omd_gv",
            "omd_gr",
        ]:
            outputs_dict[out].append(ts_vals[out])

    return outputs_dict
