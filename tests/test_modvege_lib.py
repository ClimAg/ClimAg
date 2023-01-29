"""test_modvege_lib.py

Tests for modvege_lib.py

coverage run -m pytest && coverage report -m
"""

import climag.modvege_lib as lm


def test_leaf_area_index():
    """
    Test the leaf area index function
    """

    ts_vals = {
        "bm_gv": 1200.0
    }
    params = {
        "sla": 0.033,
        "pct_lam": 0.68
    }

    assert lm.leaf_area_index(ts_vals=ts_vals, params=params) == 2.6928


def test_actual_evapotranspiration():
    """
    Test the actual evapotranspiration function
    """

    pet = 3.5
    ts_vals = {}

    # the value cannot exceed the potential evapotranspiration
    ts_vals["lai"] = 4.3
    assert lm.actual_evapotranspiration(pet=pet, ts_vals=ts_vals) == 3.5

    ts_vals["lai"] = 0.9
    assert lm.actual_evapotranspiration(pet=pet, ts_vals=ts_vals) == 1.05


def test_potential_growth():
    """
    Test the potential growth function
    """

    par_i = 15.0
    ts_vals = {"lai": 3.5}
    params = {"rue_max": 3.0}

    assert lm.potential_growth(
        par_i=par_i, ts_vals=ts_vals, params=params
    ) == 394.8946072861582


def test_par_function():
    """
    Test the par function
    """

    par_i = 15.0

    assert lm.par_function(par_i=par_i) == 0.5555555555555556


def test_sum_of_temperatures():
    """
    Test the sum of temperatures function
    """

    day = 5
    t_ts = [2.0, 7.0, 4.0, 5.0, 8.0, 3.0, 6.0]
    ts_vals = {}
    params = {"t_0": 4.0}

    ts_vals["st"] = 4.0
    assert lm.sum_of_temperatures(
        params=params, ts_vals=ts_vals, t_ts=t_ts, day=day
    ) == 8.0

    # resetting the sum to zero
    ts_vals["st"] = 0.0
    assert lm.sum_of_temperatures(
        params=params, ts_vals=ts_vals, t_ts=t_ts, day=day
    ) == 4.0


def test_ten_day_moving_avg_temperature():
    """
    Test the ten day moving avg temperature function
    """

    t_ts = [2.1, 7.0, 4.1, 5.0, 8.6, 3.2, 6.0, 9.5, 1.0, 3.3, 5.5, 7.2]

    # when day < 10
    day = 5
    assert lm.ten_day_moving_avg_temperature(day=day, t_ts=t_ts) == 8.6

    # when day >= 10
    day = 11
    assert lm.ten_day_moving_avg_temperature(
        day=day, t_ts=t_ts
    ) == 5.319999999999999


def test_temperature_function():
    """
    Test the temperature function
    """

    ts_vals = {}
    params = {
        "t_0": 4.0,
        "t_1": 10.0,
        "t_2": 20.0,
        "t_max": 40.0
    }

    # t_1 <= t_m10 <= t_2
    ts_vals["t_m10"] = 15.0
    assert lm.temperature_function(ts_vals=ts_vals, params=params) == 1.0

    # t_m10 >= t_max
    ts_vals["t_m10"] = 50.0
    assert lm.temperature_function(ts_vals=ts_vals, params=params) == 0.0

    # t_m10 <= t_0
    ts_vals["t_m10"] = 2.0
    assert lm.temperature_function(ts_vals=ts_vals, params=params) == 0.0

    # t_0 < t_m10 < t_1
    ts_vals["t_m10"] = 7.5
    assert lm.temperature_function(
        ts_vals=ts_vals, params=params
    ) == 0.5833333333333335

    # t_2 < t_m10 < t_max
    ts_vals["t_m10"] = 35.5
    assert lm.temperature_function(
        ts_vals=ts_vals, params=params
    ) == 0.22499999999999987
