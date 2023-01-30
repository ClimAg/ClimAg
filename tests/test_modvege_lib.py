"""test_modvege_lib.py

Tests for modvege_lib.py

coverage run -m pytest && coverage report -m
"""

import numpy as np
import climag.modvege_lib as lm


def test_leaf_area_index():
    """
    Test leaf_area_index
    """

    ts_vals = {
        "bm_gv": 1200.0,
        "bm_gr": 500.0
    }
    params = {
        "sla": 0.033,
        "pct_lam": 0.68
    }

    assert lm.leaf_area_index(
        ts_vals=ts_vals, params=params
    ) == 3.8148000000000004


def test_actual_evapotranspiration():
    """
    Test actual_evapotranspiration
    """

    pet = 3.5
    ts_vals = {}

    ts_vals["lai"] = 0.9
    assert lm.actual_evapotranspiration(pet=pet, ts_vals=ts_vals) == 1.05

    # the value cannot exceed the potential evapotranspiration
    ts_vals["lai"] = 4.3
    assert lm.actual_evapotranspiration(pet=pet, ts_vals=ts_vals) == 3.5


def test_potential_growth():
    """
    Test potential_growth
    """

    par_i = 15.0
    ts_vals = {"lai": 3.5}
    params = {"rue_max": 3.0}

    assert lm.potential_growth(
        par_i=par_i, ts_vals=ts_vals, params=params
    ) == 394.8946072861582


def test_par_function():
    """
    Test par_function
    """

    par_i = 15.0

    assert lm.par_function(par_i=par_i) == 0.5555555555555556


def test_sum_of_temperatures():
    """
    Test sum_of_temperatures
    """

    day = 5
    t_ts = [2.1, 7.0, 4.1, 5.0, 8.6, 3.2, 6.0, 9.5, 1.0, 3.3, 5.5, 7.2]
    ts_vals = {}
    params = {"t_0": 4.0}

    ts_vals["st"] = 4.0
    assert lm.sum_of_temperatures(
        params=params, ts_vals=ts_vals, t_ts=t_ts, day=day
    ) == 8.6

    # resetting the sum to zero
    ts_vals["st"] = 0.0
    assert lm.sum_of_temperatures(
        params=params, ts_vals=ts_vals, t_ts=t_ts, day=day
    ) == 4.6


def test_ten_day_moving_avg_temperature():
    """
    Test ten_day_moving_avg_temperature
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
    Test temperature_function
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


def test_seasonal_effect():
    """
    Test seasonal_effect
    """

    ts_vals = {}
    params = {
        "min_sea": 0.8,
        "max_sea": 1.2,
        "st_2": 1200.0
    }

    # st_1 <= 200.0
    params["st_1"] = 150.0
    assert lm.seasonal_effect(ts_vals=ts_vals, params=params) == 1.0

    # st <= 200.0
    params["st_1"] = 600.0
    ts_vals["st"] = 120.0
    assert lm.seasonal_effect(ts_vals=ts_vals, params=params) == 0.8

    # st >= st_2
    ts_vals["st"] = 1400.0
    assert lm.seasonal_effect(ts_vals=ts_vals, params=params) == 0.8

    # st_1 - 200.0 <= st <= st_1 - 100.0
    ts_vals["st"] = 450.0
    assert lm.seasonal_effect(ts_vals=ts_vals, params=params) == 1.2

    # 200.0 < st < st_1 - 200.0
    ts_vals["st"] = 275.0
    assert lm.seasonal_effect(
        ts_vals=ts_vals, params=params
    ) == 0.9500000000000001

    # st_1 - 100.0 < st < st_2
    ts_vals["st"] = 980.0
    assert lm.seasonal_effect(
        ts_vals=ts_vals, params=params
    ) == 0.9257142857142857


def test_water_reserves():
    """
    Test water_reserves
    """

    ts_vals = {
        "wr": 20.0,
        "aet": 20.0
    }
    params = {}

    # high water-holding capacity with precipitation
    precipitation = 20.0
    params["whc"] = 200.5
    assert lm.water_reserves(
        ts_vals=ts_vals, params=params, precipitation=precipitation
    ) == 20.0

    # low water-holding capacity
    params["whc"] = 10.5
    assert lm.water_reserves(
        ts_vals=ts_vals, params=params, precipitation=precipitation
    ) == 10.5

    # no precipitation and high water-holding capacity
    precipitation = 0.0
    params["whc"] = 200.5
    assert lm.water_reserves(
        ts_vals=ts_vals, params=params, precipitation=precipitation
    ) == 0.0


def test_water_stress():
    """
    Test water_stress
    """

    ts_vals = {}
    params = {}

    ts_vals["wr"] = 150.0
    params["whc"] = 200.0
    assert lm.water_stress(ts_vals=ts_vals, params=params) == 0.75

    # high water reserves
    ts_vals["wr"] = 200.0
    assert lm.water_stress(ts_vals=ts_vals, params=params) == 1.0


def test_water_stress_function():
    """
    Test water_stress_function
    """

    ts_vals = {}

    # pet < 3.8
    pet = 2.5
    ts_vals["w"] = 0.1
    assert lm.water_stress_function(ts_vals=ts_vals, pet=pet) == 0.4
    ts_vals["w"] = 0.3
    assert lm.water_stress_function(ts_vals=ts_vals, pet=pet) == 0.875
    ts_vals["w"] = 0.5
    assert lm.water_stress_function(ts_vals=ts_vals, pet=pet) == 0.975
    ts_vals["w"] = 0.7
    assert lm.water_stress_function(ts_vals=ts_vals, pet=pet) == 1.0

    # pet <= 6.5
    pet = 6.25
    ts_vals["w"] = 0.1
    assert lm.water_stress_function(ts_vals=ts_vals, pet=pet) == 0.2
    ts_vals["w"] = 0.3
    assert lm.water_stress_function(
        ts_vals=ts_vals, pet=pet
    ) == 0.5499999999999999
    ts_vals["w"] = 0.5
    assert lm.water_stress_function(
        ts_vals=ts_vals, pet=pet
    ) == 0.7999999999999999
    ts_vals["w"] = 0.65
    assert lm.water_stress_function(ts_vals=ts_vals, pet=pet) == 0.925
    ts_vals["w"] = 0.95
    assert lm.water_stress_function(ts_vals=ts_vals, pet=pet) == 1.0

    # pet > 6.5
    pet = 8.3
    assert lm.water_stress_function(ts_vals=ts_vals, pet=pet) == 0.95


def test_reproductive_function():
    """
    Test reproductive_function
    """

    ts_vals = {}
    params = {
        "st_1": 85.0,
        "st_2": 2250.0,
        "ni": 0.75
    }

    # before reproductive period
    ts_vals["st"] = 34.7
    assert lm.reproductive_function(params=params, ts_vals=ts_vals) == 0.0

    # after reproductive period
    ts_vals["st"] = 2500.5
    assert lm.reproductive_function(params=params, ts_vals=ts_vals) == 0.0

    # with grazing and harvesting
    ts_vals["st"] = 400.0
    ts_vals["i_bm"] = 13.2
    ts_vals["h_bm"] = 10.2
    assert lm.reproductive_function(params=params, ts_vals=ts_vals) == 0.0

    # with grazing and without harvesting
    ts_vals["h_bm"] = 0.0
    assert lm.reproductive_function(params=params, ts_vals=ts_vals) == 0.0

    # without grazing and harvesting
    ts_vals["i_bm"] = 0.0
    assert lm.reproductive_function(
        params=params, ts_vals=ts_vals
    ) == 0.7115384615384616

    # with harvesting and without grazing
    ts_vals["h_bm"] = 10.2
    assert lm.reproductive_function(params=params, ts_vals=ts_vals) == 0.0


def test_environmental_limitation():
    """
    Test environmental_limitation
    """

    par_i = 15.0
    ts_vals = {
        "f_t": 0.5,
        "f_w": 0.75
    }
    params = {"ni": 0.9}

    assert lm.environmental_limitation(
        ts_vals=ts_vals, params=params, par_i=par_i
    ) == 0.1875


def test_total_growth():
    """
    Test total_growth
    """

    ts_vals = {
        "pgro": 150.0,
        "env": 0.25,
        "sea": 0.95
    }

    assert lm.total_growth(ts_vals=ts_vals) == 35.625
