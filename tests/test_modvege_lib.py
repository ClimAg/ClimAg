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
    ts_vals = {"lai": 4.3}

    assert lm.actual_evapotranspiration(pet=pet, ts_vals=ts_vals) == 3.5

    ts_vals["lai"] = 0.9
    assert lm.actual_evapotranspiration(pet=pet, ts_vals=ts_vals) == 1.05
