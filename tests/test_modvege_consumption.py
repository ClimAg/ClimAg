"""test_modvege_consumption.py

Tests for modvege_consumption.py

coverage run -m pytest && coverage report -m
"""

import numpy as np
import climag.modvege_consumption as cm


def test_stocking_rate():
    """
    Test the stocking rate function
    """

    params = {"lu": 2689.7}

    # with stocking rate
    params["area"] = 831.5
    assert cm.stocking_rate(params=params) == 3.2347564642212867

    # without stocking rate
    params["area"] = 0.0
    # zero division error exception handling
    assert cm.stocking_rate(params=params) == 0.0


def test_organic_matter_digestibility():
    """
    Test the organic matter digestibility function
    """

    ts_vals = {}
    params = {
        "max_omd_gv": 0.9,
        "min_omd_gv": 0.75,
        "max_omd_gr": 0.9,
        "min_omd_gr": 0.65,
        "lls": 500.0,
        "st_1": 120.0,
        "st_2": 2300.0
    }

    ts_vals["age_gv"] = 100.0
    ts_vals["age_gr"] = 500.0
    cm.organic_matter_digestibility(ts_vals=ts_vals, params=params)
    assert ts_vals["omd_gv"] == 0.87
    assert ts_vals["omd_gr"] == 0.8426605504587156

    # ensuring the calculated values don't go below the minimum values
    ts_vals["age_gv"] = 1000.0
    ts_vals["age_gr"] = 3000.0
    cm.organic_matter_digestibility(ts_vals=ts_vals, params=params)
    assert ts_vals["omd_gv"] == 0.75
    assert ts_vals["omd_gr"] == 0.65


def test_biomass_ingestion():
    """
    Test the biomass ingestion function
    """

    ts_vals = {
        "bm_gv": 700.5,
        "bm_gr": 254.7,
        "bm_dv": 307.2,
        "bm_dr": 50.3,
        "omd_gv": 0.79,
        "omd_gr": 0.67,
    }
    params = {
        "bd_gv": 850.0,
        "bd_gr": 300.0,
        "bd_dv": 500.0,
        "bd_dr": 150.0,
        "i_bm_lu": 13.0,
        "omd_dv": 0.45,
        "omd_dr": 0.4,
        "st_g1": 180.0,
        "st_g2": 2300.0
    }

    # with stocking rate
    ts_vals["st"] = 460.5
    ts_vals["i_bm"] = 0.0
    params["sr"] = 2.5
    params["h_grass"] = 0.05
    cm.biomass_ingestion(ts_vals=ts_vals, params=params)
    assert ts_vals["i_bm"] == 43.75541125541125

    # without stocking rate
    ts_vals["i_bm"] = 0.0
    params["sr"] = 0.0
    cm.biomass_ingestion(ts_vals=ts_vals, params=params)
    assert ts_vals["i_bm"] == 0.0

    # without residual grass height
    params["sr"] = 2.5
    params["h_grass"] = np.nan
    cm.biomass_ingestion(ts_vals=ts_vals, params=params)
    assert ts_vals["i_bm"] == 0.0

    # before the grazing season
    ts_vals["st"] = 100.0
    params["h_grass"] = 0.05
    cm.biomass_ingestion(ts_vals=ts_vals, params=params)
    assert ts_vals["i_bm"] == 0.0

    # after the grazing season
    ts_vals["st"] = 2400.0
    cm.biomass_ingestion(ts_vals=ts_vals, params=params)
    assert ts_vals["i_bm"] == 0.0


def test_biomass_harvest():
    """
    Test the biomass ingestion function
    """

    ts_vals = {
        "bm_gv": 700.5,
        "bm_gr": 254.7,
        "bm_dv": 307.2,
        "bm_dr": 50.3,
    }
    params = {
        "bd_gv": 850.0,
        "bd_gr": 300.0,
        "bd_dv": 500.0,
        "bd_dr": 150.0,
        "st_h1": 2100.0,
        "st_g2": 2300.0
    }

    # during the harvest season
    ts_vals["st"] = 2200.0
    ts_vals["h_bm"] = 0.0
    params["h_grass"] = 0.05
    cm.biomass_harvest(ts_vals=ts_vals, params=params)
    assert ts_vals["h_bm"] == 393.65999999999997

    # before the harvest season
    ts_vals["h_bm"] = 0.0
    ts_vals["st"] = 1500.0
    cm.biomass_harvest(ts_vals=ts_vals, params=params)
    assert ts_vals["h_bm"] == 0.0

    # after the harvest season
    ts_vals["st"] = 2500.0
    cm.biomass_harvest(ts_vals=ts_vals, params=params)
    assert ts_vals["h_bm"] == 0.0

    # without residual grass height
    ts_vals["st"] = 2200.0
    params["h_grass"] = np.nan
    cm.biomass_harvest(ts_vals=ts_vals, params=params)
    assert ts_vals["h_bm"] == 0.0