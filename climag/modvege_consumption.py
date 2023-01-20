"""modvege_consumption.py

Functions for harvested and ingested biomass.
"""

import numpy as np

np.seterr("raise")


def stocking_rate(params: dict[str, float]) -> float:
    """
    Calculate the stocking rate

    Parameters
    ----------
    params : A dictionary containing these model parameters:
        - lu : Total number of livestock units [LU]
        - area : Total grassland area available for grazing [ha]

    Returns
    -------
    - An updated `params` dictionary with:
        - sr: Stocking rate [LU ha⁻¹]
    """

    try:
        val = params["lu"] / params["area"]
    except ZeroDivisionError:
        val = 0.0
    return val


def organic_matter_digestibility(
    ts_vals: dict[str, float], params: dict[str, float]
):
    """
    Organic matter digestibility of the green vegetative (GV) and green
    reproductive (GR) compartments.
    See Equation (9) in Jouven et al. (2006a)

    - Digestibility varies among plant parts, with leaves usually being more
      digestible than stems
    - The differing digestibility of plant parts may explain why they are
      grazed selectively
    - Digestibility of green compartments decreases linearly with compartment
      age

    Parameters
    ----------

    ts_vals : A dictionary with intermediate time series values for:
        - age_gv: Age of the GV compartment [°C d]
        - age_gr: Age of the GR compartment [°C d]
    params : A dictionary containing these model parameters:
        - min_omd_gv: Minimum OMD of the GV compartment; default is 0.75
            [dimensionless]
        - max_omd_gv: Maximum OMD of the GV compartment; default is 0.9
            [dimensionless]
        - lls: Leaf lifespan; default is 500 [°C d]
        - min_omd_gr: Minimum OMD of the GR compartment; default is 0.65
            [dimensionless]
        - max_omd_gr: Maximum OMD of the GR compartment; default is 0.9
            [dimensionless]
        - st_1: Sum of temperatures at the beginning of the reproductive
            period; default is 600 [°C d]
        - st_2: Sum of temperatures at the end of the reproductive period;
            default is 1200 [°C d]

    Returns
    -------
    - An updated `ts_vals` dictionary with:
        - omd_gv: Organic matter digestibility of the GV compartment
          [dimensionless]
        - omd_gr: Organic matter digestibility of the GR compartment
          [dimensionless]
    """

    ts_vals["omd_gv"] = max(
        params["max_omd_gv"] -
        (ts_vals["age_gv"] * (params["max_omd_gv"] - params["min_omd_gv"])) /
        params["lls"],
        params["min_omd_gv"]
    )

    ts_vals["omd_gr"] = max(
        params["max_omd_gr"] -
        (ts_vals["age_gr"] * (params["max_omd_gr"] - params["min_omd_gr"])) /
        (params["st_2"] - params["st_1"]),
        params["min_omd_gr"]
    )


def biomass_ingestion(ts_vals: dict[str, float], params: dict[str, float]):
    """
    Biomass ingestion

    - The maximum amount of biomass available for grazing and/or harvesting
      for each structural compartment.
        - Maintain the height of the residual biomass after harvest to the
          minimum residual grass height
        - The bulk density is used to convert this height to the equivalent
          biomass amount. See Jouven et al. (2006a), sec. "Harvested biomass",
          Equation (19)
    - The maximum amount of biomass ingestible based on the stocking rate
        - Ingestion takes precedence over harvesting
        - The amount of biomass ingested per livestock unit per day is
          13 kg DM LU⁻¹ based on average consumption data for dairy cows from
          Teagasc's Dairy Manual (2016) and the value used by Broad and Hough
          (1993)
        - One livestock unit (LU) is equivalent to one dairy cow, based on the
          Eurostat Glossary (2022)
    - The amount of biomass ingested in total for each structural compartment
        - This is weighted according to the organic matter digestibility.
        - Digestibility varies among plant parts, with leaves usually being
          more digestible than stems
        - The differing digestibility of plant parts may explain why they are
          grazed selectively
    - Assumption: when grazed/harvested, 10% of the available biomass in each
      structural component is lost

    Parameters
    ----------
    ts_vals : A dictionary with intermediate time series values for:
        - temperature_sum: Sum of temperatures [°C d]
        - bm_gv: Standing biomass of the green vegetative compartment
          [kg DM ha⁻¹]
        - bm_gr: Standing biomass of the green reproductive compartment
          [kg DM ha⁻¹]
        - bm_dv: Standing biomass of the dead vegetative compartment
          [kg DM ha⁻¹]
        - bm_dr: Standing biomass of the dead reproductive compartment
          [kg DM ha⁻¹]
        - omd_gv: Organic matter digestibility of the green vegetative
          compartment [dimensionless]
        - omd_gr: Organic matter digestibility of the green reproductive
          compartment [dimensionless]
    params : A dictionary containing these model parameters:
        - sr: Stocking rate [LU ha⁻¹]
        - st_1: Sum of temperatures at the beginning of the reproductive
          period; default is 600 [°C d]
        - st_2: Sum of temperatures at the end of the reproductive period;
          default is 1200 [°C d]
        - h_grass: Minimum residual grass height; default is 0.05 [m]
        - bd_gv: Bulk density of the green vegetative compartment; default is
          850 [g DM m⁻³]
        - bd_gr: Bulk density of the green reproductive compartment; default is
          300 [g DM m⁻³]
        - bd_dv: Bulk density of the dead vegetative compartment; default is
          500 [g DM m⁻³]
        - bd_dr: Bulk density of the dead reproductive compartment; default is
          150 [g DM m⁻³]
        - i_bm_lu: Amount of biomass ingested per livestock unit; default is
          13 [kg DM LU⁻¹]
        - omd_dv: Organic matter digestibility of the dead vegetative
          compartment; default is 0.45 [dimensionless]
        - omd_dr: Organic matter digestibility of the dead reproductive
          compartment; default is 0.4 [dimensionless]

    Returns
    -------
    - An updated `ts_vals` dictionary with:
        - i_bm: The total ingested biomass amount [kg DM ha⁻¹]
        - bm_gv: Updated standing biomass of the green vegetative compartment
          [kg DM ha⁻¹]
        - bm_gr: Updated standing biomass of the green reproductive compartment
          [kg DM ha⁻¹]
        - bm_dv: Updated standing biomass of the dead vegetative compartment
          [kg DM ha⁻¹]
        - bm_dr: Updated standing biomass of the dead reproductive compartment
          [kg DM ha⁻¹]
    """

    if (
        params["sr"] > 0.0 and
        params["st_2"] > ts_vals["temperature_sum"] > params["st_1"] + 50.0
    ):

        # max available compartmental biomass
        available_biomass = {}
        for key in ["gv", "gr", "dv", "dr"]:
            residual_biomass = params["h_grass"] * params[f"bd_{key}"] * 10
            if residual_biomass < ts_vals[f"bm_{key}"]:
                available_biomass[f"bm_{key}"] = (
                    (ts_vals[f"bm_{key}"] - residual_biomass) * .9
                )
            else:
                available_biomass[f"bm_{key}"] = 0.0

        # max ingested biomass
        max_ingested_biomass = params["sr"] * params["i_bm_lu"]

        # ingested compartmental biomass
        weights = {}
        ingested = {}
        available = {}

        weights_total = (
            ts_vals["omd_gv"] + ts_vals["omd_gr"] +
            params["omd_dv"] + params["omd_dr"]
        )
        weights["bm_gv"] = ts_vals["omd_gv"] / weights_total
        weights["bm_gr"] = ts_vals["omd_gr"] / weights_total
        weights["bm_dv"] = params["omd_dv"] / weights_total
        weights["bm_dr"] = params["omd_dr"] / weights_total
        # sort by weight
        weights = dict(sorted(weights.items(), key=lambda item: item[1]))

        ingested["bm_gv"] = max_ingested_biomass * weights["bm_gv"]
        ingested["bm_gr"] = max_ingested_biomass * weights["bm_gr"]
        ingested["bm_dv"] = max_ingested_biomass * weights["bm_dv"]
        ingested["bm_dr"] = max_ingested_biomass * weights["bm_dr"]

        available["bm_gv"] = available_biomass["bm_gv"]
        available["bm_gr"] = available_biomass["bm_gr"]
        available["bm_dv"] = available_biomass["bm_dv"]
        available["bm_dr"] = available_biomass["bm_dr"]

        needed = 0.0

        for key in weights:
            ingested[key] += needed
            if available[key] < ingested[key]:
                needed = ingested[key] - available[key]
                ingested[key] = available[key]

        # total ingestion
        ts_vals["i_bm"] += (
            ingested["bm_gv"] + ingested["bm_gr"] +
            ingested["bm_dv"] + ingested["bm_dr"]
        )

        # update biomass compartments
        # 10% of biomass is lost during ingestion
        ts_vals["bm_gv"] -= ingested["bm_gv"] / 0.9
        ts_vals["bm_gr"] -= ingested["bm_gr"] / 0.9
        ts_vals["bm_dv"] -= ingested["bm_dv"] / 0.9
        ts_vals["bm_dr"] -= ingested["bm_dr"] / 0.9


def biomass_harvest(ts_vals: dict[str, float], params: dict[str, float]):
    """
    Harvest biomass through a cutting event at the end of the reproductive
    period.

    Maintain the height of the residual biomass after harvest to the minimum
    cut height.
    This height is calculated by using the bulk density.

    See Jouven et al. (2006a), sec. "Harvested biomass", Equation (19).

    Assumption: during harvest, 10% of the harvestable biomass in each
    structural component is lost.

    Parameters
    ----------
    ts_vals : A dictionary with intermediate time series values for:
        - temperature_sum: Sum of temperatures [°C d]
        - bm_gv: Standing biomass of the green vegetative compartment
          [kg DM ha⁻¹]
        - bm_gr: Standing biomass of the green reproductive compartment
          [kg DM ha⁻¹]
        - bm_dv: Standing biomass of the dead vegetative compartment
          [kg DM ha⁻¹]
        - bm_dr: Standing biomass of the dead reproductive compartment
          [kg DM ha⁻¹]
    params : A dictionary containing these model parameters:
        - st_2: Sum of temperatures at the end of the reproductive period;
          default is 1200 [°C d]
        - h_grass: Minimum residual grass height; default is 0.05 [m]
        - bd_gv: Bulk density of the green vegetative compartment; default is
          850 [g DM m⁻³]
        - bd_gr: Bulk density of the green reproductive compartment; default is
          300 [g DM m⁻³]
        - bd_dv: Bulk density of the dead vegetative compartment; default is
          500 [g DM m⁻³]
        - bd_dr: Bulk density of the dead reproductive compartment; default is
          150 [g DM m⁻³]

    Returns
    -------
    - An updated `ts_vals` dictionary with:
        - h_bm: The total harvested biomass amount [kg DM ha⁻¹]
        - bm_gv: Updated standing biomass of the green vegetative compartment
          [kg DM ha⁻¹]
        - bm_gr: Updated standing biomass of the green reproductive compartment
          [kg DM ha⁻¹]
        - bm_dv: Updated standing biomass of the dead vegetative compartment
          [kg DM ha⁻¹]
        - bm_dr: Updated standing biomass of the dead reproductive compartment
          [kg DM ha⁻¹]
    """

    if (
        params["h_grass"] > 0 and
        params["st_2"] + 50.0 >= ts_vals["temperature_sum"] >= params["st_2"]
    ):
        harvested_biomass = {}
        for key in ["gv", "gr", "dv", "dr"]:
            # harvested biomass for each compartment
            residual_biomass = params["h_grass"] * params[f"bd_{key}"] * 10
            harvested_biomass[f"bm_{key}"] = 0.0
            if residual_biomass < ts_vals[f"bm_{key}"]:
                harvested_biomass[f"bm_{key}"] += (
                    (ts_vals[f"bm_{key}"] - residual_biomass) * .9
                )

        # total harvested biomass
        ts_vals["h_bm"] += (
            harvested_biomass["bm_gv"] + harvested_biomass["bm_gr"] +
            harvested_biomass["bm_dv"] + harvested_biomass["bm_dr"]
        )

        # update biomass compartments
        # 10% of biomass is lost during harvest
        ts_vals["bm_gv"] -= harvested_biomass["bm_gv"] / 0.9
        ts_vals["bm_gr"] -= harvested_biomass["bm_gr"] / 0.9
        ts_vals["bm_dv"] -= harvested_biomass["bm_dv"] / 0.9
        ts_vals["bm_dr"] -= harvested_biomass["bm_dr"] / 0.9
