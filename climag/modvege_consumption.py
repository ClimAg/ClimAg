"""Functions for computing harvested and ingested biomass

References
----------
.. [#Jouven] Jouven, M., Carrère, P., and Baumont, R. (2006). ‘Model
    predicting dynamics of biomass, structure and digestibility of herbage in
    managed permanent pastures. 1. Model description’, Grass and Forage
    Science, 61(2), pp. 112–124.
    https://doi.org/10.1111/j.1365-2494.2006.00515.x.
.. [#Kavanagh] Kavanagh, S. (2016). ‘Feeding the Dairy Cow’, in Teagasc Dairy
    Manual. Teagasc, pp. 201–212. Available at:
    https://www.teagasc.ie/publications/2016/teagasc-dairy-manual.php
    (Accessed: 2 November 2022).
.. [#Broad] Broad, H. J. and Hough, M. N. (1993). ‘The growing and grazing
    season in the United Kingdom’, Grass and Forage Science, 48(1), pp. 26–37.
    https://doi.org/10.1111/j.1365-2494.1993.tb01833.x.
.. [#Eurostat] Eurostat (2022). Glossary:Livestock unit (LSU). Available at:
    https://ec.europa.eu/eurostat/statistics-explained/index.php?title=Glossary:Livestock_unit_(LSU)
    (Accessed: 2 November 2022).
"""

import numpy as np

np.seterr("raise")


def organic_matter_digestibility(
    ts_vals: dict[str, float], params: dict[str, float]
):
    """Organic matter digestibility

    Parameters
    ----------
    ts_vals : dict
        Dictionary with intermediate time series values
    params : dict
        Dictionary of model parameters

    Returns
    -------
    dict
        Updated `ts_vals` dictionary

    Notes
    -----
    Organic matter digestibility of the green vegetative (GV) and green
    reproductive (GR) compartments. See Equations (9) and (10) in [#Jouven]_.

    Digestibility varies among plant parts, with leaves usually being more
    digestible than stems. The differing digestibility of plant parts may
    explain why they are grazed selectively. Digestibility of green
    compartments decreases linearly with compartment age.

    This function returns an updated `ts_vals` dictionary with:

    -   `omd_gv`: Organic matter digestibility of the GV compartment
        [dimensionless]
    -   `omd_gr`: Organic matter digestibility of the GR compartment
        [dimensionless]
    """
    ts_vals["omd_gv"] = max(
        params["max_omd_gv"]
        - (ts_vals["age_gv"] * (params["max_omd_gv"] - params["min_omd_gv"]))
        / params["lls"],
        params["min_omd_gv"],
    )

    ts_vals["omd_gr"] = max(
        params["max_omd_gr"]
        - (ts_vals["age_gr"] * (params["max_omd_gr"] - params["min_omd_gr"]))
        / (params["st_2"] - params["st_1"]),
        params["min_omd_gr"],
    )


def biomass_ingestion(ts_vals: dict[str, float], params: dict[str, float]):
    """Biomass ingestion through grazing

    Parameters
    ----------
    ts_vals : dict
        Dictionary of intermediate time series values
    params : dict
        Dictionary containing of model parameters

    Returns
    -------
    dict
        Updated `ts_vals` dictionary

    Notes
    -----
    The maximum amount of biomass available for grazing and/or harvesting
    for each structural compartment.

    -   Maintain the height of the residual biomass after harvest to the
        minimum residual grass height
    -   The bulk density is used to convert this height to the equivalent
        biomass amount. See [#Jouven]_, sec. "Harvested biomass",
        Equation (19)

    The maximum amount of biomass ingestible based on the stocking rate

    -   Ingestion takes precedence over harvesting
    -   The amount of biomass ingested per livestock unit per day is
        13 kg DM LU⁻¹ based on average consumption data for dairy cows from
        [#Kavanagh]_ and the value used by [#Broad]_
    -   One livestock unit (LU) is equivalent to one dairy cow, based on
        [#Eurostat]_

    The amount of biomass ingested in total for each structural compartment

    -   This is weighted according to the organic matter digestibility.
    -   Digestibility varies among plant parts, with leaves usually being
        more digestible than stems
    -   The differing digestibility of plant parts may explain why they are
        grazed selectively

    **Assumption**: when grazed/harvested, 10% of the available biomass in each
    structural component is lost.

    This function returns an updated `ts_vals` dictionary with:

    -   `i_bm`: The total ingested biomass amount [kg DM ha⁻¹]
    -   `bm_gv`: Updated standing biomass of the green vegetative
        compartment [kg DM ha⁻¹]
    -   `bm_gr`: Updated standing biomass of the green reproductive
        compartment [kg DM ha⁻¹]
    -   `bm_dv`: Updated standing biomass of the dead vegetative
        compartment [kg DM ha⁻¹]
    -   `bm_dr`: Updated standing biomass of the dead reproductive
        compartment [kg DM ha⁻¹]
    """
    if (
        params["sr"] > 0.0
        and params["h_grass"] >= 0.0
        and params["st_g1"] <= ts_vals["st"] <= params["st_g2"]
    ):
        # max available compartmental biomass
        available_biomass = {}
        for key in ["gv", "gr", "dv", "dr"]:
            residual_biomass = params["h_grass"] * params[f"bd_{key}"] * 10.0
            if residual_biomass < ts_vals[f"bm_{key}"]:
                available_biomass[f"bm_{key}"] = (
                    ts_vals[f"bm_{key}"] - residual_biomass
                ) * 0.9
            else:
                available_biomass[f"bm_{key}"] = 0.0

        # max ingested biomass
        max_ingested_biomass = params["sr"] * params["i_bm_lu"]

        # ingested compartmental biomass
        weights = {}
        ingested = {}

        weights_total = (
            ts_vals["omd_gv"]
            + ts_vals["omd_gr"]
            + params["omd_dv"]
            + params["omd_dr"]
        )

        weights["bm_gv"] = ts_vals["omd_gv"] / weights_total
        weights["bm_gr"] = ts_vals["omd_gr"] / weights_total
        weights["bm_dv"] = params["omd_dv"] / weights_total
        weights["bm_dr"] = params["omd_dr"] / weights_total

        # sort by weight
        weights = dict(sorted(weights.items(), key=lambda item: item[1]))

        needed = 0.0

        for key in weights:
            ingested[key] = max_ingested_biomass * weights[key]
            ingested[key] += needed
            if available_biomass[key] < ingested[key]:
                needed = ingested[key] - available_biomass[key]
                ingested[key] = available_biomass[key]
            else:
                needed = 0.0
            # update biomass compartments
            # 10% of biomass is lost during ingestion
            ts_vals[key] -= ingested[key] / 0.9

        # total ingestion
        ts_vals["i_bm"] += (
            ingested["bm_gv"]
            + ingested["bm_gr"]
            + ingested["bm_dv"]
            + ingested["bm_dr"]
        )

        # daily values
        ts_vals["c_bm"] = (
            ingested["bm_gv"]
            + ingested["bm_gr"]
            + ingested["bm_dv"]
            + ingested["bm_dr"]
        )

    else:
        ts_vals["c_bm"] = 0.0


def biomass_harvest(ts_vals: dict[str, float], params: dict[str, float]):
    """Harvest biomass through a cutting event

    Parameters
    ----------
    ts_vals : dict
        A dictionary with intermediate time series values.
    params : dict
        A dictionary containing of model parameters

    Returns
    -------
    dict
        An updated `ts_vals` dictionary

    Notes
    -----
    Harvest biomass through a cutting event at the end of the reproductive
    period.
    Maintain the height of the residual biomass after harvest to the minimum
    cut height. This height is calculated by using the bulk density.
    See [#Jouven]_, sec. "Harvested biomass", Equation (19).

    **Assumption**: during harvest, 10% of the harvestable biomass in each
    structural component is lost.

    This function updated `ts_vals` dictionary with:

    -   `h_bm`: The total harvested biomass amount [kg DM ha⁻¹]
    -   `bm_gv`: Updated standing biomass of the green vegetative
        compartment [kg DM ha⁻¹]
    -   `bm_gr`: Updated standing biomass of the green reproductive
        compartment [kg DM ha⁻¹]
    -   `bm_dv`: Updated standing biomass of the dead vegetative
        compartment [kg DM ha⁻¹]
    -   `bm_dr`: Updated standing biomass of the dead reproductive
        compartment [kg DM ha⁻¹]
    """
    if (
        params["h_grass"] >= 0.0
        and params["st_h1"] <= ts_vals["st"] <= params["st_g2"]
    ):
        harvested_biomass = {}
        for key in ["gv", "gr", "dv", "dr"]:
            # harvested biomass for each compartment
            residual_biomass = params["h_grass"] * params[f"bd_{key}"] * 10.0
            harvested_biomass[f"bm_{key}"] = 0.0
            if residual_biomass < ts_vals[f"bm_{key}"]:
                harvested_biomass[f"bm_{key}"] += (
                    ts_vals[f"bm_{key}"] - residual_biomass
                ) * 0.9
            # update biomass compartments
            # 10% of biomass is lost during harvest
            ts_vals[f"bm_{key}"] -= harvested_biomass[f"bm_{key}"] / 0.9

        # total harvested biomass
        ts_vals["h_bm"] += (
            harvested_biomass["bm_gv"]
            + harvested_biomass["bm_gr"]
            + harvested_biomass["bm_dv"]
            + harvested_biomass["bm_dr"]
        )
