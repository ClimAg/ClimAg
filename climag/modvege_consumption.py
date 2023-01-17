"""modvege_consumption.py

Functions for harvested and ingested biomass.
"""

from dataclasses import dataclass
import numpy as np

np.seterr("raise")


@dataclass
class OrganicMatterDigestibilityGV:
    """
    Organic matter digestibility of the green vegetative (GV) compartment.
    See Equation (9) in Jouven et al. (2006)

    - Digestibility varies among plant parts, with leaves usually being more
      digestible than stems
    - The differing digestibility of plant parts may explain why they are
      grazed selectively
    - Digestibility of green compartments decreases linearly with compartment
      age

    Parameters
    ----------
    age_gv : Age of the GV compartment [°C d]
    min_omd_gv : Minimum OMD of the GV compartment; default is 0.75
        [dimensionless]
    max_omd_gv : Maximum OMD of the GV compartment; default is 0.9
        [dimensionless]
    lls : Leaf lifespan; default is 500 [°C d]

    Returns
    -------
    - Organic matter digestibility of the GV compartment [dimensionless]
    """

    age_gv: float
    min_omd_gv: float = 0.75
    max_omd_gv: float = 0.9
    lls: float = 500.0

    def __call__(self) -> float:
        return max(
            self.max_omd_gv -
            (self.age_gv * (self.max_omd_gv - self.min_omd_gv)) / self.lls,
            self.min_omd_gv
        )


@dataclass
class OrganicMatterDigestibilityGR:
    """
    Organic matter digestibility of the green reproductive (GR) compartment.
    See Equation (9) in Jouven et al. (2006)

    - Digestibility varies among plant parts, with leaves usually being more
      digestible than stems
    - The differing digestibility of plant parts may explain why they are
      grazed selectively
    - Digestibility of green compartments decreases linearly with compartment
      age

    Parameters
    ----------
    age_gr : Age of the GR compartment [°C d]
    min_omd_gr : Minimum OMD of the GR compartment; default is 0.65
        [dimensionless]
    max_omd_gr : Maximum OMD of the GR compartment; default is 0.9
        [dimensionless]
    st_1 : Sum of temperatures at the beginning of the reproductive period;
        default is 600 [°C d]
    st_2 : Sum of temperatures at the end of the reproductive period; default
        is 1200 [°C d]

    Returns
    -------
    - Organic matter digestibility of the GR compartment [dimensionless]
    """

    age_gr: float
    min_omd_gr: float = 0.65
    max_omd_gr: float = 0.9
    st_1: float = 600.0
    st_2: float = 1200.0

    def __call__(self) -> float:
        return max(
            self.max_omd_gr -
            (self.age_gr * (self.max_omd_gr - self.min_omd_gr)) /
            (self.st_2 - self.st_1),
            self.min_omd_gr
        )


@dataclass
class MaximumAvailableBiomass:
    """
    The maximum amount of biomass available for grazing and/or harvesting for
    each structural compartment.

    Maintain the height of the residual biomass after harvest to the minimum
    cut height. The bulk density is used to convert this height to the
    equivalent biomass amount.

    See Jouven et al. (2006), sec. "Harvested biomass", Equation (19).

    Assumption: when grazed/harvested, 10% of the available biomass in each
    structural component is lost.

    Parameters
    ----------
    cut_height : Average height after the cut [m]
    bulk_density : Bulk density [g DM m⁻³]
    standing_biomass : Biomass available [kg DM ha⁻¹]

    Returns
    -------
    - Maximum available biomass [kg DM ha⁻¹]
    - Residual biomass [kg DM ha⁻¹]
    """

    bulk_density: float
    standing_biomass: float
    cut_height: float = 0.05

    def __call__(self) -> float:
        residual_biomass = self.cut_height * self.bulk_density * 10
        if residual_biomass < self.standing_biomass:
            available_biomass = (
                (self.standing_biomass - residual_biomass) * .9
            )
        else:
            available_biomass = 0.0
        return available_biomass


@dataclass
class StockingRate:
    """
    Calculate the stocking rate

    Parameters
    ----------
    livestock_units : Total number of livestock units [LU]
    grazing_area : Total grazing area (i.e. grassland available for grazing)
        [ha]

    Returns
    -------
    - Stocking rate [LU ha⁻¹]
    """

    livestock_units: float
    grazing_area: float

    def __call__(self) -> float:
        try:
            val = self.livestock_units / self.grazing_area
        except ZeroDivisionError:
            val = 0.0
        return val


@dataclass
class MaximumIngestedBiomass:
    """
    The maximum amount of biomass that can be ingested by all livestock units
    based on the stocking rate

    Ingestion takes precedence over harvesting

    Parameters
    ----------
    stocking_rate : Stocking rate [LU ha⁻¹]
    max_ingestion_per_lu : Maximum amount of grass that can be ingested by a
        livestock unit; using the average ingestion of grass by a livestock
        unit (e.g. dairy cow); default is 13 based on Teagasc data
        [kg DM LU⁻¹]

    Returns
    -------
    - Maximum amount of biomass that can be ingested in total by all livestock
      units [kg DM ha⁻¹]

    Notes
    -----
    - Grass10: https://www.teagasc.ie/crops/grassland/grass10/
        - average ingestion of grass for dairy cows is 13.41 kg DM LU⁻¹
        - average ingestion of supplements (meal, concentrate, silage) is
          4.25 kg DM LU⁻¹
    - Teagasc Dairy Manual:
      https://www.teagasc.ie/publications/2016/teagasc-dairy-manual.php
        - 8-13 kg DM grass per cow in the spring
        - increase of 0.75-1.0 kg DM until peak intake is reached
        - peak intake of 16-18 kg DM
        - average intake is 13 kg DM ((8 + 18) / 2)
    - one dairy cow is equivalent to one livestock unit; see
      https://cap-calculators.apps.rhos.agriculture.gov.ie/stocking-rate
    """

    stocking_rate: float
    max_ingestion_per_lu: float = 13.0

    def __call__(self) -> float:
        return self.stocking_rate * self.max_ingestion_per_lu


@dataclass
class Ingestion:
    """
    Calculate the amount of biomass ingested in total for each structural
    compartment. This is weighted according to the organic matter
    digestibility.

    - Digestibility varies among plant parts, with leaves usually being more
      digestible than stems
    - The differing digestibility of plant parts may explain why they are
      grazed selectively

    Parameters
    ----------
    max_ingested_biomass : Maximum amount of biomass that can be ingested by
        all livestock units [kg DM ha⁻¹]
    omd_gv : OMD of GV [dimensionless]
    omd_gr : OMD of GR [dimensionless]
    omd_dv : OMD of DV; default is 0.45 [dimensionless]
    omd_dr : OMD of DR; default is 0.40 [dimensionless]
    bm_gv_av : Available GV biomass [kg DM ha⁻¹]
    bm_gr_av : Available GR biomass [kg DM ha⁻¹]
    bm_dv_av : Available DV biomass [kg DM ha⁻¹]
    bm_dr_av : Available DR biomass [kg DM ha⁻¹]

    Returns
    -------
    - Biomass ingested by biomass compartment, returned as a dict
    """

    bm_gv_av: float
    bm_gr_av: float
    bm_dv_av: float
    bm_dr_av: float
    max_ingested_biomass: float
    omd_gv: float
    omd_gr: float
    omd_dv: float = 0.45
    omd_dr: float = 0.40

    def __call__(self) -> float:
        weights = {}
        ingested = {}
        available = {}

        weights_total = self.omd_gv + self.omd_gr + self.omd_dv + self.omd_dr
        weights["gv"] = self.omd_gv / weights_total
        weights["gr"] = self.omd_gr / weights_total
        weights["dv"] = self.omd_dv / weights_total
        weights["dr"] = self.omd_dr / weights_total
        # sort by weight
        weights = dict(sorted(weights.items(), key=lambda item: item[1]))

        ingested["gv"] = self.max_ingested_biomass * weights["gv"]
        ingested["gr"] = self.max_ingested_biomass * weights["gr"]
        ingested["dv"] = self.max_ingested_biomass * weights["dv"]
        ingested["dr"] = self.max_ingested_biomass * weights["dr"]

        available["gv"] = self.bm_gv_av
        available["gr"] = self.bm_gr_av
        available["dv"] = self.bm_dv_av
        available["dr"] = self.bm_dr_av

        needed = 0.0

        for key in weights:
            ingested[key] += needed
            if available[key] < ingested[key]:
                needed = ingested[key] - available[key]
                ingested[key] = available[key]
        return ingested


@dataclass
class HarvestedBiomass:
    """
    Harvest biomass through cuts.

    Maintain the height of the residual biomass after harvest to the minimum
    cut height.
    This height is calculated by using the bulk density.

    See Jouven et al. (2006), sec. "Harvested biomass", Equation (19).

    Assumption: during harvest, 10% of the harvestable biomass in each
    structural component is lost.

    Parameters
    ----------
    cut_height : minimum residual grass height to be maintained; default
        is 0.05 [m]
    bulk_density : Bulk density of the biomass compartment [g DM m⁻³]
    standing_biomass : Standing biomass [kg DM ha⁻¹]
    t_sum : Sum of temperatures [°C d]
    st_2 : Sum of temperatures at the end of the reproductive period; default
        is 1200 (ST₂) [°C d]

    Returns
    -------
    - Harvested biomass [kg DM ha⁻¹]
    """

    bulk_density: float
    standing_biomass: float
    t_sum: float
    cut_height: float = 0.05
    st_2: float = 1200.0

    def __call__(self) -> float:
        harvested_biomass = 0.0
        if (
            self.cut_height > 0 and
            self.st_2 >= self.t_sum >= self.st_2 - 50.0
        ):
            residual_biomass = self.cut_height * self.bulk_density * 10
            if residual_biomass < self.standing_biomass:
                harvested_biomass += (
                    (self.standing_biomass - residual_biomass) * .9
                )
        return harvested_biomass
