"""modvege_consumption.py

Functions for harvested and ingested biomass.
"""

from dataclasses import dataclass
import numpy as np
import climag.modvege_lib as lm

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
        return (
            self.max_omd_gv -
            (self.age_gv * (self.max_omd_gv - self.min_omd_gv)) / self.lls
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
        return (
            self.max_omd_gr -
            (self.age_gr * (self.max_omd_gr - self.min_omd_gr)) /
            (self.st_2 - self.st_1)
        )


@dataclass
class MaximumAvailableBiomass:
    """
    The maximum amount of biomass available for grazing and/or harvesting for
    each structural compartment.

    Maintain the height of the residual biomass after harvest to the minimum
    cut height. This height is calculated using the bulk density.

    See Jouven et al. (2006), sec. "Harvested biomass", Equation (19).

    Assumption: when grazed/harvested, 10% of the available biomass in each
    structural component is lost.

    Parameters
    ----------
    cut_height : Average height after the cut [m]
    bulk_density : Bulk density [g DM m⁻³]
    biomass : Biomass available [kg DM ha⁻¹]

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
            available_biomass = (self.standing_biomass - residual_biomass) * .9
        else:
            available_biomass = 0
            residual_biomass = self.standing_biomass
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
            val = 0
        return val


@dataclass
class MaximumIngestedBiomass:
    """
    The maximum amount of biomass that can be ingested by all livestock units

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
      units
    """

    stocking_rate: float
    max_ingestion_per_lu: float = 13.0

    def __call__(self) -> float:
        return self.stocking_rate * self.max_ingestion_per_lu


@dataclass
class Ingestion:
    """
    Calculate the amount of biomass ingested in total per structural
    compartment.
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
    cut_height : Average height after the cut [m]
    bulk_density : Bulk density [g DM m⁻³]
    biomass : Biomass available [kg DM ha⁻¹]

    Returns
    -------
    - Harvested biomass [kg DM ha⁻¹]
    - Residual biomass [kg DM ha⁻¹]
    """

    bulk_density: float
    standing_biomass: float
    cut_height: float = 0.05

    def __call__(self) -> float:
        residual_biomass = self.cut_height * self.bulk_density * 10
        if residual_biomass < self.standing_biomass:
            harvested_biomass = (self.standing_biomass - residual_biomass) * .9
        else:
            harvested_biomass = 0
            residual_biomass = self.standing_biomass
        return harvested_biomass


# ########################################################
@dataclass
class MaximumCompartmentalIngestion:
    """
    Maximum biomass ingestion by structural compartment. This is weighted
    according to the organic matter digestibility.

    - Digestibility varies among plant parts, with leaves usually being more
      digestible than stems
    - The differing digestibility of plant parts may explain why they are
      grazed selectively

    Parameters
    ----------
    max_ingested_biomass : Maximum amount of biomass that can be ingested by
        all livestock units
    omd_gv : OMD of GV [dimensionless]
    omd_gr : OMD of GR [dimensionless]
    omd_dv : OMD of DV; default is 0.45 [dimensionless]
    omd_dr : OMD of DR; default is 0.40 [dimensionless]

    Returns
    -------
    - Maximum amount of compartmental biomass that can be ingested by the
      livestock units
    """

    max_ingested_biomass: float
    omd_gv: float
    omd_gr: float
    omd_dv: float = 0.45
    omd_dr: float = 0.40

    def __call__(self) -> float:
        total = self.omd_gv + self.omd_gr + self.omd_dv + self.omd_dr
        return (
            self.max_ingested_biomass / (self.omd_gv / total),
            self.max_ingested_biomass / (self.omd_gr / total),
            self.max_ingested_biomass / (self.omd_dv / total),
            self.max_ingested_biomass / (self.omd_dr / total)
        )


@dataclass
class IngestedBiomassOld:
    """
    Return the amount of biomass ingested by the livestock units based on the
    stocking rate.

    Assume that the animals only feed by grazing on the pasture.

    Parameters
    ----------
    livestock_units : total number of livestock units [LU]
    grazing_area : total grazing area (i.e. grassland available for grazing)
        [ha]
    ingestion_per_livestock_unit : average ingestion of grass by a livestock
        unit (e.g. dairy cow); default is 13 based on Teagasc data
        [kg DM LU⁻¹]
    min_cut_height : minimum residual grass height to be maintained; default
        is 0.05 [m]
    bulk_density : bulk density of the biomass compartment [g DM m⁻³]

    Returns
    -------
    - total grass ingestion by livestock per hectare [kg DM ha⁻¹]

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

    livestock_units: float
    grazing_area: float
    bulk_density: float
    ingestion_per_lu: float = 13.0
    cut_height: float = 0.05

    def __call__(self) -> float:
        total_ingestion = (
            StockingRate(
                livestock_units=self.livestock_units,
                grazing_area=self.grazing_area
            )() * self.ingestion_per_lu
        )
        # if the total ingestion exceeds the minimum cut height of
        # 5 cm, reduce it
        if total_ingestion > self.cut_height * self.bulk_density * 10:
            total_ingestion = (
                total_ingestion - self.cut_height * self.bulk_density * 10
            )
        return total_ingestion


def avDefoliationBiomass(biomass, cutHeight, bulkDensity):
    """
    Estimate biomass available for ingestion by livestock.

    Parameters
    ----------
    biomass : Biomass [kg DM ha⁻¹]
    cutHeight : Default cut height in cut event [m]
    bulkDensity : Bulk density [g DM m⁻³]

    Returns
    -------
    - Biomass for ingestion [kg DM ha⁻¹]
    """

    biomassAfterCut = cutHeight * bulkDensity * 10
    return max(0, biomass - biomassAfterCut)


# def exeDefoliation(biomass, cut_biomass, area):
#     """
#     Defoliation method.
#     ** THIS FUNCTION IS UNUSED!

#     Parameters
#     ----------
#     biomass : Biomass [kg DM ha⁻¹] ?
#     cut_biomass : Biomass removed by cut [kg DM] ?
#     area : Area studied (i.e. pixel area) [?]

#     Returns
#     -------
#     - Updated biomass after defoliation [kg DM ha⁻¹] ?
#     """

#     biomass = biomass - cut_biomass / area
#     if biomass < 0 | np.isnan(biomass):
#         biomass = 0
#     return biomass


def harvest_biomass(bulkDensity, biomass, cutHeight=0.05):
    """
    Realise a cut so that the average height is under cutHeight.
    This height is calculated by using the bulkDensity given in parameter.
    See Jouven et al. (2006), sec. "Harvested biomass", Equation (19).
    Assumption: during harvest, 10% of the harvestable biomass in each
    structural component is lost.

    Parameters
    ----------
    cutHeight : Average height after the cut [m]
    bulkDensity : Bulk density [g DM m⁻³]
    biomass : Biomass available [kg DM ha⁻¹]

    Returns
    -------
    - Harvested biomass [kg DM ha⁻¹]
    - Residual biomass [kg DM ha⁻¹]
    """

    residual_biomass = cutHeight * bulkDensity * 10
    if residual_biomass < biomass:
        harvested_biomass = (biomass - residual_biomass) * .9
    else:
        harvested_biomass = 0
        residual_biomass = biomass
    return (harvested_biomass, residual_biomass)


def exeDefoliationByBiomass(biomass, biomassToIngest):
    """
    Realise a defoliation by providing the biomass to remove.

    Parameters
    ----------
    biomass : Biomass [kg DM ha⁻¹]
    biomassToIngest : Biomass to remove [kg DM ha⁻¹]

    Returns
    -------
    - Biomass left [kg DM ha⁻¹]
    """

    biomass -= biomassToIngest
    return biomass


# DEAD VEGETATIVE FUNCTION
def dv_update(
    gv_gamma, gv_senescent_biomass, temperature, dv_biomass, dv_avg_age
):
    """
    Update dead vegetative compartment.
    See Equation (3) in Jouven et al. (2006).

    Parameters
    ----------
    gv_gamma : Respiratory C loss during senescence (DV)
        (1 - gv_gamma = dv_gamma)
    gv_senescent_biomass : Senescence of GV compartment
    lls : Leaf lifespan (LLS) [500 °C d]
    kldv : Basic abscission rate for the dead vegetative compartment
        (Kl_DV) [0.001]
    temperature : Mean daily temperature (*T*) [°C]
    dv_biomass : DV biomass (BM_DV) [kg DM ha⁻¹]
    dv_avg_age : Average DV age (AGE_DV) [°C d]

    Returns
    -------
    - Dead vegetative biomass
    - Average DV age
    """

    abscissionBiomass = lm.AbscissionDV(
        temperature=temperature, age_dv=dv_avg_age, bm_dv=dv_biomass
    )()
    dv_biomass -= abscissionBiomass
    # at this point the biomass include cut, ingestion, and abscission, not
    # growth
    growthBiomass = (1.0 - gv_gamma) * gv_senescent_biomass
    if dv_biomass + growthBiomass > 0:
        dv_avg_age = (max(0, temperature) + dv_avg_age) * (
            dv_biomass / (dv_biomass + growthBiomass)
        )
    else:
        dv_avg_age = 0
    dv_biomass += growthBiomass
    return (dv_biomass, dv_avg_age)


# DEAD REPRODUCTIVE FUNCTION
def dr_update(
    gr_gamma, gr_senescent_biomass, temperature, dr_biomass, dr_avg_age
):
    """
    Update dead reproductive compartment.

    Parameters
    ----------
    gr_gamma : Parameter to compute growth of DR (1 - gr_gamma = dr_gamma)
    gr_senescent_biomass : Senescence of GR, computed in GR
    st1 : Sum of temperatures at the beginning of the reproductive period
        (ST₁) [600 °C d]
    st2 : Sum of temperatures at the end of the reproductive period
        (ST₂) [1200 °C d]
    temperature : Mean daily temperature (*T*) [°C]
    kldr : Basic abscission rate for the dead reproductive compartment
        (Kl_DR) [0.0005]
    dr_biomass : DR biomass
    dr_avg_age : Average DR age [°C d]

    Returns
    -------
    - Updated biomass for DR
    - Average DR age
    """

    abscissionBiomass = lm.AbscissionDR(
        bm_dr=dr_biomass, temperature=temperature, age_dr=dr_avg_age
    )()
    dr_biomass -= abscissionBiomass
    # at this point the biomass include cut, ingestion and abscission, not
    # growth
    growthBiomass = (1 - gr_gamma) * gr_senescent_biomass
    if dr_biomass + growthBiomass > 0:
        dr_avg_age = max(0, temperature) + dr_avg_age * dr_biomass / (
            dr_biomass + growthBiomass
        )
    else:
        dr_avg_age = 0
    dr_biomass += growthBiomass
    return (dr_biomass, dr_avg_age)


# GREEN VEGETATIVE FUNCTION
def gv_update(gro, a2r, temperature, t0, gv_biomass, gv_avg_age):
    """
    Update green vegetative compartment.

    Parameters
    ----------
    gro : total biomass growth (GRO) [kg DM ha⁻¹]
    a2r : Allocate to reproductive
        (REP in Jouven et al. (2006), reproductive function)
    lls : Leaf lifespan (LLS) [500 °C d]
    temperature : Mean daily temperature (*T*) [°C]
    kgv : Senescence coefficient GV [°C d]
    t0 : Minimum temperature for growth (*T*₀) [4 °C]
    gv_biomass : Updated biomass [kg DM ha⁻¹]
    gv_avg_age : Average GV age [°C day]

    Returns
    -------
    - GV biomass
    - GV average age
    - Senescent biomass
    """

    senescentBiomass = lm.SenescenceGV(
        bm_gv=gv_biomass, temperature=temperature, age_gv=gv_avg_age
    )()
    gv_biomass -= senescentBiomass
    # at this point the biomass include cut, ingestion and senescence, not
    # growth
    if temperature > t0:
        growthBiomass = gro * (1 - a2r)
    else:
        growthBiomass = 0
    if gv_biomass + growthBiomass > 0:
        gv_avg_age = (max(0, temperature) + gv_avg_age) * (
            gv_biomass / (gv_biomass + growthBiomass)
        )
    else:
        gv_avg_age = 0
    gv_biomass += growthBiomass
    return (gv_biomass, gv_avg_age, senescentBiomass)


# GREEN REPRODUCTIVE FUNCTION
def gr_update(temperature, a2r, gro, t0, gr_biomass, gr_avg_age):
    """
    Update green reproductive compartment.

    Parameters
    ----------
    temperature : Mean daily temperature (*T*) [°C]
    a2r : Allocate to reproductive
        (REP in Jouven et al. (2006), reproductive function)
    gro : Total growth (GRO) [kg DM ha⁻¹]
    st1 : Onset of reproductive growth (ST₁) [°C d]
    st2 : End of reproductive growth (ST₂) [°C d]
    kgr : Basic rates of senescence in compartment GR (K_GR)
    t0 : Minimum temperature for growth (*T*₀) [4 °C]
    gr_biomass : GR biomass (BM_GR) [kg DM ha⁻¹]
    gr_avg_age : Average GR age (AGE_GR) [°C d]
    lls : Leaf lifespan (LLS) [°C d] (** UNUSED ARGUMENT!)
    rhogr : Bulk density of GR [g m⁻³] (** UNUSED ARGUMENT!)

    Returns
    -------
    - Updated GR biomass
    - Average GR age
    - Senescent biomass
    """

    senescentBiomass = lm.SenescenceGR(
        bm_gr=gr_biomass, age_gr=gr_avg_age, temperature=temperature
    )()
    gr_biomass -= senescentBiomass
    # at this point the biomass include cut, ingestion and senescence, not
    # growth
    if temperature > t0:
        growthBiomass = gro * a2r
    else:
        growthBiomass = 0
    if gr_biomass + growthBiomass > 0:
        gr_avg_age = (max(0, temperature) + gr_avg_age) * (
            gr_biomass / (gr_biomass + growthBiomass)
        )
    else:
        gr_avg_age = 0
    gr_biomass += growthBiomass
    return (gr_biomass, gr_avg_age, senescentBiomass)


# def cut(
#     cutHeight, rhogv, rhodv, rhogr, rhodr, gvb, dvb, grb, drb
# ):
#     """
#     Realise the harvest on each compartment.

#     Parameters
#     ----------
#     cutHeight : the height of the cut [m]
#     rhogv : bulk density of green vegetative [g m⁻³]
#     rhodv : bulk density of dead vegetative [g m⁻³]
#     rhogr : bulk density of green reproductive [g m⁻³]
#     rhodr : bulk density of dead reproductive [g m⁻³]
#     gvb : standing biomass of green vegetative [kg DM ha⁻¹]
#     dvb : standing biomass of dead vegetative [kg DM ha⁻¹]
#     grb : standing biomass of green reproductive [kg DM ha⁻¹]
#     drb : standing biomass of dead reproductive [kg DM ha⁻¹]

#     Returns
#     -------
#     - the total amount of biomass cut [kg DM ha⁻¹]
#     - the amount of GV biomass cut [kg DM ha⁻¹]
#     - the amount of DV biomass cut [kg DM ha⁻¹]
#     - the amount of GR biomass cut [kg DM ha⁻¹]
#     - the amount of DR biomass cut [kg DM ha⁻¹]
#     """

#     gv_h, gv_b = harvest_biomass(
#         bulkDensity=rhogv, cutHeight=cutHeight, biomass=gvb)
#     dv_h, dv_b = harvest_biomass(
#         bulkDensity=rhodv, cutHeight=cutHeight, biomass=dvb)
#     gr_h, gr_b = harvest_biomass(
#         bulkDensity=rhogr, cutHeight=cutHeight, biomass=grb)
#     dr_h, dr_b = harvest_biomass(
#         bulkDensity=rhodr, cutHeight=cutHeight, biomass=drb)
#     sumBiomassHarvested = gv_h + dv_h + gr_h + dr_h
#     return (sumBiomassHarvested, gv_b, dv_b, gr_b, dr_b)


def getAvailableBiomassForCut(
    gv_biomass, dv_biomass, gr_biomass, dr_biomass, cutHeight, rhogv, rhodv,
    rhogr, rhodr
):
    """
    Return the amount of biomass available for cut

    Parameters
    ----------
    gv_biomass : Biomass of green vegetative
    dv_biomass : Biomass of dead vegetative
    gr_biomass : Biomass of green reproductive
    dr_biomass : Biomass of dead reproductive
    cutHeight : Height of the cut [m]
    rhogv : Volume VV [g m⁻³]
    rhodv : Volume DV [g m⁻³]
    rhogr : Volume GR [g m⁻³]
    rhodr : Volume DR [g m⁻³]

    Returns
    -------
    - Amount of biomass available for cut
    """

    avDefoliationBiomassGV = avDefoliationBiomass(
        biomass=gv_biomass, cutHeight=cutHeight, bulkDensity=rhogv
    )
    avDefoliationBiomassDV = avDefoliationBiomass(
        biomass=dv_biomass, cutHeight=cutHeight, bulkDensity=rhodv
    )
    avDefoliationBiomassGR = avDefoliationBiomass(
        biomass=gr_biomass, cutHeight=cutHeight, bulkDensity=rhogr
    )
    avDefoliationBiomassDR = avDefoliationBiomass(
        biomass=dr_biomass, cutHeight=cutHeight, bulkDensity=rhodr
    )
    return (
        avDefoliationBiomassGV
        + avDefoliationBiomassDV
        + avDefoliationBiomassGR
        + avDefoliationBiomassDR
    )


def defoliation(
    gv_biomass, dv_biomass, gr_biomass, dr_biomass, cutHeight, rhogv, rhodv,
    rhogr, rhodr, maxAmountToIngest=9999999
):
    """
    Defoliation method

    Parameters
    ----------
    gv_biomass : biomass of green vegetative
    dv_biomass : biomass of dead vegetative
    gr_biomass : biomass of green reproductive
    dr_biomass : biomass of dead reproductive
    cutHeight : height of the cut [m]
    rhogv : Volume VV [g m⁻³]
    rhodv : Volume DV [g m⁻³]
    rhogr : Volume GR [g m⁻³]
    rhodr : Volume DR [g m⁻³]
    maxAmountToIngest : The maximum amount of biomass to ingest
        (** USING A DEFAULT OF 9999999 - NEED TO CHECK!)

    Returns
    -------
    - the sum of ingested biomass
    """

    avDefoliationBiomassGV = avDefoliationBiomass(
        biomass=gv_biomass, cutHeight=cutHeight, bulkDensity=rhogv
    )
    avDefoliationBiomassDV = avDefoliationBiomass(
        biomass=dv_biomass, cutHeight=cutHeight, bulkDensity=rhodv
    )
    avDefoliationBiomassGR = avDefoliationBiomass(
        biomass=gr_biomass, cutHeight=cutHeight, bulkDensity=rhogr
    )
    avDefoliationBiomassDR = avDefoliationBiomass(
        biomass=dr_biomass, cutHeight=cutHeight, bulkDensity=rhodr
    )
    sumAvailable = (
        avDefoliationBiomassGV
        + avDefoliationBiomassDV
        + avDefoliationBiomassGR
        + avDefoliationBiomassDR
    )
    sumBiomassIngested = 0
    if sumAvailable > 0:
        if sumAvailable <= maxAmountToIngest:
            sumBiomassIngested = (
                exeDefoliationByBiomass(
                    biomass=avDefoliationBiomassGV,
                    biomassToIngest=maxAmountToIngest
                )
                + exeDefoliationByBiomass(
                    biomass=avDefoliationBiomassDV,
                    biomassToIngest=maxAmountToIngest
                )
                + exeDefoliationByBiomass(
                    biomass=avDefoliationBiomassGR,
                    biomassToIngest=maxAmountToIngest
                )
                + exeDefoliationByBiomass(
                    biomass=avDefoliationBiomassDR,
                    biomassToIngest=maxAmountToIngest
                )
            )
        else:
            sumBiomassIngested = (
                exeDefoliationByBiomass(
                    biomass=(
                        maxAmountToIngest
                        * avDefoliationBiomassGV / sumAvailable
                    ),
                    biomassToIngest=maxAmountToIngest
                )
                + exeDefoliationByBiomass(
                    biomass=(
                        maxAmountToIngest
                        * avDefoliationBiomassDV / sumAvailable
                    ),
                    biomassToIngest=maxAmountToIngest
                )
                + exeDefoliationByBiomass(
                    biomass=(
                        maxAmountToIngest
                        * avDefoliationBiomassGR / sumAvailable
                    ),
                    biomassToIngest=maxAmountToIngest
                )
                + exeDefoliationByBiomass(
                    biomass=(
                        maxAmountToIngest
                        * avDefoliationBiomassDR / sumAvailable
                    ),
                    biomassToIngest=maxAmountToIngest
                )
            )
    return sumBiomassIngested


# TO-DO: This set of functions are either not used or not useful
# def addNI(ni, amountToIncrease):
#     return max(0, min(amountToIncrease + ni, 1.2))


def stocking_rate(livestock_units, grazing_area):
    """
    Calculate the stocking rate

    Parameters
    ----------
    livestock_units : total number of livestock units [LU]
    grazing_area : total grazing area (i.e. grassland available for grazing)
        [ha]

    Returns
    -------
    stocking rate [LU ha⁻¹]
    """

    try:
        stocking_rate_ha = livestock_units / grazing_area
    except ZeroDivisionError:
        stocking_rate_ha = 0
    return stocking_rate_ha


def ingested_biomass(
    livestock_units, grazing_area, bulk_density, min_cut_height=0.05,
    ingestion_per_livestock_unit=13
):
    """
    Return the amount of biomass ingested by the livestock units based on the
    stocking rate.

    Parameters
    ----------
    livestock_units : total number of livestock units [LU]
    grazing_area : total grazing area (i.e. grassland available for grazing)
        [ha]
    ingestion_per_livestock_unit : average ingestion of grass by a livestock
        unit (e.g. dairy cow); default is 13 based on Teagasc data
        [kg DM LU⁻¹]
    min_cut_height : minimum residual grass height to be maintained; default
        is 0.05 [m]
    bulk_density : bulk density of the biomass compartment [g DM m⁻³]

    Returns
    -------
    - total grass ingestion by livestock per hectare [kg DM ha⁻¹]

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

    total_ingestion = (
        stocking_rate(
            livestock_units=livestock_units, grazing_area=grazing_area
        ) * ingestion_per_livestock_unit
    )
    # if the total ingestion exceeds the minimum cut height of 5 cm, reduce it
    if total_ingestion > min_cut_height * bulk_density * 10:
        total_ingestion = total_ingestion - min_cut_height * bulk_density * 10
    return total_ingestion
