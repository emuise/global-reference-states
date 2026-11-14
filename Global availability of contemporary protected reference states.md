# Global availability of contemporary protected reference states

Working design document · 13 November 2026

## 1. Objective

Estimate, at 1 km resolution and global scale, how much unprotected land has a suitable contemporary reference state: a protected, natural, low-modification pixel in a comparable environment. Results are reported as the proportion of land with a suitable reference, by country and by ecosystem.

The conceptual basis is the contemporary reference state framework of McNellie et al. (2020, *Global Change Biology*, doi:10.1111/gcb.15383). It defines references from current ecological patterns rather than historical conditions, so that they are measurable and falsifiable. A literature search did not turn up an existing global-scale availability assessment. Related work addresses protected area representativeness (Sayre et al. 2020; the GEO BON PARC indices) rather than reference availability. That search was not exhaustive and should be repeated before submission.

## 2. Design decisions

| Decision | Choice |
| --- | --- |
| Focal set (pixels needing a reference) | All unprotected land, natural or modified |
| Reference pool (candidate matches) | Protected pixels (IUCN Ia to IV, plus OECMs) with natural land cover and low human modification |
| Human modification | Hard filter on the pool (gHM ≤ 0.10), not a matching covariate |
| Matching method | Nearest neighbour on a whitened climate and topography distance (equivalent to Mahalanobis), with a caliper in SD units |
| Caliper | One global value, calibrated from the reference pool itself (section 6) |
| Exact-match strata, natural focal pixels | Ecoregion × modal natural cover class (like-for-like) |
| Exact-match strata, modified focal pixels | Ecoregion only (environment-only) |
| Search beyond the ecoregion | None. A focal pixel with no match in its ecoregion is unmatched |
| Domestic availability | A separate country-restricted run (strata additionally split by country) |
| Tie-break | Smallest geographic distance among candidates within the caliper |
| Reference year | 2021 |
| Outputs | Area-weighted proportions of land with a suitable reference, by country, ecoregion and biome, split by focal group |

## 3. Changes from the original plan

The original plan was caliper-based propensity score matching with exact matching on ecosystem and land cover. Four things changed.

- **Propensity scores replaced by environmental distance.** A propensity score balances covariates; it does not find environmental analogs. Two pixels with the same score can be far apart in climate and topography space, and protection is strongly biased toward high, remote, cold land. A direct distance also scales: kd-tree searches per stratum are feasible where a global matching package is not.
- **Land cover is no longer an exact-match variable for all pixels.** The focal set includes cropland and built-up land, which could only match protected cropland. That is not a reference.
- **The ecosystem stratum does not depend on current land use.** WWF ecoregions are independent of current cover, unlike ecosystem maps built partly from land cover.
- **Land cover fractions define the pool and the strata, not the distance.** Fractions are compositional and zero-inflated, which suits Mahalanobis distance poorly.

## 4. Data

All layers are aligned to a single 30 arc-second grid (about 1 km) from 56° S to 84° N, the extent of the topography layer. The pipeline downloads every layer itself.

| Layer | Source | Notes |
| --- | --- | --- |
| Climate | WorldClim 2.1, 30 arc-sec, all 19 bioclim variables | 1970 to 2000 normals; the mismatch with 2021 does not matter for long-term climate |
| Topography | EarthEnv (Amatulli et al. 2018), 1 km | Elevation, slope, roughness and TPI (median of GMTED2010 median) |
| Countries | Natural Earth 1:10m admin 0 | Coded by ADM0_A3; ISO_A3 is "-99" for several countries |
| Ecoregions | WWF / RESOLVE Ecoregions 2017 | Exact-match stratum; biome and realm used for reporting |
| Land cover | ESA WorldCover 2021 (v200), 10 m | Aggregated to 1 km as class fractions of valid pixels |
| Human modification | gHM v3 (Theobald et al. 2025), 2020, all threats combined, 300 m | Area-averaged to 1 km. 2020 is the nearest year to 2021 in the consistent time series; the 2022 static layer is the alternative |
| Protected areas | WDPA and WD-OECM (Protected Planet), current release filtered to STATUS_YR ≤ 2021 | See below |

**Protected areas.** Polygons only; point-only records are dropped. Sites must be designated, inscribed or established, not marine, and not UNESCO MAB reserves. Sites with an unknown status year (0) are kept. The reference pool uses IUCN categories Ia to IV plus all OECMs (the "strict" set); a sensitivity run uses every category plus OECMs (the "all" set). Archived WDPA releases are not openly available, so the 2021 snapshot is approximated by filtering a current release on status year. That misses sites removed or downgraded since 2021. The download date of the release is recorded alongside the data.

**Rasterization.** Per 10° tile:

- **Protected fraction:** the share of the cell covered by the union of PAs, with polygons dissolved first so that overlaps are not double counted. A cell is protected when this is at least 0.50.
- **PA identity:** the PA covering the largest share of the cell. Ties go to the smaller PA. This replaces the earlier centre-point rule and neighbourhood fill.

## 5. Definitions

**Land cover classes.** WorldCover natural classes are tree cover, shrubland, grassland, herbaceous wetland, mangroves, moss and lichen, and bare or sparse vegetation. Modified classes are cropland and built-up. Pixels with more than half permanent water or snow and ice are excluded from both the pool and the focal set.

**Collapsed natural classes.** Five: wetland (herbaceous wetland plus mangroves), tree, shrub, grass and sparse (bare or sparse vegetation plus moss and lichen).

**Natural pixel.** Natural-class fractions sum to at least X = 0.90.

**Modal natural class.** The collapsed class with the largest share of the pixel. There is no minimum share, so no natural pixel is "mixed". Ties are broken by a fixed order: wetland, tree, shrub, grass, sparse. Wetland comes first because it is the most distinctive class and the hardest to match. The rest run from most to least vegetation structure. Exact ties are rare because fractions come from 10,000 sub-pixels.

**Modified pixel.** Any non-excluded pixel that fails the natural test. These are matched environment-only.

**Reference pool.** Protected (strict set), natural, gHM ≤ 0.10, and assigned a PA identity.

**Focal groups.**

| Group | Run | Exact-match strata | Distance on |
| --- | --- | --- | --- |
| Natural, unprotected (like-for-like) | `lfl` | Ecoregion × modal natural class | Climate and topography |
| Natural, unprotected (environment-only, for comparison) | `env`, group 1 | Ecoregion | Climate and topography |
| Modified, unprotected | `env`, group 2 | Ecoregion | Climate and topography |

Comparing the first two rows shows how much availability depends on matching cover type. Each run has a country-restricted twin (`_cty`) whose strata are additionally split by country.

## 6. Matching procedure

1. **Environmental space.** Run a PCA on a regular global sample of the 19 bioclim and 4 topography variables, keep the components that explain 95% of the variance, and scale each to unit variance. Euclidean distance in this space is Mahalanobis distance in the retained subspace. The transformation is estimated once and held fixed.
2. **Caliper calibration.** Draw a sample of about 200,000 pool pixels, proportional to pool area. For each, find the environmental distance to the nearest pool pixel in the same ecoregion that lies in a different PA. The caliper is the 0.90 quantile of these distances, applied globally.
   - Pixels in the same PA are excluded because adjacent cells are near-identical and would collapse the caliper towards zero.
   - In practice, the closest neighbours of a pixel in a large PA can all lie within that PA. Those pixels are searched directly against the pool outside their PA, so their larger distances are not lost.
   - Pixels in ecoregions with only one PA have no comparison and are left out.

   The rule reads as: *a focal pixel has a suitable reference if it is at least as close to one as 90% of reference pixels are to a reference in another protected area.*
3. **Candidate search.** Within each stratum, search the reference pool for up to k = 100 neighbours of each focal pixel within the caliper.
4. **Tie-break.** Among the candidates, choose the one with the smallest geographic distance, computed on unit-sphere coordinates.
5. **Records per focal pixel.** The matched pool pixel, its country and its PA, the environmental distance, the geographic distance, the number of candidates, and the number of distinct PAs among them.
6. **Failure reason.** Recorded separately: no reference pixels in the stratum, or none within the caliper. No search is made outside the ecoregion.

Processing is done per ecoregion, with focal pixels handled in chunks, because the global table is too large for memory on a typical workstation.

## 7. Protected area importance and cross-country reliance

**Protected area importance.** Each pool cell carries the ID of the PA with the largest coverage of the cell (section 4). Three measures are reported per PA.

| Measure | Definition | Note |
| --- | --- | --- |
| Chosen reliance | Focal area whose chosen (nearest) reference lies in the PA | Credit goes to the nearest PA, so it partly reflects geography |
| Shared reliance | Each focal pixel's area split equally across the distinct PAs among its candidates | Removes the tie-break effect |
| Sole-supplier area | Focal area for which the PA is the only supplier within the caliper | Closest available measure of irreplaceability; counted only where the candidate list was not capped at k |

Each measure is also divided by the PA's qualifying pool area, so large PAs do not rank highly through size alone.

**Reliance on other countries.** Each run is paired with its country-restricted twin. The two are joined by cell, giving three measures:

- **Domestic.** The focal pixel has a reference within the caliper in its own country (matched in the country-restricted run).
- **Only foreign.** Matched in the main run, but not in the country-restricted run. This is the strict measure of dependence. Because the country run searches the whole domestic pool, it is exact and not limited by k.
- **Chosen reference is foreign.** What the nearest-match rule produces. It overstates dependence, because the nearest match can lie across a border even when a domestic candidate exists.

Foreign chosen matches are binned by geographic distance (50 km bins) to separate trivial cross-border matches from long-distance reliance. Summing foreign-chosen area per PA shows which PAs serve other countries. Continent-level reliance is not reported; ecoregions already capture the relevant geography.

**Caveats for this section.**

- Ecoregions span borders, so small countries will show high foreign reliance partly by construction.
- Shared and sole-supplier tallies use candidate lists capped at k = 100, so they are approximate where the cap is reached.
- Cells without a country (coastal and border edge effects) are excluded from the country statistics and from the country-restricted runs.

## 8. Reporting

- Area-weighted proportion of unprotected land with a suitable reference, by country, ecoregion and biome, for each focal group and overall.
- Per stratum: number of focal and reference pixels, their areas, the number of distinct PAs, and the number of independent protected patches. Patches are 8-connected components of the baseline pool. A stratum with a handful of reference cells should not be read like one with thousands.
- Reuse ratio: focal pixels served per reference pixel used.
- Environmental distance of accepted matches, as a fraction of the caliper, and the area failing for each reason.
- Matched area by distance from the focal pixel to the nearest protected cell (bins at 1, 5, 10, 25 and 50 km). Proximity may reflect spillover or shared disturbance rather than an independent benchmark.
- Cross-country reliance and per-protected-area importance (section 7).

## 9. Sensitivity analyses

One setting is changed from the baseline at a time. Each variant runs like-for-like, environment-only, and both country-restricted twins.

| Variant | Change |
| --- | --- |
| `ghm_0.4` | gHM threshold 0.40 |
| `ghm_none` | No modification filter (pixels without gHM kept) |
| `calq_0.75`, `calq_0.95` | Caliper quantile 0.75 or 0.95 |
| `nat_0.80`, `nat_0.95` | Natural threshold X = 0.80 or 0.95 |
| `pa_all` | All IUCN categories plus OECMs |

The caliper value is held at its baseline calibration in every variant except the caliper ones. That way "environmentally similar" means the same thing when the pool changes. Patch counts and distance to PA always use the baseline pool and the strict PA set. The gap between the filtered and unfiltered pool is likely to be a headline result. Changing the PCA variance or the variable set requires re-running the tabulation.

## 10. Caveats and risks

- WorldCover cannot separate plantation from natural forest, or pasture from natural grassland. The human modification filter carries that distinction.
- Low human modification correlates with remoteness, elevation and cold climate, which is the protected area bias being measured. Expect the filtered pool to represent warm, low, productive environments poorly.
- Many ecoregion × class strata will have focal pixels and no reference pixels. This is a result, not an error, but it is reported separately from the absence of an environmental analog.
- Matching with replacement lets a few protected pixels serve very many focal pixels.
- Thresholds are researcher choices and the proportions will move with them.
- Land cover (2021), human modification (2020) and protected areas (current release filtered to 2021) differ slightly in reference year. Climate normals are 1970 to 2000.
- The global caliper is calibrated mostly from data-rich ecoregions, where reference pixels are close to each other. In sparse ecoregions it may be strict.
- Whitening gives the lowest retained components the same weight as PC1.

## 11. Remaining open items

- Verify, once data are downloaded:
  - the WDPA attribute schema (field names changed in recent releases; the script accepts both old and new names);
  - that ECO_ID 0 is "Rock and Ice" in Ecoregions 2017, since the script drops it;
  - the gHM value range.
- Repeat the literature search before submission.

## 12. Pipeline outline

The R script follows these steps; each caches its output and is skipped on re-run.

0. Download all layers.
1. Template grid and layer alignment (climate, topography, gHM, ecoregions, countries).
2. Protected area rasterization: union fraction and largest-coverage PA identity, for the strict and all sets.
3. WorldCover tiles fetched one at a time, aggregated to 1 km class fractions, then deleted.
4. Natural fraction and modal natural class.
5. Environmental PCA and whitening.
6. Independent patches of the baseline pool.
7. Tabulate cells to disk, partitioned by ecoregion.
8. Distance from each cell to the nearest protected cell.
9. Caliper calibration, then per-stratum matching for each run.
10. Area-weighted summaries by country, ecoregion and biome, cross-country reliance, and per-PA importance.
11. Baseline and one-at-a-time sensitivity runs.
