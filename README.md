# BeetleBodySizeVariation

Quantifying intraspecific and community body size variation in carabid beetles across NEON sites, using elytra length measured from specimen images.

> **Status: work in progress.** This repository is under active development and is not yet associated with a published manuscript. Scripts are research code: they run top to bottom in an interactive session rather than as a packaged workflow, and several paths are hard-coded to my local machine (see [Configuration](#configuration)).

## Overview

I combine elytra length measurements from four sources into a single harmonized dataset (`BodysizeCombinedClean.csv`), then ask three questions of it:

1. **What shape are body size distributions?** I characterize skewness, kurtosis, and modality (Hartigan's dip test) of species-level size distributions.
2. **How does body size variance partition across spatial scales?** I compute CV² as a percentage of the mean (`100 * var / mean²`) for each species at the species, domain, site, and plot level, and compare nested scales.
3. **How much do co-occurring species overlap in size?** I compute community trait overlap statistics (O-statistics) at both site and plot scales, and relate them to species richness and environment.
4. **Is the observed size structure more (or less) structured than chance?** I test each community against a set of custom null models drawn from its regional pool, using five complementary metrics (overlap, niche range, and species-spacing statistics). See [Null-model framework](#null-model-framework).

Because the trait densities are estimated from an *augmented* sample (sparse groups topped up to n = 20), observation counts no longer track true abundance. Overlap statistics are therefore weighted by effort-scaled abundances taken from the standardized NEON pitfall data, not by how many specimens were imaged. See [Abundance weighting](#abundance-weighting).

Body size is measured as maximum elytra length in cm. Analyses are run on `log10` elytra length where distributions are compared or overlap is computed.

## Main pipeline

The scripts below are the current analysis. They are not sourced by a master script; I run them in this order.

| Order | Script | What it does | Key outputs |
|---|---|---|---|
| 1 | `CombineAndCleanDatasets.R` | Merges all four measurement sources, applies QC flags (z-score, robust z-score, image consistency, species extremes), integrates manual review, and writes the harmonized dataset. | `./Data/BodysizeCombinedClean.csv` |
| 2 | `BodySizeQuantification.R` | Reads `BodysizeCombinedClean.csv`. Computes distribution shape (skew, kurtosis, dip test) and CV² across nested spatial scales. Produces latitude and Bergmann's rule figures, and writes the per-scale CV² summary the overlap scripts read for augmentation. | `./Figures/BodySizeQuantification/*.png`, `./Outputs/CVpctSummary.csv` |
| 3 | `CalculateAbundance.R` | Independently of steps 1–2, reads the standardized NEON carabid counts from `neonDivData` and writes effort-scaled abundances per site and per plot (2018–2019). These weight the overlap statistics; image counts are not used as abundances because they are skewed by curation and collection. | `./Data/site_abund.csv`, `./Data/plot_abund.csv` |
| 4 | `Overlap_Site.R` | Site-level O-statistics under four sampling treatments (see [Sampling treatments](#sampling-treatments)). | `./Outputs/Site_Ostats_{unedited,20plus,augmented,augmented3to19}.csv`, `./Figures/Overlap/*.png` |
| 5 | `Overlap_Plot.R` | Plot-level O-statistics under the same four treatments, plus per-site overlap panels, overlap maps, and effect size maps. | `./Outputs/Plot_Ostats_*.csv`, `./Figures/Overlap/Plots/`, `./Figures/Overlap/Plots/Maps/` |
| 6 | `Overlap_CustomNulls.R` | Tests each community against three custom null models × five metrics, at a focal scale and regional pool set by two variables at the top of the script. Sources `community_overlap_weighted.R`. See [Null-model framework](#null-model-framework). | `./Outputs/{level}_by_{pool}_PoolNull.csv`, `./Outputs/{level}_by_{pool}_IndividualNull.csv`, `./Outputs/{level}_SwapMeansNull.csv` |
| 7 | `Environmental.R` | Assembles plot-level environmental covariates: WorldClim bioclim variables plus a PCA, MODIS NPP, and LiDAR-derived canopy structure and rugosity. | `./Outputs/BeetlePlotswEnvData.csv`, `./Outputs/BETplot_Rugosity` |
| 8 | `OverlapRichness.R` | Joins O-statistics to estimated species richness and site environment. This is the newest and least developed script. | Exploratory only |

Most analysis scripts read `Data/BodysizeCombinedClean.csv` directly, but the dependencies are no longer completely flat:

- Steps 4–6 read `Outputs/CVpctSummary.csv` to set the per-scale augmentation CV, so **step 2 must run before them**.
- Step 6 additionally reads the abundance files from **step 3** to weight its overlaps.
- Step 3 is self-contained (it pulls straight from `neonDivData`) and can run any time.

### Sampling treatments

Species-site and species-plot combinations vary enormously in sample size, and O-statistics are sensitive to small n. I therefore compute overlap four ways and compare:

| Treatment | Definition |
|---|---|
| `unedited` | All observations, no filtering |
| `20plus` | Only species-site (or species-plot) combinations with n ≥ 20 |
| `augmented` | All combinations with n < 20 padded up to n = 20 with simulated observations |
| `augmented3to19` | Only combinations with 3 ≤ n < 20 padded up to n = 20; singletons and doubletons dropped |

Augmentation draws simulated lengths from a lognormal distribution parameterized by the observed group mean and a typical CV² for that spatial scale. The CV² is no longer hard-coded: each overlap script reads `Outputs/CVpctSummary.csv` (written by `BodySizeQuantification.R`) and picks the row matching its scale, so plot-level augmentation uses the plot-level CV², site-level uses the site-level CV², and so on. The lognormal keeps simulated lengths positive and preserves observed mean-variance scaling. Simulation uses `set.seed(42)`; O-statistics use `random_seed = 517`.

### Null-model framework

`Overlap_CustomNulls.R` asks whether an observed community's body-size structure departs from what a random draw of the regional pool would produce. It is written as a flat, stepwise script: run the setup sections once, then run either null on its own. Two variables at the top set the scale and pool, and the script is meant to be re-run across pairings to trace the nested-scale trend:

- `LEVEL` — the focal community: `"plot"` or `"site"`.
- `POOL` — the regional pool the nulls draw from: `"site"` or `"domain"`. Typical pairings are plot → site and site → domain.

Three null models are evaluated, each against five metrics (observed value, null 2.5/97.5% CI, SES, and a direction flag per metric):

| Null | What it randomizes | Question it answers |
|---|---|---|
| **Pool null** (§5) | Random species assemblage from the regional pool, holding richness and total N constant | The core trait-assembly test — is overlap lower (overdispersion) or higher (clustering) than a random pool draw? |
| **Individual null** (§6) | Keeps the observed species and abundances, but redraws each species' individuals from the pool | Does a different regional sample of the *same* species change the result? (the individual-level-data test) |
| **Swap-means null** (§7) | Permutes species' community means within the focal community (`Ostats` `swap_means`); needs no pool | A within-community reference; informative for the overlap metrics only |

The five metrics are `overlap_norm` (densities normalized to area 1), `overlap_unnorm` (densities scaled by abundance), `niche_range` (2.5–97.5% span of occupied trait space), `sdnnd` (SD of nearest-neighbour distances between species means — low means even spacing), and `min_logratio` (smallest adjacent gap between species means on the log10 axis, compared to log10(1.3) as a Hutchinsonian reference). The two spacing metrics are invariant under the swap-means null by construction and correctly report a neutral result there.

Focal units whose pool holds fewer than `MIN_POOL_UNITS` (default 2) focal units are dropped, so single-plot sites and single-site domains are excluded from the pool draws.

### Abundance weighting

`community_overlap_weighted.R` is a drop-in replacement for `Ostats::community_overlap()` that weights by an external abundance vector rather than by the number of observations in each group. It is needed because augmentation tops sparse species up to n = 20, so observation counts no longer reflect true abundance — and `community_overlap()` uses those counts both for the area under each species' density (when unnormalized) and for the harmonic-mean weights on the pairwise overlaps. This version takes the effort-scaled abundances from `CalculateAbundance.R` and uses them in both places; densities are still estimated from the augmented trait values, only the weighting changes. It reproduces `community_overlap()` to numerical precision when the supplied weights equal the observation counts. Species with no true abundance are dropped from every metric so each focal unit describes one consistent community, and that coverage is reported when the script runs.

## Supporting and exploratory scripts

These are kept for provenance and for figure generation. They are not part of the main analysis path.

| Script | Purpose |
|---|---|
| `community_overlap_weighted.R` | Abundance-weighted replacement for `Ostats::community_overlap()`. Not run on its own; sourced by `Overlap_CustomNulls.R`. See [Abundance weighting](#abundance-weighting). |
| `TPDexample.R` | Conceptual trait probability density figures using D01 species. Illustration for talks and manuscript concept figures, not analysis. |
| `Summary.Rmd` | Sample size and coverage report comparing imaged specimens to the full NEON pinned collection (species imaged, species remaining, occurrences per species). Written for OSC and points at `/fs/ess/PAS2136/CarabidImaging/`. |
| `Overlap_explore.R` | Early O-statistics exploration on the BeetlePalooza data alone. Superseded by `Overlap_Site.R` and `Overlap_Plot.R`. |
| `FirstCVOutputs.R` | First look at raw computer vision keypoint output. Parses `pred_coords_px` into length and width. Superseded by the cleaned length files. |

## Data sources

| Source | File / product | Notes |
|---|---|---|
| NEON Biorepository CV measurements | `Data/beetle_lengths_cm_reviewed_clean.csv`, `Data/allIndividuals.csv` | Merged on `beetle_id` / `individualID`. Flagged records dropped; duplicates reduced to first record per beetle. |
| BeetlePalooza annotations | `Data/BeetleMeasurements.csv` | Subset to `user_name == "IsaFluck"` and `structure == "ElytraLength"` for annotator consistency. |
| Hawaii (PUUM) annotations | `Data/trait_annotations.csv` | Includes scalebar columns (`px_scalebar`, `cm_scalebar`). |
| NEON carabid pitfall data | `DP1.10022.001` via `neonUtilities` | Parataxonomist and expert taxonomist tables combined, expert IDs taking precedence. Released and provisional records both pulled, then cached to `Data/NEON_ExpertParaCombined.csv` and `Data/NEON_ExpertParaCombined_Prelim.csv`. |
| NEON effort-scaled abundances | `data_beetle` via `neonDivData` | Standardized carabid counts, 2018–2019, summed to species per site and per plot. Used to weight the overlap statistics (`CalculateAbundance.R` → `site_abund.csv`, `plot_abund.csv`). Image-derived counts are deliberately not used here because they are skewed by curation and collection. |
| NEON site metadata | `Data/NEON_Field_Site_Metadata_20260130.csv` | Committed to this repository. Downloaded 2026-01-30. |
| WorldClim bioclim | `worldclim_global`, 0.5 arc-min, USA | 19 bioclim variables, extracted at BET plot centroids, log-transformed for precipitation, then PCA (`princomp`, correlation matrix). |
| MODIS NPP | `NEON_MODIS_NPP_2018_2019.csv` | Derived in Google Earth Engine. Script link is in `Environmental.R`. |
| NEON LiDAR | `DP1.30003.001` (discrete return), `DP3.30024.001` (DTM) | 2018 preferred, falling back to 2019 then 2017. Clipped to 40 m buffers around plot centroids; 0.5 m voxel density used to derive rugosity, mean height, 98th percentile height, and height SD. |
| NEON plot geometry | `All_NEON_TOS_Plot_Polygons_V11_BET.shp`, `All_NEON_TOS_Plot_Centroids_V11_BET.shp` | Expected one level above the repository root. |
| Species richness estimates | `../BeetleBiodiversity/{Site,plot}_annualVarWeightedMean_EstimatedSppRichness.csv` | Produced by a separate repository. Required only by `OverlapRichness.R`. |

Harmonization detail worth knowing: scientific names are reduced to genus plus species by stripping parenthetical subgenera, collapsing whitespace, and cutting anything after a slash. Records identified only to genus (`sp.`) are excluded from distribution-shape and variance analyses. Non-positive lengths are dropped.

### Data availability

<!-- TODO: Complete before making this repository public. -->

**Placeholder.** Raw measurement files are not currently committed to this repository. Statements about where each dataset lives, under what license, and how to obtain it will be added here prior to publication.

## Requirements

R, plus:

**Core:** `dplyr`, `tidyr`, `stringr`, `ggplot2`, `ggpubr`, `ggrepel`

**Analysis:** `Ostats`, `moments`, `diptest`, `psych`, `parameters`, `sfsmisc`, `matrixStats` (the last two are called directly by `community_overlap_weighted.R`)

**Spatial and remote sensing:** `sf`, `terra`, `maps`, `geodata`, `lidR`, `gstat`, `purrr`

**NEON:** `neonUtilities`, `neonOS`, `neonDivData`

**Reporting:** `readxl`, `tibble`, `gridExtra`, `corrplot`

## Configuration

Before running anything:

1. **Working directory.** Every script calls `setwd()` with an absolute path (`/home/aly/Beetles/BeetleBodySizeVariation`; `Summary.Rmd` uses an OSC path). Change these to your own path.
2. **NEON API token.** Scripts read a token from `~/NEON_TOKEN` as a headerless single-value file. Request a token at the NEON data portal.
3. **Large data paths.** `Environmental.R` writes WorldClim and LiDAR downloads to an external drive (`/media/aly/Penobscot/...`). Change to a local path with sufficient space.
4. **Directory structure.** Scripts assume these exist and will not create them:

```
BeetleBodySizeVariation/
├── Data/
├── Outputs/
└── Figures/
    ├── BodySizeQuantification/
    └── Overlap/
        └── Plots/
            └── Maps/
```

5. **Cached NEON pulls.** The NEON download blocks are guarded by `file.exists()`. Delete the cached CSVs in `Data/` to force a fresh pull.

## Known rough edges

<!-- TODO: Resolve or remove before making this repository public. -->

- Absolute paths are hard-coded throughout rather than being relative or configured once.
- `Environmental.R` writes `./Outputs/BETplot_Rugosity` without a `.csv` extension.
- LiDAR downloads are gated behind `p <- NA` so they do not fire accidentally. Set `p` to `1` or `2` to enable.
- Several NEON sites download AOP data under a neighboring site code (DCFS/WOOD, KONA/KONZ, TREE/STEI, STEI/CHEQ). `Environmental.R` remaps these explicitly.
- The Biorepo source has no pixel or scalebar columns, so `px_scalebar`, `cm_scalebar`, and `px_elytra_max_length` are set to `NA` during harmonization and flagged in the code for a future update.
- `Overlap_Site.R` and `Overlap_Plot.R` select the CV² for augmentation by hard-coded row number (`cv_reslults[3, ]` / `cv_reslults[4, ]`), so they silently break if the row order of `CVpctSummary.csv` changes. `Overlap_CustomNulls.R` already does this the safe way, matching on the `scale` column; the two older scripts should be brought in line (and the `cv_reslults` misspelling cleaned up).
- `BodySizeQuantification.R` uses a single user-defined `cutoff <- 50` for minimum observations per group, but several figures apply a different threshold (n ≥ 20) inline.

## Citation

<!-- TODO: Add manuscript citation and archived release DOI when available. -->

**Placeholder.** No associated manuscript yet.

## Acknowledgements

<!-- TODO: Add funding sources, collaborators, and data provider acknowledgements. -->

**Placeholder.**

This work uses data from the National Ecological Observatory Network (NEON), a program sponsored by the U.S. National Science Foundation and operated under cooperative agreement by Battelle. NEON data are provided under their own terms of use, which are separate from the license on this code.

## License

Code in this repository is released under the MIT License. See [LICENSE](LICENSE).

## Contact

Alyson East, PhD candidate in Quantitative Ecology, University of Maine.
