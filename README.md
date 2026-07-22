# BeetleBodySizeVariation

Quantifying intraspecific and community body size variation in carabid beetles across NEON sites, using elytra length measured from specimen images.

> **Status: work in progress.** This repository is under active development and is not yet associated with a published manuscript. Scripts are research code: they run top to bottom in an interactive session rather than as a packaged workflow, and several paths are hard-coded to my local machine (see [Configuration](#configuration)).

## Overview

I combine elytra length measurements from four sources into a single harmonized table of individual beetles (`all_elytra`), then ask three questions of it:

1. **What shape are body size distributions?** I characterize skewness, kurtosis, and modality (Hartigan's dip test) of species-level size distributions.
2. **How does body size variance partition across spatial scales?** I compute CV² as a percentage of the mean (`100 * var / mean²`) for each species at the species, domain, site, and plot level, and compare nested scales.
3. **How much do co-occurring species overlap in size?** I compute community trait overlap statistics (O-statistics) at both site and plot scales, and relate them to species richness and environment.

Body size is measured as maximum elytra length in cm. Analyses are run on `log10` elytra length where distributions are compared or overlap is computed.

## Main pipeline

The scripts below are the current analysis. They are not sourced by a master script; I run them in this order.

| Order | Script | What it does | Key outputs |
|---|---|---|---|
| 1 | `BodySizeQuantification.R` | Harmonizes all four data sources into `all_elytra`. Computes distribution shape (skew, kurtosis, dip test) and CV² across nested spatial scales. Produces latitude and Bergmann's rule figures. | `./Figures/BodySizeQuantification/*.png` |
| 2 | `Overlap_Site.R` | Site-level O-statistics under four sampling treatments (see [Sampling treatments](#sampling-treatments)). | `./Outputs/Site_Ostats_{unedited,20plus,augmented,augmented3to19}.csv`, `./Figures/Overlap/*.png` |
| 3 | `Overlap_Plot.R` | Plot-level O-statistics under the same four treatments, plus per-site overlap panels, overlap maps, and effect size maps. | `./Outputs/Plot_Ostats_*.csv`, `./Figures/Overlap/Plots/`, `./Figures/Overlap/Plots/Maps/` |
| 4 | `Environmental.R` | Assembles plot-level environmental covariates: WorldClim bioclim variables plus a PCA, MODIS NPP, and LiDAR-derived canopy structure and rugosity. | `./Outputs/BeetlePlotswEnvData.csv`, `./Outputs/BETplot_Rugosity` |
| 5 | `OverlapRichness.R` | Joins O-statistics to estimated species richness and site environment. This is the newest and least developed script. | Exploratory only |

`Overlap_Site.R` and `Overlap_Plot.R` each repeat the data harmonization block from `BodySizeQuantification.R` rather than sourcing it, so they can be run standalone.

### Sampling treatments

Species-site and species-plot combinations vary enormously in sample size, and O-statistics are sensitive to small n. I therefore compute overlap four ways and compare:

| Treatment | Definition |
|---|---|
| `unedited` | All observations, no filtering |
| `20plus` | Only species-site (or species-plot) combinations with n ≥ 20 |
| `augmented` | All combinations with n < 20 padded up to n = 20 with simulated observations |
| `augmented3to19` | Only combinations with 3 ≤ n < 20 padded up to n = 20; singletons and doubletons dropped |

Augmentation draws simulated lengths from a lognormal distribution parameterized by the observed group mean and a typical CV² estimated from well-sampled groups (`typical_cvpct <- 0.47`). The lognormal keeps simulated lengths positive and preserves observed mean-variance scaling. Simulation uses `set.seed(42)`; O-statistics use `random_seed = 517`.

## Supporting and exploratory scripts

These are kept for provenance and for figure generation. They are not part of the main analysis path.

| Script | Purpose |
|---|---|
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
| NEON site metadata | `NEON_Field_Site_Metadata_20260130.csv` | Committed to this repository. Downloaded 2026-01-30. |
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

**Analysis:** `Ostats`, `moments`, `diptest`, `psych`, `parameters`

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
- The LiDAR metrics loop in `Environmental.R` starts at `i = 45` as a resume point from an interrupted run. Reset to `1` for a full run.
- LiDAR downloads are gated behind `p <- NA` so they do not fire accidentally. Set `p` to `1` or `2` to enable.
- Several NEON sites download AOP data under a neighboring site code (DCFS/WOOD, KONA/KONZ, TREE/STEI, STEI/CHEQ). `Environmental.R` remaps these explicitly.
- The Biorepo source has no pixel or scalebar columns, so `px_scalebar`, `cm_scalebar`, and `px_elytra_max_length` are set to `NA` during harmonization and flagged in the code for a future update.
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
