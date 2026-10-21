# BeetleBodySizeVariation

Quantifying intraspecific and community body size variation in carabid beetles (Coleoptera: Carabidae) across NEON sites, and testing how body-size niche space links environmental gradients to species richness.

> **Status: work in progress.** This repository supports a manuscript in preparation; there is no published paper or archived release yet. Scripts are research code: they are flat and stepwise, written to run top to bottom in an interactive R session (or section by section), and several paths are hard-coded to my machine (see [Configuration](#configuration)).

## Overview

I measure maximum elytra length (cm) from NEON carabid specimen images, harmonize measurements from four sources into one dataset (`Data/BodysizeCombinedClean.csv`), and analyze `log10` elytra length across the NEON spatial hierarchy (individual, plot, site, domain). The analysis covers 2018 and 2019 collections. <!-- TODO: confirm final counts (individuals, species, sites) against BodysizeCombinedClean.csv before release. -->

I ask four questions:

1. **How is body size distributed and partitioned?** Distribution shape (skew, kurtosis, Hartigan's dip test) and coefficient of variation across nested spatial scales.
2. **Does intraspecific body size vary in space, and why?** Tukey HSD tests among sites and among plots within sites, then per-species mixed models relating local body size to environment and competition (character displacement terms).
3. **Is community body-size structure more or less structured than chance?** Seven community-shape metrics tested against three custom null models, at plot and site scales and across regional pool extents.
4. **Does body-size niche space mediate the environment-richness relationship?** Structural equation models (lavaan) in which climate, productivity, climate velocity, and geodiversity act on richness directly and through niche breadth, niche depth, and spacing evenness.

## Repository structure

```
BeetleBodySizeVariation/
├── 1_CombineAndCleanDatasets.R        # numbered scripts = main pipeline
├── 2a_ ... 9b_*.R
├── 5_run_byYear.sh                    # shell driver for the null-model grid
├── community_depth.R                  # helper, sourced by 6a
├── community_overlap_weighted.R       # legacy helper (pairwise overlap)
├── Archive/                           # superseded scripts, kept for provenance
├── SKILL.md                           # documentation conventions for this repo
├── LICENSE
└── README.md
```

Scripts create and expect `Data/`, `Outputs/`, `Figures/`, and `logs/` in the repository root. These folders hold data and results and are not committed.

## Getting started

### 1. Clone

```bash
git clone https://github.com/AlysonEast/BeetleBodySizeVariation.git
cd BeetleBodySizeVariation
```

### 2. Install dependencies

<!-- TODO: add renv.lock (renv::snapshot()) and replace this note with `renv::restore()`. -->

**Placeholder.** I have not yet committed a machine-readable requirements file. Until I do, the packages each script needs are its `library()` calls at the top of the file. The analysis relies mainly on `dplyr`, `ggplot2`, `ggpubr`, `lme4`/`lmerTest`, `lavaan`, `ape`, `picante`, `neonUtilities`, `neonDivData`, `sf`, `terra`, and `lidR`.

### 3. Configuration

Before running anything:

1. **Working directory.** Every R script calls `setwd("/home/aly/Beetles/BeetleBodySizeVariation")`, and `5_run_byYear.sh` sets `REPO` to the same path. Change these to your clone.
2. **NEON API token.** Scripts that pull NEON data read a token from `~/NEON_TOKEN` (headerless, single value). Request one from the NEON data portal.
3. **External drive paths.** `4_Environmental.R` and `9a_SiteOverlapRichness_LDGPaths.R` read and write large rasters under `/media/aly/Penobscot/...` (WorldClim, climate velocity, geodiversity, LiDAR). Point these at a local path with enough space.
4. **Files one level above the repo.** Plot geometry shapefiles, the NEON sampling boundary shapefile, and the MODIS NPP tables are read from `../`. The richness estimates are read from the sibling repository `../BeetleBiodiversity/`.
5. **Cached NEON pulls.** NEON downloads in `1_CombineAndCleanDatasets.R` are guarded by `file.exists()`. Delete the cached CSVs in `Data/` to force a fresh pull. LiDAR downloads in `4_Environmental.R` are gated behind `p <- NA`; set `p` to enable them.

## Main pipeline

The number prefix gives the order. Scripts are not called by a master script; I run them individually.

| Step | Script | What it does | Key outputs |
|---|---|---|---|
| 1 | `1_CombineAndCleanDatasets.R` | Merges the four measurement sources, applies QC flags, and integrates manual review decisions. | `Data/BodysizeCombinedClean.csv` |
| 2a | `2a_BodySizeQuantification.R` | Distribution shape and CV by spatial scale; latitude and community-weighted-mean figures. | `Outputs/CVpctSummary.csv`, `Outputs/SpeciesList.csv`, `Figures/BodySizeQuantification/` |
| 2b | `2b_BodySizeSpatialTests.R` | One-way ANOVA + Tukey HSD of intraspecific body size among sites and among plots within sites, BH-corrected across species; nested density figures with compact letter displays. | `Outputs/Top*Spp_MeanDiffs_*.docx`, `Figures/BodySizeQuantification/SpatialTests/` |
| 2c | `2c_BodySizeAbundantSpecies.R` | Per-species models (2018) of individual body size on environment and competition terms, with and without a site random effect, plus variance partitioning. | `Outputs/IntraspecificShift_TopSpp_*`, `Figures/BodySizeQuantification/IntraspecificShift/` |
| 3 | `3_CalculateAbundance_ByYear.R` | Effort-scaled NEON pitfall abundances from `neonDivData`, per species and in total, per site and plot, per year. | `Data/{site,plot}_abund_{YEAR}.csv`, `Data/{site,plot}Total_abund_{YEAR}.csv` |
| 4 | `4_Environmental.R` | Plot-level covariates: WorldClim PCA, SoilGrids pH, climate velocity, and NEON AOP LiDAR canopy structure, rugosity, and terrain variability. | `Outputs/BeetlePlotswEnvData.csv`, `Outputs/BETplot_Rugosity.csv`, `Outputs/BETplot_Soils.csv`, `Outputs/Beetle{Plots,Sites}wVelocity.csv` |
| 5 | `5_run_byYear.sh` | Runs step 6a across every year and level/pool pairing in parallel. | one set of 6a outputs per run, `logs/` |
| 6a | `6a_CommunityShape_CustomNulls_ByYear.R` | Community-shape metrics under three null models for one `LEVEL`, `POOL`, and `YEAR`. Sources `community_depth.R`. | `Outputs/{LEVEL}_by_{POOL}_{YEAR}_{Pool,Individual,SwapMeans}Null.csv` |
| 6b | `6b_PhyloDistance_ByYear.R` | PD (with SES via `picante::ses.pd`), MPD, and MNTD for every plot-year and site-year community. | `Outputs/{plot,site}_PhyloDistance.csv` |
| 7 | `7_CommunityShape_CustomNulls_Compare.R` | Does the answer change with the null, year, pool extent, or scale? Pairwise comparison figures and an agreement table. | `Outputs/NullComparison_Agreement_CommunityShape.csv`, `Figures/NullComparisons/` |
| 8 | `8_CommunityShape_CustomNulls_DirPatterns.R` | Significance-flag composition, null agreement, and observed and SES maps. | `Outputs/DirPatterns_*.csv`, `Figures/DirPatterns/` |
| 9a | `9a_SiteOverlapRichness_LDGPaths.R` | Site-level SEM of richness. | `Outputs/finalsiteDF.csv`, `Figures/SEMs/sites*.png` |
| 9b | `9b_PlotOverlapRichness_LDGPaths.R` | Plot-level SEM of richness (bootstrap SEs clustered by site). | `Outputs/finalplotDF.csv`, `Figures/SEMs/plots*.png` |

**Dependencies that do not follow the numbering.** Steps 3 and 4 are self-contained and can run first. Step 2c needs outputs from steps 3 and 4, the phylogeny, and the `BeetleBiodiversity` richness estimates. Steps 5 and 6b need step 3. Steps 7 and 8 need the full 6a grid. Steps 9a and 9b need steps 3, 4, 6a, and 6b plus the `BeetleBiodiversity` richness estimates.

### Example: running the null-model grid

```bash
# full grid (2 years x 5 level/pool pairings), 16 jobs at a time
./5_run_byYear.sh

# fewer cores, one year
JOBS=4 YEARS="2018" ./5_run_byYear.sh

# or a single run directly
Rscript 6a_CommunityShape_CustomNulls_ByYear.R site all 2018
```

## Reference

Compact definitions. Each script's header comment holds the full rationale.

### Community-shape metrics and null models (`6a`)

`LEVEL` sets the focal community (`plot` or `site`); `POOL` sets the regional pool the nulls draw from (`site`, `domain`, or `all`). Valid pairings are plot by site, domain, or all, and site by domain or all.

| Null | What it randomizes |
|---|---|
| **Pool** | A random assemblage from the regional pool, holding richness and total N constant. The core assembly test. |
| **Individual** | Keeps observed species and abundances; redraws each species' individuals from its pool-level individuals. The individual-level-data test. |
| **Swap-means** | Keeps within-species deviations but permutes species onto each other's community means. Needs no pool. |

Each null reports, per metric, the observed value, the null 2.5% and 97.5% bounds, the SES, and a direction flag (`lower`, `neutral`, `higher`). Metrics: `overlap_depth`, `niche_range`, `cwm`, `cw_variance`, `cw_skew`, `cw_kurtosis`, and `sdnnd`. Steps 7 and 8 focus on `niche_range` and `overlap_depth`.

`overlap_depth` is computed by `community_depth.R`: peak-scaled species densities on a common grid, summed and averaged over the occupied range. It replaces the pairwise O-statistic, which becomes zero-inflated in species-rich communities. All metrics run on observed individuals with no augmentation; a species enters a community only if it has a true effort-scaled NEON abundance.

### Richness SEMs (`9a`, `9b`)

Exogenous predictors are mean annual temperature, MODIS NPP, climate velocity, and geodiversity. Three community mediators carry the body-size signal: breadth (`niche_range_obs`), depth (`overlap_depth_obs`), and evenness (`sdnnd_obs`, inverted so higher values mean more even spacing). Each script fits the model with and without phylogenetic distance (MPD) as an additional mediator. Observed metrics, not SES, enter the models. Column choices sit in the `CONFIG` block (section 1) of each script.

## Supporting scripts

| Script | Purpose |
|---|---|
| `community_depth.R` | Defines `community_depth()`. Sourced by 6a; not run on its own. |
| `community_overlap_weighted.R` | Abundance-weighted replacement for `Ostats::community_overlap()`. No longer used by the active pipeline since I dropped augmentation and the pairwise metric. |
| `Archive/` | Superseded scripts (earlier O-statistic, augmentation, and candidate-model SEM versions), kept for provenance. |

## Data sources

| Source | File / product | Used for |
|---|---|---|
| NEON Biorepository computer vision measurements | `Data/beetle_lengths_cm_clean.csv`, `Data/allIndividuals.csv` | Elytra length (main source) |
| BeetlePalooza annotations | `Data/BeetleMeasurements.csv` | Elytra length, one annotator (`IsaFluck`) for consistency |
| Hawaii (PUUM) annotations | `Data/trait_annotations.csv` | Elytra length, Hawaii site |
| NEON ground beetle pitfall data | `DP1.10022.001` via `neonUtilities` | Taxonomy and collection metadata (expert IDs take precedence) |
| NEON effort-scaled abundances | `data_beetle` via `neonDivData` | Community membership and abundance weights |
| NEON field site metadata | `Data/NEON_Field_Site_Metadata_20260130.csv` (downloaded 2026-01-30) | Site attributes |
| NEON plot geometry | `All_NEON_TOS_Plot_{Polygons,Centroids}_V11_BET.shp`, `NEON_terrestrialSamplingBoundariesDissolve.shp` | Plot and site extraction |
| WorldClim 2.1 bioclim | `geodata::worldclim_global`, 0.5 arc-min | Plot-level climate and climate PCA |
| SoilGrids | pH (H2O), 0 to 5 cm | Plot-level soil covariate |
| Climate velocity | Dryad `doi:10.5061/dryad.b13j1` | Climate velocity predictor <!-- TODO: confirm and cite the source paper. --> |
| Geodiversity | Kapsar et al. (EDI package `edi.2320.1`) | Site-level geodiversity <!-- TODO: add full citation. --> |
| MODIS NPP | `NEON_MODIS_NPP_2018_2019.csv`, `NEONSites_MODIS_NPP_2018_2019.csv` | Productivity, derived in Google Earth Engine (script links in 2b, 9a, 9b) |
| NEON AOP LiDAR | `DP1.30003.001` (discrete return), `DP3.30024.001` (DTM) | Canopy rugosity, height metrics, terrain variability within 40 m plot buffers |
| Carabid phylogeny | `Data/Tree.DNA.droppedseqs.newick` | Phylogenetic distance (6b) and nearest-relative terms (2c) <!-- TODO: add source and citation. --> |
| Species richness estimates | `../BeetleBiodiversity/{site,plot}_annual_EstimatedSppRichness.csv` | SEM response; produced by the [BeetleBiodiversity](https://github.com/AlysonEast/BeetleBiodiversity) repository |

### Input format

Downstream scripts read `Data/BodysizeCombinedClean.csv`, one row per measured individual. The columns they rely on are `individualID`, `imageID`, `scientificName_Species` (genus + species, subgenera stripped), `cm_elytra_max_length`, `siteID`, `plotID`, `yearCollected`, and `datasource`. Records identified only to genus (`sp.`) are excluded from species-level analyses.

### Data availability

<!-- TODO: Complete before making this repository public. -->

**Placeholder.** Raw measurement files are not committed to this repository. The image datasets are published on Hugging Face ([imageomics/2018-NEON-beetles](https://huggingface.co/datasets/imageomics/2018-NEON-beetles), [imageomics/Hawaii-beetles](https://huggingface.co/datasets/imageomics/Hawaii-beetles)). I will add where each processed file lives, its license, and how to obtain it before publication.

## Known rough edges

<!-- TODO: Resolve or remove before making this repository public. -->

- Absolute paths are hard-coded throughout rather than set once.
- `5_run_byYear.sh` calls `6_CommunityShape_CustomNulls_ByYear.R`, but the script is now `6a_...`, so the driver fails until that line is updated. Header comments in 3, 5, and 6a also still use pre-renumbering script names.
- Headers in 9a and 9b still describe seven candidate models and `overlap_unnorm_obs`; the scripts now fit one hypothesis-driven model (with and without phylogeny) using the three mediators above.
- 9a and 9b write their model summaries (`site_report*.txt`, `plot_report*.txt`) to the repository root instead of `Outputs/`.
- `4_Environmental.R` labels the LiDAR voxel step "1 m voxel" but uses `res = 0.5`, and labels 0.5 arc-min WorldClim as "~900 km" (it is roughly 1 km).
- Site geodiversity uses the `bio01_sq` column of the Kapsar et al. footprint file. <!-- TODO: confirm this is the intended geodiversity metric. -->
- The Biorepository source has no scalebar columns, so `px_scalebar`, `cm_scalebar`, and `px_elytra_max_length` are `NA` for those records.
- `2a_BodySizeQuantification.R` sets `cutoff <- 50` for minimum observations, while several figures apply n >= 20 inline.

## Citation

<!-- TODO: Add CITATION.cff, manuscript BibTeX, and the version-agnostic Zenodo DOI at first release. -->

**Placeholder.** No associated manuscript or archived release yet. When available, please cite both the paper and this repository.

## Acknowledgements

<!-- TODO: Add funding (Imageomics, NSF Award #2118240, if applicable), collaborators, and data provider acknowledgements. -->

**Placeholder.**

This work uses data from the National Ecological Observatory Network (NEON), a program sponsored by the U.S. National Science Foundation and operated under cooperative agreement by Battelle. NEON data are provided under their own terms of use, separate from the license on this code.

## License

Code in this repository is released under the MIT License. See [LICENSE](LICENSE).

## Contact

Alyson East, PhD candidate in Quantitative Ecology, University of Maine.
