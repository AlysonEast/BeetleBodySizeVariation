---
name: beetle-project-docs
description: Write or update the README and other documentation for the NEON carabid body-size project (repos AlysonEast/BeetleBodySizeVariation and AlysonEast/BeetleBiodiversity). Use this whenever asked to review repo status and update the README, document a new script, or refresh the docs after changes. It carries the project's standing context and vocabulary so the README can stay a lean, scannable orientation instead of re-explaining every mechanism inline. Trigger this even when the user just says "update the README" or "the docs are out of date" for these repos.
---

# Beetle project docs

Documentation for the NEON carabid (ground beetle) body-size project. The point of this skill is to keep the README **lean**: a reader should be able to scan it in two minutes and know what the repo is, what each script does, and how to run it. Deep mechanism and rationale live where they belong — in the scripts, which are already heavily commented — not duplicated in the README.

## The project in three lines

Elytra length is measured from NEON specimen images and harmonized across four sources into `BodysizeCombinedClean.csv`. The analysis asks how beetle body size varies within species, partitions across spatial scales (individual → plot → site → domain), and overlaps among co-occurring species, and tests that overlap against null models. Work in progress; no manuscript yet; MIT license. Author: Alyson East, University of Maine.

## The lean-README rule

**Explain each mechanism once, in the place it lives.** The scripts carry full block comments explaining *why* each step is done. The README's job is to *name* the mechanism and *point* to the script — not to reproduce its reasoning.

When editing the README, apply this test to every paragraph: **does this sentence help someone decide what to run or where to look, or is it re-deriving something the code comment already explains?** If the latter, cut it to a pointer.

Concrete cuts to make:

- **Rationale paragraphs → one clause + pointer.** The "why" of a design choice belongs in the script header. The README states the choice and where to read more.
- **Algorithm walk-throughs → the name and the output.** Don't narrate what a function does step by step; name it, say what it produces, link the script.
- **Repeated definitions → define once.** Terms in the [Standing reference](#standing-reference) (treatments, metrics, nulls) get one compact definition, not a re-explanation each time they appear.
- **Parameter listings → only the ones a user must set.** `setwd()`, tokens, and the `LEVEL`/`POOL` knobs matter. Internal seeds and grid constants do not need prose.

**Example — abundance weighting.**

Over-explained (what to avoid):
> `community_overlap_weighted.R` is a drop-in replacement for `Ostats::community_overlap()` that weights by an external abundance vector rather than by the number of observations in each group. It is needed because augmentation tops sparse species up to n = 20, so observation counts no longer reflect true abundance — and `community_overlap()` uses those counts both for the area under each species' density (when unnormalized) and for the harmonic-mean weights on the pairwise overlaps. This version takes the effort-scaled abundances... [continues for several more sentences]

Lean (what to write):
> Overlap is weighted by effort-scaled NEON abundances, not by imaged-specimen counts, because augmentation inflates those counts. `community_overlap_weighted.R` (sourced by `Overlap_CustomNulls.R`) handles this; see its header comment for the details.

The lean version loses nothing a scanner needs — the rationale is one clause, and the full reasoning is one file-open away.

## README structure

Keep these sections, in this order. Add only what a section genuinely needs.

1. **Title + one-line description**
2. **Status** — one blockquote: work in progress, no manuscript, research code with hard-coded paths.
3. **Overview** — the questions the project asks (currently four), 1–2 lines each. No mechanism detail here.
4. **Main pipeline** — a table: order, script, one-line "what it does", key outputs. Below it, a short note on non-obvious run dependencies only.
5. **Standing reference subsections** — treatments, null framework, abundance weighting. Compact; each is a definition + pointer, per the rule above.
6. **Supporting scripts** — table of non-pipeline scripts.
7. **Data sources** — table.
8. **Requirements / Configuration** — packages; the paths and tokens a user must change.
9. **Known rough edges** — honest TODO list.
10. **Citation / Acknowledgements / License / Contact** — placeholders where not yet available.

## Standing reference

The canonical facts the README should name but not re-explain at length. This is also the context to write *from* — you don't need to re-read every script to update a section.

**Repos.** `BeetleBodySizeVariation` (R analysis; the active repo) and `BeetleBiodiversity` (richness/rarefaction; supplies richness estimates to `OverlapRichness.R`).

**Spatial hierarchy.** individual → plot → site → domain (NEON structure).

**Body size.** Maximum elytra length in cm; analyses run on `log10` length.

**Pipeline scripts (BeetleBodySizeVariation).**

| Script | Role | Key output |
|---|---|---|
| `CombineAndCleanDatasets.R` | Merge 4 sources, QC, manual review | `Data/BodysizeCombinedClean.csv` |
| `BodySizeQuantification.R` | Distribution shape + CV² by scale | figures, `Outputs/CVpctSummary.csv` |
| `CalculateAbundance.R` | Effort-scaled abundances from `neonDivData` | `Data/{site,plot}_abund.csv` |
| `Overlap_Site.R` / `Overlap_Plot.R` | O-statistics, 4 treatments, `swap_means` null | `Outputs/{Site,Plot}_Ostats_*.csv` |
| `Overlap_CustomNulls.R` | Null-model framework (below); sources the helper | `Outputs/{level}_by_{pool}_*Null.csv` |
| `community_overlap_weighted.R` | Abundance-weighted `community_overlap()`; sourced, not run | — |
| `Environmental.R` | WorldClim PCA, MODIS NPP, LiDAR structure | `Outputs/BeetlePlotswEnvData.csv` |
| `OverlapRichness.R` | Join overlap × richness × environment (exploratory) | — |

**Non-obvious dependencies.** Overlap scripts read `CVpctSummary.csv`, so `BodySizeQuantification.R` runs before them. `Overlap_CustomNulls.R` also needs the abundance files from `CalculateAbundance.R`. `CalculateAbundance.R` is self-contained.

**Sampling treatments.** `unedited` (all), `20plus` (n ≥ 20 only), `augmented` (n < 20 padded to 20), `augmented3to19` (3 ≤ n < 20 padded; singletons/doubletons dropped). Augmentation is lognormal from the group mean and the per-scale CV² read from `CVpctSummary.csv` (no longer hard-coded). Seeds: `set.seed(42)` for simulation, `517` for O-statistics.

**Null-model framework** (`Overlap_CustomNulls.R`). Set `LEVEL` (plot/site) and `POOL` (site/domain) at the top; re-run across pairings for the nested-scale trend. Three nulls: **pool** (random assemblage from the regional pool — the core assembly test), **individual** (redraw each species' individuals from the pool — the individual-level-data test), **swap-means** (permute community means within the community; needs no pool). Five metrics per null (observed, null CI, SES, direction): `overlap_norm`, `overlap_unnorm`, `niche_range`, `sdnnd`, `min_logratio`. The two spacing metrics are invariant under swap-means by construction.

**Abundance weighting.** Overlaps are weighted by effort-scaled NEON pitfall abundances, not imaged-specimen counts (which augmentation inflates). Handled by `community_overlap_weighted.R`.

**Data sources.** NEON Biorepository CV measurements, BeetlePalooza annotations (IsaFluck / ElytraLength), Hawaii (PUUM) annotations, NEON pitfall data (`DP1.10022.001` via `neonUtilities`), NEON abundances (`data_beetle` via `neonDivData`), NEON site metadata, WorldClim, MODIS NPP, NEON LiDAR (`DP1.30003.001`, `DP3.30024.001`), plot geometry, and richness estimates from the sibling repo.

**Config a user must change.** Absolute `setwd()` paths; `~/NEON_TOKEN`; large-download paths in `Environmental.R`; the expected `Data/ Outputs/ Figures/` directory tree. NEON pulls are cached and guarded by `file.exists()`.

## Voice and conventions

- First person, matter-of-fact, honest about rough edges. Keep the existing tone.
- Tables for anything enumerable (scripts, treatments, data sources). Prose only for connective tissue.
- Mark pre-publication gaps with `<!-- TODO: ... -->` and a **Placeholder** note, as the current README does.
- The scripts are flat, inline, and stepwise by deliberate preference; describe them that way and don't imply a master-script/package structure that isn't there.
- Match the author's own concise, section-numbered style — the README should feel like the scripts.

## Update workflow

When asked to refresh the docs:

1. Find what changed since the README last moved: `git log --name-only --since=<date-of-last-README-commit>`, or diff against the commit that last touched `README.md`.
2. Read only the new/changed scripts — enough to name what they do and their outputs. Their header comments carry the rationale; don't restate it.
3. Update only the affected sections. Keep edits surgical; preserve unchanged wording so the diff stays reviewable.
4. Apply the lean rule to anything you add: name + output + pointer, not a walk-through.
5. Check that any in-README anchor links still resolve, and that the pipeline table's dependency note is still correct.
6. The repo is edited from a clone; leave committing/pushing to the user unless asked.
