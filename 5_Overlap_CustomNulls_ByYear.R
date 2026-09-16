#### CUSTOM NULL MODELS FOR OVERLAP  (flat / stepwise version) ####
# -------------------------------------------------------------------------
# This script tests body-size community assembly against two custom null models,
# each evaluated with a panel of metrics. It is written to be run top-to-bottom,
# but also stepwise: run sections 0-4 once (setup), then run either null section
# (5 or 6) on its own. Set the level and pool at the top and re-run for each
# focal scale you want.
#
# LEVEL / POOL (edit these):
#   LEVEL = "plot" -> focal community = plotID, augmentation grouped by plot
#   LEVEL = "site" -> focal community = siteID, augmentation grouped by site
#   POOL  = "site" or "domain" -> the regional pool the nulls draw from
#   Typical pairings: plot->site, site->domain. To look at the nested-scale
#   trend, also run plot->domain etc. and compare.
#
# TWO NULL MODELS:
#   Section 5  POOL null       Random species assemblage drawn from the regional
#                              pool, holding richness and total N constant. The
#                              workhorse trait-based-assembly test. Direction of
#                              the deviation is the interpretation: overlap lower
#                              than pool = overdispersion (competition-consistent),
#                              higher = clustering (filtering-consistent).
#   Section 6  INDIVIDUAL null Keep the observed species and abundances, but redraw
#                              each species' individuals from that species' pool-
#                              level individuals. Tests whether a different regional
#                              sample of the SAME species changes density / overlap
#                              (the individual-level-data test).
#
# METRICS (each: observed value, null lower/upper CI, SES, direction flag).
# Two families that differ in HOW they treat sparse species:
#
#   -- density family (require >= 2 individuals per species; the density helpers
#      enforce this themselves, so the singleton exclusion touches ONLY these) --
#   overlap_norm    overlap, each density normalized to area 1
#   overlap_unnorm  overlap, density scaled by true abundance
#   overlap_depth   mean peak-scaled species co-occupancy over the occupied range
#                   (community-wide overlap; presence-based, not abundance-weighted)
#
#   -- community-wide / mean family (NO singleton exclusion: every species with a
#      true abundance contributes, singletons included) --
#   niche_range     width of occupied trait space (2.5-97.5% span; unweighted)
#   cwm             community weighted mean (abundance-weighted community trait mean)
#   cw_variance     community-wide variance (abundance-weighted; full distribution)
#   cw_skew         community-wide skewness (abundance-weighted)
#   cw_kurtosis     community-wide excess kurtosis (abundance-weighted; 0 = Gaussian)
#   sdnnd           SD of nearest-neighbour distances between species means
#                   (LOW = even spacing = limiting-similarity signature)
#
# min_logratio was dropped in this revision.
#
# WHY THE SPLIT: augmentation pads sparse species x focal cells to n = 20, and the
# old global `table(sp) >= 2` filter meant a singleton was DROPPED without
# augmentation but INCLUDED (as 20 simulated rows) with it. That asymmetry is the
# main way augmentation moved niche_range and sdnnd. Restricting the >= 2 rule to
# the density metrics removes it: singletons now enter the community-wide and
# mean-based metrics identically with and without augmentation.
# -------------------------------------------------------------------------

library(Ostats)   # community_overlap()
library(dplyr)
library(tidyr)
library(dplyr)

setwd("/home/aly/Beetles/BeetleBodySizeVariation")

# abundance-weighted overlap (weights by true, effort-scaled abundance instead of
# the augmented observation counts). Kept in its own script for readability.
source("./community_overlap_weighted.R")
source("./community_depth.R")

#### 0. SETTINGS ####
# Defaults are for interactive runs. When driven by run_ByYear.sh these are
# overridden by command-line args, in this order:
#   Rscript Overlap_CustomNulls_ByYear.R <LEVEL> <POOL> <YEAR>
LEVEL   <- "site"     # "plot" or "site"
POOL    <- "all"      # "site" / "domain" / "all"
YEAR    <- 2018       # 2018 or 2019
AUGMENT <- TRUE       # TRUE: pad sparse species x focal cells to n = 20 (lognormal)
# FALSE: observed individuals only (no simulated rows)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) LEVEL   <- args[1]
if (length(args) >= 2) POOL    <- args[2]
if (length(args) >= 3) YEAR    <- as.numeric(args[3])
if (length(args) >= 4) AUGMENT <- as.logical(args[4])

NPERM   <- 99
NULLQS  <- c(0.025, 0.975)
SEED    <- 517
MIN_POOL_UNITS <- 2          # drop focal units whose pool holds fewer focal units
# (site->domain: drops single-site domains;
#  plot->site:  drops single-plot sites)

metric_names <- c("overlap_norm", "overlap_unnorm", "overlap_depth",
                  "niche_range", "cwm", "cw_variance", "cw_skew", "cw_kurtosis",
                  "sdnnd")

# translate the LEVEL / POOL choices into column names and the CV scale to use
if (LEVEL == "plot") { FOCAL_COL <- "plotID"; CV_SCALE <- "Plot Level" }
if (LEVEL == "site") { FOCAL_COL <- "siteID"; CV_SCALE <- "Site Level" }
POOL_COL   <- if (POOL == "site") "siteID" else 
  if (POOL == "domain") "domainID" else "all"
AUG_TAG    <- if (AUGMENT) "aug" else "noaug"
OUT_PREFIX <- paste0(LEVEL, "_by_", POOL, "_", AUG_TAG, "_", YEAR)   # e.g. "site_by_all_aug_2018"

#### 1. READ CLEAN DATA ####
all_elytra <- read.csv("./Data/BodysizeCombinedClean.csv")
all_elytra <- subset(all_elytra, yearCollected == YEAR)
message("Year: ", YEAR, "  |  specimens: ", nrow(all_elytra))
all_elytra$all<-"all"

# Drop specimens with no species label (NA or blank ""). A blank name is not a
# real species: it survives %in% filters and silently becomes a pseudo-species in
# every metric, and it breaks name-indexing in the swap null (x[""] returns NA in
# R). Must happen BEFORE augmentation, or the blank cell gets augmented to n = 20.
blank_sp <- is.na(all_elytra$scientificName_Species) | all_elytra$scientificName_Species == ""
if (any(blank_sp)) message("Dropping ", sum(blank_sp), " specimen(s) with missing species label")
all_elytra <- all_elytra[!blank_sp, ]

all_elytra$log_dist_cm <- log10(all_elytra$cm_elytra_max_length)

# crosswalks from the observed rows, used to re-attach pool ids to augmented rows
plot_to_site   <- unique(all_elytra[, c("plotID", "siteID")])
site_to_domain <- unique(all_elytra[, c("siteID", "domainID")])

# average CV^2 (%) for this scale, from the summary file (not hardcoded)
cv_results    <- read.csv("./Outputs/CVpctSummary.csv")
typical_cvpct <- cv_results$mean_cvpct[cv_results$scale == CV_SCALE]
cv2           <- typical_cvpct / 100
message("Level: ", LEVEL, "  |  Pool: ", POOL, "  |  CV scale: ", CV_SCALE,
        " (", round(typical_cvpct, 4), ")")


#### 2. (OPTIONALLY) AUGMENT SPARSE SPECIES x FOCAL CELLS TO n >= 20 ####
# Start from observed individuals; if AUGMENT, pad each sparse species x focal
# cell up to n = 20 with lognormal draws (same augmentation as Overlap_*.R). With
# AUGMENT = FALSE only observed individuals enter. Augmentation feeds ALL metrics;
# what differs by metric is the eligibility filtering applied in community_metrics().
aug <- all_elytra[, c("scientificName_Species", FOCAL_COL, "cm_elytra_max_length")]

if (AUGMENT) {
  low_n <- all_elytra %>%
    group_by(scientificName_Species, across(all_of(FOCAL_COL))) %>%
    summarise(n_obs = n(),
              mean_dist = mean(cm_elytra_max_length, na.rm = TRUE),
              .groups = "drop") %>%
    filter(n_obs < 20)
  
  set.seed(42)
  sim_low_n <- low_n %>%
    rowwise() %>%
    mutate(n_to_add = 20 - n_obs,
           sdlog    = sqrt(log(1 + cv2)),
           meanlog  = log(mean_dist) - (sdlog^2 / 2),
           sim_vals = list(rlnorm(n = n_to_add, meanlog = meanlog, sdlog = sdlog))) %>%
    unnest(cols = sim_vals) %>%
    dplyr::rename(cm_elytra_max_length = sim_vals) %>%
    select(all_of(c("scientificName_Species", FOCAL_COL)), cm_elytra_max_length) %>%
    ungroup()
  
  aug <- rbind(aug, as.data.frame(sim_low_n))
  message("Augmentation ON: added ", nrow(sim_low_n), " simulated individuals")
} else {
  message("Augmentation OFF: observed individuals only (", nrow(aug), " individuals)")
}

aug$log_dist_cm <- log10(aug$cm_elytra_max_length)

#### 3. ATTACH REGIONAL POOL IDS ####
# augmented rows carry only the focal id, so re-attach the higher-level ids
if (LEVEL == "plot") {
  aug <- merge(aug, plot_to_site,   by = "plotID", all.x = TRUE)   # + siteID
  aug <- merge(aug, site_to_domain, by = "siteID", all.x = TRUE)   # + domainID
} else {
  aug <- merge(aug, site_to_domain, by = "siteID", all.x = TRUE)   # + domainID
}
if (POOL == "all") { 
  aug$all<-"all"
}
aug$FOCAL <- aug[[FOCAL_COL]]     # generic working columns used below
aug$POOL  <- aug[[POOL_COL]]


#### 4. BUILD FOCAL LIST, DROP UNITS WITH NO REGIONAL POOL, GET LATITUDE ####
focal_pool <- unique(aug[, c("FOCAL", "POOL")])
pool_size  <- table(focal_pool$POOL)
focal_pool$n_in_pool <- as.integer(pool_size[focal_pool$POOL])

excluded    <- sort(focal_pool$FOCAL[focal_pool$n_in_pool < MIN_POOL_UNITS])
if (length(excluded))
  message("Excluding ", length(excluded), " focal unit(s) with pool < ",
          MIN_POOL_UNITS, ": ", paste(excluded, collapse = ", "))
focal_units <- sort(setdiff(unique(aug$FOCAL), excluded))

# mean latitude per focal unit (observed rows) for the latitudinal framing
all_elytra$FOCAL <- all_elytra[[FOCAL_COL]]
lat <- aggregate(latitude ~ FOCAL, data = all_elytra, FUN = mean)


#### 4b. LOAD TRUE (EFFORT-SCALED) ABUNDANCES ####
# One row per focal-unit x species: columns <FOCAL_COL>, scientificName_Species, abund.
# Keyed as "focal|species" so a community's weight vector is a single lookup.
abund_tab <- read.csv(sprintf("./Data/%s_abund_%d.csv",
                              if (LEVEL == "plot") "plot" else "site", YEAR))
abund_lookup <- setNames(abund_tab$abund,
                         paste(abund_tab[[FOCAL_COL]], abund_tab$scientificName_Species, sep = "|"))

# helper: observed true-abundance vector (species -> abund) for one focal unit,
# keeping only species that actually have a true abundance
abund_for <- function(f, species) {
  species <- unique(as.character(species))
  a <- abund_lookup[paste(f, species, sep = "|")]
  a <- setNames(as.numeric(a), species)
  a[is.finite(a)]
}

# coverage report: observed species with no true abundance (dropped from every metric)
obs_keys     <- unique(paste(aug$FOCAL, aug$scientificName_Species, sep = "|"))
missing_keys <- obs_keys[is.na(abund_lookup[obs_keys])]
if (length(missing_keys))
  message("NOTE: ", length(missing_keys), " observed (focal|species) combos have no true abundance ",
          "and drop from the metrics. e.g. ", paste(missing_keys, collapse = "; "))

#### 4c. SPECIES THAT ENTER THE OVERLAP, PER FOCAL UNIT ####
# Count of species that actually contribute to the DENSITY metrics at each focal
# unit: the eligibility rule the density helpers use (finite trait, a true
# abundance, >= 2 individuals) = the number of density curves community_overlap
# builds. A property of the OBSERVED community, so it does not depend on the
# null -- computed once here and merged into each null's output, like `lat`.
n_sp_tab <- data.frame(FOCAL = focal_units, n_overlap_sp = NA_integer_,
                       stringsAsFactors = FALSE)
for (r in seq_along(focal_units)) {
  f        <- focal_units[r]
  in_focal <- aug$FOCAL == f
  sp_f     <- aug$scientificName_Species[in_focal]
  tr_f     <- aug$log_dist_cm[in_focal]
  abund_f  <- abund_for(f, sp_f)
  
  ok  <- is.finite(tr_f) & !is.na(sp_f) & sp_f != "" & sp_f %in% names(abund_f)
  spp <- sp_f[ok]
  n_sp_tab$n_overlap_sp[r] <- length(names(which(table(spp) >= 2)))
}
message("Computable focal units this year (>= 2 overlap species): ",
        sum(n_sp_tab$n_overlap_sp >= 2, na.rm = TRUE), " of ", nrow(n_sp_tab))

#### HELPER: abundance-weighted community moments ####
# Moments of the community trait distribution, weighting each species by its TRUE
# (effort-scaled) relative abundance rather than by its (augmented) observation
# count. Each species' TOTAL weight is pinned to its relative abundance and its
# individuals are weighted equally within it, so:
#   - the mean reduces to the classic community weighted mean, sum_s p_s * mean_s;
#   - adding augmented rows to a species does not change that species' influence,
#     only refines its internal shape -> the between-species structure is
#     augmentation-invariant (the within-species spread still reflects the cv2
#     lognormal once a cell is padded).
# Returns mean, variance, skewness, excess kurtosis. NA moments where undefined.
#
# ALTERNATIVES if you want a different definition (localized change, this helper
# only): for UNWEIGHTED pooled-individual moments, set w <- rep(1/length(traits),
# length(traits)); for BETWEEN-SPECIES-MEANS moments (fully augmentation-robust,
# ignores intraspecific spread), collapse to species means first and weight those
# by p_s. For RAW (non-excess) kurtosis, drop the "- 3".
cw_moments <- function(traits, sp, abund) {
  out <- c(mean = NA_real_, var = NA_real_, skew = NA_real_, kurt = NA_real_)
  traits <- as.numeric(traits); sp <- as.character(sp)
  ok <- is.finite(traits) & !is.na(sp) & sp %in% names(abund)
  traits <- traits[ok]; sp <- sp[ok]
  if (length(traits) < 1 || length(unique(sp)) < 2) return(out)

  spp <- unique(sp)
  p   <- abund[spp]
  if (!all(is.finite(p)) || sum(p) <= 0) return(out)
  p   <- p / sum(p)                                    # relative abundance, sums to 1

  n_s <- as.numeric(table(sp)[sp])                     # individuals in each obs's species
  w   <- as.numeric(p[sp]) / n_s                       # individual weight: species total = p_s
  w   <- w / sum(w)                                    # guard: renormalize to 1

  mu  <- sum(w * traits)
  m2  <- sum(w * (traits - mu)^2)
  m3  <- sum(w * (traits - mu)^3)
  m4  <- sum(w * (traits - mu)^4)

  out["mean"] <- mu
  out["var"]  <- m2
  if (m2 > 0) {
    out["skew"] <- m3 / m2^1.5
    out["kurt"] <- m4 / m2^2 - 3                       # excess kurtosis (0 = Gaussian)
  }
  out
}

#### HELPER: the metric panel for one community ####
# The only community-metric function in the script. It is called once for the
# observed community and once for every null draw, so the metric definitions live
# in exactly one place. Returns a named vector in metric_names order; NA where a
# metric is undefined.
#
# `abund` is a named vector (species -> true, effort-scaled abundance). Species
# without a true abundance are dropped from ALL metrics so every focal unit
# describes one consistent community.
#
# FILTERING, by design:
#   base:            finite trait + real species label + has a true abundance.
#   density metrics: overlap_norm / overlap_unnorm / overlap_depth additionally
#                    need >= 2 individuals per species. That rule is applied
#                    INSIDE the helpers (community_overlap_weighted, community_depth),
#                    so it touches only these three metrics -- singletons drop out
#                    of the overlaps but nowhere else.
#   community-wide / mean metrics: niche_range, cwm, cw_variance, cw_skew,
#                    cw_kurtosis, sdnnd use every base-eligible species, singletons
#                    included. No >= 2 filter here, so augmentation no longer flips
#                    singletons in and out of these metrics.
community_metrics <- function(traits, sp, abund) {
  traits <- as.numeric(traits); sp <- as.character(sp)
  ok <- is.finite(traits) & !is.na(sp) & sp != "" & sp %in% names(abund)   # base: trait + true abundance
  traits <- traits[ok]; sp <- sp[ok]

  out <- c(overlap_norm = NA, overlap_unnorm = NA, overlap_depth = NA,
           niche_range = NA, cwm = NA, cw_variance = NA, cw_skew = NA,
           cw_kurtosis = NA, sdnnd = NA)
  if (length(unique(sp)) < 2) return(out)

  ## ---- density family: helpers self-filter to species with >= 2 individuals ----
  # overlaps weighted by TRUE abundance (not the augmented observation counts)
  out["overlap_norm"]   <- community_overlap_weighted(traits, sp, abund, normal = TRUE,  output = "mean")
  out["overlap_unnorm"] <- community_overlap_weighted(traits, sp, abund, normal = FALSE, output = "mean")
  # community-wide overlap depth (presence-based; not abundance-weighted by design)
  out["overlap_depth"]  <- community_depth(traits, sp)

  ## ---- community-wide distribution family: NO >= 2 filter (singletons included) ----
  # width of occupied trait space (robust 2.5-97.5% span; unweighted, unchanged)
  out["niche_range"] <- diff(quantile(traits, c(0.025, 0.975)))

  # abundance-weighted community moments (mean / variance / skew / excess kurtosis)
  moms <- cw_moments(traits, sp, abund)
  out["cwm"]         <- moms["mean"]
  out["cw_variance"] <- moms["var"]
  out["cw_skew"]     <- moms["skew"]
  out["cw_kurtosis"] <- moms["kurt"]

  ## ---- spacing family (mean-based): NO >= 2 filter, means from every species ----
  means <- sort(tapply(traits, sp, mean))
  if (length(means) >= 2) {
    gaps <- diff(means)                              # adjacent gaps = log10 size ratios
    nn   <- pmin(c(gaps, Inf), c(Inf, gaps))         # nearest-neighbour distance per species
    out["sdnnd"] <- sd(nn)
  }
  out
}


#### 5. NULL MODEL 1: RANDOM ASSEMBLAGE FROM THE REGIONAL POOL ####
set.seed(SEED)

pool_results <- data.frame(FOCAL = focal_units, stringsAsFactors = FALSE)
pool_results$POOL <- focal_pool$POOL[match(pool_results$FOCAL, focal_pool$FOCAL)]
for (m in metric_names) for (s in c("_obs","_lower","_upper","_ses")) pool_results[[paste0(m, s)]] <- NA_real_
for (m in metric_names) pool_results[[paste0(m, "_dir")]] <- NA_character_

for (r in seq_along(focal_units)) {
  f  <- focal_units[r]
  in_focal <- aug$FOCAL == f
  in_pool  <- aug$POOL  == pool_results$POOL[r]
  
  traits_obs  <- aug$log_dist_cm[in_focal];  sp_obs  <- aug$scientificName_Species[in_focal]
  traits_pool <- aug$log_dist_cm[in_pool];   sp_pool <- aug$scientificName_Species[in_pool]
  
  # observed metrics (weighted by this focal unit's true abundances)
  abund_f <- abund_for(f, sp_obs)
  obs <- community_metrics(traits_obs, sp_obs, abund_f)
  
  # null draws. Restrict the observed community to species that have a true
  # abundance, then PRESERVE two paired vectors and re-label them to the drawn
  # species: n_ind (individuals sampled -> density shape) and w_obs (true
  # abundance -> overlap weight). The augmented sample gives the shape; the
  # NEON div abundance gives the weight.
  null_mat <- matrix(NA, nrow = NPERM, ncol = length(metric_names),
                     dimnames = list(NULL, metric_names))
  obs_tab      <- table(sp_obs[sp_obs %in% names(abund_f)])
  obs_species  <- names(obs_tab)
  n_ind        <- as.numeric(obs_tab)
  w_obs        <- as.numeric(abund_f[obs_species])
  pool_species <- unique(sp_pool)
  for (i in 1:NPERM) {
    # ---- random assemblage from the pool: same richness, preserved (count, abundance) vectors ----
    drawn <- sample(pool_species, length(obs_species))
    traits_null <- numeric(0); sp_null <- character(0)
    for (k in seq_along(drawn)) {
      pool_k      <- traits_pool[sp_pool == drawn[k]]
      traits_null <- c(traits_null, sample(pool_k, n_ind[k], replace = TRUE))
      sp_null     <- c(sp_null, rep(drawn[k], n_ind[k]))
    }
    abund_null <- setNames(w_obs, drawn)          # observed true-abundance vector, re-labelled
    null_mat[i, ] <- community_metrics(traits_null, sp_null, abund_null)
  }
  
  # summarise observed vs null per metric; if a metric is invariant under this
  # null (sd ~ 0) report CI = obs, ses = NA
  for (m in metric_names) {
    o  <- obs[m]
    nd <- null_mat[, m]; nd <- nd[is.finite(nd)]
    pool_results[r, paste0(m, "_obs")] <- o
    if (length(nd) >= 2 && is.finite(o)) {
      if (sd(nd) > 1e-9) {
        lo <- as.numeric(quantile(nd, NULLQS[1])); hi <- as.numeric(quantile(nd, NULLQS[2]))
        pool_results[r, paste0(m, "_lower")] <- lo
        pool_results[r, paste0(m, "_upper")] <- hi
        pool_results[r, paste0(m, "_ses")]   <- (o - mean(nd)) / sd(nd)
        pool_results[r, paste0(m, "_dir")]   <- if (o < lo) "lower" else if (o > hi) "higher" else "neutral"
      } else {
        pool_results[r, paste0(m, "_lower")] <- o
        pool_results[r, paste0(m, "_upper")] <- o
        pool_results[r, paste0(m, "_ses")]   <- NA
        pool_results[r, paste0(m, "_dir")]   <- "neutral"
      }
    }
  }
}

names(pool_results)[1:2] <- c(FOCAL_COL, POOL_COL)
pool_results <- merge(pool_results, lat, by.x = FOCAL_COL, by.y = "FOCAL", all.x = TRUE)
pool_results <- merge(pool_results, n_sp_tab, by.x = FOCAL_COL, by.y = "FOCAL", all.x = TRUE)   # <- add
write.csv(pool_results, paste0("./Outputs/", OUT_PREFIX, "_PoolNull.csv"), row.names = FALSE)
message("wrote ", OUT_PREFIX, "_PoolNull.csv  (", nrow(pool_results), " focal units)")


#### 6. NULL MODEL 2: REGIONAL RESAMPLE OF INDIVIDUALS WITHIN SPECIES ####
set.seed(SEED)

indiv_results <- data.frame(FOCAL = focal_units, stringsAsFactors = FALSE)
indiv_results$POOL <- focal_pool$POOL[match(indiv_results$FOCAL, focal_pool$FOCAL)]
for (m in metric_names) for (s in c("_obs","_lower","_upper","_ses")) indiv_results[[paste0(m, s)]] <- NA_real_
for (m in metric_names) indiv_results[[paste0(m, "_dir")]] <- NA_character_

for (r in seq_along(focal_units)) {
  f  <- focal_units[r]
  in_focal <- aug$FOCAL == f
  in_pool  <- aug$POOL  == indiv_results$POOL[r]
  
  traits_obs  <- aug$log_dist_cm[in_focal];  sp_obs  <- aug$scientificName_Species[in_focal]
  traits_pool <- aug$log_dist_cm[in_pool];   sp_pool <- aug$scientificName_Species[in_pool]
  
  # species (and thus true abundances) are preserved, so use the observed weights
  abund_f <- abund_for(f, sp_obs)
  obs <- community_metrics(traits_obs, sp_obs, abund_f)
  
  null_mat <- matrix(NA, nrow = NPERM, ncol = length(metric_names),
                     dimnames = list(NULL, metric_names))
  obs_species <- unique(sp_obs)
  for (i in 1:NPERM) {
    # ---- keep species & abundances; redraw each species' individuals from the pool ----
    traits_null <- numeric(0); sp_null <- character(0)
    for (s in obs_species) {
      n_s    <- sum(sp_obs == s)
      pool_s <- traits_pool[sp_pool == s]
      draw_s <- if (length(pool_s) <= n_s) pool_s else sample(pool_s, n_s)   # no replacement
      traits_null <- c(traits_null, draw_s)
      sp_null     <- c(sp_null, rep(s, n_s))
    }
    null_mat[i, ] <- community_metrics(traits_null, sp_null, abund_f)
  }
  
  # summarise observed vs null per metric; if a metric is invariant under this
  # null (sd ~ 0, e.g. spacing metrics under swap_means) report CI = obs, ses = NA
  for (m in metric_names) {
    o  <- obs[m]
    nd <- null_mat[, m]; nd <- nd[is.finite(nd)]
    indiv_results[r, paste0(m, "_obs")] <- o
    if (length(nd) >= 2 && is.finite(o)) {
      if (sd(nd) > 1e-9) {
        lo <- as.numeric(quantile(nd, NULLQS[1])); hi <- as.numeric(quantile(nd, NULLQS[2]))
        indiv_results[r, paste0(m, "_lower")] <- lo
        indiv_results[r, paste0(m, "_upper")] <- hi
        indiv_results[r, paste0(m, "_ses")]   <- (o - mean(nd)) / sd(nd)
        indiv_results[r, paste0(m, "_dir")]   <- if (o < lo) "lower" else if (o > hi) "higher" else "neutral"
      } else {
        indiv_results[r, paste0(m, "_lower")] <- o
        indiv_results[r, paste0(m, "_upper")] <- o
        indiv_results[r, paste0(m, "_ses")]   <- NA
        indiv_results[r, paste0(m, "_dir")]   <- "neutral"
      }
    }
  }
}

names(indiv_results)[1:2] <- c(FOCAL_COL, POOL_COL)
indiv_results <- merge(indiv_results, lat, by.x = FOCAL_COL, by.y = "FOCAL", all.x = TRUE)
indiv_results <- merge(indiv_results, n_sp_tab, by.x = FOCAL_COL, by.y = "FOCAL", all.x = TRUE)   # <- add
write.csv(indiv_results, paste0("./Outputs/", OUT_PREFIX, "_IndividualNull.csv"), row.names = FALSE)
message("wrote ", OUT_PREFIX, "_IndividualNull.csv  (", nrow(indiv_results), " focal units)")


# #### 7. NULL MODEL 3: SWAP MEANS (within-community, out-of-the-box Ostats null) ####
# # The mean-swap null from Ostats::Ostats(swap_means = TRUE): within each focal
# # community, keep every species' abundance and its within-species deviations, but
# # relocate each species onto a randomly permuted community mean. This needs NO
# # regional pool (it is a within-community null), so it runs on every focal unit
# # and its result does not depend on POOL -- the output is named by LEVEL only.
# #
# # NOTE: because the SET of species means is only permuted (never changed), the
# # spacing metric is invariant by construction -- sdnnd has null CI = observed,
# # ses = NA, dir = "neutral". That is expected, not a bug: this null is only
# # informative for the overlap metrics (and weakly niche_range / cw moments). We
# # run it on the same focal_units as sections 5-6 so the three files line up; to
# # also cover the pool-excluded units, loop over sort(unique(aug$FOCAL)) instead.
# set.seed(SEED)
# 
# swap_results <- data.frame(FOCAL = focal_units, stringsAsFactors = FALSE)
# swap_results$POOL <- focal_pool$POOL[match(swap_results$FOCAL, focal_pool$FOCAL)]
# for (m in metric_names) for (s in c("_obs","_lower","_upper","_ses")) swap_results[[paste0(m, s)]] <- NA_real_
# for (m in metric_names) swap_results[[paste0(m, "_dir")]] <- NA_character_
# 
# for (r in seq_along(focal_units)) {
#   f <- focal_units[r]
#   in_focal <- aug$FOCAL == f
#   traits_obs <- aug$log_dist_cm[in_focal];  sp_obs <- aug$scientificName_Species[in_focal]
# 
#   # species (and thus true abundances) are preserved, so use the observed weights
#   abund_f <- abund_for(f, sp_obs)
#   obs <- community_metrics(traits_obs, sp_obs, abund_f)
# 
#   null_mat <- matrix(NA, nrow = NPERM, ncol = length(metric_names),
#                      dimnames = list(NULL, metric_names))
# 
#   # swap operates on exactly the community the overlaps use: finite traits, species
#   # with >= 2 individuals AND a true abundance. This keeps non-finite / singleton /
#   # unweightable means from leaking a bad value onto a real species when permuted.
#   keep <- is.finite(traits_obs) & !is.na(sp_obs) & sp_obs %in% names(abund_f)
#   tr   <- traits_obs[keep]; spp <- sp_obs[keep]
#   elig <- names(which(table(spp) >= 2))
#   tr   <- tr[spp %in% elig]; spp <- spp[spp %in% elig]
# 
#   if (length(unique(spp)) >= 2) {
#     sp_f  <- factor(spp)                             # index species by position, not by name
#     means <- as.numeric(tapply(tr, sp_f, mean))      # one mean per level, in level order
#     codes <- as.integer(sp_f)                        # each individual's species code
#     devs  <- tr - means[codes]                       # each individual's deviation from its own mean
#     for (i in 1:NPERM) {
#       # ---- permute the community means across species; keep identity/abundance/shape ----
#       means_swapped <- sample(means)
#       traits_null   <- devs + means_swapped[codes]
#       null_mat[i, ] <- community_metrics(traits_null, spp, abund_f)
#     }
#   }
# 
#   # summarise observed vs null per metric; sdnnd is invariant under swap_means
#   # (sd ~ 0) so it reports CI = obs, ses = NA (expected, not a bug)
#   for (m in metric_names) {
#     o  <- obs[m]
#     nd <- null_mat[, m]; nd <- nd[is.finite(nd)]
#     swap_results[r, paste0(m, "_obs")] <- o
#     if (length(nd) >= 2 && is.finite(o)) {
#       if (sd(nd) > 1e-9) {
#         lo <- as.numeric(quantile(nd, NULLQS[1])); hi <- as.numeric(quantile(nd, NULLQS[2]))
#         swap_results[r, paste0(m, "_lower")] <- lo
#         swap_results[r, paste0(m, "_upper")] <- hi
#         swap_results[r, paste0(m, "_ses")]   <- (o - mean(nd)) / sd(nd)
#         swap_results[r, paste0(m, "_dir")]   <- if (o < lo) "lower" else if (o > hi) "higher" else "neutral"
#       } else {
#         swap_results[r, paste0(m, "_lower")] <- o
#         swap_results[r, paste0(m, "_upper")] <- o
#         swap_results[r, paste0(m, "_ses")]   <- NA
#         swap_results[r, paste0(m, "_dir")]   <- "neutral"
#       }
#     }
#   }
# }
# 
# names(swap_results)[1:2] <- c(FOCAL_COL, POOL_COL)
# swap_results <- merge(swap_results, lat, by.x = FOCAL_COL, by.y = "FOCAL", all.x = TRUE)
# swap_results <- merge(swap_results, n_sp_tab, by.x = FOCAL_COL, by.y = "FOCAL", all.x = TRUE)   # <- add
# write.csv(swap_results, paste0("./Outputs/", LEVEL, "_SwapMeansNull.csv"), row.names = FALSE)
# message("wrote ", LEVEL, "_SwapMeansNull.csv  (", nrow(swap_results), " focal units)")
