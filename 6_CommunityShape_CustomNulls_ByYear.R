#### COMMUNITY-SHAPE CUSTOM NULL MODELS  (no augmentation; flat / stepwise) ####
# -------------------------------------------------------------------------
# This script tests body-size community structure against three custom null
# models, each evaluated on a panel of community-shape metrics. Unlike the
# overlap version, it uses NO augmentation: every metric runs on observed
# individuals, and no species is dropped for having fewer than two observations.
# The pairwise overlap metrics (overlap_norm / overlap_unnorm) are gone; depth
# and the distribution/mean metrics carry the analysis.
#
# It is written to run top-to-bottom, but also stepwise: run sections 0-4 once
# (setup), then run any null section (5, 6, or 7) on its own. Set the level and
# pool at the top and re-run for each focal scale you want.
#
# LEVEL / POOL (edit these):
#   LEVEL = "plot" -> focal community = plotID
#   LEVEL = "site" -> focal community = siteID
#   POOL  = "site" / "domain" / "all" -> the regional pool the nulls draw from
#   Typical pairings: plot->site, site->domain.
#
# THREE NULL MODELS:
#   Section 5  POOL null       Random species assemblage drawn from the regional
#                              pool, holding richness and total N constant. The
#                              workhorse assembly test. Direction of the deviation
#                              is the interpretation: depth higher than pool =
#                              clustering (filtering-consistent), lower =
#                              overdispersion (competition-consistent).
#   Section 6  INDIVIDUAL null Keep the observed species and abundances, but redraw
#                              each species' individuals from that species' pool-
#                              level individuals (the individual-level-data test).
#   Section 7  SWAP-MEANS null Within each focal community, keep every species'
#                              within-species deviations but relocate each species
#                              onto a permuted community mean. Needs no pool.
#
# METRICS (each: observed value, null lower/upper CI, SES, direction flag). All
# run on observed individuals, every species with a true abundance included
# (singletons kept):
#   overlap_depth   mean peak-scaled species co-occupancy over the occupied range
#                   (fixed KDE bandwidth; presence-based, not abundance-weighted)
#   niche_range     width of occupied trait space (2.5-97.5% span; unweighted)
#   cwm             community weighted mean (abundance-weighted trait mean)
#   cw_variance     community-wide variance (abundance-weighted; full distribution)
#   cw_skew         community-wide skewness (abundance-weighted)
#   cw_kurtosis     community-wide excess kurtosis (abundance-weighted; 0 = Gaussian)
#   sdnnd           SD of nearest-neighbour distances between species means
#                   (LOW = even spacing = limiting-similarity signature)
#
# Abundance filtering is kept as the community definition: a species enters the
# metrics only if it has a true effort-scaled abundance, which harmonizes the
# species list across data sources. It does NOT drop singletons.
# -------------------------------------------------------------------------

setwd("/home/aly/Beetles/BeetleBodySizeVariation")

source("./community_depth.R")

#### 0. SETTINGS ####
# Defaults are for interactive runs. When driven by a shell wrapper these are
# overridden by command-line args, in this order:
#   Rscript CommunityShape_CustomNulls_ByYear.R <LEVEL> <POOL> <YEAR>
LEVEL   <- "site"     # "plot" or "site"
POOL    <- "all"      # "site" / "domain" / "all"
YEAR    <- 2018       # 2018 or 2019

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) LEVEL <- args[1]
if (length(args) >= 2) POOL  <- args[2]
if (length(args) >= 3) YEAR  <- as.numeric(args[3])

NPERM   <- 99
NULLQS  <- c(0.025, 0.975)
SEED    <- 517
MIN_POOL_UNITS <- 2          # drop focal units whose pool holds fewer focal units
# (site->domain: drops single-site domains;
#  plot->site:  drops single-plot sites)

# Fixed KDE bandwidth for community_depth(), on the log10 elytra-length axis.
# Anchored to the median bw.nrd0 across species x focal cells with n >= 20
# (median 0.00953875 over that set), so sparse species and singletons receive
# the same smoothing as well-sampled ones. Rerun depth at 0.5x / 2x this value
# as a sensitivity check on the SES sign.
DEPTH_BW <- 0.009539

metric_names <- c("overlap_depth", "niche_range", "cwm", "cw_variance",
                  "cw_skew", "cw_kurtosis", "sdnnd")

# translate the LEVEL / POOL choices into column names
if (LEVEL == "plot") FOCAL_COL <- "plotID"
if (LEVEL == "site") FOCAL_COL <- "siteID"
POOL_COL   <- if (POOL == "site") "siteID" else
  if (POOL == "domain") "domainID" else "all"
OUT_PREFIX <- paste0(LEVEL, "_by_", POOL, "_", YEAR)   # e.g. "site_by_all_2018"

#### 1. READ CLEAN DATA ####
all_elytra <- read.csv("./Data/BodysizeCombinedClean.csv")
all_elytra <- subset(all_elytra, yearCollected == YEAR)
message("Year: ", YEAR, "  |  specimens: ", nrow(all_elytra))
all_elytra$all <- "all"

# Drop specimens with no species label (NA or blank ""). A blank name is not a
# real species: it survives %in% filters and silently becomes a pseudo-species in
# every metric, and it breaks name-indexing in the swap null (x[""] returns NA).
blank_sp <- is.na(all_elytra$scientificName_Species) | all_elytra$scientificName_Species == ""
if (any(blank_sp)) message("Dropping ", sum(blank_sp), " specimen(s) with missing species label")
all_elytra <- all_elytra[!blank_sp, ]

all_elytra$log_dist_cm <- log10(all_elytra$cm_elytra_max_length)

# crosswalks used to re-attach pool ids
plot_to_site   <- unique(all_elytra[, c("plotID", "siteID")])
site_to_domain <- unique(all_elytra[, c("siteID", "domainID")])
message("Level: ", LEVEL, "  |  Pool: ", POOL, "  |  depth bw: ", DEPTH_BW)


#### 2. INDIVIDUAL-LEVEL WORKING DATA (NO AUGMENTATION) ####
# Observed individuals only. No simulated rows anywhere in this script.
dat <- all_elytra[, c("scientificName_Species", FOCAL_COL, "cm_elytra_max_length")]
dat$log_dist_cm <- log10(dat$cm_elytra_max_length)
message("No augmentation: ", nrow(dat), " observed individuals")


#### 3. ATTACH REGIONAL POOL IDS ####
if (LEVEL == "plot") {
  dat <- merge(dat, plot_to_site,   by = "plotID", all.x = TRUE)   # + siteID
  dat <- merge(dat, site_to_domain, by = "siteID", all.x = TRUE)   # + domainID
} else {
  dat <- merge(dat, site_to_domain, by = "siteID", all.x = TRUE)   # + domainID
}
if (POOL == "all") dat$all <- "all"
dat$FOCAL <- dat[[FOCAL_COL]]     # generic working columns used below
dat$POOL  <- dat[[POOL_COL]]


#### 4. BUILD FOCAL LIST, DROP UNITS WITH NO REGIONAL POOL, GET LATITUDE ####
focal_pool <- unique(dat[, c("FOCAL", "POOL")])
pool_size  <- table(focal_pool$POOL)
focal_pool$n_in_pool <- as.integer(pool_size[focal_pool$POOL])

excluded    <- sort(focal_pool$FOCAL[focal_pool$n_in_pool < MIN_POOL_UNITS])
if (length(excluded))
  message("Excluding ", length(excluded), " focal unit(s) with pool < ",
          MIN_POOL_UNITS, ": ", paste(excluded, collapse = ", "))
focal_units <- sort(setdiff(unique(dat$FOCAL), excluded))

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
obs_keys     <- unique(paste(dat$FOCAL, dat$scientificName_Species, sep = "|"))
missing_keys <- obs_keys[is.na(abund_lookup[obs_keys])]
if (length(missing_keys))
  message("NOTE: ", length(missing_keys), " observed (focal|species) combos have no true abundance ",
          "and drop from the metrics. e.g. ", paste(missing_keys, collapse = "; "))

#### 4c. SPECIES THAT ENTER THE METRICS, PER FOCAL UNIT ####
# Count of species contributing to the metrics at each focal unit: finite trait
# and a true abundance. NO >= 2 filter -- singletons count. A property of the
# OBSERVED community, so it does not depend on the null; computed once and merged
# into each null's output, like `lat`.
n_sp_tab <- data.frame(FOCAL = focal_units, n_comm_sp = NA_integer_,
                       stringsAsFactors = FALSE)
for (r in seq_along(focal_units)) {
  f        <- focal_units[r]
  in_focal <- dat$FOCAL == f
  sp_f     <- dat$scientificName_Species[in_focal]
  tr_f     <- dat$log_dist_cm[in_focal]
  abund_f  <- abund_for(f, sp_f)

  ok  <- is.finite(tr_f) & !is.na(sp_f) & sp_f != "" & sp_f %in% names(abund_f)
  n_sp_tab$n_comm_sp[r] <- length(unique(sp_f[ok]))
}
message("Computable focal units this year (>= 2 community species): ",
        sum(n_sp_tab$n_comm_sp >= 2, na.rm = TRUE), " of ", nrow(n_sp_tab))

#### HELPER: abundance-weighted community moments ####
# Moments of the community trait distribution, weighting each species by its TRUE
# (effort-scaled) relative abundance rather than by its observation count. Each
# species' TOTAL weight is pinned to its relative abundance and its individuals
# are weighted equally within it, so the mean reduces to the classic community
# weighted mean, sum_s p_s * mean_s. With no augmentation the higher moments now
# reflect the REAL observed within-species spread (noisier where species are
# thin), which is the intended behaviour here.
# Returns mean, variance, skewness, excess kurtosis. NA moments where undefined.
#
# ALTERNATIVES (localized change, this helper only): UNWEIGHTED pooled-individual
# moments -> w <- rep(1/length(traits), length(traits)); BETWEEN-SPECIES-MEANS
# moments -> collapse to species means first, weight those by p_s; RAW (non-excess)
# kurtosis -> drop the "- 3".
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
# Called once for the observed community and once for every null draw, so the
# metric definitions live in exactly one place. Returns a named vector in
# metric_names order; NA where a metric is undefined.
#
# `abund` is a named vector (species -> true, effort-scaled abundance). Species
# without a true abundance are dropped from ALL metrics so every focal unit
# describes one consistent community. NO species is dropped for having a single
# observation: depth uses the fixed DEPTH_BW so singletons keep a kernel, and the
# distribution / mean metrics include them directly.
community_metrics <- function(traits, sp, abund) {
  traits <- as.numeric(traits); sp <- as.character(sp)
  ok <- is.finite(traits) & !is.na(sp) & sp != "" & sp %in% names(abund)   # community = species with a true abundance
  traits <- traits[ok]; sp <- sp[ok]

  out <- c(overlap_depth = NA, niche_range = NA, cwm = NA, cw_variance = NA,
           cw_skew = NA, cw_kurtosis = NA, sdnnd = NA)
  if (length(unique(sp)) < 2) return(out)

  # community-wide overlap depth (fixed bandwidth; presence-based, all species)
  out["overlap_depth"] <- community_depth(traits, sp, density_args = list(bw = DEPTH_BW))

  # width of occupied trait space (robust 2.5-97.5% span; unweighted)
  out["niche_range"] <- diff(quantile(traits, c(0.025, 0.975)))

  # abundance-weighted community moments (mean / variance / skew / excess kurtosis)
  moms <- cw_moments(traits, sp, abund)
  out["cwm"]         <- moms["mean"]
  out["cw_variance"] <- moms["var"]
  out["cw_skew"]     <- moms["skew"]
  out["cw_kurtosis"] <- moms["kurt"]

  # spacing of species means on the log10 axis (means from every species)
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
  in_focal <- dat$FOCAL == f
  in_pool  <- dat$POOL  == pool_results$POOL[r]

  traits_obs  <- dat$log_dist_cm[in_focal];  sp_obs  <- dat$scientificName_Species[in_focal]
  traits_pool <- dat$log_dist_cm[in_pool];   sp_pool <- dat$scientificName_Species[in_pool]

  # observed metrics (weighted by this focal unit's true abundances)
  abund_f <- abund_for(f, sp_obs)
  obs <- community_metrics(traits_obs, sp_obs, abund_f)

  # null draws: random assemblage from the pool, same richness, preserved
  # (count, abundance) vectors. Restrict to species with a true abundance, keep
  # n_ind (individuals sampled -> shape) and w_obs (true abundance -> weight).
  null_mat <- matrix(NA, nrow = NPERM, ncol = length(metric_names),
                     dimnames = list(NULL, metric_names))
  obs_tab      <- table(sp_obs[sp_obs %in% names(abund_f)])
  obs_species  <- names(obs_tab)
  n_ind        <- as.numeric(obs_tab)
  w_obs        <- as.numeric(abund_f[obs_species])
  pool_species <- unique(sp_pool)
  for (i in 1:NPERM) {
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
pool_results <- merge(pool_results, n_sp_tab, by.x = FOCAL_COL, by.y = "FOCAL", all.x = TRUE)
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
  in_focal <- dat$FOCAL == f
  in_pool  <- dat$POOL  == indiv_results$POOL[r]

  traits_obs  <- dat$log_dist_cm[in_focal];  sp_obs  <- dat$scientificName_Species[in_focal]
  traits_pool <- dat$log_dist_cm[in_pool];   sp_pool <- dat$scientificName_Species[in_pool]

  # species (and thus true abundances) are preserved, so use the observed weights
  abund_f <- abund_for(f, sp_obs)
  obs <- community_metrics(traits_obs, sp_obs, abund_f)

  null_mat <- matrix(NA, nrow = NPERM, ncol = length(metric_names),
                     dimnames = list(NULL, metric_names))
  obs_species <- unique(sp_obs)
  for (i in 1:NPERM) {
    # keep species & counts; redraw each species' individuals from the pool
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
indiv_results <- merge(indiv_results, n_sp_tab, by.x = FOCAL_COL, by.y = "FOCAL", all.x = TRUE)
write.csv(indiv_results, paste0("./Outputs/", OUT_PREFIX, "_IndividualNull.csv"), row.names = FALSE)
message("wrote ", OUT_PREFIX, "_IndividualNull.csv  (", nrow(indiv_results), " focal units)")


#### 7. NULL MODEL 3: SWAP MEANS (within-community) ####
# Within each focal community, keep every species' abundance and its within-
# species deviations, but relocate each species onto a randomly permuted
# community mean. Needs NO regional pool, so it does not depend on POOL except
# through which focal units survive the MIN_POOL_UNITS screen; we run it on the
# same focal_units as sections 5-6 so the three files line up.
#
# Singletons are KEPT (no >= 2 filter): a one-specimen species has a within-
# species deviation of 0 and is simply relocated to a permuted mean. sdnnd is
# invariant under this null by construction (the SET of means is only permuted),
# so it reports CI = obs, ses = NA -- expected, not a bug.
set.seed(SEED)

swap_results <- data.frame(FOCAL = focal_units, stringsAsFactors = FALSE)
swap_results$POOL <- focal_pool$POOL[match(swap_results$FOCAL, focal_pool$FOCAL)]
for (m in metric_names) for (s in c("_obs","_lower","_upper","_ses")) swap_results[[paste0(m, s)]] <- NA_real_
for (m in metric_names) swap_results[[paste0(m, "_dir")]] <- NA_character_

for (r in seq_along(focal_units)) {
  f <- focal_units[r]
  in_focal <- dat$FOCAL == f
  traits_obs <- dat$log_dist_cm[in_focal];  sp_obs <- dat$scientificName_Species[in_focal]

  # species (and thus true abundances) are preserved, so use the observed weights
  abund_f <- abund_for(f, sp_obs)
  obs <- community_metrics(traits_obs, sp_obs, abund_f)

  null_mat <- matrix(NA, nrow = NPERM, ncol = length(metric_names),
                     dimnames = list(NULL, metric_names))

  # swap operates on the community the metrics use: finite traits + a true
  # abundance. NO >= 2 filter, so singletons participate (deviation 0).
  keep <- is.finite(traits_obs) & !is.na(sp_obs) & sp_obs %in% names(abund_f)
  tr   <- traits_obs[keep]; spp <- sp_obs[keep]

  if (length(unique(spp)) >= 2) {
    sp_f  <- factor(spp)                             # index species by position, not by name
    means <- as.numeric(tapply(tr, sp_f, mean))      # one mean per level, in level order
    codes <- as.integer(sp_f)                        # each individual's species code
    devs  <- tr - means[codes]                       # deviation from own mean (0 for singletons)
    for (i in 1:NPERM) {
      means_swapped <- sample(means)                 # permute community means across species
      traits_null   <- devs + means_swapped[codes]
      null_mat[i, ] <- community_metrics(traits_null, spp, abund_f)
    }
  }

  for (m in metric_names) {
    o  <- obs[m]
    nd <- null_mat[, m]; nd <- nd[is.finite(nd)]
    swap_results[r, paste0(m, "_obs")] <- o
    if (length(nd) >= 2 && is.finite(o)) {
      if (sd(nd) > 1e-9) {
        lo <- as.numeric(quantile(nd, NULLQS[1])); hi <- as.numeric(quantile(nd, NULLQS[2]))
        swap_results[r, paste0(m, "_lower")] <- lo
        swap_results[r, paste0(m, "_upper")] <- hi
        swap_results[r, paste0(m, "_ses")]   <- (o - mean(nd)) / sd(nd)
        swap_results[r, paste0(m, "_dir")]   <- if (o < lo) "lower" else if (o > hi) "higher" else "neutral"
      } else {
        swap_results[r, paste0(m, "_lower")] <- o
        swap_results[r, paste0(m, "_upper")] <- o
        swap_results[r, paste0(m, "_ses")]   <- NA
        swap_results[r, paste0(m, "_dir")]   <- "neutral"
      }
    }
  }
}

names(swap_results)[1:2] <- c(FOCAL_COL, POOL_COL)
swap_results <- merge(swap_results, lat, by.x = FOCAL_COL, by.y = "FOCAL", all.x = TRUE)
swap_results <- merge(swap_results, n_sp_tab, by.x = FOCAL_COL, by.y = "FOCAL", all.x = TRUE)
write.csv(swap_results, paste0("./Outputs/", OUT_PREFIX, "_SwapMeansNull.csv"), row.names = FALSE)
message("wrote ", OUT_PREFIX, "_SwapMeansNull.csv  (", nrow(swap_results), " focal units)")
