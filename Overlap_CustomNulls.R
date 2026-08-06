#### CUSTOM NULL MODELS FOR OVERLAP  (flat / stepwise version) ####
# -------------------------------------------------------------------------
# This script tests body-size community assembly against two custom null models,
# each evaluated with five metrics. It is written to be run top-to-bottom, but
# also stepwise: run sections 0-4 once (setup), then run either null section
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
# FIVE METRICS (each: observed value, null lower/upper CI, SES, direction flag):
#   overlap_norm    overlap, each density normalized to area 1
#   overlap_unnorm  overlap, density scaled by abundance
#   niche_range     width of occupied trait space (2.5-97.5% span; NOT overlap)
#   sdnnd           SD of nearest-neighbour distances between species means
#                   (LOW = even spacing = limiting-similarity signature)
#   min_logratio    smallest adjacent gap between species means on the log10 axis
#                   (= log10 of the tightest size ratio; compare to log10(1.3)=0.114)
# -------------------------------------------------------------------------

library(Ostats)   # community_overlap()
library(dplyr)
library(tidyr)
library(dplyr)

setwd("/home/aly/Beetles/BeetleBodySizeVariation")

# abundance-weighted overlap (weights by true, effort-scaled abundance instead of
# the augmented observation counts). Kept in its own script for readability.
source("./community_overlap_weighted.R")


#### 0. SETTINGS ####
LEVEL <- "site"     # "plot" or "site"
POOL  <- "all"     # "site" or "domain" or "all"

NPERM   <- 99
NULLQS  <- c(0.025, 0.975)
HUTCH   <- 1.3               # Hutchinsonian ratio, for reference on min_logratio
SEED    <- 517
MIN_POOL_UNITS <- 2          # drop focal units whose pool holds fewer focal units
# (site->domain: drops single-site domains;
#  plot->site:  drops single-plot sites)

metric_names <- c("overlap_norm", "overlap_unnorm", "niche_range", "sdnnd", "min_logratio")

# translate the LEVEL / POOL choices into column names and the CV scale to use
if (LEVEL == "plot") { FOCAL_COL <- "plotID"; CV_SCALE <- "Plot Level" }
if (LEVEL == "site") { FOCAL_COL <- "siteID"; CV_SCALE <- "Site Level" }
POOL_COL   <- if (POOL == "site") "siteID" else 
  if (POOL == "domain") "domainID" else "all"
OUT_PREFIX <- paste0(LEVEL, "_by_", POOL)     # e.g. "plot_by_site"


#### 1. READ CLEAN DATA ####
all_elytra <- read.csv("./Data/BodysizeCombinedClean.csv")
all_elytra <- subset(all_elytra, yearCollected == 2018 | yearCollected == 2019)
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


#### 2. AUGMENT SPARSE SPECIES x FOCAL CELLS TO n >= 20 ####
# Same lognormal augmentation as Overlap_Plot.R / Overlap_Site.R, grouped by the
# focal unit for this level.
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

aug <- rbind(all_elytra[, c("scientificName_Species", FOCAL_COL, "cm_elytra_max_length")],
             as.data.frame(sim_low_n))
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
abund_tab <- read.csv(if (LEVEL == "plot") "./Data/plot_abund.csv" else "./Data/site_abund.csv")
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

# coverage report: observed species with no true abundance (dropped from the weighted overlap)
obs_keys     <- unique(paste(aug$FOCAL, aug$scientificName_Species, sep = "|"))
missing_keys <- obs_keys[is.na(abund_lookup[obs_keys])]
if (length(missing_keys))
  message("NOTE: ", length(missing_keys), " observed (focal|species) combos have no true abundance ",
          "and drop from the weighted metrics. e.g. ", paste(missing_keys, collapse = "; "))


#### HELPER: the five metrics for one community ####
# The only function in the script. It is called once for the observed community
# and once for every null draw, so the metric definitions live in exactly one
# place. Returns a named vector; NA if fewer than two eligible species.
# `abund` is a named vector (species -> true, effort-scaled abundance) used to
# weight the overlaps. Species without a true abundance are dropped from ALL
# metrics so every focal unit describes one consistent community (coverage is
# reported up front in section 3b).
community_metrics <- function(traits, sp, abund) {
  traits <- as.numeric(traits); sp <- as.character(sp)
  ok <- is.finite(traits) & !is.na(sp) & sp != "" & sp %in% names(abund)   # need a trait & a true abundance
  traits <- traits[ok]; sp <- sp[ok]
  eligible <- names(which(table(sp) >= 2))          # >=2 individuals, as community_overlap needs
  traits <- traits[sp %in% eligible]; sp <- sp[sp %in% eligible]
  
  out <- c(overlap_norm = NA, overlap_unnorm = NA, niche_range = NA,
           sdnnd = NA, min_logratio = NA)
  if (length(unique(sp)) < 2) return(out)
  
  # overlaps weighted by TRUE abundance (not the augmented observation counts)
  out["overlap_norm"]   <- community_overlap_weighted(traits, sp, abund, normal = TRUE)
  out["overlap_unnorm"] <- community_overlap_weighted(traits, sp, abund, normal = FALSE)
  
  # width of occupied trait space (robust 2.5-97.5% span)
  out["niche_range"] <- diff(quantile(traits, c(0.025, 0.975)))
  
  # spacing of species means on the log10 axis
  means <- sort(tapply(traits, sp, mean))
  gaps  <- diff(means)                              # adjacent gaps = log10 size ratios
  nn    <- pmin(c(gaps, Inf), c(Inf, gaps))         # nearest-neighbour distance per species
  out["sdnnd"]        <- sd(nn)
  out["min_logratio"] <- min(gaps)
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
  # null (sd ~ 0, e.g. spacing metrics under swap_means) report CI = obs, ses = NA
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
write.csv(indiv_results, paste0("./Outputs/", OUT_PREFIX, "_IndividualNull.csv"), row.names = FALSE)
message("wrote ", OUT_PREFIX, "_IndividualNull.csv  (", nrow(indiv_results), " focal units)")


# #### 7. NULL MODEL 3: SWAP MEANS (within-community, out-of-the-box Ostats null) ####
# # The mean-swap null from Ostats::Ostats(swap_means = TRUE): within each focal
# # community, keep every species' abundance and its within-species deviations, but
# # relocate each species onto a randomly permuted community mean. This needs NO
# # regional pool (it is a within-community null), so it runs on every focal unit
# # and its result does not depend on POOL -- the output is named by LEVEL only.
# #
# # NOTE: because the SET of species means is only permuted (never changed), the two
# # spacing metrics are invariant by construction -- sdnnd and min_logratio have
# # null CI = observed, ses = NA, dir = "neutral". That is expected, not a bug: this
# # null is only informative for the overlap metrics (and weakly niche_range). We
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
#   # swap operates on exactly the community the metrics use: finite traits, species
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
#   # summarise observed vs null per metric; sdnnd and min_logratio are invariant
#   # under swap_means (sd ~ 0) so they report CI = obs, ses = NA (expected, not a bug)
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
# write.csv(swap_results, paste0("./Outputs/", LEVEL, "_SwapMeansNull.csv"), row.names = FALSE)
# message("wrote ", LEVEL, "_SwapMeansNull.csv  (", nrow(swap_results), " focal units)")
