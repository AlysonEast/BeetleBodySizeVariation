#### SIGNIFICANCE-FLAG (_dir) PATTERNS ACROSS NULLS, POOL EXTENT, AND SCALE ####
# -----------------------------------------------------------------------------
# Companion to Overlap_CustomNulls_ByYear.R. Each null-output file carries, per
# metric, a direction flag <metric>_dir taking "lower" / "neutral" / "higher":
# whether the observed value falls below, inside, or above the null 2.5-97.5% CI.
# "neutral" = indistinguishable from the null (not significant); "lower"/"higher"
# = a significant departure, whose sign is the ecological signal (e.g. for
# overlap: lower = overdispersion, higher = clustering).
#
# This script ignores the effect sizes and asks only where the SIGNIFICANCE lands:
#   - do the two null models (PoolNull vs IndividualNull) flag the same units?
#   - does the flag distribution shift as the regional pool widens (site -> domain
#     -> all)?
#   - does it differ between scales (plot vs site)?
#
# Flags are AGGREGATED ACROSS YEARS: both years are pooled, so each observation is
# a focal-unit-year and cell counts are the summed site-years. min_logratio is
# dropped per request; four metrics are kept.
#
# A SEPARATE set of figures + tables is written for augmented and unaugmented data
# (looped over AUGS); the two are never combined into one figure. Filenames embed
# the augmentation state and the pooled-year tag.
#
# Flat / stepwise: run 0-2 once (setup + load), then section 3 or 4 on its own.
# Each producing section loops AUGS internally; set AUG <- "aug" to run one block.
# -----------------------------------------------------------------------------

library(ggplot2)
library(ggpubr)     # theme_pubr
library(svglite)    # svg device
library(dplyr)
library(tidyr)

setwd("/home/aly/Beetles/BeetleBodySizeVariation")


#### 0. SETTINGS ####
OUT_DIR <- "./Outputs"
FIG_DIR <- "./Figures/DirPatterns"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# min_logratio intentionally excluded
METRICS <- c("overlap_norm", "overlap_unnorm", "niche_range", "sdnnd")
NULLS   <- c("PoolNull", "IndividualNull")

USE_YEARS <- c(2018, 2019)              # flags are pooled (aggregated) across these
AUGS      <- c("aug", "noaug")          # one full figure+table set is written per aug
YEAR_TAG  <- paste(USE_YEARS, collapse = "-")

# pool valid per scale; ordered narrow -> wide for the x-axis
POOLS_FOR   <- list(plot = c("site", "domain", "all"), site = c("domain", "all"))
POOL_LEVELS <- c("site", "domain", "all")

DIR_LEVELS <- c("lower", "neutral", "higher")
DIR_COLS   <- c(lower = "#2C7FB8", neutral = "grey80", higher = "#D95F02")
DIR_LABS   <- c(lower = "lower than null", neutral = "neutral (n.s.)", higher = "higher than null")
SCALE_LABS <- c(plot = "Plot", site = "Site")

SAVE_SVG <- FALSE
DPI      <- 300
JITTER   <- c(plot = 1.0, site = 0.6)   # degrees of lon/lat jitter per scale (maps); plot spread most


#### 1. LOADER  (mirrors read_null_long in Overlap_CustomNulls_Compare.R) ####
# Reads one null-output file and returns it LONG: one row per focal unit x metric,
# carrying the direction flag, the SES, and the raw observed value (SES and obs
# feed the maps in sections 5-6).
read_dir_long <- function(level, pool, aug, year, null) {
  f <- file.path(OUT_DIR, sprintf("%s_by_%s_%s_%s_%s.csv", level, pool, aug, year, null))
  if (!file.exists(f)) { message("  [skip] missing ", basename(f)); return(NULL) }
  d <- read.csv(f, stringsAsFactors = FALSE)
  focal_col <- if (level == "plot") "plotID" else "siteID"
  long <- do.call(rbind, lapply(METRICS, function(m) data.frame(
    focal  = d[[focal_col]],
    metric = m,
    dir    = d[[paste0(m, "_dir")]],
    ses    = d[[paste0(m, "_ses")]],
    obs    = d[[paste0(m, "_obs")]],
    stringsAsFactors = FALSE)))
  long$level <- level; long$pool <- pool; long$aug <- aug
  long$year  <- as.character(year); long$null <- null
  long
}


#### 2. LOAD BOTH AUGMENTATIONS x BOTH YEARS INTO ONE LONG TABLE ####
dat_list <- list()
for (aug in AUGS) for (yr in USE_YEARS) for (nl in NULLS)
  for (level in names(POOLS_FOR)) for (pool in POOLS_FOR[[level]])
    dat_list[[length(dat_list) + 1]] <- read_dir_long(level, pool, aug, yr, nl)
dat <- do.call(rbind, dat_list)

# keep only computable flags, and lock factor orders for consistent plots/tables
dat <- subset(dat, !is.na(dir) & dir %in% DIR_LEVELS)
dat$dir    <- factor(dat$dir,    levels = DIR_LEVELS)
dat$metric <- factor(dat$metric, levels = METRICS)
dat$pool   <- factor(dat$pool,   levels = POOL_LEVELS)
dat$null   <- factor(dat$null,   levels = NULLS)
dat$level  <- factor(dat$level,  levels = c("plot", "site"))
message("loaded ", nrow(dat), " focal x metric flags  (years ", YEAR_TAG,
        ", aug: ", paste(AUGS, collapse = "/"), ")")


#### 3. FLAG PROPORTIONS + SIGNIFICANT FRACTION  (one set per augmentation) ####
# summ / props table / stacked composition / significant-fraction lines, all
# pooled over both years. Loops AUGS -> two independent output sets.
for (AUG in AUGS) {
  d <- subset(dat, aug == AUG)
  if (nrow(d) == 0) { message("[skip] no data for aug = ", AUG); next }

  ## per-cell flag counts and proportions (years pooled: n = summed site-years)
  summ <- d %>%
    count(metric, null, level, pool, dir, .drop = FALSE) %>%
    group_by(metric, null, level, pool) %>%
    mutate(n_total = sum(n), prop = ifelse(n_total > 0, n / n_total, NA)) %>%
    ungroup() %>%
    filter(n_total > 0)                   # drop degenerate scale x pool cells (e.g. site x site)

  props_wide <- summ %>%
    select(metric, null, level, pool, n_total, dir, n, prop) %>%
    pivot_wider(names_from = dir, values_from = c(n, prop), values_fill = 0) %>%
    mutate(prop_sig = prop_lower + prop_higher) %>%
    arrange(metric, level, null, pool)
  write.csv(props_wide,
            file.path(OUT_DIR, sprintf("DirPatterns_Proportions_%s_%s.csv", AUG, YEAR_TAG)),
            row.names = FALSE)

  ## 3a. stacked flag proportions: metric x (scale x null), bars across pool.
  ## free_x drops the pools a scale doesn't have; n (site-years) above each bar.
  n_lab <- distinct(summ, metric, null, level, pool, n_total)
  p <- ggplot(summ, aes(pool, prop, fill = dir)) +
    geom_col(width = 0.85, position = position_stack(reverse = TRUE)) +
    geom_text(data = n_lab, aes(pool, 1.03, label = n_total),
              inherit.aes = FALSE, size = 2.5, vjust = 0, colour = "grey30") +
    facet_grid(metric ~ level + null, scales = "free_x", space = "free_x",
               labeller = labeller(level = SCALE_LABS)) +
    scale_fill_manual(values = DIR_COLS, labels = DIR_LABS, name = "Observed vs null") +
    scale_y_continuous(limits = c(0, 1.1), breaks = c(0, .5, 1), expand = c(0, 0)) +
    labs(x = "Regional pool", y = "Proportion of focal-years",
         title = sprintf("Significance-flag composition  (years %s, %s)", YEAR_TAG, AUG)) +
    theme_pubr(base_size = 11, legend = "right") +
    theme(strip.background = element_blank(), strip.text = element_text(size = 9),
          panel.spacing.x = unit(4, "pt"))
  ggsave(file.path(FIG_DIR, sprintf("DirProps_stacked_%s_%s.png", AUG, YEAR_TAG)),
         p, width = 10, height = 8, dpi = DPI)
  if (SAVE_SVG) ggsave(file.path(FIG_DIR, sprintf("DirProps_stacked_%s_%s.svg", AUG, YEAR_TAG)),
                       p, width = 10, height = 8)

  ## 3b. significant fraction (lower + higher) across pool: two nulls, both scales
  sig <- summ %>%
    group_by(metric, null, level, pool) %>%
    summarise(prop_sig = sum(prop[dir != "neutral"]), n_total = first(n_total), .groups = "drop")
  p2 <- ggplot(sig, aes(pool, prop_sig, colour = null,
                        group = interaction(null, level), shape = level, linetype = level)) +
    geom_line() +
    geom_point(size = 2.4) +
    facet_wrap(~ metric, nrow = 1) +
    scale_shape_manual(values = c(plot = 16, site = 17), labels = SCALE_LABS, name = "Scale") +
    scale_linetype_manual(values = c(plot = 1, site = 2), labels = SCALE_LABS, name = "Scale") +
    scale_colour_manual(values = c(PoolNull = "#1B9E77", IndividualNull = "#7570B3"), name = "Null") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Regional pool (narrow -> wide)", y = "Fraction significant (lower + higher)",
         title = sprintf("Significant fraction across pool extent  (years %s, %s)", YEAR_TAG, AUG)) +
    theme_pubr(base_size = 11, legend = "right") +
    theme(strip.background = element_blank())
  ggsave(file.path(FIG_DIR, sprintf("DirSigFraction_%s_%s.png", AUG, YEAR_TAG)),
         p2, width = 11, height = 3.8, dpi = DPI)
  if (SAVE_SVG) ggsave(file.path(FIG_DIR, sprintf("DirSigFraction_%s_%s.svg", AUG, YEAR_TAG)),
                       p2, width = 11, height = 3.8)

  message("  [", AUG, "] wrote proportions CSV + stacked + sig-fraction figures")
}


#### 4. DO THE TWO NULLS AGREE?  (paired within focal-unit-year; one set per aug) ####
# PoolNull and IndividualNull are evaluated on the SAME focal units, so this is a
# PAIRED comparison. Join the two flags per focal-unit-YEAR (year in the key so a
# unit's 2018 flags pair with its 2018 flags), then aggregate over years:
#   (a) agreement per (metric, scale, pool); (b) flag x flag confusion per metric.
for (AUG in AUGS) {
  d <- subset(dat, aug == AUG)
  if (nrow(d) == 0) { message("[skip] no data for aug = ", AUG); next }

  paired <- inner_join(
    d %>% filter(null == "PoolNull")       %>% select(level, pool, year, focal, metric, dir_pool = dir),
    d %>% filter(null == "IndividualNull") %>% select(level, pool, year, focal, metric, dir_indiv = dir),
    by = c("level", "pool", "year", "focal", "metric"))

  ## (a) agreement table (pooled over years)
  null_agree <- paired %>%
    group_by(metric, level, pool) %>%
    summarise(n_paired = n(),
              prop_same      = mean(dir_pool == dir_indiv),
              prop_pool_sig  = mean(dir_pool  != "neutral"),
              prop_indiv_sig = mean(dir_indiv != "neutral"),
              .groups = "drop") %>%
    mutate(across(c(prop_same, prop_pool_sig, prop_indiv_sig), ~ round(.x, 3))) %>%
    arrange(metric, level, pool)
  write.csv(null_agree,
            file.path(OUT_DIR, sprintf("DirPatterns_NullAgreement_%s_%s.csv", AUG, YEAR_TAG)),
            row.names = FALSE)
  message("  [", AUG, "] null agreement:"); print(null_agree)

  ## (b) confusion of the two nulls' flags, per metric (counts pooled over scale x pool x year)
  conf <- paired %>%
    count(metric, dir_pool, dir_indiv, .drop = FALSE)   # factor grouping fills 0-count cells
  p3 <- ggplot(conf, aes(dir_indiv, dir_pool, fill = n)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = n), size = 3) +
    facet_wrap(~ metric, nrow = 1) +
    scale_fill_gradient(low = "grey95", high = "#2C7FB8", name = "focal-years") +
    scale_x_discrete(labels = DIR_LABS) +
    scale_y_discrete(labels = DIR_LABS) +
    labs(x = "IndividualNull flag", y = "PoolNull flag",
         title = sprintf("Null-model flag agreement  (years %s, %s; pooled over scale x pool)",
                         YEAR_TAG, AUG)) +
    theme_pubr(base_size = 11, legend = "right") +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank())
  ggsave(file.path(FIG_DIR, sprintf("DirNullConfusion_%s_%s.png", AUG, YEAR_TAG)),
         p3, width = 11, height = 3.6, dpi = DPI)
  if (SAVE_SVG) ggsave(file.path(FIG_DIR, sprintf("DirNullConfusion_%s_%s.svg", AUG, YEAR_TAG)),
                       p3, width = 11, height = 3.6)

  message("  [", AUG, "] wrote null-agreement CSV + confusion figure")
}


#### 5-6. MAPS SETUP: coordinates from neonDivData + shared map builder ####
# Coordinates are appended from neonDivData::neon_location beetle ("bet") sampling
# locations, as in the other scripts; site points are the centroid of their beetle
# plots. Each point is a focal-unit-YEAR -- both years are shown (consistent with
# the pooled-years treatment above) rather than averaged -- with shape = the
# direction flag (circle = neutral/n.s., triangles = significant lower/higher).
# Heavy spatial overlap (plots within a site, the two years at a location) is
# spread with jitter, larger for the denser plot scale (JITTER in settings).
library(neonDivData)
library(maps)

## coordinate lookups from neonDivData (beetle "bet" locations only)
bet_loc <- subset(neon_location,
                  substr(location_id, nchar(location_id) - 2, nchar(location_id)) == "bet")
plot_xy <- unique(bet_loc[, c("plotID", "longitude", "latitude")])
plot_xy <- plot_xy[!is.na(plot_xy$plotID) & !is.na(plot_xy$longitude), ]
plot_xy$siteID <- sub("_.*", "", plot_xy$plotID)                 # HARV_001 -> HARV
site_xy <- aggregate(cbind(longitude, latitude) ~ siteID, plot_xy, mean)   # site = plot centroid

attach_xy <- function(df, level) {
  if (level == "plot")
    merge(df, plot_xy[, c("plotID", "longitude", "latitude")],
          by.x = "focal", by.y = "plotID", all.x = TRUE)
  else
    merge(df, site_xy, by.x = "focal", by.y = "siteID", all.x = TRUE)
}

states_map <- map_data("state")
CONUS      <- list(xlim = c(-125, -66), ylim = c(24, 50))
DIR_SHAPES <- c(lower = 25, neutral = 21, higher = 24)   # triangle-down / circle / triangle-up

# Shared builder: base CONUS polygons + jittered points (fill = `fillvar`, shape =
# dir) + the caller's fill scale and facet. `fillvar` is set by the caller so the
# same function serves the SES maps (section 5) and the raw-value maps (section 6).
# Reused across both sections, so extracted per the one-helper-per-reuse rule.
make_map <- function(d, jit, facet, fill_scale, title) {
  set.seed(517)                                          # reproducible jitter
  ggplot() +
    geom_polygon(data = states_map, aes(long, lat, group = group),
                 fill = "grey95", colour = "grey75", linewidth = 0.2) +
    geom_jitter(data = d, aes(longitude, latitude, fill = fillvar, shape = dir),
                width = jit, height = jit, size = 2, stroke = 0.2,
                colour = "grey30", alpha = 0.9) +
    facet +
    fill_scale +
    scale_shape_manual(values = DIR_SHAPES, labels = DIR_LABS, name = "Observed vs null") +
    guides(shape = guide_legend(override.aes = list(fill = "grey50"))) +
    coord_quickmap(xlim = CONUS$xlim, ylim = CONUS$ylim) +
    labs(title = title, x = NULL, y = NULL) +
    theme_pubr(base_size = 11, legend = "right") +
    theme(strip.background = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.spacing = unit(6, "pt"))
}

# focal-unit-year table for one scale x aug, with coordinates attached and
# non-CONUS points reported. LV (not `level`) so it can't shadow the column.
map_frame <- function(AUG, LV) {
  d <- subset(dat, aug == AUG & level == LV & !is.na(ses))
  d <- attach_xy(d, LV)
  d <- d[!is.na(d$longitude), ]
  if (nrow(d)) {
    n_out <- sum(d$longitude < CONUS$xlim[1] | d$longitude > CONUS$xlim[2] |
                 d$latitude  < CONUS$ylim[1] | d$latitude  > CONUS$ylim[2])
    if (n_out > 0)
      message("  [", LV, " ", AUG, "] ", n_out,
              " focal-year points outside the CONUS frame (AK/HI/PR) not drawn")
  }
  d
}


#### 5. SES MAPS: fill = SES (diverging), facets = null x metric ####
# SES is standardized and comparable across metrics, so all four share one
# diverging fill scale in a single null x metric figure per scale x aug.
SES_CLIP <- 4                                            # symmetric fill limit; |SES| beyond squished
for (AUG in AUGS) for (LV in names(POOLS_FOR)) {
  d <- map_frame(AUG, LV)
  if (nrow(d) == 0) { message("[skip] no mappable rows: ", LV, " ", AUG); next }
  d$fillvar <- d$ses
  p <- make_map(d, JITTER[[LV]], facet_grid(null ~ metric),
                scale_fill_gradient2(low = "#2C7FB8", mid = "grey92", high = "#D95F02",
                                     midpoint = 0, limits = c(-SES_CLIP, SES_CLIP),
                                     oob = scales::squish, name = "SES"),
                sprintf("%s-level SES map  (years %s, %s)", SCALE_LABS[[LV]], YEAR_TAG, AUG))
  ggsave(file.path(FIG_DIR, sprintf("DirMapSES_%s_%s_%s.png", LV, AUG, YEAR_TAG)),
         p, width = 12, height = 6, dpi = DPI)
  if (SAVE_SVG) ggsave(file.path(FIG_DIR, sprintf("DirMapSES_%s_%s_%s.svg", LV, AUG, YEAR_TAG)),
                       p, width = 12, height = 6)
  message("  [", LV, " ", AUG, "] wrote SES map (", nrow(d), " focal-year points)")
}


#### 6. RAW-VALUE MAPS: fill = observed value (sequential), facets = null ####
# Raw values are NOT comparable across metrics (overlap in [0,1], sdnnd tiny), so
# each metric gets its own figure and its own sequential fill scale spanning that
# metric's range; the two null panels sit side by side. One figure per
# metric x scale x aug.
for (AUG in AUGS) for (LV in names(POOLS_FOR)) {
  base <- map_frame(AUG, LV)
  if (nrow(base) == 0) { message("[skip] no mappable rows: ", LV, " ", AUG); next }
  for (m in METRICS) {
    d <- subset(base, metric == m & !is.na(obs))
    if (nrow(d) == 0) { message("  [skip] no obs for ", m, " ", LV, " ", AUG); next }
    d$fillvar <- d$obs
    p <- make_map(d, JITTER[[LV]], facet_wrap(~ null, nrow = 1),
                  scale_fill_viridis_c(name = m, option = "viridis"),
                  sprintf("%s  |  %s-level raw value  (years %s, %s)",
                          m, SCALE_LABS[[LV]], YEAR_TAG, AUG))
    ggsave(file.path(FIG_DIR, sprintf("DirMapRaw_%s_%s_%s_%s.png", m, LV, AUG, YEAR_TAG)),
           p, width = 9, height = 4.2, dpi = DPI)
    if (SAVE_SVG) ggsave(file.path(FIG_DIR, sprintf("DirMapRaw_%s_%s_%s_%s.svg", m, LV, AUG, YEAR_TAG)),
                         p, width = 9, height = 4.2)
  }
  message("  [", LV, " ", AUG, "] wrote ", length(METRICS), " raw-value maps")
}
