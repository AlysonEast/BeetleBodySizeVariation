#### SIGNIFICANCE FLAGS, SES, AND MAPS ACROSS NULL MODELS  (niche_range + overlap_depth) ####
# -----------------------------------------------------------------------------
# Companion to 6a_CommunityShape_CustomNulls_ByYear.R, which writes one file per
#   <LEVEL>_by_<POOL>_<YEAR>_<PoolNull|IndividualNull|SwapMeansNull>.csv
# Each carries, per metric, <metric>_obs / _ses / _dir, where _dir is
# "lower" / "neutral" / "higher": whether the observed value falls below, inside,
# or above the null 2.5-97.5% CI. "neutral" = not significant; "lower"/"higher"
# = a significant departure, whose sign is the ecological signal:
#   overlap_depth  higher = clustering (filtering-consistent),
#                  lower  = overdispersion (competition-consistent)
#   niche_range    higher = wider occupied trait space than the null,
#                  lower  = contracted trait space
#
# Only the two analysis metrics are kept (niche_range, overlap_depth). There is
# no augmentation dimension any more, so every section writes ONE set of outputs.
#
# YEARS: for the flag and SES summaries (sections 3-4), both years are pooled, so
# each observation is a focal-unit-year. For the maps (sections 5-6), the two
# years are averaged to one point per focal unit (see section 5 setup).
#
#   3  Flag composition (stacked bars) + SES distributions, by null x pool x scale
#   4  Do the three nulls agree? pairwise agreement table + confusion tiles
#   5  Observed-value maps (null- and pool-invariant: one map per metric x scale)
#   6  SES maps (rows = null, cols = metric) at the headline pool per scale
#
# Swap-means needs no pool, so its bars/boxes change across pool only through
# the MIN_POOL_UNITS screen in 6a. It is a weak test for niche_range (permuting
# species means barely moves the 2.5-97.5% span), so expect mostly neutral there.
#
# Flat / stepwise: run 0-2 once (setup + load), then any later section on its own.
# -----------------------------------------------------------------------------

library(ggplot2)
library(ggpubr)     # theme_pubr
library(svglite)    # svg device
library(patchwork)  # side-by-side observed-value maps
library(dplyr)
library(tidyr)

setwd("/home/aly/Beetles/BeetleBodySizeVariation")


#### 0. SETTINGS ####
OUT_DIR <- "./Outputs"
FIG_DIR <- "./Figures/DirPatterns"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

METRICS     <- c("niche_range", "overlap_depth")
METRIC_LABS <- c(niche_range = "Niche range", overlap_depth = "Overlap depth")
NULLS       <- c("PoolNull", "IndividualNull", "SwapMeansNull")
NULL_LABS   <- c(PoolNull = "Pool null", IndividualNull = "Individual null",
                 SwapMeansNull = "Swap-means null")
NULL_COLS   <- c(PoolNull = "#1B9E77", IndividualNull = "#7570B3", SwapMeansNull = "#E7298A")

USE_YEARS <- c(2018, 2019)
YEAR_TAG  <- paste(USE_YEARS, collapse = "-")

# pools valid per scale; ordered narrow -> wide for the x-axis
POOLS_FOR   <- list(plot = c("site", "domain", "all"), site = c("domain", "all"))
POOL_LEVELS <- c("site", "domain", "all")

# headline pool per scale (the typical pairings), used for the confusion tiles
# (section 4) and the maps (sections 5-6)
FOCUS_POOL <- c(plot = "site", site = "domain")

DIR_LEVELS <- c("lower", "neutral", "higher")
DIR_COLS   <- c(lower = "#2C7FB8", neutral = "grey80", higher = "#D95F02")
DIR_LABS   <- c(lower = "lower than null", neutral = "neutral (n.s.)", higher = "higher than null")
SCALE_LABS <- c(plot = "Plot", site = "Site")

SES_CLIP <- 4                          # symmetric SES fill limit on maps; beyond is squished
JITTER   <- c(plot = 0.6, site = 0)    # degrees of lon/lat jitter (plots within a site overlap)
SAVE_SVG <- FALSE
DPI      <- 300


#### 1. HELPERS ####

# read one null-output file and return it LONG: one row per focal unit x metric,
# with the direction flag, SES, and observed value
read_dir_long <- function(level, pool, year, null) {
  f <- file.path(OUT_DIR, sprintf("%s_by_%s_%s_%s.csv", level, pool, year, null))
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
  long$level <- level; long$pool <- pool
  long$year  <- as.character(year); long$null <- null
  long
}

# save a figure as png (+ svg if SAVE_SVG)
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = DPI)
  if (SAVE_SVG) ggsave(file.path(FIG_DIR, paste0(name, ".svg")), p, width = w, height = h)
  message("  wrote ", name)
}


#### 2. LOAD BOTH YEARS INTO ONE LONG TABLE ####
dat_list <- list()
for (yr in USE_YEARS) for (nl in NULLS)
  for (level in names(POOLS_FOR)) for (pool in POOLS_FOR[[level]])
    dat_list[[length(dat_list) + 1]] <- read_dir_long(level, pool, yr, nl)
dat <- do.call(rbind, dat_list)

# keep only computable flags (single-species units carry NA), lock factor orders
dat <- subset(dat, !is.na(dir) & dir %in% DIR_LEVELS)
dat$dir    <- factor(dat$dir,    levels = DIR_LEVELS)
dat$metric <- factor(dat$metric, levels = METRICS)
dat$pool   <- factor(dat$pool,   levels = POOL_LEVELS)
dat$null   <- factor(dat$null,   levels = NULLS)
dat$level  <- factor(dat$level,  levels = c("plot", "site"))
message("loaded ", nrow(dat), " focal-year x metric flags  (years ", YEAR_TAG, ")")

facet_labs <- labeller(level = SCALE_LABS, null = NULL_LABS, metric = METRIC_LABS)


#### 3. FLAG COMPOSITION + SES DISTRIBUTIONS ####
## per-cell flag counts and proportions (years pooled: n = summed focal-years)
summ <- dat %>%
  count(metric, null, level, pool, dir, .drop = FALSE) %>%
  group_by(metric, null, level, pool) %>%
  mutate(n_total = sum(n), prop = ifelse(n_total > 0, n / n_total, NA)) %>%
  ungroup() %>%
  filter(n_total > 0)                     # drop degenerate scale x pool cells (site x site)

## summary table: flag proportions + median SES per cell
med_ses <- dat %>%
  group_by(metric, null, level, pool) %>%
  summarise(median_ses = median(ses, na.rm = TRUE), .groups = "drop")
summary_tab <- summ %>%
  select(metric, null, level, pool, n_total, dir, prop) %>%
  pivot_wider(names_from = dir, names_prefix = "prop_", values_from = prop, values_fill = 0) %>%
  mutate(prop_sig = prop_lower + prop_higher) %>%
  left_join(med_ses, by = c("metric", "null", "level", "pool")) %>%
  mutate(across(starts_with("prop_") | matches("median_ses"), ~ round(.x, 3))) %>%
  arrange(metric, level, null, pool)
write.csv(summary_tab, file.path(OUT_DIR, sprintf("DirPatterns_Summary_%s.csv", YEAR_TAG)),
          row.names = FALSE)
print(summary_tab)

## 3a. stacked flag composition: rows = metric, cols = scale x null, bars = pool.
## n (focal-years) above each bar; free_x drops the pools a scale doesn't have.
n_lab <- distinct(summ, metric, null, level, pool, n_total)
p <- ggplot(summ, aes(pool, prop, fill = dir)) +
  geom_col(width = 0.85, position = position_stack(reverse = TRUE)) +
  geom_text(data = n_lab, aes(pool, 1.03, label = n_total),
            inherit.aes = FALSE, size = 2.5, vjust = 0, colour = "grey30") +
  facet_grid(metric ~ level + null, scales = "free_x", space = "free_x", labeller = facet_labs) +
  scale_fill_manual(values = DIR_COLS, labels = DIR_LABS, name = "Observed vs null") +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, .5, 1), expand = c(0, 0)) +
  labs(x = "Regional pool", y = "Proportion of focal-years",
       title = sprintf("Significance-flag composition  (years %s pooled)", YEAR_TAG)) +
  theme_pubr(base_size = 11, legend = "right") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 9),
        panel.spacing.x = unit(4, "pt"))
save_fig(p, sprintf("DirProps_stacked_%s", YEAR_TAG), 12, 5.5)

## 3b. SES distributions: magnitude AND sign of the departure, which the flags
## alone hide. Dotted lines at +/-1.96 are a visual guide only; the flags use the
## empirical null 2.5-97.5% quantiles.
d_ses <- subset(dat, is.finite(ses))
p <- ggplot(d_ses, aes(pool, ses, fill = null)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 3, colour = "grey60") +
  geom_boxplot(position = position_dodge(0.8), width = 0.7,
               outlier.size = 0.6, outlier.alpha = 0.5, linewidth = 0.3) +
  facet_grid(metric ~ level, scales = "free", space = "free_x", labeller = facet_labs) +
  scale_fill_manual(values = NULL_COLS, labels = NULL_LABS, name = "Null") +
  labs(x = "Regional pool", y = "SES (observed vs null)",
       title = sprintf("Departure from null  (years %s pooled)", YEAR_TAG)) +
  theme_pubr(base_size = 11, legend = "right") +
  theme(strip.background = element_blank())
save_fig(p, sprintf("SES_distributions_%s", YEAR_TAG), 9, 6)


#### 4. DO THE THREE NULLS AGREE?  (paired within focal-unit-year) ####
# All nulls are evaluated on the SAME focal units, so each pair of nulls is a
# PAIRED comparison. Join flags per focal-unit-YEAR, then aggregate over years.
NULL_PAIRS <- combn(NULLS, 2, simplify = FALSE)
paired <- do.call(rbind, lapply(NULL_PAIRS, function(np) {
  a <- dat %>% filter(null == np[1]) %>% select(level, pool, year, focal, metric, dir_a = dir)
  b <- dat %>% filter(null == np[2]) %>% select(level, pool, year, focal, metric, dir_b = dir)
  j <- inner_join(a, b, by = c("level", "pool", "year", "focal", "metric"))
  j$pair <- sprintf("%s (y) vs\n%s (x)", NULL_LABS[[np[1]]], NULL_LABS[[np[2]]])
  j
}))
paired$pair <- factor(paired$pair, levels = unique(paired$pair))

## (a) agreement table: every metric x scale x pool x null pair (pooled over years)
null_agree <- paired %>%
  group_by(pair, metric, level, pool) %>%
  summarise(n_paired  = n(),
            prop_same = mean(dir_a == dir_b),
            prop_sig_a = mean(dir_a != "neutral"),
            prop_sig_b = mean(dir_b != "neutral"),
            .groups = "drop") %>%
  mutate(pair = gsub("\n", " ", pair),
         across(starts_with("prop_"), ~ round(.x, 3))) %>%
  arrange(metric, level, pool, pair)
write.csv(null_agree, file.path(OUT_DIR, sprintf("DirPatterns_NullAgreement_%s.csv", YEAR_TAG)),
          row.names = FALSE)
print(null_agree)

## (b) confusion tiles at the headline pool per scale: rows = null pair,
## cols = scale x metric. Counts are focal-years.
conf <- paired %>%
  filter(as.character(pool) == FOCUS_POOL[as.character(level)]) %>%
  count(pair, level, metric, dir_a, dir_b, .drop = FALSE)
p <- ggplot(conf, aes(dir_b, dir_a, fill = n)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n), size = 3) +
  facet_grid(pair ~ level + metric, labeller = labeller(level = SCALE_LABS, metric = METRIC_LABS)) +
  scale_fill_gradient(low = "grey95", high = "#2C7FB8", name = "focal-years") +
  scale_x_discrete(labels = c(lower = "lower", neutral = "n.s.", higher = "higher")) +
  scale_y_discrete(labels = c(lower = "lower", neutral = "n.s.", higher = "higher")) +
  labs(x = "Flag under second null", y = "Flag under first null",
       title = sprintf("Null-model flag agreement  (plot by %s, site by %s; years %s pooled)",
                       FOCUS_POOL[["plot"]], FOCUS_POOL[["site"]], YEAR_TAG)) +
  theme_pubr(base_size = 11, legend = "right") +
  theme(strip.background = element_blank(), strip.text.y = element_text(size = 8),
        panel.grid = element_blank())
save_fig(p, sprintf("DirNullConfusion_%s", YEAR_TAG), 11, 7.5)


#### 5-6. MAPS SETUP: coordinates, year averaging, shared map builder ####
# Coordinates come from neonDivData::neon_location beetle ("bet") locations; a
# site point is the centroid of its beetle plots. Maps use the headline pool per
# scale (FOCUS_POOL) and average the two years to ONE point per focal unit:
#   obs, ses  mean over available years
#   dir       "lower"/"higher" only if the flag agreed in EVERY available year,
#             otherwise "neutral" (n.s. or mixed). Conservative by design.
library(neonDivData)
library(maps)

bet_loc <- subset(neon_location,
                  substr(location_id, nchar(location_id) - 2, nchar(location_id)) == "bet")
plot_xy <- unique(bet_loc[, c("plotID", "longitude", "latitude")])
plot_xy <- plot_xy[!is.na(plot_xy$plotID) & !is.na(plot_xy$longitude), ]
plot_xy$siteID <- sub("_.*", "", plot_xy$plotID)                            # HARV_001 -> HARV
site_xy <- aggregate(cbind(longitude, latitude) ~ siteID, plot_xy, mean)   # site = plot centroid

states_map <- map_data("state")
CONUS      <- list(xlim = c(-125, -66), ylim = c(24, 50))
DIR_SHAPES <- c(lower = 25, neutral = 21, higher = 24)   # triangle-down / circle / triangle-up
MAP_DIR_LABS <- c(lower = "lower (all years)", neutral = "n.s. or mixed",
                  higher = "higher (all years)")

# one row per null x metric x focal unit at the headline pool, years averaged,
# with coordinates attached and non-CONUS points reported
map_frame <- function(LV) {
  d <- dat %>%
    filter(level == LV, pool == FOCUS_POOL[[LV]]) %>%
    group_by(null, metric, focal) %>%
    summarise(obs     = mean(obs, na.rm = TRUE),
              ses     = if (all(is.na(ses))) NA_real_ else mean(ses, na.rm = TRUE),
              n_years = n(),
              dir     = if (all(dir == "lower")) "lower" else
                        if (all(dir == "higher")) "higher" else "neutral",
              .groups = "drop") %>%
    mutate(dir = factor(dir, levels = DIR_LEVELS))
  xy <- if (LV == "plot") plot_xy[, c("plotID", "longitude", "latitude")] else site_xy
  d  <- merge(d, xy, by.x = "focal", by.y = names(xy)[1])
  n_out <- length(unique(d$focal[d$longitude < CONUS$xlim[1] | d$longitude > CONUS$xlim[2] |
                                 d$latitude  < CONUS$ylim[1] | d$latitude  > CONUS$ylim[2]]))
  if (n_out > 0) message("  [", LV, "] ", n_out, " focal units outside CONUS (AK/HI/PR) not drawn")
  d
}

# base CONUS map + points. fill = `fillvar` (set by the caller); with use_shape,
# shape = consensus direction flag. Jitter is seeded in position_jitter so it is
# reproducible at render time.
make_map <- function(d, jit, fill_scale, title = NULL, facet = NULL,
                     use_shape = FALSE, legend = "right") {
  pos <- position_jitter(width = jit, height = jit, seed = 517)
  pts <- if (use_shape)
    geom_point(data = d, aes(longitude, latitude, fill = fillvar, shape = dir),
               position = pos, size = 2, stroke = 0.2, colour = "grey30", alpha = 0.9)
  else
    geom_point(data = d, aes(longitude, latitude, fill = fillvar), shape = 21,
               position = pos, size = 2, stroke = 0.2, colour = "grey30", alpha = 0.9)
  p <- ggplot() +
    geom_polygon(data = states_map, aes(long, lat, group = group),
                 fill = "grey95", colour = "grey75", linewidth = 0.2) +
    pts + facet + fill_scale +
    coord_quickmap(xlim = CONUS$xlim, ylim = CONUS$ylim) +
    labs(title = title, x = NULL, y = NULL) +
    theme_pubr(base_size = 11, legend = legend) +
    theme(strip.background = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(), panel.spacing = unit(6, "pt"))
  if (use_shape)
    p <- p + scale_shape_manual(values = DIR_SHAPES, labels = MAP_DIR_LABS, drop = FALSE,
                                name = "Observed vs null") +
      guides(shape = guide_legend(override.aes = list(fill = "grey50")))
  p
}


#### 5. OBSERVED-VALUE MAPS ####
# The observed value is the same under every null and pool, so one copy is
# mapped (PoolNull rows). The two metrics live on different scales, so each gets
# its own map and fill scale, placed side by side.
for (LV in names(POOLS_FOR)) {
  base <- subset(map_frame(LV), null == "PoolNull" & is.finite(obs))
  if (nrow(base) == 0) { message("[skip] no mappable obs: ", LV); next }
  maps_m <- lapply(METRICS, function(m) {
    d <- subset(base, metric == m); d$fillvar <- d$obs
    make_map(d, JITTER[[LV]], scale_fill_viridis_c(name = NULL),
             title = METRIC_LABS[[m]], legend = "bottom")
  })
  p <- wrap_plots(maps_m, nrow = 1) +
    plot_annotation(title = sprintf("%s-level observed values  (mean of %s)",
                                    SCALE_LABS[[LV]], YEAR_TAG))
  save_fig(p, sprintf("MapObs_%s_%s", LV, YEAR_TAG), 11, 4.2)
}


#### 6. SES MAPS ####
# SES is standardized, so both metrics share one diverging fill scale. Rows =
# null, cols = metric; shape = consensus flag across years.
for (LV in names(POOLS_FOR)) {
  d <- subset(map_frame(LV), is.finite(ses))
  if (nrow(d) == 0) { message("[skip] no mappable SES: ", LV); next }
  d$fillvar <- d$ses
  p <- make_map(d, JITTER[[LV]],
                scale_fill_gradient2(low = "#2C7FB8", mid = "grey92", high = "#D95F02",
                                     midpoint = 0, limits = c(-SES_CLIP, SES_CLIP),
                                     oob = scales::squish, name = "SES"),
                title = sprintf("%s-level SES  (pool = %s; mean of %s)",
                                SCALE_LABS[[LV]], FOCUS_POOL[[LV]], YEAR_TAG),
                facet = facet_grid(null ~ metric, labeller = facet_labs),
                use_shape = TRUE)
  save_fig(p, sprintf("MapSES_%s_by_%s_%s", LV, FOCUS_POOL[[LV]], YEAR_TAG), 10, 8)
}
