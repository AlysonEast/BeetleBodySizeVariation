#### COMPARE OVERLAP METHODS & NULL MODELS  (paired scatters + departure-from-null) ####
# -----------------------------------------------------------------------------
# Companion to Overlap_CustomNulls_ByYear.R. That script writes one file per
#   <LEVEL>_by_<POOL>_<aug|noaug>_<YEAR>_<PoolNull|IndividualNull>.csv
# holding, for each of five metrics, the observed value plus its null CI, SES,
# and a direction flag:
#   <metric>_obs   <metric>_lower  <metric>_upper  <metric>_ses  <metric>_dir
# metrics: overlap_norm, overlap_unnorm, niche_range, sdnnd, min_logratio
# plus columns: <focal id> (plotID/siteID), <pool id>, latitude, n_overlap_sp.
#
# This script does NOT recompute anything. It reads those outputs and makes
# PAIRED comparisons, holding everything else constant, to answer "does this
# methodological choice change my answer?" For each comparison it shows BOTH:
#   (1) the observed metric value on each axis (does the raw number move?), and
#   (2) the DEPARTURE FROM THE NULL MODEL as the SES (does the inference move?).
# A point is a focal unit; the dashed line is 1:1. If choices don't matter,
# points sit on the line.
#
# Five comparisons (each holds the other axes at the DEFAULT context below):
#   A  Augmented vs unaugmented         (aug   vs noaug)
#   B  Interannual                      (2018  vs 2019)
#   C  Norm vs unnorm overlap           (is abundance-weighting important?)
#   D  Pool extent
#        plot level: site vs domain, domain vs all
#        site level: domain vs all
#   E  Plot vs site                     (plot means aggregated to site vs site)
#
# Flat / stepwise: run 0-2 once (setup + load), then any comparison section on
# its own. Sweep a different context by editing the DEFAULT block and re-running.
# -----------------------------------------------------------------------------

library(ggplot2)
library(ggpubr)     # theme_pubr
library(svglite)    # svg device

setwd("/home/aly/Beetles/BeetleBodySizeVariation")


#### 0. SETTINGS ####
OUT_DIR <- "./Outputs"
FIG_DIR <- "./Figures/NullComparisons"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

METRICS <- c("overlap_norm", "overlap_unnorm", "niche_range", "sdnnd", "min_logratio")
NULLS   <- c("PoolNull", "IndividualNull")   # loop the comparisons over each null

# Context held constant in every comparison EXCEPT the axis that comparison varies.
# Edit these and re-run to sweep (e.g. LEVEL_DEF <- "plot").
LEVEL_DEF <- "site"     # "site" or "plot"
POOL_DEF  <- "all"      # a pool valid for LEVEL_DEF (site: domain/all; plot: site/domain/all)
YEAR_DEF  <- 2018       # 2018 or 2019
AUG_DEF   <- "aug"      # "aug" or "noaug"

SAVE_SVG  <- FALSE       # also write an .svg next to each .png
FIG_W     <- 11         # inches; height scales with the facet grid
DPI       <- 300

# pairings that actually exist in the output grid (a pool must sit above the focal level)
POOLS_FOR <- list(plot = c("site", "domain", "all"), site = c("domain", "all"))


#### 1. HELPERS  (each reused across >= 2 comparisons) ####

# plot -> site crosswalk, so plot-level outputs can be aggregated up to site (comparison E)
.clean <- read.csv("./Data/BodysizeCombinedClean.csv", stringsAsFactors = FALSE)
plot_to_site <- unique(.clean[, c("plotID", "siteID")])
rm(.clean)

# read one null-output file and return it LONG: one row per focal unit x metric.
# Missing files are skipped with a message (partial grids don't abort the script).
read_null_long <- function(level, pool, aug, year, null) {
  f <- file.path(OUT_DIR, sprintf("%s_by_%s_%s_%s_%s.csv", level, pool, aug, year, null))
  if (!file.exists(f)) { message("  [skip] missing ", basename(f)); return(NULL) }
  d <- read.csv(f, stringsAsFactors = FALSE)
  focal_col <- if (level == "plot") "plotID" else "siteID"
  long <- do.call(rbind, lapply(METRICS, function(m) data.frame(
    focal        = d[[focal_col]],
    metric       = m,
    obs          = d[[paste0(m, "_obs")]],
    ses          = d[[paste0(m, "_ses")]],
    dir          = d[[paste0(m, "_dir")]],
    latitude     = d$latitude,
    n_overlap_sp = d$n_overlap_sp,
    stringsAsFactors = FALSE)))
  long$siteID <- if (level == "plot")
    plot_to_site$siteID[match(long$focal, plot_to_site$plotID)] else long$focal
  long$level <- level; long$pool <- pool; long$aug <- aug
  long$year  <- as.character(year); long$null <- null
  long
}

# build a paired wide frame: split `df` on `split_col` into values a/b and join on
# `key` + metric, so obs/ses/dir land as *_a / *_b columns for a paired scatter.
pair_wide <- function(df, key, split_col, a, b) {
  keep <- c(key, "metric", "obs", "ses", "dir")
  da <- df[df[[split_col]] == a, keep]
  db <- df[df[[split_col]] == b, keep]
  merge(da, db, by = c(key, "metric"), suffixes = c("_a", "_b"))
}

# draw a paired scatter. `wide` has metric, obs_a/obs_b, ses_a/ses_b. Rows = the two
# quantities (observed value; departure-from-null SES), cols = metric (unless
# metric_facet = FALSE, for the single metric-pair comparison C). Writes png (+svg).
plot_paired <- function(wide, lab_a, lab_b, title, file, metric_facet = TRUE) {
  if (is.null(wide) || nrow(wide) == 0) { message("  [skip] no paired rows for: ", title); return(invisible()) }
  long <- rbind(
    data.frame(metric = wide$metric, quantity = "Observed value",
               x = wide$obs_a, y = wide$obs_b, stringsAsFactors = FALSE),
    data.frame(metric = wide$metric, quantity = "Departure from null (SES)",
               x = wide$ses_a, y = wide$ses_b, stringsAsFactors = FALSE))
  long <- long[is.finite(long$x) & is.finite(long$y), ]
  if (nrow(long) == 0) { message("  [skip] all-NA paired rows for: ", title); return(invisible()) }
  long$quantity <- factor(long$quantity, levels = c("Observed value", "Departure from null (SES)"))
  if (metric_facet) long$metric <- factor(long$metric, levels = METRICS)

  # per-facet Pearson r + n, placed top-left of each facet
  grp <- if (metric_facet) list(long$metric, long$quantity) else list(long$quantity)
  stats <- do.call(rbind, lapply(split(long, grp, drop = TRUE), function(s) data.frame(
    metric = s$metric[1], quantity = s$quantity[1],
    r = if (nrow(s) >= 3) cor(s$x, s$y) else NA_real_, n = nrow(s),
    stringsAsFactors = FALSE)))
  stats$lab <- ifelse(is.na(stats$r), sprintf("n=%d", stats$n),
                      sprintf("r=%.2f  n=%d", stats$r, stats$n))

  p <- ggplot(long, aes(x, y)) +
    geom_hline(yintercept = 0, colour = "grey88") +
    geom_vline(xintercept = 0, colour = "grey88") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_text(data = stats, aes(x = -Inf, y = Inf, label = lab),
              inherit.aes = FALSE, hjust = -0.08, vjust = 1.4, size = 3, colour = "grey25") +
    labs(x = lab_a, y = lab_b, title = title) +
    theme_pubr(base_size = 11) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 9, lineheight = 0.9),
          plot.title = element_text(size = 11))
  # facet_wrap (not facet_grid) so EVERY panel gets its own x and y scale. Under
  # facet_grid, scales = "free" still shares x down each column, forcing the
  # Observed-value and SES panels of a metric onto one x-axis; the wide SES range
  # then squashes the [0,1]-ish overlap values. Independent panels keep both
  # rows readable. ncol = #metrics preserves the quantity-row x metric-column grid.
  facet <- if (metric_facet)
    facet_wrap(vars(quantity, metric), scales = "free", ncol = length(METRICS),
               labeller = labeller(.multi_line = TRUE))
  else
    facet_wrap(~ quantity, scales = "free", nrow = 1)
  p <- p + facet

  h <- if (metric_facet) 6 else 3.6
  w <- if (metric_facet) FIG_W else 7
  ggsave(file, p, width = w, height = h, dpi = DPI)
  if (SAVE_SVG) ggsave(sub("\\.png$", ".svg", file), p, width = w, height = h)
  message("  wrote ", basename(file))
  invisible(p)
}

# summarise a paired wide frame: agreement of the DEPARTURE-FROM-NULL inference
# across the two conditions, per metric. Collected into one table at the end.
summarise_pair <- function(wide, comparison, context) {
  if (is.null(wide) || nrow(wide) == 0) return(NULL)
  do.call(rbind, lapply(split(wide, wide$metric), function(s) {
    both <- is.finite(s$ses_a) & is.finite(s$ses_b)
    dboth <- !is.na(s$dir_a) & !is.na(s$dir_b)
    data.frame(
      comparison = comparison, context = context, metric = s$metric[1],
      n_paired          = sum(both),
      pearson_ses       = if (sum(both)  >= 3) cor(s$ses_a[both], s$ses_b[both]) else NA,
      spearman_ses      = if (sum(both)  >= 3) cor(s$ses_a[both], s$ses_b[both], method = "spearman") else NA,
      mean_absdiff_ses  = if (any(both))  mean(abs(s$ses_a[both] - s$ses_b[both])) else NA,
      prop_dir_agree    = if (any(dboth)) mean(s$dir_a[dboth] == s$dir_b[dboth]) else NA,
      stringsAsFactors = FALSE)
  }))
}


#### 2. LOAD ALL AVAILABLE OUTPUTS INTO ONE LONG TABLE ####
grid <- expand.grid(level = c("plot", "site"), aug = c("aug", "noaug"),
                    year = c(2018, 2019), null = NULLS, stringsAsFactors = FALSE)
dat_list <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  for (pool in POOLS_FOR[[g$level]])
    dat_list[[length(dat_list) + 1]] <-
      read_null_long(g$level, pool, g$aug, g$year, g$null)
}
dat <- do.call(rbind, dat_list)
message("loaded ", nrow(dat), " focal x metric rows across ",
        length(unique(paste(dat$level, dat$pool, dat$aug, dat$year, dat$null))), " files")

agree <- list()   # collector for the summary table (filled by each section)


#### 3. COMPARISON A: AUGMENTED vs UNAUGMENTED ####
# Holds LEVEL_DEF / POOL_DEF / YEAR_DEF; varies aug. Does padding sparse cells to
# n = 20 move the observed overlap and/or the inference?
for (nl in NULLS) {
  ctx <- sprintf("%s, pool=%s, %d", LEVEL_DEF, POOL_DEF, YEAR_DEF)
  df  <- subset(dat, level == LEVEL_DEF & pool == POOL_DEF & year == as.character(YEAR_DEF) & null == nl)
  w   <- pair_wide(df, key = "focal", split_col = "aug", a = "aug", b = "noaug")
  plot_paired(w, "Augmented", "Unaugmented",
              sprintf("Augmented vs unaugmented  |  %s  |  %s", ctx, nl),
              file.path(FIG_DIR, sprintf("A_aug_v_noaug_%s_by_%s_%d_%s.png",
                                         LEVEL_DEF, POOL_DEF, YEAR_DEF, nl)))
  agree[[length(agree) + 1]] <- summarise_pair(w, "aug_vs_noaug", paste(ctx, nl, sep = " | "))
}


#### 4. COMPARISON B: INTERANNUAL (2018 vs 2019) ####
# Holds LEVEL_DEF / POOL_DEF / AUG_DEF; varies year.
for (nl in NULLS) {
  ctx <- sprintf("%s, pool=%s, %s", LEVEL_DEF, POOL_DEF, AUG_DEF)
  df  <- subset(dat, level == LEVEL_DEF & pool == POOL_DEF & aug == AUG_DEF & null == nl)
  w   <- pair_wide(df, key = "focal", split_col = "year", a = "2018", b = "2019")
  plot_paired(w, "2018", "2019",
              sprintf("Interannual 2018 vs 2019  |  %s  |  %s", ctx, nl),
              file.path(FIG_DIR, sprintf("B_2018_v_2019_%s_by_%s_%s_%s.png",
                                         LEVEL_DEF, POOL_DEF, AUG_DEF, nl)))
  agree[[length(agree) + 1]] <- summarise_pair(w, "2018_vs_2019", paste(ctx, nl, sep = " | "))
}


#### 5. COMPARISON C: NORM vs UNNORM OVERLAP (is abundance important?) ####
# Both metrics live in the SAME file, so this pairs metric-to-metric on the focal
# unit (not a config axis). Single metric-pair -> facet by quantity only.
for (nl in NULLS) {
  ctx <- sprintf("%s, pool=%s, %d, %s", LEVEL_DEF, POOL_DEF, YEAR_DEF, AUG_DEF)
  df  <- subset(dat, level == LEVEL_DEF & pool == POOL_DEF & year == as.character(YEAR_DEF) &
                     aug == AUG_DEF & null == nl & metric %in% c("overlap_norm", "overlap_unnorm"))
  wa <- df[df$metric == "overlap_norm",   c("focal", "obs", "ses", "dir")]
  wb <- df[df$metric == "overlap_unnorm", c("focal", "obs", "ses", "dir")]
  w  <- merge(wa, wb, by = "focal", suffixes = c("_a", "_b"))
  w$metric <- "overlap"   # single pseudo-metric so plot_paired's columns line up
  plot_paired(w, "overlap_norm (unweighted)", "overlap_unnorm (abundance-weighted)",
              sprintf("Norm vs unnorm overlap  |  %s  |  %s", ctx, nl),
              file.path(FIG_DIR, sprintf("C_norm_v_unnorm_%s_by_%s_%d_%s_%s.png",
                                         LEVEL_DEF, POOL_DEF, YEAR_DEF, AUG_DEF, nl)),
              metric_facet = FALSE)
  s <- summarise_pair(w, "norm_vs_unnorm", paste(ctx, nl, sep = " | "))
  if (!is.null(s)) s$metric <- "overlap_norm_vs_unnorm"
  agree[[length(agree) + 1]] <- s
}


#### 6. COMPARISON D: POOL EXTENT ####
# Same focal units evaluated against different regional pools. Does widening the
# pool (site -> domain -> all) change the departure from null?
# 6a. Plot level: site vs domain, then domain vs all.
for (nl in NULLS) {
  for (pp in list(c("site", "domain"), c("domain", "all"))) {
    ctx <- sprintf("plot, %s, %d", AUG_DEF, YEAR_DEF)
    df  <- subset(dat, level == "plot" & aug == AUG_DEF & year == as.character(YEAR_DEF) &
                       null == nl & pool %in% pp)
    w   <- pair_wide(df, key = "focal", split_col = "pool", a = pp[1], b = pp[2])
    plot_paired(w, paste0("pool = ", pp[1]), paste0("pool = ", pp[2]),
                sprintf("Plot pool extent: %s vs %s  |  %s  |  %s", pp[1], pp[2], ctx, nl),
                file.path(FIG_DIR, sprintf("D_pool_plot_%s_v_%s_%s_%d_%s.png",
                                           pp[1], pp[2], AUG_DEF, YEAR_DEF, nl)))
    agree[[length(agree) + 1]] <- summarise_pair(
      w, sprintf("pool_plot_%s_vs_%s", pp[1], pp[2]), paste(ctx, nl, sep = " | "))
  }
}
# 6b. Site level: domain vs all.
for (nl in NULLS) {
  ctx <- sprintf("site, %s, %d", AUG_DEF, YEAR_DEF)
  df  <- subset(dat, level == "site" & aug == AUG_DEF & year == as.character(YEAR_DEF) &
                     null == nl & pool %in% c("domain", "all"))
  w   <- pair_wide(df, key = "focal", split_col = "pool", a = "domain", b = "all")
  plot_paired(w, "pool = domain", "pool = all",
              sprintf("Site pool extent: domain vs all  |  %s  |  %s", ctx, nl),
              file.path(FIG_DIR, sprintf("D_pool_site_domain_v_all_%s_%d_%s.png",
                                         AUG_DEF, YEAR_DEF, nl)))
  agree[[length(agree) + 1]] <- summarise_pair(w, "pool_site_domain_vs_all", paste(ctx, nl, sep = " | "))
}


#### 7. COMPARISON E: PLOT vs SITE ####
# Focal units differ (plotID vs siteID), so aggregate plot results to the site
# (mean over the plots in a site) and pair against the site-level result for the
# SAME pool. Matched pools available at both levels: domain and all.
for (nl in NULLS) {
  for (pl in c("domain", "all")) {
    ctx   <- sprintf("pool=%s, %s, %d", pl, AUG_DEF, YEAR_DEF)
    plotd <- subset(dat, level == "plot" & pool == pl & aug == AUG_DEF &
                         year == as.character(YEAR_DEF) & null == nl)
    sited <- subset(dat, level == "site" & pool == pl & aug == AUG_DEF &
                         year == as.character(YEAR_DEF) & null == nl)
    if (nrow(plotd) == 0 || nrow(sited) == 0) {
      message("  [skip] plot-vs-site pool=", pl, " ", nl, " (missing a level)"); next }
    # aggregate obs and ses to site x metric (na.pass so all-NA cells stay NA, not dropped)
    ao <- aggregate(obs ~ siteID + metric, plotd, mean, na.rm = TRUE, na.action = na.pass)
    as_ <- aggregate(ses ~ siteID + metric, plotd, mean, na.rm = TRUE, na.action = na.pass)
    agg <- merge(ao, as_, by = c("siteID", "metric"))
    agg$dir_a <- NA_character_   # aggregated direction is undefined; agreement skipped for E
    w <- merge(agg, sited[, c("focal", "metric", "obs", "ses", "dir")],
               by.x = c("siteID", "metric"), by.y = c("focal", "metric"),
               suffixes = c("_a", "_b"))
    names(w)[names(w) == "dir"] <- "dir_b"
    plot_paired(w, "Plot (mean over plots in site)", "Site",
                sprintf("Plot vs site  |  %s  |  %s", ctx, nl),
                file.path(FIG_DIR, sprintf("E_plot_v_site_pool_%s_%s_%d_%s.png",
                                           pl, AUG_DEF, YEAR_DEF, nl)))
    agree[[length(agree) + 1]] <- summarise_pair(w, "plot_vs_site", paste(ctx, nl, sep = " | "))
  }
}


#### 8. AGREEMENT SUMMARY TABLE ####
# One row per (comparison x context x metric): how well do the two conditions
# agree on the departure-from-null inference? pearson/spearman on SES, mean
# absolute SES difference, and the fraction of focal units whose direction flag
# (lower/higher/neutral) matches. High agreement = the choice doesn't matter.
summary_tab <- do.call(rbind, agree)
num <- c("pearson_ses", "spearman_ses", "mean_absdiff_ses", "prop_dir_agree")
summary_tab[num] <- lapply(summary_tab[num], function(x) round(x, 3))
write.csv(summary_tab, file.path(OUT_DIR, "NullComparison_Agreement.csv"), row.names = FALSE)
message("wrote NullComparison_Agreement.csv  (", nrow(summary_tab), " rows)")
print(summary_tab)
