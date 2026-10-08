#### COMPARE NULL MODELS, YEARS, POOLS AND SCALES  (niche_range + overlap_depth) ####
# -----------------------------------------------------------------------------
# Companion to 6a_CommunityShape_CustomNulls_ByYear.R, which writes one file per
#   <LEVEL>_by_<POOL>_<YEAR>_<PoolNull|IndividualNull|SwapMeansNull>.csv
# holding, per metric, <metric>_obs / _lower / _upper / _ses / _dir, plus the
# focal id (plotID/siteID), the pool id, latitude, and n_comm_sp.
#
# 6a writes seven metrics; this script compares only the two that carry the
# analysis:
#   niche_range    width of occupied trait space (2.5-97.5% span)
#   overlap_depth  mean peak-scaled species co-occupancy over the occupied range
# There is no augmentation dimension any more.
#
# Nothing is recomputed. Each comparison pairs focal units across ONE choice,
# holding everything else at the DEFAULT context, and asks "does this choice
# change my answer?" A point is a focal unit; the dashed line is 1:1.
#
# What is worth plotting depends on the comparison. The observed value is
# computed from the focal community alone, so it is IDENTICAL across nulls and
# across pools (section 2 checks this). Those comparisons therefore show SES
# only; the year and scale comparisons show the observed value AND the SES.
#
#   A  Null vs null   Pool vs Individual, Pool vs SwapMeans,       (SES)
#                     Individual vs SwapMeans
#   B  Interannual    2018 vs 2019                                 (obs + SES)
#   C  Pool extent    plot: site vs domain, domain vs all          (SES)
#                     site: domain vs all
#                     (SwapMeans skipped: it needs no pool)
#   D  Plot vs site   plot values averaged to site vs site         (obs + SES)
#
# Every figure has one column per metric and one row per quantity (observed
# value, or SES under a given null / null pair).
#
# Flat / stepwise: run 0-2 once (setup + load), then any comparison on its own.
# Sweep a different context by editing the DEFAULT block and re-running.
# -----------------------------------------------------------------------------

library(ggplot2)
library(ggpubr)     # theme_pubr
library(svglite)    # svg device

setwd("/home/aly/Beetles/BeetleBodySizeVariation")


#### 0. SETTINGS ####
OUT_DIR <- "./Outputs"
FIG_DIR <- "./Figures/NullComparisons"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

METRICS     <- c("niche_range", "overlap_depth")
METRIC_LABS <- c(niche_range = "Niche range", overlap_depth = "Overlap depth")
NULLS       <- c("PoolNull", "IndividualNull", "SwapMeansNull")
NULL_LABS   <- c(PoolNull = "Pool null", IndividualNull = "Individual null",
                 SwapMeansNull = "Swap-means null")
YEARS       <- c(2018, 2019)

# Context held constant in every comparison EXCEPT the axis that comparison varies.
LEVEL_DEF <- "site"     # "site" or "plot"
POOL_DEF  <- "all"      # a pool valid for LEVEL_DEF (site: domain/all; plot: site/domain/all)
YEAR_DEF  <- 2018       # 2018 or 2019

SAVE_SVG <- FALSE       # also write an .svg next to each .png
DPI      <- 300

# pairings that exist in the output grid (a pool must sit above the focal level)
POOLS_FOR <- list(plot = c("site", "domain", "all"), site = c("domain", "all"))

OBS_ROW <- "Observed value"


#### 1. HELPERS  (each reused across >= 2 comparisons) ####

# plot -> site crosswalk, so plot-level outputs can be aggregated to site (D)
.clean <- read.csv("./Data/BodysizeCombinedClean.csv", stringsAsFactors = FALSE)
plot_to_site <- unique(.clean[, c("plotID", "siteID")])
rm(.clean)

# read one null-output file and return it LONG: one row per focal unit x metric.
# Missing files are skipped with a message (partial grids don't abort the script).
read_null_long <- function(level, pool, year, null) {
  f <- file.path(OUT_DIR, sprintf("%s_by_%s_%s_%s.csv", level, pool, year, null))
  if (!file.exists(f)) { message("  [skip] missing ", basename(f)); return(NULL) }
  d <- read.csv(f, stringsAsFactors = FALSE)
  focal_col <- if (level == "plot") "plotID" else "siteID"
  long <- do.call(rbind, lapply(METRICS, function(m) data.frame(
    focal     = d[[focal_col]],
    metric    = m,
    obs       = d[[paste0(m, "_obs")]],
    ses       = d[[paste0(m, "_ses")]],
    dir       = d[[paste0(m, "_dir")]],
    latitude  = d$latitude,
    n_comm_sp = d$n_comm_sp,
    stringsAsFactors = FALSE)))
  long$siteID <- if (level == "plot")
    plot_to_site$siteID[match(long$focal, plot_to_site$plotID)] else long$focal
  long$level <- level; long$pool <- pool
  long$year  <- as.character(year); long$null <- null
  long
}

# paired wide frame: split `df` on `split_col` into values a/b and join on
# `key` + metric, so obs/ses/dir land as *_a / *_b columns.
pair_wide <- function(df, key, split_col, a, b) {
  keep <- c(key, "metric", "obs", "ses", "dir")
  da <- df[df[[split_col]] == a, keep]
  db <- df[df[[split_col]] == b, keep]
  merge(da, db, by = c(key, "metric"), suffixes = c("_a", "_b"))
}

# one block of plotting rows from a paired wide frame: quantity `qty` ("obs" or
# "ses") on x (condition a) and y (condition b), labelled with facet row `row`.
stack_pair <- function(wide, qty, row) {
  if (is.null(wide) || nrow(wide) == 0) return(NULL)
  data.frame(metric = wide$metric, row = row,
             x = wide[[paste0(qty, "_a")]], y = wide[[paste0(qty, "_b")]],
             stringsAsFactors = FALSE)
}

# paired scatter: rows = the stacked blocks (in the order given), cols = metric.
# facet_wrap with free scales so every panel gets its own x and y (observed
# values and SES live on very different ranges). Writes png (+svg).
plot_paired <- function(long, lab_a, lab_b, title, file) {
  if (!is.null(long)) long <- long[is.finite(long$x) & is.finite(long$y), ]
  if (is.null(long) || nrow(long) == 0) { message("  [skip] no paired rows for: ", title); return(invisible()) }
  long$row    <- factor(long$row, levels = unique(long$row))
  long$metric <- factor(long$metric, levels = METRICS, labels = METRIC_LABS[METRICS])

  # per-panel Pearson r + n, placed top-left
  stats <- do.call(rbind, lapply(split(long, list(long$row, long$metric), drop = TRUE),
    function(s) data.frame(row = s$row[1], metric = s$metric[1],
                           r = if (nrow(s) >= 3) cor(s$x, s$y) else NA_real_, n = nrow(s))))
  stats$lab <- ifelse(is.na(stats$r), sprintf("n=%d", stats$n),
                      sprintf("r=%.2f  n=%d", stats$r, stats$n))

  p <- ggplot(long, aes(x, y)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_text(data = stats, aes(x = -Inf, y = Inf, label = lab),
              inherit.aes = FALSE, hjust = -0.08, vjust = 1.4, size = 3, colour = "grey25") +
    facet_wrap(vars(row, metric), scales = "free", ncol = length(METRICS)) +
    labs(x = lab_a, y = lab_b, title = title) +
    theme_pubr(base_size = 11) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 9, lineheight = 0.9),
          plot.title = element_text(size = 11))

  h <- 0.8 + 2.6 * nlevels(long$row)
  ggsave(file, p, width = 7, height = h, dpi = DPI)
  if (SAVE_SVG) ggsave(sub("\\.png$", ".svg", file), p, width = 7, height = h)
  message("  wrote ", basename(file))
  invisible(p)
}

# agreement of the departure-from-null inference across the two conditions,
# per metric. Collected into one table at the end (section 7).
summarise_pair <- function(wide, comparison, context) {
  if (is.null(wide) || nrow(wide) == 0) return(NULL)
  do.call(rbind, lapply(split(wide, wide$metric), function(s) {
    both  <- is.finite(s$ses_a) & is.finite(s$ses_b)
    dboth <- !is.na(s$dir_a) & !is.na(s$dir_b)
    data.frame(
      comparison = comparison, context = context, metric = s$metric[1],
      n_paired         = sum(both),
      pearson_ses      = if (sum(both) >= 3) cor(s$ses_a[both], s$ses_b[both]) else NA,
      spearman_ses     = if (sum(both) >= 3) cor(s$ses_a[both], s$ses_b[both], method = "spearman") else NA,
      mean_absdiff_ses = if (any(both))  mean(abs(s$ses_a[both] - s$ses_b[both])) else NA,
      prop_dir_agree   = if (any(dboth)) mean(s$dir_a[dboth] == s$dir_b[dboth]) else NA,
      stringsAsFactors = FALSE)
  }))
}


#### 2. LOAD ALL AVAILABLE OUTPUTS INTO ONE LONG TABLE ####
dat_list <- list()
for (level in names(POOLS_FOR)) for (pool in POOLS_FOR[[level]])
  for (yr in YEARS) for (nl in NULLS)
    dat_list[[length(dat_list) + 1]] <- read_null_long(level, pool, yr, nl)
dat <- do.call(rbind, dat_list)
message("loaded ", nrow(dat), " focal x metric rows across ",
        length(unique(paste(dat$level, dat$pool, dat$year, dat$null))), " files")

# Sanity check: the observed value should not depend on the null. Compare each
# null's obs against PoolNull's for the same level / pool / year / focal / metric.
ref <- dat[dat$null == "PoolNull", c("level", "pool", "year", "focal", "metric", "obs")]
chk <- merge(dat[dat$null != "PoolNull", c("level", "pool", "year", "focal", "metric", "obs")],
             ref, by = c("level", "pool", "year", "focal", "metric"), suffixes = c("", "_ref"))
if (nrow(chk)) {
  max_diff <- max(abs(chk$obs - chk$obs_ref), na.rm = TRUE)
  message("max |obs difference| across nulls: ", signif(max_diff, 3),
          if (max_diff > 1e-8) "  <- CHECK: obs should be null-invariant" else "  (null-invariant, as expected)")
}
rm(ref, chk)

agree <- list()   # collector for the summary table (filled by each section)


#### 3. COMPARISON A: NULL vs NULL (SES only) ####
# Same focal units, same observed value, different reference distribution. Do
# the three nulls rank communities the same way? Held at the DEFAULT context.
NULL_PAIRS <- combn(NULLS, 2, simplify = FALSE)
ctx  <- sprintf("%s, pool=%s, %d", LEVEL_DEF, POOL_DEF, YEAR_DEF)
df   <- subset(dat, level == LEVEL_DEF & pool == POOL_DEF & year == as.character(YEAR_DEF))
rows <- list()
for (np in NULL_PAIRS) {
  w <- pair_wide(df, key = "focal", split_col = "null", a = np[1], b = np[2])
  rows[[length(rows) + 1]] <- stack_pair(w, "ses",
    sprintf("%s (x) vs %s (y)", NULL_LABS[[np[1]]], NULL_LABS[[np[2]]]))
  agree[[length(agree) + 1]] <- summarise_pair(
    w, sprintf("%s_vs_%s", np[1], np[2]), ctx)
}
plot_paired(do.call(rbind, rows), "SES under first null", "SES under second null",
            sprintf("Null vs null SES  |  %s", ctx),
            file.path(FIG_DIR, sprintf("A_null_v_null_%s_by_%s_%d.png",
                                       LEVEL_DEF, POOL_DEF, YEAR_DEF)))


#### 4. COMPARISON B: INTERANNUAL (2018 vs 2019) ####
# Holds LEVEL_DEF / POOL_DEF; varies year. Top row = observed value (taken from
# PoolNull; identical in every null), then one SES row per null.
ctx  <- sprintf("%s, pool=%s", LEVEL_DEF, POOL_DEF)
df   <- subset(dat, level == LEVEL_DEF & pool == POOL_DEF)
rows <- list()
for (nl in NULLS) {
  w <- pair_wide(df[df$null == nl, ], key = "focal", split_col = "year", a = "2018", b = "2019")
  if (nl == "PoolNull") rows[[length(rows) + 1]] <- stack_pair(w, "obs", OBS_ROW)
  rows[[length(rows) + 1]] <- stack_pair(w, "ses", paste("SES:", NULL_LABS[[nl]]))
  agree[[length(agree) + 1]] <- summarise_pair(w, "2018_vs_2019", paste(ctx, nl, sep = " | "))
}
plot_paired(do.call(rbind, rows), "2018", "2019",
            sprintf("Interannual 2018 vs 2019  |  %s", ctx),
            file.path(FIG_DIR, sprintf("B_2018_v_2019_%s_by_%s.png", LEVEL_DEF, POOL_DEF)))


#### 5. COMPARISON C: POOL EXTENT (SES only) ####
# Same focal units evaluated against a narrower vs wider regional pool. The
# observed value cannot move, so only SES is shown. SwapMeans is skipped: it
# draws from no pool, so any difference is only the MIN_POOL_UNITS screen.
POOL_PAIRS <- list(plot = list(c("site", "domain"), c("domain", "all")),
                   site = list(c("domain", "all")))
for (LV in names(POOL_PAIRS)) for (pp in POOL_PAIRS[[LV]]) {
  ctx  <- sprintf("%s, %d", LV, YEAR_DEF)
  df   <- subset(dat, level == LV & year == as.character(YEAR_DEF) & pool %in% pp)
  rows <- list()
  for (nl in c("PoolNull", "IndividualNull")) {
    w <- pair_wide(df[df$null == nl, ], key = "focal", split_col = "pool", a = pp[1], b = pp[2])
    rows[[length(rows) + 1]] <- stack_pair(w, "ses", paste("SES:", NULL_LABS[[nl]]))
    agree[[length(agree) + 1]] <- summarise_pair(
      w, sprintf("pool_%s_%s_vs_%s", LV, pp[1], pp[2]), paste(ctx, nl, sep = " | "))
  }
  plot_paired(do.call(rbind, rows), paste("pool =", pp[1]), paste("pool =", pp[2]),
              sprintf("%s pool extent: %s vs %s  |  %d", tools::toTitleCase(LV), pp[1], pp[2], YEAR_DEF),
              file.path(FIG_DIR, sprintf("C_pool_%s_%s_v_%s_%d.png", LV, pp[1], pp[2], YEAR_DEF)))
}


#### 6. COMPARISON D: PLOT vs SITE ####
# Focal units differ, so plot results are averaged to their site and paired
# against the site-level result for the SAME pool (domain and all exist at both
# levels). The aggregated direction flag is undefined, so dir agreement is NA.
# Note niche_range at plot level is bounded by the site range, so points are
# expected below 1:1; read the ranking (r), not the offset.
for (pl in c("domain", "all")) {
  ctx  <- sprintf("pool=%s, %d", pl, YEAR_DEF)
  rows <- list()
  for (nl in NULLS) {
    plotd <- subset(dat, level == "plot" & pool == pl & year == as.character(YEAR_DEF) & null == nl)
    sited <- subset(dat, level == "site" & pool == pl & year == as.character(YEAR_DEF) & null == nl)
    if (nrow(plotd) == 0 || nrow(sited) == 0) {
      message("  [skip] plot-vs-site pool=", pl, " ", nl, " (missing a level)"); next }
    # mean over plots in a site (na.pass so all-NA cells stay NA rather than dropping)
    agg <- aggregate(cbind(obs, ses) ~ siteID + metric, plotd,
                     function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE),
                     na.action = na.pass)
    w <- merge(agg, sited[, c("focal", "metric", "obs", "ses", "dir")],
               by.x = c("siteID", "metric"), by.y = c("focal", "metric"),
               suffixes = c("_a", "_b"))
    names(w)[names(w) == "dir"] <- "dir_b"
    w$dir_a <- NA_character_
    if (nl == "PoolNull") rows[[length(rows) + 1]] <- stack_pair(w, "obs", OBS_ROW)
    rows[[length(rows) + 1]] <- stack_pair(w, "ses", paste("SES:", NULL_LABS[[nl]]))
    agree[[length(agree) + 1]] <- summarise_pair(w, "plot_vs_site", paste(ctx, nl, sep = " | "))
  }
  plot_paired(do.call(rbind, rows), "Plot (mean over plots in site)", "Site",
              sprintf("Plot vs site  |  %s", ctx),
              file.path(FIG_DIR, sprintf("D_plot_v_site_pool_%s_%d.png", pl, YEAR_DEF)))
}


#### 7. AGREEMENT SUMMARY TABLE ####
# One row per (comparison x context x metric): how well do the two conditions
# agree on the departure-from-null inference? Pearson / Spearman on SES, mean
# absolute SES difference, and the fraction of focal units whose direction flag
# matches. High agreement = the choice doesn't change the answer.
summary_tab <- do.call(rbind, agree)
num <- c("pearson_ses", "spearman_ses", "mean_absdiff_ses", "prop_dir_agree")
summary_tab[num] <- lapply(summary_tab[num], function(x) round(x, 3))
write.csv(summary_tab, file.path(OUT_DIR, "NullComparison_Agreement_CommunityShape.csv"),
          row.names = FALSE)
message("wrote NullComparison_Agreement_CommunityShape.csv  (", nrow(summary_tab), " rows)")
print(summary_tab)
