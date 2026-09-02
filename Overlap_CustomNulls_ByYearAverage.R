#### AVERAGE BY-YEAR OVERLAP ACROSS YEARS + YEAR-TO-YEAR STABILITY ####
# -----------------------------------------------------------------------------
# Companion to Overlap_CustomNulls_ByYear.R. That script writes one set of null
# outputs per year (<prefix>_2018_*Null.csv, <prefix>_2019_*Null.csv). Overlap is
# computed per year on purpose: richness is a median across annual rarefaction
# estimates (rarefaction assumes no within-window turnover), so overlap and
# richness only line up on a shared per-year basis. This script:
#
#   (1) averages each metric's observed value and SES across the two years, per
#       focal unit  -> the per-year-then-average basis that matches median-annual
#       rarefied richness. n_overlap_sp is now a PER-YEAR species count, so its
#       mean should sit at/below the median-annual richness estimate -- the check
#       that closes the original "more overlap species than richness" mismatch.
#   (2) reports how stable each metric is from 2018 to 2019 (a sanity check that
#       the average combines like with like, rather than papering over a swing).
#
# Writes averaged null files for the SEM and stability figures/tables. Does NOT
# overwrite the pooled or per-year outputs.
# -----------------------------------------------------------------------------

library(ggplot2)

setwd("/home/aly/Beetles/BeetleBodySizeVariation")


#### 0. SETTINGS  (match the by-year run you want to combine) ####
# Defaults are for interactive runs. When driven by run_ByYear.sh these are
# overridden by command-line args, in this order:
#   Rscript Overlap_CustomNulls_ByYear.R <LEVEL> <POOL> <YEAR>
LEVEL <- "site"     # "plot" or "site"
POOL  <- "all"      # "site" / "domain" / "all"
NULLS <- c("PoolNull", "IndividualNull")     # null types the by-year script writes
YEARS <- c(2018, 2019)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) LEVEL <- args[1]
if (length(args) >= 2) POOL  <- args[2]

metric_names <- c("overlap_norm", "overlap_unnorm", "niche_range", "sdnnd", "min_logratio")
FOCAL_COL <- if (LEVEL == "plot") "plotID" else "siteID"
PREFIX    <- paste0(LEVEL, "_by_", POOL)     # per-year files append _<YEAR>
FIG_DIR   <- "./Figures/ByYearStability"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)


#### 1. LOOP NULL TYPES: LOAD BOTH YEARS, AVERAGE, MEASURE STABILITY ####
stability_all <- list()      # collects one row per (null x metric)

for (nl in NULLS) {

  # ---- read the two per-year files ----
  y1 <- read.csv(sprintf("./Outputs/%s_%d_%s.csv", PREFIX, YEARS[1], nl))
  y2 <- read.csv(sprintf("./Outputs/%s_%d_%s.csv", PREFIX, YEARS[2], nl))

  # keep only what we combine: focal id, per-metric obs/ses/dir, lat, n_overlap_sp
  keep <- c(FOCAL_COL,
            paste0(rep(metric_names, each = 3), c("_obs", "_ses", "_dir")),
            "latitude", "n_overlap_sp")
  y1 <- y1[, intersect(keep, names(y1))]
  y2 <- y2[, intersect(keep, names(y2))]

  # full join on focal unit so single-year units are kept (flagged by n_years below)
  m <- merge(y1, y2, by = FOCAL_COL, all = TRUE, suffixes = c("_y1", "_y2"))

  # ---- build the averaged output ----
  out <- data.frame(m[FOCAL_COL], stringsAsFactors = FALSE)
  out$latitude        <- ifelse(is.na(m$latitude_y1),      m$latitude_y2,      m$latitude_y1)
  out$n_overlap_sp_y1 <- m$n_overlap_sp_y1
  out$n_overlap_sp_y2 <- m$n_overlap_sp_y2
  out$n_overlap_sp_mean <- rowMeans(cbind(m$n_overlap_sp_y1, m$n_overlap_sp_y2), na.rm = TRUE)

  for (metric in metric_names) {
    o1 <- m[[paste0(metric, "_obs_y1")]]; o2 <- m[[paste0(metric, "_obs_y2")]]
    s1 <- m[[paste0(metric, "_ses_y1")]]; s2 <- m[[paste0(metric, "_ses_y2")]]
    d1 <- m[[paste0(metric, "_dir_y1")]]; d2 <- m[[paste0(metric, "_dir_y2")]]

    # per-year values retained for transparency
    out[[paste0(metric, "_obs_y1")]] <- o1
    out[[paste0(metric, "_obs_y2")]] <- o2
    out[[paste0(metric, "_ses_y1")]] <- s1
    out[[paste0(metric, "_ses_y2")]] <- s2

    # mean across AVAILABLE years, and how many years contributed (0/1/2)
    out[[paste0(metric, "_obs_mean")]] <- rowMeans(cbind(o1, o2), na.rm = TRUE)
    out[[paste0(metric, "_ses_mean")]] <- rowMeans(cbind(s1, s2), na.rm = TRUE)
    out[[paste0(metric, "_n_years")]]  <- (!is.na(o1)) + (!is.na(o2))

    # direction agreement (TRUE only when both years present and the flags match)
    out[[paste0(metric, "_dir_agree")]] <- ifelse(is.na(d1) | is.na(d2), NA, d1 == d2)
  }

  # rowMeans over an all-NA row returns NaN -> restore NA
  for (cc in grep("_mean$", names(out), value = TRUE)) out[[cc]][is.nan(out[[cc]])] <- NA

  write.csv(out, sprintf("./Outputs/%s_ByYearAvg_%s.csv", PREFIX, nl), row.names = FALSE)
  message("wrote ", PREFIX, "_ByYearAvg_", nl, ".csv  (", nrow(out), " focal units; ",
          sum(out[[paste0(metric_names[1], "_n_years")]] == 2, na.rm = TRUE), " with both years)")

  # ---- stability per metric: 2018 vs 2019 ----
  for (metric in metric_names) {
    o1   <- out[[paste0(metric, "_obs_y1")]]; o2 <- out[[paste0(metric, "_obs_y2")]]
    both <- is.finite(o1) & is.finite(o2)
    da   <- out[[paste0(metric, "_dir_agree")]]
    stability_all[[length(stability_all) + 1]] <- data.frame(
      null           = nl,
      metric         = metric,
      n_both         = sum(both),
      pearson_obs    = if (sum(both) >= 3) cor(o1[both], o2[both]) else NA,
      spearman_obs   = if (sum(both) >= 3) cor(o1[both], o2[both], method = "spearman") else NA,
      mean_absdiff   = if (any(both))  mean(abs(o1[both] - o2[both])) else NA,
      prop_dir_agree = if (any(!is.na(da))) mean(da, na.rm = TRUE) else NA,
      stringsAsFactors = FALSE
    )
  }

  # ---- 1:1 scatter of observed overlap, 2018 vs 2019, faceted by metric ----
  long <- do.call(rbind, lapply(metric_names, function(metric)
    data.frame(metric = metric,
               y1 = out[[paste0(metric, "_obs_y1")]],
               y2 = out[[paste0(metric, "_obs_y2")]])))
  p <- ggplot(subset(long, is.finite(y1) & is.finite(y2)), aes(y1, y2)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60") +
    geom_point(alpha = 0.6) +
    facet_wrap(~ metric, scales = "free") +
    labs(x = paste0(YEARS[1], " observed"), y = paste0(YEARS[2], " observed"),
         title = paste0(nl, "  --  year-to-year stability (", LEVEL, ", pool = ", POOL, ")")) +
    theme_bw()
  ggsave(sprintf("%s/%s_%s_obs_2018v2019.png", FIG_DIR, PREFIX, nl),
         p, width = 9, height = 6, dpi = 300)
}


#### 2. STABILITY SUMMARY TABLE ####
stability <- do.call(rbind, stability_all)
num <- c("pearson_obs", "spearman_obs", "mean_absdiff", "prop_dir_agree")
stability[num] <- round(stability[num], 3)
write.csv(stability, sprintf("./Outputs/%s_ByYearStability.csv", PREFIX), row.names = FALSE)
print(stability)
