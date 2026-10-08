################################################################################
# 2b_BodySizeSpatialTests.R
#
# Tests whether intraspecific elytra length differs (1) among sites and
# (2) among plots within each site, for every species with enough data, then
# draws the nested site/plot density figure (the SpeciesNestedWrappedDensityPlots
# layout from 2_BodySizeQuantification.R) with compact letter displays.
#
# Site test: lmer(length ~ siteID + (1 | plotID)). Plots are the unit of
#   replication for a site difference, so beetles are not treated as
#   independent draws from a site. Omnibus test is a likelihood ratio test
#   against the intercept-only model; pairwise contrasts are Tukey-adjusted
#   emmeans.
# Plot test: lm(length ~ plotID), fit separately within each site. Omnibus F
#   test; pairwise contrasts are Tukey-adjusted emmeans.
# Letters: groups sharing a letter do not differ at ALPHA.
# Omnibus p-values are BH-corrected across species within each test type, and
#   letters are only drawn for tests that survive the correction (unless
#   SHOW_NS_LETTERS is TRUE).
################################################################################

library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(lme4)
library(emmeans)
# multcomp and multcompView must be installed for multcomp::cld(). Not loaded
# with library() because multcomp attaches MASS, which masks dplyr::select.

#### 0. User defined variables ####
MIN_N_SITE      <- 20     # min individuals for a site to enter the site test
MIN_N_PLOT      <- 20     # min individuals for a plot to enter the plot test
ALPHA           <- 0.05
SHOW_NS_LETTERS <- FALSE  # draw letters even when the omnibus test is n.s. after FDR
DOMAINS         <- NULL   # e.g. "D07" to restrict; NULL = all domains

emm_options(lmer.df = "satterthwaite")  # needs lmerTest installed (not loaded)

setwd("/home/aly/Beetles/BeetleBodySizeVariation")
fig_dir <- "./Figures/BodySizeQuantification/SpatialTests/"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
param_tag <- paste0("minSite", MIN_N_SITE, "_minPlot", MIN_N_PLOT,
                    if (!is.null(DOMAINS)) paste0("_", paste(DOMAINS, collapse = "-")))

#### 1. Load data (same exclusions as 2_BodySizeQuantification.R) ####
all_elytra <- read.csv("./Data/BodysizeCombinedClean.csv")
all_elytra <- subset(all_elytra, individualID != "NEON.BET.D10.016322")
all_elytra <- subset(all_elytra, individualID != "NEON.BET.D08.002280")
all_elytra <- subset(all_elytra, imageID != "MLBS_009.S.20180522.jpg")
all_elytra <- subset(all_elytra, imageID != "MLBS_009.E.20180522.CARABIDS.01.jpg")
all_elytra <- subset(all_elytra, individualID != "NEON.BET.D07.003475")

all_elytra <- all_elytra %>%
  filter(!grepl("sp\\.", scientificName_Species),
         !is.na(scientificName_Species),
         scientificName_Species != "",
         !is.na(cm_elytra_max_length))

if (!is.null(DOMAINS)) all_elytra <- all_elytra %>% filter(domainID %in% DOMAINS)

all_elytra<-subset(all_elytra, yearCollected==2018)

#### 2. Identify species with enough data ####
site_n <- all_elytra %>%
  count(scientificName_Species, domainID, siteID, name = "n")
plot_n <- all_elytra %>%
  count(scientificName_Species, domainID, siteID, plotID, name = "n")

sites_ok <- site_n %>% filter(n >= MIN_N_SITE)
plots_ok <- plot_n %>% filter(n >= MIN_N_PLOT)

# Site test needs >= 2 qualifying sites
site_test_spp <- sites_ok %>%
  count(scientificName_Species, name = "n_sites") %>%
  filter(n_sites >= 2)

# Plot test needs >= 2 qualifying plots within a site
plot_test_sites <- plots_ok %>%
  count(scientificName_Species, siteID, name = "n_plots") %>%
  filter(n_plots >= 2)

species_list <- sort(union(site_test_spp$scientificName_Species,
                           plot_test_sites$scientificName_Species))
length(species_list)
species_list

#### 3. Compact letter display helper (used by both the site and plot tests) ####
get_letters <- function(model, group) {
  emm <- emmeans(model, specs = group)
  cl  <- as.data.frame(multcomp::cld(emm, Letters = letters,
                                     alpha = ALPHA, adjust = "tukey"))
  data.frame(group  = as.character(cl[[group]]),
             letter = trimws(cl$.group))
}

#### 4. Run site and plot tests for each species ####
results <- list()

for (sp in species_list) {

  sp_dat <- all_elytra %>% filter(scientificName_Species == sp)

  ## 4a. Among-site test ##
  sp_sites <- sites_ok %>% filter(scientificName_Species == sp) %>% pull(siteID)
  d_site   <- sp_dat %>% filter(siteID %in% sp_sites)

  site_res <- d_site %>%
    group_by(siteID) %>%
    summarise(n       = n(),
              mean_cm = mean(cm_elytra_max_length),
              n_plots = n_distinct(plotID),
              n_years = n_distinct(yearCollected),
              .groups = "drop") %>%
    mutate(species = sp, test = "site", test_unit = "all sites",
           group = siteID, letter = NA_character_,
           stat = NA_real_, df = NA_real_, p = NA_real_)

  if (length(sp_sites) >= 2) {
    m1  <- lmer(cm_elytra_max_length ~ siteID + (1 | plotID), data = d_site, REML = FALSE)
    m0  <- lmer(cm_elytra_max_length ~ 1 + (1 | plotID), data = d_site, REML = FALSE)
    lrt <- anova(m0, m1)

    site_letters <- get_letters(update(m1, REML = TRUE), "siteID")
    site_res <- site_res %>%
      select(-letter) %>%
      left_join(site_letters, by = "group") %>%
      mutate(stat = lrt$Chisq[2], df = lrt$Df[2], p = lrt$`Pr(>Chisq)`[2])
  }
  results[[length(results) + 1]] <- site_res

  ## 4b. Among-plot test within each site ##
  test_sites <- plot_test_sites %>%
    filter(scientificName_Species == sp) %>%
    pull(siteID)

  for (s in test_sites) {
    s_plots <- plots_ok %>%
      filter(scientificName_Species == sp, siteID == s) %>%
      pull(plotID)
    d_s <- sp_dat %>% filter(plotID %in% s_plots)

    m <- lm(cm_elytra_max_length ~ plotID, data = d_s)
    a <- anova(m)

    plot_res <- d_s %>%
      group_by(siteID, plotID) %>%
      summarise(n       = n(),
                mean_cm = mean(cm_elytra_max_length),
                n_years = n_distinct(yearCollected),
                .groups = "drop") %>%
      mutate(species = sp, test = "plot", test_unit = s, group = plotID) %>%
      left_join(get_letters(m, "plotID"), by = "group") %>%
      mutate(stat = a$`F value`[1], df = a$Df[1], p = a$`Pr(>F)`[1])

    results[[length(results) + 1]] <- plot_res
  }
}

#### 5. Combine results, FDR-correct across species, write ####
results_df <- bind_rows(results)

fdr <- results_df %>%
  filter(!is.na(p)) %>%
  distinct(test, species, test_unit, p) %>%
  group_by(test) %>%
  mutate(p_fdr = p.adjust(p, method = "BH")) %>%
  ungroup()

results_df <- results_df %>%
  left_join(fdr, by = c("test", "species", "test_unit", "p")) %>%
  select(any_of(c("species", "test", "test_unit", "siteID", "plotID", "group", "n",
                  "n_plots", "n_years", "mean_cm", "letter", "stat", "df", "p", "p_fdr")))

write.csv(results_df,
          paste0("./Outputs/BodySizeSpatialTests_", param_tag, ".csv"),
          row.names = FALSE)

# How many tests are significant before and after correction
fdr %>%
  group_by(test) %>%
  summarise(n_tests   = n(),
            n_sig_raw = sum(p < ALPHA),
            n_sig_fdr = sum(p_fdr < ALPHA))

# Flag sites with a single plot: site and plot effects are confounded there
results_df %>% filter(test == "site", !is.na(p), n_plots == 1)

length(species_list)
table(results_df$test)
site_diff_results_df<-subset(results_df, test=="site")
table(site_diff_results_df$test_unit)
min(site_diff_results_df$n)

site_diff_results_df<-subset(site_diff_results_df, !is.na(stat))
hist(site_diff_results_df$p)
table(ifelse(site_diff_results_df$p>0.05,"same","different"))
site_diff_results_df$diffClass<-ifelse(site_diff_results_df$p>0.05,"same","different")
table(ifelse(site_diff_results_df$p>0.05,"same","different"))
dim(table(subset(site_diff_results_df, diffClass=="same")$species))
dim(table(subset(site_diff_results_df, diffClass=="different")$species))

plot_diff_results_df<-subset(results_df, test=="plot")
plot_diff_results_df<-subset(plot_diff_results_df, !is.na(stat))
hist(plot_diff_results_df$p)
plot_diff_results_df$diffClass<-ifelse(plot_diff_results_df$p>0.05,"same","different")
table(ifelse(plot_diff_results_df$p>0.05,"same","different"))
dim(table(subset(plot_diff_results_df, diffClass=="same")$species))
dim(table(subset(plot_diff_results_df, diffClass=="different")$species))

#### 6. Figures ####
tab10 <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
           "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC")

fmt_p <- function(p) format.pval(p, digits = 2, eps = 0.001)

for (sp in species_list) {

  sp_dat   <- all_elytra %>% filter(scientificName_Species == sp)
  site_res <- results_df %>% filter(species == sp, test == "site")
  plot_res <- results_df %>% filter(species == sp, test == "plot")

  # Blank letters for tests that are not significant after FDR
  if (!SHOW_NS_LETTERS) {
    site_res <- site_res %>% mutate(letter = ifelse(!is.na(p_fdr) & p_fdr < ALPHA, letter, NA))
    plot_res <- plot_res %>% mutate(letter = ifelse(!is.na(p_fdr) & p_fdr < ALPHA, letter, NA))
  }

  # Colors: one hue per site, plots are shades of their site's hue
  all_sites <- sort(unique(c(site_res$siteID, plot_res$siteID)))
  site_cols <- if (length(all_sites) <= 10) tab10[seq_along(all_sites)] else colorRampPalette(tab10)(length(all_sites))
  names(site_cols) <- all_sites

  plot_cols <- c()
  for (s in unique(plot_res$siteID)) {
    s_plots <- sort(plot_res$plotID[plot_res$siteID == s])
    k <- length(s_plots)
    shades <- colorRampPalette(c("white", site_cols[[s]], "black"))(k + 4)[3:(k + 2)]
    plot_cols <- c(plot_cols, setNames(shades, s_plots))
  }

  n_col <- min(4, max(nrow(site_res), length(unique(plot_res$siteID))))

  ## 6a. Site panel ##
  d_site <- sp_dat %>% filter(siteID %in% site_res$siteID)

  site_lab <- d_site %>%
    group_by(siteID) %>%
    summarise(ymax = max(density(cm_elytra_max_length)$y), .groups = "drop") %>%
    left_join(site_res %>% select(siteID, mean_cm, letter), by = "siteID")

  site_sub <- if (any(!is.na(site_res$p))) {
    sprintf("Sites: LRT chi2(%d) = %.2f, p = %s, FDR p = %s",
            as.integer(site_res$df[1]), site_res$stat[1],
            fmt_p(site_res$p[1]), fmt_p(site_res$p_fdr[1]))
  } else {
    "Sites: fewer than 2 sites meet MIN_N_SITE, not tested"
  }

  p_site <- ggplot(d_site, aes(x = cm_elytra_max_length, fill = siteID)) +
    geom_density(alpha = 0.5) +
    geom_rug(aes(colour = siteID)) +
    geom_vline(data = site_res, aes(xintercept = mean_cm, colour = siteID)) +
    geom_text(data = site_lab, aes(x = mean_cm, y = ymax * 1.12, label = letter),
              inherit.aes = FALSE, fontface = "bold", size = 5, na.rm = TRUE) +
    scale_fill_manual(values = site_cols) +
    scale_colour_manual(values = site_cols) +
    labs(x = NULL, y = "Site Density", title = sp, subtitle = site_sub) +
    theme_pubr() +
    theme(legend.position = "none") +
    facet_wrap(~ siteID, ncol = n_col)

  fig        <- p_site
  rows_total <- ceiling(nrow(site_res) / n_col)

  ## 6b. Plot panel (only sites with >= 2 qualifying plots) ##
  if (nrow(plot_res) > 0) {
    d_plot <- sp_dat %>% filter(plotID %in% plot_res$plotID)

    # Stack plot letters above the tallest density in each site, ordered by mean
    plot_lab <- d_plot %>%
      group_by(siteID, plotID) %>%
      summarise(dmax = max(density(cm_elytra_max_length)$y), .groups = "drop") %>%
      group_by(siteID) %>%
      mutate(ymax = max(dmax)) %>%
      ungroup() %>%
      left_join(plot_res %>% select(plotID, mean_cm, letter), by = "plotID") %>%
      group_by(siteID) %>%
      arrange(mean_cm, .by_group = TRUE) %>%
      mutate(y = ymax * (1.08 + 0.10 * (row_number() - 1)),
             label = ifelse(is.na(letter), NA, paste0(sub(".*_", "", plotID), " ", letter))) %>%
      ungroup()

    plot_p <- plot_res %>%
      distinct(siteID, stat, df, p, p_fdr) %>%
      mutate(label = sprintf("Plots: F(%d) = %.2f\np = %s, FDR p = %s",
                             as.integer(df), stat, fmt_p(p), fmt_p(p_fdr)))

    p_plot <- ggplot(d_plot, aes(x = cm_elytra_max_length, fill = plotID)) +
      geom_density(alpha = 0.2) +
      geom_rug(aes(colour = plotID)) +
      geom_vline(data = plot_res, aes(xintercept = mean_cm, colour = plotID)) +
      geom_text(data = plot_lab, aes(x = mean_cm, y = y, label = label, colour = plotID),
                inherit.aes = FALSE, fontface = "bold", size = 3, na.rm = TRUE) +
      geom_text(data = plot_p, aes(x = -Inf, y = Inf, label = label),
                inherit.aes = FALSE, hjust = -0.05, vjust = 1.2, size = 2.8) +
      scale_fill_manual(values = plot_cols) +
      scale_colour_manual(values = plot_cols) +
      labs(x = "Elytra Length (cm)", y = "Plot Density") +
      theme_pubr() +
      theme(legend.position = "none") +
      facet_wrap(~ siteID, ncol = n_col)

    fig        <- p_site / p_plot
    rows_total <- rows_total + ceiling(length(unique(plot_res$siteID)) / n_col)
  }

  ggsave(paste0(fig_dir, gsub(" ", "_", sp), "_SitePlotTests_", param_tag, ".png"),
         fig, width = 2.5 * n_col + 1, height = 2.8 * rows_total + 1,
         units = "in", dpi = 300)
}

df<-all_elytra %>% filter(scientificName_Species %in% species_list)
head(df)

# df <- df %>%
#   group_by(siteID, plotID, scientificName_Species) %>%
#   summarise(n       = n(),
#             mean_cm = mean(cm_elytra_max_length),
#             .groups = "drop")

df<-subset(df, n>=20)

plot_richness<-read.csv("../BeetleBiodiversity/plot_annual_EstimatedSppRichness.csv")
head(plot_richness)
plot_richness<-subset(plot_richness, Year==2018)
plot_richness$X<-NULL
plotDF<-merge(df, plot_richness, by.x="plotID", by.y = "PlotID", all.x = TRUE, all.y = FALSE)
head(plotDF)

phylo<-read.csv("./Outputs/plot_PhyloDistance.csv")
phylo<-subset(phylo, Year==2018)
plotDF<-merge(plotDF, phylo, by="plotID")
head(plotDF)

struc<-read.csv("./Outputs/BETplot_Rugosity.csv")
struc$X<-NULL
env<-read.csv("./Outputs/BeetlePlotswEnvData.csv")
NPP<-read.csv("../NEON_MODIS_NPP_2018_2019.csv") #from https://code.earthengine.google.com/b41a55076352b2d9e21ac5e74bf337bc
velocity<-read.csv("./Outputs/BeetlePlotswVelocity.csv")
head(velocity)
velocity<-velocity[,c("plotID","Velocity")]

plotDF<-merge(plotDF, struc, by="plotID")
plotDF<-merge(plotDF, env, by="plotID")
plotDF<-merge(plotDF, NPP[,c("Npp","Gpp","plotID")], by="plotID")
plotDF<-merge(plotDF, velocity, by="plotID")
head(plotDF)

m<-glm(cm_elytra_max_length~scientificName_Species+Estimator+mntd+rugosity_RC+bio_1+bio_12+geodiv, data = plotDF)
#m<-glm(mean_cm~scientificName_Species*bio1+Estimator, data = plotDF)

plot(m)
summary(m)

library(jtools)
effect_plot(m, pred = bio_1)
effect_plot(m, pred = bio_12)
effect_plot(m, pred = rugosity_RC)
effect_plot(m, pred = Estimator)
effect_plot(m, pred = rugosity_RC)
