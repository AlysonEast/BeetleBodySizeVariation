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
YEARS           <- c(2018, 2019)  # each year is analyzed as its own sample unit

emm_options(lmer.df = "satterthwaite")  # needs lmerTest installed (not loaded)

setwd("/home/aly/Beetles/BeetleBodySizeVariation")
fig_dir <- "./Figures/BodySizeQuantification/SpatialTests/"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
param_tag <- paste0("minSite", MIN_N_SITE, "_minPlot", MIN_N_PLOT,
                    "_yr", paste(YEARS, collapse = "-"),
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

all_elytra <- all_elytra %>% filter(yearCollected %in% YEARS)

#### 2. Identify species x year combinations with enough data ####
site_n <- all_elytra %>%
  count(scientificName_Species, yearCollected, domainID, siteID, name = "n")
plot_n <- all_elytra %>%
  count(scientificName_Species, yearCollected, domainID, siteID, plotID, name = "n")

sites_ok <- site_n %>% filter(n >= MIN_N_SITE)
plots_ok <- plot_n %>% filter(n >= MIN_N_PLOT)

# Site test needs >= 2 qualifying sites in the same year
site_test_spp <- sites_ok %>%
  count(scientificName_Species, yearCollected, name = "n_sites") %>%
  filter(n_sites >= 2)

# Plot test needs >= 2 qualifying plots within a site in the same year
plot_test_sites <- plots_ok %>%
  count(scientificName_Species, yearCollected, siteID, name = "n_plots") %>%
  filter(n_plots >= 2)

# Species x year combinations with at least one test
spp_years <- bind_rows(site_test_spp %>% select(scientificName_Species, yearCollected),
                       plot_test_sites %>% select(scientificName_Species, yearCollected)) %>%
  distinct() %>%
  arrange(scientificName_Species, yearCollected)

species_list <- sort(unique(spp_years$scientificName_Species))  # kept for code after section 6
length(species_list)
table(spp_years$yearCollected)

#### 3. Compact letter display helper (used by both the site and plot tests) ####
get_letters <- function(model, group) {
  emm <- emmeans(model, specs = group)
  cl  <- as.data.frame(multcomp::cld(emm, Letters = letters,
                                     alpha = ALPHA, adjust = "tukey"))
  data.frame(group  = as.character(cl[[group]]),
             letter = trimws(cl$.group))
}

#### 4. Run site and plot tests for each species x year ####
results <- list()

for (i in seq_len(nrow(spp_years))) {
  
  sp <- spp_years$scientificName_Species[i]
  yr <- spp_years$yearCollected[i]
  sp_dat <- all_elytra %>% filter(scientificName_Species == sp, yearCollected == yr)
  
  ## 4a. Among-site test ##
  sp_sites <- sites_ok %>%
    filter(scientificName_Species == sp, yearCollected == yr) %>%
    pull(siteID)
  d_site <- sp_dat %>% filter(siteID %in% sp_sites)
  
  site_res <- d_site %>%
    group_by(siteID) %>%
    summarise(n       = n(),
              mean_cm = mean(cm_elytra_max_length),
              n_plots = n_distinct(plotID),
              .groups = "drop") %>%
    mutate(species = sp, year = yr, test = "site", test_unit = "all sites",
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
    filter(scientificName_Species == sp, yearCollected == yr) %>%
    pull(siteID)
  
  for (s in test_sites) {
    s_plots <- plots_ok %>%
      filter(scientificName_Species == sp, yearCollected == yr, siteID == s) %>%
      pull(plotID)
    d_s <- sp_dat %>% filter(plotID %in% s_plots)
    
    m <- lm(cm_elytra_max_length ~ plotID, data = d_s)
    a <- anova(m)
    
    plot_res <- d_s %>%
      group_by(siteID, plotID) %>%
      summarise(n       = n(),
                mean_cm = mean(cm_elytra_max_length),
                .groups = "drop") %>%
      mutate(species = sp, year = yr, test = "plot", test_unit = s, group = plotID) %>%
      left_join(get_letters(m, "plotID"), by = "group") %>%
      mutate(stat = a$`F value`[1], df = a$Df[1], p = a$`Pr(>F)`[1])
    
    results[[length(results) + 1]] <- plot_res
  }
}

#### 5. Combine results, FDR-correct across species, write ####
results_df <- bind_rows(results)

fdr <- results_df %>%
  filter(!is.na(p)) %>%
  distinct(test, species, year, test_unit, p) %>%
  group_by(test) %>%
  mutate(p_fdr = p.adjust(p, method = "BH")) %>%
  ungroup()

results_df <- results_df %>%
  left_join(fdr, by = c("test", "species", "year", "test_unit", "p")) %>%
  select(any_of(c("species", "year", "test", "test_unit", "siteID", "plotID", "group", "n",
                  "n_plots", "mean_cm", "letter", "stat", "df", "p", "p_fdr")))

write.csv(results_df,
          paste0("./Outputs/BodySizeSpatialTests_", param_tag, ".csv"),
          row.names = FALSE)

#### 5b. How many species vary between sites and between plots ####
results_df

site_results_df<-subset(results_df, test=="site")
site_results_df<-subset(site_results_df,!is.na(p))
dim(table(site_results_df$species))
site_results_df

site_results_summary <- site_results_df %>%
  group_by(species, year) %>%
  summarise(n_sites = n_distinct(siteID),
            n_groups = n_distinct(letter),
            .groups = "drop")
site_results_summary$sig<-ifelse(site_results_summary$n_groups>1, 1, 0)

table(site_results_summary$n_sites,site_results_summary$sig)
site_results_summary %>%
  group_by(sig) %>%
  summarise(median(n_sites))

site_results_summary <- site_results_summary %>%
  group_by(species) %>%
  summarise(n_years = n_distinct(year),
            min_sites = min(n_sites),
            max_sites = max(n_sites),
            min_groups = min(n_groups),
            max_groups = max(n_groups),
            n_sig = sum(sig),
            pct_sig = n_sig/n_years,
            .groups = "drop")


site_results_summary
table(site_results_summary$n_years)

site_results_summary_1year<-subset(site_results_summary, n_years==1)
dim(site_results_summary_1year)
table(site_results_summary_1year$pct_sig)
table(site_results_summary_1year$pct_sig,site_results_summary_1year$max_sites)

site_results_summary_2year<-subset(site_results_summary, n_years==2)
dim(site_results_summary_2year)
table(site_results_summary_2year$pct_sig)
table(site_results_summary_2year$pct_sig,site_results_summary_2year$max_sites)


plot_results_df<-subset(results_df, test=="plot")
plot_results_df<-subset(plot_results_df,!is.na(p))
dim(table(plot_results_df$species))
plot_results_df

plot_results_summary <- plot_results_df %>%
  group_by(species, year) %>%
  summarise(n_plots = n_distinct(plotID),
            n_groups = n_distinct(letter),
            .groups = "drop")
plot_results_summary$sig<-ifelse(plot_results_summary$n_groups>1, 1, 0)

table(plot_results_summary$n_plots,plot_results_summary$sig)
plot_results_summary %>%
  group_by(sig) %>%
  summarise(median(n_plots))

plot_results_summary <- plot_results_summary %>%
  group_by(species) %>%
  summarise(n_years = n_distinct(year),
            min_plots = min(n_plots),
            max_plots = max(n_plots),
            min_groups = min(n_groups),
            max_groups = max(n_groups),
            n_sig = sum(sig),
            pct_sig = n_sig/n_years,
            .groups = "drop")


plot_results_summary
table(plot_results_summary$n_years)

plot_results_summary_1year<-subset(plot_results_summary, n_years==1)
dim(plot_results_summary_1year)
table(plot_results_summary_1year$pct_sig)
table(plot_results_summary_1year$pct_sig,plot_results_summary_1year$max_plots)

plot_results_summary_2year<-subset(plot_results_summary, n_years==2)
dim(plot_results_summary_2year)
table(plot_results_summary_2year$pct_sig)
table(plot_results_summary_2year$pct_sig,plot_results_summary_2year$max_plots)


results_effect_sizes <- results_df %>%
  group_by(species, year, test, test_unit) %>%
  summarise(n_obs = n(),
            diff = max(mean_cm)-min(mean_cm),
            p = p,
            n_groups = n_distinct(letter),
            sig = ifelse(n_groups>1, 1, 0),
            .groups = "drop")

results_effect_sizes %>% 
  group_by(test, sig) %>%
  summarise(mean(diff),
            .groups = "drop")

SIG <- "sig_fdr"  # or "sig_raw" for uncorrected p < ALPHA

# One row per test (results_df has one row per site or plot)
tests <- results_df %>%
  filter(!is.na(p)) %>%
  distinct(test, species, year, test_unit, stat, df, p, p_fdr) %>%
  mutate(sig_raw = p < ALPHA,
         sig_fdr = p_fdr < ALPHA)

# Tests and significant tests by test type and year
tests %>%
  group_by(test, year) %>%
  summarise(n_tests   = n(),
            n_species = n_distinct(species),
            n_sig_raw = sum(sig_raw),
            n_sig_fdr = sum(sig_fdr),
            .groups = "drop")

# Replication across years for each test unit
# (site test: unit = species; plot test: unit = species x site)
unit_consistency <- tests %>%
  group_by(test, species, test_unit) %>%
  summarise(n_years_tested = n(),
            n_years_sig    = sum(.data[[SIG]]),
            .groups = "drop") %>%
  mutate(class = case_when(n_years_sig == 0               ~ "never differs",
                           n_years_tested == 1            ~ "differs, 1 year tested",
                           n_years_sig == n_years_tested  ~ "differs every year",
                           TRUE                           ~ "differs some years"))
table(unit_consistency$test, unit_consistency$class)

# One class per species: the strongest evidence across its units
spp_class <- unit_consistency %>%
  group_by(test, species) %>%
  summarise(n_units          = n(),
            n_units_sig      = sum(n_years_sig > 0),
            n_units_every    = sum(class == "differs every year"),
            max_years_tested = max(n_years_tested),
            .groups = "drop") %>%
  mutate(class = case_when(n_units_sig == 0      ~ "never differs",
                           n_units_every > 0     ~ "differs every year (>=1 unit)",
                           max_years_tested == 1 ~ "differs, 1 year tested",
                           TRUE                  ~ "differs some years only"))
table(spp_class$test, spp_class$class)

#### 6. Figures (one per species, years overlaid) ####
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
    s_plots <- sort(unique(plot_res$plotID[plot_res$siteID == s]))
    k <- length(s_plots)
    shades <- colorRampPalette(c("white", site_cols[[s]], "black"))(k + 4)[3:(k + 2)]
    plot_cols <- c(plot_cols, setNames(shades, s_plots))
  }
  
  n_col <- min(4, max(n_distinct(site_res$siteID), n_distinct(plot_res$siteID)))
  
  ## 6a. Site panel ##
  # Densities pool the years in which each site entered the site test
  d_site <- sp_dat %>%
    semi_join(site_res, by = c("siteID", "yearCollected" = "year"))
  
  # One letter per site x year, stacked by year above the tallest density
  site_lab <- d_site %>%
    group_by(siteID) %>%
    summarise(ymax = max(density(cm_elytra_max_length)$y), .groups = "drop") %>%
    left_join(site_res %>% select(siteID, year, mean_cm, letter), by = "siteID") %>%
    mutate(y   = ymax * (1.10 + 0.12 * (match(year, YEARS) - 1)),
           col = year_cols[as.character(year)])
  
  site_tests <- site_res %>%
    filter(!is.na(p)) %>%
    distinct(year, stat, df, p, p_fdr) %>%
    arrange(year)
  
  site_sub <- if (nrow(site_tests) > 0) {
    paste(sprintf("%s sites: LRT chi2(%d) = %.2f, p = %s, FDR p = %s",
                  site_tests$year, as.integer(site_tests$df), site_tests$stat,
                  fmt_p(site_tests$p), fmt_p(site_tests$p_fdr)),
          collapse = "\n")
  } else {
    "Sites: fewer than 2 sites meet MIN_N_SITE in any year, not tested"
  }
  
  p_site <- ggplot(d_site, aes(x = cm_elytra_max_length, fill = siteID)) +
    geom_density(alpha = 0.5) +
    geom_rug(aes(colour = siteID)) +
    geom_vline(data = site_res,
               aes(xintercept = mean_cm, colour = siteID, linetype = as.character(year))) +
    geom_text(data = site_lab, aes(x = mean_cm, y = y, label = letter),
              colour = site_lab$col, inherit.aes = FALSE,
              fontface = "bold", size = 5, na.rm = TRUE) +
    scale_fill_manual(values = site_cols) +
    scale_colour_manual(values = site_cols) +
    scale_linetype_manual(values = year_lines) +
    labs(x = NULL, y = "Site Density", title = sp, subtitle = site_sub) +
    theme_pubr() +
    theme(legend.position = "none") +
    facet_wrap(~ siteID, ncol = n_col)
  
  fig        <- p_site
  rows_total <- ceiling(n_distinct(site_res$siteID) / n_col)
  
  ## 6b. Plot panel (only sites with >= 2 qualifying plots in a year) ##
  if (nrow(plot_res) > 0) {
    d_plot <- sp_dat %>%
      semi_join(plot_res, by = c("plotID", "yearCollected" = "year"))
    
    # Letters stacked above the tallest density in each site: first year,
    # then second, and within a year ordered by plot mean (low to high)
    plot_lab <- d_plot %>%
      group_by(siteID, plotID) %>%
      summarise(dmax = max(density(cm_elytra_max_length)$y), .groups = "drop") %>%
      group_by(siteID) %>%
      summarise(ymax = max(dmax), .groups = "drop") %>%
      left_join(plot_res %>% select(siteID, plotID, year, mean_cm, letter), by = "siteID") %>%
      group_by(siteID) %>%
      arrange(year, mean_cm, .by_group = TRUE) %>%
      mutate(y = ymax * (1.08 + 0.10 * (row_number() - 1))) %>%
      ungroup() %>%
      mutate(col = year_cols[as.character(year)])
    
    # One stat line per site x year, coloured by year
    plot_p <- plot_res %>%
      distinct(siteID, year, stat, df, p, p_fdr) %>%
      arrange(siteID, year) %>%
      group_by(siteID) %>%
      mutate(vj = 1.2 + 1.3 * (row_number() - 1)) %>%
      ungroup() %>%
      mutate(label = sprintf("%s plots: F(%d) = %.2f, p = %s, FDR p = %s",
                             year, as.integer(df), stat, fmt_p(p), fmt_p(p_fdr)),
             col   = year_cols[as.character(year)])
    
    p_plot <- ggplot(d_plot, aes(x = cm_elytra_max_length, fill = plotID)) +
      geom_density(alpha = 0.2) +
      geom_rug(aes(colour = plotID)) +
      geom_vline(data = plot_res,
                 aes(xintercept = mean_cm, colour = plotID, linetype = as.character(year))) +
      geom_text(data = plot_lab, aes(x = mean_cm, y = y, label = letter),
                colour = plot_lab$col, inherit.aes = FALSE,
                fontface = "bold", size = 3.5, na.rm = TRUE) +
      geom_text(data = plot_p, aes(x = -Inf, y = Inf, label = label, vjust = vj),
                colour = plot_p$col, inherit.aes = FALSE, hjust = -0.05, size = 2.5) +
      scale_fill_manual(values = plot_cols) +
      scale_colour_manual(values = plot_cols) +
      scale_linetype_manual(values = year_lines) +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +  # headroom for stat text
      labs(x = "Elytra Length (cm)", y = "Plot Density") +
      theme_pubr() +
      theme(legend.position = "none") +
      facet_wrap(~ siteID, ncol = n_col)
    
    fig        <- p_site / p_plot
    rows_total <- rows_total + ceiling(n_distinct(plot_res$siteID) / n_col)
  }
  
  fig <- fig + plot_annotation(
    caption = paste0("Letters: ", paste(names(year_cols), year_cols, sep = " = ", collapse = ", "),
                     ". Mean lines: ", paste(names(year_lines), year_lines, sep = " = ", collapse = ", "),
                     ". Densities pool the years tested."))
  
  ggsave(paste0(fig_dir, gsub(" ", "_", sp), "_SitePlotTests_", param_tag, ".png"),
         fig, width = 2.5 * n_col + 1, height = 2.8 * rows_total + 1,
         units = "in", dpi = 300)
}

#### 7. Single figure: one species, one year, letters only ####
sp <- "Pasimachus elongatus"
yr <- 2018

sp_dat   <- all_elytra %>% filter(scientificName_Species == sp, yearCollected == yr)
site_res <- results_df %>% filter(species == sp, year == yr, test == "site")
plot_res <- results_df %>% filter(species == sp, year == yr, test == "plot")

# Blank letters for tests that are not significant after FDR
if (!SHOW_NS_LETTERS) {
  site_res <- site_res %>% mutate(letter = ifelse(!is.na(p_fdr) & p_fdr < ALPHA, letter, NA))
  plot_res <- plot_res %>% mutate(letter = ifelse(!is.na(p_fdr) & p_fdr < ALPHA, letter, NA))
}

# Facet order: sites from smallest to largest mean elytra length
site_order <- site_res %>% arrange(mean_cm) %>% pull(siteID)
# Sites that only appear in the plot test go at the end
site_order <- c(site_order, setdiff(sort(unique(plot_res$siteID)), site_order))


# Colors: one hue per site, plots are shades of their site's hue
all_sites <- sort(unique(c(site_res$siteID, plot_res$siteID)))
site_cols <- if (length(all_sites) <= 10) tab10[seq_along(all_sites)] else colorRampPalette(tab10)(length(all_sites))
names(site_cols) <- all_sites

plot_cols <- c()
for (s in unique(plot_res$siteID)) {
  s_plots <- sort(unique(plot_res$plotID[plot_res$siteID == s]))
  k <- length(s_plots)
  shades <- colorRampPalette(c("white", site_cols[[s]], "black"))(k + 4)[3:(k + 2)]
  plot_cols <- c(plot_cols, setNames(shades, s_plots))
}

n_col <- min(4, max(n_distinct(site_res$siteID), n_distinct(plot_res$siteID)))

## 7a. Site panel ##
d_site <- sp_dat %>% filter(siteID %in% site_res$siteID)

d_site   <- d_site   %>% mutate(siteID = factor(siteID, levels = site_order))
site_res <- site_res %>% mutate(siteID = factor(siteID, levels = site_order))
site_lab <- site_lab %>% mutate(siteID = factor(siteID, levels = site_order))

site_lab <- d_site %>%
  group_by(siteID) %>%
  summarise(ymax = max(density(cm_elytra_max_length)$y), .groups = "drop") %>%
  left_join(site_res %>% select(siteID, mean_cm, letter), by = "siteID") %>%
  mutate(y = ymax * 1.10)

p_site <- ggplot(d_site, aes(x = cm_elytra_max_length, fill = siteID)) +
  geom_density(alpha = 0.5) +
  geom_rug(aes(colour = siteID)) +
  geom_vline(data = site_res, aes(xintercept = mean_cm, colour = siteID)) +
  geom_text(data = site_lab, aes(x = mean_cm, y = y, label = letter),
            colour = year_cols[as.character(yr)], inherit.aes = FALSE,
            fontface = "bold", size = 5, na.rm = TRUE) +
  scale_fill_manual(values = site_cols) +
  scale_colour_manual(values = site_cols) +
  labs(x = NULL, y = "Site Density", title = paste0(sp, " (", yr, ")")) +
  theme_pubr() +
  # facet_wrap(~ siteID, ncol = n_col) +
  # theme(legend.position = "none")
  theme(legend.position = c(0.2, 0.8))

fig        <- p_site
rows_total <- ceiling(n_distinct(site_res$siteID) / n_col)

## 7b. Plot panel (only sites with >= 2 qualifying plots) ##
if (nrow(plot_res) > 0) {
  d_plot <- sp_dat %>% filter(plotID %in% plot_res$plotID)
  
  # Letters stacked above the tallest density in each site, ordered by plot mean
  plot_lab <- d_plot %>%
    group_by(siteID, plotID) %>%
    summarise(dmax = max(density(cm_elytra_max_length)$y), .groups = "drop") %>%
    group_by(siteID) %>%
    summarise(ymax = max(dmax), .groups = "drop") %>%
    left_join(plot_res %>% select(siteID, plotID, mean_cm, letter), by = "siteID") %>%
    group_by(siteID) %>%
    arrange(mean_cm, .by_group = TRUE) %>%
    mutate(y = ymax * (1.08 + 0.10 * (row_number() - 1))) %>%
    ungroup()
  d_plot   <- d_plot   %>% mutate(siteID = factor(siteID, levels = site_order))
  plot_res <- plot_res %>% mutate(siteID = factor(siteID, levels = site_order))
  plot_lab <- plot_lab %>% mutate(siteID = factor(siteID, levels = site_order))
  
  
  p_plot <- ggplot(d_plot, aes(x = cm_elytra_max_length, fill = plotID)) +
    geom_density(alpha = 0.2) +
    geom_rug(aes(colour = plotID)) +
    geom_vline(data = plot_res, aes(xintercept = mean_cm, colour = plotID)) +
    geom_text(data = plot_lab, aes(x = mean_cm, y = y, label = letter),
              colour = year_cols[as.character(yr)], inherit.aes = FALSE,
              fontface = "bold", size = 3.5, na.rm = TRUE) +
    scale_fill_manual(values = plot_cols) +
    scale_colour_manual(values = plot_cols) +
    labs(x = "Elytra Length (cm)", y = "Plot Density") +
    theme_pubr() +
    theme(legend.position = "none") +
    facet_wrap(~ siteID, ncol = n_col)
  
  fig        <- p_site / p_plot
  rows_total <- rows_total + ceiling(n_distinct(plot_res$siteID) / n_col)
}

fig

ggsave(paste0(fig_dir, gsub(" ", "_", sp), "_", yr, "_SitePlotTests_LettersOnly_", param_tag, ".png"),
       fig, width = 3 * n_col + 1, height = 2 * rows_total + 1,
       units = "in", dpi = 300)

df<-all_elytra %>% filter(scientificName_Species %in% species_list)
head(df)

sort(table(df$scientificName_Species), decreasing = TRUE)[1:10]
table(subset(df, scientificName_Species=="Amara quenseli")$siteID)
table(subset(df, scientificName_Species=="Harpalus pensylvanicus")$siteID)
table(subset(df, scientificName_Species=="Chlaenius aestivus")$siteID)
table(subset(df, scientificName_Species=="Pterostichus melanarius")$siteID)
table(subset(df, scientificName_Species=="Cyclotrachelus torvus")$siteID)
table(subset(df, scientificName_Species=="Pasimachus elongatus")$siteID)
table(subset(df, scientificName_Species=="Pterostichus lachrymosus")$siteID)
table(subset(df, scientificName_Species=="Calathus advena")$siteID)
table(subset(df, scientificName_Species=="Synuchus impunctatus")$siteID)
table(subset(df, scientificName_Species=="Carabus goryi")$siteID)

# df <- df %>%
#   group_by(siteID, plotID, scientificName_Species) %>%
#   summarise(n       = n(),
#             mean_cm = mean(cm_elytra_max_length),
#             .groups = "drop")

# df<-subset(df, n>=20)

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

# m<-glm(cm_elytra_max_length~scientificName_Species+Estimator+mntd+rugosity_RC+bio_1+bio_12+geodiv, data = plotDF)
# #m<-glm(mean_cm~scientificName_Species*bio1+Estimator, data = plotDF)
# 
# plot(m)
# summary(m)
# 
# library(jtools)
# effect_plot(m, pred = bio_1)
# effect_plot(m, pred = bio_12)
# effect_plot(m, pred = rugosity_RC)
# effect_plot(m, pred = Estimator)
# effect_plot(m, pred = rugosity_RC)

df_sub<-subset(plotDF, scientificName_Species=="Carabus goryi")

m<-glm(cm_elytra_max_length~Estimator+mntd+rugosity_RC+bio_1+bio_12+geodiv, data = plotDF)

