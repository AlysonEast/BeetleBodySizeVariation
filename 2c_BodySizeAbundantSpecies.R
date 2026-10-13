################################################################################
# 2c_IntraspecificShiftDrivers.R
#
# What explains shifts in a species' local mean body size away from its
# species-wide mean? One cross-species mixed model of individual elytra length.
#
# RESPONSE
#   log_len_c = log(elytra length) minus the species' mean log length.
#   Coefficients are proportional shifts: 100 * (exp(b) - 1) = % change in
#   local mean size per 1 SD of the predictor.
#
# PREDICTORS (all z-scored at the level they vary, not across beetles)
#   Climate (TEMP_COL), NPP, rugosity : per plot, split into site mean (_site)
#     and the plot's deviation from it (_dev), so among-site and within-site
#     effects are estimated separately.
#   Richness : per plot-year, split the same way.
#   Nearest relative : per focal species x plot-year.
#     nn_dist = phylogenetic distance to the closest co-occurring species in
#       the tree (community = species with an effort-scaled abundance in that
#       plot-year, matching 6b_PhyloDistance_ByYear.R).
#     nn_diff = relative's mean log size minus the focal species' mean log
#       size. The relative's size leaves out the focal plot, so it is never
#       measured on the same beetles' neighbors.
#     Character displacement predicts a negative nn_diff effect (shift away
#       from the relative's size) that is strongest when nn_dist is small,
#       i.e. a positive nn_dist:nn_diff interaction.
#
# RANDOM EFFECTS
#   (1 | siteID) + (1 | plotID)      community-wide size differences
#   (1 | species:siteID) + (1 | species:plotID)
#                                    species-specific local means; species x
#                                    plot is the real unit of replication for
#                                    plot-level predictors, not the beetle
#   (0 + slopes || species)          species-specific responses. A random
#                                    slope variance > 0 means species differ
#                                    in whether and how much they shift.
#
# MODELS
#   m_env : climate, NPP, rugosity, richness, year
#   m_nn  : m_env + nn_dist * nn_diff, on rows with nearest-relative data
#
# OUTPUTS (filenames carry param_tag)
#   ./Outputs/IntraspecificShift_Data_<tag>.csv
#   ./Outputs/IntraspecificShift_Fixed_<tag>.csv
#   ./Outputs/IntraspecificShift_SpeciesSlopes_<tag>.csv
#   ./Figures/BodySizeQuantification/IntraspecificShift/*.png
################################################################################

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ape)
library(lme4)
library(lmerTest)   # Satterthwaite p-values and ranova(); masks lme4::lmer

#### 0. User defined variables ####
YEARS            <- c(2018,2019)
TEMP_COL         <- "bio_1"       # or "Comp.1" for climate PC1 if temp and NPP are collinear
MIN_N_SPECIES    <- 50            # min individuals for a species to enter
MIN_SITES        <- 3             # min sites per species (needed for among-site slopes)
MIN_PLOTS        <- 5             # min plots per species
MIN_COMPLETENESS <- 0             # drop plot-years with richness completeness below this (0 = keep all)
#RANDOM_SLOPES    <- c("temp_site", "rich_dev")  # start small; add terms if they converge

TREE_FILE <- "./Data/Tree.DNA.droppedseqs.newick"

setwd("/home/aly/Beetles/BeetleBodySizeVariation")
fig_dir <- "./Figures/BodySizeQuantification/IntraspecificShift/"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
param_tag <- paste0(TEMP_COL, "_minN", MIN_N_SPECIES, "_minSites", MIN_SITES,
                    "_minPlots", MIN_PLOTS, "_comp", MIN_COMPLETENESS)

# z-score helper, applied to plot, plot-year, and species x plot-year tables
z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

#### 1. Individuals (same exclusions as 2_BodySizeQuantification.R) ####
all_elytra <- read.csv("./Data/BodysizeCombinedClean.csv")
all_elytra <- subset(all_elytra, individualID != "NEON.BET.D10.016322")
all_elytra <- subset(all_elytra, individualID != "NEON.BET.D08.002280")
all_elytra <- subset(all_elytra, imageID != "MLBS_009.S.20180522.jpg")
all_elytra <- subset(all_elytra, imageID != "MLBS_009.E.20180522.CARABIDS.01.jpg")
all_elytra <- subset(all_elytra, individualID != "NEON.BET.D07.003475")

# Same name fix as 6b so species match tree tips
name_fix <- c("Amara crassipina" = "Amara crassispina")
hit <- all_elytra$scientificName_Species %in% names(name_fix)
all_elytra$scientificName_Species[hit] <- name_fix[all_elytra$scientificName_Species[hit]]

ind <- all_elytra %>%
  filter(yearCollected %in% YEARS,
         !grepl("sp\\.", scientificName_Species),
         !is.na(scientificName_Species),
         scientificName_Species != "",
         is.finite(cm_elytra_max_length),
         cm_elytra_max_length > 0) %>%
  transmute(individualID, species = scientificName_Species, domainID, siteID, plotID,
            year = yearCollected, Assemblage = paste0(plotID, "_", yearCollected),
            log_len = log(cm_elytra_max_length))
head(ind)

#### 2. Plot-level environment (static per plot) ####
env  <- read.csv("./Outputs/BeetlePlotswEnvData.csv")
struc <- read.csv("./Outputs/BETplot_Rugosity.csv")
NPP  <- read.csv("../NEON_MODIS_NPP_2018_2019.csv")

plot_env <- ind %>%
  distinct(siteID, plotID) %>%
  left_join(env %>% distinct(plotID, .keep_all = TRUE) %>%
              transmute(plotID, temp = .data[[TEMP_COL]]), by = "plotID") %>%
  left_join(env %>% distinct(plotID, .keep_all = TRUE) %>%
              transmute(plotID, precip = bio_12), by = "plotID") %>%
  left_join(env %>% distinct(plotID, .keep_all = TRUE) %>%
              transmute(plotID, Comp.1 = Comp.1), by = "plotID") %>%
  left_join(env %>% distinct(plotID, .keep_all = TRUE) %>%
              transmute(plotID, Comp.2 = Comp.2), by = "plotID") %>%
  left_join(NPP %>% distinct(plotID, .keep_all = TRUE) %>%
              transmute(plotID, npp = Npp), by = "plotID") %>%
  left_join(struc %>% distinct(plotID, .keep_all = TRUE) %>%
              transmute(plotID, rug = rugosity_RC), by = "plotID")%>%
  left_join(struc %>% distinct(plotID, .keep_all = TRUE) %>%
              transmute(plotID, geodiv = geodiv), by = "plotID")

# Plots missing any environmental variable (e.g. no LiDAR)
plot_env %>% filter(is.na(temp) | is.na(npp) | is.na(rug))
plot_env <- plot_env %>% filter(!is.na(temp), !is.na(npp), !is.na(rug))

#### 3. Plot-year richness ####
plot_rich <- read.csv("../BeetleBiodiversity/plot_annual_EstimatedSppRichness.csv")
plot_rich$X <- NULL

rich <- ind %>%
  distinct(siteID, plotID, year, Assemblage) %>%
  left_join(plot_rich %>% transmute(Assemblage, rich = Estimator, completeness),
            by = "Assemblage")

#### 4. Nearest co-occurring relative (species x plot-year) ####
tree <- read.tree(TREE_FILE)
tree$tip.label <- gsub("_", " ", sub("_\\d+$", "", tree$tip.label))
bad_tip <- !grepl("^[A-Z][a-z]+ [a-z]", tree$tip.label)
if (any(bad_tip)) tree <- drop.tip(tree, which(bad_tip))
dup <- duplicated(tree$tip.label)
if (any(dup)) tree <- drop.tip(tree, which(dup))
cophen <- cophenetic(tree)

# Community per plot-year: species with an effort-scaled abundance (as in 6b)
comm <- bind_rows(lapply(YEARS, function(yr) {
  read.csv(sprintf("./Data/plot_abund_%d.csv", yr)) %>%
    mutate(year = yr)
})) %>%
  mutate(scientificName_Species = ifelse(scientificName_Species %in% names(name_fix),
                                         name_fix[scientificName_Species],
                                         scientificName_Species),
         Assemblage = paste0(plotID, "_", year)) %>%
  filter(scientificName_Species %in% tree$tip.label)

# Species-wide mean log size, plus per-plot sums for the leave-plot-out size
sp_mean <- ind %>% group_by(species) %>% summarise(sp_mean_log = mean(log_len), .groups = "drop")
sp_tot  <- ind %>% group_by(species) %>% summarise(sum_all = sum(log_len), n_all = n(), .groups = "drop")
sp_plot <- ind %>% group_by(species, plotID) %>% summarise(sum_plot = sum(log_len), n_plot = n(), .groups = "drop")

focal <- ind %>% distinct(species, plotID, Assemblage)

nn_rows <- vector("list", nrow(focal))
for (i in seq_len(nrow(focal))) {
  f  <- focal$species[i]
  pl <- focal$plotID[i]
  others <- setdiff(comm$scientificName_Species[comm$Assemblage == focal$Assemblage[i]], f)
  
  nn <- NA_character_; nn_dist <- NA_real_; nn_size <- NA_real_
  if (f %in% rownames(cophen) && length(others) > 0) {
    d       <- cophen[f, others]
    nn      <- others[which.min(d)]
    nn_dist <- min(d)
    
    # Relative's mean log size, excluding its beetles from the focal plot
    tot <- sp_tot %>% filter(species == nn)
    pp  <- sp_plot %>% filter(species == nn, plotID == pl)
    s   <- tot$sum_all - sum(pp$sum_plot)
    n   <- tot$n_all - sum(pp$n_plot)
    if (length(n) == 1 && n > 0) nn_size <- s / n
  }
  nn_rows[[i]] <- data.frame(species = f, Assemblage = focal$Assemblage[i],
                             nn = nn, nn_dist = nn_dist, nn_size = nn_size)
}

nn_tab <- bind_rows(nn_rows) %>%
  left_join(sp_mean, by = "species") %>%
  mutate(nn_diff = nn_size - sp_mean_log,
         nn_dist_z = z(nn_dist),
         nn_diff_z = z(nn_diff)) %>%
  select(species, Assemblage, nn, nn_dist, nn_diff, nn_dist_z, nn_diff_z)

# Coverage: how many focal species x plot-years have a sized relative
nn_tab %>% summarise(n = n(), has_nn = sum(!is.na(nn)), has_size = sum(!is.na(nn_diff)))

#### Nearest sized species in community ####
# Size distance to the nearest other species at the same plot (plot-level species means)
sp_plot_means <- ind %>%
  group_by(plotID, species) %>%
  summarise(plot_mean_log = mean(log_len), .groups = "drop")

nn_size <- sp_plot_means %>%
  group_by(plotID) %>%
  mutate(nn_size_diff = sapply(seq_along(plot_mean_log), function(i) {
    d <- plot_mean_log[-i] - plot_mean_log[i]
    if (length(d) == 0) NA_real_ else d[which.min(abs(d))]
  }),
  nn_species = sapply(seq_along(species), function(i) {
    others <- species[-i]
    if (length(others) == 0) NA_character_ else others[which.min(abs(plot_mean_log[-i] - plot_mean_log[i]))]
  })) %>%
  ungroup()


#### 5. Assemble modeling table ####
dat_full <- ind %>%
  inner_join(plot_env %>% select(-siteID), by = "plotID") %>%
  inner_join(rich %>% select(Assemblage, rich, completeness), by = "Assemblage") %>%
  left_join(nn_tab, by = c("species", "Assemblage")) %>%
  left_join(sp_mean, by = "species") %>%
  left_join(nn_size, by = c("plotID", "species")) %>% 
  mutate(log_len_c = log_len - sp_mean_log,
         year = factor(year))
dim(dat_full)
head(dat_full)

cat("Individuals:", nrow(dat_full), " Species:", n_distinct(dat_full$species),
    " Plots:", n_distinct(dat_full$plotID), " Sites:", n_distinct(dat_full$siteID),
    " Species x plots:", n_distinct(paste(dat_full$species, dat_full$plotID)), "\n")

#### 6. Per-species models: settings ####
library(lmerTest)
library(MuMIn)          # r.squaredGLMM
library(broom.mixed)    # tidy() fixed effects with Satterthwaite p and Wald CIs
library(ggeffects)      # predict_response()
library(patchwork)
library(flextable)
library(tidyr)
library(svglite)

N_TOP_SPECIES <- 10
MODEL_TERMS   <- c("Comp.1", "Comp.2", "geodiv", "rug", "rich",
                   "nn_dist", "nn_diff", "nn_size_diff")
RANDOM_EFFECTS <- c(between = "(1 | year)",
                    within  = "(1 | siteID) + (1 | year)")
MODEL_LABELS   <- c(between = "Across sites", within = "Within sites")

pairs.panels(dat_full[,MODEL_TERMS])
dat_full$rug<-log10(dat_full$rug)
dat_full$geodiv<-log10(dat_full$geodiv)
dat_full$rich<-log10(dat_full$rich)
pairs.panels(dat_full[,MODEL_TERMS])

# Display labels and predictor groups, in table/figure order
term_info <- data.frame(
  term  = c("(Intercept)", "Comp.1", "Comp.2", "geodiv", "rug", "rich",
            "nn_dist", "nn_diff", "nn_dist:nn_diff", "nn_size_diff"),
  label = c("Intercept", "Climate PC1", "Climate PC2", "Log 10 Geodiversity",
            "Log 10 Canopy rugosity", "Log 10 Species richness",
            "Phylo. distance to nearest relative",
            "Size difference from nearest relative",
            "Phylo. distance x size difference",
            "Size difference from nearest-sized species"),
  group = c("", "Climate", "Climate", "Habitat", "Habitat", "Competition",
            "Competition", "Competition", "Competition", "Competition"))

TAB_COLS <- c("#4E79A7", "#F28E2B", "#E15759")   # Tableau 10

table_dir <- "./Outputs/IntraspecificShift_BySpecies/"
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

top_spp <- names(sort(table(dat_full$species), decreasing = TRUE)[1:N_TOP_SPECIES])
top_spp

#### 6a. Helpers (reused for every species x model) ####
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
fmt_p_eq <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
wrap_lab <- function(x) paste(strwrap(x, 32), collapse = "\n")

re_sd <- function(m, grp) {
  vc <- as.data.frame(VarCorr(m))
  x  <- vc$sdcor[vc$grp == grp]
  if (length(x) == 1) x else NA_real_
}

# Main effect of one predictor, back-transformed to original units and cm.
# Other predictors held at their mean (0 after scaling); population-level (no REs).
plot_effect <- function(tm, m, dat, scl, fx) {
  mu <- scl$mean[scl$term == tm]; s <- scl$sd[scl$term == tm]
  pr <- as.data.frame(predict_response(m, terms = tm, margin = "mean_reference")) %>%
    mutate(x = x * s + mu, across(c(predicted, conf.low, conf.high), exp))
  b  <- fx[fx$term == tm, ]
  ggplot(pr, aes(x, predicted)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "grey85") +
    geom_line(linewidth = 0.8, linetype = ifelse(b$p.value < 0.05, "solid", "dashed")) +
    geom_rug(data = data.frame(x = dat[[tm]] * s + mu), aes(x = x),
             inherit.aes = FALSE, sides = "b", alpha = 0.15) +
    labs(x = wrap_lab(term_info$label[term_info$term == tm]), y = "Elytra length (cm)",
         subtitle = sprintf("b = %.3f, %s", b$estimate, fmt_p_eq(b$p.value))) +
    theme_pubr(base_size = 9)
}

# nn_diff effect at low / mean / high phylogenetic distance to the nearest relative
plot_interaction <- function(m, dat, scl, fx) {
  mu_d <- scl$mean[scl$term == "nn_diff"]; s_d <- scl$sd[scl$term == "nn_diff"]
  mu_p <- scl$mean[scl$term == "nn_dist"]; s_p <- scl$sd[scl$term == "nn_dist"]
  pr <- as.data.frame(predict_response(m, terms = c("nn_diff", "nn_dist [-1, 0, 1]"),
                                       margin = "mean_reference")) %>%
    mutate(x = x * s_d + mu_d,
           across(c(predicted, conf.low, conf.high), exp),
           group = factor(group, levels = c("-1", "0", "1"),
                          labels = sprintf(c("-1 SD (%.2f)", "Mean (%.2f)", "+1 SD (%.2f)"),
                                           mu_p + c(-1, 0, 1) * s_p)))
  b <- fx[fx$term == "nn_dist:nn_diff", ]
  ggplot(pr, aes(x, predicted, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = TAB_COLS, name = "Phylo. distance") +
    scale_fill_manual(values = TAB_COLS, name = "Phylo. distance") +
    labs(x = wrap_lab("Size difference from nearest relative (log ratio)"), y = "Elytra length (cm)",
         subtitle = sprintf("Interaction b = %.3f, %s", b$estimate, fmt_p_eq(b$p.value))) +
    theme_pubr(base_size = 9, legend = "bottom") +
    guides(colour = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))
}

#### 6b. Loop over species ####
fixed_all <- list()   # tidy fixed effects, all species x models
fit_all   <- list()   # fit statistics, all species x models
tables    <- list()   # one flextable per species

for (sp in top_spp) {
  sp_tag <- gsub(" ", "_", sp)
  cat("\n====", sp, "====\n")
  
  # Complete cases first, so scaling and both models use the same rows
  dat <- dat_full %>%
    filter(species == sp) %>%
    filter(if_all(all_of(MODEL_TERMS), ~ !is.na(.x)))
  cat("Complete-case individuals:", nrow(dat), "\n")
  if (nrow(dat) < 30) { cat("Skipping: too few complete cases\n"); next }
  
  # Drop predictors with no variation for this species
  keep <- MODEL_TERMS[sapply(dat[MODEL_TERMS], function(x) isTRUE(sd(x) > 0))]
  if (length(keep) < length(MODEL_TERMS))
    cat("Dropped (no variation):", setdiff(MODEL_TERMS, keep), "\n")
  
  scl <- data.frame(term = keep,
                    mean = sapply(dat[keep], mean),
                    sd   = sapply(dat[keep], sd))
  dat[keep] <- lapply(dat[keep], z)
  
  fixed <- paste(keep, collapse = " + ")
  if (all(c("nn_dist", "nn_diff") %in% keep)) fixed <- paste(fixed, "+ nn_dist:nn_diff")
  
  fx_sp <- list(); fit_sp <- list()
  
  for (mod in names(RANDOM_EFFECTS)) {
    f <- as.formula(paste("log_len ~", fixed, "+", RANDOM_EFFECTS[[mod]]))
    m <- tryCatch(lmer(f, data = dat, REML = TRUE),
                  error = function(e) { cat("FAILED", mod, ":", conditionMessage(e), "\n"); NULL })
    if (is.null(m)) next
    print(summary(m))
    
    fx <- tidy(m, effects = "fixed", conf.int = TRUE) %>%
      mutate(species = sp, model = mod, pct_change = 100 * (exp(estimate) - 1))
    r2 <- suppressWarnings(r.squaredGLMM(m))
    fs <- data.frame(species = sp, model = mod, n = nobs(m),
                     n_sites = n_distinct(dat$siteID), n_plots = n_distinct(dat$plotID),
                     R2m = r2[1, "R2m"], R2c = r2[1, "R2c"],
                     sd_site = re_sd(m, "siteID"), sd_year = re_sd(m, "year"),
                     sd_resid = sigma(m), singular = isSingular(m),
                     dropped = paste(setdiff(MODEL_TERMS, keep), collapse = ", "),
                     row.names = NULL)
    
    fx_sp[[mod]] <- fx; fit_sp[[mod]] <- fs
    fixed_all[[paste(sp, mod)]] <- fx; fit_all[[paste(sp, mod)]] <- fs
    
    # Effect and interaction panels
    main_terms <- intersect(term_info$term, keep)
    panels <- lapply(main_terms, plot_effect, m = m, dat = dat, scl = scl, fx = fx)
    if ("nn_dist:nn_diff" %in% fx$term)
      panels <- c(panels, list(plot_interaction(m, dat, scl, fx)))
    
    fig <- wrap_plots(panels, ncol = 3) +
      plot_annotation(title = bquote(italic(.(sp)) * ":" ~ .(MODEL_LABELS[[mod]])),
                      tag_levels = "a")
    fn <- paste0(fig_dir, "Effects_", sp_tag, "_", mod, "_", param_tag)
    h  <- 2.8 * ceiling(length(panels) / 3)
    ggsave(paste0(fn, ".png"), fig, width = 10, height = h, dpi = 300)
    ggsave(paste0(fn, ".svg"), fig, width = 10, height = h)
  }
  if (length(fx_sp) == 0) next
  
  # Publication table: both models side by side
  mods <- names(fx_sp)
  coef_rows <- bind_rows(fx_sp) %>%
    left_join(term_info, by = "term") %>%
    transmute(group, label, model,
              est = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high),
              p   = fmt_p(p.value)) %>%
    pivot_wider(names_from = model, values_from = c(est, p)) %>%
    arrange(match(label, term_info$label))
  
  fit_vals <- function(fs) {
    v <- c(as.character(c(fs$n, fs$n_sites, fs$n_plots)),
           sprintf("%.3f", c(fs$R2m, fs$R2c, fs$sd_site, fs$sd_year, fs$sd_resid)),
           ifelse(fs$singular, "yes", "no"))
    replace(v, v == "NA", "")
  }
  fit_rows <- data.frame(group = "Model fit",
                         label = c("N individuals", "Sites", "Plots",
                                   "Marginal R2", "Conditional R2",
                                   "SD site", "SD year", "SD residual", "Singular fit"))
  for (mod in mods) {
    fit_rows[[paste0("est_", mod)]] <- fit_vals(fit_sp[[mod]])
    fit_rows[[paste0("p_", mod)]]   <- ""
  }
  
  tab <- bind_rows(coef_rows, fit_rows)
  tab <- tab[, c("group", "label", as.vector(rbind(paste0("est_", mods), paste0("p_", mods))))]
  write.csv(tab, paste0(table_dir, "Table_", sp_tag, "_", param_tag, ".csv"), row.names = FALSE)
  
  hdr <- setNames(as.list(c("", "Predictor", rep(c("Estimate (95% CI)", "p"), length(mods)))),
                  names(tab))
  tables[[sp]] <- flextable(tab) %>%
    set_header_labels(values = hdr) %>%
    add_header_row(values = c("", "", MODEL_LABELS[mods]),
                   colwidths = c(1, 1, rep(2, length(mods)))) %>%
    merge_v(j = "group") %>%
    valign(j = "group", valign = "top") %>%
    bold(j = "group") %>%
    align(i = 1, part = "header", align = "center") %>%
    theme_booktabs() %>%
    hline(i = nrow(coef_rows), border = officer::fp_border(width = 0.75)) %>%
    fontsize(size = 9, part = "all") %>%
    set_caption(as_paragraph(
      "Linear mixed models of log elytra length for ", as_i(sp),
      ". Predictors z-scored within species; estimates are change in log length per 1 SD ",
      "(100 x (exp(b) - 1) = % change). Satterthwaite p-values, Wald 95% CIs.")) %>%
    autofit()
}

#### 6c. Combined outputs ####
fixed_df <- bind_rows(fixed_all) %>% left_join(term_info, by = "term")
fit_df   <- bind_rows(fit_all)
# write.csv(fixed_df, paste0("./Outputs/IntraspecificShift_TopSpp_Fixed_", param_tag, ".csv"), row.names = FALSE)
# write.csv(fit_df,   paste0("./Outputs/IntraspecificShift_TopSpp_Fit_", param_tag, ".csv"), row.names = FALSE)

save_as_docx(values = setNames(tables, paste0("Table S", seq_along(tables), ". ", names(tables))),
             path = paste0(table_dir, "Tables_TopSpp_", param_tag, ".docx"))

# Cross-species coefficient plot: climate vs habitat vs competition, both scales
coef_plot <- fixed_df %>%
  filter(term != "(Intercept)") %>%
  mutate(label   = factor(label, levels = term_info$label),
         species = factor(species, levels = rev(top_spp)),
         model   = factor(MODEL_LABELS[model], levels = MODEL_LABELS)) %>%
  ggplot(aes(estimate, species, colour = model)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                  position = position_dodge(width = 0.6), size = 0.2) +
  facet_wrap(~ label, ncol = 3, labeller = label_wrap_gen(30)) +
  scale_colour_manual(values = TAB_COLS[1:2], name = NULL) +
  labs(x = "Effect on log elytra length (per 1 SD)", y = NULL) +
  theme_pubr(base_size = 9) +
  theme(axis.text.y = element_text(face = "italic"),
        panel.spacing.x = unit(1.2, "lines")) +
  xlim(-0.15,0.15)

ggsave(paste0(fig_dir, "Coefficients_TopSpp_", param_tag, ".png"), coef_plot, width = 9, height = 9, dpi = 300)
ggsave(paste0(fig_dir, "Coefficients_TopSpp_", param_tag, ".svg"), coef_plot, width = 9, height = 9)
coef_plot
