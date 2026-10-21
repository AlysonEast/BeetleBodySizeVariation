################################################################################
# 2c_BodySizeAbundantSpecies.R
#
# What explains shifts in a species' local body size? Per-species linear mixed
# models of individual log elytra length for the N_TOP_SPECIES most abundant
# species, in a single year (YEARS; default 2018) so no year term is needed.
#
# RESPONSE
#   log_len = natural log of elytra length. Coefficients are change in log
#   length per 1 SD of the predictor: 100 * (exp(b) - 1) = % change.
#
# PREDICTORS (MODEL_TERMS, z-scored within species; groups set in term_info)
#   Environment : climate PC1-3, log10 geodiversity, log10 canopy rugosity
#   Competition : log10 richness, phylo. distance to nearest relative (nn_dist),
#                 size difference from that relative (nn_diff), nn_dist:nn_diff,
#                 size difference from nearest-sized species (nn_size_diff),
#                 log10 abundance share of that species (nn_size_abund),
#                 nn_size_diff:nn_size_abund
#   All size differences are leave-plot-out (section 4 / 4b).
#   Character displacement predicts a negative size-difference effect that is
#   strongest when the neighbour is closely related (positive nn_dist:nn_diff)
#   or abundant (negative nn_size_diff:nn_size_abund).
#
# MODELS (RANDOM_EFFECTS)
#   between : no random effect, lm()     drivers across sites
#   within  : (1 | siteID), lmer()       drivers within sites
#
# VARIANCE PARTITION (Borcard et al. 1992, on the marginal R2 scale)
#   Each model is refit with only the Environment terms and only the Competition
#   terms (same rows, same random effects). Unique env = R2(full) - R2(comp),
#   unique comp = R2(full) - R2(env), shared = R2(env) + R2(comp) - R2(full).
#   Random effects = R2c - R2m; residual = 1 - R2c.
#
# OUTPUTS (filenames carry param_tag, built in section 6 from N_TOP_SPECIES and
# MODEL_TERMS, so a new predictor set never overwrites an old run)
#   ./Outputs/IntraspecificShift_TopSpp_{Fixed,Fit,DiagFlags}_<tag>.csv
#   ./Outputs/IntraspecificShift_TopSpp_ManuscriptTable_<tag>.{csv,docx}
#   ./Outputs/IntraspecificShift_BySpecies/Table_<species>_<tag>.csv, Tables_TopSpp_<tag>.docx
#   ./Figures/BodySizeQuantification/IntraspecificShift/
#       Effects_<species>_<model>_<tag>.{png,svg}, Coefficients_TopSpp_<tag>.{png,svg},
#       Diagnostics_TopSpp_<tag>.pdf (one 2x2 page per species x model)
################################################################################

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ape)
library(lme4)
library(lmerTest)   # Satterthwaite p-values and ranova(); masks lme4::lmer

#### 0. User defined variables ####
YEARS            <- 2018         # single year: no year term in the models
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
# param_tag (output filename tag) is built in section 6, once MODEL_TERMS is set

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
soil <- read.csv("./Outputs/BETplot_Soils.csv")


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
  left_join(env %>% distinct(plotID, .keep_all = TRUE) %>%
              transmute(plotID, Comp.3 = Comp.3), by = "plotID") %>%
  # left_join(soil %>% distinct(plotID, .keep_all = TRUE) %>%
  #             transmute(plotID, soil_ph_0_5cm = phh2o), by = "plotID") %>%
  left_join(NPP %>% distinct(plotID, .keep_all = TRUE) %>%
              transmute(plotID, npp = Npp), by = "plotID") %>%
  left_join(struc %>% distinct(plotID, .keep_all = TRUE) %>%
              transmute(plotID, rug = rugosity_RC), by = "plotID") %>%
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

# Community per plot-year: species with an effort-scaled abundance (as in 6b).
# comm_all keeps every species-level taxon (used for the nearest-sized species
# and its abundance share); comm keeps tree tips only (nearest relative).
comm_all <- bind_rows(lapply(YEARS, function(yr) {
  read.csv(sprintf("./Data/plot_abund_%d.csv", yr)) %>%
    mutate(year = yr)
})) %>%
  mutate(scientificName_Species = ifelse(scientificName_Species %in% names(name_fix),
                                         name_fix[scientificName_Species],
                                         scientificName_Species),
         Assemblage = paste0(plotID, "_", year)) %>%
  filter(!is.na(scientificName_Species), scientificName_Species != "", abund > 0) %>%
  group_by(Assemblage) %>%
  mutate(abund_share = abund / sum(abund)) %>%   # share of the plot-year community
  ungroup()
comm <- comm_all %>% filter(scientificName_Species %in% tree$tip.label)

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

#### 4b. Nearest-sized co-occurring species (species x plot-year) ####
# For each focal species x plot-year, the co-occurring species (any in comm_all
# for that plot-year) whose mean size is closest to the focal species'.
# Leave-plot-out on BOTH sides, as for nn_diff: the focal species' size comes
# from its beetles in other plots, so the predictor is never built from the
# individuals in the response (the old plot-mean version was, which made a
# negative nn_size_diff effect partly mechanical). The neighbour's size also
# leaves out the focal plot, so a shared good or bad plot cannot drive both.
#   nn_size_diff  = neighbour mean log size - focal mean log size (log ratio)
#   nn_size_abund = log10 of the neighbour's abundance share in the plot-year
# Displacement predicts a negative nn_size_diff effect that strengthens with
# the neighbour's abundance, i.e. a negative nn_size_diff:nn_size_abund term.
lpo_size <- function(sp, pl) {          # leave-plot-out mean log size
  tot <- sp_tot[sp_tot$species %in% sp, ]
  pp  <- sp_plot[sp_plot$species %in% sp & sp_plot$plotID == pl, ]
  ps  <- pp$sum_plot[match(tot$species, pp$species)]; ps[is.na(ps)] <- 0
  pn  <- pp$n_plot[match(tot$species, pp$species)];   pn[is.na(pn)] <- 0
  n_  <- tot$n_all - pn
  setNames(ifelse(n_ > 0, (tot$sum_all - ps) / n_, NA_real_), tot$species)
}

ns_rows <- vector("list", nrow(focal))
for (i in seq_len(nrow(focal))) {
  f  <- focal$species[i]
  pl <- focal$plotID[i]
  cm <- comm_all[comm_all$Assemblage == focal$Assemblage[i] &
                   comm_all$scientificName_Species != f, ]
  
  out <- data.frame(species = f, Assemblage = focal$Assemblage[i], nn_size_sp = NA_character_,
                    nn_size_diff = NA_real_, nn_size_abund = NA_real_)
  f_size <- lpo_size(f, pl)
  o_size <- lpo_size(cm$scientificName_Species, pl)
  o_size <- o_size[is.finite(o_size)]
  if (length(f_size) == 1 && is.finite(f_size) && length(o_size) > 0) {
    d   <- o_size - f_size
    k   <- which.min(abs(d))
    out$nn_size_sp    <- names(o_size)[k]
    out$nn_size_diff  <- d[[k]]
    out$nn_size_abund <- log10(cm$abund_share[cm$scientificName_Species == names(o_size)[k]][1])
  }
  ns_rows[[i]] <- out
}
nn_sized <- bind_rows(ns_rows)

# Coverage: focal species x plot-years with a sized neighbour and its abundance
nn_sized %>% summarise(n = n(), has_nn = sum(!is.na(nn_size_diff)),
                       has_abund = sum(is.finite(nn_size_abund)))


#### 5. Assemble modeling table ####
dat_full <- ind %>%
  inner_join(plot_env %>% select(-siteID), by = "plotID") %>%
  inner_join(rich %>% select(Assemblage, rich, completeness), by = "Assemblage") %>%
  left_join(nn_tab, by = c("species", "Assemblage")) %>%
  left_join(sp_mean, by = "species") %>%
  left_join(nn_sized, by = c("species", "Assemblage")) %>%
  mutate(log_len_c = log_len - sp_mean_log)
dim(dat_full)
head(dat_full)

cat("Individuals:", nrow(dat_full), " Species:", n_distinct(dat_full$species),
    " Plots:", n_distinct(dat_full$plotID), " Sites:", n_distinct(dat_full$siteID),
    " Species x plots:", n_distinct(paste(dat_full$species, dat_full$plotID)), "\n")

#### 6. Per-species models: settings ####
library(lmerTest)
library(MuMIn)          # r.squaredGLMM
library(broom)          # tidy() for the OLS (no random effect) model
library(broom.mixed)    # tidy() fixed effects with Satterthwaite p and Wald CIs
library(ggeffects)      # predict_response()
library(patchwork)
library(flextable)
library(tidyr)
library(svglite)
library(psych)

N_TOP_SPECIES <- 10
# Comp.3 and rich dropped: species have at most 12 plots in 2018 and every
# predictor is plot-level, so the model must stay under 12 coefficients
MODEL_TERMS   <- c("Comp.1", "Comp.2", "geodiv", "rug", #"Comp.3", "rich", "soil_ph_0_5cm",
                   "nn_dist", "nn_diff", "nn_size_diff", "nn_size_abund")
# "" = no random effect, fit with lm(); otherwise lmer()
RANDOM_EFFECTS <- c(between = "",
                    within  = "(1 | siteID)")
# Interactions added when both terms survive. x = predictor on the plot's x axis,
# by = moderator shown at -1 SD / mean / +1 SD, term = coefficient name.
INTERACTIONS   <- list(
  phylo = c(x = "nn_diff",      by = "nn_dist",       term = "nn_dist:nn_diff"),
  abund = c(x = "nn_size_diff", by = "nn_size_abund", term = "nn_size_diff:nn_size_abund"))
MODEL_LABELS   <- c(between = "Across sites", within = "Within sites")
SIG_LEVEL      <- 0.05    # predictors listed in the manuscript table
SHOW_SIGN      <- FALSE   # TRUE appends (+) / (-) to each listed predictor

# Filename tag: everything that changes the per-species models
param_tag <- paste0("y", paste(YEARS, collapse = "-"), "_top", N_TOP_SPECIES, "_",
                    gsub("[._]", "", paste(MODEL_TERMS, collapse = "-")))
param_tag


pairs.panels(dat_full[,MODEL_TERMS])
dat_full$rug<-log10(dat_full$rug + 0.01)
dat_full$geodiv<-log10(dat_full$geodiv + 0.01)
dat_full$rich<-log10(dat_full$rich + 0.01)
pairs.panels(dat_full[,MODEL_TERMS])

top_spp <- names(sort(table(dat_full$species), decreasing = TRUE)[1:N_TOP_SPECIES])
top_spp
top_spp<-c("Amara quenseli",
       "Calathus advena",
       "Carabus goryi",
       "Chlaenius aestivus",
       "Cyclotrachelus torvus",
       "Harpalus pensylvanicus",
       "Pasimachus elongatus",
       "Pterostichus lachrymosus",
       "Pterostichus melanarius",
       "Synuchus impunctatus")
#spp<-names(sort(table(dat_full$species), decreasing = TRUE)[1:N_TOP_SPECIES])
pairs.panels(subset(dat_full, species==top_spp[1])[,MODEL_TERMS])
pairs.panels(subset(dat_full, species==top_spp[2])[,MODEL_TERMS])
pairs.panels(subset(dat_full, species==top_spp[3])[,MODEL_TERMS])
pairs.panels(subset(dat_full, species==top_spp[4])[,MODEL_TERMS])
pairs.panels(subset(dat_full, species==top_spp[5])[,MODEL_TERMS])
pairs.panels(subset(dat_full, species==top_spp[6])[,MODEL_TERMS])
pairs.panels(subset(dat_full, species==top_spp[7])[,MODEL_TERMS])
pairs.panels(subset(dat_full, species==top_spp[8])[,MODEL_TERMS])
pairs.panels(subset(dat_full, species==top_spp[9])[,MODEL_TERMS])
pairs.panels(subset(dat_full, species==top_spp[10])[,MODEL_TERMS])

# Display labels and predictor groups, in table/figure order
term_info <- data.frame(
  term  = c("(Intercept)", "Comp.1", "Comp.2", "Comp.3", "geodiv", "rug", #"soil_ph_0_5cm",
            "rich", "nn_dist", "nn_diff", "nn_dist:nn_diff",
            "nn_size_diff", "nn_size_abund", "nn_size_diff:nn_size_abund"),
  label = c("Intercept", "Climate PC1", "Climate PC2", "Climate PC3", "Log 10 Geodiversity",
            "Log 10 Canopy rugosity", #"Soil pH", 
            "Log 10 Species richness",
            "Phylo. distance to nearest relative",
            "Size difference from nearest relative",
            "Phylo. distance x size difference",
            "Size difference from nearest-sized species",
            "Log 10 abundance share of nearest-sized species",
            "Size difference x abundance of nearest-sized species"),
  group = c("", "Environment", "Environment", "Environment", "Environment", "Environment", #"Environment",
            "Competition", "Competition", "Competition", "Competition", "Competition",
            "Competition", "Competition"))

TAB_COLS <- c("#4E79A7", "#F28E2B", "#E15759")   # Tableau 10

table_dir <- "./Outputs/IntraspecificShift_BySpecies/"
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

#### 6a. Helpers (reused for every species x model) ####
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
fmt_p_eq <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
wrap_lab <- function(x) paste(strwrap(x, 32), collapse = "\n")

# One fitting path for both models: lm() when the random-effect string is
# empty, lmer() otherwise. R2 for lm: marginal = conditional = ordinary R2.
fit_model <- function(fixed, re, dat) {
  if (re == "") return(lm(as.formula(paste("log_len ~", fixed)), data = dat))
  lmer(as.formula(paste("log_len ~", fixed, "+", re)), data = dat, REML = TRUE)
}
r2_of <- function(m) {
  if (inherits(m, "merMod")) return(suppressWarnings(r.squaredGLMM(m)))
  r2 <- summary(m)$r.squared
  cbind(R2m = r2, R2c = r2)
}

# Coefficients that could not be estimated (rank-deficient design). With every
# predictor measured per plot, this happens when a species has about as many
# plots as the model has coefficients.
aliased_terms <- function(m) {
  if (inherits(m, "merMod")) names(attr(getME(m, "X"), "col.dropped"))
  else names(which(is.na(coef(m))))
}

re_sd <- function(m, grp) {
  if (!inherits(m, "merMod")) return(NA_real_)
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

# Effect of ix["x"] at -1 SD / mean / +1 SD of the moderator ix["by"]
plot_interaction <- function(ix, m, dat, scl, fx) {
  xt <- ix[["x"]]; bt <- ix[["by"]]
  mu_x <- scl$mean[scl$term == xt]; s_x <- scl$sd[scl$term == xt]
  mu_b <- scl$mean[scl$term == bt]; s_b <- scl$sd[scl$term == bt]
  by_lab <- term_info$label[term_info$term == bt]
  pr <- as.data.frame(predict_response(m, terms = c(xt, paste0(bt, " [-1, 0, 1]")),
                                       margin = "mean_reference")) %>%
    mutate(x = x * s_x + mu_x,
           across(c(predicted, conf.low, conf.high), exp),
           group = factor(group, levels = c("-1", "0", "1"),
                          labels = sprintf(c("-1 SD (%.2f)", "Mean (%.2f)", "+1 SD (%.2f)"),
                                           mu_b + c(-1, 0, 1) * s_b)))
  b <- fx[fx$term == ix[["term"]], ]
  ggplot(pr, aes(x, predicted, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = TAB_COLS, name = wrap_lab(by_lab)) +
    scale_fill_manual(values = TAB_COLS, name = wrap_lab(by_lab)) +
    labs(x = wrap_lab(paste(term_info$label[term_info$term == xt], "(log ratio)")),
         y = "Elytra length (cm)",
         subtitle = sprintf("Interaction b = %.3f, %s", b$estimate, fmt_p_eq(b$p.value))) +
    theme_pubr(base_size = 9, legend = "bottom") +
    guides(colour = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))
}

# Marginal R2 of a refit with only `terms` as fixed effects (same rows, same
# random effects). Used for the Borcard partition. 0 if the group is empty.
r2m_subset <- function(terms, re, dat) {
  if (length(terms) == 0) return(0)
  m_sub <- tryCatch(fit_model(paste(terms, collapse = " + "), re, dat),
                    error = function(e) NULL)
  if (is.null(m_sub)) return(NA_real_)
  r2_of(m_sub)[1, "R2m"]
}

# 2x2 diagnostics in the style of par(mfrow = c(2, 2)); plot(lm).
# For lmer, residuals and fitted values are conditional (include the random
# effects), as in plot(lmerMod); for lm this reproduces plot(lm). Standardized residuals = resid / (sigma * sqrt(1 - h)), with
# h from the mixed-model hat matrix, so Cook's distance is an approximation.
# Points are labelled by row number in `dat`; the three most extreme per panel
# are returned with their individualID so they can be looked up afterwards.
# Leverage-1 rows (e.g. the only beetle in a plot when every predictor is
# plot-level) have no finite standardized residual; they are left off the plots
# and returned with flag = "leverage 1".
diag_plot <- function(m, dat, main) {
  fit <- fitted(m); res <- resid(m); h <- hatvalues(m)
  rs   <- res / (sigma(m) * sqrt(1 - h))
  p    <- if (inherits(m, "merMod")) length(fixef(m)) else m$rank   # estimable coefficients
  cook <- rs^2 * h / (p * (1 - h))
  lev1 <- which(!is.finite(rs) | h > 1 - 1e-8)
  ok   <- setdiff(seq_along(rs), lev1)
  fit <- fit[ok]; res <- res[ok]; h <- h[ok]; rs <- rs[ok]; cook <- cook[ok]
  top_r <- order(abs(rs), decreasing = TRUE)[1:3]
  top_c <- order(cook,    decreasing = TRUE)[1:3]
  pt <- grey(0, 0.35)

  op <- par(mfrow = c(2, 2), oma = c(0, 0, 2, 0), mar = c(4, 4, 2, 1))
  on.exit(par(op))

  # 1. Residuals vs Fitted: linearity
  plot(fit, res, col = pt, pch = 16, cex = 0.6,
       xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
  abline(h = 0, lty = 3, col = "grey50")
  lines(lowess(fit, res), col = "red")
  text(fit[top_r], res[top_r], ok[top_r], pos = 4, cex = 0.7)

  # 2. Normal Q-Q: residual normality
  qq <- qqnorm(rs, col = pt, pch = 16, cex = 0.6,
               ylab = "Standardized residuals", main = "Normal Q-Q")
  qqline(rs, lty = 3, col = "grey50")
  text(qq$x[top_r], qq$y[top_r], ok[top_r], pos = 4, cex = 0.7)

  # 3. Scale-Location: homoskedasticity
  sl <- sqrt(abs(rs))
  plot(fit, sl, col = pt, pch = 16, cex = 0.6, xlab = "Fitted values",
       ylab = expression(sqrt("|Standardized residuals|")), main = "Scale-Location")
  lines(lowess(fit, sl), col = "red")
  text(fit[top_r], sl[top_r], ok[top_r], pos = 4, cex = 0.7)

  # 4. Residuals vs Leverage with Cook's distance 0.5 / 1 contours: influence
  plot(h, rs, col = pt, pch = 16, cex = 0.6, xlim = c(0, max(h)),
       xlab = "Leverage", ylab = "Standardized residuals", main = "Residuals vs Leverage")
  abline(h = 0, v = 0, lty = 3, col = "grey50")
  lines(lowess(h, rs), col = "red")
  hh <- seq(max(h) / 100, max(h), length.out = 100)
  for (cd in c(0.5, 1)) {
    lines(hh,  sqrt(cd * p * (1 - hh) / hh), lty = 2, col = "red")
    lines(hh, -sqrt(cd * p * (1 - hh) / hh), lty = 2, col = "red")
  }
  text(h[top_c], rs[top_c], ok[top_c], pos = 4, cex = 0.7)

  mtext(main, outer = TRUE, font = 3)

  # Rows are reported by their position in `dat`, which is also the plot label
  flag <- sort(unique(c(top_r, top_c)))
  out <- data.frame(row = ok[flag], individualID = dat$individualID[ok[flag]],
                    plotID = dat$plotID[ok[flag]], std_resid = rs[flag], leverage = h[flag],
                    cooks = cook[flag], flag = "extreme", row.names = NULL)
  if (length(lev1) > 0)
    out <- rbind(out, data.frame(row = lev1, individualID = dat$individualID[lev1],
                                 plotID = dat$plotID[lev1], std_resid = NA, leverage = 1,
                                 cooks = NA, flag = "leverage 1"))
  out
}

#### 6b. Loop over species ####
fixed_all <- list()   # tidy fixed effects, all species x models
fit_all   <- list()   # fit statistics, all species x models
tables    <- list()   # one flextable per species
diag_flags <- list()  # most extreme residual / Cook's rows, all species x models

# One diagnostics page per species x model. If the loop errors out, run
# dev.off(diag_dev) before re-running so the PDF is closed.
pdf(paste0(fig_dir, "Diagnostics_TopSpp_", param_tag, ".pdf"), width = 8, height = 8)
diag_dev <- dev.cur()

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
  
  fixed_terms <- keep
  for (ix in INTERACTIONS)
    if (all(c(ix[["x"]], ix[["by"]]) %in% keep)) fixed_terms <- c(fixed_terms, ix[["term"]])
  fixed <- paste(fixed_terms, collapse = " + ")
  # Every predictor varies at the plot level, so plots (not beetles) bound how
  # many coefficients the data can support
  cat(sprintf("Plots: %d, sites: %d, coefficients: %d (%.1f plots per coefficient)\n",
              n_distinct(dat$plotID), n_distinct(dat$siteID), length(fixed_terms) + 1,
              n_distinct(dat$plotID) / (length(fixed_terms) + 1)))
  
  # Predictor groups for the variance partition (from term_info)
  env_terms  <- intersect(fixed_terms, term_info$term[term_info$group == "Environment"])
  comp_terms <- intersect(fixed_terms, term_info$term[term_info$group == "Competition"])
  
  fx_sp <- list(); fit_sp <- list()
  
  for (mod in names(RANDOM_EFFECTS)) {
    m <- tryCatch(fit_model(fixed, RANDOM_EFFECTS[[mod]], dat),
                  error = function(e) { cat("FAILED", mod, ":", conditionMessage(e), "\n"); NULL })
    if (is.null(m)) next
    print(summary(m))
    aliased <- aliased_terms(m)
    if (length(aliased) > 0) cat("RANK DEFICIENT", mod, ": not estimable:", aliased, "\n")
    
    fx <- tidy(m, effects = "fixed", conf.int = TRUE) %>%
      filter(is.finite(estimate)) %>%   # lm() keeps aliased terms as NA rows
      mutate(species = sp, model = mod, pct_change = 100 * (exp(estimate) - 1))
    r2 <- r2_of(m)
    
    # Borcard partition of the marginal R2 into environment vs competition
    r2_env  <- r2m_subset(env_terms,  RANDOM_EFFECTS[[mod]], dat)
    r2_comp <- r2m_subset(comp_terms, RANDOM_EFFECTS[[mod]], dat)
    fs <- data.frame(species = sp, model = mod, n = nobs(m),
                     n_sites = n_distinct(dat$siteID), n_plots = n_distinct(dat$plotID),
                     R2m = r2[1, "R2m"], R2c = r2[1, "R2c"],
                     sd_site = re_sd(m, "siteID"),
                     sd_resid = sigma(m), singular = inherits(m, "merMod") && isSingular(m),
                     R2m_env = r2_env, R2m_comp = r2_comp,
                     env_unique  = r2[1, "R2m"] - r2_comp,
                     comp_unique = r2[1, "R2m"] - r2_env,
                     shared      = r2_env + r2_comp - r2[1, "R2m"],
                     random      = r2[1, "R2c"] - r2[1, "R2m"],
                     residual    = 1 - r2[1, "R2c"],
                     dropped = paste(setdiff(MODEL_TERMS, keep), collapse = ", "),
                     n_coef = length(fixed_terms) + 1,
                     aliased = paste(aliased, collapse = ", "),
                     row.names = NULL)
    
    fx_sp[[mod]] <- fx; fit_sp[[mod]] <- fs
    
    # Diagnostics page
    dev.set(diag_dev)
    diag_flags[[paste(sp, mod)]] <- tryCatch(
      diag_plot(m, dat, paste0(sp, ": ", MODEL_LABELS[[mod]])) %>%
        mutate(species = sp, model = mod, .before = 1),
      error = function(e) { cat("DIAGNOSTICS FAILED", mod, ":", conditionMessage(e), "\n"); NULL })
    fixed_all[[paste(sp, mod)]] <- fx; fit_all[[paste(sp, mod)]] <- fs
    
    # Effect and interaction panels
    main_terms <- intersect(intersect(term_info$term, keep), fx$term)   # skips aliased terms
    panels <- lapply(main_terms, plot_effect, m = m, dat = dat, scl = scl, fx = fx)
    for (ix in INTERACTIONS)
      if (ix[["term"]] %in% fx$term)
        panels <- c(panels, list(plot_interaction(ix, m, dat, scl, fx)))
    
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
           sprintf("%.3f", c(fs$R2m, fs$R2c, fs$env_unique, fs$comp_unique, fs$shared,
                             fs$sd_site, fs$sd_resid)),
           ifelse(fs$singular, "yes", "no"), fs$aliased)
    replace(v, v == "NA", "")
  }
  fit_rows <- data.frame(group = "Model fit",
                         label = c("N individuals", "Sites", "Plots",
                                   "Marginal R2", "Conditional R2",
                                   "Environment (unique)", "Competition (unique)", "Shared",
                                   "SD site", "SD residual", "Singular fit", "Not estimable"))
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
      "(100 x (exp(b) - 1) = % change). Across sites: OLS, t-test p-values; within sites: ",
      "site random intercept, Satterthwaite p-values. Wald 95% CIs.")) %>%
    autofit()
}

dev.off(diag_dev)

#### 6c. Combined outputs ####
fixed_df <- bind_rows(fixed_all) %>% left_join(term_info, by = "term")
fit_df   <- bind_rows(fit_all)
write.csv(fixed_df, paste0("./Outputs/IntraspecificShift_TopSpp_Fixed_", param_tag, ".csv"), row.names = FALSE)
write.csv(fit_df,   paste0("./Outputs/IntraspecificShift_TopSpp_Fit_", param_tag, ".csv"), row.names = FALSE)
write.csv(bind_rows(diag_flags), paste0("./Outputs/IntraspecificShift_TopSpp_DiagFlags_", param_tag, ".csv"),
          row.names = FALSE)

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
  coord_cartesian(xlim = c(-0.15, 0.15))   # zooms; xlim() silently dropped wide CIs

ggsave(paste0(fig_dir, "Coefficients_TopSpp_", param_tag, ".png"), coef_plot, width = 9, height = 9, dpi = 300)
ggsave(paste0(fig_dir, "Coefficients_TopSpp_", param_tag, ".svg"), coef_plot, width = 9, height = 9)
coef_plot

#### 6d. Manuscript table: model variance + significant predictors ####
# One row per species, matching the "Body Size Mixed Model" columns of the
# intraspecific variation table: Model variance | Predictors, for each model.
# Cells carry line breaks, so in the docx each item sits on its own line.
stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))

var_cells <- fit_df %>%
  transmute(species, model,
            variance = paste0(sprintf("Marginal R2 %.3f\nConditional R2 %.3f\nEnvironment %.3f\nCompetition %.3f\nShared %.3f",
                                      R2m, R2c, env_unique, comp_unique, shared),
                              ifelse(singular, "\nSingular fit", ""),
                              ifelse(aliased != "", "\nRank deficient", "")))

pred_cells <- fixed_df %>%
  filter(term != "(Intercept)", p.value < SIG_LEVEL) %>%
  arrange(species, model, match(term, term_info$term)) %>%
  mutate(item = paste0(label, if (SHOW_SIGN) ifelse(estimate > 0, " (+)", " (-)") else "",
                       " ", stars(p.value))) %>%
  group_by(species, model) %>%
  summarise(predictors = paste(item, collapse = "\n"), .groups = "drop")

ms_tab <- var_cells %>%
  left_join(pred_cells, by = c("species", "model")) %>%
  mutate(predictors = ifelse(is.na(predictors), "None", predictors)) %>%
  pivot_wider(names_from = model, values_from = c(variance, predictors)) %>%
  arrange(match(species, top_spp))
ms_tab <- ms_tab[, c("species", as.vector(rbind(paste0("variance_", names(RANDOM_EFFECTS)),
                                                paste0("predictors_", names(RANDOM_EFFECTS)))))]
ms_tab

write.csv(ms_tab, paste0("./Outputs/IntraspecificShift_TopSpp_ManuscriptTable_", param_tag, ".csv"),
          row.names = FALSE)

ms_ft <- flextable(ms_tab) %>%
  set_header_labels(values = setNames(as.list(c("Species", rep(c("Model variance", "Predictors"),
                                                               length(RANDOM_EFFECTS)))),
                                      names(ms_tab))) %>%
  add_header_row(values = c("", MODEL_LABELS[names(RANDOM_EFFECTS)]),
                 colwidths = c(1, rep(2, length(RANDOM_EFFECTS)))) %>%
  italic(j = "species", part = "body") %>%
  valign(valign = "top", part = "body") %>%
  align(i = 1, part = "header", align = "center") %>%
  theme_booktabs() %>%
  fontsize(size = 9, part = "all") %>%
  set_caption(paste0(
    "Model variance: marginal and conditional R2, with the marginal R2 partitioned into ",
    "variance unique to environment, unique to competition, and shared (Borcard partition; ",
    "shared can be negative when predictor groups suppress each other). Predictors: fixed effects with ",
    "p < ", SIG_LEVEL, " (OLS t-test across sites, Satterthwaite within sites)", "; * p < 0.05, ** p < 0.01, *** p < 0.001.")) %>%
  width(j = 1, width = 1.1) %>%
  width(j = grep("^variance_",   names(ms_tab)), width = 1.5) %>%
  width(j = grep("^predictors_", names(ms_tab)), width = 2.4)
save_as_docx(ms_ft, path = paste0("./Outputs/IntraspecificShift_TopSpp_ManuscriptTable_", param_tag, ".docx"),
             pr_section = officer::prop_section(page_size = officer::page_size(orient = "landscape")))
