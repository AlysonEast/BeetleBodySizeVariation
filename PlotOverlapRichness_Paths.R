##############################################################################
## plotOverlapRichness_Paths.R
## plot-level path analysis (SEM) of species richness.
## Builds on OverlapRichness.R.
##
## Construct mapping (from OverlapRichness.R "####Paths####" block):
##   Climate      : Tmean (bio01_mean) + Precip (bio12_mean)
##   Productivity : NPP   (added below from NEONplotNPP.csv)
##   Heterogeneity: Complexity (an srtm_* surface metric)  <-- CONFIRM WHICH COLUMN
##   Interaction  : Overlap   (trait overlap)              <-- overlap_unnorm_obs per instruction
##   Response     : Richness (median_richness)
##
## Full path model + 7 candidate models for selection (from sketch).
##############################################################################
library(ggplot2)
library(ggpubr)
library(neonDivData)
library(lavaan)
library(dplyr)
library(psych)

setwd("/home/aly/Beetles/BeetleBodySizeVariation")

## ============================================================ ##
## 0. CONFIG -- edit these, everything downstream is parameterized
## ============================================================ ##

#bio_1","bio_12","rugosity_RC","Npp","overlap_unnorm_obs
Overlap_COL    <- "overlap_norm_obs"   # NOTE: project default elsewhere is overlap_norm_obs
Complexity_COL <- "log_rugosity_RC"              # heterogeneity proxy; alternatives: srtm_sdq, srtm_sq ... SRTM excludes AK plots... 
RICH_COL   <- "median_richness"
TMEAN_COL  <- "bio_11" #Mean daily mean temperature of coldest quarter
PPT_COL    <- "log_bio_17" #Mean monthly precipitation of the driest quarter
NPP_COL  <- "log_Npp"                    

## Transforms (applied before standardizing)
LOG_Overlap        <- TRUE               # overlap spans many orders of magnitude -> log
RICH_TRANSFORM <- "log"             # "none", "sqrt", or "log"
STANDARDIZE    <- TRUE               # z-score all model vars (coeffs in SD units)

## ============================================================ ##
## 1. Assemble plot data: start from plotDF, add NPP
## ============================================================ ##
#Read in and merge overlap and richness data
plot_overlap<-read.csv("./Outputs/plot_by_all_IndividualNull.csv") #use plot_by_all becuase there are no exclusions due to domains with 1 site
plot_overlap$latitude<-NULL
plot_richness<-read.csv("../BeetleBiodiversity/plot_annualVarWeightedMean_EstimatedSppRichness.csv")
plot_richness$X<-NULL
plotDF<-merge(plot_overlap, plot_richness, by.x = "plotID", by.y = "PlotID")
head(plotDF)

struc<-read.csv("./Outputs/BETplot_Rugosity.csv")
struc$X<-NULL
env<-read.csv("./Outputs/BeetlePlotswEnvData.csv")
NPP<-read.csv("../NEON_MODIS_NPP_2018_2019.csv") #from https://code.earthengine.google.com/b41a55076352b2d9e21ac5e74bf337bc

plotDF<-merge(plotDF, struc, by="plotID")
plotDF<-merge(plotDF, env, by="plotID")
plotDF<-merge(plotDF, NPP[,c("Npp","Gpp","plotID")], by="plotID")
head(plotDF)

pairs.panels(plotDF[,c("bio_11","bio_17","rugosity_RC","Npp","overlap_norm_obs","median_richness")])
plotDF$log_rugosity_RC<-log10((plotDF$rugosity_RC+0.01))
plotDF$log_overlap_unnorm_obs<-log10(plotDF$overlap_norm_obs)
plotDF$log_median_richness<-log10(plotDF$median_richness)
plotDF$sqrt_median_richness<-sqrt(plotDF$median_richness)
plotDF$log_Npp<-log10(plotDF$Npp)
plotDF$log_bio_17<-log10(plotDF$bio_17)
pairs.panels(plotDF[,c("bio_11","log_bio_17","log_rugosity_RC","log_Npp","log_overlap_unnorm_obs","sqrt_median_richness","log_median_richness")])

## ============================================================ ##
## 2. Build modeling frame: select, rename, transform, complete-case, scale
## ============================================================ ##

dat <- data.frame(
  plotID = plotDF$plotID,
  tmean  = plotDF[[TMEAN_COL]],
  ppt    = plotDF[[PPT_COL]],
  npp    = plotDF[[NPP_COL]],
  Complexity = plotDF[[Complexity_COL]],
  Overlap    = plotDF[[Overlap_COL]],
  rich   = plotDF[[RICH_COL]]
)

## transforms
if (LOG_Overlap) {
  if (any(dat$Overlap <= 0, na.rm = TRUE))
    stop("Non-positive Overlap values present; choose log1p or a small constant before logging.")
  dat$Overlap <- log10(dat$Overlap)
}
if (RICH_TRANSFORM == "sqrt") dat$rich <- sqrt(dat$rich)
if (RICH_TRANSFORM == "log")  dat$rich <- log(dat$rich)

## complete-case across ALL model variables so every candidate model is fit on
## identical rows (required for valid AIC/BIC comparison). With the current
## Complexity column all 47 sites should be retained -- verify in the printout.
model_vars <- c("tmean", "ppt", "npp", "Complexity", "Overlap", "rich")
cc <- complete.cases(dat[, model_vars])

cat("\n--- complete-case summary ---\n")
cat("N total sites :", nrow(dat), "\n")
cat("N used (cc)   :", sum(cc), "\n")
cat("Dropped sites :", paste(dat$plotID[!cc], collapse = ", "), "\n\n")

dat <- dat[cc, ]

## standardize (keep raw copy in case you want unscaled effects later)
dat_raw <- dat
if (STANDARDIZE) {
  dat[, model_vars] <- scale(dat[, model_vars])
}

## ============================================================ ##
## 3. FULL path model: fit + direct / indirect / total effects on richness
## ============================================================ ##
## The FULL model = every theory operating at once. Reading the richness paths:
##   rich ~ tmean   ambient-energy / kinetic (MTE): temperature acts directly
##   rich ~ ppt     water-energy: precipitation acts directly
##   rich ~ npp     species-energy ("more individuals")
##   rich ~ Complexity  habitat heterogeneity
##   rich ~ Overlap     niche packing (focal mechanism)
## Backbone: climate -> NPP (structural); environment -> Overlap (trait space).
## The indirect/total blocks split each driver into its direct effect vs the
## parts routed through productivity (NPP) and through trait space (Overlap) --
## i.e. they partition each variable's action across the competing theories.

m1_full <- '
  npp  ~ a1*tmean + a2*ppt
  Overlap  ~ b1*tmean + b2*ppt + b3*npp + b4*Complexity
  rich ~ c1*tmean + c2*ppt + c3*npp + c4*Complexity + d*Overlap

  # indirect paths to richness
  ind_tmean_Overlap     := b1*d
  ind_ppt_Overlap       := b2*d
  ind_npp_Overlap       := b3*d
  ind_Complexity_Overlap    := b4*d
  ind_tmean_npp     := a1*c3
  ind_ppt_npp       := a2*c3
  ind_tmean_npp_Overlap := a1*b3*d
  ind_ppt_npp_Overlap   := a2*b3*d

  # total effects on richness
  tot_tmean  := c1 + b1*d + a1*c3 + a1*b3*d
  tot_ppt    := c2 + b2*d + a2*c3 + a2*b3*d
  tot_npp    := c3 + b3*d
  tot_Complexity := c4 + b4*d
'
## ============================================================ ##
## 4. Candidate models = competing theories of the latitudinal gradient
## ============================================================ ##
## SHARED BACKBONE (identical in every model, so all candidates share the same
## six variables and N -> AIC/BIC valid across the WHOLE set):
##   npp ~ tmean + ppt            climate drives productivity (structural)
##   Overlap ~ tmean+ppt+npp+Complexity   environment shapes trait space
## Models differ ONLY in the RICHNESS equation: which direct-to-richness paths
## are free vs fixed to zero. Each choice IS a theory. The reduced models are
## over-identified (df > 0), so CFI/RMSEA/chisq are informative again -- a
## good-fitting reduced model means the omitted direct paths were not needed.

## richness equation per theory (Overlap retained except in model 8)
m2_dropGeo <- '  
  npp  ~ a1*tmean + a2*ppt
  Overlap  ~ b1*tmean + b2*ppt + b3*npp 
  rich ~ c1*tmean + c2*ppt + c3*npp + d*Overlap

  # indirect paths to richness
  ind_tmean_Overlap     := b1*d
  ind_ppt_Overlap       := b2*d
  ind_npp_Overlap       := b3*d
  ind_tmean_npp     := a1*c3
  ind_ppt_npp       := a2*c3
  ind_tmean_npp_Overlap := a1*b3*d
  ind_ppt_npp_Overlap   := a2*b3*d

  # total effects on richness
  tot_tmean  := c1 + b1*d + a1*c3 + a1*b3*d
  tot_ppt    := c2 + b2*d + a2*c3 + a2*b3*d
  tot_npp    := c3 + b3*d
'
m3_dropNPP <- '  
  Overlap  ~ b1*tmean + b2*ppt + b4*Complexity
  rich ~ c1*tmean + c2*ppt + c4*Complexity + d*Overlap

  # indirect paths to richness
  ind_tmean_Overlap     := b1*d
  ind_ppt_Overlap       := b2*d
  ind_Complexity_Overlap    := b4*d

  # total effects on richness
  tot_tmean  := c1 + b1*d 
  tot_ppt    := c2 + b2*d
  tot_Complexity := c4 + b4*d
'
m4_climOverlap <- '
  Overlap  ~ b1*tmean + b2*ppt 
  rich ~ c1*tmean + c2*ppt + d*Overlap

  # indirect paths to richness
  ind_tmean_Overlap     := b1*d
  ind_ppt_Overlap       := b2*d

  # total effects on richness
  tot_tmean  := c1 + b1*d
  tot_ppt    := c2 + b2*d
'
m5_DropClim <- '
  Overlap  ~ b3*npp + b4*Complexity
  rich ~ c3*npp + c4*Complexity + d*Overlap

  # indirect paths to richness
  ind_npp_Overlap       := b3*d
  ind_Complexity_Overlap    := b4*d
  
  # total effects on richness
  tot_npp    := c3 + b3*d
  tot_Complexity := c4 + b4*d
'
m6_NPP <- '
  Overlap  ~ b3*npp
  rich ~ c3*npp + d*Overlap

  # indirect paths to richness
  ind_npp_Overlap       := b3*d

  # total effects on richness
  tot_npp    := c3 + b3*d
'
m7_Geo <- '
  Overlap  ~ b4*Complexity
  rich ~ c4*Complexity + d*Overlap

  # indirect paths to richness
  ind_Complexity_Overlap    := b4*d
  
  # total effects on richness
  tot_Complexity := c4 + b4*d
'
models <- list(
  "1_Full"              = m1_full,
  "2_dropGeo"           = m2_dropGeo,
  "3_dropNPP"           = m3_dropNPP,
  "4_climOverlap"       = m4_climOverlap,
  "5_DropClimate"       = m5_DropClim,
  "6_Npp"               = m6_NPP,
  "7_Geo"               = m7_Geo
)

fits <- lapply(models, function(spec) sem(spec, data = dat, estimator = "ML"))

## fit stats + response set (a guard: 'endog' should be identical for every row)
get_fit <- function(fit) {
  m <- fitMeasures(fit, c("npar", "df", "chisq", "pvalue",
                          "cfi", "rmsea", "aic", "bic"))
  data.frame(as.list(round(m, 3)),
             endog = paste(sort(lavNames(fit, "ov.y")), collapse = "+"))
}

fit_table <- do.call(rbind, lapply(fits, get_fit))
fit_table$model <- names(models)
fit_table$dAIC  <- round(fit_table$aic - min(fit_table$aic), 2)
fit_table$wAIC  <- round(exp(-0.5 * fit_table$dAIC) /
                           sum(exp(-0.5 * fit_table$dAIC)), 3)
fit_table <- fit_table[order(fit_table$aic),
                       c("model", "npar", "df", "chisq", "pvalue",
                         "cfi", "rmsea", "aic", "dAIC", "wAIC", "endog")]
cat("\n--- competing-theory path models (AIC valid across all; 'endog' must match) ---\n")
print(fit_table, row.names = FALSE)

## ============================================================ ##
## 5. Targeted nested tests (the two questions that matter)
## ============================================================ ##

## ============================================================ ##
## 6. (optional) visualize best/full model
## ============================================================ ##
library(lavaanPlot)

lavaanPlot(model = fits$`3_dropNPP`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`5_DropClimate`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`4_climOverlap`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

library(tidySEM)

# Create a default graph from the fitted model
graph_data <- graph_sem(fits$`3_dropNPP`)

# Plot the graph
plot(graph_data)


(ggplot(dat, aes(y=Complexity, x=Overlap)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm") +
  labs(y = "Complexity",
    x = "Overlap",
    title = "Complexity → Overlap") +
  theme_pubr() |
    ggplot(dat, aes(y=Complexity, x=rich)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm") +
    labs(y = "Complexity",
         x = "Richness",
         title = "Complexity → Richness") +
    theme_pubr() |
    ggplot(dat, aes(y=Overlap, x=rich)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm") +
    labs(x = "Richness",
         y = "Overlap",
         title = "Overlap → Richness") +
    theme_pubr() |
    ggplot(dat, aes(y=npp, x=Overlap)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm") +
    labs(x = "Overlap",
         y = "NPP",
         title = "NPP → Overlap") +
    theme_pubr() |
    ggplot(dat, aes(y=npp, x=rich)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm") +
    labs(x = "Richness",
         y = "NPP",
         title = "NPP → Richness") +
    theme_pubr()
)

lavaanPlot(model = fits$`3_dropNPP`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`7_Geo`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions


# Standardized coefficients
std <- standardizedSolution(fits$`5_DropClimate`)

std %>%
  filter(op == "~") %>%
  select(lhs, rhs, est.std, pvalue, ci.lower, ci.upper)
dat_plot <- dat %>%
  mutate(
    npp_resid = residuals(lm(npp ~ Complexity + Overlap, data = .)),
    rich_resid = residuals(lm(rich ~ Complexity + Overlap, data = .))
  )

ggplot(dat_plot, aes(npp_resid, rich_resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "NPP (adjusted for Complexityersity and Overlap)",
    y = "Richness (adjusted for Complexityersity and Overlap)"
  ) +
  theme_classic()

m5_DropClim_fit <- sem(
  m5_DropClim,
  data = dat,
  estimator = "ML")
library(dplyr)

sem_coefs <- standardizedSolution(m5_DropClim_fit) %>%
  filter(op == "~") %>%
  select(lhs, rhs, est.std, pvalue, ci.lower, ci.upper)

sem_coefs

library(ggplot2)

p1 <- ggplot(dat, aes(npp, Overlap)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "NPP",
    y = "Overlap",
    title = "NPP → Overlap"
  ) +
  theme_classic()

p2 <- ggplot(dat, aes(Complexity, Overlap)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Complexityersity",
    y = "Overlap",
    title = "Complexityersity → Overlap"
  ) +
  theme_classic()

p3 <- ggplot(dat, aes(npp, rich)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "NPP",
    y = "Species richness",
    title = "NPP → Richness"
  ) +
  theme_classic()

p4 <- ggplot(dat, aes(Complexity, rich)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Complexityersity",
    y = "Species richness",
    title = "Complexityersity → Richness"
  ) +
  theme_classic()

p5 <- ggplot(dat, aes(Overlap, rich)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Overlap",
    y = "Species richness",
    title = "Overlap → Richness"
  ) +
  theme_classic()

library(patchwork)

(p1 | p2) /
  (p3 | p4) /
  p5
