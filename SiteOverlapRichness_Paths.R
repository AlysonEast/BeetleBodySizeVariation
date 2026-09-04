##############################################################################
## SiteOverlapRichness_Paths.R
## Site-level path analysis (SEM) of species richness.
## Builds on OverlapRichness.R.
##
## Construct mapping (from OverlapRichness.R "####Paths####" block):
##   Climate      : Tmean (bio01_mean) + Precip (bio12_mean)
##   Productivity : NPP   (added below from NEONSiteNPP.csv)
##   Heterogeneity: Geodiv (an srtm_* surface metric)  <-- CONFIRM WHICH COLUMN
##   Interaction  : ITV   (trait overlap)              <-- overlap_unnorm_obs per instruction
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
geodiv_dir<-"/media/aly/Penobscot/NEON/Geodiversity/edi.2320.1/"

## ============================================================ ##
## 0. CONFIG -- edit these, everything downstream is parameterized
## ============================================================ ##

ITV_COL    <- "sqrt_overlap_unnorm_obs"   
GEODIV_COL <- "log_bio01_sq"              # heterogeneity proxy; alternatives: srtm_sdq, srtm_sq ... SRTM excludes AK sites... 
RICH_COL   <- "richness"
TMEAN_COL  <- "bio01_mean"                # This is second order.. 
PPT_COL    <- "log_bio12_mean"
NPP_COL  <- "log_Npp" 

## Transforms (applied before standardizing)
STANDARDIZE    <- TRUE               # z-score all model vars (coeffs in SD units)
EXCLUDE_ISLANDS <- TRUE

## ============================================================ ##
## 1. Assemble site data: start from siteDF, add NPP
## ============================================================ ##
#Read in and merge overlap and richness data
# site_overlap<-read.csv("./Outputs/site_by_all_noaug_ByYearAvg_IndividualNull.csv") #use site_by_all becuase there are no exclusions due to domains with 1 site
#Read in overlap data
site_2018<-read.csv("./Outputs/site_by_all_aug_2018_IndividualNull.csv")
site_2019<-read.csv("./Outputs/site_by_all_aug_2019_IndividualNull.csv")
head(site_2018)
site_2018$Year<-2018
site_2019$Year<-2019

site_overlap<-rbind(site_2018, site_2019)
site_overlap$latitude<-NULL
site_overlap$Assemblage<-paste0(site_overlap$siteID,"_",site_overlap$Year)

site_richness<-read.csv("../BeetleBiodiversity/site_annual_EstimatedSppRichness.csv")
head(site_richness)
site_richness$X<-NULL
site_richness$Assemblage<-paste0(site_richness$Assemblage,"_",site_richness$Year)
head(site_richness)
siteDF<-merge(site_overlap, site_richness, by = "Assemblage", all.x = TRUE, all.y = FALSE)
head(siteDF)

site_abund2018<-read.csv("./Data/siteTotal_abund_2018.csv")
site_abund2019<-read.csv("./Data/siteTotal_abund_2019.csv")
head(site_abund2018)
site_abund2018$Assemblage<-paste0(site_abund2018$siteID,"_2018")
site_abund2019$Assemblage<-paste0(site_abund2019$siteID,"_2019")
site_abund<-rbind(site_abund2018, site_abund2019)
head(site_abund)

siteDF<-merge(siteDF, site_abund, by="Assemblage")
head(siteDF)

#How stable is overlap from year to year
siteDF2018<-subset(siteDF, Year.x==2018)
siteDF2019<-subset(siteDF, Year.x==2019)
pair<-merge(siteDF2018, siteDF2019, by="siteID.x", all=TRUE)
head(pair)

plot(pair$n_overlap_sp.x~pair$n_overlap_sp.y)
abline(a=0, b=1)

plot(pair$overlap_unnorm_obs.x~pair$overlap_unnorm_obs.y)
abline(a=0, b=1)
plot(sqrt(pair$overlap_unnorm_obs.x)~sqrt(pair$overlap_unnorm_obs.y))
abline(a=0, b=1)

plot(pair$niche_range_obs.x~pair$niche_range_obs.y)
abline(a=0, b=1)

siteDF$richness<-siteDF$Estimator
#What overlap values need to be removed?
plot(siteDF$richness~siteDF$n_overlap_sp)
abline(a=0, b=1)
plot(siteDF$Observed~siteDF$n_overlap_sp)
abline(a=0, b=1)

siteDF$diff<-siteDF$Observed-siteDF$n_overlap_sp
hist(siteDF$diff)

siteDF$diffpct<-((siteDF$Observed-siteDF$n_overlap_sp)/siteDF$Observed)
# siteDF$diffpct<-as.numeric(ifelse(siteDF$diffpct<0, paste0(NA), siteDF$diffpct))

table(siteDF$diffpct, useNA = "ifany")
hist(siteDF$diffpct)
siteDF$diffdouble<-ifelse(siteDF$diffpct>.5, paste0(1), paste0(0))
siteDF$diffthird<-ifelse(siteDF$diffpct>(2/3), paste0(1), paste0(0))


ggplot(siteDF, aes(x=richness, y=n_overlap_sp, colour = log(overlap_norm_obs))) +
  geom_point(alpha=0.5) +
  geom_errorbar(aes(xmin = LCL, xmax=UCL), alpha=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradient(low = "purple", high = "orange")

ggplot(siteDF, aes(x=richness, y=n_overlap_sp, colour = diffpct, shape = diffdouble)) +
  geom_point(alpha=0.5, size=3) +
  geom_errorbar(aes(xmin = LCL, xmax=UCL), alpha=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradient(low = "purple", high = "orange")

table(siteDF$diffdouble)
table(siteDF$diffsig)

#Evaluate validity of richness estimates
ggplot(siteDF, aes(x=richness, y=n_overlap_sp, colour = completeness, shape = diffdouble)) +
  geom_point(alpha=0.5, size=3) +
  geom_errorbar(aes(xmin = LCL, xmax=UCL), alpha=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradient(low = "purple", high = "orange")

hist(siteDF$completeness)

siteDF$poorRichnessEstimate<-ifelse(siteDF$completeness<.5, paste0(1), paste0(0))
table(siteDF$poorRichnessEstimate)
table(siteDF$poorRichnessEstimate, siteDF$diffdouble)
table(siteDF$poorRichnessEstimate, siteDF$diffthird)
table(siteDF$poorRichnessEstimate, siteDF$siteID.x)


ggplot(siteDF, aes(x=richness, y=n_overlap_sp, colour = poorRichnessEstimate, shape = diffthird)) +
  geom_point(alpha=0.5, size=3) +
  geom_errorbar(aes(xmin = LCL, xmax=UCL), alpha=0.5) +
  geom_abline(intercept = 0, slope = 1) 

#### Exclusion ####
preExclusion<-siteDF
if (EXCLUDE_ISLANDS) siteDF <- siteDF %>% 
  filter(!siteID.x %in% c("PUUM","LAJA","GUAN"))
siteDF<-subset(siteDF, completeness>=.5)
siteDF<-subset(siteDF, diffpct<=(2/3))
siteDF<-subset(siteDF, !is.na(overlap_norm_obs))

dim(preExclusion)
dim(siteDF)
dim(preExclusion)[1]-dim(siteDF)[1]

symdiff(levels(as.factor(preExclusion$siteID.x)),levels(as.factor(siteDF$siteID.x)))
dim(table(siteDF$siteID.x))

#Env Variaibles
siteDF$domainID<-NULL
neonDivData::neon_sites
siteDF$siteID<-siteDF$siteID.x
siteDF$siteID.x<-NULL
siteDF$siteID.y<-NULL
siteDF<-merge(neonDivData::neon_sites, siteDF,  by = "siteID")

geodiv_dir<-"/media/aly/Penobscot/NEON/Geodiversity/edi.2320.1/"
geodiv<-read.csv(paste0(geodiv_dir,"NEON_site_footprint_elev30m.csv"))
head(geodiv)
geodiv$domainID<-NULL

siteDF<-merge(siteDF, geodiv, by="siteID")
head(siteDF)

NPP<-read.csv("../NEONSites_MODIS_NPP_2018_2019.csv") #from https://code.earthengine.google.com/aab72ea17a0eb8cecae913e2ee839254
siteDF<-merge(siteDF, NPP[,c("Npp","Gpp","siteID")], by="siteID")
head(siteDF)

#### Pair site#

pairs.panels(siteDF[,c("bio01_mean","bio12_mean","bio01_sq","Npp","overlap_unnorm_obs","richness")])
siteDF$log_bio12_mean<-log10(siteDF$bio12_mean)
siteDF$log_bio01_sq<-log10((siteDF$bio01_sq+0.01))
siteDF$log_Npp<-log10(siteDF$Npp)

siteDF$log_abund<-log10(siteDF$abund)
siteDF$log_overlap_unnorm_obs<-log10(siteDF$overlap_unnorm_obs)
siteDF$sqrt_overlap_unnorm_obs<-sqrt(siteDF$overlap_unnorm_obs)
siteDF$log_Npp<-log10(siteDF$Npp)
siteDF$log_bio1_mean<-log10((siteDF$bio01_mean+0.001))

pairs.panels(siteDF[,c("bio01_mean","log_bio1_mean","bio12_mean","log_bio12_mean",
                       "bio01_sq","log_bio01_sq","Npp","log_Npp",
                       "overlap_unnorm_obs","log_overlap_unnorm_obs","sqrt_overlap_unnorm_obs",
                       "richness")])

pairs.panels(siteDF[,c("bio01_sq","bio12_sq","srtm_sq",
                       "richness")])


## ============================================================ ##
## 2. Build modeling frame: select, rename, transform, complete-case, scale
## ============================================================ ##

dat <- data.frame(
  siteID = siteDF$siteID,
  tmean  = siteDF[[TMEAN_COL]],
  ppt    = siteDF[[PPT_COL]],
  npp    = siteDF[[NPP_COL]],
  geodiv = siteDF[[GEODIV_COL]],
  itv    = siteDF[[ITV_COL]],
  rich   = siteDF[[RICH_COL]]
)

## complete-case across ALL model variables so every candidate model is fit on
## identical rows (required for valid AIC/BIC comparison). With the current
## Geodiv column all 47 sites should be retained -- verify in the printout.
model_vars <- c("tmean", "ppt", "npp", "geodiv", "itv", "rich")
cc <- complete.cases(dat[, model_vars])

cat("\n--- complete-case summary ---\n")
cat("N total sites :", nrow(dat), "\n")
cat("N used (cc)   :", sum(cc), "\n")
cat("Dropped sites :", paste(dat$siteID[!cc], collapse = ", "), "\n\n")

dat <- dat[cc, ]

## standardize (keep raw copy in case you want unscaled effects later)
dat_raw <- dat
if (STANDARDIZE) {
  dat[, model_vars] <- scale(dat[, model_vars])
}

## quadratic temperature term. Built AFTER standardizing, so tmean is already
## mean-centered: squaring it puts the vertex at the mean temperature, minimizes
## collinearity with the linear term, and keeps the marginal slopes below exact.
## Only bio_1 gets this -- precip curvature is handled by its log transform.
if (STANDARDIZE) {
  dat$tmean_sq <- dat$tmean^2
} else {
  dat$tmean_sq <- (dat$tmean - mean(dat$tmean))^2
}
cat("cor(tmean, tmean_sq) =", round(cor(dat$tmean, dat$tmean_sq), 3),
    "  (large |r| => temp is skewed; consider poly(tmean,2))\n")

## ============================================================ ##
## 2.5 Quadratic pretest: does temp curvature occur INSIDE the data range?
## ============================================================ ##
fit_lin  <- lm(rich ~ tmean, data = dat)
fit_quad <- lm(rich ~ tmean + tmean_sq, data = dat)
cat("\n--- quadratic temperature pretest (plain lm) ---\n")
print(anova(fit_lin, fit_quad))                 # F-test on the added quadratic
b1q <- coef(fit_quad)[["tmean"]]; b2q <- coef(fit_quad)[["tmean_sq"]]
vertex <- -b1q / (2 * b2q)
cat(sprintf("curvature = %.3f (%s); vertex at std tmean = %.2f; range = [%.2f, %.2f]\n",
            b2q, ifelse(b2q < 0, "concave/hump", "convex/U"),
            vertex, min(dat$tmean), max(dat$tmean)))
cat(ifelse(vertex >= min(dat$tmean) & vertex <= max(dat$tmean),
           ">> turnover is inside the data: quadratic supported.\n",
           ">> vertex outside the data: you sample one limb only -- a monotone fit\n   is more honest; reconsider the quadratic before adding it to the SEM.\n"))

## ============================================================ ##
## 3. FULL path model: fit + direct / indirect / total effects on richness
## ============================================================ ##
## The FULL model = every theory operating at once. Reading the richness paths:
##   rich ~ tmean   ambient-energy / kinetic (MTE): temperature acts directly
##   rich ~ ppt     water-energy: precipitation acts directly
##   rich ~ npp     species-energy ("more individuals")
##   rich ~ geodiv  habitat heterogeneity
##   rich ~ itv     niche packing (focal mechanism)
## Backbone: climate -> NPP (structural); environment -> ITV (trait space).
## The indirect/total blocks split each driver into its direct effect vs the
## parts routed through productivity (NPP) and through trait space (ITV) --
## i.e. they partition each variable's action across the competing theories.

m1_full <- '
  npp  ~ a1*tmean + a2*ppt
  itv  ~ b1*tmean + b2*ppt + b3*npp + b4*geodiv
  rich ~ c1*tmean + q1*tmean_sq + c2*ppt + c3*npp + c4*geodiv + d*itv

  ind_tmean_itv     := b1*d
  ind_ppt_itv       := b2*d
  ind_npp_itv       := b3*d
  ind_geodiv_itv    := b4*d
  ind_tmean_npp     := a1*c3
  ind_ppt_npp       := a2*c3
  ind_tmean_npp_itv := a1*b3*d
  ind_ppt_npp_itv   := a2*b3*d

  # temperature curvature -- the LDG test (expect q1 < 0: thermal optimum)
  curv_tmean       := q1
  # marginal dRich/dTmean at cold / mean / warm sites (std temp = -1, 0, +1)
  slope_tmean_cold := c1 + 2*q1*(-1)
  slope_tmean_mean := c1
  slope_tmean_warm := c1 + 2*q1*(1)

  # LINEAR-component totals; temperature total is level-dependent (see slopes)
  tot_tmean  := c1 + b1*d + a1*c3 + a1*b3*d
  tot_ppt    := c2 + b2*d + a2*c3 + a2*b3*d
  tot_npp    := c3 + b3*d
  tot_geodiv := c4 + b4*d
'

## ============================================================ ##
## 4. Candidate models = competing theories of the latitudinal gradient
## ============================================================ ##
## SHARED BACKBONE (identical in every model, so all candidates share the same
## six variables and N -> AIC/BIC valid across the WHOLE set):
##   npp ~ tmean + ppt            climate drives productivity (structural)
##   itv ~ tmean+ppt+npp+geodiv   environment shapes trait space
## Models differ ONLY in the RICHNESS equation: which direct-to-richness paths
## are free vs fixed to zero. Each choice IS a theory. The reduced models are
## over-identified (df > 0), so CFI/RMSEA/chisq are informative again -- a
## good-fitting reduced model means the omitted direct paths were not needed.

## richness equation per theory (ITV retained except in model 8)
m2_dropGeo <- '  
  npp  ~ a1*tmean + a2*ppt
  itv  ~ b1*tmean + b2*ppt + b3*npp 
  rich ~ c1*tmean + q1*tmean_sq + c2*ppt + c3*npp + d*itv

  # indirect paths to richness
  ind_tmean_itv     := b1*d
  ind_ppt_itv       := b2*d
  ind_npp_itv       := b3*d
  ind_tmean_npp     := a1*c3
  ind_ppt_npp       := a2*c3
  ind_tmean_npp_itv := a1*b3*d
  ind_ppt_npp_itv   := a2*b3*d

  # temperature curvature -- LDG test (expect q1 < 0: thermal optimum)
  curv_tmean       := q1
  slope_tmean_cold := c1 + 2*q1*(-1)
  slope_tmean_mean := c1
  slope_tmean_warm := c1 + 2*q1*(1)

  # LINEAR-component totals; temperature total is level-dependent (see slopes)
  tot_tmean  := c1 + b1*d + a1*c3 + a1*b3*d
  tot_ppt    := c2 + b2*d + a2*c3 + a2*b3*d
  tot_npp    := c3 + b3*d
'

m3_dropNPP <- '  
  itv  ~ b1*tmean + b2*ppt + b4*geodiv
  rich ~ c1*tmean + q1*tmean_sq + c2*ppt + c4*geodiv + d*itv

  # indirect paths to richness
  ind_tmean_itv     := b1*d
  ind_ppt_itv       := b2*d
  ind_geodiv_itv    := b4*d

  # temperature curvature -- LDG test (expect q1 < 0: thermal optimum)
  curv_tmean       := q1
  slope_tmean_cold := c1 + 2*q1*(-1)
  slope_tmean_mean := c1
  slope_tmean_warm := c1 + 2*q1*(1)

  # LINEAR-component totals; temperature total is level-dependent (see slopes)
  tot_tmean  := c1 + b1*d 
  tot_ppt    := c2 + b2*d
  tot_geodiv := c4 + b4*d
'

m4_climITV <- '
  itv  ~ b1*tmean + b2*ppt 
  rich ~ c1*tmean + q1*tmean_sq + c2*ppt + d*itv

  # indirect paths to richness
  ind_tmean_itv     := b1*d
  ind_ppt_itv       := b2*d

  # temperature curvature -- LDG test (expect q1 < 0: thermal optimum)
  curv_tmean       := q1
  slope_tmean_cold := c1 + 2*q1*(-1)
  slope_tmean_mean := c1
  slope_tmean_warm := c1 + 2*q1*(1)

  # LINEAR-component totals; temperature total is level-dependent (see slopes)
  tot_tmean  := c1 + b1*d
  tot_ppt    := c2 + b2*d
'

m5_DropClim <- '
  itv  ~ b3*npp + b4*geodiv
  rich ~ c3*npp + c4*geodiv + d*itv

  # indirect paths to richness
  ind_npp_itv       := b3*d
  ind_geodiv_itv    := b4*d
  
  # total effects on richness
  tot_npp    := c3 + b3*d
  tot_geodiv := c4 + b4*d
'
m6_NPP <- '
  itv  ~ b3*npp
  rich ~ c3*npp + d*itv

  # indirect paths to richness
  ind_npp_itv       := b3*d

  # total effects on richness
  tot_npp    := c3 + b3*d
'
m7_Geo <- '
  itv  ~ b4*geodiv
  rich ~ c4*geodiv + d*itv

  # indirect paths to richness
  ind_geodiv_itv    := b4*d
  
  # total effects on richness
  tot_geodiv := c4 + b4*d
'
models <- list(
  "1_Full"              = m1_full,
  "2_dropGeo"           = m2_dropGeo,
  "3_dropNPP"           = m3_dropNPP,
  "4_climITV"           = m4_climITV,
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
## 6. visualize models
## ============================================================ ##
library(lavaanPlot)
lavaanPlot(model = fits$`1_Full`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

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

lavaanPlot(model = fits$`4_climITV`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`6_Npp`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

library(tidySEM)

# Create a default graph from the fitted model
graph_sem(fits$`3_dropNPP`)

ggplot(dat, aes(x=rich, y=itv)) +
  geom_point() +
  theme_pubr()

ggplot(dat, aes(x=itv, y=geodiv)) +
  geom_point() +
  theme_pubr()

ggplot(dat, aes(x=rich, y=geodiv)) +
  geom_point() +
  theme_pubr()

## ============================================================ ##
## 6. visualize TOP model
## ============================================================ ##

library(lavaan); library(ggplot2); library(ggpubr); library(dplyr)

USE_BOOT <- TRUE
N_BOOT   <- 2000   # bump to 5000 for the final figure

fit3 <- sem(m3_dropNPP, data = dat, estimator = "ML",
            se = if (USE_BOOT) "bootstrap" else "standard",
            bootstrap = N_BOOT, iseed = 42)

## R^2 for the two responses (how much of ITV and richness the model explains)
cat("\n--- R^2 (endogenous) ---\n"); print(round(lavInspect(fit3, "rsquare"), 3))

## all standardized paths + defined effects, with CIs
pe <- parameterEstimates(fit3, standardized = TRUE, ci = TRUE)

## structural paths (b* = env->ITV, c*/q1 = ->richness, d = ITV->richness)
paths <- subset(pe, op == "~",
                c("lhs","rhs","label","est","ci.lower","ci.upper","pvalue","std.all"))
cat("\n--- structural paths (std.all = fully standardized) ---\n")
print(paths, row.names = FALSE, digits = 3)

## effect decomposition on richness (from the := lines in m3)
eff <- subset(pe, op == ":=",
              c("label","est","ci.lower","ci.upper","pvalue"))
cat("\n--- effects on richness: direct via ITV (ind_*), totals (tot_*), curvature ---\n")
print(eff, row.names = FALSE, digits = 3)

## relative contribution ranking: |standardized total effect| on richness.
## Temperature is split: linear-route total + curvature (its total is
## level-dependent, so read curv_tmean and the slopes alongside).
rank_tbl <- data.frame(
  driver = c("temperature (linear route)", "temperature (curvature)",
             "precip", "heterogeneity", "ITV (direct)"),
  effect = c(eff$est[eff$label=="tot_tmean"],  eff$est[eff$label=="curv_tmean"],
             eff$est[eff$label=="tot_ppt"],    eff$est[eff$label=="tot_geodiv"],
             paths$std.all[paths$label=="d"])
)
rank_tbl <- rank_tbl[order(-abs(rank_tbl$effect)), ]
cat("\n--- relative contribution (|standardized effect on richness|) ---\n")
print(rank_tbl, row.names = FALSE, digits = 3)

library(lavaanPlot)
lavaanPlot(model = fit3, coefs = TRUE, stand = TRUE, sig = 0.05,
           stars = c("regress"),
           graph_options = list(rankdir = "LR"))

## tidySEM alternative with an explicit layout (tmean_sq sits beside tmean)
library(tidySEM)
lay <- get_layout(
  "tmean", "tmean_sq", "ppt", "geodiv",
  NA,      "itv",      NA,    NA,
  NA,      "rich",     NA,    NA,
  rows = 3)
graph_sem(fit3, layout = lay)

gb  <- function(l) pe$est[pe$label == l]           # grab a labeled coef
c1<-gb("c1"); q1<-gb("q1"); c2<-gb("c2"); c4<-gb("c4"); d<-gb("d")
b1<-gb("b1"); b2<-gb("b2"); b4<-gb("b4")
mu <- function(v) mean(dat_raw[[v]]); sdv <- function(v) sd(dat_raw[[v]])

## (i) temperature -> richness, with the hump
tz <- seq(min(dat$tmean), max(dat$tmean), length.out = 250)
grid <- data.frame(
  tmean = tz * sdv("tmean") + mu("tmean"),
  total  = ((c1 + d*b1)*tz + q1*tz^2) * sdv("rich") + mu("rich"),
  direct = (c1*tz + q1*tz^2)          * sdv("rich") + mu("rich"))
vz  <- -(c1 + d*b1) / (2*q1)                        # total-effect vertex (std)
vx  <- vz * sdv("tmean") + mu("tmean")
pts <- data.frame(tmean = dat_raw$tmean, rich = dat_raw$rich)

p_temp <- ggplot() +
  geom_point(data = pts, aes(tmean, rich), alpha = .55) +
  geom_line(data = grid, aes(tmean, total),  linewidth = 1.1, colour = "#c1440e") +
  geom_line(data = grid, aes(tmean, direct), linewidth = .8, linetype = 2, colour = "grey45") +
  {if (vz >= min(tz) & vz <= max(tz)) geom_vline(xintercept = vx, linetype = 3)} +
  labs(x = "Mean annual temperature (bio_1)", y = "Estimated richness",
       title = "Temperature–richness (solid = total, dashed = direct)") +
  theme_pubr()

## ============================================================ ##
## 8.3 Model-implied trends: TOTAL vs DIRECT for each env predictor
##     total  = direct-to-richness + the part routed through ITV (d * b_k)
##     direct = the coefficient(s) straight into the richness equation
## In m3 ITV is the only mediator, so total = direct + d*b_k exactly.
## Each panel holds the other predictors at their means, so both lines
## pivot on the bivariate mean; the gap between them IS the mediated share.
## ============================================================ ##
gb  <- function(l) pe$est[pe$label == l]
c1<-gb("c1"); q1<-gb("q1"); c2<-gb("c2"); c4<-gb("c4"); d<-gb("d")
b1<-gb("b1"); b2<-gb("b2"); b4<-gb("b4")
mu  <- function(v) mean(dat_raw[[v]]); sdv <- function(v) sd(dat_raw[[v]])

trend_panel <- function(v, direct_slope, total_slope, xlab,
                        quad = 0, mark_vertex = FALSE) {
  z  <- seq(min(dat[[v]]), max(dat[[v]]), length.out = 250)
  bt <- function(slope) (slope*z + quad*z^2) * sdv("rich") + mu("rich")
  df <- rbind(
    data.frame(x = z*sdv(v)+mu(v), rich = bt(total_slope),  path = "total"),
    data.frame(x = z*sdv(v)+mu(v), rich = bt(direct_slope), path = "direct"))
  p <- ggplot() +
    geom_point(data = data.frame(x = dat_raw[[v]], rich = dat_raw$rich),
               aes(x, rich), alpha = .5, colour = "grey40") +
    geom_line(data = df, aes(x, rich, colour = path, linetype = path),
              linewidth = 1) +
    scale_colour_manual(values = c(total = "#c1440e", direct = "grey35")) +
    scale_linetype_manual(values = c(total = 1, direct = 2)) +
    labs(x = xlab, y = "Estimated richness", colour = NULL, linetype = NULL) +
    theme_pubr()
  if (mark_vertex && quad != 0) {                 # only meaningful for temp
    vz <- -total_slope / (2*quad)
    if (vz >= min(z) & vz <= max(z))
      p <- p + geom_vline(xintercept = vz*sdv(v)+mu(v),
                          linetype = 3, colour = "grey60")
  }
  p
}

p_temp   <- trend_panel("tmean",  direct_slope = c1,
                        total_slope = c1 + d*b1,
                        xlab = "Mean annual temp (bio_1)",
                        quad = q1, mark_vertex = TRUE)
p_ppt    <- trend_panel("ppt",    direct_slope = c2,
                        total_slope = c2 + d*b2,
                        xlab = "log precip (bio_12)")
p_geodiv <- trend_panel("geodiv", direct_slope = c4,
                        total_slope = c4 + d*b4,
                        xlab = "log thermal heterogeneity")

ggarrange(p_temp, p_ppt, p_geodiv, ncol = 3,
          common.legend = TRUE, legend = "bottom", labels = "AUTO")

## (ii) ITV -> richness (the focal mechanism, slope d)
iz <- seq(min(dat$itv), max(dat$itv), length.out = 100)
p_itv <- ggplot() +
  geom_point(data = data.frame(itv = dat_raw$itv, rich = dat_raw$rich),
             aes(itv, rich), alpha = .55) +
  geom_line(data = data.frame(itv = iz*sdv("itv")+mu("itv"),
                              rich = (d*iz)*sdv("rich")+mu("rich")),
            aes(itv, rich), linewidth = 1.1, colour = "#1f6f6f") +
  labs(x = "ITV (sqrt overlap)", y = "Estimated richness",
       title = "Trait overlap -> richness") + theme_pubr()

## (iii) environment -> ITV (the mediator's drivers: b1, b2, b4)
env_panel <- function(v, coef, xlab) {
  z <- seq(min(dat[[v]]), max(dat[[v]]), length.out = 100)
  ggplot() +
    geom_point(data = data.frame(x = dat_raw[[v]], itv = dat_raw$itv),
               aes(x, itv), alpha = .55) +
    geom_line(data = data.frame(x = z*sdv(v)+mu(v),
                                itv = (coef*z)*sdv("itv")+mu("itv")),
              aes(x, itv), linewidth = 1, colour = "#555599") +
    labs(x = xlab, y = "ITV (sqrt overlap)") + theme_pubr()
}
p_e1 <- env_panel("tmean",  b1, "Temperature (bio_1)")
p_e2 <- env_panel("ppt",    b2, "log precip (bio_12)")
p_e3 <- env_panel("geodiv", b4, "log thermal heterogeneity")

ggarrange(p_temp, p_itv, p_e1, p_e2, p_e3, ncol = 2, nrow = 3, labels = "AUTO")

fp <- subset(pe, label %in% c("tot_tmean","curv_tmean","slope_tmean_cold",
                              "slope_tmean_warm","tot_ppt","tot_geodiv",
                              "ind_tmean_itv","ind_ppt_itv","ind_geodiv_itv","d"),
             c("label","est","ci.lower","ci.upper"))
fp$label <- factor(fp$label, levels = rev(fp$label))

ggplot(fp, aes(est, label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = ci.lower, xmax = ci.upper)) +
  labs(x = "Standardized effect on richness (bootstrap CI)", y = NULL,
       title = "m3_dropNPP: effect decomposition") + theme_pubr()
