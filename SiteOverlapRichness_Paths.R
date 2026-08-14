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

ITV_COL    <- "overlap_norm_obs"   # NOTE: project default elsewhere is overlap_norm_obs
GEODIV_COL <- "log_bio01_sq"              # heterogeneity proxy; alternatives: srtm_sdq, srtm_sq ... SRTM excludes AK sites... 
RICH_COL   <- "median_richness"
TMEAN_COL  <- "bio01_mean"
PPT_COL    <- "log_bio12_mean"
NPP_COL  <- "Npp" 

## Transforms (applied before standardizing)
LOG_ITV        <- TRUE               # overlap spans many orders of magnitude -> log
RICH_TRANSFORM <- "none"             # "none", "sqrt", or "log"
STANDARDIZE    <- TRUE               # z-score all model vars (coeffs in SD units)

EXCLUDE_ISLANDS <- TRUE

## ============================================================ ##
## 1. Assemble site data: start from siteDF, add NPP
## ============================================================ ##
#Read in and merge overlap and richness data
site_overlap<-read.csv("./Outputs/site_by_all_IndividualNull.csv") #use site_by_all because there are no exclusions due to domains with 1 site
site_overlap$latitude<-NULL
site_richness<-read.csv("../BeetleBiodiversity/Site_annualVarWeightedMean_EstimatedSppRichness.csv")
site_richness$X<-NULL
siteDF<-merge(site_overlap, site_richness, by.x = "siteID", by.y = "Assemblage")
head(siteDF)
siteDF$domainID<-NULL
neonDivData::neon_sites
siteDF<-merge(neonDivData::neon_sites, siteDF,  by = "siteID")

geodiv<-read.csv(paste0(geodiv_dir,"NEON_site_footprint_elev30m.csv"))
head(geodiv)
geodiv$domainID<-NULL

siteDF<-merge(siteDF, geodiv, by="siteID")
head(siteDF)

NPP<-read.csv("../NEONSites_MODIS_NPP_2018_2019.csv") #from https://code.earthengine.google.com/aab72ea17a0eb8cecae913e2ee839254
siteDF<-merge(siteDF, NPP[,c("Npp","Gpp","siteID")], by="siteID")
head(siteDF)

plot(siteDF$median_richness~siteDF$n_overlap_sp)
abline(a=0, b=1)
siteDF$diff<-siteDF$median_richness-siteDF$n_overlap_sp
hist(siteDF$diff)
siteDF$diffpct<-(siteDF$median_richness-siteDF$n_overlap_sp)/siteDF$median_richness
siteDF$diffpct<-as.numeric(ifelse(siteDF$diffpct<0, paste0(NA), siteDF$diffpct))
table(siteDF$diffpct, useNA = "ifany")
hist(siteDF$diffpct)
15/nrow(siteDF)
siteDF$diffsig<-ifelse(siteDF$n_overlap_sp<siteDF$median_lcl | siteDF$n_overlap_sp>siteDF$median_ucl, paste0(1), paste0(0))
siteDF$diffdouble<-ifelse(siteDF$diffpct>.5, paste0(1), paste0(0))

ggplot(siteDF, aes(x=median_richness, y=n_overlap_sp, colour = log(overlap_norm_obs))) +
  geom_point(alpha=0.5) +
  geom_errorbar(aes(xmin = median_lcl, xmax=median_ucl), alpha=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradient(low = "purple", high = "orange")

ggplot(siteDF, aes(x=median_richness, y=n_overlap_sp, colour = diffpct)) +
  geom_point(alpha=0.5) +
  geom_errorbar(aes(xmin = median_lcl, xmax=median_ucl), alpha=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradient(low = "purple", high = "orange")

(ggplot(siteDF, aes(x= median_richness, y=n_overlap_sp, colour = diffdouble)) +
    geom_point(alpha=0.5) +
    geom_errorbar(aes(xmin = median_lcl, xmax=median_ucl), alpha=0.3) +
    geom_abline(intercept = 0, slope = 1) |
    ggplot(siteDF, aes(x=median_richness, y=n_overlap_sp, colour = diffsig)) +
    geom_point() +
    geom_errorbar(aes(xmin = median_lcl, xmax=median_ucl), alpha=0.3) +
    geom_abline(intercept = 0, slope = 1) )

table(siteDF$diffdouble)
table(siteDF$diffsig)

#siteDF<-subset(siteDF, diffdouble==0)


if (EXCLUDE_ISLANDS) dat <- dat %>% 
  filter(!siteID %in% c("PUUM","LAJA","GUAN"))


pairs.panels(siteDF[,c("Latitude","bio01_mean","bio12_mean","bio01_sq","Npp","overlap_unnorm_obs","median_richness")])
siteDF$log_bio12_mean<-log10(siteDF$bio12_mean)
siteDF$log_bio01_sq<-log10((siteDF$bio01_sq+0.01))
siteDF$log_Npp<-log10(siteDF$Npp)
siteDF$log_median_richness<-log10(siteDF$median_richness)
pairs.panels(siteDF[,c("Latitude","bio01_mean","log_bio12_mean","log_bio01_sq","log_Npp","overlap_unnorm_obs","log_median_richness")])



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

## transforms
if (LOG_ITV) {
  if (any(dat$itv <= 0, na.rm = TRUE))
    stop("Non-positive ITV values present; choose log1p or a small constant before logging.")
  dat$itv <- log(dat$itv)
}
if (RICH_TRANSFORM == "sqrt") dat$rich <- sqrt(dat$rich)
if (RICH_TRANSFORM == "log")  dat$rich <- log(dat$rich)

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
  rich ~ c1*tmean + c2*ppt + c3*npp + c4*geodiv + d*itv

  # indirect paths to richness
  ind_tmean_itv     := b1*d
  ind_ppt_itv       := b2*d
  ind_npp_itv       := b3*d
  ind_geodiv_itv    := b4*d
  ind_tmean_npp     := a1*c3
  ind_ppt_npp       := a2*c3
  ind_tmean_npp_itv := a1*b3*d
  ind_ppt_npp_itv   := a2*b3*d

  # total effects on richness
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
  rich ~ c1*tmean + c2*ppt + c3*npp + d*itv

  # indirect paths to richness
  ind_tmean_itv     := b1*d
  ind_ppt_itv       := b2*d
  ind_npp_itv       := b3*d
  ind_tmean_npp     := a1*c3
  ind_ppt_npp       := a2*c3
  ind_tmean_npp_itv := a1*b3*d
  ind_ppt_npp_itv   := a2*b3*d

  # total effects on richness
  tot_tmean  := c1 + b1*d + a1*c3 + a1*b3*d
  tot_ppt    := c2 + b2*d + a2*c3 + a2*b3*d
  tot_npp    := c3 + b3*d
'
m3_dropNPP <- '  
  itv  ~ b1*tmean + b2*ppt + b4*geodiv
  rich ~ c1*tmean + c2*ppt + c4*geodiv + d*itv

  # indirect paths to richness
  ind_tmean_itv     := b1*d
  ind_ppt_itv       := b2*d
  ind_geodiv_itv    := b4*d

  # total effects on richness
  tot_tmean  := c1 + b1*d 
  tot_ppt    := c2 + b2*d
  tot_geodiv := c4 + b4*d
'
m4_climITV <- '
  itv  ~ b1*tmean + b2*ppt 
  rich ~ c1*tmean + c2*ppt + d*itv

  # indirect paths to richness
  ind_tmean_itv     := b1*d
  ind_ppt_itv       := b2*d

  # total effects on richness
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
## 5. Targeted nested tests (the two questions that matter)
## ============================================================ ##

## ============================================================ ##
## 6. (optional) visualize best/full model
## ============================================================ ##
library(lavaanPlot)

lavaanPlot(model = fits$`7_Geo`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`5_DropClimate`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`3_dropNPP`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`4_climITV`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions
library(tidySEM)

# Create a default graph from the fitted model
graph_sem(fits$`3_dropNPP`)
graph_sem(fits$`7_Geo`)
graph_sem(fits$`4_climITV`)
graph_sem(fits$`5_DropClimate`)

ggplot(siteDF, aes(x=median_richness, y=overlap_norm_obs)) +
  geom_point(aes(colour = overlap_norm_dir))
