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

siteDF$diffpct<-(siteDF$Observed-siteDF$n_overlap_sp)/siteDF$Observed
siteDF$diffpct<-as.numeric(ifelse(siteDF$diffpct<0, paste0(NA), siteDF$diffpct))

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
if (EXCLUDE_ISLANDS) siteDF <- siteDF %>% 
  filter(!siteID %in% c("PUUM","LAJA","GUAN"))

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
