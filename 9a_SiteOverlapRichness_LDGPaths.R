##############################################################################
## SiteOverlapRichness_Paths.R
## Site-level path analysis (SEM) of species richness.
## Builds on OverlapRichness.R.
##
## Construct mapping (from OverlapRichness.R "####Paths####" block):
##   Climate      : Tmean (bio01_mean) + Precip (bio12_mean)
##   Productivity : NPP   (added below from NEONSiteNPP.csv)
##   Heterogeneity: Geodiv (an srtm_* surface metric)  <-- CONFIRM WHICH COLUMN
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
geodiv_dir<-"/media/aly/Penobscot/NEON/Geodiversity/edi.2320.1/"


## ============================================================ ##
## 0. Assemble site data: start from siteDF, add NPP
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


ggplot(siteDF, aes(x=richness, y=n_overlap_sp, colour = overlap_norm_obs)) +
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

EXCLUDE_ISLANDS <- TRUE
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

ggplot(preExclusion, aes(x=richness, y=n_overlap_sp)) +
  geom_point(alpha=0.5, size=2, col="grey") +
  geom_errorbar(aes(xmin = LCL, xmax=UCL), alpha=0.5, col="grey") +
  geom_point(data = siteDF, alpha=0.5, size=2, col="black") +
  geom_errorbar(data = siteDF, aes(xmin = LCL, xmax=UCL), alpha=0.5, col="black") +
  geom_abline(intercept = 0, slope = 1) +
  theme_pubr()

png("./Figures/SiteRichnessMetrics.png", res = 300, height = 8, width = 8, units = "in")
ggplot(siteDF, aes(x=niche_range_obs, y=overlap_depth_obs, colour = richness)) +
  geom_point(size=4) +
  scale_colour_gradient(low = "purple", high = "orange") +
  theme_pubr() +
  xlab("Niche Space") +
  ylab("Average Co-occurance") +
  labs(colour = "Observed \n Richness") +
  annotate(geom = "text", x = 1.2, y = 3.8, label = "Highest Potential \n Richness", size = 5)+
  annotate(geom = "text", x = .25, y = .4, label = "Lowest Potential \n Richness", size = 5)+
  annotate(geom = "text", x = .25, y = 3.8, label = "Lowest Total \n Partitioning", size = 5)+
  annotate(geom = "text", x = 1.2, y = .4, label = "Highest Total \n Partitioning", size = 5)+
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.7))
dev.off()
png("./Figures/SiteRichnessMetrics_blank.png", res = 300, height = 8, width = 8, units = "in")
ggplot(siteDF, aes(x=niche_range_obs, y=overlap_depth_obs)) +
  geom_point(size=4, colour="white") +
  theme_pubr() +
  xlab("Niche Space") +
  ylab("Average Co-occurance") +
  labs(colour = "Observed \n Richness") +
  annotate(geom = "text", x = 1.2, y = 3.8, label = "Highest Potential \n Richness", size = 5)+
  annotate(geom = "text", x = .25, y = .4, label = "Lowest Potential \n Richness", size = 5)+
  annotate(geom = "text", x = .25, y = 3.8, label = "Lowest Total \n Partitioning", size = 5)+
  annotate(geom = "text", x = 1.2, y = .4, label = "Highest Total \n Partitioning", size = 5) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.9))
dev.off()


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

vel<-read.csv("./Outputs/BeetleSiteswVelocity.csv")
siteDF<-merge(siteDF, vel, by="siteID")
#### Pair site#

pairs.panels(siteDF[,c("bio01_mean","Npp","Velocity","bio01_sq","overlap_unnorm_obs","niche_range_obs","overlap_depth_obs","richness")])
siteDF$log_richness<-log10(siteDF$richness)
siteDF$log_bio01_sq<-log10(siteDF$bio01_sq+.001)
siteDF$log_bio01_sq
pairs.panels(siteDF[,c("bio01_mean","Npp","Velocity","log_bio01_sq","overlap_unnorm_obs","niche_range_obs","overlap_depth_obs","richness","log_richness")])


## ============================================================ ##
## 1. CONFIG -- edit these, everything downstream is parameterized
## ============================================================ ##
RANGE_COL     <-"niche_range_obs"
COOCCURANCE_COL     <-"overlap_depth_obs"
RICH_COL   <- "log_richness"
TMEAN_COL  <- "bio01_mean" 
NPP_COL  <- "Npp"       
VELOCITY_COL  <- "Velocity"       
GEODIV_COL  <- "log_bio01_sq"       


## Transforms (applied before standardizing)
STANDARDIZE    <- TRUE               # z-score all model vars (coeffs in SD units)

## ============================================================ ##
## 2. Build modeling frame: select, rename, transform, complete-case, scale
## ============================================================ ##

dat <- data.frame(
  siteID = siteDF$siteID,
  tmean  = siteDF[[TMEAN_COL]],
  npp    = siteDF[[NPP_COL]],
  velocity    = siteDF[[VELOCITY_COL]],
  SpatialHet    = siteDF[[GEODIV_COL]],
  range      = siteDF[[RANGE_COL]],
  cooccurrence = siteDF[[COOCCURANCE_COL]],
  rich   = siteDF[[RICH_COL]]
)
## complete-case across ALL model variables so every candidate model is fit on
## identical rows (required for valid AIC/BIC comparison). With the current
## Geodiv column all 47 sites should be retained -- verify in the printout.
model_vars <- c("tmean", "npp", "cooccurrence", "range", "rich","velocity","SpatialHet")
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
## 3. Candidate Models
## ============================================================ ##
#___________________Env Only___________________
m_env_direct <- '
  rich ~ c1*tmean + c2*npp + c3*velocity + c4*SpatialHet
'
sem_env_direct<-sem(m_env_direct, data = dat, estimator = "MLR")
lavaanPlot(model = sem_env_direct,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions

#___________________Env and Range___________________
m_env_range <- '
  range ~ r1*tmean + r2*npp + r3*velocity + r4*SpatialHet

  rich ~ c1*tmean + c2*npp + c3*velocity + c4*SpatialHet + 
         d1*range

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_npp_range := r2*d1
  ind_velocity_range := r3*d1
  ind_spatial_range := r4*d1
  
  # total effects on richness
  tot_tmean := c1 + r1*d1
  tot_npp := c2 + r2*d1
  tot_velocity := c3 + r3*d1
  tot_spatial := c4 + r4*d1

'
sem_env_range<-sem(m_env_range, data = dat, estimator = "MLR")
lavaanPlot(model = sem_env_range,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions

m1a_env_range <- '
  range ~ r1*tmean + r3*velocity + r4*SpatialHet

  rich ~ c1*tmean + c3*velocity + c4*SpatialHet + 
         d1*range

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_velocity_range := r3*d1
  ind_spatial_range := r4*d1
  
  # total effects on richness
  tot_tmean := c1 + r1*d1
  tot_velocity := c3 + r3*d1
  tot_spatial := c4 + r4*d1
'
sem1a_env_range<-sem(m1a_env_range, data = dat, estimator = "MLR")


m1b_env_range <- '
  range ~ r1*tmean + r3*velocity

  rich ~ c1*tmean + c3*velocity + 
         d1*range

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_velocity_range := r3*d1
  
  # total effects on richness
  tot_tmean := c1 + r1*d1
  tot_velocity := c3 + r3*d1
'
sem1_env_range<-sem(m1b_env_range, data = dat, estimator = "MLR")
summary(sem1_env_range)

m2_env_range <- '
  range ~ r3*velocity

  rich ~ c3*velocity + 
         d1*range

  # indirect paths to richness
  ind_velocity_range := r3*d1
  
  # total effects on richness
  tot_velocity := c3 + r3*d1
'
sem2_env_range<-sem(m2_env_range, data = dat, estimator = "MLR")

m3_env_range <- '
  range ~ r3*velocity

  rich ~ c1*tmean + c3*velocity + 
         d1*range

  # indirect paths to richness
  ind_velocity_range := r3*d1
  
  # total effects on richness
  tot_tmean := c1
  tot_velocity := c3 + r3*d1
'
sem3_env_range<-sem(m3_env_range, data = dat, estimator = "MLR")

#___________________Env and Depth___________________
m_env_depth <- '
  cooccurrence ~ o1*tmean + o2*npp + o3*velocity + o4*SpatialHet
            
  rich ~ c1*tmean + c2*npp + c3*velocity + c4*SpatialHet +
         d2*cooccurrence

  # indirect paths to richness
  ind_temp_co := o1*d2
  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  ind_spatial_co := o4*d2
  
  # total effects on richness
  tot_tmean := c1 + o1*d2
  tot_npp := c2 +  o2*d2
  tot_velocity := c3 + o3*d2
  tot_spatial := c4 + o4*d2
'

sem_env_depth<-sem(m_env_depth, data = dat, estimator = "MLR")
lavaanPlot(model = sem_env_depth,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions
summary(sem_env_depth)

m1_env_depth <- '
  cooccurrence ~ o2*npp + o3*velocity
            
  rich ~ c1*tmean + c2*npp + c3*velocity + 
         d2*cooccurrence

  # indirect paths to richness
  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  
  # total effects on richness
  tot_tmean := c1
  tot_npp := c2 +  o2*d2
  tot_velocity := c3 + o3*d2
'
sem1_env_depth<-sem(m1_env_depth, data = dat, estimator = "MLR")
lavaanPlot(model = sem1_env_depth,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions

m2_env_depth <-'
  cooccurrence ~ o2*npp + o3*velocity
            
  rich ~ c2*npp + c3*velocity + 
         d2*cooccurrence

  # indirect paths to richness
  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  
  # total effects on richness
  tot_npp := c2 +  o2*d2
  tot_velocity := c3 + o3*d2
'
sem2_env_depth<-sem(m2_env_depth, data = dat, estimator = "MLR")

#___________________Env Range and Depth___________________
m_env_range_depth <- '
  range ~ r1*tmean + r2*npp + r3*velocity + r4*SpatialHet

  cooccurrence ~ o1*tmean + o2*npp + o3*velocity + o4*SpatialHet
            
  rich ~ c1*tmean + c2*npp + c3*velocity + c4*SpatialHet +
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_npp_range := r2*d1
  ind_velocity_range := r3*d1
  ind_spatial_range := r4*d1

  ind_temp_co := o1*d2
  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  ind_spatial_co := o4*d2

  # total effects on richness
  tot_tmean := c1 + r1*d1 + o1*d2
  tot_npp := c2 + r2*d1 + o2*d2
  tot_velocity := c3 + r3*d1 + o3*d2
  tot_spatial := c4 + r4*d1 + o4*d2
'
sem_env_range_depth<-sem(m_env_range_depth, data = dat, estimator = "MLR")
lavaanPlot(model = sem_env_range_depth,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions
summary(sem_env_range_depth)

m1_env_range_depth <-  '
  range ~ r1*tmean + r3*velocity

  cooccurrence ~ o1*tmean + o2*npp + o3*velocity
            
  rich ~ c1*tmean + c2*npp + c3*velocity + 
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_velocity_range := r3*d1
  
  ind_temp_co := o1*d2
  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  
  # total effects on richness
  tot_tmean := c1 + r1*d1 + o1*d2
  tot_npp := c2 +  o2*d2
  tot_velocity := c3 + r3*d1 + o3*d2
'
sem1_env_range_depth<-sem(m1_env_range_depth, data = dat, estimator = "MLR")

m2_env_range_depth <- '
  range ~ r1*tmean + r3*velocity

  cooccurrence ~ o2*npp + o3*velocity
            
  rich ~ c1*tmean + c2*npp + c3*velocity + 
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_velocity_range := r3*d1
  
  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  
  # total effects on richness
  tot_tmean := c1 + r1*d1
  tot_npp := c2 +  o2*d2
  tot_velocity := c3 + r3*d1 + o3*d2
'
sem2_env_range_depth<-sem(m2_env_range_depth, data = dat, estimator = "MLR")
summary(sem2_env_range_depth)

m3_env_range_depth <- '
  range ~ r1*tmean + r3*velocity

  cooccurrence ~ o3*velocity
            
  rich ~ c1*tmean + c3*velocity + 
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_velocity_range := r3*d1
  
  ind_velocity_co := o3*d2
  
  # total effects on richness
  tot_tmean := c1 + r1*d1
  tot_velocity := c3 + r3*d1 + o3*d2
'
sem3_env_range_depth<-sem(m3_env_range_depth, data = dat, estimator = "MLR")
summary(sem3_env_range_depth)

## ============================================================ ##
## 3. evaluate Models
## ============================================================ ##

models <- list(
  "Full_env_only"              = m_env_direct,
  "Full_Env_Depth"             = m_env_depth,
  "Full_Env_Range"             = m_env_range,
  "Full_Env_Range_Depth"       = m_env_range_depth,
  "m1a_Env_Range"               = m1a_env_range,
  "m1b_Env_Range"               = m1b_env_range,
  "m2_Env_Range"               = m2_env_range,
  "m3_Env_Range"               = m3_env_range,
  "m1_Env_Depth"               = m1_env_depth,
  "m2_Env_Depth"               = m2_env_depth,
  "m1_Env_Range_Depth"         = m1_env_range_depth,
  "m2_Env_Range_Depth"         = m2_env_range_depth,
  "m3_Env_Range_Depth"         = m3_env_range_depth
)

fits <- lapply(models, function(spec) sem(spec, data = dat, estimator = "MLR"))

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

lavaanPlot(model = fits$`Full_env_only`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`m1_Env_Depth`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`m1b_Env_Range`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`m2_Env_Range_Depth`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions
lavaanPlot(model = fits$`m3_Env_Range_Depth`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions


library(tidySEM)
lay <- get_layout(
  "velocity", NA, "tmean", NA, "npp",
  NA, "range", NA, "cooccurrence", NA,
  NA, NA, "rich", NA, NA,
  rows = 3)
lay <- get_layout(
  "velocity","range", NA,
  "tmean", NA,  "rich",
  "npp", "cooccurrence", NA,
  rows = 3)

make_sem_graph <- function(model, layout, scale = 5) {
  g <- prepare_graph(model = model)
  # Standardized path coefficients
  g$edges$linewidth <- abs(as.numeric(g$edges$est_std)) * scale
  graph_sem(model, layout = layout)
}

p1 <- make_sem_graph(fits$`Full_env_only`, lay)
p2 <- make_sem_graph(fits$`m1b_Env_Range`, lay)
p3 <- make_sem_graph(fits$`m1_Env_Depth`, lay)
p4 <- make_sem_graph(fits$`m2_Env_Range_Depth`, lay)

library(patchwork)
png("./Figures/SEMs/sitesSEMsLDG.png", res = 300, height = 10, width = 13, units = "in")
(p1 | p2) /
  (p3 | p4)
dev.off()

summary(sem2_env_range_depth)

p4b <- make_sem_graph(fits$`m3_Env_Range_Depth`, lay)
png("./Figures/SEMs/sitesSEMsLDG_EnvRangeDepthEquivelent.png", res = 300, height = 10, width = 13, units = "in")
(p4 | p4b)
dev.off()
