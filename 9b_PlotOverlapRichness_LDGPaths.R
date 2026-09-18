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
library(tidySEM)
library(lavaanPlot)


setwd("/home/aly/Beetles/BeetleBodySizeVariation")

## ============================================================ ##
## 0. Assemble plot data and visulally inspect to choose parameters
## ============================================================ ##
#Read in and merge overlap and richness data
# plot_overlap<-read.csv("./Outputs/plot_by_all_noaug_ByYearAvg_IndividualNull.csv") #use plot_by_all becuase there are no exclusions due to domains with 1 site
#Read in overlap data
plot_2018<-read.csv("./Outputs/plot_by_site_aug_2018_PoolNull.csv")
plot_2019<-read.csv("./Outputs/plot_by_site_aug_2019_PoolNull.csv")
head(plot_2018)
plot_2018$Year<-2018
plot_2019$Year<-2019

plot_overlap<-rbind(plot_2018, plot_2019)
plot_overlap$latitude<-NULL
plot_overlap$Assemblage<-paste0(plot_overlap$plotID,"_",plot_overlap$Year)

plot_richness<-read.csv("../BeetleBiodiversity/plot_annual_EstimatedSppRichness.csv")
plot_richness$X<-NULL
head(plot_richness)
plotDF<-merge(plot_overlap, plot_richness, by = "Assemblage", all.x = TRUE, all.y = FALSE)
head(plotDF)

plot_abund2018<-read.csv("./Data/plotTotal_abund_2018.csv")
plot_abund2019<-read.csv("./Data/plotTotal_abund_2019.csv")
head(plot_abund2018)
plot_abund2018$Assemblage<-paste0(plot_abund2018$plotID,"_2018")
plot_abund2019$Assemblage<-paste0(plot_abund2019$plotID,"_2019")
plot_abund<-rbind(plot_abund2018, plot_abund2019)
head(plot_abund)

plotDF<-merge(plotDF, plot_abund, by="Assemblage")
head(plotDF)

#How stable is overlap from year to year
plotDF2018<-subset(plotDF, Year.x==2018)
plotDF2019<-subset(plotDF, Year.x==2019)
pair<-merge(plotDF2018, plotDF2019, by="plotID.x", all=TRUE)
head(pair)

plot(pair$n_overlap_sp.x~pair$n_overlap_sp.y)
abline(a=0, b=1)

plot(pair$overlap_unnorm_obs.x~pair$overlap_unnorm_obs.y)
abline(a=0, b=1)
plot(sqrt(pair$overlap_unnorm_obs.x)~sqrt(pair$overlap_unnorm_obs.y))
abline(a=0, b=1)

plot(pair$niche_range_obs.x~pair$niche_range_obs.y)
abline(a=0, b=1)

plot(pair$overlap_depth_obs.x~pair$overlap_depth_obs.y)
abline(a=0, b=1)


plotDF$richness<-plotDF$Estimator
#What overlap values need to be removed?
plot(plotDF$richness~plotDF$n_overlap_sp)
abline(a=0, b=1)
plot(plotDF$Observed~plotDF$n_overlap_sp)
abline(a=0, b=1)

plotDF$diff<-plotDF$Observed-plotDF$n_overlap_sp
hist(plotDF$diff)

plotDF$diffpct<-((plotDF$Estimator-plotDF$n_overlap_sp)/plotDF$Estimator)
# plotDF$diffpct<-as.numeric(ifelse(plotDF$diffpct<0, paste0(NA), plotDF$diffpct))

table(plotDF$diffpct, useNA = "ifany")
hist(plotDF$diffpct)
plotDF$diffdouble<-ifelse(plotDF$diffpct>.5, paste0(1), paste0(0))
plotDF$diffthird<-ifelse(plotDF$diffpct>(2/3), paste0(1), paste0(0))


ggplot(plotDF, aes(x=richness, y=n_overlap_sp, colour = overlap_unnorm_obs)) +
  geom_point(alpha=0.5) +
  geom_errorbar(aes(xmin = LCL, xmax=UCL), alpha=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradient(low = "purple", high = "orange")

ggplot(plotDF, aes(x=richness, y=n_overlap_sp, colour = diffpct, shape = diffdouble)) +
  geom_point(alpha=0.5, size=3) +
  geom_errorbar(aes(xmin = LCL, xmax=UCL), alpha=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradient(low = "purple", high = "orange")

table(plotDF$diffdouble)

#Evaluate validity of richness estimates
ggplot(plotDF, aes(x=richness, y=n_overlap_sp, colour = completeness, shape = diffdouble)) +
  geom_point(alpha=0.5, size=3) +
  geom_errorbar(aes(xmin = LCL, xmax=UCL), alpha=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradient(low = "purple", high = "orange")

hist(plotDF$completeness)

plotDF$poorRichnessEstimate<-ifelse(plotDF$completeness<.5, paste0(1), paste0(0))
table(plotDF$poorRichnessEstimate)
table(plotDF$poorRichnessEstimate, plotDF$diffdouble)
table(plotDF$poorRichnessEstimate, plotDF$diffthird)
table(plotDF$poorRichnessEstimate, plotDF$plotID.x)


ggplot(plotDF, aes(x=richness, y=n_overlap_sp, colour = poorRichnessEstimate, shape = diffthird)) +
  geom_point(alpha=0.5, size=3) +
  geom_errorbar(aes(xmin = LCL, xmax=UCL), alpha=0.5) +
  geom_abline(intercept = 0, slope = 1) 

#### Exclusion ####
preExclusion<-plotDF

EXCLUDE_ISLANDS <- TRUE
if (EXCLUDE_ISLANDS) plotDF<-plotDF %>% 
  filter(!grepl('PUUM', plotDF$plotID.x),
         !grepl('LAJA', plotDF$plotID.x),
         !grepl('GUAN', plotDF$plotID.x)) #c("PUUM","LAJA","GUAN"))

plotDF<-subset(plotDF, completeness>=.5)
plotDF<-subset(plotDF, diffpct<.5 & n_overlap_sp<=2 | 
                 n_overlap_sp>2 & diffpct<=(2/3))

dim(preExclusion)
plotDF<-subset(plotDF, !is.na(overlap_unnorm_obs))
dim(plotDF)
dim(preExclusion)[1]-dim(plotDF)[1]

symdiff(levels(as.factor(preExclusion$plotID.x)),levels(as.factor(plotDF$plotID.x)))
length(symdiff(levels(as.factor(preExclusion$plotID.x)),levels(as.factor(plotDF$plotID.x))))
dim(table(plotDF$plotID.x))

symdiff(levels(as.factor(preExclusion$SiteID)),levels(as.factor(plotDF$SiteID)))
#Evaluate validity of richness estimates
ggplot(preExclusion, aes(x=richness, y=n_overlap_sp)) +
  geom_point(alpha=0.5, size=2, col="grey") +
  geom_errorbar(aes(xmin = LCL, xmax=UCL), alpha=0.5, col="grey") +
  geom_point(data = plotDF, alpha=0.5, size=2, col="black") +
  geom_errorbar(data = plotDF, aes(xmin = LCL, xmax=UCL), alpha=0.5, col="black") +
  geom_abline(intercept = 0, slope = 1) +
  theme_pubr()

ggplot(preExclusion, aes(x=richness)) +
  geom_histogram(fill="grey") +
  geom_histogram(data = plotDF, alpha=0.5, col="black") +
  theme_pubr()


plotDF2018<-subset(plotDF, Year.x==2018)
plotDF2019<-subset(plotDF, Year.x==2019)
pair<-merge(plotDF2018, plotDF2019, by="plotID.x", all=TRUE)
head(pair)

plot(pair$n_overlap_sp.x~pair$n_overlap_sp.y)
abline(a=0, b=1)

plot(pair$overlap_unnorm_obs.x~pair$overlap_unnorm_obs.y)
abline(a=0, b=1)
plot(sqrt(pair$overlap_unnorm_obs.x)~sqrt(pair$overlap_unnorm_obs.y))
abline(a=0, b=1)

plot(pair$niche_range_obs.x~pair$niche_range_obs.y)
abline(a=0, b=1)

plot(plotDF$richness~plotDF$overlap_unnorm_obs)
plot(plotDF$n_overlap_sp~plotDF$overlap_unnorm_obs)

plot(plotDF$richness~plotDF$overlap_depth_obs)

plot(plotDF$overlap_norm_obs~plotDF$overlap_depth_obs)
plot(plotDF$niche_range_obs~plotDF$overlap_depth_obs)
plot(plotDF$overlap_depth_obs~plotDF$niche_range_obs)

png("./Figures/PlotRichnessMetrics.png", res = 300, height = 8, width = 8, units = "in")
ggplot(plotDF, aes(x=niche_range_obs, y=overlap_depth_obs, colour = richness)) +
  geom_point(size=4) +
  scale_colour_gradient(low = "purple", high = "orange") +
  theme_pubr() +
  xlab("Niche Space") +
  ylab("Average Co-occurance") +
  labs(colour = "Observed \n Richness") +
  annotate(geom = "text", x = 1.2, y = 3.1, label = "Highest Potential \n Richness", size = 5)+
  annotate(geom = "text", x = .17, y = .3, label = "Lowest Potential \n Richness", size = 5)+
  annotate(geom = "text", x = .17, y = 3, label = "Lowest Total \n Partitioning", size = 5)+
  annotate(geom = "text", x = 1.2, y = .3, label = "Highest Total \n Partitioning", size = 5)+
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.7))
dev.off()
png("./Figures/PlotRichnessMetrics_blank.png", res = 300, height = 8, width = 8, units = "in")
ggplot(plotDF, aes(x=niche_range_obs, y=overlap_depth_obs)) +
  geom_point(size=4, colour="white") +
  theme_pubr() +
  xlab("Niche Space") +
  ylab("Average Co-occurance") +
  labs(colour = "Observed \n Richness") +
  annotate(geom = "text", x = 1.2, y = 3.1, label = "Highest Potential \n Richness", size = 5)+
  annotate(geom = "text", x = .17, y = .3, label = "Lowest Potential \n Richness", size = 5)+
  annotate(geom = "text", x = .17, y = 3, label = "Lowest Total \n Partitioning", size = 5)+
  annotate(geom = "text", x = 1.2, y = .3, label = "Highest Total \n Partitioning", size = 5)
dev.off()

#Env Variaibles
struc<-read.csv("./Outputs/BETplot_Rugosity.csv")
struc$X<-NULL
env<-read.csv("./Outputs/BeetlePlotswEnvData.csv")
NPP<-read.csv("../NEON_MODIS_NPP_2018_2019.csv") #from https://code.earthengine.google.com/b41a55076352b2d9e21ac5e74bf337bc
plotDF$plotID<-plotDF$plotID.x
plotDF$plotID.x<-NULL
plotDF$plotID.y<-NULL
velocity<-read.csv("./Outputs/BeetlePlotswVelocity.csv")
head(velocity)
velocity<-velocity[,c("plotID","Velocity")]

plotDF<-merge(plotDF, struc, by="plotID")
plotDF<-merge(plotDF, env, by="plotID")
plotDF<-merge(plotDF, NPP[,c("Npp","Gpp","plotID")], by="plotID")
plotDF<-merge(plotDF, velocity, by="plotID")
head(plotDF)

#### Pair plot#
plotDF$log_richness<-log10(plotDF$richness)
head(plotDF)
plotDF$log_rugosity<-log10(plotDF$rugosity_RC)
plotDF$log_geodiv<-log10(plotDF$geodiv)

pairs.panels(plotDF[,c("bio_1","Npp","Velocity","rugosity_RC","geodiv","niche_range_obs","overlap_depth_obs","overlap_unnorm_obs","richness")])
pairs.panels(plotDF[,c("bio_1","Npp","Velocity","log_rugosity","log_geodiv","niche_range_obs","overlap_depth_obs","overlap_unnorm_obs","richness","log_richness")])

## ============================================================ ##
## 1. CONFIG -- edit these, everything downstream is parameterized
## ============================================================ ##

RANGE_COL     <-"niche_range_obs"
COOCCURANCE_COL     <-"overlap_depth_obs"
RICH_COL   <- "log_richness"
TMEAN_COL  <- "bio_1" 
NPP_COL  <- "Npp"       
VELOCITY_COL  <- "Velocity"       
GEODIV_COL <- "log_geodiv"


## Transforms (applied before standardizing)
STANDARDIZE    <- TRUE               # z-score all model vars (coeffs in SD units)

## ============================================================ ##
## 2. Build modeling frame: select, rename, transform, complete-case, scale
## ============================================================ ##

dat <- data.frame(
  plotID = plotDF$plotID,
  tmean  = plotDF[[TMEAN_COL]],
  npp    = plotDF[[NPP_COL]],
  velocity    = plotDF[[VELOCITY_COL]],
  geodiv = plotDF[[GEODIV_COL]],
  range      = plotDF[[RANGE_COL]],
  cooccurrence = plotDF[[COOCCURANCE_COL]],
  rich   = plotDF[[RICH_COL]]
)

## complete-case across ALL model variables so every candidate model is fit on
## identical rows (required for valid AIC/BIC comparison). With the current
## Complexity column all 47 sites should be retained -- verify in the printout.
model_vars <- c("tmean", "npp","geodiv", "cooccurrence", "range", "rich","velocity")
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
head(dat)
hist(dat$rich)

pairs.panels(dat[,c(2:ncol(dat))])


## ============================================================ ##
## 3. Candidate Models
## ============================================================ ##
#___________________Env Only___________________
m_env_direct <- '
  rich ~ c1*tmean + c2*npp + c3*velocity + c4*geodiv
'
sem_env_direct<-sem(m_env_direct, data = dat, estimator = "MLR")
lavaanPlot(model = sem_env_direct,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions

#___________________Env and Range___________________
m_env_range <- '
  range ~ r1*tmean + r2*npp + r3*velocity + r4*geodiv

  rich ~ c1*tmean + c2*npp + c3*velocity + c4*geodiv +
         d1*range

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_npp_range := r2*d1
  ind_velocity_range := r3*d1
  ind_geodiv_range := r4*d1
  
  # total effects on richness
  tot_tmean := c1 + r1*d1
  tot_npp := c2 + r2*d1
  tot_velocity := c3 + r3*d1
  tot_geodiv := c4 + r4*d1
'
sem_env_range<-sem(m_env_range, data = dat, estimator = "MLR")
lavaanPlot(model = sem_env_range,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions

m1_env_range <-  '
  range ~ r1*tmean + r3*velocity + r4*geodiv

  rich ~ c1*tmean + c2*npp + c3*velocity + c4*geodiv +
         d1*range

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_velocity_range := r3*d1
  ind_geodiv_range := r4*d1
  
  # total effects on richness
  tot_tmean := c1 + r1*d1
  tot_npp := c2
  tot_velocity := c3 + r3*d1
  tot_geodiv := c4 + r4*d1
'
sem1_env_range<-sem(m1_env_range, data = dat, estimator = "MLR")

m2_env_range <-  '
  range ~ r1*tmean + r3*velocity

  rich ~ c1*tmean + c2*npp + c3*velocity + c4*geodiv +
         d1*range

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_velocity_range := r3*d1

  # total effects on richness
  tot_tmean := c1 + r1*d1
  tot_npp := c2
  tot_velocity := c3 + r3*d1
  tot_geodiv := c4
'
sem2_env_range<-sem(m2_env_range, data = dat, estimator = "MLR")

m3_env_range <- '
  range ~ r3*velocity

  rich ~ c2*npp + c3*velocity + c4*geodiv +
         d1*range

  # indirect paths to richness
  ind_velocity_range := r3*d1

  # total effects on richness
  tot_npp := c2
  tot_velocity := c3 + r3*d1
  tot_geodiv := c4
'
sem3_env_range<-sem(m3_env_range, data = dat, estimator = "MLR")

#___________________Env and Depth___________________
m_env_depth <- '
  cooccurrence ~ o1*tmean + o2*npp + o3*velocity + o4*geodiv
            
  rich ~ c1*tmean + c2*npp + c3*velocity + c4*geodiv +
         d2*cooccurrence

  # indirect paths to richness
  ind_temp_co := o1*d2
  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  ind_geodiv_co := o4*d2
  
  # total effects on richness
  tot_tmean := c1 + o1*d2
  tot_npp := c2 +  o2*d2
  tot_velocity := c3 + o3*d2
  tot_velocity := c4 + o4*d2
'
sem_env_depth<-sem(m_env_depth, data = dat, estimator = "MLR")
lavaanPlot(model = sem_env_depth,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions


m1_env_depth <- '
  cooccurrence ~ o1*tmean + o2*npp + o4*geodiv
            
  rich ~ c1*tmean + c2*npp + c3*velocity + c4*geodiv +
         d2*cooccurrence

  # indirect paths to richness
  ind_temp_co := o1*d2
  ind_npp_co := o2*d2
  ind_geodiv_co := o4*d2
  
  # total effects on richness
  tot_tmean := c1 + o1*d2
  tot_npp := c2 +  o2*d2
  tot_velocity := c3
  tot_geodiv := c4 +  o4*d2
'
sem1_env_depth<-sem(m1_env_depth, data = dat, estimator = "MLR")
lavaanPlot(model = sem1_env_depth,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions

m2_env_depth <- '
  cooccurrence ~ o2*npp + o3*velocity + o4*geodiv
            
  rich ~ c2*npp + c3*velocity + c4*geodiv +
         d2*cooccurrence

  # indirect paths to richness
  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  ind_geodiv_co := o4*d2
  
  # total effects on richness
  tot_npp := c2 +  o2*d2
  tot_velocity := c3 + o3*d2
  tot_geodiv := c4 +  o4*d2
'
sem2_env_depth<-sem(m2_env_depth, data = dat, estimator = "MLR")

m3_env_depth <- '
  cooccurrence ~ o2*npp + o4*geodiv
            
  rich ~ c2*npp + c3*velocity + c4*geodiv +
         d2*cooccurrence

  # indirect paths to richness
  ind_npp_co := o2*d2
  ind_geodiv_co := o4*d2
  
  # total effects on richness
  tot_npp := c2 +  o2*d2
  tot_velocity := c3 
  tot_geodiv := c4 +  o4*d2
'
sem3_env_depth<-sem(m3_env_depth, data = dat, estimator = "MLR")


#___________________Env Range and Depth___________________
m_env_range_depth <- '
  range ~ r1*tmean + r2*npp + r3*velocity + r4*geodiv

  cooccurrence ~ o1*tmean + o2*npp + o3*velocity + o4*geodiv
            
  rich ~ c1*tmean + c2*npp + c3*velocity + c4*geodiv +
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_npp_range := r2*d1
  ind_velocity_range := r3*d1
  ind_geodiv_range := r4*d1
  
  ind_temp_co := o1*d2
  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  ind_geodiv_co := o4*d2
  
  # total effects on richness
  tot_tmean := c1 + r1*d1 + o1*d2
  tot_npp := c2 + r2*d1 + o2*d2
  tot_velocity := c3 + r3*d1 + o3*d2
  tot_geodiv := c4 + r4*d1 + o4*d2
'
sem_env_range_depth<-sem(m_env_range_depth, data = dat, estimator = "MLR")
lavaanPlot(model = sem_env_range_depth,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions
summary(sem_env_range_depth)

m1_env_range_depth <- '
  range ~ r1*tmean + r3*velocity + r4*geodiv

  cooccurrence ~ o1*tmean + o2*npp + o3*velocity + o4*geodiv
            
  rich ~ c1*tmean + c2*npp + c3*velocity + c4*geodiv +
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_velocity_range := r3*d1
  ind_geodiv_range := r4*d1
  
  ind_temp_co := o1*d2
  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  ind_geodiv_co := o4*d2
  
  # total effects on richness
  tot_tmean := c1 + r1*d1 + o1*d2
  tot_npp := c2 + o2*d2
  tot_velocity := c3 + r3*d1 + o3*d2
  tot_geodiv := c4 + r4*d1 + o4*d2
'
sem1_env_range_depth<-sem(m1_env_range_depth, data = dat, estimator = "MLR")

m2_env_range_depth <- '
  range ~ r1*tmean + r3*velocity

  cooccurrence ~ o1*tmean + o2*npp + o3*velocity + o4*geodiv
            
  rich ~ c1*tmean + c2*npp + c3*velocity + c4*geodiv +
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_velocity_range := r3*d1
  
  ind_temp_co := o1*d2
  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  ind_geodiv_co := o4*d2
  
  # total effects on richness
  tot_tmean := c1 + r1*d1 + o1*d2
  tot_npp := c2 + o2*d2
  tot_velocity := c3 + r3*d1 + o3*d2
  tot_geodiv := c4 + o4*d2
'
sem2_env_range_depth<-sem(m2_env_range_depth, data = dat, estimator = "MLR")
summary(sem2_env_range_depth)

m3_env_range_depth <- '
  range ~ r3*velocity

  cooccurrence ~ o2*npp + o3*velocity + o4*geodiv
            
  rich ~ c2*npp + c3*velocity + c4*geodiv +
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_velocity_range := r3*d1
  
  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  ind_geodiv_co := o4*d2
  
  # total effects on richness
  tot_npp := c2 + o2*d2
  tot_velocity := c3 + r3*d1 + o3*d2
  tot_geodiv := c4 + o4*d2
'
sem3_env_range_depth<-sem(m3_env_range_depth, data = dat, estimator = "MLR")
summary(sem3_env_range_depth)

m4_env_range_depth <- '
  range ~ r1*tmean + r3*velocity

  cooccurrence ~ o1*tmean + o2*npp + o4*geodiv
            
  rich ~ c1*tmean + c2*npp + c3*velocity + c4*geodiv +
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_velocity_range := r3*d1
  
  ind_temp_co := o1*d2
  ind_npp_co := o2*d2
  ind_geodiv_co := o4*d2
  
  # total effects on richness
  tot_tmean := c1 + r1*d1 + o1*d2
  tot_npp := c2 + o2*d2
  tot_velocity := c3 + r3*d1
  tot_geodiv := c4 + o4*d2
'
sem4_env_range_depth<-sem(m4_env_range_depth, data = dat, estimator = "MLR")
summary(sem4_env_range_depth)

m5_env_range_depth <- '
  range ~ r3*velocity

  cooccurrence ~ o2*npp + o4*geodiv
            
  rich ~ c2*npp + c3*velocity + c4*geodiv +
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_velocity_range := r3*d1
  
  ind_npp_co := o2*d2
  ind_geodiv_co := o4*d2
  
  # total effects on richness
  tot_npp := c2 + o2*d2
  tot_velocity := c3 + r3*d1
  tot_geodiv := c4 + o4*d2
'
sem5_env_range_depth<-sem(m5_env_range_depth, data = dat, estimator = "MLR")
summary(sem5_env_range_depth)


m6_env_range_depth <- '
  range ~ r3*velocity

  cooccurrence ~ o2*npp 
            
  rich ~ c2*npp + c3*velocity +
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_velocity_range := r3*d1
  
  ind_npp_co := o2*d2
  
  # total effects on richness
  tot_npp := c2 + o2*d2
  tot_velocity := c3 + r3*d1
'
sem6_env_range_depth<-sem(m6_env_range_depth, data = dat, estimator = "MLR")
summary(sem6_env_range_depth)

m7_env_range_depth <- '
  range ~ r1*tmean + r2*npp + r4*geodiv

  cooccurrence ~ o1*tmean + o2*npp + o4*geodiv
            
  rich ~ c1*tmean + c2*npp + c4*geodiv +
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_npp_range := r2*d1
  ind_geodiv_range := r4*d1
  
  ind_temp_co := o1*d2
  ind_npp_co := o2*d2
  ind_geodiv_co := o4*d2
  
  # total effects on richness
  tot_tmean := c1 + r1*d1 + o1*d2
  tot_npp := c2 + r2*d1 + o2*d2
  tot_geodiv := c4 + r4*d1 + o4*d2
'
sem7_env_range_depth<-sem(m7_env_range_depth, data = dat, estimator = "MLR")
lavaanPlot(model = sem7_env_range_depth,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions
summary(sem7_env_range_depth)

m8_env_range_depth <- '
  range ~ r2*npp + r3*velocity + r4*geodiv

  cooccurrence ~ o2*npp + o3*velocity + o4*geodiv
            
  rich ~ c2*npp + c3*velocity + c4*geodiv +
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_npp_range := r2*d1
  ind_velocity_range := r3*d1
  ind_geodiv_range := r4*d1
  
  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  ind_geodiv_co := o4*d2
  
  # total effects on richness
  tot_npp := c2 + r2*d1 + o2*d2
  tot_velocity := c3 + r3*d1 + o3*d2
  tot_geodiv := c4 + r4*d1 + o4*d2
'
sem8_env_range_depth<-sem(m8_env_range_depth, data = dat, estimator = "MLR")
lavaanPlot(model = sem8_env_range_depth,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions
summary(sem8_env_range_depth)

## ============================================================ ##
## 3. evaluate Models
## ============================================================ ##

models <- list(
  "Full_env_only"              = m_env_direct,
  "Full_Env_Depth"             = m_env_depth,
  "Full_Env_Range"             = m_env_range,
  "Full_Env_Range_Depth"       = m_env_range_depth,
  "m1_Env_Range"               = m1_env_range,
  "m2_Env_Range"               = m2_env_range,
  "m3_Env_Range"               = m3_env_range,
  "m1_Env_Depth"               = m1_env_depth,
  "m2_Env_Depth"               = m2_env_depth,
  "m3_Env_Depth"               = m3_env_depth,
  "m1_Env_Range_Depth"         = m1_env_range_depth,
  "m2_Env_Range_Depth"         = m2_env_range_depth,
  "m3_Env_Range_Depth"         = m3_env_range_depth,
  "m4_Env_Range_Depth"         = m4_env_range_depth,
  "m5_Env_Range_Depth"         = m5_env_range_depth,
  "m6_Env_Range_Depth"         = m6_env_range_depth,
  "m7_Env_Range_Depth"         = m7_env_range_depth,
  "m8_Env_Range_Depth"         = m8_env_range_depth
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

lavaanPlot(model = fits$`m2_Env_Range`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`m2_Env_Depth`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`m2_Env_Range_Depth`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions
lavaanPlot(model = fits$`m1_Env_Range_Depth`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions
lavaanPlot(model = fits$`m4_Env_Range_Depth`,
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
  "geodiv", NA, NA,
  rows = 4)

make_sem_graph <- function(model, layout, scale = 5) {
  g <- prepare_graph(model = model)
  # Standardized path coefficients
  g$edges$linewidth <- abs(as.numeric(g$edges$est_std)) * scale
  graph_sem(model, layout = layout)
}

p1 <- make_sem_graph(fits$`Full_env_only`, lay)
p2 <- make_sem_graph(fits$`m2_Env_Range`, lay)
p3 <- make_sem_graph(fits$`m2_Env_Depth`, lay)
p4 <- make_sem_graph(fits$`m4_Env_Range_Depth`, lay)

library(patchwork)
png("./Figures/SEMs/plotsSEMsLDG.png", res = 300, height = 10, width = 13, units = "in")
(p1 | p2) /
  (p3 | p4)
dev.off()

p4b <- make_sem_graph(fits$`m2_Env_Range_Depth`, lay)
p4c <- make_sem_graph(fits$`m3_Env_Range_Depth`, lay)
png("./Figures/SEMs/plotsSEMsLDG_EnvRangeDepthEquivelant.png", res = 300, height = 10, width = 13, units = "in")
(p4 | p4b | p4c)
dev.off()


graph_sem(fits$`Full_env_only`, layout = lay)
graph_sem(fits$`m1_Env_Range`, layout = lay)
graph_sem(fits$`m1_Env_Depth`, layout = lay)
graph_sem(fits$`m1_Env_Range_Depth`, layout = lay)

graph_data <- prepare_graph(model = fits$`m1_Env_Range_Depth`)
graph_data$edges$linewidth <- abs(as.numeric(graph_data$edges$est)) * 5
plot(graph_data)


library(semPlot)
par(mfrow=c(2,2))
semPaths(fits$`Full_env_only`, 
         what = "std",          # Proportional thickness based on standardized paths
         layout = "tree",
         fade = FALSE,
         # --- FONT & LABEL SIZE CUSTOMIZATION ---
         edge.color = "black",   # Consistent line color
         edge.label.cex = 3,  # Enlarges the path coefficient numbers (Default is 1.0)
         sizeLat = 10,          # Enlarges the text/box size for Latent variables
         sizeMan = 14,          # Enlarges the text/box size for Manifest/observed variables
         label.cex = 1.2)       # Globally scales up node text size inside the boxes
semPaths(fits$`m1_Env_Range`, 
         what = "std",          # Proportional thickness based on standardized paths
         layout = "tree",
         fade = FALSE,
         # --- FONT & LABEL SIZE CUSTOMIZATION ---
         edge.color = "black",   # Consistent line color
         edge.label.cex = 3,  # Enlarges the path coefficient numbers (Default is 1.0)
         sizeLat = 10,          # Enlarges the text/box size for Latent variables
         sizeMan = 14,          # Enlarges the text/box size for Manifest/observed variables
         label.cex = 1.2)       # Globally scales up node text size inside the boxes
semPaths(fits$`m1_Env_Depth`, 
         what = "std",          # Proportional thickness based on standardized paths
         layout = "tree",
         fade = FALSE,
         # --- FONT & LABEL SIZE CUSTOMIZATION ---
         edge.color = "black",   # Consistent line color
         edge.label.cex = 3,  # Enlarges the path coefficient numbers (Default is 1.0)
         sizeLat = 10,          # Enlarges the text/box size for Latent variables
         sizeMan = 14,          # Enlarges the text/box size for Manifest/observed variables
         label.cex = 1.2)       # Globally scales up node text size inside the boxes
semPaths(fits$`m1_Env_Range_Depth`, 
         what = "std",          # Proportional thickness based on standardized paths
         layout = "tree",
         fade = FALSE,
         # --- FONT & LABEL SIZE CUSTOMIZATION ---
         edge.color = "black",   # Consistent line color
         edge.label.cex = 3,  # Enlarges the path coefficient numbers (Default is 1.0)
         sizeLat = 10,          # Enlarges the text/box size for Latent variables
         sizeMan = 14,          # Enlarges the text/box size for Manifest/observed variables
         label.cex = 1.2)       # Globally scales up node text size inside the boxes
