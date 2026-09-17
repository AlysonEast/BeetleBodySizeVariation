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

pairs.panels(plotDF[,c("bio_1","bio_12","Npp","Velocity","niche_range_obs","overlap_depth_obs","overlap_unnorm_obs","richness")])
pairs.panels(plotDF[,c("bio_1","bio_12","Npp","Velocity","niche_range_obs","overlap_depth_obs","overlap_unnorm_obs","richness","log_richness")])

## ============================================================ ##
## 1. CONFIG -- edit these, everything downstream is parameterized
## ============================================================ ##

RANGE_COL     <-"niche_range_obs"
COOCCURANCE_COL     <-"overlap_depth_obs"
RICH_COL   <- "log_richness"
TMEAN_COL  <- "bio_1" 
NPP_COL  <- "Npp"       
VELOCITY_COL  <- "Velocity"       


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
  range      = plotDF[[RANGE_COL]],
  cooccurrence = plotDF[[COOCCURANCE_COL]],
  rich   = plotDF[[RICH_COL]]
)

## complete-case across ALL model variables so every candidate model is fit on
## identical rows (required for valid AIC/BIC comparison). With the current
## Complexity column all 47 sites should be retained -- verify in the printout.
model_vars <- c("tmean", "npp", "cooccurrence", "range", "rich","velocity")
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
  rich ~ c1*tmean + c2*npp + c3*velocity
'
sem_env_direct<-sem(m_env_direct, data = dat, estimator = "MLR")
lavaanPlot(model = sem_env_direct,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions

#___________________Env and Range___________________
m_env_range <- '
  range ~ r1*tmean + r2*npp + r3*velocity

  rich ~ c1*tmean + c2*npp + c3*velocity + 
         d1*range

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_npp_range := r2*d1
  ind_velocity_range := r3*d1
  
  # total effects on richness
  tot_tmean := c1 + r1*d1
  tot_npp := c2 + r2*d1
  tot_velocity := c3 + r3*d1
'
sem_env_range<-sem(m_env_range, data = dat, estimator = "MLR")
lavaanPlot(model = sem_env_range,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions

m1_env_range <- '
  range ~ r1*tmean + r3*velocity

  rich ~ c1*tmean + c2*npp + c3*velocity + 
         d1*range

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_velocity_range := r3*d1
  
  # total effects on richness
  tot_tmean := c1 + r1*d1
  tot_npp := c2 
  tot_velocity := c3 + r3*d1
'
sem1_env_range<-sem(m1_env_range, data = dat, estimator = "MLR")

m2_env_range <- '
  range ~ r2*npp + r3*velocity

  rich ~  c2*npp + c3*velocity + 
         d1*range

  # indirect paths to richness
  ind_npp_range := r2*d1
  ind_velocity_range := r3*d1
  
  # total effects on richness
  tot_npp := c2 + r2*d1
  tot_velocity := c3 + r3*d1
'
sem2_env_range<-sem(m2_env_range, data = dat, estimator = "MLR")

m3_env_range <- '
  range ~ r3*velocity

  rich ~  c2*npp + c3*velocity + 
         d1*range

  # indirect paths to richness
  ind_velocity_range := r3*d1
  
  # total effects on richness
  tot_npp := c2 
  tot_velocity := c3 + r3*d1
'
sem3_env_range<-sem(m3_env_range, data = dat, estimator = "MLR")

#___________________Env and Depth___________________
m_env_depth <- '
  cooccurrence ~ o1*tmean + o2*npp + o3*velocity
            
  rich ~ c1*tmean + c2*npp + c3*velocity + 
         d2*cooccurrence

  # indirect paths to richness
  ind_temp_co := o1*d2
  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  
  # total effects on richness
  tot_tmean := c1 + o1*d2
  tot_npp := c2 +  o2*d2
  tot_velocity := c3 + o3*d2
'
sem_env_depth<-sem(m_env_depth, data = dat, estimator = "MLR")
lavaanPlot(model = sem_env_depth,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions


m1_env_depth <- '
  cooccurrence ~ o1*tmean + o2*npp
            
  rich ~ c1*tmean + c2*npp + c3*velocity + 
         d2*cooccurrence

  # indirect paths to richness
  ind_temp_co := o1*d2
  ind_npp_co := o2*d2
  
  # total effects on richness
  tot_tmean := c1 + o1*d2
  tot_npp := c2 +  o2*d2
  tot_velocity := c3
'
sem1_env_depth<-sem(m1_env_depth, data = dat, estimator = "MLR")
lavaanPlot(model = sem1_env_depth,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions

m2_env_depth <- '
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

m3_env_depth <- '
  cooccurrence ~ o2*npp
            
  rich ~ c2*npp + c3*velocity + 
         d2*cooccurrence

  # indirect paths to richness
  ind_npp_co := o2*d2

  # total effects on richness
  tot_npp := c2 +  o2*d2
  tot_velocity := c3
'
sem3_env_depth<-sem(m3_env_depth, data = dat, estimator = "MLR")


#___________________Env Range and Depth___________________
m_env_range_depth <- '
  range ~ r1*tmean + r2*npp + r3*velocity

  cooccurrence ~ o1*tmean + o2*npp + o3*velocity
            
  rich ~ c1*tmean + c2*npp + c3*velocity + 
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_npp_range := r2*d1
  ind_velocity_range := r3*d1
  
  ind_temp_co := o1*d2
  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  
  # total effects on richness
  tot_tmean := c1 + r1*d1 + o1*d2
  tot_npp := c2 + r2*d1 + o2*d2
  tot_velocity := c3 + r3*d1 + o3*d2
'
sem_env_range_depth<-sem(m_env_range_depth, data = dat, estimator = "MLR")
lavaanPlot(model = sem_env_range_depth,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions
summary(sem_env_range_depth)

m1_env_range_depth <- '
  range ~ r1*tmean + r3*velocity

  cooccurrence ~ o1*tmean + o2*npp 
            
  rich ~ c1*tmean + c2*npp + c3*velocity + 
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_velocity_range := r3*d1
  
  ind_temp_co := o1*d2
  ind_npp_co := o2*d2

  # total effects on richness
  tot_tmean := c1 + r1*d1 + o1*d2
  tot_npp := c2 + o2*d2
  tot_velocity := c3 + r3*d1
'
sem1_env_range_depth<-sem(m1_env_range_depth, data = dat, estimator = "MLR")

m2_env_range_depth <- '
  range ~ r1*tmean + r3*velocity

  cooccurrence ~ o2*npp 
            
  rich ~ c1*tmean + c2*npp + c3*velocity + 
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_temp_range := r1*d1
  ind_velocity_range := r3*d1
  
  ind_npp_co := o2*d2

  # total effects on richness
  tot_tmean := c1 + r1*d1 
  tot_npp := c2 + o2*d2
  tot_velocity := c3 + r3*d1
'
sem2_env_range_depth<-sem(m2_env_range_depth, data = dat, estimator = "MLR")
summary(sem2_env_range_depth)

m3_env_range_depth <- '
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
  "m1_Env_Range"               = m1_env_range,
  "m2_Env_Range"               = m2_env_range,
  "m3_Env_Range"               = m3_env_range,
  "m1_Env_Depth"               = m1_env_depth,
  "m2_Env_Depth"               = m2_env_depth,
  "m3_Env_Depth"               = m3_env_depth,
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

lavaanPlot(model = fits$`m1_Env_Range`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`m1_Env_Depth`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`m1_Env_Range_Depth`,
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
p2 <- make_sem_graph(fits$`m1_Env_Range`, lay)
p3 <- make_sem_graph(fits$`m1_Env_Depth`, lay)
p4 <- make_sem_graph(fits$`m1_Env_Range_Depth`, lay)

library(patchwork)
png("./Figures/SEMs/plotsSEMsLDG.png", res = 300, height = 7, width = 10, units = "in")
(p1 | p2) /
  (p3 | p4)
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

# =====================================================================
# X. Effect decomposition for top model (m2): direct / indirect / total
# =====================================================================
# Parameters ----------------------------------------------------------
m<-sem(m2, data = dat, estimator = "MLR")

FIT      <- m                 # point at your fitted sem() object for m2
STD      <- TRUE                   # TRUE = standardized (est.std); FALSE = raw
CI_LEVEL <- 0.95
FIG_DIR  <- "./Figures/PathEffects"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

library(lavaan)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(svglite)
library(dplyr)
library(forcats)

# Predictor colours (replace with the Tableau palette from TPDexample.R) ---
pred_cols <- c(Temperature   = "#4E79A7",
               Precipitation = "#59A14F",
               Complexity    = "#F28E2B",
               Mediator      = "#BAB0AC")

# 1. Pull the parameter table (standardized or unstandardized) --------
pe <- if (STD) {
  standardizedSolution(FIT, level = CI_LEVEL) |> rename(est = est.std)
} else {
  parameterEstimates(FIT, level = CI_LEVEL, standardized = FALSE)
}
# columns used downstream: lhs, op, rhs, est, ci.lower, ci.upper

# 2. Direct effects on richness (all rich ~ paths) --------------------
direct_df <- pe %>%
  filter(op == "~", lhs == "rich") %>%
  mutate(
    label = recode(rhs,
                   tmean      = "Temperature",
                   tmean_sq   = "Temperature\u00B2 (curv.)",
                   ppt        = "Precipitation",
                   Complexity = "Complexity",
                   Overlap    = "Overlap \u2192 Rich",
                   Range      = "Range \u2192 Rich"),
    predictor = case_when(
      rhs %in% c("tmean", "tmean_sq") ~ "Temperature",
      rhs == "ppt"                    ~ "Precipitation",
      rhs == "Complexity"             ~ "Complexity",
      TRUE                            ~ "Mediator"),
    sig = ci.lower > 0 | ci.upper < 0)

# 3. Indirect effects (the ind_* defined parameters) ------------------
indirect_df <- pe %>%
  filter(op == ":=", grepl("^ind_", lhs)) %>%
  mutate(
    label = recode(lhs,
                   ind_tmean_Overlap            = "Temp \u2192 Overlap",
                   ind_ppt_Overlap              = "Precip \u2192 Overlap",
                   ind_Complexity_Overlap       = "Complexity \u2192 Overlap",
                   ind_tmean_Range              = "Temp \u2192 Range",
                   ind_Complexity_Range         = "Complexity \u2192 Range",
                   ind_tmean_Range_Overlap      = "Temp \u2192 Range \u2192 Overlap",
                   ind_Complexity_Range_Overlap = "Complexity \u2192 Range \u2192 Overlap"),
    predictor = case_when(
      grepl("tmean", lhs)      ~ "Temperature",
      grepl("ppt", lhs)        ~ "Precipitation",
      grepl("Complexity", lhs) ~ "Complexity",
      TRUE                     ~ "Mediator"),
    sig = ci.lower > 0 | ci.upper < 0)

# 4. Total effects (the tot_* defined parameters) ---------------------
total_df <- pe %>%
  filter(op == ":=", grepl("^tot_", lhs)) %>%
  mutate(
    label = recode(lhs,
                   tot_tmean      = "Temperature",
                   tot_ppt        = "Precipitation",
                   tot_Complexity = "Complexity"),
    predictor = case_when(
      lhs == "tot_tmean"      ~ "Temperature",
      lhs == "tot_ppt"        ~ "Precipitation",
      lhs == "tot_Complexity" ~ "Complexity",
      TRUE                    ~ "Mediator"),
    sig = ci.lower > 0 | ci.upper < 0)

# 5. Shared plotting helper (used across all three panels) ------------
effect_plot <- function(df, title) {
  ggplot(df, aes(x = est, y = fct_reorder(label, est),
                 colour = predictor, alpha = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(xmin = ci.lower, xmax = ci.upper),
                    fatten = 3, linewidth = 0.6) +
    scale_colour_manual(values = pred_cols, drop = FALSE, name = "Predictor") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
    labs(x = if (STD) "Standardized effect on richness" else "Effect on richness",
         y = NULL, title = title) +
    theme_pubr(legend = "right") +
    theme(plot.title = element_text(face = "bold", size = 11))
}

# 6. Build and combine -------------------------------------------------
p_direct   <- effect_plot(direct_df,   "Direct effects on richness")
p_indirect <- effect_plot(indirect_df, "Indirect effects on richness")
p_total    <- effect_plot(total_df,    "Total effects on richness")

combined <- (p_direct / p_indirect / p_total) +
  plot_layout(guides = "collect", heights = c(1, 1, 0.6)) +
  plot_annotation(
    title    = "Model m2: decomposition of effects on richness",
    subtitle = if (STD) "Standardized paths, 95% CI (faded = CI spans 0)"
    else       "Unstandardized paths, 95% CI (faded = CI spans 0)",
    tag_levels = "A")

combined
# 7. Save (filename embeds the STD parameter) -------------------------
ggsave(file.path(FIG_DIR, sprintf("m2_effects_%s.svg", if (STD) "std" else "raw")),
       combined, width = 8, height = 10, device = svglite::svglite)

## ============================================================ ##
## 6. Visualize TOP model (m2): two mediators, Range -> Overlap
## ============================================================ ##
library(lavaan); library(ggplot2); library(ggpubr); library(dplyr)

USE_BOOT <- TRUE
N_BOOT   <- 2000   # bump to 5000 for the final figure

fit_m2 <- sem(m2, data = dat, estimator = "MLR",
              se = if (USE_BOOT) "bootstrap" else "standard",
              bootstrap = N_BOOT, iseed = 42)

## R^2 for the three endogenous responses (Range, Overlap, rich)
cat("\n--- R^2 (endogenous) ---\n"); print(round(lavInspect(fit_m2, "rsquare"), 3))

pe <- parameterEstimates(fit_m2, standardized = TRUE, ci = TRUE)

## structural paths (r* = env->Range, b* = env/Range->Overlap,
##                    c*/q1 = ->rich, d1 = Overlap->rich, d2 = Range->rich)
paths <- subset(pe, op == "~",
                c("lhs","rhs","label","est","ci.lower","ci.upper","pvalue","std.all"))
cat("\n--- structural paths (std.all = fully standardized) ---\n")
print(paths, row.names = FALSE, digits = 3)

## effect decomposition on richness (the := lines in m2)
eff <- subset(pe, op == ":=",
              c("label","est","ci.lower","ci.upper","pvalue"))
cat("\n--- effects on richness: indirect (ind_*), totals (tot_*), curvature ---\n")
print(eff, row.names = FALSE, digits = 3)

## relative contribution ranking: |standardized effect on richness|.
## Temperature splits into a linear-route total plus curvature (level-dependent,
## so read curv_tmean and the marginal slopes alongside it).
rank_tbl <- data.frame(
  driver = c("temperature (linear route)", "temperature (curvature)",
             "precipitation", "complexity", "Overlap (direct)", "Range (direct)"),
  effect = c(eff$est[eff$label=="tot_tmean"],  eff$est[eff$label=="curv_tmean"],
             eff$est[eff$label=="tot_ppt"],    eff$est[eff$label=="tot_Complexity"],
             paths$std.all[paths$label=="d1"], paths$std.all[paths$label=="d2"]))
rank_tbl <- rank_tbl[order(-abs(rank_tbl$effect)), ]
cat("\n--- relative contribution (|standardized effect on richness|) ---\n")
print(rank_tbl, row.names = FALSE, digits = 3)

## path diagrams -------------------------------------------------------
library(lavaanPlot)
lavaanPlot(model = fit_m2, coefs = TRUE, stand = TRUE, sig = 0.05,
           stars = c("regress"), graph_options = list(rankdir = "LR"))

library(tidySEM)
lay <- get_layout(
  "tmean", "tmean_sq", "ppt",     "Complexity",
  NA,      "Range",   "Overlap",  NA,
  NA,       NA,       "rich",     NA,
  rows = 3)
graph_sem(fit_m2, layout = lay)

## coefficient grabber + back-transform helpers -----------------------
gb  <- function(l) pe$est[pe$label == l]           # labeled path OR := effect
mu  <- function(v) mean(dat_raw[[v]]); sdv <- function(v) sd(dat_raw[[v]])
c1<-gb("c1"); c2<-gb("c2"); c4<-gb("c4"); q1<-gb("q1")
d1<-gb("d1"); d2<-gb("d2")
r1<-gb("r1"); r4<-gb("r4"); b1<-gb("b1"); b2<-gb("b2"); b4<-gb("b4"); b5<-gb("b5")

## ============================================================ ##
## 6.1 Model-implied trends on richness: TOTAL vs DIRECT
##     direct = coefficient straight into the richness equation
##     total  = tot_* from the := block (all mediated routes summed)
##     Pulling total from the defined effect keeps the line and the model
##     in lockstep: tmean/Complexity route through Overlap, Range, and
##     Range->Overlap; ppt routes through Overlap only.
## ============================================================ ##
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
  if (mark_vertex && quad != 0) {
    vz <- -total_slope / (2*quad)
    if (vz >= min(z) & vz <= max(z))
      p <- p + geom_vline(xintercept = vz*sdv(v)+mu(v),
                          linetype = 3, colour = "grey60")
  }
  p
}

p_temp <- trend_panel("tmean", direct_slope = c1,
                      total_slope = gb("tot_tmean"),
                      xlab = "Mean annual temp (bio_1)",
                      quad = q1, mark_vertex = TRUE)
p_ppt  <- trend_panel("ppt", direct_slope = c2,
                      total_slope = gb("tot_ppt"), xlab = "Precipitation (bio_12)")
p_comp <- trend_panel("Complexity", direct_slope = c4,
                      total_slope = gb("tot_Complexity"), xlab = "Geodiversity / complexity")

ggarrange(p_temp, p_ppt, p_comp, ncol = 3,
          common.legend = TRUE, legend = "bottom", labels = "AUTO")

## ============================================================ ##
## 6.2 Focal mechanism: each mediator -> richness (slopes d1, d2)
## ============================================================ ##
mech_panel <- function(med, coef, xlab, col) {
  z <- seq(min(dat[[med]]), max(dat[[med]]), length.out = 100)
  ggplot() +
    geom_point(data = data.frame(m = dat_raw[[med]], rich = dat_raw$rich),
               aes(m, rich), alpha = .55) +
    geom_line(data = data.frame(m = z*sdv(med)+mu(med),
                                rich = (coef*z)*sdv("rich")+mu("rich")),
              aes(m, rich), linewidth = 1.1, colour = col) +
    labs(x = xlab, y = "Estimated richness") + theme_pubr()
}
p_ov <- mech_panel("Overlap", d1, "Body-size overlap", "#4576b5")
p_rg <- mech_panel("Range",   d2, "Body-size range",   "#1f6f6f")

## ============================================================ ##
## 6.3 Mediator drivers: env -> Range, env -> Overlap, Range -> Overlap
## ============================================================ ##
driver_panel <- function(v, coef, xlab, med, col = "#555599") {
  z <- seq(min(dat[[v]]), max(dat[[v]]), length.out = 100)
  ggplot() +
    geom_point(data = data.frame(x = dat_raw[[v]], m = dat_raw[[med]]),
               aes(x, m), alpha = .55) +
    geom_line(data = data.frame(x = z*sdv(v)+mu(v),
                                m = (coef*z)*sdv(med)+mu(med)),
              aes(x, m), linewidth = 1, colour = col) +
    labs(x = xlab, y = med) + theme_pubr()
}
## env -> Range
p_rg_t <- driver_panel("tmean",      r1, "Mean annual temp (bio_1)", "Range")
p_rg_c <- driver_panel("Complexity", r4, "Geodiversity / complexity", "Range")
## env -> Overlap (+ Range -> Overlap, the cross-mediator link b5)
p_ov_t <- driver_panel("tmean",      b1, "Mean annual temp (bio_1)", "Overlap")
p_ov_p <- driver_panel("ppt",        b2, "Precipitation (bio_12)",   "Overlap")
p_ov_c <- driver_panel("Complexity", b4, "Geodiversity / complexity", "Overlap")
p_ov_r <- driver_panel("Range",      b5, "Body-size range",           "Overlap", col = "#1f6f6f")

ggarrange(p_rg_t, p_rg_c, p_ov,
          p_ov_t, p_ov_p, p_ov_c,
          p_ov_r, p_rg,   NULL,
          ncol = 3, nrow = 3, labels = "AUTO")

## ============================================================ ##
## 6.4 Effect-decomposition forest plot (bootstrap CIs)
## ============================================================ ##
fp_labels <- c(
  "tot_tmean", "tot_ppt", "tot_Complexity",         # totals
  "curv_tmean", "slope_tmean_cold", "slope_tmean_warm",  # temp curvature
  "ind_tmean_Overlap", "ind_ppt_Overlap", "ind_Complexity_Overlap",   # via Overlap
  "ind_tmean_Range", "ind_Complexity_Range",                          # via Range
  "ind_tmean_Range_Overlap", "ind_Complexity_Range_Overlap",          # via Range->Overlap
  "d1", "d2")                                        # mediator direct effects
fp <- subset(pe, label %in% fp_labels, c("label","est","ci.lower","ci.upper"))
fp$label <- factor(fp$label, levels = rev(fp_labels))

ggplot(fp, aes(est, label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = ci.lower, xmax = ci.upper)) +
  labs(x = "Standardized effect on richness (bootstrap CI)", y = NULL,
       title = "m2: effect decomposition") + theme_pubr()
