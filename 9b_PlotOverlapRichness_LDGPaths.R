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
library(MVN)


setwd("/home/aly/Beetles/BeetleBodySizeVariation")

## ============================================================ ##
## 0. Assemble plot data and visulally inspect to choose parameters
## ============================================================ ##
#Read in and merge overlap and richness data
# plot_overlap<-read.csv("./Outputs/plot_by_all_noaug_ByYearAvg_IndividualNull.csv") #use plot_by_all becuase there are no exclusions due to domains with 1 site
#Read in overlap data
plot_2018<-read.csv("./Outputs/plot_by_site_2018_PoolNull.csv")
plot_2019<-read.csv("./Outputs/plot_by_site_2019_PoolNull.csv")
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

phylo<-read.csv("./Outputs/plot_PhyloDistance.csv")
plotDF<-merge(plotDF, phylo, by="Assemblage")
head(plotDF)

#How stable is overlap from year to year
plotDF2018<-subset(plotDF, Year.x==2018)
plotDF2019<-subset(plotDF, Year.x==2019)
pair<-merge(plotDF2018, plotDF2019, by="plotID.x", all=TRUE)
head(pair)

plot(pair$n_comm_sp.x~pair$n_comm_sp.y)
abline(a=0, b=1)

plot(pair$niche_range_obs.x~pair$niche_range_obs.y)
abline(a=0, b=1)

plot(pair$overlap_depth_obs.x~pair$overlap_depth_obs.y)
abline(a=0, b=1)

plot(pair$pd.x~pair$pd.y)
abline(a=0, b=1)

plotDF$richness<-plotDF$Estimator
#What overlap values need to be removed?
plot(plotDF$richness~plotDF$n_comm_sp)
abline(a=0, b=1)
plot(plotDF$Observed~plotDF$n_comm_sp)
abline(a=0, b=1)

plotDF$diff<-plotDF$Observed-plotDF$n_comm_sp
hist(plotDF$diff)

plotDF$diffpct<-((plotDF$Estimator-plotDF$n_comm_sp)/plotDF$Estimator)
# plotDF$diffpct<-as.numeric(ifelse(plotDF$diffpct<0, paste0(NA), plotDF$diffpct))

table(plotDF$diffpct, useNA = "ifany")
hist(plotDF$diffpct)
plotDF$diffdouble<-ifelse(plotDF$diffpct>.5, paste0(1), paste0(0))
plotDF$diffthird<-ifelse(plotDF$diffpct>(2/3), paste0(1), paste0(0))



ggplot(plotDF, aes(x=richness, y=n_comm_sp, colour = diffpct, shape = diffdouble)) +
  geom_point(alpha=0.5, size=3) +
  geom_errorbar(aes(xmin = LCL, xmax=UCL), alpha=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradient(low = "purple", high = "orange")

table(plotDF$diffdouble)

#Evaluate validity of richness estimates
ggplot(plotDF, aes(x=richness, y=n_comm_sp, colour = completeness, shape = diffdouble)) +
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


ggplot(plotDF, aes(x=richness, y=n_comm_sp, colour = poorRichnessEstimate, shape = diffthird)) +
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
plotDF<-subset(plotDF, diffpct<.5 & n_comm_sp<=2 | 
                 n_comm_sp>2 & diffpct<=(2/3))

dim(preExclusion)
dim(plotDF)
dim(preExclusion)[1]-dim(plotDF)[1]

symdiff(levels(as.factor(preExclusion$plotID.x)),levels(as.factor(plotDF$plotID.x)))
length(symdiff(levels(as.factor(preExclusion$plotID.x)),levels(as.factor(plotDF$plotID.x))))
dim(table(plotDF$plotID.x))

write.csv(plotDF, "./Outputs/finalplotDF.csv")

symdiff(levels(as.factor(preExclusion$SiteID)),levels(as.factor(plotDF$SiteID)))
#Evaluate validity of richness estimates
ggplot(preExclusion, aes(x=richness, y=n_comm_sp)) +
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

plot(pair$n_comm_sp.x~pair$n_comm_sp.y)
abline(a=0, b=1)

plot(pair$niche_range_obs.x~pair$niche_range_obs.y)
abline(a=0, b=1)

plot(plotDF$richness~plotDF$overlap_depth_obs)

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
  annotate(geom = "text", x = 1.15, y = 2, label = "Highest Potential \n Richness", size = 5)+
  annotate(geom = "text", x = .1, y = .15, label = "Lowest Potential \n Richness", size = 5)+
  annotate(geom = "text", x = .1, y = 2, label = "Lowest Total \n Partitioning", size = 5)+
  annotate(geom = "text", x = 1.15, y = .15, label = "Highest Total \n Partitioning", size = 5)+
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
  annotate(geom = "text", x = 1.15, y = 2, label = "Highest Potential \n Richness", size = 5)+
  annotate(geom = "text", x = .1, y = .15, label = "Lowest Potential \n Richness", size = 5)+
  annotate(geom = "text", x = .1, y = 2, label = "Lowest Total \n Partitioning", size = 5)+
  annotate(geom = "text", x = 1.15, y = .15, label = "Highest Total \n Partitioning", size = 5)
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
plotDF$log_rugosity<-log10(plotDF$rugosity_RC +0.01)
plotDF$sqrt_rugosity<-sqrt(plotDF$rugosity_RC)

plotDF$log_geodiv<-log10(plotDF$geodiv)

pairs.panels(plotDF[,c("bio_1","Npp","Velocity","rugosity_RC","geodiv","niche_range_obs","overlap_depth_obs","pd","mpd","mntd","richness")])
pairs.panels(plotDF[,c("bio_1","Npp","Velocity","log_rugosity","log_geodiv","niche_range_obs","overlap_depth_obs","pd","mpd","richness","log_richness")])

n<-mvn(data = plotDF[,c("bio_1","Npp","Velocity","log_geodiv",
                     "niche_range_obs","overlap_depth_obs","mpd",
                     "log_richness")], multivariate_outlier_method="adj")
n$multivariate_normality
head(n$multivariate_outliers)
plot(n)
#View(plotDF[n$multivariate_outliers$Observation,])
# pairs.panels(plotDF[n$multivariate_outliers$Observation,c("bio_1","Npp","Velocity","log_geodiv",
#                                                           "niche_range_obs","overlap_depth_obs","mpd",
#                                                           "log_richness")])
# pairs.panels(plotDF[,c("bio_1","Npp","Velocity","log_geodiv",
#                        "niche_range_obs","overlap_depth_obs","mpd",
#                        "log_richness")])

n<-mvn(data = plotDF[,c("bio_1","Npp","Velocity","log_geodiv",
                     "niche_range_obs","overlap_depth_obs","mpd",
                     "richness")], multivariate_outlier_method="adj")
n$multivariate_normality
head(n$multivariate_outliers)
plot(n)
# View(plotDF[n$multivariate_outliers$Observation,])
# pairs.panels(plotDF[n$multivariate_outliers$Observation,c("bio_1","Npp","Velocity","log_geodiv",
#                                                           "niche_range_obs","overlap_depth_obs","mpd",
#                                                           "richness")])
# pairs.panels(plotDF[,c("bio_1","Npp","Velocity","log_geodiv",
#                        "niche_range_obs","overlap_depth_obs","mpd",
#                        "richness")])

## ============================================================ ##
## 1. CONFIG -- edit these, everything downstream is parameterized
## ============================================================ ##

RANGE_COL     <-"niche_range_obs"
COOCCURANCE_COL     <-"overlap_depth_obs"
PHYLO_COL <- "mpd"
RICH_COL   <- "richness"
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
  siteID = plotDF$siteID,
  plotID = plotDF$plotID,
  tmean  = plotDF[[TMEAN_COL]],
  npp    = plotDF[[NPP_COL]],
  velocity    = plotDF[[VELOCITY_COL]],
  geodiv = plotDF[[GEODIV_COL]],
  range      = plotDF[[RANGE_COL]],
  cooccurrence = plotDF[[COOCCURANCE_COL]],
  phylo = plotDF[[PHYLO_COL]],
  rich   = plotDF[[RICH_COL]]
)

# dat$phylo[is.na(dat$phylo)]<-0
## complete-case across ALL model variables so every candidate model is fit on
## identical rows (required for valid AIC/BIC comparison). With the current
## Complexity column all 47 sites should be retained -- verify in the printout.
model_vars <- c("tmean", "npp","geodiv", "cooccurrence", "range", "rich","velocity","phylo")
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

png("./Figures/SEMs/plotsPairs.png", res = 300, height = 13, width = 13, units = "in")
pairs.panels(dat[,c(2:ncol(dat))])
dev.off()

## ============================================================ ##
## 3. Candidate Models
## ============================================================ ##
#___________________Env Range and Depth___________________

Hyp_driven_model <- '
  range ~ r1*tmean + r3*velocity

  cooccurrence ~ o2*npp + o3*velocity + o4*geodiv

  phylo ~ p1*tmean + p2*npp + p3*velocity + p4*geodiv

  range ~~ cooccurrence
  phylo ~~ cooccurrence
  phylo ~~ range

  rich ~ c1*tmean + c2*npp + c3*velocity + c4*geodiv +
         d1*range + d2*cooccurrence + d3*phylo

  # indirect paths to richness
  ind_temp_range := r1*d1
  #ind_npp_range := r2*d1
  ind_velocity_range := r3*d1

  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  ind_spatial_co := o4*d2
  
  ind_temp_phylo := p1*d3
  ind_npp_phylo := p2*d3
  ind_velocity_phylo := p3*d3
  ind_spatial_phylo := p4*d3


  # total effects on richness
  tot_tmean := c1 + r1*d1 + p1*d3
  tot_npp := c2 + o2*d2
  tot_velocity := c3 + r3*d1 + o3*d2 + p3*d3
  tot_spatial := c4 + o4*d2 + p4*d3
'
Hyp_driven_model_nophylo <- '
  range ~ r1*tmean + r3*velocity

  cooccurrence ~ o2*npp + o3*velocity + o4*geodiv

  range ~~ cooccurrence

  rich ~ c1*tmean + c2*npp + c3*velocity + c4*geodiv +
         d1*range + d2*cooccurrence

  # indirect paths to richness
  ind_temp_range := r1*d1
  # ind_npp_range := r2*d1
  ind_velocity_range := r3*d1

  ind_npp_co := o2*d2
  ind_velocity_co := o3*d2
  ind_spatial_co := o4*d2

  # total effects on richness
  tot_tmean := c1 + r1*d1
  tot_npp := c2 + o2*d2
  tot_velocity := c3 + r3*d1 + o3*d2
  tot_spatial := c4 + o4*d2
'

sem_Hyp_driven_model_grouped<-sem(Hyp_driven_model, 
                                  data = dat, 
                                  estimator = "ML",
                                  se = "bootstrap",
                                  cluster = "siteID")

sem_Hyp_driven_model_noPhylo_grouped<-sem(Hyp_driven_model_nophylo, 
                                          data = dat, 
                                          estimator = "ML",
                                          se = "bootstrap",
                                          cluster = "siteID")

sem_Hyp_driven_model_grouped_exclude<-sem(Hyp_driven_model, 
                                          data = subset(dat, !siteID=="BARR"), 
                                          estimator = "ML",
                                          se = "bootstrap",
                                          cluster = "siteID")

sem_Hyp_driven_model_noPhylo_grouped_exclude<-sem(Hyp_driven_model_nophylo, 
                                                  data = subset(dat, !siteID=="BARR"), 
                                                  estimator = "ML",
                                                  se = "bootstrap",
                                                  cluster = "siteID")

## ============================================================ ##
## 3. evaluate Models
## ============================================================ ##
# AIC(sem_Saturated_model_grouped, sem_Saturated_model_no_direct_grouped,
#     sem_Hyp_driven_model_grouped, sem_Hyp_driven_model_nodirect_grouped)

fitMeasures(sem_Hyp_driven_model_grouped,
            c("chisq","df","pvalue","pvalue.scaled","cfi","tli","rmsea","srmr"))
fitMeasures(sem_Hyp_driven_model_grouped_exclude,
            c("chisq","df","pvalue","pvalue.scaled","cfi","tli","rmsea","srmr"))

fitMeasures(sem_Hyp_driven_model_noPhylo_grouped,
            c("chisq","df","pvalue","pvalue.scaled","cfi","tli","rmsea","srmr"))
fitMeasures(sem_Hyp_driven_model_noPhylo_grouped_exclude,
            c("chisq","df","pvalue","pvalue.scaled","cfi","tli","rmsea","srmr"))


modindices(sem_Hyp_driven_model_grouped, sort. = TRUE, minimum.value = 3.84)
residuals(sem_Hyp_driven_model_grouped, type = "cor")

summary(sem_Hyp_driven_model_grouped)
summary(sem_Hyp_driven_model_grouped_exclude)

# Write the captured console text to a file
writeLines(capture.output(summary(sem_Hyp_driven_model_grouped)), con = "plot_report.txt")
writeLines(capture.output(summary(sem_Hyp_driven_model_grouped_exclude)), con = "plot_report_excludeBARR.txt")

writeLines(capture.output(summary(sem_Hyp_driven_model_noPhylo_grouped)), con = "plot_reportNoPhylo.txt")
writeLines(capture.output(summary(sem_Hyp_driven_model_noPhylo_grouped_exclude)), con = "plot_reportNoPhylo_excludeBARR.txt")


## ============================================================ ##
## 6. visualize best/full model
## ============================================================ ##
library(lavaanPlot)


lavaanPlot(model = sem_Hyp_driven_model_grouped,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = sem_Hyp_driven_model_grouped_exclude,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = sem_Hyp_driven_model_noPhylo_grouped,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = sem_Hyp_driven_model_noPhylo_grouped_exclude,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions


lay <- get_layout(
  "velocity","range", NA,
  "tmean", NA,  "rich",
  NA, "phylo", NA,
  "npp",NA, NA,
  "geodiv", "cooccurrence", NA,
  rows = 5)
#
lay2 <- get_layout(
  "velocity","range", NA,
  "tmean", NA,  "rich",
  NA, "NA", NA,
  "npp",NA, NA,
  "geodiv", "cooccurrence", NA,
  rows = 5)
#
make_sem_graph <- function(model, layout, scale = 5) {
  g <- prepare_graph(model = model)
  # Standardized path coefficients
  g$edges$linewidth <- abs(as.numeric(g$edges$est_std)) * scale
  graph_sem(model, layout = layout)
}

library(patchwork)
png("./Figures/SEMs/plotsSEMstars.png", res = 300, height = 10, width = 11, units = "in")
make_sem_graph(sem_Hyp_driven_model_grouped_exclude, lay)
dev.off()

png("./Figures/SEMs/plotsSEMNoPhylostars.png", res = 300, height = 10, width = 11, units = "in")
make_sem_graph(sem_Hyp_driven_model_grouped_exclude, lay2)
dev.off()


library(lavaan)
library(semPlot)

lay <- matrix(
  c(.7,  0,   # range 
    0,  1.5,   # cooccurrence 
    0, -1.5,   # phylo
    3,  0,      # rich
    -3,  1.5,   # tmean
    -3,  0.5,   # npp
    -3, -0.5,   # velocity
    -3, -1.5),   # geodiv
  ncol = 2, byrow = TRUE)

png("./Figures/SEMs/plotsSEMLDG.png", res = 300, height = 10, width = 13, units = "in")
semPaths(
  sem_Hyp_driven_model_grouped,
  # layout = lay,
  what = "std",
  whatLabels = "std",
  residuals = TRUE,
  # edge.color = c("black", "grey75"),
  exoVar = FALSE, exoCov = FALSE,
  sizeMan = 8,
  label.cex = 1.1,
  edge.label.cex = .85,
  edge.width = 2,
  fade = TRUE,
  curve = 2)
dev.off()
