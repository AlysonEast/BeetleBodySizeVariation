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
site_2018<-read.csv("./Outputs/site_by_all_2018_IndividualNull.csv")
site_2019<-read.csv("./Outputs/site_by_all_2019_IndividualNull.csv")
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

phylo<-read.csv("./Outputs/site_PhyloDistance.csv")
siteDF<-merge(siteDF, phylo, by="Assemblage")
head(siteDF)

#How stable is overlap from year to year
siteDF2018<-subset(siteDF, Year.x==2018)
siteDF2019<-subset(siteDF, Year.x==2019)
pair<-merge(siteDF2018, siteDF2019, by="siteID.x", all=TRUE)
head(pair)

plot(pair$n_comm_sp.x~pair$n_comm_sp.y)
abline(a=0, b=1)

plot(pair$niche_range_obs.x~pair$niche_range_obs.y)
abline(a=0, b=1)

plot(pair$pd.x~pair$pd.y)

siteDF$richness<-siteDF$Estimator
#What overlap values need to be removed?
plot(siteDF$richness~siteDF$n_comm_sp)
abline(a=0, b=1)
plot(siteDF$Observed~siteDF$n_comm_sp)
abline(a=0, b=1)

siteDF$diff<-siteDF$Observed-siteDF$n_comm_sp
hist(siteDF$diff)

siteDF$diffpct<-((siteDF$Observed-siteDF$n_comm_sp)/siteDF$Observed)

table(siteDF$diffpct, useNA = "ifany")
hist(siteDF$diffpct)
siteDF$diffdouble<-ifelse(siteDF$diffpct>.5, paste0(1), paste0(0))
siteDF$diffthird<-ifelse(siteDF$diffpct>(2/3), paste0(1), paste0(0))


ggplot(siteDF, aes(x=richness, y=n_comm_sp, colour = overlap_depth_obs)) +
  geom_point(alpha=0.5) +
  geom_errorbar(aes(xmin = LCL, xmax=UCL), alpha=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradient(low = "purple", high = "orange")

ggplot(siteDF, aes(x=richness, y=n_comm_sp, colour = diffpct, shape = diffdouble)) +
  geom_point(alpha=0.5, size=3) +
  geom_errorbar(aes(xmin = LCL, xmax=UCL), alpha=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_gradient(low = "purple", high = "orange")

table(siteDF$diffdouble)
table(siteDF$diffsig)

#Evaluate validity of richness estimates
ggplot(siteDF, aes(x=richness, y=n_comm_sp, colour = completeness, shape = diffdouble)) +
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


ggplot(siteDF, aes(x=richness, y=n_comm_sp, colour = poorRichnessEstimate, shape = diffthird)) +
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
#siteDF<-subset(siteDF, !is.na(overlap_depth_obs))

dim(preExclusion)
dim(siteDF)
dim(preExclusion)[1]-dim(siteDF)[1]

write.csv(siteDF, "./Outputs/finalsiteDF.csv")


symdiff(levels(as.factor(preExclusion$siteID.x)),levels(as.factor(siteDF$siteID.x)))
dim(table(siteDF$siteID.x))




ggplot(preExclusion, aes(x=richness, y=n_comm_sp)) +
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
  annotate(geom = "text", x = 1.15, y = 3.8, label = "Highest Potential \n Richness", size = 5)+
  annotate(geom = "text", x = .25, y = .35, label = "Lowest Potential \n Richness", size = 5)+
  annotate(geom = "text", x = .25, y = 3.8, label = "Lowest Total \n Partitioning", size = 5)+
  annotate(geom = "text", x = 1.15, y = .35, label = "Highest Total \n Partitioning", size = 5)+
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
  annotate(geom = "text", x = 1.15, y = 3.8, label = "Highest Potential \n Richness", size = 5)+
  annotate(geom = "text", x = .25, y = .35, label = "Lowest Potential \n Richness", size = 5)+
  annotate(geom = "text", x = .25, y = 3.8, label = "Lowest Total \n Partitioning", size = 5)+
  annotate(geom = "text", x = 1.15, y = .35, label = "Highest Total \n Partitioning", size = 5)+
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

pairs.panels(siteDF[,c("bio01_mean","Npp","Velocity","bio01_sq","niche_range_obs","overlap_depth_obs","pd","mpd","mntd","pd_ses","pd_p","richness")])
siteDF$log_richness<-log10(siteDF$richness)
siteDF$log_bio01_sq<-log10(siteDF$bio01_sq+.001)
siteDF$log_bio01_sq
siteDF$log_mpd<-sqrt(siteDF$mpd)
pairs.panels(siteDF[,c("bio01_mean","Npp","Velocity","log_bio01_sq","niche_range_obs","overlap_depth_obs","pd","log_mpd","mpd","richness","log_richness")])

n<-mvn(data = siteDF[,c("bio01_mean","Npp","Velocity","log_bio01_sq",
                        "niche_range_obs","overlap_depth_obs","mpd",
                        "log_richness")], multivariate_outlier_method="adj")
n$multivariate_normality
head(n$multivariate_outliers)
plot(n)

n<-mvn(data = siteDF[,c("bio01_mean","Npp","Velocity","log_bio01_sq",
                        "niche_range_obs","overlap_depth_obs","mpd",
                        "richness")], multivariate_outlier_method="adj")
n$multivariate_normality
head(n$multivariate_outliers)
plot(n)
View(siteDF[n$multivariate_outliers$Observation,])
## ============================================================ ##
## 1. CONFIG -- edit these, everything downstream is parameterized
## ============================================================ ##
RANGE_COL     <-"niche_range_obs"
COOCCURANCE_COL     <-"overlap_depth_obs"
PHYLO_COL <- "mpd"
RICH_COL   <- "richness"
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
  geodiv    = siteDF[[GEODIV_COL]],
  range      = siteDF[[RANGE_COL]],
  cooccurrence = siteDF[[COOCCURANCE_COL]],
  phylo = siteDF[[PHYLO_COL]],
  rich   = siteDF[[RICH_COL]]
)
## complete-case across ALL model variables so every candidate model is fit on
## identical rows (required for valid AIC/BIC comparison). With the current
## Geodiv column all 47 sites should be retained -- verify in the printout.
model_vars <- c("tmean", "npp", "cooccurrence", "range", "rich","velocity","geodiv","phylo")
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


png("./Figures/SEMs/sitesPairs.png", res = 300, height = 13, width = 13, units = "in")
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
  tot_npp := c2 + o2*d2 + p2*d3
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
writeLines(capture.output(summary(sem_Hyp_driven_model_grouped)), con = "site_report.txt")
writeLines(capture.output(summary(sem_Hyp_driven_model_grouped_exclude)), con = "site_report_excludeBARR.txt")

writeLines(capture.output(summary(sem_Hyp_driven_model_noPhylo_grouped)), con = "site_reportNoPhylo.txt")
writeLines(capture.output(summary(sem_Hyp_driven_model_noPhylo_grouped_exclude)), con = "site_reportNoPhylo_excludeBARR.txt")


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
png("./Figures/SEMs/sitesSEMstars.png", res = 300, height = 10, width = 11, units = "in")
make_sem_graph(sem_Hyp_driven_model_grouped_exclude, lay)
dev.off()

png("./Figures/SEMs/sitesSEMnoPhylostars.png", res = 300, height = 10, width = 11, units = "in")
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

png("./Figures/SEMs/sitesSEMLDG.png", res = 300, height = 10, width = 13, units = "in")
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
