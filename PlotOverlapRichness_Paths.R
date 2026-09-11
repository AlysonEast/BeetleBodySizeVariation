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


#Env Variaibles
struc<-read.csv("./Outputs/BETplot_Rugosity.csv")
struc$X<-NULL
env<-read.csv("./Outputs/BeetlePlotswEnvData.csv")
NPP<-read.csv("../NEON_MODIS_NPP_2018_2019.csv") #from https://code.earthengine.google.com/b41a55076352b2d9e21ac5e74bf337bc
plotDF$plotID<-plotDF$plotID.x
plotDF$plotID.x<-NULL
plotDF$plotID.y<-NULL

plotDF<-merge(plotDF, struc, by="plotID")
plotDF<-merge(plotDF, env, by="plotID")
plotDF<-merge(plotDF, NPP[,c("Npp","Gpp","plotID")], by="plotID")
head(plotDF)

#### Pair plot#
pairs.panels(plotDF[,c("bio_1","bio_12","rugosity_RC","Npp","abund","overlap_unnorm_obs","richness")])
plotDF$log_rugosity_RC<-log10((plotDF$rugosity_RC+0.01))
plotDF$log_abund<-log10(plotDF$abund)
plotDF$log_overlap_unnorm_obs<-log10(plotDF$overlap_unnorm_obs)
plotDF$sqrt_overlap_unnorm_obs<-sqrt(plotDF$overlap_unnorm_obs)
plotDF$log_richness<-log10(plotDF$richness)
plotDF$sqrt_richness<-sqrt(plotDF$richness)
plotDF$log_Npp<-log10(plotDF$Npp)
plotDF$log_bio_12<-log10(plotDF$bio_12)
plotDF$log_bio_1<-log10(plotDF$bio_1)
plotDF$log_comp.1<-log10(plotDF$Comp.1)
plotDF$log_sdnnd_obs<-log10(plotDF$sdnnd_obs)
pairs.panels(plotDF[,c("bio_1","log_bio_12","log_rugosity_RC","log_Npp","log_abund","overlap_unnorm_obs","sqrt_richness","log_richness")])

pairs.panels(plotDF[,c("bio_1",
                       "bio_12",
                       "rugosity_RC",
                       "Npp",
                       "niche_range_obs",
                       "overlap_unnorm_obs",
                       "richness")])

pairs.panels(plotDF[,c("bio_1","log_bio_1","bio_12","log_bio_12","rugosity_RC","log_rugosity_RC",
                       "Npp","log_Npp","log_abund","abund",
                       "overlap_unnorm_obs","sqrt_overlap_unnorm_obs",
                       "richness","sqrt_richness","log_richness")])

pairs.panels(plotDF[,c("bio_1","log_bio_1","bio_12","log_bio_12",
                       "rugosity_RC","log_rugosity_RC",
                       "Npp","log_Npp","log_abund",
                       "niche_range_obs",
                       "sqrt_overlap_unnorm_obs",
                       "overlap_unnorm_obs",
                       "overlap_norm_obs",
                       "sdnnd_obs",
                       "richness",
                       "log_richness")])

plotDF$sqrt_overlap_norm_obs<-sqrt(plotDF$overlap_norm_obs)
pairs.panels(plotDF[,c("niche_range_obs",
                       "overlap_unnorm_obs","sqrt_overlap_unnorm_obs",
                       "overlap_norm_obs","sqrt_overlap_norm_obs",
                       "sdnnd_obs",
                       "richness",
                       "log_richness")])

plot(plotDF$sqrt_overlap_unnorm_obs~plotDF$richness)
ggplot(plotDF, aes(x=richness, y=overlap_unnorm_obs))+
  geom_point()

ggplot(plotDF, aes(x=richness, y=sqrt_overlap_unnorm_obs, colour = overlap_norm_ses))+
  scale_color_gradient2(low = "#2C7FB8", mid = "grey80", high = "#D95F02", midpoint = 0) +
  geom_point() +
  theme_pubr()

ggplot(plotDF, aes(x=richness, y=overlap_unnorm_obs, colour = overlap_norm_ses))+
  scale_color_gradient2(low = "#2C7FB8", mid = "grey80", high = "#D95F02", midpoint = 0) +
  geom_point() +
  theme_pubr()
  
ggplot(plotDF, aes(x=richness, y=niche_range_obs, colour = niche_range_ses))+
  scale_color_gradient2(low = "#2C7FB8", mid = "grey80", high = "#D95F02", midpoint = 0) +
  geom_point() +
  theme_pubr()

pairs.panels(plotDF[,c("Comp.1", "Comp.2","Comp.3","Comp.4","Comp.5",
                       "richness","sqrt_richness","log_richness")])

pairs.panels(plotDF[,c("bio_1","log_bio_12",
                       "log_rugosity_RC",
                       "log_Npp",
                       "niche_range_obs",
                       "overlap_unnorm_obs",
                       "sqrt_overlap_unnorm_obs",
                       "log_overlap_unnorm_obs",
                       "sdnnd_obs",
                       "richness",
                       "log_richness")])

ggarrange(
ggplot(plotDF, aes(y=richness, col=SiteID)) +
  geom_boxplot(),
ggplot(plotDF, aes(y=richness)) +
  geom_boxplot(),
ggplot(plotDF, aes(y=overlap_unnorm_obs, col=SiteID)) +
  geom_boxplot(),
ggplot(plotDF, aes(y=overlap_unnorm_obs)) +
  geom_boxplot(),
ggplot(plotDF, aes(y=niche_range_obs, col=SiteID)) +
  geom_boxplot(),
ggplot(plotDF, aes(y=niche_range_obs)) +
  geom_boxplot(),
ggplot(plotDF, aes(y=bio_1, col=SiteID)) +
  geom_boxplot(),
ggplot(plotDF, aes(y=bio_1)) +
  geom_boxplot(),
ggplot(plotDF, aes(y=log_rugosity_RC, col=SiteID)) +
  geom_boxplot(),
ggplot(plotDF, aes(y=log_rugosity_RC)) +
  geom_boxplot(),
nrow=5, ncol=2)

## ============================================================ ##
## 1. CONFIG -- edit these, everything downstream is parameterized
## ============================================================ ##

Overlap_COL    <- "sqrt_overlap_unnorm_obs" 
Range_COL     <-"niche_range_obs"
Complexity_COL <- "log_rugosity_RC"
RICH_COL   <- "richness"
TMEAN_COL  <- "bio_1" # Second Order Mean daily mean temperature of coldest quarter
PPT_COL    <- "log_bio_12" #Mean monthly precipitation of the driest quarter
NPP_COL  <- "log_Npp"       

Climate <- "Comp.1"
Abundance <- "log_abound"


## Transforms (applied before standardizing)
STANDARDIZE    <- TRUE               # z-score all model vars (coeffs in SD units)

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
  Range      = plotDF[[Range_COL]],
  rich   = plotDF[[RICH_COL]]
)

## complete-case across ALL model variables so every candidate model is fit on
## identical rows (required for valid AIC/BIC comparison). With the current
## Complexity column all 47 sites should be retained -- verify in the printout.
model_vars <- c("tmean", "ppt", "npp", "Complexity", "Overlap", "Range", "rich")
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
dat$RangeXOverlap <- dat$Range * dat$Overlap

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

# m1_full <- '
#   npp  ~ a1*tmean + a2*ppt
#   Overlap  ~ b1*tmean + b2*ppt + b3*npp + b4*Complexity
#   rich ~ c1*tmean + q1*tmean_sq + c2*ppt + c3*npp + c4*Complexity + d*Overlap
# 
#   # indirect paths to richness
#   ind_tmean_Overlap     := b1*d
#   ind_ppt_Overlap       := b2*d
#   ind_npp_Overlap       := b3*d
#   ind_Complexity_Overlap    := b4*d
#   ind_tmean_npp     := a1*c3
#   ind_ppt_npp       := a2*c3
#   ind_tmean_npp_Overlap := a1*b3*d
#   ind_ppt_npp_Overlap   := a2*b3*d
#   
#   # temperature curvature -- the LDG test (expect q1 < 0: thermal optimum)
#   curv_tmean       := q1
#   # marginal dRich/dTmean at cold / mean / warm sites (std temp = -1, 0, +1)
#   slope_tmean_cold := c1 + 2*q1*(-1)
#   slope_tmean_mean := c1
#   slope_tmean_warm := c1 + 2*q1*(1)
# 
#   # total effects on richness
#   tot_tmean  := c1 + b1*d + a1*c3 + a1*b3*d
#   tot_ppt    := c2 + b2*d + a2*c3 + a2*b3*d
#   tot_npp    := c3 + b3*d
#   tot_Complexity := c4 + b4*d
# '


m1_full <- '
  Range  ~ r1*tmean + r3*npp + r4*Complexity
  Overlap  ~ b1*tmean + b2*ppt + b3*npp + b4*Complexity + b5*Range
  rich ~ c1*tmean + q1*tmean_sq + c2*ppt + c3*npp + c4*Complexity + d1*Overlap + d2*Range

  # indirect paths to richness
  ind_tmean_Overlap     := b1*d1
  ind_ppt_Overlap       := b2*d1
  ind_npp_Overlap       := b3*d1
  ind_Complexity_Overlap    := b4*d1

  ind_tmean_Range     := r1*d2
  ind_npp_Range       := r3*d2
  ind_Complexity_Range  := r4*d2
  ind_tmean_Range_Overlap := r1*b5*d1
  ind_npp_Range_Overlap   := r3*b5*d1
  ind_Complexity_Range_Overlap   := r4*b5*d1

  # temperature curvature -- the LDG test (expect q1 < 0: thermal optimum)
  curv_tmean       := q1
  # marginal dRich/dTmean at cold / mean / warm sites (std temp = -1, 0, +1)
  slope_tmean_cold := c1 + 2*q1*(-1)
  slope_tmean_mean := c1
  slope_tmean_warm := c1 + 2*q1*(1)

  # total effects on richness
  tot_tmean  := c1 + b1*d1 + r1*d2 + r1*b5*d1
  tot_ppt    := c2 + b2*d1
  tot_npp    := c3 + b3*d1 + r3*d2 + r3*b5*d1
  tot_Complexity := c4 + b4*d1 + r4*d2 + r4*b5*d1
'
# 
# m1_int <- '
#   # Range no longer predicts Overlap (dropped the b5 path)
#   Range   ~ r1*tmean + r3*npp + r4*Complexity
#   Overlap ~ b1*tmean + b2*ppt + b3*npp + b4*Complexity
#   rich ~ c1*tmean + q1*tmean_sq + c2*ppt + c3*npp + c4*Complexity +
#          d1*Overlap + d2*Range + d3*RangeXOverlap
# 
#   # two mediators now covary rather than being causally ordered
#   # (lavaan estimates this by default for two endogenous vars; explicit here)
#   Range ~~ Overlap
# 
#   # the interaction itself -- conditional slopes are the actual test
#   slope_Overlap_loRange   := d1 + d3*(-1)
#   slope_Overlap_meanRange := d1
#   slope_Overlap_hiRange   := d1 + d3*(1)
#   slope_Range_loOverlap   := d2 + d3*(-1)
#   slope_Range_meanOverlap := d2
#   slope_Range_hiOverlap   := d2 + d3*(1)
# 
#   # indirect paths (evaluated at the mean of the OTHER mediator)
#   ind_tmean_Overlap      := b1*d1
#   ind_ppt_Overlap        := b2*d1
#   ind_npp_Overlap        := b3*d1
#   ind_Complexity_Overlap := b4*d1
#   ind_tmean_Range        := r1*d2
#   ind_npp_Range          := r3*d2
#   ind_Complexity_Range   := r4*d2
# 
#   # temperature curvature -- the LDG test (expect q1 < 0: thermal optimum)
#   curv_tmean       := q1
#   slope_tmean_cold := c1 + 2*q1*(-1)
#   slope_tmean_mean := c1
#   slope_tmean_warm := c1 + 2*q1*(1)
# 
#   # total effects on richness (at the mean of the mediators)
#   tot_tmean      := c1 + b1*d1 + r1*d2
#   tot_ppt        := c2 + b2*d1
#   tot_npp        := c3 + b3*d1 + r3*d2
#   tot_Complexity := c4 + b4*d1 + r4*d2
# '
# 
# lavaanPlot(model = sem(m1_int, data = dat, estimator = "ML"),
#            coefs = TRUE,          # Display the path coefficients
#            stand = TRUE,          # Standardize the coefficients
#            sig = 0.05,            # Only highlight significant paths
#            stars = c("regress"))  # Append significance stars to regressions
# graph_sem(sem(m1_int, data = dat, estimator = "ML"))

lavaanPlot(model = sem(m1_full, data = dat, estimator = "ML"),
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions
graph_sem(sem(m1_full, data = dat, estimator = "ML"))

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

m2 <- ' #Take out NPP
  Range  ~ r1*tmean + r4*Complexity
  Overlap  ~ b1*tmean + b2*ppt + b4*Complexity + b5*Range
  rich ~ c1*tmean + q1*tmean_sq + c2*ppt + c4*Complexity + d1*Overlap + d2*Range

  # indirect paths to richness
  ind_tmean_Overlap     := b1*d1
  ind_ppt_Overlap       := b2*d1
  ind_Complexity_Overlap    := b4*d1

  ind_tmean_Range     := r1*d2
  ind_Complexity_Range  := r4*d2
  ind_tmean_Range_Overlap := r1*b5*d1
  ind_Complexity_Range_Overlap   := r4*b5*d1

  # temperature curvature -- the LDG test (expect q1 < 0: thermal optimum)
  curv_tmean       := q1
  # marginal dRich/dTmean at cold / mean / warm sites (std temp = -1, 0, +1)
  slope_tmean_cold := c1 + 2*q1*(-1)
  slope_tmean_mean := c1
  slope_tmean_warm := c1 + 2*q1*(1)

  # total effects on richness
  tot_tmean  := c1 + b1*d1 + r1*d2 + r1*b5*d1
  tot_ppt    := c2 + b2*d1
  tot_Complexity := c4 + b4*d1 + r4*d2 + r4*b5*d1
'
lavaanPlot(model = sem(m2, data = dat, estimator = "ML"),
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           stars = c("regress"))  # Append significance stars to regressions
graph_sem(sem(m2, data = dat, estimator = "ML"))
summary(sem(m2, data = dat, estimator = "ML"), fit.measures=TRUE)
m<-sem(m2, data = dat, estimator = "ML")
summary(m, fit.measures=TRUE)


fitMeasures(m, c("chisq", "df", "pvalue", "rmsea", "cfi", "tli", "srmr"))

m3 <- ' #Take out climate
  Range  ~ r3*npp + r4*Complexity
  Overlap  ~ b3*npp + b4*Complexity + b5*Range
  rich ~ c3*npp + c4*Complexity + d1*Overlap + d2*Range

  # indirect paths to richness
  ind_npp_Overlap       := b3*d1
  ind_Complexity_Overlap    := b4*d1

  ind_npp_Range       := r3*d2
  ind_npp_Complexity  := r4*d2
  ind_npp_Range_Overlap   := r3*b5*d1
  ind_npp_Complexity_Overlap   := r4*b5*d1

  # total effects on richness
  tot_npp    := c3 + b3*d1 + r3*d2 + r3*b5*d1
  tot_Complexity := c4 + b4*d1 + r4*d2 + r4*b5*d1
'
lavaanPlot(model = sem(m3, data = dat, estimator = "ML"),
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions
graph_sem(sem(m3, data = dat, estimator = "ML"))

AIC(sem(m2, data = dat, estimator = "ML"), sem(m3, data = dat, estimator = "ML"))

m4 <-  ' #Take out range to overlap
  Range  ~ r1*tmean + r4*Complexity
  Overlap  ~ b1*tmean + b2*ppt + b4*Complexity
  rich ~ c1*tmean + q1*tmean_sq + c2*ppt + c4*Complexity + d1*Overlap + d2*Range

  # indirect paths to richness
  ind_tmean_Overlap     := b1*d1
  ind_ppt_Overlap       := b2*d1
  ind_Complexity_Overlap    := b4*d1

  ind_tmean_Range     := r1*d2
  ind_Complexity_Range  := r4*d2

  # temperature curvature -- the LDG test (expect q1 < 0: thermal optimum)
  curv_tmean       := q1
  # marginal dRich/dTmean at cold / mean / warm sites (std temp = -1, 0, +1)
  slope_tmean_cold := c1 + 2*q1*(-1)
  slope_tmean_mean := c1
  slope_tmean_warm := c1 + 2*q1*(1)

  # total effects on richness
  tot_tmean  := c1 + b1*d1 + r1*d2 
  tot_ppt    := c2 + b2*d1
  tot_Complexity := c4 + b4*d1 + r4*d2
'

lavaanPlot(model = sem(m4, data = dat, estimator = "ML"),
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions
graph_sem(sem(m4, data = dat, estimator = "ML"))

AIC(sem(m4, data = dat, estimator = "ML"), sem(m2, data = dat, estimator = "ML"))

m5 <- ' #Take out NPP & Swap direction
  Overlap  ~ b1*tmean + b2*ppt + b4*Complexity 
  Range  ~ r1*tmean + r4*Complexity + r5*Overlap
  rich ~ c1*tmean + q1*tmean_sq + c2*ppt + c4*Complexity + d1*Overlap + d2*Range

  # indirect paths to richness
  ind_tmean_Overlap     := b1*d1
  ind_ppt_Overlap       := b2*d1
  ind_Complexity_Overlap    := b4*d1

  ind_tmean_Range     := r1*d2
  ind_Complexity_Range  := r4*d2
  ind_Overlap_Range  := r5*d2
  
  ind_tmean_Overlap_Range := b1*r1*d1
  ind_Complexity_Range_Overlap   := b4*r4*d1

  # temperature curvature -- the LDG test (expect q1 < 0: thermal optimum)
  curv_tmean       := q1
  # marginal dRich/dTmean at cold / mean / warm sites (std temp = -1, 0, +1)
  slope_tmean_cold := c1 + 2*q1*(-1)
  slope_tmean_mean := c1
  slope_tmean_warm := c1 + 2*q1*(1)

  # total effects on richness
  tot_tmean  := c1 + b1*d1 + r1*d2 + b1*r1*d1
  tot_ppt    := c2 + b2*d1
  tot_Complexity := c4 + b4*d1 + r4*d2 + b4*r4*d1
'
lavaanPlot(model = sem(m5, data = dat, estimator = "ML"),
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions
graph_sem(sem(m5, data = dat, estimator = "ML"))

AIC(sem(m5, data = dat, estimator = "ML"), sem(m2, data = dat, estimator = "ML"))


models <- list(
  "1_Full"              = m1_full,
  "2_ClimComplexityOverlap"           = m2_ClimComplexityOverlap,
  "3_climOverlap"           = m3_climOverlap,
  "4_NPPComplexityOverlap"       = m4_NPPComplexityOverlap,
  "5_NPPOverlap"       = m5_NPPOverlap,
  "6_ComplexityOverlap"               = m6_ComplexityOverlap,
  "7_ClimComplex_noOverlap"               = m7_ClimComplex_noOverlap,
  "8_NPPComplexity_NoOverlap"               = m8_NPPComplexity_NoOverlap
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

lavaanPlot(model = fits$`1_Full`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`7_ClimComplex_noOverlap`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`2_ClimComplexityOverlap`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`2_ClimComplexityOverlap`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

library(tidySEM)

# Create a default graph from the fitted model
graph_sem(fits$`2_ClimComplexityOverlap`)
graph_sem(fits$`7_ClimComplex_noOverlap`)

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


# =====================================================================
# X. Effect decomposition for top model (m2): direct / indirect / total
# =====================================================================
# Parameters ----------------------------------------------------------
m<-sem(m2, data = dat, estimator = "ML")

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

fit_m2 <- sem(m2, data = dat, estimator = "ML",
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
