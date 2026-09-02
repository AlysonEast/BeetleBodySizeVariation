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

Overlap_COL    <- "sqrt_overlap_unnorm_obs" 
Complexity_COL <- "rugosity_RC"
RICH_COL   <- "median_richness"
TMEAN_COL  <- "bio_1" # Second Order Mean daily mean temperature of coldest quarter
PPT_COL    <- "log_bio_12" #Mean monthly precipitation of the driest quarter
NPP_COL  <- "log_Npp"                    
Climate <- "Comp.1"
Abundance <- "log_abound"


## Transforms (applied before standardizing)
LOG_Overlap        <- TRUE               # overlap spans many orders of magnitude -> log
RICH_TRANSFORM <- "log"             # "none", "sqrt", or "log"
STANDARDIZE    <- TRUE               # z-score all model vars (coeffs in SD units)
EXCLUDE_ISLANDS <- TRUE

## ============================================================ ##
## 1. Assemble plot data: start from plotDF, add NPP
## ============================================================ ##
#Read in and merge overlap and richness data
# plot_overlap<-read.csv("./Outputs/plot_by_all_noaug_ByYearAvg_IndividualNull.csv") #use plot_by_all becuase there are no exclusions due to domains with 1 site
#Read in overlap data
plot_2018<-read.csv("./Outputs/plot_by_all_aug_2018_IndividualNull.csv")
plot_2019<-read.csv("./Outputs/plot_by_all_aug_2019_IndividualNull.csv")
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

plotDF$diffpct<-(plotDF$Observed-plotDF$n_overlap_sp)/plotDF$Observed
plotDF$diffpct<-as.numeric(ifelse(plotDF$diffpct<0, paste0(NA), plotDF$diffpct))

table(plotDF$diffpct, useNA = "ifany")
hist(plotDF$diffpct)
plotDF$diffdouble<-ifelse(plotDF$diffpct>.5, paste0(1), paste0(0))
plotDF$diffthird<-ifelse(plotDF$diffpct>(2/3), paste0(1), paste0(0))


ggplot(plotDF, aes(x=richness, y=n_overlap_sp, colour = log(overlap_norm_obs))) +
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
plotDF<-subset(plotDF, completeness>=.5)
plotDF<-subset(plotDF, diffpct<.5 & n_overlap_sp<=2 | 
                 n_overlap_sp>2 & diffpct<=(2/3))

dim(preExclusion)
dim(plotDF)
dim(preExclusion)[1]-dim(plotDF)[1]
plotDF<-subset(plotDF, !is.na(overlap_norm_obs))

symdiff(levels(as.factor(preExclusion$plotID.x)),levels(as.factor(plotDF$plotID.x)))
length(symdiff(levels(as.factor(preExclusion$plotID.x)),levels(as.factor(plotDF$plotID.x))))
dim(table(plotDF$plotID.x))

symdiff(levels(as.factor(preExclusion$SiteID)),levels(as.factor(plotDF$SiteID)))

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

plot(plotDF$richness~plotDF$overlap_norm_obs)
plot(plotDF$n_overlap_sp~plotDF$overlap_norm_obs)


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
pairs.panels(plotDF[,c("bio_1","log_bio_12","log_rugosity_RC","log_Npp","log_abund","log_overlap_unnorm_obs","sqrt_richness","log_richness")])


if (EXCLUDE_ISLANDS) plotDF<-plotDF %>% 
  filter(!grepl('PUUM', plotDF$plotID),
         !grepl('LAJA', plotDF$plotID),
         !grepl('GUAN', plotDF$plotID)) #c("PUUM","LAJA","GUAN"))


pairs.panels(plotDF[,c("Comp.1", "log_comp.1","bio_1","log_bio_1","bio_12","log_bio_12","rugosity_RC","log_rugosity_RC",
                       "Npp","log_Npp","log_abund","abund",
                       "overlap_unnorm_obs","log_overlap_unnorm_obs","sqrt_overlap_unnorm_obs",
                       "richness","sqrt_richness","log_richness")])

pairs.panels(plotDF[,c("Comp.1", "Comp.2","Comp.3","Comp.4","Comp.5",
                       "richness","sqrt_richness","log_richness")])


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

m2_ClimComplexityOverlap <- '  
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
m3_climOverlap <- '
  Overlap  ~ b1*tmean + b2*ppt 
  rich ~ c1*tmean + c2*ppt + d*Overlap

  # indirect paths to richness
  ind_tmean_Overlap     := b1*d
  ind_ppt_Overlap       := b2*d

  # total effects on richness
  tot_tmean  := c1 + b1*d
  tot_ppt    := c2 + b2*d
'
m4_NPPComplexityOverlap <- '
  Overlap  ~ b3*npp + b4*Complexity
  rich ~ c3*npp + c4*Complexity + d*Overlap

  # indirect paths to richness
  ind_npp_Overlap       := b3*d
  ind_Complexity_Overlap    := b4*d
  
  # total effects on richness
  tot_npp    := c3 + b3*d
  tot_Complexity := c4 + b4*d
'
m5_NPPOverlap <- '
  Overlap  ~ b3*npp
  rich ~ c3*npp + d*Overlap

  # indirect paths to richness
  ind_npp_Overlap       := b3*d

  # total effects on richness
  tot_npp    := c3 + b3*d
'
m6_ComplexityOverlap <- '
  Overlap  ~ b4*Complexity
  rich ~ c4*Complexity + d*Overlap

  # indirect paths to richness
  ind_Complexity_Overlap    := b4*d
  
  # total effects on richness
  tot_Complexity := c4 + b4*d
'
m7_ClimComplex_noOverlap <- '  
  rich ~ c1*tmean + c2*ppt + c4*Complexity + d*Overlap

  # total effects on richness
  tot_tmean  := c1 
  tot_ppt    := c2 
  tot_Complexity := c4 
  tot_Overlap := d 
'
m8_NPPComplexity_NoOverlap <- '
  rich ~ c3*npp + c4*Complexity + d*Overlap
  
  # total effects on richness
  tot_npp    := c3 
  tot_Complexity := c4 
  tot_Overlap := d 
'

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

lavaanPlot(model = fits$`8_NPPComplexity_NoOverlap`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`7_ClimComplex_noOverlap`,
           coefs = TRUE,          # Display the path coefficients
           stand = TRUE,          # Standardize the coefficients
           sig = 0.05,            # Only highlight significant paths
           stars = c("regress"))  # Append significance stars to regressions

lavaanPlot(model = fits$`4_NPPComplexityOverlap`,
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
graph_sem(fits$`4_NPPComplexityOverlap`)
graph_sem(fits$`8_NPPComplexity_NoOverlap`)
graph_sem(fits$`2_ClimComplexityOverlap`)

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
