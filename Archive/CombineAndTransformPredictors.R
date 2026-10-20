library(ggplot2)
library(ggpubr)
library(neonDivData)
library(lavaan)
library(dplyr)
library(psych)

setwd("/home/aly/Beetles/BeetleBodySizeVariation")

## PLOT ####

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
pairs.panels(plotDF[,c("bio_1","bio_12","rugosity_RC","Npp","abund","overlap_norm_obs","richness")])
plotDF$log_rugosity_RC<-log10((plotDF$rugosity_RC+0.01))
plotDF$log_abund<-log10(plotDF$abund)
plotDF$log_overlap_unnorm_obs<-log10(plotDF$overlap_norm_obs)
plotDF$sqrt_overlap_unnorm_obs<-sqrt(plotDF$overlap_norm_obs)
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


#columns to use
Temp <-"bio_1 second order"
Precip <- "log_bio_12"
Climate <- "Comp.1"
Rugosity <- "rugosity_RC"
NPP <- "log_Npp"
Abundance <- "log_abound"
Overlap <- "sqrt_overlap_unnorm_obs"


## SITE ####


