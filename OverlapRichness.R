library(ggplot2)
library(sf)
library(maps)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(neonDivData)

setwd("/home/aly/Beetles/BeetleBodySizeVariation")

site_overlap<-read.csv("./Outputs/Site_Ostats_augmented.csv")
plot_overlap<-read.csv("./Outputs/Plot_Ostats_augmented.csv")

site_richness<-read.csv("../BeetleBiodiversity/Site_annualVarWeightedMean_EstimatedSppRichness.csv")
site_richness$X<-NULL
plot_richness<-read.csv("../BeetleBiodiversity/plot_annualVarWeightedMean_EstimatedSppRichness.csv")
plot_richness$X<-NULL

siteDF<-merge(site_overlap, site_richness, by.x = "site_id", by.y = "Assemblage")
plotDF<-merge(plot_overlap, plot_richness, by.x = "plotID", by.y = "PlotID")

#Site####
head(siteDF)
neonDivData::neon_sites
siteDF<-merge(neonDivData::neon_sites, siteDF,  by.x = "siteID", by.y = "site_id")
siteDF$mean_annual_temperature_C<-as.numeric(substr(siteDF$mean_annual_temperature_C, 1, (nchar(siteDF$mean_annual_temperature_C)-2)))

library(psych)
pairs.panels(siteDF[,c("Latitude","mean_annual_temperature_C","mean_evelation_m","mean_canopy_height_m","mean_annual_precipitation_mm","log_dist_cm","median_richness")])
siteDF$overlap<-siteDF$log_dist_cm
siteDF$logOverlap<-log10(siteDF$log_dist_cm)
pairs.panels(siteDF[,c("Latitude","mean_annual_temperature_C","mean_evelation_m","mean_canopy_height_m","mean_annual_precipitation_mm","logOverlap","median_richness")])

#Heterogeneity: Spatial (geodiversity), temporal (growing degree days)
#Productivity: NPP
#Climate:
## Averages: Temp, Precipitation
## Extremes: Min temp cold quarter, max temp warm quarter, ppt driest quarter
#Interaction: ITV

