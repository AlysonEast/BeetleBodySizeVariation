library(ggplot2)
library(sf)
library(maps)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(neonDivData)

setwd("/home/aly/Beetles/BeetleBodySizeVariation")

geodiv_dir<-"/media/aly/Penobscot/NEON/Geodiversity/edi.2320.1/"

#Read in and merge overlap and richness data
site_overlap<-read.csv("./Outputs/site_by_domain_PoolNull.csv")
site_overlap$latitude<-NULL
plot_overlap<-read.csv("./Outputs/plot_by_site_PoolNull.csv")
plot_overlap$latitude<-NULL

site_richness<-read.csv("../BeetleBiodiversity/Site_annualVarWeightedMean_EstimatedSppRichness.csv")
site_richness$X<-NULL
plot_richness<-read.csv("../BeetleBiodiversity/plot_annualVarWeightedMean_EstimatedSppRichness.csv")
plot_richness$X<-NULL

siteDF<-merge(site_overlap, site_richness, by.x = "siteID", by.y = "Assemblage")
plotDF<-merge(plot_overlap, plot_richness, by.x = "plotID", by.y = "PlotID")

#Site####
head(siteDF)
siteDF$domainID<-NULL
neonDivData::neon_sites
siteDF<-merge(neonDivData::neon_sites, siteDF,  by = "siteID")

geodiv<-read.csv(paste0(geodiv_dir,"NEON_site_footprint_elev30m.csv"))
head(geodiv)
geodiv$domainID<-NULL

siteDF<-merge(siteDF, geodiv, by="siteID")
head(siteDF)

library(psych)
pairs.panels(siteDF[,c("Latitude","bio01_mean","srtm_mean","bio12_mean","overlap_unnorm_obs","median_richness")])

#### Plot ####
#Enviromental Data
struc<-read.csv("./Outputs/BETplot_Rugosity.csv")
struc$X<-NULL
env<-read.csv("./Outputs/BeetlePlotswEnvData.csv")
NPP<-read.csv("../NEON_MODIS_NPP_2018_2019.csv") #from https://code.earthengine.google.com/b41a55076352b2d9e21ac5e74bf337bc


plotDF<-merge(plotDF, struc, by="plotID")
plotDF<-merge(plotDF, env, by="plotID")
plotDF<-merge(plotDF, NPP[,c("Npp","Gpp","plotID")], by="plotID")
head(plotDF)

####Paths####
#Heterogeneity: Spatial (geodiversity), temporal (growing degree days)
#Productivity: NPP
#Climate:
## Averages: Temp, Precipitation
## Extremes: Min temp cold quarter, max temp warm quarter, ppt driest quarter
#Interaction: ITV

