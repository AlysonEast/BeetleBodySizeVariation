library(ggplot2)
library(sf)
library(maps)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(neonDivData)
library(sf)
library(psych)
library(geodata)
library(terra)

setwd("/home/aly/Beetles/BeetleBodySizeVariation")

BETplots<-read_sf("../All_NEON_TOS_Plot_Polygons_V11_BET.shp")
BETpts<-read_sf("../All_NEON_TOS_Plot_Centroids_V11_BET.shp")

plot(BETplots)
head(BETplots)
plot(BETpts)
head(BETpts)



# Download WorldClim bioclim variables
bio <- worldclim_global(country = "USA",
                         var = "bio",
                         res = 0.5,      # 0.5 arc-min (~900 km)
                         path = "/media/aly/Penobscot/Data/WorldClim/")

names(bio)

BETpts<-st_transform(BETpts, st_crs(bio))

climate <- extract(bio, BETpts)
head(climate)
colnames(climate)<-substr(colnames(climate),11,20)
colnames(climate)[1]<-"ID"
colnames(climate)

plots_climate <- cbind(BETpts, climate)

# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (×100)
# BIO4 = Temperature Seasonality (standard deviation ×100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter

climate_matrix<-as.data.frame(plots_climate[,c("bio_1","bio_2","bio_4","bio_5","bio_6","bio_12","bio_15","bio_16","bio_17","bio_18")])
climate_matrix$geometry<-NULL

climate_matrix<-as.data.frame(plots_climate[,c(44:63)])
climate_matrix$geometry<-NULL

pairs.panels(climate_matrix)

# climate_matrix$bio_12<-log10(climate_matrix$bio_12)
# climate_matrix$bio_16<-log10(climate_matrix$bio_16)
# climate_matrix$bio_17<-log10(climate_matrix$bio_17)

pairs.panels(climate_matrix) 

clim_pca <- princomp(climate_matrix, cor = TRUE)
summary(clim_pca)
(clim_pca$sdev^2 / sum(clim_pca$sdev^2))

plot(clim_pca)

loadings(clim_pca)
loadings(clim_pca)[,1:2]
scores(clim_pca)

library("corrplot")
corrplot(clim_pca$loadings, is.corr=FALSE, col.lim=c(-1,1))

pts_pca<-cbind(plots_climate, clim_pca$scores)
plot(pts_pca, max.plot=100)

NPP<-read.csv("../NEON_MODIS_NPP_2018_2019.csv") #from https://code.earthengine.google.com/b41a55076352b2d9e21ac5e74bf337bc
head(NPP)

pts_gpp<-merge(pts_pca, NPP[,c("plotID","Npp")], by="plotID")

write.csv(pts_pca, "./Outputs/BeetlePlotswEnvData.csv", row.names = FALSE)

### Plot Level LiDAR data ####
library(neonUtilities)
library(neonOS)

#Get LiDAR data as needed ####
NEON_TOKEN<-read.delim("~/NEON_TOKEN",header = FALSE)[1,1]

product<-c("DP3.30024.001","DP1.30003.001")
p<-NA #Make this 1 or 2 to download
sites<-unique(BETpts$siteID)
for (i in seq_along(sites)) {
  tryCatch({
    byTileAOP(
      dpID = product[p], # "DP1.30003.001"
      site = sites[i],
      easting = subset(BETpts, siteID == sites[i])$easting,
      northing = subset(BETpts, siteID == sites[i])$northing,
      year = 2018,
      token = NEON_TOKEN,
      savepath = "/media/aly/Penobscot/NEON/LiDAR/",
      check.size = FALSE,
      progress = TRUE
    )
  }, error = function(e) {
    if (grepl("There are no data at the selected site and year", e$message)) {
      message("No LiDAR for ", sites[i], " in 2018; skipping.")} 
    else {
      stop(e)}
  })
}

files<-list.files(paste0("/media/aly/Penobscot/NEON/LiDAR/",product[p],"/neon-aop-products"), recursive = TRUE)
sites_done <- sapply(strsplit(files, "/"), `[`, 4)
sites_done <- sapply(strsplit(sites_done, "_"), `[`, 2)
sites_done <- unique(sites_done)
length(sort(sites_done))
sort(sites_done)

sites2<-sort(symdiff(sites_done, sites))
length(sites2)
sites2

for (i in seq_along(sites2)) {
  tryCatch({
    byTileAOP(
      dpID = product[p],
      site = sites2[i],
      easting = subset(BETpts, siteID == sites2[i])$easting,
      northing = subset(BETpts, siteID == sites2[i])$northing,
      year = 2019,
      token = NEON_TOKEN,
      savepath = "/media/aly/Penobscot/NEON/LiDAR/",
      check.size = FALSE,
      progress = TRUE
    )
  }, error = function(e) {
    if (grepl("There are no data at the selected site and year", e$message)) {
      message("No LiDAR for ", sites[i], " in 2019; skipping.")} 
    else {
      stop(e)}
  })
}

files<-list.files("/media/aly/Penobscot/NEON/LiDAR/DP1.30003.001/neon-aop-products", recursive = TRUE)
sites_done <- sapply(strsplit(files, "/"), `[`, 4)
sites_done <- sapply(strsplit(sites_done, "_"), `[`, 2)
sites_done <- unique(sites_done)
length(sort(sites_done))
sort(sites_done)

sites3<-sort(symdiff(sites_done, sites))
length(sites3)
sites3

for (i in seq_along(sites3)) {
  tryCatch({
    byTileAOP(
      dpID = product[p],
      site = sites3[i],
      easting = subset(BETpts, siteID == sites3[i])$easting,
      northing = subset(BETpts, siteID == sites3[i])$northing,
      year = 2017,
      token = NEON_TOKEN,
      savepath = "/media/aly/Penobscot/NEON/LiDAR/",
      check.size = FALSE,
      progress = TRUE
    )
  }, error = function(e) {
    if (grepl("There are no data at the selected site and year", e$message)) {
      message("No LiDAR for ", sites[i], " in 2017; skipping.")} 
    else {
      stop(e)}
  })
}

files<-list.files("/media/aly/Penobscot/NEON/LiDAR/DP1.30003.001/neon-aop-products", recursive = TRUE)
sites_done <- sapply(strsplit(files, "/"), `[`, 4)
sites_done <- sapply(strsplit(sites_done, "_"), `[`, 2)
sites_done <- unique(sites_done)
length(sort(sites_done))
sort(sites_done)
#DCFS is part of the flight box for WOOD. Downloading data from WOOD
#KONA is part of the flight box for KONZ. Downloading data from KONZ
#TREE is part of the flight box for STEI. Downloading data from STEI

filesDF <- as.data.frame(sapply(strsplit(files, "/"), `[`, 4))
colnames(filesDF)<-"siteYearCode"
filesDF <- cbind(filesDF, sapply(strsplit(files, "/"), `[`, 1))
filesDF$site <- substr(filesDF$siteYearCode, 6, 9)
table(filesDF$`sapply(strsplit(files, "/"), \`[\`, 1)`, filesDF$site)


# calculate metrics ####
library(lidR)
library(gstat)
library(purrr)
files<-list.files("/media/aly/Penobscot/NEON/LiDAR/DP1.30003.001/neon-aop-products", recursive = TRUE, full.names = TRUE)
DTMfiles<-list.files("/media/aly/Penobscot/NEON/LiDAR/DP3.30024.001/neon-aop-products", recursive = TRUE, full.names = TRUE, pattern="DTM.tif$")

rugosity_results <- list()

sites_list<-sort(unique(BETpts$siteID))
AOP_list<-sites_list
AOP_list[AOP_list == "DCFS"] <- "WOOD"
AOP_list[AOP_list == "KONA"] <- "KONZ"
AOP_list[AOP_list == "STEI"] <- "CHEQ" #Why is STEI downloading CHEQ AOP without any notice?
AOP_list[AOP_list == "TREE"] <- "STEI" #And TREE downloads as STEI (with warning)

sites_list
AOP_list

for (i in 1:length(sites_list)) {#
  # Get LiDAR files for site
  site_files <- grep(
    AOP_list[i],
    files,
    value = TRUE)
  
  ctg <- readLAScatalog(site_files)
  
  site_dtm <- grep(
    AOP_list[i],
    DTMfiles,
    value=TRUE)
  
  dtm <- vrt(site_dtm)

  # Plot centroids for site
  pts<-as.data.frame(subset(BETpts, siteID==sites_list[i]))
  
  # UTM projection
  zone <- as.numeric(gsub("N","", unique(pts$utmZone)))
  epsg <- 32600 + zone
  pts_utm <- st_as_sf(pts,
                      coords = c("easting","northing"),
                      crs = as.numeric(epsg))
  if(sites_list[i]=="BLAN") {
    pts_utm<-st_transform(pts_utm, crs = "epsg:32617")
  }
  
  buffers <- st_buffer(pts_utm, 40)

  #-----------------------------
  # Loop through plots
  #-----------------------------
  
  site_results <- map_dfr(
    1:nrow(buffers),
    function(j){
      print(paste("Processing  Plot", pts$plotID[j]))
      # extract lidar only in 40m radius
      las <- clip_roi(ctg,
                      buffers[j,])
      if(is.empty(las)){
        return(
          data.frame(
            plotID=pts$plotID[j],
            siteID=sites_list[i],
            rugosity_RC=NA,
            mean_height=NA,
            sd_height=NA
          )
        )
      }
      
      dtm_plot <- crop(dtm,
                       vect(buffers[j,]))

    # Remove noise
    las <- filter_poi(las, Classification != 7)
    las_norm <- normalize_height(las, dtm_plot)

    # 1 m voxel structural density
    vox <- voxel_metrics(las,
                         ~length(Z),
                         res = 0.5)
    
    names(vox)[names(vox)=="V1"] <- "density"
    
    # Calculate RC
    # SD of vertical density profile
    # within each x column
    column_sd <- vox %>%
      group_by(X,Y) %>%
      summarise(sd_z = 
                  if(n() > 1)
                    sd(density)
                else
                  0,
                .groups="drop")
    
    rugosity_RC <- sd(column_sd$sd_z,
                      na.rm=TRUE)
    

    # Additional structure metrics
    mean_height <- mean(las_norm@data$Z, na.rm=TRUE)
    height_98 <- quantile(las_norm@data$Z, na.rm=TRUE, probs = 0.98)
    
    sd_height <- sd(las_norm@data$Z,na.rm=TRUE)
    
    data.frame(plotID = pts$plotID[j],
               siteID = sites_list[i],
               rugosity_RC = rugosity_RC,
               height_98 = height_98,
               mean_height = mean_height,
               sd_height = sd_height,
               n_points = npoints(las_norm))
})
  rugosity_results[[i]] <- site_results
  }

rugosity_results <- bind_rows(rugosity_results)
  
write.csv(rugosity_results, "./Outputs/BETplot_Rugosity")
