library(ggplot2)
library(dplyr)
library(ggplot2)
library(sf)
library(maps)
library(dplyr)
library(ggrepel)

#### READ IN DATA####
setwd("/home/aly/Beetles/BeetleBodySizeVariation")

#### Load NEON Token and Data Product ####
# Read NEON token from file
neon_token <- read.delim("~/NEON_TOKEN", header = FALSE)[1, 1]
Beetle_dpID <- "DP1.10022.001"

#### READ IN CLEAN DATA ####
all_elytra<-read.csv("./Data/BodysizeCombinedClean.csv")
all_elytra<-subset(all_elytra, yearCollected==2018 | yearCollected==2019)

#### Explore the data ####
# Get unique siteID levels ordered by latitude
site_order <- unique(all_elytra[order(all_elytra$latitude, decreasing = TRUE), "siteID"])

# Convert siteID into a factor with the correct order
all_elytra$siteID <- factor(all_elytra$siteID, levels = site_order)
all_elytra$log_dist_cm<-log10(all_elytra$cm_elytra_max_length)

library(ggpubr)
table(subset(all_elytra, siteID=="STER")$plotID)
sort(table(subset(all_elytra, siteID=="STER")$scientificName_Species))
table(subset(all_elytra, siteID=="STER")$scientificName_Species, subset(all_elytra, siteID=="STER")$plotID)

# ggarrange(nrow=4,
#           ggplot(data = subset(all_elytra, siteID=="STER"), aes(x=log_dist_cm, fill = scientificName_Species)) +
#             geom_density(alpha=0.5) +
#             #  theme(legend.position="none") +
#             facet_wrap(.~siteID, scales = "free_y"),
#           ggplot(data = subset(all_elytra, siteID=="STER"), aes(x=log_dist_cm, fill = scientificName_Species)) +
#             geom_density(alpha=0.5) +
#             #  theme(legend.position="0") +
#             facet_wrap(.~plotID, scales = "free_y"),
#           ggplot(data = subset(all_elytra_site50, siteID=="STER"), aes(x=log_dist_cm, fill = scientificName_Species)) +
#             geom_density(alpha=0.5) +
#             #  theme(legend.position="none") +
#             facet_wrap(.~siteID, scales = "free_y"),
#           ggplot(data = subset(all_elytra_plot50, siteID=="STER"), aes(x=log_dist_cm, fill = scientificName_Species)) +
#             geom_density(alpha=0.5) +
#             #  theme(legend.position="0") +
#             facet_wrap(.~plotID, scales = "free_y")
# )

#### OStats for the site level
#### Overlap Stats ####
library(Ostats)
library(neonDivData)
#### All Obs, unedited####
elytraDF<-all_elytra[,c("domainID","siteID","plotID","scientificName_Species","log_dist_cm")]

Ostats_unedited <- Ostats(traits = as.matrix(elytraDF[,'log_dist_cm', drop = FALSE]),
                        sp = factor(elytraDF$scientificName_Species),
                        plots = factor(elytraDF$plotID),
                        random_seed = 517)

as.data.frame(Ostats_unedited$overlaps_norm)
Ostats_unedited$overlaps_norm_ses

neon_location_BET<-subset(neon_location, substr(location_id, (nchar(location_id)-2), nchar(location_id))=="bet")
head(as.data.frame(Ostats_unedited))
Ostats_unedited_df<-as.data.frame(Ostats_unedited)
colnames(Ostats_unedited_df)<-c("overlaps_norm","overlaps_unnorm",
                                "overlaps_norm_ses","overlaps_norm_ses_lower","overlaps_norm_ses_upper",
                                "overlaps_norm_raw_lower","overlaps_norm_raw_upper",
                                "overlaps_unnorm_ses","overlaps_unnorm_ses_lower","overlaps_unnorm_ses_upper",
                                "overlaps_unnorm_raw_lower","overlaps_unnorm_raw_upper")
colnames(Ostats_unedited_df)
Ostats_unedited_map<-merge(neon_location_BET, Ostats_unedited_df, by.x="plotID", by.y = "row.names")
Ostats_unedited_map <- Ostats_unedited_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

#### 20 plus obs ####
all_elytra_plot50 <- all_elytra %>%
  group_by(plotID, scientificName_Species) %>%
  filter(n() >= 20) %>%
  ungroup()
all_elytra_site50 <- all_elytra %>%
  group_by(siteID, scientificName_Species) %>%
  filter(n() >= 20) %>%
  ungroup()

plot20plus_elytraDF<-all_elytra_plot50[,c("domainID","siteID","plotID","scientificName_Species","log_dist_cm")]

Ostats_20plus <- Ostats(traits = as.matrix(plot20plus_elytraDF[,'log_dist_cm', drop = FALSE]),
                         sp = factor(plot20plus_elytraDF$scientificName_Species),
                         plots = factor(plot20plus_elytraDF$plotID),
                         random_seed = 517)

as.data.frame(Ostats_20plus$overlaps_norm)
Ostats_20plus$overlaps_norm_ses

plots<-sort(unique(plot20plus_elytraDF$plotID))
plots
subplots<-as.data.frame(plots)
subplots$siteID<-substr(subplots$plots, 1, 4)
sites<-unique(subplots$siteID)
table(subplots$siteID)

head(as.data.frame(Ostats_20plus))
Ostats_20plus_df<-as.data.frame(Ostats_20plus)
colnames(Ostats_20plus_df)<-c("overlaps_norm","overlaps_unnorm",
                                "overlaps_norm_ses","overlaps_norm_ses_lower","overlaps_norm_ses_upper",
                                "overlaps_norm_raw_lower","overlaps_norm_raw_upper",
                                "overlaps_unnorm_ses","overlaps_unnorm_ses_lower","overlaps_unnorm_ses_upper",
                                "overlaps_unnorm_raw_lower","overlaps_unnorm_raw_upper")
colnames(Ostats_20plus_df)
Ostats_20plus_map<-merge(neon_location_BET, Ostats_20plus_df, by.x="plotID", by.y = "row.names")
Ostats_20plus_map <- Ostats_20plus_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

#### Average distribution shape ####
# -------------------------------------------------------------------------
# SIMULATE ADDITIONAL OBSERVATIONS FOR LOW-SAMPLE SPECIES-SITE COMBINATIONS
#
# Goal:
#   Increase all species-site combinations to a minimum sample size of 20
#   using an empirically-derived variance relationship.
#
# Assumptions:
#   1. Well-sampled species-site combinations provide a representative
#      estimate of within-group variability.
#
#   2. Typical variability is described by:
#
#         cv2_pct = 100 * variance / mean^2
#
#      where cv2_pct is the average observed value across high-sample groups.
#
#   3. Elytra lengths are positive continuous measurements, so simulated
#      values are drawn from a lognormal distribution rather than a normal
#      distribution.
#
#      This prevents impossible negative lengths and preserves the
#      observed mean-variance scaling.
# -------------------------------------------------------------------------
# Average CV² (%) estimated from well-sampled species-site combinations
typical_cvpct <- 0.487

# Convert percentage to proportion
cv2 <- typical_cvpct / 100

# Identify low-sample species-site combinations

low_n <- all_elytra %>%
  group_by(scientificName_Species, plotID) %>%
  summarise(
    n_obs = n(),
    mean_dist = mean(cm_elytra_max_length, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_obs < 20)

# -------------------------------------------------------------------------
# Simulate observations
# -------------------------------------------------------------------------
#
# For a lognormal distribution:
#
#   CV² = exp(sdlog²) - 1
#
# therefore:
#
#   sdlog = sqrt(log(1 + CV²))
#
# and
#
#   meanlog = log(mean) - sdlog² / 2
#
# This parameterization ensures that the simulated observations have:
#
#   E[X] = observed mean
#
# and approximately:
#
#   Var[X] / E[X]² = CV²
#
# -------------------------------------------------------------------------

set.seed(42)

sim_low_n <- low_n %>%
  rowwise() %>%
  mutate(n_to_add = 20 - n_obs, # Number of observations needed to reach n = 20
         sdlog = sqrt(log(1 + cv2)), # Lognormal parameters implied by the typical CV² relationship
         meanlog = log(mean_dist) - (sdlog^2 / 2),
         # Simulate additional observations
         sim_vals = list(rlnorm(n = n_to_add,
                                meanlog = meanlog,
                                sdlog = sdlog))) %>%
  unnest(cols = sim_vals) %>%
  rename(cm_elytra_max_length = sim_vals) %>%
  select(scientificName_Species,
         plotID,
         cm_elytra_max_length) %>%
  ungroup()

sim_low_n$siteID<-NA
sim_low_n$domainID<-NA

all_elytraDF<-all_elytra[,c("domainID","siteID","plotID","scientificName_Species","cm_elytra_max_length")]
head(all_elytraDF)
head(sim_low_n)

all_elytra_aug <- rbind(all_elytraDF, sim_low_n[,colnames(all_elytraDF)])
all_elytra_aug$log_dist_cm<-log10(all_elytra_aug$cm_elytra_max_length)

Ostats_aug <- Ostats(traits = as.matrix(all_elytra_aug[,'log_dist_cm', drop = FALSE]),
                         sp = factor(all_elytra_aug$scientificName_Species),
                         plots = factor(all_elytra_aug$plotID),
                         random_seed = 517)

hist(Ostats_20plus$overlaps_norm)
hist(Ostats_aug$overlaps_norm)

head(as.data.frame(Ostats_aug))
Ostats_aug_df<-as.data.frame(Ostats_aug)
colnames(Ostats_aug_df)<-c("overlaps_norm","overlaps_unnorm",
                              "overlaps_norm_ses","overlaps_norm_ses_lower","overlaps_norm_ses_upper",
                              "overlaps_norm_raw_lower","overlaps_norm_raw_upper",
                              "overlaps_unnorm_ses","overlaps_unnorm_ses_lower","overlaps_unnorm_ses_upper",
                              "overlaps_unnorm_raw_lower","overlaps_unnorm_raw_upper")
colnames(Ostats_aug_df)
Ostats_aug_map<-merge(neon_location_BET, Ostats_aug_df, by.x="plotID", by.y = "row.names")
Ostats_aug_map <- Ostats_aug_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)


#### Average distribution shape 3-19 augment ####
low_n <- all_elytra %>%
  group_by(scientificName_Species, plotID) %>%
  summarise(
    n_obs = n(),
    mean_dist = mean(cm_elytra_max_length, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_obs < 20 & 
           n_obs>=3)

set.seed(42)

sim_low_n <- low_n %>%
  rowwise() %>%
  mutate(n_to_add = 20 - n_obs, # Number of observations needed to reach n = 20
         sdlog = sqrt(log(1 + cv2)), # Lognormal parameters implied by the typical CV² relationship
         meanlog = log(mean_dist) - (sdlog^2 / 2),
         # Simulate additional observations
         sim_vals = list(rlnorm(n = n_to_add,
                                meanlog = meanlog,
                                sdlog = sdlog))) %>%
  unnest(cols = sim_vals) %>%
  rename(cm_elytra_max_length = sim_vals) %>%
  select(scientificName_Species,
         plotID,
         cm_elytra_max_length) %>%
  ungroup()

sim_low_n$siteID<-NA
sim_low_n$domainID<-NA

all_elytra_aug3to19 <- rbind(all_elytraDF, sim_low_n[,colnames(all_elytraDF)])
all_elytra_aug3to19$log_dist_cm<-log10(all_elytra_aug3to19$cm_elytra_max_length)

Ostats_aug3to19 <- Ostats(traits = as.matrix(all_elytra_aug3to19[,'log_dist_cm', drop = FALSE]),
                     sp = factor(all_elytra_aug3to19$scientificName_Species),
                     plots = factor(all_elytra_aug3to19$plotID),
                     random_seed = 517)


head(as.data.frame(Ostats_aug3to19))
Ostats_aug3to19_df<-as.data.frame(Ostats_aug3to19)
colnames(Ostats_aug3to19_df)<-c("overlaps_norm","overlaps_unnorm",
                              "overlaps_norm_ses","overlaps_norm_ses_lower","overlaps_norm_ses_upper",
                              "overlaps_norm_raw_lower","overlaps_norm_raw_upper",
                              "overlaps_unnorm_ses","overlaps_unnorm_ses_lower","overlaps_unnorm_ses_upper",
                              "overlaps_unnorm_raw_lower","overlaps_unnorm_raw_upper")
colnames(Ostats_aug3to19_df)
Ostats_aug3to19_map<-merge(neon_location_BET, Ostats_aug3to19_df, by.x="plotID", by.y = "row.names")
Ostats_aug3to19_map <- Ostats_aug3to19_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

####Effect size Null model difference
Ostats_unedited_map$ses_direction<-ifelse(Ostats_unedited_map$overlaps_norm_ses<Ostats_unedited_map$overlaps_norm_ses_lower, 
                                          paste0("lower"), ifelse(
                                            Ostats_unedited_map$overlaps_norm_ses>Ostats_unedited_map$overlaps_norm_ses_upper, 
                                            paste0("higher"), paste0("neutral")))
table(Ostats_unedited_map$ses_direction, exclude = NULL)
Ostats_unedited_map$ses_diff<-(Ostats_unedited_map$overlaps_norm_ses-Ostats_unedited_map$overlaps_norm_ses_lower)
hist(Ostats_unedited_map$ses_diff)


Ostats_aug_map$ses_direction<-ifelse(Ostats_aug_map$overlaps_norm_ses<Ostats_aug_map$overlaps_norm_ses_lower, 
                                     paste0("lower"), ifelse(
                                       Ostats_aug_map$overlaps_norm_ses>Ostats_aug_map$overlaps_norm_ses_upper, 
                                       paste0("higher"), paste0("neutral")))
table(Ostats_aug_map$ses_direction, exclude = NULL)
Ostats_aug_map$ses_diff<-(Ostats_aug_map$overlaps_norm_ses-Ostats_aug_map$overlaps_norm_ses_lower)
hist(Ostats_aug_map$ses_diff)

Ostats_aug3to19_map$ses_direction<-ifelse(Ostats_aug3to19_map$overlaps_norm_ses<Ostats_aug3to19_map$overlaps_norm_ses_lower, 
                                          paste0("lower"), ifelse(
                                            Ostats_aug3to19_map$overlaps_norm_ses>Ostats_aug3to19_map$overlaps_norm_ses_upper, 
                                            paste0("higher"), paste0("neutral")))
table(Ostats_aug3to19_map$ses_direction, exclude = NULL)
Ostats_aug3to19_map$ses_diff<-(Ostats_aug3to19_map$overlaps_norm_ses-Ostats_aug3to19_map$overlaps_norm_ses_lower)
hist(Ostats_aug3to19_map$ses_diff)

Ostats_20plus_map$ses_direction<-ifelse(Ostats_20plus_map$overlaps_norm_ses<Ostats_20plus_map$overlaps_norm_ses_lower, 
                                        paste0("lower"), ifelse(
                                          Ostats_20plus_map$overlaps_norm_ses>Ostats_20plus_map$overlaps_norm_ses_upper, 
                                          paste0("higher"), paste0("neutral")))
table(Ostats_20plus_map$ses_direction, exclude = NULL)
Ostats_20plus_map$ses_diff<-(Ostats_20plus_map$overlaps_norm_ses-Ostats_20plus_map$overlaps_norm_ses_lower)
hist(Ostats_20plus_map$ses_diff)

#### Plots ####

plots<-sort(unique(all_elytra_aug$plotID))
plots
subplots<-as.data.frame(plots)
subplots$siteID<-substr(subplots$plots, 1, 4)
sites<-unique(subplots$siteID)
table(subplots$siteID)


for (i in 1:length(sites)) {
  png(paste0("./Figures/Overlap/Plots/PlotOverlap",sites[i],"Unedited.png"), units = "in", width = 9, height = 5, res=300)
  print(Ostats_plot(plots = all_elytra$plotID,
                    sp = all_elytra$scientificName_Species,
                    traits = all_elytra$log_dist_cm,
                    overlap_dat = Ostats_unedited,
                    use_plots = subset(subplots, siteID==sites[i])$plots,
                    name_x = 'log10(Elytra Length (cm))',
                    means = FALSE,
                    n_col = 5,
                    scale = "free_y"))
  dev.off()
}

for (i in 1:length(sites)) {
  png(paste0("./Figures/Overlap/Plots/PlotOverlap",sites[i],"20Plus.png"), units = "in", width = 9, height = 5, res=300)
  print(Ostats_plot(plots = plot20plus_elytraDF$plotID,
                    sp = plot20plus_elytraDF$scientificName_Species,
                    traits = plot20plus_elytraDF$log_dist_cm,
                    overlap_dat = Ostats_20plus,
                    use_plots = subset(subplots, siteID==sites[i])$plots,
                    name_x = 'log10(Elytra Length (cm))',
                    means = FALSE,
                    n_col = 5,
                    scale = "free_y"))
  dev.off()
}

for (i in 1:length(sites)){
  png(paste0("./Figures/Overlap/Plots/PlotOverlap",sites[i],"Augmented.png"), units = "in", width = 9, height = 5, res=300)
  print(Ostats_plot(plots = all_elytra_aug$plotID,
                    sp = all_elytra_aug$scientificName_Species,
                    traits = all_elytra_aug$log_dist_cm,
                    overlap_dat = Ostats_aug,
                    use_plots = subset(subplots, siteID==sites[i])$plots,
                    name_x = 'log10(Elytra Length (cm))',
                    means = FALSE,
                    n_col = 5,
                    scale = "free_y"))
  dev.off()
}

for (i in 1:length(sites)){
  png(paste0("./Figures/Overlap/Plots/PlotOverlap",sites[i],"Augmented3to19.png"), units = "in", width = 9, height = 5, res=300)
  print(Ostats_plot(plots = all_elytra_aug3to19$plotID,
                    sp = all_elytra_aug3to19$scientificName_Species,
                    traits = all_elytra_aug3to19$log_dist_cm,
                    overlap_dat = Ostats_aug3to19,
                    use_plots = subset(subplots, siteID==sites[i])$plots,
                    name_x = 'log10(Elytra Length (cm))',
                    means = FALSE,
                    n_col = 5,
                    scale = "free_y"))
  dev.off()
}


# Convert state map to dataframe
library(tigris)
states_map <- states(cb = TRUE, year = 2024)

dataDescrip<-c("Unedited Data","Augmented 1 to 19", "Augmented 3 to 19", "Subset >20")
dataframes<-c("Ostats_unedited_map","Ostats_aug_map","Ostats_aug3to19_map","Ostats_20plus")
datacodes<-c("Unedited","Aug1to19","Aug3to19","20plus")


i<-3
plot_dat <- cbind(Ostats_aug3to19_map,
                  st_coordinates(Ostats_aug3to19_map))

xlim <- range(plot_dat$X, na.rm = TRUE) + c(-11, 5)
ylim <- range(plot_dat$Y, na.rm = TRUE) + c(-5, 5)

png(paste0("./Figures/Overlap/Plots/Maps/",datacodes[i],"OverlapMap.png"), units = "in", width = 11, height = 6, res=300)
ggplot() +
  geom_sf(data = states_map,
          fill = "gray95",
          color = "gray70",
          linewidth = 0.2) +
  geom_point(data = plot_dat,
             aes(X, Y, color = overlaps_norm),
             position = position_jitter(width = 1.2, height = 1.2),
             size = 4, shape = 1, stroke = 1) +
  scale_color_viridis_c(name = "Overlap", option = "magma") +
  coord_sf(
    xlim = xlim,
    ylim = ylim,
    expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        legend.position = "right") +
  labs(title = "Plot Locations",
       subtitle = paste0("Points colored by overlap (",dataDescrip[i],")"),
       x = NULL,
       y = NULL)
dev.off()

png(paste0("./Figures/Overlap/Plots/Maps/",datacodes[i],"EffectSizeMap.png"), units = "in", width = 11, height = 6, res=300)
ggplot() +
  geom_sf(data = states_map,
          fill = "gray95",
          color = "gray70",
          linewidth = 0.2) +
  geom_point(data = subset(plot_dat, ses_diff>=0),
             aes(X, Y), fill = "pink",
             position = position_jitter(width = 1.2, height = 1,2),
             size = 4, shape = 21) +
  geom_point(data = subset(plot_dat, ses_diff<0),
             aes(X, Y, color = ses_diff),
             position = position_jitter(width = 1.2, height = 1,2),
             size = 4, shape = 1, stroke = .8) +
  geom_point(data = subset(plot_dat, is.na(ses_diff)),
             aes(X, Y), colour = "gray60",
             position = position_jitter(width = 1.2, height = 1,2),
             size = 2, shape = 16) +
  coord_sf(xlim = xlim,
           ylim = ylim,
           expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        legend.position = "right") +
  labs(title = "Plot Locations",
       subtitle = paste0("Points colored by Effect Size - CI Lower Bound \n (",dataDescrip[i],")"),
       x = NULL,
       y = NULL,
       color = "Effect Size \n Difference")
dev.off()

#### Write out data 
write.csv(Ostats_unedited_map, "./Outputs/Plot_Ostats_unedited.csv", row.names = FALSE)
write.csv(Ostats_aug_map, "./Outputs/Plot_Ostats_augmented.csv", row.names = FALSE)
write.csv(Ostats_aug3to19_map, "./Outputs/Plot_Ostats_augmented3to19.csv", row.names = FALSE)
write.csv(Ostats_20plus_map, "./Outputs/Plot_Ostats_20plus.csv", row.names = FALSE)

#Comparisions
Ostats_merge<-merge(Ostats_unedited_map, as.data.frame(Ostats_20plus$overlaps_norm), by.x="plotID", by.y=0, all.x=TRUE)
dim(Ostats_merge)
colnames(Ostats_merge)<-c(colnames(Ostats_merge)[1:20],"Overlap_20plus", "geometry")
colnames(Ostats_merge)
Ostats_merge<-merge(Ostats_merge, as.data.frame(Ostats_aug$overlaps_norm), by.x="plotID", by.y=0, all.x=TRUE)
Ostats_merge<-merge(Ostats_merge, as.data.frame(Ostats_aug3to19$overlaps_norm), by.x="plotID", by.y=0, all.x=TRUE)
colnames(Ostats_merge)
colnames(Ostats_merge)<-c(colnames(Ostats_merge)[1:21], "Overlap_aug1to19", "Overlap_aug3to19", "geometry")

plot(Ostats_merge$overlaps_norm~Ostats_merge$Overlap_20plus)
plot(Ostats_merge$overlaps_norm~Ostats_merge$Overlap_aug3to19)
plot(Ostats_merge$overlaps_norm~Ostats_merge$Overlap_aug1to19)

plot(Ostats_merge$Overlap_20plus~Ostats_merge$overlaps_norm)
plot(Ostats_merge$Overlap_20plus~Ostats_merge$Overlap_aug1to19)

png(paste0("./Figures/Overlap/SubsetComparison.png"), units = "in", width = 10, height = 5, res=300)
ggarrange(ggplot(data=Ostats_merge, aes(x=Overlap_20plus, y=overlaps_norm)) +
            geom_point() +
            labs(x="Overlap for 20 plus subset",
                 y="Overlap for Unedited Data") +
            xlim(0,1) +
            ylim(0,1) +
            geom_abline(slope = 1, intercept = 0),
          ggplot(data=Ostats_merge, aes(x=Overlap_20plus, y=Overlap_aug1to19)) +
            geom_point() +
            labs(x="Overlap for 20 plus subset",
                 y="Overlap for Augmented (1-19) Data") +
            xlim(0,1) +
            ylim(0,1) +
            geom_abline(slope = 1, intercept = 0),
          nrow=1)
dev.off()

png(paste0("./Figures/Overlap/UneditedComparison.png"), units = "in", width = 15, height = 5, res=300)
ggarrange(ggplot(data=Ostats_merge, aes(x=overlaps_norm, y=Overlap_20plus)) +
            geom_point() +
            labs(y="Overlap for 20 plus subset",
                 x="Overlap for Unedited Data") +
            xlim(0,1) +
            ylim(0,1) +
            geom_abline(slope = 1, intercept = 0),
          ggplot(data=Ostats_merge, aes(x=overlaps_norm, y=Overlap_aug3to19)) +
            geom_point() +
            labs(x="Overlap for Unedited Data",
                 y="Overlap for Augmented (3-19) Data") +
            xlim(0,1) +
            ylim(0,1) +
            geom_abline(slope = 1, intercept = 0),
          ggplot(data=Ostats_merge, aes(x=overlaps_norm, y=Overlap_aug1to19)) +
            geom_point() +
            labs(x="Overlap for Unedited Data",
                 y="Overlap for Augmented (1-19) Data") +
            xlim(0,1) +
            ylim(0,1) +
            geom_abline(slope = 1, intercept = 0),
          nrow=1)
dev.off()

plot(Ostats_merge$Overlap_aug1to19~Ostats_merge$Overlap_aug3to19)

table(is.na(Ostats_merge$Overlap_20plus))
table(is.na(Ostats_merge$overlaps_norm))
table(is.na(Ostats_merge$Overlap_aug1to19))
table(is.na(Ostats_merge$Overlap_aug3to19))

plot(log(Ostats_merge$overlaps_norm)~log(Ostats_merge$Overlap_20plus))
plot(log(Ostats_merge$overlaps_norm)~log(Ostats_merge$Overlap_aug1to19))
plot(log(Ostats_merge$overlaps_norm)~log(Ostats_merge$Overlap_aug3to19))
plot(log(Ostats_merge$Overlap_aug1to19)~log(Ostats_merge$Overlap_aug3to19))

