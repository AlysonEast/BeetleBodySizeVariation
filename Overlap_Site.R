library(ggplot2)
library(dplyr)
library(ggplot2)
library(sf)
library(maps)
library(dplyr)
library(ggrepel)

#### READ IN DATA####
setwd("/home/aly/Beetles/BeetleBodySizeVariation")

df<-read.csv("./Data/beetle_lengths_cm_reviewed_clean.csv")
allInd<-read.csv("./Data/allIndividuals.csv")

head(df)
str(df)
biorepo<-merge(df, allInd, by.x="beetle_id", by.y = "individualID", all.x = TRUE)
table(biorepo$flag)
biorepo<-subset(biorepo, is.na(flag))

#Update with new data, this is data specific. 
biorepo <- biorepo %>% 
  group_by(beetle_id) %>% 
  slice_head(n = 1) %>% 
  ungroup()

biorepo$cm_elytra_max_length<-biorepo$length_cm

#

biorepo$individualID<-biorepo$beetle_id
biorepo$datasource<-"Biorepo"

sp_count <- biorepo %>%
  group_by(scientific_name) %>%
  summarise(n = n()) %>%
  ungroup()

#Add in Beetlepalooza annotations
paloozadf<-read.csv("./Data/BeetleMeasurements.csv")
meta<-read.csv("./Data/NEON_Field_Site_Metadata_20260130.csv")

paloozadf<-merge(paloozadf, meta, by.x="siteID", by.y="site_id", all.x=TRUE)
df_sub<-subset(paloozadf, user_name=="IsaFluck")

ElytraLength<-subset(df_sub, structure=="ElytraLength")
colnames(ElytraLength)
ElytraLength$cm_elytra_max_length<-ElytraLength$dist_cm
ElytraLength$datasource<-"Beetle Palooza"
head(ElytraLength)

library(tidyr)
library(dplyr)
library(stringr)
library(ggpubr)

ElytraLength <- ElytraLength %>%
  mutate(
    sample_date = str_extract(NEON_sampleID, "\\d{8}"),
    collectDate = as.Date(sample_date, format = "%Y%m%d"),
    yearCollected = as.integer(substr(sample_date, 1, 4))
  )
#Add in Hawaii data ####
HI<-read.csv("./Data/trait_annotations.csv")
colnames(HI)
HI$datasource<-"PUUM"

#### Load NEON Token and Data Product ####
# Read NEON token from file
neon_token <- read.delim("~/NEON_TOKEN", header = FALSE)[1, 1]
Beetle_dpID <- "DP1.10022.001"

# If already combined data exists, load it, else fetch and combine
if (file.exists("./Data/NEON_ExpertParaCombined.csv")) {
  combined_data <- read.csv("./Data/NEON_ExpertParaCombined.csv")
} else {
  neon_df <- neonUtilities::loadByProduct(
    dpID = Beetle_dpID,
    token = neon_token,
    site = "PUUM",
    include.provisional = FALSE,
    check.size = FALSE
  )
  
  neon_para <- neon_df$bet_parataxonomistID
  neon_expert <- neon_df$bet_expertTaxonomistIDProcessed
  
  neon_para_clean <- neon_para %>%
    filter(!(individualID %in% neon_expert$individualID))
  
  common_cols <- intersect(names(neon_para_clean), names(neon_expert))
  neon_para_common <- neon_para_clean[, common_cols]
  neon_expert_common <- neon_expert[, common_cols]
  
  neon_para_common$ID_status <- "Para"
  neon_expert_common$ID_status <- "Expert"
  
  combined_data <- bind_rows(neon_para_common, neon_expert_common)
  combined_data$numbericID <- as.numeric(substr(combined_data$individualID, 
                                                (nchar(combined_data$individualID) - 5), 
                                                nchar(combined_data$individualID)))
  
  write.csv(combined_data, "./Data/NEON_ExpertParaCombined.csv", row.names = FALSE)
}

if (file.exists("./Data/NEON_ExpertParaCombined_Prelim.csv")) {
  combined_data_prelim <- read.csv("./Data/NEON_ExpertParaCombined_Prelim.csv")
} else {
  neon_df <- neonUtilities::loadByProduct(
    dpID = Beetle_dpID,
    token = neon_token,
    include.provisional = TRUE,
    site = "PUUM",
    check.size = FALSE
  )
  
  neon_para <- neon_df$bet_parataxonomistID
  neon_expert <- neon_df$bet_expertTaxonomistIDProcessed
  
  neon_para_clean <- neon_para %>%
    filter(!(individualID %in% neon_expert$individualID))
  
  common_cols <- intersect(names(neon_para_clean), names(neon_expert))
  neon_para_common <- neon_para_clean[, common_cols]
  neon_expert_common <- neon_expert[, common_cols]
  
  neon_para_common$ID_status <- "Para"
  neon_expert_common$ID_status <- "Expert"
  
  combined_data_prelim <- bind_rows(neon_para_common, neon_expert_common)
  combined_data_prelim$numbericID <- as.numeric(substr(combined_data_prelim$individualID, 
                                                       (nchar(combined_data_prelim$individualID) - 5), 
                                                       nchar(combined_data_prelim$individualID)))
  
  write.csv(combined_data_prelim, "./Data/NEON_ExpertParaCombined_Prelim.csv", row.names = FALSE)
}

combined_data_prelim<-subset(combined_data_prelim, release=="PROVISIONAL")

dim(combined_data)
combined_data<-rbind(combined_data, combined_data_prelim)
dim(combined_data)

# Add year to combined_data
combined_data$yearCollected <- as.numeric(substr(combined_data$collectDate, 1, 4))

# Extract genus and species only from scientific name
combined_data$scientificName_Species<-gsub(r"{\s*\([^\)]+\)}","",as.character(combined_data$scientificName))
combined_data$scientificName_Species<-gsub(" {2,}", " ", combined_data$scientificName_Species)
combined_data$scientificName_Species<-sub("^(\\S*\\s+\\S+).*", "\\1", combined_data$scientificName_Species)

combined_data$scientificName_Species <-
  gsub("/.*$", "", combined_data$scientificName_Species)

#Add metadata to PUUM dataset
HI_meta<-merge(HI, combined_data, by="individualID", all.x = TRUE)

# Merge all data
EL_harm <- ElytraLength %>%
  transmute(
    datasource,
    siteID,
    plotID,
    domainID = domain_id,
    individualID = combinedID,   # probably better than NEON_sampleID
    sampleID = NEON_sampleID,
    scientificName,
    collectDate,
    yearCollected,
    imageID = pictureID,
    cm_elytra_max_length,
    latitude
  )

HI_harm <- HI_meta %>%
  transmute(
    datasource,
    siteID,
    plotID,
    domainID,
    individualID,
    scientificName,
    imageID = groupImageFilePath,
    px_scalebar,
    cm_scalebar,
    px_elytra_max_length,
    cm_elytra_max_length,
    yearCollected
  )

biorepo_harm  <- biorepo %>%
  transmute(
    datasource,
    siteID,
    plotID,
    domainID,
    individualID,
    scientificName,
    imageID,
    px_scalebar = NA, ##Update!!!!
    cm_scalebar = NA, ##Update!!!!
    px_elytra_max_length = NA, ##Update!!!!
    cm_elytra_max_length,
    yearCollected
  )

#Add in biorepo
all_elytra <- bind_rows(EL_harm, HI_harm, biorepo_harm)

#### Clean the data####
all_elytra<-subset(all_elytra, !cm_elytra_max_length<=0)

# Extract genus and species only from scientific name
all_elytra$scientificName_Species<-gsub(r"{\s*\([^\)]+\)}","",as.character(all_elytra$scientificName))
all_elytra$scientificName_Species<-gsub(" {2,}", " ", all_elytra$scientificName_Species)
all_elytra$scientificName_Species<-sub("^(\\S*\\s+\\S+).*", "\\1", all_elytra$scientificName_Species)

all_elytra$scientificName_Species <-
  gsub("/.*$", "", all_elytra$scientificName_Species)
all_elytra<-subset(all_elytra, !is.na(scientificName_Species))

all_elytra <- all_elytra %>%
  filter(!grepl("sp\\.", scientificName_Species))

all_elytra<-subset(all_elytra, yearCollected == 2018 | yearCollected == 2019)


#Remove afer manual review
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D09.000975") # Wrong Spp, prob individudula ID error
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D07.013169") # Bad Measure
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D10.016322") # Check in Nathan outputs
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D08.002280") # Check in Nathan outputs
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D10.015494") # Wrong Spp, prob individudula ID error
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D09.002980") # Check in Nathan outputs
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D01.005349") # Check in Nathan outputs
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D13.001261") # Wrong Spp, prob individudula ID error
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D13.001271") # Wrong Spp, prob individudula ID error
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D07.003475") # Check in Nathan outputs
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D09.002671") # Check in Nathan outputs
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D09.002187") # Check in Nathan outputs
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D05.003147") # Wrong Spp, prob individudula ID error

all_elytra<-subset(all_elytra, individualID!="NEON.BET.D18.000854") # Check in Nathan outputs
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D18.001484") # Check in Nathan outputs
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D18.002816") # Check in Nathan outputs

all_elytra<-subset(all_elytra, individualID!="NEON.BET.D07.012343") # Check in Nathan outputs
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D07.012065") # Check in Nathan outputs
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D08.003658") # Check in Nathan outputs

all_elytra<-subset(all_elytra, individualID!="NEON.BET.D07.012011") # #Damaged Specimen
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D07.007024") # #Damaged Specimen

all_elytra$cm_elytra_max_length<-ifelse(all_elytra$imageID=="group_images/IMG_0510.png", all_elytra$cm_elytra_max_length*2, all_elytra$cm_elytra_max_length)   # Scale Bar Doubled


all_elytra<-subset(all_elytra, imageID!="MLBS_009.S.20180522.jpg") # Beetlepalooza data, one image, all outliers
all_elytra<-subset(all_elytra, imageID!="MLBS_009.E.20180522.CARABIDS.01.jpg") # Beetlepalooza data, one image, all outliers
all_elytra<-subset(all_elytra, imageID!="Cicindela_punctulata-Btray-Y2022-NEON.BET.D13.001489-NEON.BET.D13.001511.png") # all indivID yield Calathus advena

#### Explore the data ####
# Get unique siteID levels ordered by latitude
site_order <- unique(all_elytra[order(all_elytra$latitude, decreasing = TRUE), "siteID"])

# Convert siteID into a factor with the correct order
all_elytra$siteID <- factor(all_elytra$siteID, levels = site_order)
all_elytra$log_dist_cm<-log10(all_elytra$cm_elytra_max_length)

ggplot(data = all_elytra, aes(x=log_dist_cm, fill = scientificName_Species)) +
  geom_density(alpha=0.5) +
  theme(legend.position="none") +
  facet_wrap(.~siteID, scales = "free_y")

library(ggpubr)
ggarrange(nrow=2,
  ggplot(data = subset(all_elytra, siteID=="STER"), aes(x=log_dist_cm, fill = scientificName_Species)) +
    geom_density(alpha=0.5) +
  #  theme(legend.position="none") +
    facet_wrap(.~siteID, scales = "free_y"),
  ggplot(data = subset(all_elytra, siteID=="STER"), aes(x=log_dist_cm, fill = scientificName_Species)) +
    geom_density(alpha=0.5) +
  #  theme(legend.position="0") +
    facet_wrap(.~plotID, scales = "free_y")
)
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
#### All Obs, unedited####
elytraDF<-all_elytra[,c("domainID","siteID","plotID","scientificName_Species","log_dist_cm")]

Ostats_unedited <- Ostats(traits = as.matrix(elytraDF[,'log_dist_cm', drop = FALSE]),
                        sp = factor(elytraDF$scientificName_Species),
                        plots = factor(elytraDF$siteID),
                        random_seed = 517)

as.data.frame(Ostats_unedited$overlaps_norm)
Ostats_unedited$overlaps_norm_ses

meta
Ostats_unedited_df<-as.data.frame(Ostats_unedited)
str(Ostats_unedited)
colnames(Ostats_unedited_df)<-c("overlaps_norm","overlaps_unnorm",
                                "overlaps_norm_ses","overlaps_norm_ses_lower","overlaps_norm_ses_upper",
                                "overlaps_norm_raw_lower","overlaps_norm_raw_upper",
                                "overlaps_unnorm_ses","overlaps_unnorm_ses_lower","overlaps_unnorm_ses_upper",
                                "overlaps_unnorm_raw_lower","overlaps_unnorm_raw_upper")
colnames(Ostats_unedited_df)
Ostats_unedited_map<-merge(meta, Ostats_unedited_df, by.x="site_id", by.y = "row.names")
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

site20plus_elytraDF<-all_elytra_site50[,c("domainID","siteID","plotID","scientificName_Species","log_dist_cm")]

Ostats_20plus <- Ostats(traits = as.matrix(site20plus_elytraDF[,'log_dist_cm', drop = FALSE]),
                         sp = factor(site20plus_elytraDF$scientificName_Species),
                         plots = factor(site20plus_elytraDF$siteID),
                         random_seed = 517)

as.data.frame(Ostats_20plus$overlaps_norm)
Ostats_20plus$overlaps_norm_ses

sites<-levels(site20plus_elytraDF$siteID)


head(as.data.frame(Ostats_20plus))
Ostats_20plus_df<-as.data.frame(Ostats_20plus)
colnames(Ostats_20plus_df)<-c("overlaps_norm","overlaps_unnorm",
                                "overlaps_norm_ses","overlaps_norm_ses_lower","overlaps_norm_ses_upper",
                                "overlaps_norm_raw_lower","overlaps_norm_raw_upper",
                                "overlaps_unnorm_ses","overlaps_unnorm_ses_lower","overlaps_unnorm_ses_upper",
                                "overlaps_unnorm_raw_lower","overlaps_unnorm_raw_upper")
colnames(Ostats_20plus_df)
Ostats_20plus_map<-merge(meta, Ostats_20plus_df, by.x="site_id", by.y = "row.names")
Ostats_20plus_map <- Ostats_20plus_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
#### Average distribution shape 1-19 augment ####
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
typical_cvpct <- 0.506

# Convert percentage to proportion
cv2 <- typical_cvpct / 100

# Identify low-sample species-site combinations

low_n <- all_elytra %>%
  group_by(scientificName_Species, siteID) %>%
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
         siteID,
         cm_elytra_max_length) %>%
  ungroup()

sim_low_n$plotID<-NA
sim_low_n$domainID<-NA

all_elytraDF<-all_elytra[,c("domainID","siteID","plotID","scientificName_Species","cm_elytra_max_length")]
head(all_elytraDF)
head(sim_low_n)

all_elytra_aug <- rbind(all_elytraDF, sim_low_n[,colnames(all_elytraDF)])
all_elytra_aug$log_dist_cm<-log10(all_elytra_aug$cm_elytra_max_length)

Ostats_aug <- Ostats(traits = as.matrix(all_elytra_aug[,'log_dist_cm', drop = FALSE]),
                         sp = factor(all_elytra_aug$scientificName_Species),
                         plots = factor(all_elytra_aug$siteID),
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
Ostats_aug_map<-merge(meta, Ostats_aug_df, by.x="site_id", by.y = "row.names")
Ostats_aug_map <- Ostats_aug_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

#### Average distribution shape 3-19 augment ####
low_n <- all_elytra %>%
  group_by(scientificName_Species, siteID) %>%
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
         siteID,
         cm_elytra_max_length) %>%
  ungroup()

sim_low_n$plotID<-NA
sim_low_n$domainID<-NA

all_elytra_aug3to19 <- rbind(all_elytraDF, sim_low_n[,colnames(all_elytraDF)])
all_elytra_aug3to19$log_dist_cm<-log10(all_elytra_aug3to19$cm_elytra_max_length)

Ostats_aug3to19 <- Ostats(traits = as.matrix(all_elytra_aug3to19[,'log_dist_cm', drop = FALSE]),
                     sp = factor(all_elytra_aug3to19$scientificName_Species),
                     plots = factor(all_elytra_aug3to19$siteID),
                     random_seed = 517)


head(as.data.frame(Ostats_aug3to19))
Ostats_aug3to19_df<-as.data.frame(Ostats_aug3to19)
colnames(Ostats_aug3to19_df)<-c("overlaps_norm","overlaps_unnorm",
                              "overlaps_norm_ses","overlaps_norm_ses_lower","overlaps_norm_ses_upper",
                              "overlaps_norm_raw_lower","overlaps_norm_raw_upper",
                              "overlaps_unnorm_ses","overlaps_unnorm_ses_lower","overlaps_unnorm_ses_upper",
                              "overlaps_unnorm_raw_lower","overlaps_unnorm_raw_upper")
colnames(Ostats_aug3to19_df)
Ostats_aug3to19_map<-merge(meta, Ostats_aug3to19_df, by.x="site_id", by.y = "row.names")
Ostats_aug3to19_map <- Ostats_aug3to19_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

####Effect size Null model difference
Ostats_unedited_map$norm_ses_direction<-ifelse(Ostats_unedited_map$overlaps_norm_ses<Ostats_unedited_map$overlaps_norm_ses_lower, 
                                          paste0("lower"), ifelse(
                                            Ostats_unedited_map$overlaps_norm_ses>Ostats_unedited_map$overlaps_norm_ses_upper, 
                                            paste0("higher"), paste0("neutral")))
table(Ostats_unedited_map$norm_ses_direction, exclude = NULL)
Ostats_unedited_map$norm_ses_diff<-(Ostats_unedited_map$overlaps_norm_ses-Ostats_unedited_map$overlaps_norm_ses_lower)
hist(Ostats_unedited_map$norm_ses_diff)

Ostats_unedited_map$unnorm_ses_direction<-ifelse(Ostats_unedited_map$overlaps_unnorm_ses<Ostats_unedited_map$overlaps_unnorm_ses_lower, 
                                               paste0("lower"), ifelse(
                                                 Ostats_unedited_map$overlaps_unnorm_ses>Ostats_unedited_map$overlaps_unnorm_ses_upper, 
                                                 paste0("higher"), paste0("neutral")))
table(Ostats_unedited_map$unnorm_ses_direction, exclude = NULL)
Ostats_unedited_map$unnorm_ses_diff<-(Ostats_unedited_map$overlaps_unnorm_ses-Ostats_unedited_map$overlaps_unnorm_ses_lower)
hist(Ostats_unedited_map$unnorm_ses_diff)



Ostats_aug_map$norm_ses_direction<-ifelse(Ostats_aug_map$overlaps_norm_ses<Ostats_aug_map$overlaps_norm_ses_lower, 
                                     paste0("lower"), ifelse(
                                       Ostats_aug_map$overlaps_norm_ses>Ostats_aug_map$overlaps_norm_ses_upper, 
                                       paste0("higher"), paste0("neutral")))
table(Ostats_aug_map$norm_ses_direction, exclude = NULL)
Ostats_aug_map$norm_ses_diff<-(Ostats_aug_map$overlaps_norm_ses-Ostats_aug_map$overlaps_norm_ses_lower)
hist(Ostats_aug_map$norm_ses_diff)

Ostats_aug_map$unnorm_ses_direction<-ifelse(Ostats_aug_map$overlaps_unnorm_ses<Ostats_aug_map$overlaps_unnorm_ses_lower, 
                                                 paste0("lower"), ifelse(
                                                   Ostats_aug_map$overlaps_unnorm_ses>Ostats_aug_map$overlaps_unnorm_ses_upper, 
                                                   paste0("higher"), paste0("neutral")))
table(Ostats_aug_map$unnorm_ses_direction, exclude = NULL)
Ostats_aug_map$unnorm_ses_diff<-(Ostats_aug_map$overlaps_unnorm_ses-Ostats_aug_map$overlaps_unnorm_ses_lower)
hist(Ostats_aug_map$unnorm_ses_diff)

Ostats_aug3to19_map$norm_ses_direction<-ifelse(Ostats_aug3to19_map$overlaps_norm_ses<Ostats_aug3to19_map$overlaps_norm_ses_lower, 
                                          paste0("lower"), ifelse(
                                            Ostats_aug3to19_map$overlaps_norm_ses>Ostats_aug3to19_map$overlaps_norm_ses_upper, 
                                            paste0("higher"), paste0("neutral")))
table(Ostats_aug3to19_map$norm_ses_direction, exclude = NULL)
Ostats_aug3to19_map$norm_ses_diff<-(Ostats_aug3to19_map$overlaps_norm_ses-Ostats_aug3to19_map$overlaps_norm_ses_lower)
hist(Ostats_aug3to19_map$norm_ses_diff)

Ostats_aug3to19_map$unnorm_ses_direction<-ifelse(Ostats_aug3to19_map$overlaps_unnorm_ses<Ostats_aug3to19_map$overlaps_unnorm_ses_lower, 
                                            paste0("lower"), ifelse(
                                              Ostats_aug3to19_map$overlaps_unnorm_ses>Ostats_aug3to19_map$overlaps_unnorm_ses_upper, 
                                              paste0("higher"), paste0("neutral")))
table(Ostats_aug3to19_map$unnorm_ses_direction, exclude = NULL)
Ostats_aug3to19_map$unnorm_ses_diff<-(Ostats_aug3to19_map$overlaps_unnorm_ses-Ostats_aug3to19_map$overlaps_unnorm_ses_lower)
hist(Ostats_aug3to19_map$unnorm_ses_diff)

Ostats_20plus_map$norm_ses_direction<-ifelse(Ostats_20plus_map$overlaps_norm_ses<Ostats_20plus_map$overlaps_norm_ses_lower, 
                                        paste0("lower"), ifelse(
                                          Ostats_20plus_map$overlaps_norm_ses>Ostats_20plus_map$overlaps_norm_ses_upper, 
                                          paste0("higher"), paste0("neutral")))
table(Ostats_20plus_map$norm_ses_direction, exclude = NULL)
Ostats_20plus_map$norm_ses_diff<-(Ostats_20plus_map$overlaps_norm_ses-Ostats_20plus_map$overlaps_norm_ses_lower)
hist(Ostats_20plus_map$norm_ses_diff)

Ostats_20plus_map$unnorm_ses_direction<-ifelse(Ostats_20plus_map$overlaps_unnorm_ses<Ostats_20plus_map$overlaps_unnorm_ses_lower, 
                                                 paste0("lower"), ifelse(
                                                   Ostats_20plus_map$overlaps_unnorm_ses>Ostats_20plus_map$overlaps_unnorm_ses_upper, 
                                                   paste0("higher"), paste0("neutral")))
table(Ostats_20plus_map$unnorm_ses_direction, exclude = NULL)
Ostats_20plus_map$unnorm_ses_diff<-(Ostats_20plus_map$overlaps_unnorm_ses-Ostats_20plus_map$overlaps_unnorm_ses_lower)
hist(Ostats_20plus_map$unnorm_ses_diff)

#### Plots ####

#Single Site
png("./Figures/Overlap/GRSM_OstatsUnedited.png", units = "in", width = 20, height = 11, res=300)
Ostats_plot(plots = elytraDF$siteID, 
            sp = elytraDF$scientificName_Species, 
            traits = elytraDF$log_dist_cm, 
            use_plots = "GRSM", 
            name_x = 'log10(Elytra Length (cm))', 
            means = FALSE,
            legend = TRUE,
            n_col = 6,
            scale = "free_y")
dev.off()

WREF

png("./Figures/Overlap/GRSMplots_OstatsUnedited.png", units = "in", width = 20, height = 11, res=300)
Ostats_plot(plots = elytraDF$plotID, 
            sp = elytraDF$scientificName_Species, 
            traits = elytraDF$log_dist_cm, 
            use_plots = c("GRSM_001", "GRSM_006", "GRSM_008", "GRSM_012", "GRSM_013",
                          "GRSM_014", "GRSM_020" ,"GRSM_021", "GRSM_022", "GRSM_024"), 
            name_x = 'log10(Elytra Length (cm))', 
            means = FALSE,
            legend = FALSE,
            scale = "free_y",
            n_col = 5)
dev.off()

png("./Figures/Overlap/WREF_009_OstatsUnedited.png", units = "in", width = 20, height = 11, res=300)
Ostats_plot(plots = elytraDF$plotID, 
            sp = elytraDF$scientificName_Species, 
            traits = elytraDF$log_dist_cm, 
            use_plots = c("WREF_009"), 
            name_x = 'log10(Elytra Length (cm))', 
            means = FALSE,
            legend = FALSE,
            scale = "free_y",
            n_col = 5)
dev.off()

png("./Figures/Overlap/GRSMplots_OstatsUnedited.png", units = "in", width = 20, height = 11, res=300)
Ostats_plot(plots = elytraDF$siteID, 
            sp = elytraDF$scientificName_Species, 
            traits = elytraDF$log_dist_cm, 
            use_plots = c("BARR","OSBS"), 
            name_x = 'log10(Elytra Length (cm))', 
            means = FALSE,
            legend = FALSE,
            scale = "free_y",
            limits_x = c(0.5, 1),
            n_col = 5)
dev.off()

png("./Figures/Overlap/OstatsUnedited.png", units = "in", width = 20, height = 11, res=300)
Ostats_plot(plots = elytraDF$siteID, 
            sp = elytraDF$scientificName_Species, 
            traits = elytraDF$log_dist_cm, 
            overlap_dat = Ostats_unedited, 
            use_plots = sites, 
            name_x = 'Elytra Length (cm)', 
            means = FALSE,
            n_col = 6,
            scale = "free_y",
            limits_x = c(0.5, 0.75)
) 
dev.off()

png("./Figures/Overlap/Ostats20Plus.png", units = "in", width = 20, height = 11, res=300)
Ostats_plot(plots = site20plus_elytraDF$siteID, 
            sp = site20plus_elytraDF$scientificName_Species, 
            traits = site20plus_elytraDF$log_dist_cm, 
            overlap_dat = Ostats_20plus, 
            use_plots = sites, 
            name_x = 'Elytra Length (cm)', 
            means = FALSE,
            n_col = 6,
            scale = "free_y",
            limits_x = c(0.5, 0.75)
            ) 
dev.off()

png("./Figures/Overlap/OstatsAugmented.png", units = "in", width = 20, height = 11, res=300)
Ostats_plot(plots = all_elytra_aug$siteID, 
            sp = all_elytra_aug$scientificName_Species, 
            traits = all_elytra_aug$log_dist_cm, 
            overlap_dat = Ostats_aug, 
            use_plots = sort(sites)[c(1:12,14,16:47)], 
            name_x = 'Elytra Length (cm)', 
            means = FALSE,
            n_col = 6,
            scale = "free_y",
            limits_x = c(0.5, 0.75)
) 
dev.off()

png("./Figures/Overlap/OstatsAugmented3tp19.png", units = "in", width = 20, height = 11, res=300)
Ostats_plot(plots = all_elytra_aug3to19$siteID, 
            sp = all_elytra_aug3to19$scientificName_Species, 
            traits = all_elytra_aug3to19$log_dist_cm, 
            overlap_dat = Ostats_aug3to19, 
            use_plots = sort(sites)[c(1:12,14,16:47)], 
            name_x = 'Elytra Length (cm)', 
            means = FALSE,
            n_col = 6,
            scale = "free_y",
            limits_x = c(0.5, 0.75)
            ) 
dev.off()

# par(mfrow=c(1,2))
# png("./Figures/OverlapExample_20plus.png", height = 3, width = 9, units = "in", res = 300)
Ostats_plot(plots = all_elytraDF$siteID,
            sp = all_elytraDF$scientificName_Species,
            traits = all_elytraDF$log_dist_cm,
            overlap_dat = Ostats_20plus,
            use_plots = c("SERC","SRER"),
            name_x = 'Elytra Length (cm)',
            means = TRUE,
            n_col = 6,
            scale = "free",
            limits_x = c(0,1))
# dev.off()
# 
# png("./Figures/OverlapExample_augmented.png", height = 3, width = 9, units = "in", res = 300)
Ostats_plot(plots = all_elytra_aug$siteID,
            sp = all_elytra_aug$scientificName_Species,
            traits = all_elytra_aug$log_dist_cm,
            overlap_dat = Ostats_aug,
            use_plots = c("SERC","SRER"),
            name_x = 'Elytra Length (cm)',
            means = TRUE,
            n_col = 6,
            scale = "free",
            limits_x = c(0,1))
# dev.off()

plot_dat <- Ostats_unedited_map %>%
  arrange(overlaps_unnorm_ses) %>%
  mutate(site = factor(row_number(), levels = row_number()))

ggplot(plot_dat,
       aes(y = reorder(site, overlaps_unnorm_ses),
           x = overlaps_unnorm_ses)) +
  
  geom_errorbarh(aes(xmin = overlaps_unnorm_ses_lower,
                     xmax = overlaps_unnorm_ses_upper),
                 height = 0) +
  
  geom_point(size = 2,
             colour = "red")

#### Map of overlaps ####
# Convert state map to dataframe
states_map <- map_data("state")

ggplot() + geom_polygon(data = states_map, aes(x = long, y = lat, group = group),   # Base map
                        fill = "gray95",
                        color = "gray70",
                        linewidth = 0.2) +
  geom_sf(data = Ostats_unedited_map, aes(color = overlaps_norm),   # Site points
          size = 6, alpha = 0.9, inherit.aes = FALSE) +
  geom_text_repel(data = Ostats_unedited_map, aes(label = site_id, geometry = geometry), # Site labels
                  stat = "sf_coordinates", size = 3,
                  min.segment.length = 0, segment.color = NA,
                  seed = 42, inherit.aes = FALSE) +
  scale_color_viridis_c(name = "Overlap", option = "magma") +
  coord_sf() +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        legend.position = "right") +
  labs(title = "Site Locations",
       subtitle = "Points colored by Unedited overlap",
       x = NULL,
       y = NULL)

ggplot() + geom_polygon(data = states_map, aes(x = long, y = lat, group = group),   # Base map
                        fill = "gray95",
                        color = "gray70",
                        linewidth = 0.2) +
  geom_sf(data = subset(Ostats_aug_map, site_id!="ONAQ"), aes(color = overlaps_norm),   # Site points
          size = 6, alpha = 0.9, inherit.aes = FALSE) +
  geom_text_repel(data = subset(Ostats_aug_map), aes(label = site_id, geometry = geometry), # Site labels
                  stat = "sf_coordinates", size = 3,
                  min.segment.length = 0, segment.color = NA,
                  seed = 42, inherit.aes = FALSE) +
  scale_color_viridis_c(name = "Overlap", option = "magma") +
  coord_sf() +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        legend.position = "right") +
  labs(title = "Site Locations",
       subtitle = "Points colored by Augmented overlap",
       x = NULL,
       y = NULL)

# Convert state map to dataframe
library(tigris)
states_map <- states(cb = TRUE, year = 2024)

dataDescrip<-c("Unedited Data","Augmented 1 to 19", "Augmented 3 to 19", "Subset >20")
dataframes<-c("Ostats_unedited_map","Ostats_aug_map","Ostats_aug3to19_map","Ostats_20plus_map")
datacodes<-c("Unedited","Aug1to19","Aug3to19","20plus")


overlap_max<-max(Ostats_unedited_map$overlaps_norm, Ostats_aug_map$overlaps_norm, Ostats_aug3to19_map$overlaps_norm, Ostats_20plus_map$overlaps_norm, na.rm = TRUE)

i<-4
plot_dat <- cbind(Ostats_20plus_map,
                  st_coordinates(Ostats_20plus_map))

xlim <- range(plot_dat$X, na.rm = TRUE) + c(-11, 5)
ylim <- range(plot_dat$Y, na.rm = TRUE) + c(-5, 5)

png(paste0("./Figures/Overlap/Sites/Maps/",datacodes[i],"OverlapMap.png"), units = "in", width = 11, height = 6, res=300)
ggplot() +
  geom_sf(data = states_map,
          fill = "gray95",
          color = "gray70",
          linewidth = 0.2) +
  geom_point(data = plot_dat,
             aes(X, Y, color = overlaps_norm),
             # position = position_jitter(width = 1.2, height = 1.2),
             size = 4) +
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

png(paste0("./Figures/Overlap/Sites/Maps/",datacodes[i],"EffectSizeMap.png"), units = "in", width = 11, height = 6, res=300)
ggplot() +
  geom_sf(data = states_map,
          fill = "gray95",
          color = "gray70",
          linewidth = 0.2) +
  geom_point(data = subset(plot_dat, ses_diff<0),
             aes(X, Y, color = ses_diff),
             # position = position_jitter(width = 1.2, height = 1,2),
             size = 4) +
  geom_point(data = subset(plot_dat, is.na(ses_diff)),
             aes(X, Y), colour = "gray60",
             # position = position_jitter(width = 1.2, height = 1,2),
             ) +
  geom_point(data = subset(plot_dat, ses_diff>=0),
             aes(X, Y), color = "pink",
             # position = position_jitter(width = 1.2, height = 1,2),
             size = 4) +
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
write.csv(Ostats_unedited_map, "./Outputs/Site_Ostats_unedited.csv", row.names = FALSE)
write.csv(Ostats_aug_map, "./Outputs/Site_Ostats_augmented.csv", row.names = FALSE)
write.csv(Ostats_aug3to19_map, "./Outputs/Site_Ostats_augmented3to19.csv", row.names = FALSE)
write.csv(Ostats_20plus_map, "./Outputs/Site_Ostats_20plus.csv", row.names = FALSE)

Ostats_unedited_map<-read.csv("./Outputs/Site_Ostats_unedited.csv")

dim(Ostats_unedited_map)
dim(Ostats_aug_map)
dim(Ostats_aug3to19_map)
dim(Ostats_20plus_map)

Ostats_merge<-merge(Ostats_unedited_map, as.data.frame(Ostats_20plus$overlaps_norm), by.x="site_id", by.y=0, all.x=TRUE)
dim(Ostats_merge)
colnames(Ostats_merge)<-c(colnames(Ostats_merge)[1:52], "Overlap_20plus", "geometry")
colnames(Ostats_merge)
Ostats_merge<-merge(Ostats_merge, as.data.frame(Ostats_aug$overlaps_norm), by.x="site_id", by.y=0, all.x=TRUE)
Ostats_merge<-merge(Ostats_merge, as.data.frame(Ostats_aug3to19$overlaps_norm), by.x="site_id", by.y=0, all.x=TRUE)
colnames(Ostats_merge)
colnames(Ostats_merge)<-c(colnames(Ostats_merge)[1:53], "Overlap_aug1to19", "Overlap_aug3to19", "geometry")

plot(Ostats_merge$overlaps_norm~Ostats_merge$Overlap_20plus)
plot(Ostats_merge$overlaps_norm~Ostats_merge$Overlap_aug1to19)
plot(Ostats_merge$overlaps_norm~Ostats_merge$Overlap_aug3to19)
plot(Ostats_merge$Overlap_aug1to19~Ostats_merge$Overlap_aug3to19)


plot(log(Ostats_merge$overlaps_norm)~log(Ostats_merge$Overlap_20plus))
plot(log(Ostats_merge$overlaps_norm)~log(Ostats_merge$Overlap_aug1to19))
plot(log(Ostats_merge$overlaps_norm)~log(Ostats_merge$Overlap_aug3to19))
plot(log(Ostats_merge$Overlap_aug1to19)~log(Ostats_merge$Overlap_aug3to19))

Ostats_merge$diff_og_20plus<-Ostats_merge$overlaps_norm-Ostats_merge$Overlap_20plus
Ostats_merge$diff_og_aug1to19<-Ostats_merge$overlaps_norm-Ostats_merge$Overlap_aug1to19
Ostats_merge$diff_og_aug3to19<-Ostats_merge$overlaps_norm-Ostats_merge$Overlap_aug3to19
Ostats_merge$diff_aug1to19_aug3to19<-Ostats_merge$Overlap_aug1to19-Ostats_merge$Overlap_aug3to19

hist(Ostats_merge$diff_og_20plus)
hist(Ostats_merge$diff_og_aug1to19)
hist(Ostats_merge$diff_og_aug3to19)
