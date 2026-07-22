library(dplyr)
library(tidyr)
library(stringr)
library(parameters)
library(diptest)
library(moments)


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
    order = individual,
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
    order = BeetlePosition,
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
    order = Order,
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
no_measurment<-subset(all_elytra, cm_elytra_max_length<=0)
all_elytra<-subset(all_elytra, !cm_elytra_max_length<=0)

# Extract genus and species only from scientific name
all_elytra$scientificName_Species<-gsub(r"{\s*\([^\)]+\)}","",as.character(all_elytra$scientificName))
all_elytra$scientificName_Species<-gsub(" {2,}", " ", all_elytra$scientificName_Species)
all_elytra$scientificName_Species<-sub("^(\\S*\\s+\\S+).*", "\\1", all_elytra$scientificName_Species)

all_elytra$scientificName_Species <-
  gsub("/.*$", "", all_elytra$scientificName_Species)

all_elytra<-subset(all_elytra, !is.na(scientificName_Species))

#### Generate individuals for manual review
#---------------------------------------------------------
# 1. Species-level robust statistics
#---------------------------------------------------------

all_elytra_flagged <- all_elytra %>%
  group_by(scientificName_Species) %>%
  mutate(
    species_median = median(cm_elytra_max_length, na.rm = TRUE),
    species_mad = mad(cm_elytra_max_length,
                      constant = 1.4826,
                      na.rm = TRUE),
    
    robust_z_species =
      (cm_elytra_max_length - species_median) / species_mad,
    
    flag_species_z3 = abs(robust_z_species) > 3,
    flag_species_z4 = abs(robust_z_species) > 4,
    flag_species_z5 = abs(robust_z_species) > 5
  ) %>%
  ungroup()

#---------------------------------------------------------
# 2. Species x Site robust statistics
#---------------------------------------------------------

all_elytra_flagged <- all_elytra_flagged %>%
  group_by(scientificName_Species, siteID) %>%
  mutate(
    
    site_n = n(),
    
    site_median = median(cm_elytra_max_length, na.rm = TRUE),
    site_mad = mad(cm_elytra_max_length,
                   constant = 1.4826,
                   na.rm = TRUE),
    
    robust_z_site =
      (cm_elytra_max_length - site_median) / site_mad,
    
    flag_site = site_n >= 5 & abs(robust_z_site) > 3
  ) %>%
  ungroup()

#---------------------------------------------------------
# 3. Image-level consistency
#---------------------------------------------------------

all_elytra_flagged <- all_elytra_flagged %>%
  group_by(imageID) %>%
  mutate(
    
    image_n = n(),
    
    image_mean = mean(cm_elytra_max_length, na.rm = TRUE),
    image_sd   = sd(cm_elytra_max_length, na.rm = TRUE),
    
    image_z =
      ifelse(image_sd > 0,
             (cm_elytra_max_length - image_mean)/image_sd,
             0),
    
    flag_image = image_n >= 3 & abs(image_z) > 2
  ) %>%
  ungroup()

#---------------------------------------------------------
# 4. Species extremes (top/bottom 1%)
#---------------------------------------------------------

all_elytra_flagged <- all_elytra_flagged %>%
  group_by(scientificName_Species) %>%
  mutate(
    
    lower1 = quantile(cm_elytra_max_length,
                      0.01,
                      na.rm = TRUE),
    
    upper99 = quantile(cm_elytra_max_length,
                       0.99,
                       na.rm = TRUE),
    
    flag_extreme =
      cm_elytra_max_length < lower1 |
      cm_elytra_max_length > upper99
    
  ) %>%
  ungroup()

#---------------------------------------------------------
# 5. Overall review score
#---------------------------------------------------------

all_elytra_flagged <- all_elytra_flagged %>%
  mutate(
    
    review_score =
      1*flag_species_z3 +
      2*flag_species_z4 +
      3*flag_species_z5 +
      2*flag_site +
      1*flag_image +
      1*flag_extreme,
    
    review_flag = review_score >= 3
  )

review_list <- all_elytra_flagged %>%
  filter(review_flag) %>%
  arrange(desc(review_score),
          desc(abs(robust_z_species)))
dim(table(review_list$imageID))
#write.csv(review_list, "./Outputs/BodySizeMeasuremntsForManualReview.csv", row.names = FALSE)

#### Reconcile after manual review
reviewed<-read.csv("./Outputs/BodySizeMeasuremntsForManualReview_Reviewed.csv")
dim(reviewed)
table(reviewed$status)
552/(735-115)
7/(735-115)
5/(735-115)
26/(735-115)
30/(735-115)
dim(table(reviewed$imageID))
dim(table(subset(reviewed,status=="Reorder")$imageID))
(dim(table(subset(reviewed,status=="Reorder")$imageID))/dim(table(reviewed$imageID)))

#Deal with bad annotations
Need_annotation<-subset(reviewed, status=="Bad")
write.csv(Need_annotation, "./Outputs/BodySizeMeasuremntsForManualAnnotation.csv")

exclusion<-subset(reviewed, status=="Bad" | status=="Damaged" | status=="Missing")
exclusion<-unique(exclusion$individualID)
#remove bad annotations from dataset
all_elytra<-all_elytra %>%
  filter(!individualID %in% exclusion)
#add in manual corrections
# manual_corrections<-read.csv("./")
# cbind(all_elytra, manual_corrections)

#deal with reorders
imageExclusion<-subset(reviewed, status=="Reorder")$imageID
write.csv(imageExclusion, "./Outputs/BodySizeMeasuremntsForManualReorder.csv")
#remove bad annotations from dataset
all_elytra<-all_elytra %>%
  filter(!imageID %in% imageExclusion)
#add in manual corrections
# manual_reorder<-read.csv("./")
# cbind(all_elytra, manual_reorder)

#Remove afer manual review
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D09.000975") # Wrong Spp, prob individudula ID error
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D10.015494") # Wrong Spp, prob individudula ID error
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D13.001261") # Wrong Spp, prob individudula ID error
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D13.001271") # Wrong Spp, prob individudula ID error
all_elytra<-subset(all_elytra, individualID!="NEON.BET.D05.003147") # Wrong Spp, prob individudula ID error

#Hawaii scale bar double sized
all_elytra$cm_elytra_max_length<-ifelse(all_elytra$imageID=="group_images/IMG_0510.png", all_elytra$cm_elytra_max_length*2, all_elytra$cm_elytra_max_length)   # Scale Bar Doubled


#all_elytra<-subset(all_elytra, imageID!="MLBS_009.S.20180522.jpg") # Beetlepalooza data, one image, all outliers
#all_elytra<-subset(all_elytra, imageID!="MLBS_009.E.20180522.CARABIDS.01.jpg") # Beetlepalooza data, one image, all outliers
all_elytra<-subset(all_elytra, imageID!="Cicindela_punctulata-Btray-Y2022-NEON.BET.D13.001489-NEON.BET.D13.001511.png") # all indivID yield Calathus advena

write.csv(all_elytra, "./Data/BodysizeCombinedClean.csv", row.names = FALSE)