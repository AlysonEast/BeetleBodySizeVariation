library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(parameters)
library(diptest)
library(moments)

#User defined variables#
cutoff<-50
#

setwd("/home/aly/Beetles/BeetleBodySizeVariation")

df<-read.csv("./Data/beetle_lengths_cm_reviewed_clean.csv")
allInd<-read.csv("./Data/allIndividuals.csv")

head(df)

# df <- df %>%
#   # Split the single column into 8 numeric columns
#   separate(pred_coords_px,
#            into = paste0("v", 1:8),
#            sep = "[, ]+",   # handles comma or space separated
#            convert = TRUE) %>%
#   # Compute distance using first four values
#   mutate(
#     width = sqrt((v3 - v1)^2 + (v4 - v2)^2)
#   ) %>%
#   # Compute distance using second four values
#   mutate(
#     length = sqrt((v7 - v5)^2 + (v8 - v6)^2)
#   )

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

#Distribution of Mean Body Size in Carabids####
ElytraSummary<- all_elytra %>%
  group_by(scientificName_Species) %>%
  summarise(
    n_obs = n(),
    mean_dist = mean(cm_elytra_max_length, na.rm = TRUE),
    var_dist = var(cm_elytra_max_length, na.rm = TRUE),
    skew = skewness(cm_elytra_max_length, na.rm = TRUE),
    kurtosis = datawizard::kurtosis(cm_elytra_max_length, na.rm = TRUE)[["Kurtosis"]],
    kurtosisSE = datawizard::kurtosis(cm_elytra_max_length, na.rm = TRUE)[["SE"]],
    dstat = dip.test(cm_elytra_max_length)[1],
    dpval = dip.test(cm_elytra_max_length)[2],
    .groups = "drop"
  )

ElytraSummary <- ElytraSummary %>%
  filter(!grepl("sp\\.", scientificName_Species))

table(ElytraSummary$scientificName_Species)
str(ElytraSummary)

#png("./Figures/BodySizeQuantification/AllBodySizeDist.png", units = "in", width = 6, height = 4, res=300)
ggarrange(ggplot(data = ElytraSummary, aes(x=mean_dist)) +
            geom_histogram() +
            theme_pubr() +
            xlab("Species Mean Bodysize (cm)") +
            theme(legend.position="none"),
          ggplot(data = ElytraSummary, aes(x=log10(mean_dist))) +
            geom_histogram() +
            theme_pubr() +
            xlab("log10(Species Mean Bodysize (cm))") +
            theme(legend.position="none"),
          nrow=1)
dev.off()
#ITV###
#What is the shape of distributions?####
ElytraSummary_n<-subset(ElytraSummary, n_obs>=cutoff)
#Skewness#
#png("./Figures/BodySizeQuantification/SpeciesDistribution_Skew.png", units = "in", width = 6, height = 5, res=300)
ggplot(data = ElytraSummary_n, aes(skew)) +
  geom_histogram() + 
  theme(legend.position="none") +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = -1) +
  geom_vline(xintercept = -2) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = 2) +
  annotate("text", x = 0, y = 20, label = "Normal", angle = 90, vjust = -0.5) +
  annotate("text", x = -1, y = 20, label = "left-skewed", angle = 90, vjust = -0.5) +
  annotate("text", x = 1, y = 20, label = "right-skewed", angle = -90, vjust = -0.5) +
  annotate("text", x = 2, y = 20, label = "Exponential right-skewed", angle = -90, vjust = -0.5)+
  annotate("text", x = -2, y = 20, label = "Exponential left-skewed", angle = 90, vjust = -0.5) +
  theme_pubr() +
  labs(title = "Skew of Body Size Distribution",
       subtitle = "Species Level",
       x = "Skew",
       y = "Count of Species",
       caption = paste0("Subset to obersvations of species with >",cutoff," individuals"))
dev.off()

print(paste0("left-skewed: n = ", nrow(subset(ElytraSummary_n,skew > 1))))

print(paste0("right-skewed: n = ", nrow(subset(ElytraSummary_n,skew < -1))))

#Kurtosis#
dat <- ElytraSummary_n %>%
  arrange(kurtosis) %>%
  mutate(rank = row_number())
dat$kurtosisLower<-dat$kurtosis - 1.96*as.numeric(dat$kurtosisSE)
dat$sig<-ifelse((dat$kurtosis - 1.96*dat$kurtosisSE)>0 | (dat$kurtosis + 1.96*dat$kurtosisSE)<0,
                paste0("*"),NA)

#png("./Figures/BodySizeQuantification/SpeciesDistribution_Kurtosis.png", units = "in", width = 6, height = 5, res=300)
ggplot(dat, aes(kurtosis, rank, colour = sig)) +
  geom_errorbarh(aes(xmin = kurtosis - 1.96*kurtosisSE,
                     xmax = kurtosis + 1.96*kurtosisSE),
                 height = 0) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "Kurtosis", y = NULL) +
  # annotate("text", x = 0, y = 15, label = "Peaked", angle = -90, vjust = -0.5) +
  # annotate("text", x = 0, y = 15, label = "Flattened", angle = 90, vjust = -0.5) +
  theme_pubr() +
  labs(title = "Kurtosis of Body Size Distribution",
       subtitle = "Species Level",
       x = "Kurtosis",
       y = "Rank",
       caption = paste0("Subset to obersvations of species with >",cutoff," individuals"))+
  theme(legend.position = "None")
dev.off()

#png("./Figures/BodySizeQuantification/SpeciesDistribution_KurtosisHigh.png", units = "in", width = 3, height = 3, res=300, bg = "transparent")
ggplot(subset(all_elytra,scientificName_Species=="Amara conflata"), aes(cm_elytra_max_length)) +
  geom_histogram() +  
  geom_density() +
  theme_pubr() +
  labs(title = "Amara conflata",
       x = "Body Size (cm)",
       y = "Count of Indivudals")
dev.off()

#png("./Figures/BodySizeQuantification/SpeciesDistribution_KurtosisLow.png", units = "in", width = 3, height = 3, res=300, bg = "transparent")
ggplot(subset(all_elytra,scientificName_Species=="Cyclotrachelus constrictus"), aes(cm_elytra_max_length)) +
  geom_histogram() +  
  geom_density() +
  theme_pubr() +
  labs(title = "Cyclotrachelus constrictus",
       x = "Body Size (cm)",
       y = "Count of Indivudals")
dev.off()



#Modality
#png("./Figures/BodySizeQuantification/SpeciesDistribution_Modality.png", units = "in", width = 6, height = 5, res=300)
ggplot(data = ElytraSummary_n, aes(as.numeric(dpval))) +
  geom_histogram() +
  geom_vline(xintercept = 0.05) +
  annotate("text", x = 0.05, y = 15, label = "Multimodal", angle = 90, vjust = -0.5)+
  annotate("text", x = 0.05, y = 15, label = "Unimodial", angle = -90, vjust = -0.5) +
  theme_pubr() +
  labs(title = "Modality of Body Size Distribution",
       subtitle = "Species Level",
       x = "Dip Test P-value",
       y = "Count of Species",
       caption = paste0("Subset to obersvations of species with >",cutoff," individuals"))
dev.off()


print(paste0("Unimodal Distributions: n = ", nrow(subset(ElytraSummary_n,dpval>=0.05))))
  
print(paste0("Multimodal Distributions: n = ", nrow(subset(ElytraSummary_n,dpval<=0.05))))
subset(ElytraSummary_n,dpval<=0.05)

#Variance across Scales ####
head(all_elytra)

species_stats <- all_elytra %>%
  filter(!grepl("sp\\.", scientificName_Species)) %>%
  group_by(scientificName_Species) %>%
  summarise(
    n = n(),
    mean_cm_elytra_max_length = mean(cm_elytra_max_length, na.rm = TRUE),
    var_cm_elytra_max_length  = var(cm_elytra_max_length, na.rm = TRUE),
    sd_cm_elytra_max_length   = sd(cm_elytra_max_length, na.rm = TRUE),
    skew = skewness(cm_elytra_max_length, na.rm = TRUE),
    kurtosis = datawizard::kurtosis(cm_elytra_max_length, na.rm = TRUE)[["Kurtosis"]],
    kurtosisSE = datawizard::kurtosis(cm_elytra_max_length, na.rm = TRUE)[["SE"]],
    dstat = dip.test(cm_elytra_max_length)[1],
    dpval = dip.test(cm_elytra_max_length)[2]
  ) %>%
  ungroup()

species_plot_stats <- all_elytra %>%
  filter(!grepl("sp\\.", scientificName_Species)) %>%
  group_by(plotID, scientificName_Species) %>%
  summarise(
    n = n(),
    mean_cm_elytra_max_length = mean(cm_elytra_max_length, na.rm = TRUE),
    var_cm_elytra_max_length  = var(cm_elytra_max_length, na.rm = TRUE),
    sd_cm_elytra_max_length   = sd(cm_elytra_max_length, na.rm = TRUE),
    skew = skewness(cm_elytra_max_length, na.rm = TRUE),
    kurtosis = datawizard::kurtosis(cm_elytra_max_length, na.rm = TRUE)[["Kurtosis"]],
    kurtosisSE = datawizard::kurtosis(cm_elytra_max_length, na.rm = TRUE)[["SE"]],
    dstat = dip.test(cm_elytra_max_length)[1],
    dpval = dip.test(cm_elytra_max_length)[2]
  ) %>%
  ungroup()

species_site_stats <- all_elytra %>%
  filter(!grepl("sp\\.", scientificName_Species)) %>%
  group_by(siteID, scientificName_Species) %>%
  summarise(
    n = n(),
    mean_cm_elytra_max_length = mean(cm_elytra_max_length, na.rm = TRUE),
    var_cm_elytra_max_length  = var(cm_elytra_max_length, na.rm = TRUE),
    sd_cm_elytra_max_length   = sd(cm_elytra_max_length, na.rm = TRUE),
    skew = skewness(cm_elytra_max_length, na.rm = TRUE),
    kurtosis = datawizard::kurtosis(cm_elytra_max_length, na.rm = TRUE)[["Kurtosis"]],
    kurtosisSE = datawizard::kurtosis(cm_elytra_max_length, na.rm = TRUE)[["SE"]],
    dstat = dip.test(cm_elytra_max_length)[1],
    dpval = dip.test(cm_elytra_max_length)[2]
  ) %>%
  ungroup()

species_domain_stats <- all_elytra %>%
  filter(!grepl("sp\\.", scientificName_Species)) %>%
  group_by(domainID, scientificName_Species) %>%
  summarise(
    n = n(),
    mean_cm_elytra_max_length = mean(cm_elytra_max_length, na.rm = TRUE),
    var_cm_elytra_max_length  = var(cm_elytra_max_length, na.rm = TRUE),
    sd_cm_elytra_max_length   = sd(cm_elytra_max_length, na.rm = TRUE),
    skew = skewness(cm_elytra_max_length, na.rm = TRUE),
    kurtosis = datawizard::kurtosis(cm_elytra_max_length, na.rm = TRUE)[["Kurtosis"]],
    kurtosisSE = datawizard::kurtosis(cm_elytra_max_length, na.rm = TRUE)[["SE"]],
    dstat = dip.test(cm_elytra_max_length)[1],
    dpval = dip.test(cm_elytra_max_length)[2]
  ) %>%
  ungroup()


# png("./Figures/BodySizeQuantification/plotDistribution_Skew.png", units = "in", width = 6, height = 5, res=300)
ggplot(data = subset(species_plot_stats, n>=cutoff), aes(skew)) +
  geom_histogram() + 
  theme(legend.position="none") +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = -1) +
  geom_vline(xintercept = -2) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = 2) +
  annotate("text", x = 0, y = 40, label = "Normal", angle = 90, vjust = -0.5) +
  annotate("text", x = -1, y = 40, label = "left-skewed", angle = 90, vjust = -0.5) +
  annotate("text", x = 1, y = 40, label = "right-skewed", angle = -90, vjust = -0.5) +
  annotate("text", x = 2, y = 40, label = "Exponential right-skewed", angle = -90, vjust = -0.5)+
  annotate("text", x = -2, y = 40, label = "Exponential left-skewed", angle = 90, vjust = -0.5) +
  theme_pubr() +
  labs(title = "Skew of Body Size Distribution",
       subtitle = "Species Level",
       x = "Skew",
       y = "Count of Species",
       caption = paste0("Subset to obersvations of species with >",cutoff," individuals"))
# dev.off()

dat <- species_plot_stats %>%
  arrange(kurtosis) %>%
  mutate(rank = row_number())
dat$sig<-ifelse((dat$kurtosis - 1.96*dat$kurtosisSE)>0 | (dat$kurtosis + 1.96*dat$kurtosisSE)<0,
                paste0("*"),NA)
#png("./Figures/BodySizeQuantification/plotDistribution_Kurtosis.png", units = "in", width = 6, height = 5, res=300)
ggplot(subset(dat, n>=cutoff), aes(kurtosis, rank, colour = sig)) +
  geom_errorbarh(aes(xmin = kurtosis - 1.96*kurtosisSE,
                     xmax = kurtosis + 1.96*kurtosisSE),
                 height = 0) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "Kurtosis", y = NULL) +
  # annotate("text", x = 0, y = 15, label = "Peaked", angle = -90, vjust = -0.5) +
  # annotate("text", x = 0, y = 15, label = "Flattened", angle = 90, vjust = -0.5) +
  theme_pubr() +
  labs(title = "Kurtosis of Body Size Distribution",
       subtitle = "Plot Level",
       x = "Kurtosis",
       y = "Rank",
       caption = paste0("Subset to obersvations of species with >",cutoff," individuals"))+
  theme(legend.position = "None")
dev.off()
dat <- species_site_stats %>%
  arrange(kurtosis) %>%
  mutate(rank = row_number())
dat$sig<-ifelse((dat$kurtosis - 1.96*dat$kurtosisSE)>0 | (dat$kurtosis + 1.96*dat$kurtosisSE)<0,
                paste0("*"),NA)
#png("./Figures/BodySizeQuantification/siteDistribution_Kurtosis.png", units = "in", width = 6, height = 5, res=300)
ggplot(subset(dat, n>=cutoff), aes(kurtosis, rank, colour = sig)) +
  geom_errorbarh(aes(xmin = kurtosis - 1.96*kurtosisSE,
                     xmax = kurtosis + 1.96*kurtosisSE),
                 height = 0) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "Kurtosis", y = NULL) +
  # annotate("text", x = 0, y = 15, label = "Peaked", angle = -90, vjust = -0.5) +
  # annotate("text", x = 0, y = 15, label = "Flattened", angle = 90, vjust = -0.5) +
  theme_pubr() +
  labs(title = "Kurtosis of Body Size Distribution",
       subtitle = "Site Level",
       x = "Kurtosis",
       y = "Rank",
       caption = paste0("Subset to obersvations of species with >",cutoff," individuals"))+
  theme(legend.position = "None")
dev.off()


#Modality
#png("./Figures/BodySizeQuantification/plotDistribution_Modality.png", units = "in", width = 6, height = 5, res=300)
ggplot(data = subset(species_plot_stats, n>=cutoff), aes(as.numeric(dpval))) +
  geom_histogram() +
  geom_vline(xintercept = 0.05) +
  annotate("text", x = 0.05, y = 15, label = "Multimodal", angle = 90, vjust = -0.5)+
  annotate("text", x = 0.05, y = 15, label = "Unimodial", angle = -90, vjust = -0.5) +
  theme_pubr() +
  labs(title = "Modality of Body Size Distribution",
       subtitle = "Species Level",
       x = "Dip Test P-value",
       y = "Count of Species",
       caption = paste0("Subset to obersvations of species with >",cutoff," individuals"))
dev.off()

print(paste0("Multimodal Distributions: n = ", nrow(subset(subset(species_plot_stats, n>=cutoff),dpval<=0.05))))
subset(subset(species_plot_stats, n>=cutoff),dpval<=0.05)

ggplot(all_elytra, aes(scientificName_Species, log10(cm_elytra_max_length))) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(species_site_stats,
       aes(scientificName_Species, mean_cm_elytra_max_length, color = siteID)) +
  geom_point() +
  geom_errorbar(
    aes(
      ymin = mean_cm_elytra_max_length - sd_cm_elytra_max_length,
      ymax = mean_cm_elytra_max_length + sd_cm_elytra_max_length
    ),
    width = 0
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

species_occurrence <- all_elytra %>%
  group_by(scientificName_Species) %>%
  summarise(
    n_sites = n_distinct(siteID),
    n_plots = n_distinct(plotID),
    n_obs   = n()
  )
# species_keep <- species_occurrence %>%
#   filter(n_sites > 1, n_plots > 1) %>%
#   pull(scientificName_Species)
# 
# all_elytraDF_multi <- all_elytra %>%
#   filter(scientificName_Species %in% species_keep)

species_var <- all_elytra %>%
  filter(!grepl("sp\\.", scientificName_Species)) %>%
  group_by(scientificName_Species) %>%
  summarise(
    n = n(),
    mean_cm = mean(cm_elytra_max_length, na.rm = TRUE),
    var_cm  = var(cm_elytra_max_length, na.rm = TRUE),
    cv2_pct = 100 * var_cm / mean_cm^2,
    ID = "Spp",
    Above = "None",
    siteID = "All",
    domainID = "All",
    plotID = "All"
  ) %>%
  ungroup()

plot_species_var <- all_elytra %>%
  filter(!grepl("sp\\.", scientificName_Species)) %>%
  group_by(plotID, scientificName_Species) %>%
  summarise(
    n = n(),
    mean_cm = mean(cm_elytra_max_length, na.rm = TRUE),
    var_cm  = var(cm_elytra_max_length, na.rm = TRUE),
    cv2_pct = 100 * var_cm / mean_cm^2,
    Above = unique(siteID),
    siteID = unique(siteID),
    domainID = unique(domainID),
    plotID = unique(plotID)
  ) %>% 
  ungroup()
plot_species_var$ID<-plot_species_var$siteID
# colnames(plot_species_var)<-c("ID",colnames(plot_species_var)[2:length(colnames(plot_species_var))])

site_species_var <- all_elytra %>%
  filter(!grepl("sp\\.", scientificName_Species)) %>%
  group_by(siteID, scientificName_Species) %>%
  summarise(
    n = n(),
    mean_cm = mean(cm_elytra_max_length, na.rm = TRUE),
    var_cm  = var(cm_elytra_max_length, na.rm = TRUE),
    cv2_pct = 100 * var_cm / mean_cm^2,
    Above = unique(domainID),
    siteID = unique(siteID),
    domainID = unique(domainID),
    plotID = "All"
  ) %>%
  ungroup()
site_species_var$ID<-site_species_var$siteID
# colnames(site_species_var)<-c("ID",colnames(site_species_var)[2:length(colnames(site_species_var))])


domain_species_var <- all_elytra %>%
  filter(!grepl("sp\\.", scientificName_Species)) %>%
  group_by(domainID, scientificName_Species) %>%
  summarise(
    n = n(),
    mean_cm = mean(cm_elytra_max_length, na.rm = TRUE),
    var_cm  = var(cm_elytra_max_length, na.rm = TRUE),
    cv2_pct = 100 * var_cm / mean_cm^2,
    Above = "Spp",
    siteID = "All",
    domainID = unique(domainID),
    plotID = "All"
  ) %>%
  ungroup()
domain_species_var$ID<-domain_species_var$domainID
# colnames(domain_species_var)<-c("ID",colnames(domain_species_var)[2:length(colnames(domain_species_var))])


plot_species_var$scale <- "Plot Level"
species_var$scale <- "Species Level"
site_species_var$scale <- "Site Level"
domain_species_var$scale <- "Domain Level"

list<-plot_species_var %>% filter(n >= cutoff)
list<-list %>%
  group_by(scientificName_Species) %>%
  summarise(
    n=n()
  ) %>%
  ungroup()
list<-subset(list, n>1)
list<-unique(list$scientificName_Species)

var_all_scales <- bind_rows(
  species_var %>% filter(n >= cutoff),
  plot_species_var %>%
    filter(n >= cutoff),
  site_species_var %>%
    filter(n >= cutoff),
  domain_species_var %>%
    filter(n >= cutoff)
)

var_all_scales$scale<-factor(
  var_all_scales$scale,
  ordered = TRUE,
  levels = c("Species Level", "Domain Level", "Site Level", "Plot Level"))

dim(var_all_scales)
var_all_scales <- var_all_scales %>% filter(scientificName_Species %in% list)
dim(var_all_scales)

head(var_all_scales)

summary<-var_all_scales %>%
  group_by(scale) %>%
  summarise(
    n = n(),
    mean_cvpct = mean(cv2_pct, na.rm = TRUE),
    lower_ci = t.test(cv2_pct)$conf.int[1],
    upper_ci = t.test(cv2_pct)$conf.int[2]
  ) %>%
  ungroup()

summary
#png("./Figures/BodySizeQuantification/CV_Nested.png", units = "in", width = 7, height = 6, res=300)
ggplot(var_all_scales,
       aes(cv2_pct)) +
  theme_pubr() + 
  geom_histogram() +
  geom_vline(data = summary, aes(xintercept = mean_cvpct)) +
  geom_vline(data = summary, aes(xintercept = lower_ci), linetype="dashed") +
  geom_vline(data = summary, aes(xintercept = upper_ci), linetype="dashed") +
  scale_x_continuous(expand = c(0,0), limits = c(0,3), breaks = c(seq(0,3,0.5))) +
  xlab("Variance as % of mean² (CV² × 100)") +
  facet_wrap(.~scale, scale="free_y", ncol=1) +
  theme(legend.position = "None") +
  labs(title = "Variance in bodysize as % of mean across nested spatial scales",
       caption = paste0("Subset to obersvations of species with >",cutoff," individuals"))
dev.off()

#png("./Figures/NestedCVpct.png", height = 10, width = 10, units = "in", res = 300)
ggplot(subset(var_all_scales, scientificName_Species!="Carabidae sp."),
       aes(scale, cv2_pct, group = scientificName_Species, col = siteID, shape = domainID)) +
  theme_pubr() + 
  scale_shape_manual(values=c(16,15,17:25,6:14))+
  geom_point(size = 2) +
  # scale_color_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
  #                               "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
  #                               "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
  #                               "#8A7C64", "#599861")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Variance as % of mean² (CV² × 100)") +
  facet_wrap(.~scientificName_Species, scales = "free_y") +
  theme(legend.position = "None")
dev.off()

#Plots for single species
ggplot(subset(var_all_scales, scientificName_Species=="Pterostichus lama" | scientificName_Species=="Pterostichus mutus" | 
                scientificName_Species=="Pterostichus rostratus"),
       aes(scale, cv2_pct, group = scientificName_Species, col = siteID, shape = domainID)) +
  theme_pubr() + 
  scale_shape_manual(values=c(16,15,17:25))+
  geom_point(size = 4,
             stroke = 2) +
  scale_color_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                "#673770", "#D3D93E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
                                "#8A7C64", "#599861")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(.~scientificName_Species) +
  ylab("Variance as % of mean² (CV² × 100)") # +
  # theme(legend.position = "None")

sort(table(var_all_scales$scientificName_Species))

species<-"Carabus goryi"
plotList<-subset(plot_species_var, scientificName_Species==species &
                   n>=20)
plotList<-unique(plotList$plotID)
siteList<-subset(site_species_var, scientificName_Species==species &
                   n>=20)
siteList<-unique(siteList$siteID)
subset(site_species_var, scientificName_Species==species &
               n>=20)
domList<-subset(domain_species_var, scientificName_Species==species &
                   n>=20)
domList<-unique(domList$domainID)
domList<-"D07"

#png("./Figures/BodySizeQuantification/SpeciesNestedDensityPlots.png", units = "in", width = 4, height = 10, res=300)
annotate_figure(
  ggarrange(
    ggplot(data = subset(all_elytra, scientificName_Species==species),
           aes(x=cm_elytra_max_length, fill = domainID)) +
      geom_density(alpha=0.5, fill="#660033") +
      geom_rug(colour = "#660033") +
      geom_vline(data = subset(domain_species_var %>%
                                 filter(domainID %in% domList), 
                               scientificName_Species==species &
                                 n>=cutoff),
                 aes(xintercept = mean_cm),
                 col="#660033") +
      ylab("Domain Density") +
      xlab(NULL) +
      labs(title=species) +
      theme_pubr()+
      theme(legend.position="right"),
      # facet_wrap(.~domainID),
    ggplot(data = subset(all_elytra, scientificName_Species==species)%>%
             filter(siteID %in% siteList) %>%
             filter(domainID %in% domList),
           aes(x=cm_elytra_max_length, fill = siteID)) +
      geom_density(alpha=.5) +
      geom_rug(aes(colour = siteID)) +
      geom_vline(data = subset(site_species_var, scientificName_Species==species &
                                 n>=cutoff)%>%
                   filter(siteID %in% siteList) %>%
                   filter(domainID %in% domList),
                 aes(xintercept = mean_cm,
                     colour = siteID)) +
      scale_fill_manual(values=c("#990000", "#006699")) +
      scale_color_manual(values=c("#990000", "#006699")) +
      xlab(NULL) +
      theme_pubr()+
      theme(legend.position = "inside",
            legend.position.inside = c(0.85, 0.85)) +
      ylab("Plot Density"),
    ggplot(data = subset(all_elytra, scientificName_Species==species) %>%
             filter(plotID %in% plotList) %>%
             filter(domainID %in% domList),
           aes(x=cm_elytra_max_length, fill= plotID)) +
      geom_density(alpha=0.2) +
      geom_rug(aes(colour = plotID)) +
      geom_vline(data = subset(plot_species_var, scientificName_Species==species &
                                 n>=cutoff) %>%
                   filter(plotID %in% plotList) %>%
                   filter(domainID %in% domList),
                 aes(xintercept = mean_cm,
                     color = plotID)) +
      scale_fill_manual(values=c("#00CCCC", "#009999", "#3399cc",
                            "#3366cc", "#0000ff", "#33ccff",
                            "#0099cc", "#0066ff", "red",
                            "#660000", "#990000", "#cc3333",
                            "#ff0000", "#cc3300", "#ff6666",
                            "#ff6600", "#cc0033", "#fc4e2a",
                            "#d23428")) +
      scale_color_manual(values=c("#00CCCC", "#009999", "#3399cc",
                                 "#3366cc", "#0000ff", "#33ccff",
                                 "#0099cc", "#0066ff", "red",
                                 "#660000", "#990000", "#cc3333",
                                 "#ff0000", "#cc3300", "#ff6666",
                                 "#ff6600", "#cc0033", "#fc4e2a",
                                 "#d23428")) +
      theme_pubr() +
      theme(legend.position = "none") +
      # theme(legend.position = "inside",
      #       legend.position.inside = c(0.85, 0.85)) +
      # guides(color = guide_legend(ncol = 5),
      #        fill = guide_legend(ncol = 5)) +
      ylab("Plot Density") +
      xlab("Elytra Length (cm)"),
      # facet_wrap(.~domainID),
    nrow = 3
  )
)
dev.off()

#png("./Figures/BodySizeQuantification/SpeciesNestedWrappedDensityPlots.png", units = "in", width = 8, height = 6, res=300)
annotate_figure(
  ggarrange(
    ggplot(data = subset(all_elytra, scientificName_Species==species)%>%
             filter(siteID %in% siteList) %>%
             filter(domainID %in% domList),
           aes(x=cm_elytra_max_length, fill = siteID)) +
      geom_density(alpha=.5) +
      geom_rug(aes(colour = siteID)) +
      geom_vline(data = subset(site_species_var, scientificName_Species==species &
                                 n>=cutoff)%>%
                   filter(siteID %in% siteList) %>%
                   filter(domainID %in% domList),
                 aes(xintercept = mean_cm,
                     colour = siteID)) +
      scale_fill_manual(values=c("#006699", "#990000")) +
      scale_color_manual(values=c("#006699", "#990000")) +
      ylab("Site Density") +
      xlab(NULL) +
      theme_pubr() +
      theme(legend.position="none") +
      facet_wrap(.~siteID),
    ggplot(data = subset(all_elytra, scientificName_Species==species) %>%
             filter(plotID %in% plotList) %>%
             filter(domainID %in% domList),
           aes(x=cm_elytra_max_length, fill= plotID)) +
      geom_density(alpha=0.2) +
      geom_rug(aes(colour = plotID)) +
      geom_vline(data = subset(plot_species_var, scientificName_Species==species &
                                 n>=cutoff) %>%
                   filter(plotID %in% plotList) %>%
                   filter(domainID %in% domList),
                 aes(xintercept = mean_cm,
                     color = plotID)) +
      scale_fill_manual(values=c("#00CCCC", "#009999", "#3399cc",
                                 "#3366cc", "#0000ff", "#33ccff",
                                 "#0099cc", "#0066ff", "#3399ff",
                                 "#660000", "#990000", "#cc3333",
                                 "#ff0000", "#cc3300", "#ff6666",
                                 "#ff6600", "#cc0033", "#fc4e2a",
                                 "#d23428")) +
      scale_color_manual(values=c("#00CCCC", "#009999", "#3399cc",
                                  "#3366cc", "#0000ff", "#33ccff",
                                  "#0099cc", "#0066ff", "#3399ff",
                                  "#660000", "#990000", "#cc3333",
                                  "#ff0000", "#cc3300", "#ff6666",
                                  "#ff6600", "#cc0033", "#fc4e2a",
                                  "#d23428")) +
      ylab("Plot Density") +
      xlab("Elytra Length (cm)") +
      theme_pubr()+
      theme(legend.position="none")+
      facet_wrap(.~siteID),
    nrow = 2
  )
)
dev.off()
#################### Questions #########################################3
#Are high CV plots just combined across multiple years?


#png("./Figures/NestedVar.png", height = 10, width = 12, units = "in", res = 300)
ggplot(subset(var_all_scales, scientificName_Species!="Carabidae sp." &  scientificName_Species!="Pterostichus coracinus"),
       aes(scale, var_cm, group = scientificName_Species, col = siteID, shape = domainID)) +
  theme_pubr() + 
  # scale_shape_manual(values=c(16,15,17:25))+
  scale_shape_manual(values=c(16,15,17:25,6:14))+
  geom_point(size = 2) +
  # scale_color_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
  #                               "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
  #                               "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
  #                               "#8A7C64", "#599861")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Variance (cm)") +
  facet_wrap(.~scientificName_Species, scale="free_y") 
dev.off()

ggplot(var_all_scales,
       aes(var_cm)) +
  theme_pubr() + 
  scale_shape_manual(values=c(16,15,17:25))+
  geom_histogram() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("Variance (cm)") +
  facet_wrap(.~scale, scale="free_y", ncol=1) 

ggplot(var_all_scales,
       aes(cv2_pct)) +
  theme_pubr() + 
  geom_histogram() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("Variance as % of mean² (CV² × 100)") +
  facet_wrap(.~scale, scale="free_y", ncol=1)

var_all_scales %>%
  group_by(scale) %>%
  summarise(
    n = n(),
    mean_cvpct = mean(cv2_pct, na.rm = TRUE),
    lower_ci = t.test(cv2_pct)$conf.int[1],
    upper_ci = t.test(cv2_pct)$conf.int[2]
  ) %>%
  ungroup()

#### Geo speatial ####
library(neonDivData)
head(var_all_scales)
neon_location_BET<-subset(neon_location, substr(location_id, (nchar(location_id)-2), nchar(location_id))=="bet")
plot_species_var_pre<-plot_species_var
plot_species_var<-merge(plot_species_var, neon_location_BET, by="plotID")
plot_species_var
site_species_var<-merge(site_species_var, neon_sites, by.x="ID", by.y="siteID")
site_species_var
domain_species_var

ggplot(site_species_var, aes(x=Latitude, y=mean_cm)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  theme_pubr() +
  theme(legend.position = "None")
ggplot(site_species_var, aes(x=Latitude, y=mean_cm, colour = scientificName_Species)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  theme_pubr() +
  theme(legend.position = "None")
ggplot(site_species_var, aes(x=Latitude, y=var_cm, colour = n)) +
  geom_point() +
  theme_pubr()
ggplot(site_species_var, aes(x=Latitude, y=cv2_pct, colour = n)) +
  geom_point() +
  theme_pubr() 
ggplot(site_species_var, aes(x=Latitude, y = n)) +
  geom_point() +
  theme_pubr()
ggplot(site_species_var, aes(x=Latitude, y = range, colour = log10(n))) +
  geom_point() +
  theme_pubr()
ggplot(site_species_var, aes(y = range, x = n)) +
  geom_point() +
  theme_pubr()

ggplot(site_species_var, aes(x=n, y=log(cv2_pct))) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  theme_pubr() +
  theme(legend.position = "None")

####HERE####
#Plot
plot_cummunity_summary<-all_elytra %>%
  filter(!grepl("sp\\.", scientificName_Species)) %>%
  group_by(plotID) %>%
  summarise(
    n_spp = length(unique(scientificName_Species)),
    max_cm = quantile(cm_elytra_max_length, probs=c(.98), na.rm = TRUE),
    min_cm= quantile(cm_elytra_max_length, probs=c(.02), na.rm = TRUE),
    range= (max_cm-min_cm),
    mean_cm = mean(cm_elytra_max_length, na.rm = TRUE),
    mean_logcm = mean(log10(cm_elytra_max_length), na.rm = TRUE),
    skew = skewness(log10(cm_elytra_max_length), na.rm = TRUE),
    kurtosis = datawizard::kurtosis(log10(cm_elytra_max_length), na.rm = TRUE)[["Kurtosis"]],
    kurtosisSE = datawizard::kurtosis(log(cm_elytra_max_length), na.rm = TRUE)[["SE"]]
  ) %>% 
  ungroup()
plot_cummunity_summary
plot_cummunity_summary<-merge(plot_cummunity_summary, neon_location_BET, by="plotID")

#Site
site_cummunity_summary<-all_elytra %>%
  filter(!grepl("sp\\.", scientificName_Species)) %>%
  group_by(siteID) %>%
  summarise(
    n_spp = length(unique(scientificName_Species)),
    max_cm = quantile(cm_elytra_max_length, probs=c(.98), na.rm = TRUE),
    min_cm= quantile(cm_elytra_max_length, probs=c(.02), na.rm = TRUE),
    range= (max_cm-min_cm),
    mean_cm = mean(cm_elytra_max_length, na.rm = TRUE),
    mean_logcm = mean(log10(cm_elytra_max_length), na.rm = TRUE),
    skew = skewness(log10(cm_elytra_max_length), na.rm = TRUE),
    kurtosis = datawizard::kurtosis(log10(cm_elytra_max_length), na.rm = TRUE)[["Kurtosis"]],
    kurtosisSE = datawizard::kurtosis(log10(cm_elytra_max_length), na.rm = TRUE)[["SE"]]
  ) %>% 
  ungroup()
site_cummunity_summary
site_cummunity_summary<-merge(site_cummunity_summary, neon_sites, by="siteID")

ggplot(site_cummunity_summary, aes(x=Latitude, y=mean_cm, colour = n_spp)) +
  geom_point() +
  geom_smooth(method = "lm")+
  theme_pubr()
ggplot(site_cummunity_summary, aes(x=Latitude, y=min_cm, colour = n_spp)) +
  geom_point() +
  geom_smooth(method = "lm")+
  theme_pubr()
ggplot(site_cummunity_summary, aes(x=Latitude, y=max_cm, colour = n_spp)) +
  geom_point() +
  geom_smooth(method = "lm")+
  theme_pubr()
ggplot(site_cummunity_summary, aes(x=Latitude, y = range)) +
  geom_point() +
  theme_pubr()
ggplot(site_cummunity_summary, aes(y = range, x = n_spp)) +
  geom_point() +
  theme_pubr()

ggplot(plot_cummunity_summary, aes(x=latitude, y=mean_cm)) +
  geom_point() +
  geom_smooth(method = "lm")+
  theme_pubr()
ggplot(plot_cummunity_summary, aes(x=latitude, y=min_cm)) +
  geom_point() +
  geom_smooth(method = "lm")+
  theme_pubr()
ggplot(plot_cummunity_summary, aes(x=latitude, y=max_cm)) +
  geom_point() +
  geom_smooth(method = "lm")+
  theme_pubr()

ggplot(plot_cummunity_summary, aes(x=latitude, y = range)) +
  geom_point() +
  geom_smooth(method = "lm")+
  theme_pubr()
ggplot(plot_cummunity_summary, aes(y = range, x = n_spp)) +
  geom_point() +
  theme_pubr()

#Latitude Graphs####
#Mean size
#communit level means
png("./Figures/BodySizeQuantification/CWM(site)_latitude.png", units = "in", width = 6, height = 5, res=300)
ggplot(site_cummunity_summary, aes(x=Latitude, y=mean_cm)) +
  geom_point(col="grey") +
  geom_smooth(method = "lm", col="grey", alpha=0.15)+
  geom_point(data=subset(site_cummunity_summary, Latitude>=25), 
             aes(x=Latitude, y=mean_cm)) +
  geom_smooth(data=subset(site_cummunity_summary, Latitude>=25), 
              aes(x=Latitude, y=mean_cm), 
              method = "lm", col="black")+
  theme_pubr() +
  labs(title = "Community Weighted Mean Bodysize by Latitude",
       subtitle = "Site Level",
       x = "Latitude",
       y = "Elytra Length (cm)")
dev.off()

ggplot(site_cummunity_summary, aes(x=Latitude, y=mean_logcm)) +
  geom_point(col="grey") +
  geom_smooth(method = "lm", col="grey", alpha=0.15)+
  geom_point(data=subset(site_cummunity_summary, Latitude>=25), 
             aes(x=Latitude, y=mean_logcm)) +
  geom_smooth(data=subset(site_cummunity_summary, Latitude>=25), 
              aes(x=Latitude, y=mean_logcm), 
              method = "lm", col="black")+
  theme_pubr() +
  labs(title = "Community Weighted Mean Bodysize by Latitude",
       subtitle = "Site Level",
       x = "Latitude",
       y = "log Elytra Length (cm)")

png("./Figures/BodySizeQuantification/CWM(plot)_latitude.png", units = "in", width = 6, height = 5, res=300)
ggplot(plot_cummunity_summary, aes(x=latitude, y=mean_cm)) +
  geom_point(col="grey") +
  geom_smooth(method = "lm", col="grey", alpha=0.15)+
  geom_point(data=subset(plot_cummunity_summary, latitude>=25), 
             aes(x=latitude, y=mean_cm)) +
  geom_smooth(data=subset(plot_cummunity_summary, latitude>=25), 
              aes(x=latitude, y=mean_cm), 
              method = "lm", col="black")+
  theme_pubr()+
  labs(title = "Community Weighted Mean Bodysize by Latitude",
       subtitle = "Plot Level",
       x = "Latitude",
       y = "Elytra Length (cm)")
dev.off()

ggplot(plot_cummunity_summary, aes(x=latitude, y=mean_logcm)) +
  geom_point(col="grey") +
  geom_smooth(method = "lm", col="grey", alpha=0.15)+
  geom_point(data=subset(plot_cummunity_summary, latitude>=25), 
             aes(x=latitude, y=mean_logcm)) +
  geom_smooth(data=subset(plot_cummunity_summary, latitude>=25), 
              aes(x=latitude, y=mean_logcm), 
              method = "lm", col="black")+
  theme_pubr()+
  labs(title = "Community Weighted Mean Bodysize by Latitude",
       subtitle = "Plot Level",
       x = "Latitude",
       y = "log(Elytra Length (cm))")


png("./Figures/BodySizeQuantification/CWStack(plot)_latitude.png", units = "in", width = 6, height = 5, res=300)
ggplot(plot_cummunity_summary, aes(x=latitude, y=min_cm)) +
  geom_point(col="lightblue") +
  geom_smooth(method = "lm", col="lightblue", alpha=0.15)+
  geom_point(data=subset(plot_cummunity_summary, latitude>=25), 
             aes(x=latitude, y=min_cm), col="navy") +
  geom_smooth(data=subset(plot_cummunity_summary, latitude>=25), 
              aes(x=latitude, y=min_cm), 
              method = "lm", col="navy")+
  geom_point(data = plot_cummunity_summary, aes(x=latitude, y=max_cm), col="#FF9966") +
  geom_smooth(data = plot_cummunity_summary, aes(x=latitude, y=max_cm),
              method = "lm", col="#FF9966", alpha=0.15)+
  geom_point(data=subset(plot_cummunity_summary, latitude>=25), 
             aes(x=latitude, y=max_cm), col="#660000") +
  geom_smooth(data=subset(plot_cummunity_summary, latitude>=25), 
              aes(x=latitude, y=max_cm), 
              method = "lm", col="#660000")+
  geom_point(col="grey", aes(x=latitude, y=mean_cm)) +
  geom_smooth(aes(x=latitude, y=mean_cm) , method = "lm", col="grey", alpha=0.15)+
  geom_point(data=subset(plot_cummunity_summary, latitude>=25), 
             aes(x=latitude, y=mean_cm)) +
  geom_smooth(data=subset(plot_cummunity_summary, latitude>=25), 
              aes(x=latitude, y=mean_cm), 
              method = "lm", col="black") +
  theme_pubr()+
  labs(title = "Community Weighted Bodysize by Latitude",
       subtitle = "Plot Level",
       x = "Latitude",
       y = "Elytra Length (cm)",
       caption = "98% (red), Mean (black), 2% (blue)")
dev.off()

#No aggrigation
ggplot(all_elytra_lat, aes(x=latitude, y=cm_elytra_max_length)) +
  geom_point(col="grey") +
  geom_smooth(method = "lm", col="grey", alpha=0.15)+
  geom_point(data=subset(all_elytra_lat, latitude>=25), 
             aes(x=latitude, y=cm_elytra_max_length)) +
  geom_smooth(data=subset(all_elytra_lat, latitude>=25), 
              aes(x=latitude, y=cm_elytra_max_length), 
              method = "lm", col="black")+
  theme_pubr()

ggplot(all_elytra_lat, aes(x=latitude, y=cm_elytra_max_length)) +
  geom_point(alpha=0.2) +
  geom_smooth(method = "lm")+
  theme_pubr() +
  theme(legend.position = "None")


sort(table(plot_species_var$scientificName_Species), decreasing = T)[1:50]
all_elytra_lat<-merge(all_elytra, neon_location_BET, by="plotID")
head(all_elytra_lat)

all_elytra_lat_subset <- all_elytra_lat %>%
  group_by(scientificName_Species) %>%
  filter(n_distinct(siteID.y) > 3) %>%
  ungroup()

#png("./Figures/BodySizeQuantification/SpeciesBodysize_latitude.png", units = "in", width = 6, height = 5, res=300)
ggplot(all_elytra_lat_subset, aes(x=latitude, y=cm_elytra_max_length, colour = scientificName_Species)) +
  geom_point(alpha=0.2) +
  geom_smooth(method = "lm")+
  theme_pubr() +
  theme(legend.position = "None")+
  labs(title = "Species Bodysize by Latitude",
       x = "Latitude",
       y = "Elytra Length (cm)",
       caption = "subset to species that appear in more than 3 sites")
dev.off()

#Variance in size
#png("./Figures/BodySizeQuantification/SpeciesCVpct(plot)_latitude.png", units = "in", width = 6, height = 5, res=300)
ggplot(subset(plot_species_var, n>20), aes(x=latitude, y=cv2_pct)) +
  geom_point(alpha=0.2) +
  geom_smooth(method = "lm")+
  theme_pubr() +
  theme(legend.position = "None")+
  labs(title = "Species CV% by Latitude",
       subtitle = "Plot Level",
       x = "Latitude",
       y = "Variance as % of mean² (CV² × 100)",
       caption = "Subset to obersvations of species with >20 individuals per plot")
dev.off()

ggplot(subset(plot_species_var, n>20), aes(x=latitude, y=cv2_pct, colour = scientificName_Species)) +
  geom_point() +
  geom_smooth(method = "lm")+
  theme_pubr() +
  theme(legend.position = "None")

#png("./Figures/BodySizeQuantification/CVpct(plot)_N.png", units = "in", width = 6, height = 5, res=300)
ggplot(plot_species_var, aes(x=log10(n), y=log10(cv2_pct))) +
  geom_point(alpha=0.2) +
  theme_pubr() +
  theme(legend.position = "None")
dev.off()

#Range of Body Size
ggplot(site_cummunity_summary, aes(x=latitude, y=range)) +
  geom_point(col="grey") +
  geom_smooth(method = "lm", col="grey", alpha=0.15)+
  geom_point(data=subset(site_cummunity_summary, latitude>=25), 
             aes(x=latitude, y=range)) +
  geom_smooth(data=subset(site_cummunity_summary, latitude>=25), 
              aes(x=latitude, y=range), 
              method = "lm", col="black")+
  theme_pubr()

ggplot(plot_cummunity_summary, aes(x=latitude, y=range)) +
  geom_point(col="grey") +
  geom_smooth(method = "lm", col="grey", alpha=0.15)+
  geom_point(data=subset(plot_cummunity_summary, latitude>=25), 
             aes(x=latitude, y=range)) +
  geom_smooth(data=subset(plot_cummunity_summary, latitude>=25), 
              aes(x=latitude, y=range), 
              method = "lm", col="black")+
  theme_pubr()

ggplot(plot_cummunity_summary, aes(x=mean_cm, y=range)) +
  geom_point(col="grey") +
  geom_smooth(method = "lm", col="grey", alpha=0.15)+
  # geom_point(data=subset(site_cummunity_summary, Latitude>=25), 
  #            aes(x=Latitude, y=mean_cm)) +
  # geom_smooth(data=subset(site_cummunity_summary, Latitude>=25), 
              # aes(x=Latitude, y=mean_cm), 
              # method = "lm", col="black")+
  theme_pubr()

#Skew
ggplot(site_cummunity_summary, aes(x=Latitude, y=skew)) +
  geom_point(col="grey") +
  geom_smooth(method = "lm", col="grey", alpha=0.15)+
  geom_point(data=subset(site_cummunity_summary, Latitude>=25), 
             aes(x=Latitude, y=skew)) +
  geom_smooth(data=subset(site_cummunity_summary, Latitude>=25), 
              aes(x=Latitude, y=skew), 
              method = "lm", col="black")+
  theme_pubr() +
  labs(title = "Community Weighted Mean Bodysize by Latitude",
       subtitle = "Site Level",
       x = "Latitude",
       y = "Skew")

ggplot(site_cummunity_summary, aes(x=Latitude, y=kurtosis)) +
  geom_point(col="grey") +
  geom_smooth(method = "lm", col="grey", alpha=0.15)+
  geom_point(data=subset(site_cummunity_summary, Latitude>=25), 
             aes(x=Latitude, y=kurtosis)) +
  geom_smooth(data=subset(site_cummunity_summary, Latitude>=25), 
              aes(x=Latitude, y=kurtosis), 
              method = "lm", col="black")+
  theme_pubr() +
  labs(title = "Community Weighted Mean Bodysize by Latitude",
       subtitle = "Site Level",
       x = "Latitude",
       y = "kurtosis")

ggplot(plot_cummunity_summary, aes(x=latitude, y=skew)) +
  geom_point(col="grey") +
  geom_smooth(method = "lm", col="grey", alpha=0.15)+
  geom_point(data=subset(plot_cummunity_summary, latitude>=25), 
             aes(x=latitude, y=skew)) +
  geom_smooth(data=subset(plot_cummunity_summary, latitude>=25), 
              aes(x=latitude, y=skew), 
              method = "lm", col="black")+
  theme_pubr() +
  labs(title = "Community Weighted Mean Bodysize by Latitude",
       subtitle = "Plot Level",
       x = "Latitude",
       y = "Skew")

ggplot(plot_cummunity_summary, aes(x=latitude, y=kurtosis)) +
  geom_point(col="grey") +
  geom_smooth(method = "lm", col="grey", alpha=0.15)+
  geom_point(data=subset(plot_cummunity_summary, latitude>=25), 
             aes(x=latitude, y=kurtosis)) +
  geom_smooth(data=subset(plot_cummunity_summary, latitude>=25), 
              aes(x=latitude, y=kurtosis), 
              method = "lm", col="black")+
  theme_pubr() +
  labs(title = "Community Weighted Mean Bodysize by Latitude",
       subtitle = "Plot Level",
       x = "Latitude",
       y = "kurtosis")


ggplot(plot_cummunity_summary, aes(x=(skew)^2, y=kurtosis)) +
  geom_point(col="grey") +
  geom_smooth(method = "lm", col="grey", alpha=0.15)+
  geom_point(data=subset(plot_cummunity_summary, latitude>=25), 
             aes(x=(skew)^2, y=kurtosis)) +
  geom_smooth(data=subset(plot_cummunity_summary, latitude>=25), 
              aes(x=(skew)^2, y=kurtosis), 
              method = "lm", col="black")+
  theme_pubr() +
  labs(title = "Community Weighted Mean Bodysize by Latitude",
       subtitle = "Plot Level",
       x = "(skew)^2",
       y = "kurtosis")

