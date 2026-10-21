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
# 
# #### READ IN CLEAN DATA ####
# all_elytra<-read.csv("./Data/BodysizeCombinedClean.csv")
# meta<-read.csv("./Data/NEON_Field_Site_Metadata_20260130.csv")
# 
# #### Explore the data ####
# # Get unique siteID levels ordered by latitude
# site_order <- unique(all_elytra[order(all_elytra$latitude, decreasing = TRUE), "siteID"])
# 
# # Convert siteID into a factor with the correct order
# all_elytra$siteID <- factor(all_elytra$siteID, levels = site_order)
# all_elytra$log_dist_cm<-log10(all_elytra$cm_elytra_max_length)
# 
# library(ggpubr)
# table(subset(all_elytra, siteID=="STER")$plotID)
# sort(table(subset(all_elytra, siteID=="STER")$scientificName_Species))
# table(subset(all_elytra, siteID=="STER")$scientificName_Species, subset(all_elytra, siteID=="STER")$plotID)
# 
# all_elytra <- read.csv("./Data/BodysizeCombinedClean.csv")
# all_elytra <- subset(all_elytra, yearCollected == YEAR)
# message("Year: ", YEAR, "  |  specimens: ", nrow(all_elytra))
# all_elytra$all<-"all"

library(Ostats)   # community_overlap()
library(dplyr)
library(tidyr)
library(dplyr)

setwd("/home/aly/Beetles/BeetleBodySizeVariation")

#### 0. SETTINGS ####
# Defaults are for interactive runs. When driven by run_ByYear.sh these are
# overridden by command-line args, in this order:
#   Rscript Overlap_CustomNulls_ByYear.R <LEVEL> <POOL> <YEAR>
LEVEL   <- "plot"     # "plot" or "site"
POOL    <- "all"      # "site" / "domain" / "all"
AUGMENT <- TRUE       # TRUE: pad sparse species x focal cells to n = 20 (lognormal)
# FALSE: observed individuals only (no simulated rows)

# translate the LEVEL / POOL choices into column names and the CV scale to use
if (LEVEL == "plot") { FOCAL_COL <- "Assemblage"; CV_SCALE <- "Plot Level" }


#### 1. READ CLEAN DATA ####
all_elytra <- read.csv("./Data/BodysizeCombinedClean.csv")
all_elytra<-subset(all_elytra, yearCollected==2018 | yearCollected==2019)

# Drop specimens with no species label (NA or blank ""). A blank name is not a
# real species: it survives %in% filters and silently becomes a pseudo-species in
# every metric, and it breaks name-indexing in the swap null (x[""] returns NA in
# R). Must happen BEFORE augmentation, or the blank cell gets augmented to n = 20.
blank_sp <- is.na(all_elytra$scientificName_Species) | all_elytra$scientificName_Species == ""
if (any(blank_sp)) message("Dropping ", sum(blank_sp), " specimen(s) with missing species label")
all_elytra <- all_elytra[!blank_sp, ]

all_elytra$log_dist_cm <- log10(all_elytra$cm_elytra_max_length)

# crosswalks from the observed rows, used to re-attach pool ids to augmented rows
plot_to_site   <- unique(all_elytra[, c("plotID", "siteID")])
site_to_domain <- unique(all_elytra[, c("siteID", "domainID")])

# average CV^2 (%) for this scale, from the summary file (not hardcoded)
cv_results    <- read.csv("./Outputs/CVpctSummary.csv")
typical_cvpct <- cv_results$mean_cvpct[cv_results$scale == CV_SCALE]
cv2           <- typical_cvpct / 100
message("Level: ", LEVEL, "  |  Pool: ", POOL, "  |  CV scale: ", CV_SCALE,
        " (", round(typical_cvpct, 4), ")")

all_elytra$Assemblage<-paste0(all_elytra$plotID,"_",all_elytra$yearCollected)

#### 2. (OPTIONALLY) AUGMENT SPARSE SPECIES x FOCAL CELLS TO n >= 20 ####
# Start from observed individuals; if AUGMENT, pad each sparse species x focal
# cell up to n = 20 with lognormal draws (same augmentation as Overlap_*.R). With
# AUGMENT = FALSE only species with >= 2 real individuals enter each overlap.
aug <- all_elytra[, c("scientificName_Species", FOCAL_COL, "cm_elytra_max_length")]

if (AUGMENT) {
  low_n <- all_elytra %>%
    group_by(scientificName_Species, across(all_of(FOCAL_COL))) %>%
    summarise(n_obs = n(),
              mean_dist = mean(cm_elytra_max_length, na.rm = TRUE),
              .groups = "drop") %>%
    filter(n_obs < 20)
  
  set.seed(42)
  sim_low_n <- low_n %>%
    rowwise() %>%
    mutate(n_to_add = 20 - n_obs,
           sdlog    = sqrt(log(1 + cv2)),
           meanlog  = log(mean_dist) - (sdlog^2 / 2),
           sim_vals = list(rlnorm(n = n_to_add, meanlog = meanlog, sdlog = sdlog))) %>%
    unnest(cols = sim_vals) %>%
    dplyr::rename(cm_elytra_max_length = sim_vals) %>%
    select(all_of(c("scientificName_Species", FOCAL_COL)), cm_elytra_max_length) %>%
    ungroup()
  
  aug <- rbind(aug, as.data.frame(sim_low_n))
  message("Augmentation ON: added ", nrow(sim_low_n), " simulated individuals")
} else {
  message("Augmentation OFF: observed individuals only (", nrow(aug), " individuals)")
}

aug$log_dist_cm <- log10(aug$cm_elytra_max_length)

aug$plotID<-substr(aug$FOCAL, 1, 8)
head(aug)
#### 3. ATTACH REGIONAL POOL IDS ####
# augmented rows carry only the focal id, so re-attach the higher-level ids
if (LEVEL == "plot") {
  aug <- merge(aug, plot_to_site,   by = "plotID", all.x = TRUE)   # + siteID
  aug <- merge(aug, site_to_domain, by = "siteID", all.x = TRUE)   # + domainID
} else {
  aug <- merge(aug, site_to_domain, by = "siteID", all.x = TRUE)   # + domainID
}
if (POOL == "all") { 
  aug$all<-"all"
}
aug$FOCAL <- aug[[FOCAL_COL]]     # generic working columns used below
aug$POOL  <- aug[[POOL_COL]]

head(aug)

#### OStats for the site level
library(neonDivData)
elytraDF<-aug[,c("Assemblage","domainID","siteID","plotID","scientificName_Species","log_dist_cm")]
head(elytraDF)

Ostats_aug <- Ostats(traits = as.matrix(elytraDF[,'log_dist_cm', drop = FALSE]),
                     sp = factor(elytraDF$scientificName_Species),
                     plots = factor(elytraDF$Assemblage),
                     output = "mean",
                     run_null_model = FALSE)

as.data.frame(Ostats_aug$overlaps_norm)

Assemblage_list<-unique(elytraDF$Assemblage)

for (i in 1:length(Assemblage_list)) {
  png(paste0("./Figures/Overlap/Plots/Overlap",Assemblage_list[i],".png"), units = "in", width = 5, height = 5, res=300)
  print(Ostats_plot(plots = elytraDF$Assemblage,
                    sp = elytraDF$scientificName_Species,
                    traits = elytraDF$log_dist_cm,
                    overlap_dat = Ostats_aug,
                    use_plots = Assemblage_list[i],
                    name_x = 'log10(Elytra Length (cm))',
                    means = FALSE,
                    legend = TRUE,
                    scale = "free_y"))
  dev.off()
}
