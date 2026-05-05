library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(parameters)
library(diptest)
library(moments)

df<-read.csv("/home/aly/Downloads/1ykzeqat.csv")

df <- df %>%
  # Split the single column into 8 numeric columns
  separate(pred_coords_px,
           into = paste0("v", 1:8),
           sep = "[, ]+",   # handles comma or space separated
           convert = TRUE) %>%
  # Compute distance using first four values
  mutate(
    width = sqrt((v3 - v1)^2 + (v4 - v2)^2)
  ) %>%
  # Compute distance using second four values
  mutate(
    length = sqrt((v7 - v5)^2 + (v8 - v6)^2)
  )

sp_count <- df %>%
  group_by(scientific_name) %>%
  summarise(n = n()) %>%
  ungroup()

#Add in Beetlepalooza annotations
getwd()
setwd("/home/aly/Beetles/BeetleBodySizeVariation")
paloozadf<-read.csv("./BeetleMeasurements.csv")
meta<-read.csv("./NEON_Field_Site_Metadata_20260130.csv")

paloozadf<-merge(paloozadf, meta, by.x="siteID", by.y="site_id", all.x=TRUE)
df_sub<-subset(paloozadf, user_name=="IsaFluck")
ElytraLength<-subset(df_sub, structure=="ElytraLength")

#Distribution of Mean Body Size in Carabids####
ElytraSummary<- ElytraLength %>%
  group_by(scientificName) %>%
  summarise(
    n_obs = n(),
    mean_dist = mean(dist_cm, na.rm = TRUE),
    var_dist = var(dist_cm, na.rm = TRUE),
    skew = skewness(dist_cm, na.rm = TRUE),
    kurtosis = kurtosis(dist_cm, na.rm = TRUE),
    dstat = dip.test(dist_cm)[1],
    dpval = dip.test(dist_cm)[2],
    .groups = "drop"
  )

ggplot(data = ElytraSummary, aes(x=log10(mean_dist))) +
  geom_histogram() + 
  theme(legend.position="none") 

#ITV###
#What is the shape of distributions?####
ElytraSummary_n<-subset(ElytraSummary, n_obs>=20)
#Skewness#
ggplot(data = ElytraSummary_n, aes(skew)) +
  geom_histogram() + 
  theme(legend.position="none") +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = -1) +
  geom_vline(xintercept = -2) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = 2) +
  annotate("text", x = 0, y = 15, label = "Normal", angle = 90, vjust = -0.5) +
  annotate("text", x = -1, y = 15, label = "left-skewed", angle = 90, vjust = -0.5) +
  annotate("text", x = 1, y = 15, label = "right-skewed", angle = -90, vjust = -0.5) +
  annotate("text", x = 2, y = 15, label = "Exponential right-skewed", angle = -90, vjust = -0.5)+
  annotate("text", x = -2, y = 15, label = "Exponential left-skewed", angle = 90, vjust = -0.5)

print(paste0("left-skewed: n = ", nrow(subset(ElytraSummary_n,skew > 1))))

print(paste0("right-skewed: n = ", nrow(subset(ElytraSummary_n,skew < -1))))

#Kurtosis#
ggplot(data = ElytraSummary_n, aes(kurtosis)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  annotate("text", x = 0, y = 15, label = "Peaked", angle = -90, vjust = -0.5) +
  annotate("text", x = 0, y = 15, label = "Flattened", angle = 90, vjust = -0.5) 

#Modality#
ggplot(data = ElytraSummary_n, aes(as.numeric(dpval))) +
  geom_histogram() +
  geom_vline(xintercept = 0.05) +
  annotate("text", x = 0.05, y = 15, label = "Multimodal", angle = 90, vjust = -0.5)+
  annotate("text", x = 0.05, y = 15, label = "Unimodial", angle = -90, vjust = -0.5)

print(paste0("Unimodal Distributions: n = ", nrow(subset(ElytraSummary_n,dpval>=0.05))))
  
print(paste0("Multimodal Distributions: n = ", nrow(subset(ElytraSummary_n,dpval<=0.05))))

#Variance across Scales ####
head(ElytraLength)
species_stats <- ElytraLength %>%
  group_by(scientificName) %>%
  summarise(
    n = n(),
    mean_dist_cm = mean(dist_cm, na.rm = TRUE),
    var_dist_cm  = var(dist_cm, na.rm = TRUE),
    sd_dist_cm   = sd(dist_cm, na.rm = TRUE)
  ) %>%
  ungroup()

species_site_stats <- ElytraLength %>%
  group_by(siteID, scientificName) %>%
  summarise(
    n = n(),
    mean_dist_cm = mean(dist_cm, na.rm = TRUE),
    var_dist_cm  = var(dist_cm, na.rm = TRUE),
    sd_dist_cm   = sd(dist_cm, na.rm = TRUE)
  ) %>%
  ungroup()

species_domain_stats <- ElytraLength %>%
  group_by(domain_id, scientificName) %>%
  summarise(
    n = n(),
    mean_dist_cm = mean(dist_cm, na.rm = TRUE),
    var_dist_cm  = var(dist_cm, na.rm = TRUE),
    sd_dist_cm   = sd(dist_cm, na.rm = TRUE)
  ) %>%
  ungroup()

ggplot(ElytraLength, aes(scientificName, log10(dist_cm))) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(species_site_stats,
       aes(scientificName, mean_dist_cm, color = siteID)) +
  geom_point() +
  geom_errorbar(
    aes(
      ymin = mean_dist_cm - sd_dist_cm,
      ymax = mean_dist_cm + sd_dist_cm
    ),
    width = 0
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


species_occurrence <- ElytraLength %>%
  group_by(scientificName) %>%
  summarise(
    n_sites = n_distinct(siteID),
    n_plots = n_distinct(plotID),
    n_obs   = n()
  )
species_keep <- species_occurrence %>%
  filter(n_sites > 1, n_plots > 1) %>%
  pull(scientificName)

ElytraLengthDF_multi <- ElytraLength %>%
  filter(scientificName %in% species_keep)

species_var <- ElytraLengthDF_multi %>%
  group_by(scientificName) %>%
  summarise(
    n = n(),
    mean_cm = mean(dist_cm, na.rm = TRUE),
    var_cm  = var(dist_cm, na.rm = TRUE),
    cv2_pct = 100 * var_cm / mean_cm^2,
    ID = "Spp",
    Above = "None",
    Site = "All",
    Domain = "All"
  ) %>%
  ungroup()

plot_species_var <- ElytraLengthDF_multi %>%
  group_by(plotID, scientificName) %>%
  summarise(
    n = n(),
    mean_cm = mean(dist_cm, na.rm = TRUE),
    var_cm  = var(dist_cm, na.rm = TRUE),
    cv2_pct = 100 * var_cm / mean_cm^2,
    Above = unique(siteID),
    Site = unique(siteID),
    Domain = unique(domain_id)
  ) %>% 
  ungroup()
colnames(plot_species_var)<-c("ID",colnames(plot_species_var)[2:length(colnames(plot_species_var))])

site_species_var <- ElytraLengthDF_multi %>%
  group_by(siteID, scientificName) %>%
  summarise(
    n = n(),
    mean_cm = mean(dist_cm, na.rm = TRUE),
    var_cm  = var(dist_cm, na.rm = TRUE),
    cv2_pct = 100 * var_cm / mean_cm^2,
    Above = unique(domain_id),
    Site = unique(siteID),
    Domain = unique(domain_id)
  ) %>%
  ungroup()
colnames(site_species_var)<-c("ID",colnames(site_species_var)[2:length(colnames(site_species_var))])


domain_species_var <- ElytraLengthDF_multi %>%
  group_by(domain_id, scientificName) %>%
  summarise(
    n = n(),
    mean_cm = mean(dist_cm, na.rm = TRUE),
    var_cm  = var(dist_cm, na.rm = TRUE),
    cv2_pct = 100 * var_cm / mean_cm^2,
    Above = "Spp",
    Site = "All",
    Domain = unique(domain_id)
  ) %>%
  ungroup()
colnames(domain_species_var)<-c("ID",colnames(domain_species_var)[2:length(colnames(domain_species_var))])


plot_species_var$scale <- "Plot Level"
species_var$scale <- "Species Level"
site_species_var$scale <- "Site Level"
domain_species_var$scale <- "Domain Level"

var_all_scales <- bind_rows(species_var,
                            plot_species_var,
                            site_species_var,
                            domain_species_var)

var_all_scales<-var_all_scales %>% filter(n >= 10)

var_all_scales$scale<-factor(
  var_all_scales$scale,
  ordered = TRUE,
  levels = c("Species Level", "Domain Level", "Site Level", "Plot Level"))

head(var_all_scales)

png("./Figures/NestedCVpct.png", height = 10, width = 10, units = "in", res = 300)
ggplot(subset(var_all_scales, scientificName!="Carabidae sp." &  scientificName!="Pterostichus coracinus"),
       aes(scale, cv2_pct, group = scientificName, col = Site, shape = Domain)) +
  theme_pubr() + 
  scale_shape_manual(values=c(16,15,17:25))+
  geom_point(size = 2) +
  scale_color_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
                                "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                                "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
                                "#8A7C64", "#599861")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Variance as % of mean² (CV² × 100)") +
  facet_wrap(.~scientificName) 
dev.off()

png("./Figures/NestedVar.png", height = 10, width = 12, units = "in", res = 300)
ggplot(subset(var_all_scales, scientificName!="Carabidae sp." &  scientificName!="Pterostichus coracinus"),
       aes(scale, var_cm, group = scientificName, col = Site, shape = Domain)) +
  theme_pubr() + 
  scale_shape_manual(values=c(16,15,17:25))+
  geom_point(size = 2) +
  scale_color_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
                                "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                                "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
                                "#8A7C64", "#599861")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Variance (cm)") +
  facet_wrap(.~scientificName, scale="free_y") 
dev.off()

ggplot(subset(var_all_scales, scientificName!="Carabidae sp." &  scientificName!="Pterostichus coracinus"),
       aes(var_cm)) +
  theme_pubr() + 
  scale_shape_manual(values=c(16,15,17:25))+
  geom_histogram() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("Variance (cm)") +
  facet_wrap(.~scale, scale="free_y", ncol=1) 

ggplot(subset(var_all_scales, scientificName!="Carabidae sp." &  scientificName!="Pterostichus coracinus"),
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
    mean = mean(cv2_pct, na.rm = TRUE),
    lower_ci = t.test(cv2_pct)$conf.int[1],
    upper_ci = t.test(cv2_pct)$conf.int[2]
  ) %>%
  ungroup()

