library(ggplot2)
library(dplyr)
library(ggplot2)
library(sf)
library(maps)
library(dplyr)
library(ggrepel)

getwd()
setwd("/home/aly/Beetles/BeetleBodySizeVariation")
df<-read.csv("./BeetleMeasurements.csv")
meta<-read.csv("./NEON_Field_Site_Metadata_20260130.csv")
puum<-read.csv("./trait_annotations.csv")

df<-merge(df, meta, by.x="siteID", by.y="site_id", all.x=TRUE)

head(df)

df_sub<-subset(df, user_name=="IsaFluck")

head(df_sub)
str(df_sub)

ElytraLength<-subset(df_sub, structure=="ElytraLength")

hist(ElytraLength$dist_cm)

site_spp<-table(ElytraLength$species, ElytraLength$siteID)
site_spp<-as.data.frame(site_spp)

ElytraLength<-subset(ElytraLength, siteID!="DELA" & siteID!="DSNY")

#### Explore the data ####
# Get unique siteID levels ordered by latitude
site_order <- unique(ElytraLength[order(ElytraLength$latitude, decreasing = TRUE), "siteID"])

# Convert siteID into a factor with the correct order
ElytraLength$siteID <- factor(ElytraLength$siteID, levels = site_order)
ElytraLength$log_dist_cm<-log10(ElytraLength$dist_cm)

ggplot(data = ElytraLength, aes(x=log_dist_cm, fill = scientificName)) +
  geom_density(alpha=0.5) +
  theme(legend.position="none") +
  facet_wrap(.~siteID, scales = "free_y")

#BLAN
ggplot(data = subset(ElytraLength, siteID=="BLAN"), aes(x=dist_cm, fill = scientificName)) +
  geom_density(alpha=0.5) +
  theme(legend.position="none") +
  facet_wrap(.~siteID, scales = "free_y")

ggplot(data = subset(ElytraLength, siteID=="BLAN"), aes(x=dist_cm, fill = scientificName)) +
  geom_density(alpha=0.5) +
  theme(legend.position="none") +
  facet_wrap(.~plotID, scales = "free_y")
table(subset(ElytraLength, siteID=="BLAN")$plotID)

#MLBS
ggplot(data = subset(ElytraLength, siteID=="MLBS"), aes(x=dist_cm, fill = scientificName)) +
  geom_density(alpha=0.5) +
  theme(legend.position="none") +
  facet_wrap(.~siteID, scales = "free_y")

ggplot(data = subset(ElytraLength, siteID=="MLBS"), aes(x=dist_cm, fill = scientificName)) +
  geom_density(alpha=0.5) +
  theme(legend.position="none") +
  facet_wrap(.~plotID, scales = "free_y")
table(subset(ElytraLength, siteID=="MLBS")$plotID)


library(ggpubr)
ggarrange(nrow=2,
  ggplot(data = subset(ElytraLength, siteID=="STER"), aes(x=log_dist_cm, fill = scientificName)) +
    geom_density(alpha=0.5) +
  #  theme(legend.position="none") +
    facet_wrap(.~siteID, scales = "free_y"),
  ggplot(data = subset(ElytraLength, siteID=="STER"), aes(x=log_dist_cm, fill = scientificName)) +
    geom_density(alpha=0.5) +
  #  theme(legend.position="0") +
    facet_wrap(.~plotID, scales = "free_y")
)
table(subset(ElytraLength, siteID=="STER")$plotID)
table(subset(ElytraLength, siteID=="STER")$scientificName)
table(subset(ElytraLength, siteID=="STER")$scientificName, subset(ElytraLength, siteID=="STER")$plotID)

ggarrange(nrow=2,
          ggplot(data = subset(ElytraLength, plotID=="STER_006"), aes(x=dist_cm, fill = scientificName)) +
            geom_density(alpha=0.5) +
            theme(legend.position="none") +
            xlim(0.4,1.6) +
            facet_wrap(.~siteID, scales = "free_y"),
          ggplot(data = subset(ElytraLength, siteID=="STER"), aes(x=dist_cm, fill = scientificName)) +
            geom_density(alpha=0.5) +
            theme(legend.position="none") +
            xlim(0.4,1.6)
)
table(subset(ElytraLength, siteID=="STER")$scientificName)
table(subset(ElytraLength, plotID=="STER_006")$scientificName)

#### Overlap Stats ####
library(Ostats)

ElytraLengthDF<-ElytraLength[,c("domain_id","siteID","plotID","scientificName","dist_cm")]

Ostats_example <- Ostats(traits = as.matrix(ElytraLengthDF[,'dist_cm', drop = FALSE]),
                         sp = factor(ElytraLengthDF$scientificName),
                         plots = factor(ElytraLengthDF$siteID),
                         random_seed = 517)

as.data.frame(Ostats_example$overlaps_norm)
Ostats_example$overlaps_norm_ses

sites<-levels(ElytraLengthDF$siteID)

Ostats_plot(plots = ElytraLengthDF$siteID, 
            sp = ElytraLengthDF$scientificName, 
            traits = ElytraLengthDF$dist_cm, 
            overlap_dat = Ostats_example, 
            #use_plots = sites[14], 
            name_x = 'Elytra Length (cm)', 
            means = FALSE,
            n_col = 6,
            scale = "free_y") 

#### Map of overlaps ####
meta
Ostats_map<-merge(meta, Ostats_example$overlaps_norm, by.x="site_id", by.y = "row.names")
Ostats_map <- Ostats_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Convert state map to dataframe
states_map <- map_data("state")

ggplot() +
  # Base map
  geom_polygon(
    data = states_map,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray70",
    linewidth = 0.2
  ) +
  
  # Site points
  geom_sf(
    data = Ostats_map,
    aes(color = dist_cm),
    size = 6,
    alpha = 0.9,
    inherit.aes = FALSE
  ) +
  
  # Site labels
  geom_text_repel(
    data = Ostats_map,
    aes(label = site_id, geometry = geometry),
    stat = "sf_coordinates",
    size = 3,
    min.segment.length = 0,
    segment.color = NA,
    seed = 42,
    inherit.aes = FALSE
  ) +
  
  # Color scale
  scale_color_viridis_c(
    name = "Overlap",
    option = "magma"
  ) +
  
  coord_sf() +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  
  labs(
    title = "Site Locations",
    subtitle = "Points colored by overlap",
    x = NULL,
    y = NULL
  )

#### Species wise density distributions ####
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

ggplot(ElytraLength, aes(scientificName, dist_cm)) +
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


table(ElytraLengthDF_multi$scientificName)

####Synuchus impunctatus####
SYIM_all<-subset(ElytraLengthDF_multi, scientificName=="Synuchus impunctatus")
SYIM_Domain<-SYIM_all
SYIM_Site<-SYIM_all
SYIM_Plot<-SYIM_all

SYIM_all$Level<-"Species"
SYIM_Domain$Level<-paste0("Domain")
SYIM_Site$Level<-paste0("Site")
SYIM_Plot$Level<-paste0("Plot")

SYIM_all$LevelID<-"Species"
SYIM_Domain$LevelID<-paste0("Domain:",SYIM_Domain$domain_id)
SYIM_Site$LevelID<-paste0("Domain:",SYIM_Site$domain_id," Site",SYIM_Site$siteID)
SYIM_Plot$LevelID<-paste0("Domain:",SYIM_Plot$domain_id," Plot",SYIM_Plot$plotID)
table(SYIM_Plot$LevelID)

SYIM<-rbind(SYIM_Domain,SYIM_Site,SYIM_Plot)

png("./Figures/SynuchusImpunctatusNestedDensity.png", height = 9, width = 7, units = "in", res = 300)
annotate_figure(
  ggarrange(
    ggplot(data = SYIM_Domain,
           aes(x=dist_cm, fill = LevelID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Domain Density") +
      facet_wrap(.~domain_id),
    ggplot(data = SYIM_Site,
           aes(x=dist_cm, fill = LevelID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Site Density") +
      facet_wrap(.~domain_id),
    ggplot(data = SYIM_Plot,
           aes(x=dist_cm, fill = LevelID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Plot Density") +
      facet_wrap(.~domain_id),
    nrow = 3
  ),
  top = text_grob("Synuchus impunctatus")
)
dev.off()

png("./Figures/SynuchusImpunctatusSitePlotNestedDensity.png", height = 5, width = 7, units = "in", res = 300)
annotate_figure(
  ggarrange(
    ggplot(data = SYIM_Site,
           aes(x=dist_cm, fill = LevelID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Site Density") +
      facet_wrap(.~siteID, nrow=1),
    ggplot(data = SYIM_Plot,
           aes(x=dist_cm, fill = LevelID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Plot Density") +
      facet_wrap(.~siteID, nrow=1),
    nrow = 2
  ),
  top = text_grob("Synuchus impunctatus")
)
dev.off()

####Calathus advena####
SYIM_all<-subset(ElytraLengthDF_multi, scientificName=="Calathus advena")
SYIM_Domain<-SYIM_all
SYIM_Site<-SYIM_all
SYIM_Plot<-SYIM_all

SYIM_all$Level<-"Species"
SYIM_Domain$Level<-paste0("Domain")
SYIM_Site$Level<-paste0("Site")
SYIM_Plot$Level<-paste0("Plot")

SYIM_all$LevelID<-"Species"
SYIM_Domain$LevelID<-paste0("Domain:",SYIM_Domain$domain_id)
SYIM_Site$LevelID<-paste0("Domain:",SYIM_Site$domain_id," Site",SYIM_Site$siteID)
SYIM_Plot$LevelID<-paste0("Domain:",SYIM_Plot$domain_id," Plot",SYIM_Plot$plotID)
table(SYIM_Plot$LevelID)

SYIM<-rbind(SYIM_Domain,SYIM_Site,SYIM_Plot)

png("./Figures/CalathusAdvenaNestedDensity.png", height = 9, width = 7, units = "in", res = 300)
annotate_figure(
  ggarrange(
    ggplot(data = SYIM_Domain,
           aes(x=dist_cm, fill = LevelID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Domain Density") +
      facet_wrap(.~domain_id),
    ggplot(data = SYIM_Site,
           aes(x=dist_cm, fill = LevelID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Site Density") +
      facet_wrap(.~domain_id),
    ggplot(data = SYIM_Plot,
           aes(x=dist_cm, fill = LevelID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Plot Density") +
      facet_wrap(.~domain_id),
    nrow = 3
  ),
  top = text_grob("Calathus advena")
)
dev.off()

####Chlaenius aestivus####
SYIM_all<-subset(ElytraLengthDF_multi, scientificName=="Chlaenius aestivus")
SYIM_Domain<-SYIM_all
SYIM_Site<-SYIM_all
SYIM_Plot<-SYIM_all

SYIM_all$Level<-"Species"
SYIM_Domain$Level<-paste0("Domain")
SYIM_Site$Level<-paste0("Site")
SYIM_Plot$Level<-paste0("Plot")

SYIM_all$LevelID<-"Species"
SYIM_Domain$LevelID<-paste0("Domain:",SYIM_Domain$domain_id)
SYIM_Site$LevelID<-paste0("Domain:",SYIM_Site$domain_id," Site",SYIM_Site$siteID)
SYIM_Plot$LevelID<-paste0("Domain:",SYIM_Plot$domain_id," Plot",SYIM_Plot$plotID)
table(SYIM_Plot$LevelID)

SYIM<-rbind(SYIM_Domain,SYIM_Site,SYIM_Plot)

png("./Figures/ChlaeniusAestivusNestedDensity.png", height = 9, width = 7, units = "in", res = 300)
annotate_figure(
  ggarrange(
    ggplot(data = SYIM_Domain,
           aes(x=dist_cm, fill = LevelID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Domain Density") +
      facet_wrap(.~domain_id),
    ggplot(data = SYIM_Site,
           aes(x=dist_cm, fill = siteID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Site Density") +
      facet_wrap(.~domain_id),
    ggplot(data = SYIM_Plot,
           aes(x=dist_cm, fill = LevelID, color = siteID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Plot Density") +
      facet_wrap(.~domain_id),
    nrow = 3
  ),
  top = text_grob("Chlaenius aestivus")
)
dev.off()

annotate_figure(
  ggarrange(
    ggplot(data = SYIM_Domain,
           aes(x=dist_cm)) +
      geom_density(alpha=0.5, aes(fill=siteID)) +
      geom_density(alpha=0.5, fill="black") +
      theme(legend.position="none") +
      ylab("Domain Density") +
      facet_wrap(.~domain_id),
    nrow = 1
  ),
  top = text_grob("Chlaenius aestivus")
)
annotate_figure(
  ggarrange(
    ggplot(data = SYIM_Domain,
           aes(x=dist_cm)) +
      geom_density(alpha=0.5, aes(fill=plotID, col=plotID)) +
      geom_density(alpha=0.5, fill="black") +
      theme(legend.position="none") +
      ylab("Domain Density") +
      facet_wrap(.~siteID),
    nrow = 1
  ),
  top = text_grob("Chlaenius aestivus")
)

annotate_figure(
  ggarrange(
    ggplot(data = SYIM_Domain,
           aes(x=dist_cm)) +
      geom_density(alpha=0.5, fill="black") +
      geom_histogram(alpha=0.5, position = "identity", aes(fill=plotID, col=plotID)) +
      theme(legend.position="none") +
      ylab("Domain Density") +
      facet_wrap(.~siteID),
    nrow = 1
  ),
  top = text_grob("Chlaenius aestivus")
)

png("./Figures/ChlaeniusAestivusSitePlotNestedDensity.png", height = 6, width = 12, units = "in", res = 300)
annotate_figure(
  ggarrange(
    ggplot(data = SYIM_Domain,
           aes(x=dist_cm)) +
      geom_density(alpha=0.5, position = "identity", aes(fill=plotID, col=plotID)) +
      geom_density(alpha=0.5, fill="black") +
      #theme(legend.position="none") +
      ylab("Domain Density") +
      facet_wrap(.~siteID),
    nrow = 1
  ),
  top = text_grob("Chlaenius aestivus")
)
dev.off()


####Carabus goryi####
SYIM_all<-subset(ElytraLengthDF_multi, scientificName=="Carabus goryi")
SYIM_Domain<-SYIM_all
SYIM_Site<-SYIM_all
SYIM_Plot<-SYIM_all

SYIM_all$Level<-"Species"
SYIM_Domain$Level<-paste0("Domain")
SYIM_Site$Level<-paste0("Site")
SYIM_Plot$Level<-paste0("Plot")

SYIM_all$LevelID<-"Species"
SYIM_Domain$LevelID<-paste0("Domain:",SYIM_Domain$domain_id)
SYIM_Site$LevelID<-paste0("Domain:",SYIM_Site$domain_id," Site",SYIM_Site$siteID)
SYIM_Plot$LevelID<-paste0("Domain:",SYIM_Plot$domain_id," Plot",SYIM_Plot$plotID)
table(SYIM_Plot$LevelID)

SYIM<-rbind(SYIM_Domain,SYIM_Site,SYIM_Plot)

png("./Figures/CarabusGoryiNestedDensity.png", height = 9, width = 7, units = "in", res = 300)
annotate_figure(
  ggarrange(
    ggplot(data = SYIM_Domain,
           aes(x=dist_cm, fill = LevelID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Domain Density") +
      facet_wrap(.~domain_id),
    ggplot(data = SYIM_Site,
           aes(x=dist_cm, fill = LevelID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Site Density") +
      facet_wrap(.~domain_id),
    ggplot(data = SYIM_Plot,
           aes(x=dist_cm, fill = LevelID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Plot Density") +
      facet_wrap(.~domain_id),
    nrow = 3
  ),
  top = text_grob("Carabus goryi")
)
dev.off()


####Cyclotrachelus torvus####
SYIM_all<-subset(ElytraLengthDF_multi, scientificName=="Cyclotrachelus torvus")
SYIM_Domain<-SYIM_all
SYIM_Site<-SYIM_all
SYIM_Plot<-SYIM_all

SYIM_all$Level<-"Species"
SYIM_Domain$Level<-paste0("Domain")
SYIM_Site$Level<-paste0("Site")
SYIM_Plot$Level<-paste0("Plot")

SYIM_all$LevelID<-"Species"
SYIM_Domain$LevelID<-paste0("Domain:",SYIM_Domain$domain_id)
SYIM_Site$LevelID<-paste0("Domain:",SYIM_Site$domain_id," Site",SYIM_Site$siteID)
SYIM_Plot$LevelID<-paste0("Domain:",SYIM_Plot$domain_id," Plot",SYIM_Plot$plotID)
table(SYIM_Plot$LevelID)

SYIM<-rbind(SYIM_Domain,SYIM_Site,SYIM_Plot)

png("./Figures/CyclotrachelusTorvusNestedDensity.png", height = 9, width = 7, units = "in", res = 300)
annotate_figure(
  ggarrange(
    ggplot(data = SYIM_Domain,
           aes(x=dist_cm, fill = LevelID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Domain Density") +
      facet_wrap(.~domain_id),
    ggplot(data = SYIM_Site,
           aes(x=dist_cm, fill = LevelID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Site Density") +
      facet_wrap(.~domain_id),
    ggplot(data = SYIM_Plot,
           aes(x=dist_cm, fill = LevelID)) +
      geom_density(alpha=0.5) +
      theme(legend.position="none") +
      ylab("Plot Density") +
      facet_wrap(.~domain_id),
    nrow = 3
  ),
  top = text_grob("Cyclotrachelus torvus")
)
dev.off()
#### Average distribution shape ####

#Compute typical distribution
# Only high-n species-plot combos
high_n <- ElytraLengthDF %>%
  group_by(scientificName, plotID) %>%
  summarise(
    n_obs = n(),
    mean_dist = mean(dist_cm, na.rm = TRUE),
    var_dist = var(dist_cm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_obs >= 20)

# Typical variance across well-sampled species
typical_var <- median(high_n$var_dist, na.rm = TRUE)
typical_sd <- sqrt(typical_var)

#Pull out low observation combinations
low_n <- ElytraLengthDF %>%
  group_by(scientificName, siteID) %>%
  summarise(n_obs = n(), mean_dist = mean(dist_cm, na.rm = TRUE), .groups = "drop") %>%
  filter(n_obs < 10)


#Generate additional simulated observations
set.seed(42)

library(tidyr)
sim_low_n <- low_n %>%
  rowwise() %>%
  mutate(
    n_to_add = 10 - n_obs,  # fill up to 10
    sim_vals = list(rnorm(n_to_add, mean = mean_dist, sd = typical_sd))
  ) %>%
  unnest(cols = c(sim_vals)) %>%
  rename(dist_cm = sim_vals) %>%
  mutate(domain_id = NA)  # will fill later from original DF

sim_low_n$plotID<-NA

head(ElytraLengthDF)
head(sim_low_n)

ElytraLength_aug <- rbind(ElytraLengthDF, sim_low_n[,colnames(ElytraLengthDF)])

Ostats_aug <- Ostats(traits = as.matrix(ElytraLength_aug[,'dist_cm', drop = FALSE]),
                         sp = factor(ElytraLength_aug$scientificName),
                         plots = factor(ElytraLength_aug$siteID),
                         random_seed = 517)

Ostats_plot(plots = ElytraLengthDF$siteID, 
            sp = ElytraLengthDF$scientificName, 
            traits = ElytraLengthDF$dist_cm, 
            overlap_dat = Ostats_example, 
            #use_plots = sites[14], 
            name_x = 'Elytra Length (cm)', 
            means = FALSE,
            n_col = 6,
            scale = "free_y") 

Ostats_plot(plots = ElytraLength_aug$siteID, 
            sp = ElytraLength_aug$scientificName, 
            traits = ElytraLength_aug$dist_cm, 
            overlap_dat = Ostats_aug, 
            #use_plots = sites[14], 
            name_x = 'Elytra Length (cm)', 
            means = FALSE,
            n_col = 6,
            scale = "free_y") 

par(mfrow=c(1,2))
png("./Figures/Overlap.png", height = 3, width = 9, units = "in", res = 300)
Ostats_plot(plots = ElytraLengthDF$siteID, 
            sp = ElytraLengthDF$scientificName, 
            traits = ElytraLengthDF$dist_cm, 
            overlap_dat = Ostats_example, 
            use_plots = c("SERC","SRER"), 
            name_x = 'Elytra Length (cm)', 
            means = TRUE,
            n_col = 6,
            scale = "free",
            limits_x = c(0,1)) 
dev.off()

png("./Figures/Overlap_augmented.png", height = 3, width = 9, units = "in", res = 300)
Ostats_plot(plots = ElytraLength_aug$siteID, 
            sp = ElytraLength_aug$scientificName, 
            traits = ElytraLength_aug$dist_cm, 
            overlap_dat = Ostats_aug, 
            use_plots = c("SERC","SRER"), 
            name_x = 'Elytra Length (cm)', 
            means = TRUE,
            n_col = 6,
            scale = "free",
            limits_x = c(0,1)) 
dev.off()
Ostats_aug_map<-merge(meta, Ostats_aug$overlaps_norm, by.x="site_id", by.y = "row.names")
Ostats_aug_map <- Ostats_aug_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

ggplot() +
  # Base map
  geom_polygon(
    data = states_map,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray70",
    linewidth = 0.2
  ) +
  
  # Site points
  geom_sf(
    data = Ostats_aug_map,
    aes(color = dist_cm),
    size = 6,
    alpha = 0.9,
    inherit.aes = FALSE
  ) +
  
  # Site labels
  geom_text_repel(
    data = Ostats_aug_map,
    aes(label = site_id, geometry = geometry),
    stat = "sf_coordinates",
    size = 3,
    min.segment.length = 0,
    segment.color = NA,
    seed = 42,
    inherit.aes = FALSE
  ) +
  
  # Color scale
  scale_color_viridis_c(
    name = "Overlap",
    option = "magma"
  ) +
  
  coord_sf() +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  
  labs(
    title = "Site Locations",
    subtitle = "Points colored by overlap",
    x = NULL,
    y = NULL
  )

#### Species means
ggplot(data = ElytraLength, aes(x=log10(dist_cm), fill = scientificName)) +
     geom_histogram() + 
     theme(legend.position="none") 

ElytraMean<- ElytraLength %>%
  group_by(scientificName) %>%
  summarise(
    n_obs = n(),
    mean_dist = mean(dist_cm, na.rm = TRUE),
    var_dist = var(dist_cm, na.rm = TRUE),
    .groups = "drop"
  )
ggplot(data = ElytraMean, aes(x=log10(mean_dist))) +
  geom_histogram() + 
  theme(legend.position="none") 

#### PUUM Explore ####
# Load the beetle data from NEON API using the 'neonUtilities' package
Beetle_dpID <- "DP1.10022.001"
# Load NEON API token for accessing the data
NEON_TOKEN <- read.delim("../../NEON_TOKEN", header = FALSE)[1, 1]  # Read in the NEON API token path anonymized

# Load the beetle data from NEON API using the 'neonUtilities' package
NEON_df <- neonUtilities::loadByProduct(
  dpID = Beetle_dpID,  # Data product ID for beetles
  site = "PUUM",  # Specify the site (PUUM is one of the NEON sites)
  token = NEON_TOKEN,  # API token for accessing the NEON data
  include.provisional = TRUE,  # Include provisional data (Must for PUUM)
  check.size = FALSE)  # Do not check the file size

# Extract taxonomist IDs for further analysis
Neon_para <- NEON_df$bet_parataxonomistID  # Parataxonomist IDs from NEON data
Neon_expert <- NEON_df$bet_expertTaxonomistIDProcessed  # Expert taxonomist IDs from NEON data

#filter out parataxonomist IDs that were later re-IDed by experts
neon_para_clean <- Neon_para %>%
  filter(!(individualID %in% Neon_expert$individualID))
# Find common columns between the two dataframes
common_cols <- intersect(names(neon_para_clean), names(Neon_expert))

# Subset dataframes to the common columns
neon_para_common <- neon_para_clean[, common_cols]
neon_expert_common <- Neon_expert[, common_cols]

#Add a column to specify the identifier type
neon_para_common$ID_status<-"para"
neon_expert_common$ID_status<-"expert"

# Row bind the two dataframes to create on harmonized reference
combined_data <- bind_rows(neon_para_common, neon_expert_common)

puum<-subset(puum, cm_elytra_max_length>0)
puum<-merge(puum, combined_data, by ="individualID")
  
ggplot(data = puum, aes(x = cm_elytra_max_length, fill = scientificName)) +
  geom_density(alpha = 0.5) +
  xlab("Elytra Length (cm)") +
  ylab("Density") + 
#  scale_fill_discrete(labels = legend_labels, name = "Species (n)") +
  theme_pubr()

table(puum$scientificName)

puum_n_site <- puum %>%
  group_by(scientificName, plotID) %>%
  summarise(
    n_obs = n(),
    mean_cm = mean(cm_elytra_max_length, na.rm = TRUE),
    var_cm  = var(cm_elytra_max_length, na.rm = TRUE),
    cv2_pct = 100 * var_cm / mean_cm^2,
    .groups = "drop"
  ) %>%
  filter(n_obs >= 19)

puum_n_site$var_mm<-puum_n_site$var_cm*100

puum_n <- puum %>%
  group_by(scientificName) %>%
  summarise(
    n_obs = n(),
    mean_cm = mean(cm_elytra_max_length, na.rm = TRUE),
    var_cm  = var(cm_elytra_max_length, na.rm = TRUE),
    cv2_pct = 100 * var_cm / mean_cm^2,
    .groups = "drop"
  ) %>%
  filter(n_obs >= 19)

puum_n$var_mm<-puum_n$var_cm*100

puum_n

#### Inter annotator Variation ####
table(df$user_name)
df_lenght<-subset(df, structure=="ElytraLength")
table(df_lenght$user_name)
library(reshape2)

