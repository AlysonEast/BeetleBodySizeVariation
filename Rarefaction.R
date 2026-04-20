library(neonDivData)
library(dplyr)
library(tidyr)
library(vegan)
library(iNEXT)
library(scales)

data(data_beetle)

df <- data_beetle

str(df)

df$taxon_name
table(df$taxon_rank)
df$Species<-ifelse(df$taxon_rank=="species", paste0(df$taxon_name),
                   ifelse(df$taxon_rank=="subspecies", sub("^(\\S+\\s+\\S+).*", "\\1", 
                                                           df$taxon_name), NA))

#### Sites ####
df_counts <- df %>%
  filter(!is.na(Species),
         variable_name == "abundance") %>%
  mutate(trappingDays = as.numeric(trappingDays),
         count_est = round(value * trappingDays))

site_list <- df_counts %>%
  group_by(siteID, Species) %>%
  summarise(abundance = sum(count_est), .groups = "drop") %>%
  group_split(siteID) %>%
  setNames(unique(df_counts$siteID)) %>%
  lapply(function(x) x$abundance[x$abundance > 0])

iNEXT_sites <- iNEXT(site_list,
                     q = 0,
                     datatype = "abundance")

ggiNEXT(iNEXT_sites, type = 1) +
  aes(x=(x/100))+
  theme_pubr() +
  xlab("Number of Individuals (x100)") +
  scale_shape_manual(values=c(rep(4,47))) +
  scale_color_manual(values=c(rep("black",47))) +
  scale_fill_manual(values=c(rep("black",47))) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1),
                     breaks = breaks_extended(n = 4)) +
  theme(legend.position = "none") +
  facet_wrap(.~Assemblage, nrow=4, scale="free")

ggiNEXT(iNEXT_sites, type = 1) +
  theme_pubr() +
  scale_shape_manual(values=c(rep(4,47))) +
  guides(colour = guide_colorbar()) +
  scale_x_continuous(limits=c(0,30000),expand = c(0,0)) +#,
  theme(legend.position = "right") 
  

site_completeness <- iNEXT_sites$AsyEst
head(site_completeness)
site_completeness<-subset(site_completeness, Diversity=="Species richness")
site_completeness$completeness<-site_completeness$Observed/site_completeness$Estimator
site_completeness <- site_completeness %>%
  arrange(desc(completeness))

plot(site_completeness$completeness~site_completeness$Observed)
plot(site_completeness$completeness~site_completeness$Estimator)
plot(site_completeness$completeness~site_completeness$s.e.)
plot(site_completeness$Observed~site_completeness$Estimator)


# Year restricted
df_counts_1819 <- df %>%
  filter(!is.na(Species),
         variable_name == "abundance",
         observation_datetime>="2018-01-01",
         observation_datetime<="2019-12-31") %>%
  mutate(trappingDays = as.numeric(trappingDays),
         count_est = round(value * trappingDays))

site_list_1819 <- df_counts_1819 %>%
  group_by(siteID, Species) %>%
  summarise(abundance = sum(count_est), .groups = "drop") %>%
  group_split(siteID) %>%
  setNames(unique(df_counts_1819$siteID)) %>%
  lapply(function(x) x$abundance[x$abundance > 0])

iNEXT_sites_1819 <- iNEXT(site_list_1819,
                     q = 0,
                     datatype = "abundance")

ggiNEXT(iNEXT_sites_1819, type = 1) +
  theme_pubr() +
  xlab("Number of Individuals (x100)") +
  scale_x_continuous(labels = scales::label_number(accuracy = 1)) +
  theme(legend.position = "none") +
  scale_shape_manual(values=c(rep(4,47))) +
  facet_wrap(.~Assemblage, nrow=4, scale="free")

ggiNEXT(iNEXT_sites, type = 1) +
  theme_pubr() +
  xlab("Number of Individuals (x100)") +
  scale_x_continuous(labels = scales::label_number(accuracy = 1)) +
  theme(legend.position = "none") +
  scale_shape_manual(values=c(rep(4,47))) +
  facet_wrap(.~Assemblage, nrow=4, scale="free")


site_completeness_1819 <- iNEXT_sites_1819$AsyEst
head(site_completeness_1819)
site_completeness_1819<-subset(site_completeness_1819, Diversity=="Species richness")
site_completeness_1819$completeness<-site_completeness_1819$Observed/site_completeness_1819$Estimator
site_completeness_1819 <- site_completeness_1819 %>%
  arrange(desc(Assemblage))
site_completeness <- site_completeness %>%
  arrange(desc(Assemblage))

site_completeness_diff<-site_completeness_1819[,c(3:8)]-site_completeness[,c(3:8)]
site_completeness_diff$Assemblage<-site_completeness$Assemblage
plot(site_completeness_1819$Observed~site_completeness_1819$Estimator)
plot(site_completeness_1819$Estimator~site_completeness$Estimator)
abline(a=0,b=1)
plot(site_completeness_1819$completeness~site_completeness$completeness)


#### Plots ####
valid_plots<- df %>%
  filter(!is.na(Species),
         variable_name == "abundance",
         observation_datetime>="2018-01-01",
         observation_datetime<="2019-12-31")
valid_plots<-unique(valid_plots$plotID)

df_counts <- df %>%
  filter(!is.na(Species),
         variable_name == "abundance",
         plotID %in% valid_plots) %>%
  mutate(trappingDays = as.numeric(trappingDays),
         count_est = round(value * trappingDays))

plot_list <- df_counts %>%
  group_by(plotID, Species) %>%
  summarise(abundance = sum(count_est), .groups = "drop") %>%
  group_split(plotID) %>%
  setNames(unique(df_counts$plotID)) %>%
  lapply(function(x) x$abundance[x$abundance > 0])

iNEXT_plots <- iNEXT(plot_list,
                     q = 0,
                     datatype = "abundance")

ggiNEXT(iNEXT_plots, type = 1) +
  theme_pubr() +
  theme(legend.position = "none")#+
  # scale_shape_manual(values=c(rep(4,47))) +
#  facet_wrap(.~Assemblage, nrow=4, scale="free")

# test<-ggiNEXT(iNEXT_plots, type = 1, color.var = "None")
# test$data$col<-substr(test$data$col, 6, 8)
# test


ggiNEXT(iNEXT_plots, type = 1) + 
  aes(colour = substr(Assemblage, 6, 8), x=(x/100))+
  theme_pubr() +
  xlab("Number of Individuals (x100)") +
  scale_shape_manual(values=c(rep(4,length(valid_plots)))) +
  scale_fill_manual(values = c(rep("grey",523)))+
  scale_x_continuous(labels = scales::label_number(accuracy = 1),
                     breaks = breaks_extended(n = 4)) +
  theme(legend.position = "none") +
  facet_wrap(.~substr(Assemblage, 1, 4), nrow=4, scale="free")

plot_completeness <- iNEXT_plots$AsyEst
head(plot_completeness)
plot_completeness<-subset(plot_completeness, Diversity=="Species richness")
plot_completeness$completeness<-plot_completeness$Observed/plot_completeness$Estimator
plot_completeness <- plot_completeness %>%
  arrange(desc(completeness))
hist(plot_completeness$completeness)

#### domains ####
df_counts <- df %>%
  filter(!is.na(Species),
         variable_name == "abundance") %>%
  mutate(trappingDays = as.numeric(trappingDays),
         count_est = round(value * trappingDays))
data(neon_sites, package = "neonDivData")

df_counts <- df_counts %>%
  left_join(neon_sites %>% select(siteID, domainID),
            by = "siteID")

domain_list <- df_counts %>%
  group_by(domainID, Species) %>%
  summarise(abundance = sum(count_est), .groups = "drop") %>%
  group_split(domainID) %>%
  setNames(unique(df_counts$domainID)) %>%
  lapply(function(x) x$abundance[x$abundance > 0])

iNEXT_domains <- iNEXT(domain_list,
                     q = 0,
                     datatype = "abundance")

ggiNEXT(iNEXT_domains, type = 1) +
  theme_pubr() +
  theme(legend.position = "none")#+
# scale_shape_manual(values=c(rep(4,47))) +
#  facet_wrap(.~Assemblage, nrow=4, scale="free")

ggiNEXT(iNEXT_domains, type = 1) + 
  aes(x=(x/100))+
  theme_pubr() +
  theme(legend.position = "none") +
  xlab("Number of Individuals (x100)") +
  scale_shape_manual(values=c(rep(4,20))) +
  scale_color_manual(values=c(rep("black",20))) +
  scale_fill_manual(values=c(rep("black",20))) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1),
                     breaks = breaks_extended(n = 4)) +
  facet_wrap(.~substr(Assemblage, 1, 4), nrow=4, scale="free")

table(subset(neon_sites, `Site Type`=="Relocatable Terrestrial"|`Site Type`=="Core Terrestrial")$domainID)
oneSiteDomain<-c("D12","D15","D20")
twoSiteDomain<-c("D01","D04","D11","D13","D14","D16","D18")
threeSiteDomain<-c("D01","D03","D05","D06","D07","D08","D09","D10",
                   "D17","D19")
domainSampling<-table(subset(neon_sites, `Site Type`=="Relocatable Terrestrial"|`Site Type`=="Core Terrestrial")$domainID)
domainSampling

ggiNEXT(iNEXT_domains, type = 1) + 
  theme_pubr() +
  theme(legend.position = "right") +
  scale_x_continuous(limits = c(0,50000), expand = c(0,0))+
  scale_shape_manual(values=c(17,15,15,17,15,15,15,15,15,15,
                              17,1,17,17,1,17,15,17,15,1))


domain_completeness <- iNEXT_domains$AsyEst
head(domain_completeness)
domain_completeness<-subset(domain_completeness, Diversity=="Species richness")
domain_completeness$completeness<-domain_completeness$Observed/domain_completeness$Estimator
domain_completeness <- domain_completeness %>%
  arrange(desc(completeness))
hist(domain_completeness$completeness)