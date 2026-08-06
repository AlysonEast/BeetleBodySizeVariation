library(neonDivData)
library(dplyr)
library(tidyr)

setwd("/home/aly/Beetles/BeetleBodySizeVariation")

data(data_beetle)

df <- data_beetle

#### Create Species Column ####
df$scientificName_Species <- ifelse(
  df$taxon_rank == "species",
  paste0(df$taxon_name),
  ifelse(
    df$taxon_rank == "subspecies",
    sub("^(\\S+\\s+\\S+).*", "\\1", df$taxon_name),
    NA
  )
)
df$yearCollected<-as.numeric(substr(df$observation_datetime, 1, 4))

df<-subset(df, yearCollected==2018 | yearCollected==2019)
table(df$scientificName_Species)
table(is.na(df$scientificName_Species))

df<-subset(df, !is.na(scientificName_Species))

site_abund <- df %>% 
  group_by(siteID, scientificName_Species) %>%
  summarise(abund = sum(value))

plot_abund <- df %>% 
  group_by(plotID, scientificName_Species) %>%
  summarise(abund = sum(value))

write.csv(site_abund, "./Data/site_abund.csv", row.names = FALSE)
write.csv(plot_abund, "./Data/plot_abund.csv", row.names = FALSE)