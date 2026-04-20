df<-read.csv("/home/aly/Downloads/1ykzeqat.csv")

head(df)
str(df)

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

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

sp_20<-subset(sp_count, n>=100)

df20<-df %>% filter(scientific_name %in% sp_20$scientific_name)
table(df20$scientific_name)

ggplot(data = df20, aes(x=length)) +
  geom_density() +
  facet_wrap(.~scientific_name, scales = "free")

sp1<-subset(df, scientific_name=="Metrius contractus")
