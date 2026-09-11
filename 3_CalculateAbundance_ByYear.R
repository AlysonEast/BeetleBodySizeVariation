#### PER-YEAR EFFORT-SCALED ABUNDANCES  (by-year companion) ####
# -----------------------------------------------------------------------------
# By-year companion to CalculateAbundance.R. Identical effort-scaled NEON
# abundances (data_beetle$value), but summed WITHIN year rather than pooled
# across 2018-2019, so Overlap_CustomNulls_ByYear.R can weight each year's
# community with that year's abundances.
#
# Writes one file per year x level:
#   ./Data/site_abund_2018.csv   ./Data/site_abund_2019.csv
#   ./Data/plot_abund_2018.csv   ./Data/plot_abund_2019.csv
# Same columns as the pooled files (<FOCAL>, scientificName_Species, abund), so
# the by-year overlap script's lookup is unchanged apart from the file name.
# Does NOT touch the pooled site_abund.csv / plot_abund.csv.
# -----------------------------------------------------------------------------

library(neonDivData)
library(dplyr)

setwd("/home/aly/Beetles/BeetleBodySizeVariation")

data(data_beetle)
df <- data_beetle


#### 1. SPECIES COLUMN + YEAR  (identical to CalculateAbundance.R) ####
df$scientificName_Species <- ifelse(
  df$taxon_rank == "species",
  paste0(df$taxon_name),
  ifelse(
    df$taxon_rank == "subspecies",
    sub("^(\\S+\\s+\\S+).*", "\\1", df$taxon_name),
    NA
  )
)
df$yearCollected <- as.numeric(substr(df$observation_datetime, 1, 4))

df <- subset(df, yearCollected == 2018 | yearCollected == 2019)
df <- subset(df, !is.na(scientificName_Species))


#### 2. SUM ABUNDANCE WITHIN YEAR, WRITE ONE FILE PER YEAR x LEVEL ####
YEARS <- c(2018, 2019)

for (yr in YEARS) {
  dfy <- subset(df, yearCollected == yr)

  site_abund <- dfy %>%
    group_by(siteID, scientificName_Species) %>%
    summarise(abund = sum(value), .groups = "drop")

  plot_abund <- dfy %>%
    group_by(plotID, scientificName_Species) %>%
    summarise(abund = sum(value), .groups = "drop")

  write.csv(site_abund, sprintf("./Data/site_abund_%d.csv", yr), row.names = FALSE)
  write.csv(plot_abund, sprintf("./Data/plot_abund_%d.csv", yr), row.names = FALSE)
  message("wrote site_abund_", yr, ".csv (", nrow(site_abund), " rows)  and  ",
          "plot_abund_", yr, ".csv (", nrow(plot_abund), " rows)")
}

for (yr in YEARS) {
  dfy <- subset(df, yearCollected == yr)
  
  site_abund <- dfy %>%
    group_by(siteID) %>%
    summarise(abund = sum(value), .groups = "drop")
  
  plot_abund <- dfy %>%
    group_by(plotID) %>%
    summarise(abund = sum(value), .groups = "drop")
  
  write.csv(site_abund, sprintf("./Data/siteTotal_abund_%d.csv", yr), row.names = FALSE)
  write.csv(plot_abund, sprintf("./Data/plotTotal_abund_%d.csv", yr), row.names = FALSE)
  message("wrote siteTotal_abund_", yr, ".csv (", nrow(site_abund), " rows)  and  ",
          "plotTotal_abund_", yr, ".csv (", nrow(plot_abund), " rows)")
}
