# REMIT data compilation ####

require(arrow)
require(dplyr)
require(rnaturalearth)
require(sf)
require(ggplot2)
require(ggsci)
require(ggthemes)
require(FNN)
require(data.table)
require(geosphere)

source("aux_funct.R")

## read data ####
data_path <- "~/Documents/elexon/data_by_year"


# test <- get_remit(year = 2024, path = data_path, end_time = "2024-01-31")

test <- get_remit(year = 2025, path = data_path, end_time = "2025-09-30")

# Data Download ####

years <- 2023:2024

remit_full <- lapply(
  years,
  \(y) get_remit(year = y, path = data_path)
)

# Data exploration
remit_df <- read_parquet(
  file.path("~/Documents/elexon/", "remit_all_2023.parquet")
) %>%
  mutate(inelexon = assetId %in% wind.bmus.alt$elexonBmUnit)

# event statust type
remit_df %>%
  ggplot() +
  geom_bar(
    aes(x = eventStatus),
    stat = "count",
    position = "dodge",
    fill = "darkblue"
  )
# fuel type
remit_df %>%
  ggplot() +
  geom_bar(
    aes(x = forcats::fct_infreq(fuelType)),
    stat = "count",
    position = "dodge",
    fill = "darkblue"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18)
  )
# not in elexon but classified as wind
remit_df %>%

  filter(!inelexon, grepl("Wind", fuelType)) %>%
  ggplot() +
  geom_bar(
    aes(x = forcats::fct_infreq(fuelType)),
    stat = "count",
    position = "dodge",
    fill = "darkblue"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18)
  )

# in elexon data
remit_df %>%
  filter(assetId %in% wind.bmus.alt$elexonBmUnit) %>%
  ggplot() +
  geom_bar(
    aes(x = forcats::fct_infreq(fuelType)),
    stat = "count",
    position = "dodge",
    fill = "darkblue"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18)
  )


# not in elexon data
remit_df %>%
  filter(assetId %in% wind.bmus.alt$elexonBmUnit) %>%
  ggplot() +
  geom_bar(
    aes(x = forcats::fct_infreq(fuelType)),
    stat = "count",
    position = "dodge",
    fill = "darkblue"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18)
  )

## joining years####

## filtering wind related ####

## included in elexon power data vs ####

# Figures ####
