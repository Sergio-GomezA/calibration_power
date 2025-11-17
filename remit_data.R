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

## joining years####
# catalog
ref_catalog_2025 <- fread(
  file.path("data/ref_catalog_wind_2025.csv.gz")
)
# Data exploration
remit_df <- lapply(
  2019:2025,
  \(y) {
    read_parquet(
      file.path("~/Documents/elexon/", sprintf("remit_all_%d.parquet", y))
    )
  }
) %>%
  bind_rows()

# get BmUnit
remit_df <- remit_df %>%
  mutate(inelexon = assetId %in% wind.bmus.alt$elexonBmUnit) %>%
  mutate(
    elexonBmUnit = case_when(
      # asset matches bmUnit
      assetId %in% wind.bmus.alt$elexonBmUnit ~ assetId,
      # asset matches with T_
      paste0("T_", assetId) %in% wind.bmus.alt$elexonBmUnit ~ paste0(
        "T_",
        assetId
      ),
      # affected unit matches
      toupper(affectedUnit) %in% wind.bmus.alt$elexonBmUnit ~ affectedUnit,
      # affected unit matches with T_
      paste0("T_", toupper(affectedUnit)) %in%
        wind.bmus.alt$elexonBmUnit ~ paste0("T_", toupper(affectedUnit)),
      # special cases
      grepl("LARYO", assetId) ~ paste0("T_", gsub("O", "W", assetId)),
      grepl("RAMP", assetId) ~ paste0("T_", gsub("RAMP", "RMPNO", assetId)),
      # else
      TRUE ~ NA
    )
  ) %>%
  left_join(
    ref_catalog_2025 %>% select(bmUnit, tech_typ),
    by = c("elexonBmUnit" = "bmUnit")
  )

# check missing
# remit_df %>%
#   filter(
#     eventStatus == "Active",
#     !inelexon,
#     # !is_asset_in_elexon,
#     grepl("Wind", fuelType)
#   ) %>%
#   pull(assetId) %>%
#   unique()

# remit_df %>%
#   filter(
#     eventStatus == "Active",
#     is.na(elexonBmUnit),
#     grepl("Wind", fuelType)
#   ) %>%
#   # pull(assetId) %>%
#   pull(affectedUnit) %>%
#   unique()

# c("RAMP", "T_RMPNO", "LARYO", "LARYW", "FALGW-1", "T_FALGW-1", )
# Figures ####
# event statust type
remit_df %>%
  ggplot() +
  geom_bar(
    aes(x = eventStatus),
    stat = "count",
    position = "dodge",
    fill = "darkblue"
  )

# active by year
remit_df %>%
  filter(eventStatus == "Active") %>%
  mutate(year = factor(year(eventStartTime))) %>%
  ggplot() +
  geom_bar(
    aes(x = year),
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


## included in elexon power data ####
# Wind related events
remit_df %>%
  filter(
    eventStatus == "Active",
    !is.na(elexonBmUnit),
    # grepl("Wind", fuelType)
  ) %>%
  mutate(year = factor(year(eventStartTime))) %>%
  ggplot() +
  geom_bar(
    aes(x = forcats::fct_infreq(tech_typ)),
    stat = "count",
    position = "dodge",
    fill = "darkblue"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18)
  )
# remit_df %>%
#   filter(
#     eventStatus == "Active",
#     !is.na(elexonBmUnit),
#     is.na(tech_typ)
#   ) %>%
#   pull(elexonBmUnit) %>%
#   unique()

# Wind related by year
remit_df %>%
  filter(
    eventStatus == "Active",
    !is.na(elexonBmUnit),
    # grepl("Wind", fuelType)
  ) %>%
  mutate(year = factor(year(eventStartTime))) %>%
  ggplot() +
  geom_bar(
    aes(x = year, fill = tech_typ),
    stat = "count",
    position = "stack"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = mypalette) +
  labs(fill = "", x = "year", y = "count")
