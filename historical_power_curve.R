require(arrow)
require(dplyr)
# require(rnaturalearth)
require(sf)
require(ggplot2)
require(ggthemes)
# require(FNN)
require(data.table)

# read Data ####
data_path <- "~/Documents/ERA5_at_wf/"

era_df <- read_parquet(
  file.path(data_path, "era5_combined.parquet")
) %>%
  mutate(
    ws100 = sqrt(u100^2 + v100^2),
    ws10 = sqrt(u10^2 + v10^2),
    wd100 = (atan2(-u100, -v100) * 180 / pi) %% 360,
    wd10 = (atan2(-u10, -v10) * 180 / pi) %% 360
  )


gen_adj <- read_parquet(
  file.path(data_path, "gen_adj.parquet")
)

ref_catalog_2025 <- fread(
  file.path("data/ref_catalog_wind_2025.csv.gz")
)

coords_tb <- read.csv("data/era5_loc_mapping.csv")
wind.bmus.alt <- read.csv("data/wind_bmu_alt.csv")

pot_summary <- read.csv("data/hist_pot2024.csv")


# filter one wind farm ###

## Offshore no curtailment ####

bmu_code <- "HOWAO-1"

pwr_curv_1wf <- gen_adj %>%
  filter(grepl(bmu_code, bmUnit)) %>%
  left_join()

## Offshore with curtailment ####
"T_SGRWO-1"

## Onshore high curtailment ####
"T_VKNGW-1"

## Onshore low curtailment ####
"T_FALGW-1"
