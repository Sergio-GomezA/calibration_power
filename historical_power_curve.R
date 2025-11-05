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


year_seq <- seq(2019, 2025, 1)

bmu_df <- lapply(year_seq, function(x) {
  file_path <- file.path(
    path,
    "data_by_year",
    paste0("wind_gen_bmu_", x, ".csv.gz")
  )
  fread(file_path)
}) %>%
  bind_rows()


ref_catalog_2025 <- fread(
  file.path("data/ref_catalog_wind_2025.csv.gz")
)

coords_tb <- read.csv("data/era5_loc_mapping.csv")
wind.bmus.alt <- read.csv("data/wind_bmu_alt.csv")
