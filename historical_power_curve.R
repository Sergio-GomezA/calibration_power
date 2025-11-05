require(arrow)
require(dplyr)
# require(rnaturalearth)
require(sf)
require(ggplot2)
require(ggthemes)
# require(FNN)
require(data.table)
require(parallel)

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

generic_pc <- fread("data/generic_powerCurves.csv.gz")
# filter one wind farm ###

## Offshore no curtailment ####
source("aux_funct.R")
bmu_code <- "HOWAO-1"
pc_temp <- power_curve_data1f(bmu_code)
file_string <- "fig/pc_comparison_%s_val%s.pdf"
ggsave(
  sprintf(file_string, bmu_code, "q0"),
  pc_temp$p_quant,
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)
ggsave(
  sprintf(file_string, bmu_code, "pot"),
  pc_temp$p_pot,
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)
## Offshore with curtailment ####
bmu_code <- "T_SGRWO-1"
pc_temp <- power_curve_data1f(bmu_code)
ggsave(
  sprintf(file_string, bmu_code, "q0"),
  pc_temp$p_quant,
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)
ggsave(
  sprintf(file_string, bmu_code, "pot"),
  pc_temp$p_pot,
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)
## Onshore high curtailment ####
bmu_code <- "T_VKNGW-1"
pc_temp <- power_curve_data1f(bmu_code)
ggsave(
  sprintf(file_string, bmu_code, "q0"),
  pc_temp$p_quant,
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)
ggsave(
  sprintf(file_string, bmu_code, "pot"),
  pc_temp$p_pot,
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)
## Onshore low curtailment ####
bmu_code <- "T_FALGW-1"
pc_temp <- power_curve_data1f(bmu_code)
ggsave(
  sprintf(file_string, bmu_code, "q0"),
  pc_temp$p_quant,
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)
ggsave(
  sprintf(file_string, bmu_code, "pot"),
  pc_temp$p_pot,
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)

# apply to generic power curve to entire data ###

pwr_curv_df <- gen_adj %>%
  mutate(
    quantity = quantity + lag(quantity),
    curtailment = curtailment + lag(curtailment),
    potential = potential + lag(potential)
  ) %>%
  filter(minute(halfHourEndTime) == 0) %>%
  left_join(
    ref_catalog_2025 %>%
      select(bmUnit, matches("lon|lat"), site_name, tech_typ, turb_class),
    by = c("bmUnit")
  ) %>%
  left_join(
    era_df %>% select(time, longitude, latitude, ws100, wd100),
    by = c(
      "halfHourEndTime" = "time",
      "era5lon" = "longitude",
      "era5lat" = "latitude"
    )
  )
pwr_curv_df$power_est0 <- mcmapply(
  generic_pow_conv,
  wind_speed = pwr_curv_df$ws100,
  turb_class = pwr_curv_df$turb_class,
  turb_capacity = pwr_curv_df$capacity,
  mc.cores = parallel::detectCores() - 2, # use all but one core
  SIMPLIFY = TRUE
)


arrow::write_parquet(
  pwr_curv_df,
  file.path(path, "power_curve.parquet")
)
