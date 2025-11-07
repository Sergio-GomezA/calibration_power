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
gen_path <- "~/Documents/elexon/"

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
  file.path(gen_path, "gen_adj.parquet")
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
t0 <- "2024-01-01"
t1 <- "2024-12-31"

pwr_curv_df <- gen_adj %>%
  filter(between(halfHourEndTime, t0, t1)) %>%
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

class_curves <- fread("data/generic_powerCurves.csv.gz") %>%
  group_by(class) %>%
  mutate(power_scaled = power_kw / ratedPower) %>%
  ungroup()
classes <- class_curves$class %>% unique()
pwr_curv_df$power_est0 <- 0
for (class in classes) {
  # print(class)
  # browser()
  class_ind <- pwr_curv_df$turb_class == class
  pwr_curv_df$power_est0[class_ind] <- approx(
    x = class_curve$wind_speed,
    y = class_curve$power_scaled,
    xout = pwr_curv_df$ws100[class_ind],
    rule = 2
  )$y *
    pwr_curv_df$capacity[class_ind]
}


pwr_curv_df %>%
  filter(grepl("T_FALGW-1", bmUnit)) %>%
  ggplot(aes(power_est0, potential)) +
  geom_point(col = "darkblue", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, col = "red", linetype = 2)
ggsave(
  "fig/falgw_conv.pdf",
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)

pwr_curv_df %>%
  filter(grepl("HOWAO", bmUnit)) %>%
  ggplot(aes(power_est0, potential)) +
  geom_point(col = "darkblue", alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, col = "red", linetype = 2)
ggsave(
  "fig/howao_conv.pdf",
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)

arrow::write_parquet(
  pwr_curv_df,
  file.path(path, "power_curve.parquet")
)


gen_adj %>%
  filter(
    grepl("HOWAO-1", bmUnit),
    between(halfHourEndTime, "2024-01-23", "2024-01-31")
  ) %>%
  ggplot() +
  geom_line(aes(halfHourEndTime, quantity, col = "gen")) +
  geom_line(aes(halfHourEndTime, curtailment, col = "curt")) +
  geom_line(aes(halfHourEndTime, potential, col = "potential")) +
  scale_color_aaas()

gen_adj %>%
  filter(
    grepl("HOWAO-1", bmUnit),
    between(halfHourEndTime, "2024-02-1", "2024-5-1")
  ) %>%
  ggplot() +
  geom_line(aes(halfHourEndTime, quantity, col = "gen")) +
  geom_line(aes(halfHourEndTime, curtailment, col = "curt")) +
  geom_line(aes(halfHourEndTime, potential, col = "potential")) +
  scale_color_aaas()

gen_adj %>%
  filter(
    grepl("HOWAO-1", bmUnit),
    between(halfHourEndTime, "2024-05-1", "2024-5-10")
  ) %>%
  View()
ggplot() +
  geom_line(aes(halfHourEndTime, quantity, col = "gen")) +
  geom_line(aes(halfHourEndTime, curtailment, col = "curt")) +
  geom_line(aes(halfHourEndTime, potential, col = "potential")) +
  scale_color_aaas()

class(gen_adj$halfHourEndTime)
lubridate::tz(gen_adj$halfHourEndTime)
gen_adj %>%
  filter(
    grepl("HOWAO-1", bmUnit),
    between(halfHourEndTime, "2024-05-09T06:00:00Z", "2024-05-09T11:00:00Z")
  ) %>%
  ggplot() +
  geom_line(aes(halfHourEndTime, quantity, col = "gen")) +
  geom_line(aes(halfHourEndTime, curtailment, col = "curt")) +
  geom_line(aes(halfHourEndTime, potential, col = "potential")) +
  scale_color_aaas()

gen_adj %>%
  filter(
    grepl("HOWAO-1", bmUnit),
    between(halfHourEndTime, "2024-5-10", "2024-06-04")
  ) %>% #View()
  ggplot() +
  geom_line(aes(halfHourEndTime, quantity, col = "gen")) +
  geom_line(aes(halfHourEndTime, curtailment, col = "curt")) +
  geom_line(aes(halfHourEndTime, potential, col = "potential")) +
  scale_color_aaas()


gen_adj %>%
  filter(
    grepl("SGRWO-1", bmUnit),
    between(halfHourEndTime, "2025-01-28", "2025-01-29")
  ) %>%
  ggplot() +
  geom_line(aes(halfHourEndTime, quantity, col = "gen")) +
  geom_line(aes(halfHourEndTime, curtailment, col = "curt")) +
  geom_line(aes(halfHourEndTime, potential, col = "potential")) +
  scale_color_aaas()
