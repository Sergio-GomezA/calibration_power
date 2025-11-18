require(arrow)
require(dplyr)
# require(rnaturalearth)
require(sf)
require(ggplot2)
require(ggthemes)
require(ggsci)
# require(FNN)
require(data.table)
require(parallel)

require(brms)

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

# historical generation with curtailment and outages list
gen_adj <- read_parquet(
  file.path(gen_path, "gen_adj_v2.parquet")
)
# wind farm catalog based on 2025 data
ref_catalog_2025 <- fread(
  file.path("data/ref_catalog_wind_2025.csv.gz")
)

# era 5 coordinates list
coords_tb <- read.csv("data/era5_loc_mapping.csv")
# elexon wind bmus with updated capacity
wind.bmus.alt <- read.csv("data/wind_bmu_alt.csv")

# potential energy historical summary
pot_summary <- read.csv("data/hist_pot2024.csv")

generic_pc <- fread("data/generic_powerCurves.csv.gz")
# filter one wind farm ###

## Offshore no curtailment ####
source("aux_funct.R")
bmu_code <- "HOWAO-1"
pc_temp <- power_curve_data1f(bmu_code)
file_string <- "fig/pc_comparison_%s_val%s.pdf"
fig_codes <- c("q0", "pot", "out", "nout")

lapply(
  2:length(pc_temp),
  \(i) {
    ggsave(
      filename = sprintf(file_string, bmu_code, fig_codes[i - 1]),
      plot = pc_temp[[i]],
      width = 6,
      height = 4,
      units = "in",
      dpi = 300
    )
  }
)

bmu_code <- "T_HOWAO-2"
pc_temp <- power_curve_data1f(bmu_code)

lapply(
  2:length(pc_temp),
  \(i) {
    ggsave(
      filename = sprintf(file_string, bmu_code, fig_codes[i - 1]),
      plot = pc_temp[[i]],
      width = 6,
      height = 4,
      units = "in",
      dpi = 300
    )
  }
)

## Offshore with curtailment ####
bmu_code <- "T_SGRWO-1"
pc_temp <- power_curve_data1f(bmu_code)
lapply(
  2:length(pc_temp),
  \(i) {
    ggsave(
      filename = sprintf(file_string, bmu_code, fig_codes[i - 1]),
      plot = pc_temp[[i]],
      width = 6,
      height = 4,
      units = "in",
      dpi = 300
    )
  }
)
## Onshore high curtailment ####
bmu_code <- "T_VKNGW-1"
pc_temp <- power_curve_data1f(bmu_code)
lapply(
  2:length(pc_temp),
  \(i) {
    ggsave(
      filename = sprintf(file_string, bmu_code, fig_codes[i - 1]),
      plot = pc_temp[[i]],
      width = 6,
      height = 4,
      units = "in",
      dpi = 300
    )
  }
)
## Onshore low curtailment ####
bmu_code <- "T_FALGW-1"
pc_temp <- power_curve_data1f(bmu_code)
lapply(
  2:length(pc_temp),
  \(i) {
    ggsave(
      filename = sprintf(file_string, bmu_code, fig_codes[i - 1]),
      plot = pc_temp[[i]],
      width = 6,
      height = 4,
      units = "in",
      dpi = 300
    )
  }
)

# apply to generic power curve to entire data ###
t0 <- "2024-01-01"
t1 <- "2024-12-31"

pwr_curv_df <- gen_adj %>%
  filter(
    between(halfHourEndTime, t0, t1),
    is.na(outageCapacity)
  ) %>%
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
  filter(!is.na(lon)) %>%
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
  class_curve <- class_curves %>% filter(class == !!class)
  class_ind <- pwr_curv_df$turb_class == class
  pwr_curv_df$power_est0[class_ind] <- approx(
    x = class_curve$wind_speed,
    y = class_curve$power_scaled,
    xout = pwr_curv_df$ws100[class_ind],
    rule = 2
  )$y *
    pwr_curv_df$capacity[class_ind]
}

bmu_codes <- c("T_FALGW-1", "T_VKNGW-1", "T_SGRWO-1", "T_HOWAO-1", "T_HOWAO-2")

lm_t_calib <- lapply(
  bmu_codes,
  \(bmu_code) {
    # filter data
    bmu_data <- pwr_curv_df %>%
      filter(grepl(bmu_code, bmUnit))

    # estimate lm calibration
    mod1_t <- brm(
      formula = potential ~ -1 + power_est0,
      data = bmu_data,
      family = student(), # <-- Student-t likelihood
      chains = 4,
      cores = 4,
      iter = 5000
    )
    list(data = bmu_data, model = mod1_t)
  }
)

saveRDS(
  lm_t_calib,
  file = file.path("~/Documents/elexon/model_objects", "lm_t_calib_v0.rds")
)
lm_t_calib[[5]]$model %>% summary()


vik_df <- lm_t_calib[[2]]$data %>%
  filter(potential > 0 | power_est0 < 0.1 * capacity)

vik_df %>%
  ggplot(aes(power_est0, potential)) +
  geom_point(col = "darkblue", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, col = "red", linetype = 2)
mod1_t <- brm(
  formula = potential ~ -1 + power_est0,
  data = vik_df,
  family = student(), # <-- Student-t likelihood
  chains = 4,
  cores = 4,
  iter = 5000
)
# mod1_t %>% summary()
lm_t_calib[[2]] <- list(data = vik_df, model = mod1_t)

file_string <- "fig/pc_scatter_%s.pdf"
lapply(
  seq_along(bmu_codes),
  \(i) {
    lm_t_calib[[i]]$data %>%
      ggplot(aes(power_est0, potential)) +
      geom_point(col = blues9[6], alpha = 0.5) +
      geom_abline(intercept = 0, slope = 1, col = "black", linetype = 2) +
      geom_abline(
        intercept = 0,
        slope = fixef(lm_t_calib[[i]]$model)[1],
        col = "darkred",
        linetype = 2
      ) +
      annotate(
        "label",
        x = mean(lm_t_calib[[i]]$data$potential),
        y = mean(lm_t_calib[[i]]$data$potential),
        label = paste0("y = x"),
        fill = "gray90",
        color = "gray20",
        hjust = 1.5
      ) +
      annotate(
        "label",
        x = mean(lm_t_calib[[i]]$data$power_est0),
        y = mean(lm_t_calib[[i]]$data$potential),
        label = paste0("y = ", round(fixef(lm_t_calib[[i]]$model)[1], 2), " x"),
        fill = "gray90",
        color = "gray20",
        hjust = -0.5
      ) +
      labs(
        x = "ERA5 derived power est. (MW)",
        y = "Elexon generation + curtailment (MW)"
      )
    ggsave(
      sprintf(file_string, bmu_codes[i]),
      width = 6,
      height = 4,
      units = "in",
      dpi = 300
    )
  }
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


gen_adj %>%
  filter(
    grepl("HOWAO-2", bmUnit),
    between(halfHourEndTime, "2019-08-08", "2019-08-10")
  ) %>% #View()
  ggplot() +
  geom_line(aes(halfHourEndTime, quantity, col = "gen")) +
  geom_line(aes(halfHourEndTime, curtailment, col = "curt")) +
  geom_line(aes(halfHourEndTime, potential, col = "potential")) +
  scale_color_aaas()
