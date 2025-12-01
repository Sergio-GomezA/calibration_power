require(arrow)
require(dplyr)
require(tidyr)
require(rnaturalearth)
require(rnaturalearthdata)
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
  file.path("data/ref_catalog_wind_2025_era.csv.gz")
)

# era 5 coordinates list
coords_tb <- read.csv("data/era5_loc_mapping.csv")
# elexon wind bmus with updated capacity
wind.bmus.alt <- read.csv("data/wind_bmu_alt.csv")

# potential energy historical summary
pot_summary <- read.csv("data/hist_pot2024.csv")

generic_pc <- fread("data/generic_powerCurves.csv.gz")
# filter one wind farm ###

class_curves <- fread("data/generic_powerCurves.csv.gz") %>%
  group_by(class) %>%
  mutate(power_scaled = power_kw / ratedPower) %>%
  ungroup()

pwr_curv_df <- read_parquet(file.path(gen_path, "power_curve_all.parquet"))

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
t0 <- "2019-01-01"
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
      select(
        bmUnit,
        matches("lon|lat"),
        site_name,
        tech_typ,
        turb_class,
        height_turb_imp
      ),
    by = c("bmUnit")
  ) %>%
  filter(!is.na(lon)) %>%
  left_join(
    era_df %>% select(time, longitude, latitude, ws100, wd100, ws10, wd10),
    by = c(
      "halfHourEndTime" = "time",
      "era5lon" = "longitude",
      "era5lat" = "latitude"
    )
  ) %>%
  # wind speed vertical interpolation
  mutate(
    ws_log = log(height_turb_imp / 10) / log(100 / 10) * (ws100 - ws10) + ws10,
    ws_h = ws100 * (height_turb_imp / 100)^(1 / 7)
  )


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
    xout = pwr_curv_df$ws_h[class_ind],
    rule = 2
  )$y *
    pwr_curv_df$capacity[class_ind]
  pwr_curv_df$power_est_w100[class_ind] <- approx(
    x = class_curve$wind_speed,
    y = class_curve$power_scaled,
    xout = pwr_curv_df$ws_h[class_ind],
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
      filter(grepl(bmu_code, bmUnit)) %>%
      filter(potential > 0 | power_est0 < 0.15 * capacity)

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
  file = file.path("~/Documents/elexon/model_objects", "lm_t_calib_v1.rds")
)

lm_t_calib[[5]]
lapply(lm_t_calib, \(x) fixef(x$model))
lapply(lm_t_calib, \(x) summary(x$model))


# vik_df <- lm_t_calib[[2]]$data %>%
#   filter(potential > 0 | power_est0 < 0.1 * capacity)

# vik_df %>%
#   ggplot(aes(power_est0, potential)) +
#   geom_point(col = "darkblue", alpha = 0.5) +
#   geom_abline(intercept = 0, slope = 1, col = "red", linetype = 2)
# mod1_t <- brm(
#   formula = potential ~ -1 + power_est0,
#   data = vik_df,
#   family = student(), # <-- Student-t likelihood
#   chains = 4,
#   cores = 4,
#   iter = 5000
# )
# # mod1_t %>% summary()
# lm_t_calib[[2]] <- list(data = vik_df, model = mod1_t)

lm_t_calib <- readRDS(
  file = file.path("~/Documents/elexon/model_objects", "lm_t_calib_v1.rds")
)
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

    # autocorrelation
    df <- lm_t_calib[[i]]$data %>%
      mutate(est_err = potential - power_est0 * fixef(lm_t_calib[[i]]$model)[1])
    acf_df <- acf(df$est_err, plot = FALSE)

    acf_tidy <- tibble(
      lag = acf_df$lag[, 1, 1],
      acf = acf_df$acf[, 1, 1]
    )
    ggplot(acf_tidy, aes(x = lag, y = acf)) +
      geom_col(width = 0.01) +
      geom_hline(yintercept = 0, color = "black") +
      theme_minimal() +
      labs(x = "Lag", y = "ACF")
    ggsave(
      sprintf("fig/lcalib_err_acf_%s.pdf", bmu_codes[i]),
      width = 6,
      height = 4,
      units = "in",
      dpi = 300
    )
  }
)


arrow::write_parquet(
  pwr_curv_df,
  file.path(gen_path, "power_curve_all.parquet")
)


## time series during curtailment ####
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


## whole data figures ####
set.seed(0)
pwr_curv_df %>%
  filter(potential > 0 | power_est0 < 0.15 * capacity) %>%
  slice_sample(n = 10000) %>%
  mutate(across(c(power_est0, potential), ~ . / capacity * 100)) %>%
  filter(power_est0 < 100, potential < 100) %>%
  ggplot(aes(power_est0, potential)) +
  geom_point(col = blues9[6], alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = 2) +
  # geom_abline(
  #   intercept = 0,
  #   slope = fixef(lm_t_calib[[i]]$model)[1],
  #   col = "darkred",
  #   linetype = 2
  # ) +
  geom_smooth(
    method = 'gam',
    formula = y ~ s(x, bs = "cs"),
    color = "darkred"
  ) +
  annotate(
    "label",
    x = mean(pwr_curv_df$potential),
    y = mean(pwr_curv_df$potential),
    label = paste0("y = x"),
    fill = "gray90",
    color = "gray20",
    hjust = 1.5
  ) +
  annotate(
    "label",
    x = 50,
    y = 25,
    label = paste0("y ~ splines(x)"),
    fill = "gray90",
    color = "gray20",
    hjust = 0
  ) +
  labs(
    x = "ERA5 derived capacity factor (%)",
    y = "Elexon capacity factor"
  )

ggsave("fig/pc_scatter_all.pdf", width = 6, height = 4, units = "in", dpi = 300)

# map
uk_ie <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(admin %in% c("United Kingdom", "Ireland")) %>%
  st_transform(crs = 4326)
errdat <- pwr_curv_df %>%
  filter(potential > 0 | power_est0 < 0.15 * capacity) %>%
  slice_sample(n = 10000) %>%
  mutate(across(c(power_est0, potential), ~ . / capacity * 100)) %>%
  filter(power_est0 < 100, potential < 100) %>%
  group_by(lon, lat, tech_typ) %>%
  summarise(
    across(c(power_est0, potential), mean, na.rm = TRUE)
  ) %>%
  mutate(mape = abs(potential - power_est0) / potential) %>%
  st_as_sf(
    coords = c("lon", "lat"),
    crs = 4326
  )

ggplot() +
  geom_sf(data = uk_ie, fill = "gray95", color = "gray70") + # map background
  geom_sf(
    data = errdat,
    aes(geometry = geometry, color = tech_typ, size = mape)
  ) +
  scale_color_manual(values = mypalette) +
  theme_map() +
  labs(col = "", size = "MAPE") +
  theme(legend.position = "right")
ggsave(
  "fig/pc_mape_map_all.pdf",
)


# comparison after agregation ####
pwr_curv_df %>% names()
GB_df <- pwr_curv_df %>%
  group_by(tech_typ, halfHourEndTime) %>%
  summarise(
    ws_h_wmean = sum(ws_h * capacity),
    across(
      c(power_est0, potential, capacity),
      sum
    ),
    across(c(ws_h), mean, .names = "{.col}_mean")
  ) %>%
  mutate(
    across(
      c(power_est0, potential),
      ~ . / capacity,
      .names = "norm_{.col}"
    ),
    ws_h_wmean = ws_h_wmean / capacity,
    month = factor(month(halfHourEndTime)),
    hour = factor(hour(halfHourEndTime))
  )


lm0 <- lm(
  formula = norm_potential ~ norm_power_est0,
  data = GB_df
)
summary(lm0)
lm1 <- lm(
  formula = norm_potential ~ tech_typ *
    norm_power_est0,
  data = GB_df
)
summary(lm1)
lm2 <- lm(
  formula = norm_potential ~ tech_typ *
    norm_power_est0 *
    month,
  data = GB_df
)
summary(lm2)
# lm_t <- brm(
#   formula = norm_potential ~ +norm_power_est0,
#   data = GB_df,
#   family = student(), # <-- Student-t likelihood
#   chains = 4,
#   cores = 4,
#   iter = 5000
# )

saveRDS(
  lm_t,
  file = file.path("~/Documents/elexon/model_objects", "lm_t_calib_all.rds")
)
lm_t <- readRDS(
  file.path("~/Documents/elexon/model_objects", "lm_t_calib_all.rds")
)
plot(lm_t)
summary(lm_t)

fixef(lm_t)[, 1]
lm0$coefficients

## scatter #####

GB_df %>%
  # slice_sample(n = 10000) %>%
  ggplot() +
  geom_point(
    aes(norm_power_est0, norm_potential),
    color = blues9[5],
    # fill = NA,
    alpha = 0.2,
    pch = 16
  ) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  geom_abline(
    aes(intercept = lm0$coefficients[1], slope = lm0$coefficients[2]),
    linetype = 2,
    col = "darkred"
  ) +
  annotate(
    "label",
    x = mean(GB_df$norm_potential),
    y = mean(GB_df$norm_potential),
    label = paste0("y = x"),
    fill = "gray90",
    color = "gray20",
    hjust = 1.5
  ) +
  annotate(
    "label",
    x = mean(GB_df$norm_power_est0),
    y = mean(GB_df$norm_power_est0) * lm0$coefficients[2] + lm0$coefficients[1],
    label = sprintf(
      "y = %0.2f + %0.2f x",
      lm0$coefficients[1],
      lm0$coefficients[2]
    ),
    fill = "gray90",
    color = "gray20",
    hjust = 0
  ) +
  labs(
    x = "ERA5 derived capacity factor (%)",
    y = "Elexon capacity factor"
  )

ggsave("fig/gb_aggr_scatter_all.png", width = 6, height = 4)

GB_df %>%
  ggplot(aes(norm_power_est0, norm_potential)) +
  stat_bin2d(aes(fill = (after_stat(count))), bins = 100) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  geom_abline(
    aes(intercept = lm0$coefficients[1], slope = lm0$coefficients[2]),
    linetype = 2,
    col = "darkred"
  ) +
  annotate(
    "label",
    x = mean(GB_df$norm_potential),
    y = mean(GB_df$norm_potential),
    label = paste0("y = x"),
    fill = "gray90",
    color = "gray20",
    hjust = 1.5
  ) +
  annotate(
    "label",
    x = mean(GB_df$norm_power_est0),
    y = mean(GB_df$norm_power_est0) * lm0$coefficients[2] + lm0$coefficients[1],
    label = sprintf(
      "y = %0.2f + %0.2f x",
      lm0$coefficients[1],
      lm0$coefficients[2]
    ),
    fill = "gray90",
    color = "gray20",
    hjust = 0
  ) +
  labs(
    x = "ERA5 derived capacity factor (%)",
    y = "Elexon capacity factor"
  ) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    breaks = c(3, 30, 300)
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.1, 0.78),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA)
  )
ggsave("fig/gb_aggr_bindens_all.png", width = 6, height = 4)

GB_df %>%
  ggplot(aes(
    pmax(0, lm2$fitted.values),
    norm_potential
  )) +
  stat_bin2d(aes(fill = (after_stat(count))), bins = 100) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  annotate(
    "label",
    x = mean(GB_df$norm_potential),
    y = mean(GB_df$norm_potential),
    label = paste0("y = x"),
    fill = "gray90",
    color = "gray20",
    hjust = 1.5
  ) +
  labs(
    x = "ERA5 calibrated CF(%)",
    y = "Elexon capacity factor"
  ) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    breaks = c(3, 30, 300)
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.1, 0.78),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA)
  )
ggsave("fig/gb_calib_bindens_all.png", width = 6, height = 4)


with(GB_df, ModelMetrics::rmse(norm_potential, lm0$fitted.values))
with(GB_df, ModelMetrics::rmse(norm_potential, lm2$fitted.values))


### temporal correlation ####
pdf("fig/lcalib_err_acf_offshore.pdf")
acf_df <- acf(
  lm2$residuals[GB_df$tech_typ == "Wind Offshore"],
  plot = TRUE,
  main = "Wind Offshore"
)
dev.off()

pdf("fig/lcalib_err_acf_onshore.pdf")
acf_df <- acf(
  lm2$residuals[GB_df$tech_typ == "Wind Onshore"],
  plot = TRUE,
  main = "Wind Onshore"
)
dev.off()


### scatter plots by type #####
line_df <- GB_df %>%
  group_by(tech_typ) %>%
  summarise(
    intercept = coef(lm(norm_potential ~ norm_power_est0))[1],
    slope = coef(lm(norm_potential ~ norm_power_est0))[2]
  ) %>%
  mutate(
    x = mean(GB_df$norm_power_est0, na.rm = TRUE),
    y = intercept + slope * x,
    label = sprintf("y = %.2f + %.2f x", intercept, slope)
  )
GB_df %>%
  # slice_sample(n = 10000) %>%
  ggplot() +
  geom_point(
    aes(norm_power_est0, norm_potential, col = tech_typ),
    # fill = NA,
    alpha = 0.2,
    pch = 16
  ) +
  facet_wrap(~tech_typ, nrow = 2) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  geom_abline(
    aes(intercept = intercept, slope = slope),
    data = line_df,
    linetype = 2,
    col = "darkred"
  ) +
  geom_label(
    data = line_df,
    aes(x = x, y = y, label = label),
    fill = "gray90",
    color = "gray20",
    hjust = 0
  ) +
  labs(
    x = "ERA5 derived capacity factor (%)",
    y = "Elexon capacity factor"
  ) +
  scale_color_manual(values = mypalette) +
  theme(legend.position = "none")
ggsave("fig/gb_aggr_scatter_typ_all.png", width = 6, height = 4)

### scatter plots by month #####
line_df <- GB_df %>%
  group_by(tech_typ, month) %>%
  summarise(
    intercept = coef(lm(norm_potential ~ norm_power_est0))[1],
    slope = coef(lm(norm_potential ~ norm_power_est0))[2]
  ) %>%
  mutate(
    x = mean(GB_df$norm_power_est0, na.rm = TRUE),
    y = intercept + slope * x,
    label = sprintf("y = %.2f + %.2f x", intercept, slope)
  )
GB_df %>%
  filter(tech_typ == "Wind Onshore") %>%
  ggplot() +
  geom_point(
    aes(norm_power_est0, norm_potential, col = tech_typ),
    # fill = NA,
    alpha = 0.2,
    pch = 16
  ) +
  facet_wrap(~month, ncol = 4) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  geom_abline(
    aes(intercept = intercept, slope = slope),
    data = line_df %>% filter(tech_typ == "Wind Onshore"),
    linetype = 2,
    col = "darkred"
  ) +
  geom_text(
    data = line_df %>% filter(tech_typ == "Wind Onshore"),
    aes(x = 0.2, y = 0.05, label = label),
    color = "gray20",
    size = 3.5,
    hjust = 0
  ) +
  labs(
    x = "ERA5 derived capacity factor (%)",
    y = "Elexon capacity factor"
  ) +
  scale_color_manual(values = mypalette[2]) +
  theme(legend.position = "none", axis.text.x = element_text(size = 7))
ggsave("fig/gb_aggr_scatter_on_month_all.png", width = 8, height = 6)

GB_df %>%
  filter(tech_typ == "Wind Offshore") %>%
  ggplot() +
  geom_point(
    aes(norm_power_est0, norm_potential, col = tech_typ),
    # fill = NA,
    alpha = 0.2,
    pch = 16
  ) +
  facet_wrap(~month, ncol = 4) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  geom_abline(
    aes(intercept = intercept, slope = slope),
    data = line_df %>% filter(tech_typ == "Wind Offshore"),
    linetype = 2,
    col = "darkred"
  ) +
  geom_text(
    data = line_df %>% filter(tech_typ == "Wind Offshore"),
    aes(x = 0.2, y = 0.05, label = label),
    color = "gray20",
    size = 3.5,
    hjust = 0
  ) +
  labs(
    x = "ERA5 derived capacity factor (%)",
    y = "Elexon capacity factor"
  ) +
  scale_color_manual(values = mypalette[1]) +
  theme(legend.position = "none", axis.text.x = element_text(size = 7))
ggsave("fig/gb_aggr_scatter_off_month_all.png", width = 8, height = 6)

### scatter plots by hour #####
line_df <- GB_df %>%
  group_by(tech_typ, hour) %>%
  summarise(
    intercept = coef(lm(norm_potential ~ norm_power_est0))[1],
    slope = coef(lm(norm_potential ~ norm_power_est0))[2]
  ) %>%
  mutate(
    x = mean(GB_df$norm_power_est0, na.rm = TRUE),
    y = intercept + slope * x,
    label = sprintf("y = %.2f + %.2f x", intercept, slope)
  )
GB_df %>%
  filter(tech_typ == "Wind Onshore") %>%
  ggplot() +
  geom_point(
    aes(norm_power_est0, norm_potential, col = tech_typ),
    # fill = NA,
    alpha = 0.2,
    pch = 16
  ) +
  facet_wrap(~hour, ncol = 6) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  geom_abline(
    aes(intercept = intercept, slope = slope),
    data = line_df %>% filter(tech_typ == "Wind Onshore"),
    linetype = 2,
    col = "darkred"
  ) +
  geom_text(
    data = line_df %>% filter(tech_typ == "Wind Onshore"),
    aes(x = 0.2, y = 0.05, label = label),
    color = "gray20",
    size = 2,
    hjust = 0
  ) +
  labs(
    x = "ERA5 derived capacity factor (%)",
    y = "Elexon capacity factor"
  ) +
  scale_color_manual(values = mypalette[2]) +
  theme(legend.position = "none", axis.text.x = element_text(size = 5))
ggsave("fig/gb_aggr_scatter_on_hour_all.png", width = 8, height = 6)

GB_df %>%
  filter(tech_typ == "Wind Offshore") %>%
  ggplot() +
  geom_point(
    aes(norm_power_est0, norm_potential, col = tech_typ),
    # fill = NA,
    alpha = 0.2,
    pch = 16
  ) +
  facet_wrap(~hour, ncol = 6) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  geom_abline(
    aes(intercept = intercept, slope = slope),
    data = line_df %>% filter(tech_typ == "Wind Offshore"),
    linetype = 2,
    col = "darkred"
  ) +
  geom_text(
    data = line_df %>% filter(tech_typ == "Wind Offshore"),
    aes(x = 0.2, y = 0.05, label = label),
    color = "gray20",
    size = 2,
    hjust = 0
  ) +
  labs(
    x = "ERA5 derived capacity factor (%)",
    y = "Elexon capacity factor"
  ) +
  scale_color_manual(values = mypalette[1]) +
  theme(legend.position = "none", axis.text.x = element_text(size = 5))
ggsave("fig/gb_aggr_scatter_off_hour_all.png", width = 8, height = 6)


## densities ####
GB_df %>%
  # slice_sample(n = 10000) %>%
  ggplot() +
  geom_density2d_filled(
    aes(norm_power_est0, norm_potential),
    # color = blues9[5],
    # alpha = 0.5
  ) +
  scale_fill_manual(values = blues9)

GB_df %>%
  # slice_sample(n = 10000) %>%
  pivot_longer(cols = matches("norm_")) %>%
  ggplot() +
  geom_density(
    aes(x = value, fill = name),
    # color = blues9[5],
    alpha = 0.5
  ) +
  facet_wrap(~tech_typ, nrow = 2) +
  scale_fill_manual(
    values = pal_lancet()(2),
    labels = c("elexon", "era5 estimate")
  ) +
  labs(x = "Capacity factor", fill = "") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.5, 0.9),
    # legend.box.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA)
  )

ggsave("fig/density_comparison_all.pdf", width = 6, height = 4)

### dens plots by season ####

GB_df %>%
  filter(tech_typ == "Wind Offshore") %>%
  # slice_sample(n = 10000) %>%
  pivot_longer(cols = matches("norm_")) %>%
  ggplot() +
  geom_density(
    aes(x = value, fill = name),
    # color = blues9[5],
    alpha = 0.5
  ) +
  facet_wrap(~month, ncol = 4) +
  scale_fill_manual(
    values = pal_lancet()(2),
    labels = c("elexon", "era5 estimate")
  ) +
  labs(x = "Capacity factor", fill = "") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.35, 0.95),
    # legend.box.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7)
  )
ggsave("fig/cf_density_comparison_off_month_all.pdf", width = 8, height = 6)
GB_df %>%
  filter(tech_typ == "Wind Onshore") %>%
  # slice_sample(n = 10000) %>%
  pivot_longer(cols = matches("norm_")) %>%
  ggplot() +
  geom_density(
    aes(x = value, fill = name),
    # color = blues9[5],
    alpha = 0.5
  ) +
  facet_wrap(~month, ncol = 4) +
  scale_fill_manual(
    values = pal_lancet()(2),
    labels = c("elexon", "era5 estimate")
  ) +
  labs(x = "Capacity factor", fill = "") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.35, 0.95),
    # legend.box.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7)
  )
ggsave("fig/cf_density_comparison_on_month_all.pdf", width = 8, height = 6)

### dens by hour ####
GB_df %>%
  filter(tech_typ == "Wind Offshore") %>%
  # slice_sample(n = 10000) %>%
  pivot_longer(cols = matches("norm_")) %>%
  ggplot() +
  geom_density(
    aes(x = value, fill = name),
    # color = blues9[5],
    alpha = 0.5
  ) +
  facet_wrap(~hour, ncol = 6) +
  scale_fill_manual(
    values = pal_lancet()(2),
    labels = c("elexon", "era5 estimate")
  ) +
  labs(x = "Capacity factor", fill = "") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.35, 0.95),
    # legend.box.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7)
  )
ggsave("fig/cf_density_comparison_off_hour_all.pdf", width = 8, height = 6)
GB_df %>%
  filter(tech_typ == "Wind Onshore") %>%
  # slice_sample(n = 10000) %>%
  pivot_longer(cols = matches("norm_")) %>%
  ggplot() +
  geom_density(
    aes(x = value, fill = name),
    # color = blues9[5],
    alpha = 0.5
  ) +
  facet_wrap(~hour, ncol = 6) +
  scale_fill_manual(
    values = pal_lancet()(2),
    labels = c("elexon", "era5 estimate")
  ) +
  labs(x = "Capacity factor", fill = "") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.35, 0.95),
    # legend.box.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7)
  )
ggsave("fig/cf_density_comparison_on_hour_all.pdf", width = 8, height = 6)

### power curve plots by season ####

source("aux_funct.R")

curve_GB <- est_pwr_curv(
  GB_df,
  avg_fun = median,
  n_bins = 50,
  quantile_bins = TRUE,
  plot = FALSE
)
speeds <- seq(2, 25, length.out = 100)
pc_GB <- data.frame(
  ws_mean = speeds
) %>%
  mutate(
    power_mean = approx(
      x = curve_GB$ws_mean,
      y = curve_GB$power_mean,
      xout = ws_mean,
      rule = 2 # flat extrapolation
    )$y
  )
GB_df %>%
  ggplot(aes(ws_h_wmean, norm_potential, col = "observed")) +
  geom_point(alpha = 0.2) +
  geom_line(
    aes(ws_mean, power_mean, col = "power curve est."),
    data = curve_GB
  ) +
  scale_color_manual(
    values = c(
      "observed" = blues9[7],
      "power curve est." = "darkred",
      "outages" = "darkorange"
    )
  ) +
  theme(
    legend.position = "bottom"
  ) +
  labs(col = "", x = "wind speed", y = "generation (MW)") +
  scale_x_continuous(breaks = 5 * (0:5))


ggsave("fig/pc_scatter_GB_valnout_all.png", width = 6, height = 4)

GB_df %>%
  ggplot(aes(ws_h_wmean, norm_potential)) +
  stat_bin2d(aes(fill = (after_stat(count))), bins = 100) +
  geom_line(
    aes(ws_mean, power_mean, col = "power curve est."),
    data = curve_GB
  ) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    breaks = c(3, 30, 300)
  ) +
  scale_color_manual(values = c("power curve est." = "darkred")) +
  labs(x = "wind speed", y = "capacity factor %", col = "") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.2, 0.7),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA)
  )
ggsave("fig/pc_bindens_GB_valnout_all.png", width = 6, height = 4)

curve_GB_typ <- lapply(
  c("Wind Offshore", "Wind Onshore"),
  \(typ) {
    curve_month <- est_pwr_curv(
      GB_df %>% filter(tech_typ == typ),
      avg_fun = median,
      n_bins = 30,
      quantile_bins = TRUE,
      plot = FALSE
    ) %>%
      mutate(tech_typ = typ)
  }
) %>%
  bind_rows()
GB_df %>%
  ggplot(aes(ws_h_wmean, norm_potential)) +
  stat_bin2d(aes(fill = (after_stat(count))), bins = 100) +
  geom_line(
    aes(ws_mean, power_mean, col = "power curve est."),
    data = curve_GB
  ) +
  geom_line(
    aes(ws_mean, power_mean, col = "type est."),
    data = curve_GB_typ
  ) +
  facet_wrap(~tech_typ, nrow = 2) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    breaks = c(3, 30, 300)
  ) +
  scale_color_manual(
    values = c("power curve est." = "darkred", "type est." = "darkorange")
  ) +
  labs(x = "wind speed", y = "capacity factor %", col = "") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.85),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    # shrink keys
    legend.key.size = unit(0.4, "lines"),
    legend.spacing.y = unit(0, "lines"), # reduce vertical gap between legends
    legend.box.spacing = unit(0, "lines"), # reduce space around the legend box
    # shrink text
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )
ggsave("fig/pc_bindens_GB_valnout_all_typ.png", width = 6, height = 6)
curve_GB_month <- lapply(
  1:12,
  \(m) {
    curve_month <- est_pwr_curv(
      GB_df %>% filter(month == m),
      avg_fun = median,
      n_bins = 30,
      quantile_bins = TRUE,
      plot = FALSE
    ) %>%
      mutate(month = m)
  }
) %>%
  bind_rows()

GB_df %>%
  ggplot(aes(ws_h_wmean, norm_potential, col = "observed")) +
  geom_point(alpha = 0.2) +
  geom_line(
    aes(ws_mean, power_mean, col = "global pc. est."),
    data = curve_GB
  ) +
  geom_line(
    aes(ws_mean, power_mean, col = "monthly est."),
    data = curve_GB_month
  ) +
  scale_color_manual(
    values = c(
      "observed" = blues9[7],
      "global pc. est." = "darkred",
      "monthly est." = "darkorange"
    )
  ) +
  facet_wrap(~month) +
  theme(
    legend.position = "bottom"
  ) +
  labs(col = "", x = "wind speed", y = "generation (MW)") +
  scale_x_continuous(breaks = 5 * (0:5))
ggsave("fig/pc_scatter_GB_valnout_month.png", width = 8, height = 6)

GB_df %>%
  ggplot(aes(ws_h_wmean, norm_potential)) +
  stat_bin2d(aes(fill = (after_stat(count))), bins = 100) +
  geom_line(
    aes(ws_mean, power_mean, col = "PC est."),
    data = curve_GB
  ) +
  geom_line(
    aes(ws_mean, power_mean, col = "monthly est."),
    data = curve_GB_month
  ) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "freq.",
    breaks = c(1, 10, 50)
  ) +
  scale_color_manual(values = c("PC est." = "darkred")) +
  facet_wrap(~month) +
  labs(x = "wind speed", y = "capacity factor %", col = "") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.2, 0.45),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    # shrink keys
    legend.key.size = unit(0.4, "lines"),
    legend.spacing.y = unit(0, "lines"), # reduce vertical gap between legends
    legend.box.spacing = unit(0, "lines"), # reduce space around the legend box
    # shrink text
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )
ggsave("fig/pc_bindens_GB_valnout_month_all.png", width = 8, height = 6)
