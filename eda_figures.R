fig_path <- "fig"

# time series thining for spatial correlation ######
n.days <- 7
time_skips <- 29
# slice data to 12 hour time steps
gb_thin <- pwr_curv_df %>%
  mutate(
    site_id = as.integer(factor(site_name))
  ) %>%
  # rename(halfHourEndTime = time) %>%
  group_by(bmUnit) %>%
  arrange(halfHourEndTime) %>%
  mutate(hour = row_number()) %>%
  # filter(hour %% (24 * n.days) == 9) %>%
  filter(hour %% time_skips == 0) %>%
  mutate(
    norm_potential = potential / capacity,
    norm_power_est0 = power_est0 / capacity,
    error0 = norm_potential - norm_power_est0
  ) %>%
  group_by(lon, lat, halfHourEndTime) %>%
  summarise(
    site_name = first(site_name),
    site_id = first(site_id),
    across(c(norm_potential, norm_power_est0, error0), sum),
    ws_h_wmean = sum(ws_h * capacity) / sum(capacity),
    capacity = sum(capacity),
    .groups = "drop"
  )

max_lag <- 20

pacf_df <- gb_thin %>%
  group_by(site_name, site_id) %>%
  summarise(
    n = sum(!is.na(norm_potential)),
    pacf = list(
      pacf(norm_potential, lag.max = max_lag, plot = FALSE)$acf
    ),
    lag = list(
      seq_along(pacf(norm_potential, lag.max = max_lag, plot = FALSE)$acf)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    ci = 1.96 / sqrt(n)
  ) %>%
  unnest(c(lag, pacf)) %>%
  mutate(
    significant = abs(pacf) >= ci
  )
ggplot(pacf_df, aes(lag, site_id, group = site_id, fill = significant)) +
  geom_tile() +
  scale_fill_manual(values = c("grey90", "red")) +
  theme_minimal()


ggplot(pacf_df, aes(lag, site_id, fill = pacf)) +
  geom_tile() +
  scale_fill_gradient2(
    midpoint = 0,
    low = "blue",
    mid = "white",
    high = "red",
    na.value = "grey90"
  ) +
  theme_minimal() +
  labs(x = "Lag", y = "Time series", fill = "PACF")


# spatial correlation by distance #####

## wind generation ######
sp_corr <- spatial_corr_by_distance_fast(
  gb_thin,
  value_col = "norm_potential",
  n_bins = 15,
  bin_type = "quantile",
  bands = c(0.05, 0.95)
)

corr_df <- sp_corr$summary

ggplot(corr_df, aes(x = dist_mean, y = corr_mean)) +
  geom_ribbon(
    aes(ymin = corr_lower, ymax = corr_upper, fill = "90% CI"),
    # alpha = 0.2,
    # fill = blues9[3]
  ) +
  geom_line(color = blues9[5]) +
  geom_point(color = blues9[7]) +
  # geom_line(aes(y = threshold_95), linetype = "dashed", color = "red") +
  # geom_line(aes(y = -threshold_95), linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Distance (km)",
    y = "Correlation",
    # title = "Spatial correlation vs distance with variability and significance"
  ) +
  theme_minimal() +
  scale_fill_manual(values = blues9[3], name = "") +
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8))

ggsave("fig/norm_pot_corr_dist_90CI.pdf", width = 6, height = 4)


outliers_df <- sp_corr$pairs %>%
  filter(n_obs_per_pair >= 50) %>%
  group_by(bin) %>%
  mutate(
    Q1 = quantile(corr, 0.25, na.rm = TRUE),
    Q3 = quantile(corr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower = Q1 - 1.5 * IQR,
    upper = Q3 + 1.5 * IQR,
    is_outlier = corr < lower | corr > upper
  ) %>%
  filter(is_outlier)

bin_centers <- sp_corr$pairs %>%
  group_by(bin) %>%
  summarise(dist_mean = mean(dist_km), .groups = "drop")

outliers_df <- outliers_df %>%
  left_join(bin_centers, by = "bin")

ggplot(corr_df, aes(x = dist_mean, y = corr_mean)) +
  geom_ribbon(aes(ymin = corr_lower, ymax = corr_upper, fill = "90% CI")) +
  geom_line(color = blues9[5]) +
  geom_point(color = blues9[7]) +
  geom_point(
    data = outliers_df,
    aes(x = dist_mean, y = corr, col = "outliers"),
    # color = "red",
    alpha = 0.6,
    size = 1.5
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Distance (km)", y = "Correlation") +
  theme_minimal() +
  scale_fill_manual(values = blues9[3], name = "") +
  scale_color_manual(values = "darkred", name = "") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8),
    legend.box.spacing = unit(0, "pt"),
    legend.spacing.y = unit(0, "pt"),
    legend.key.height = unit(8, "pt")
  )
ggsave("fig/norm_pot_corr_dist_90CI_outliers.pdf", width = 6, height = 4)

## wind generation error ######
sp_corr_e <- spatial_corr_by_distance_fast(
  gb_thin,
  value_col = "error0",
  n_bins = 15,
  bin_type = "quantile",
  bands = c(0.05, 0.95)
)

corr_df_e <- sp_corr_e$summary

ggplot(corr_df_e, aes(x = dist_mean, y = corr_mean)) +
  geom_ribbon(
    aes(ymin = corr_lower, ymax = corr_upper, fill = "90% CI"),
    # alpha = 0.2,
    # fill = blues9[3]
  ) +
  geom_line(color = blues9[5]) +
  geom_point(color = blues9[7]) +
  # geom_line(aes(y = threshold_95), linetype = "dashed", color = "red") +
  # geom_line(aes(y = -threshold_95), linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Distance (km)",
    y = "Correlation",
    # title = "Spatial correlation vs distance with variability and significance"
  ) +
  theme_minimal() +
  scale_fill_manual(values = blues9[3], name = "") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8),
    legend.box.spacing = unit(0, "pt"),
    legend.spacing.y = unit(0, "pt"),
    legend.key.height = unit(8, "pt")
  )

ggsave("fig/error0_corr_dist_90CI.pdf", width = 6, height = 4)


outliers_df_e <- sp_corr_e$pairs %>%
  filter(n_obs_per_pair >= 50) %>%
  group_by(bin) %>%
  mutate(
    Q1 = quantile(corr, 0.25, na.rm = TRUE),
    Q3 = quantile(corr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower = Q1 - 1.5 * IQR,
    upper = Q3 + 1.5 * IQR,
    is_outlier = corr < lower | corr > upper
  ) %>%
  filter(is_outlier)

bin_centers_e <- sp_corr_e$pairs %>%
  group_by(bin) %>%
  summarise(dist_mean = mean(dist_km), .groups = "drop")

outliers_df_e <- outliers_df_e %>%
  left_join(bin_centers_e, by = "bin")

ggplot(corr_df_e, aes(x = dist_mean, y = corr_mean)) +
  geom_ribbon(aes(ymin = corr_lower, ymax = corr_upper, fill = "90% CI")) +
  geom_line(color = blues9[5]) +
  geom_point(color = blues9[7]) +
  geom_point(
    data = outliers_df_e,
    aes(x = dist_mean, y = corr, col = "outliers"),
    # color = "red",
    alpha = 0.6,
    size = 1.5
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Distance (km)", y = "Correlation") +
  theme_minimal() +
  scale_fill_manual(values = blues9[3], name = "") +
  scale_color_manual(values = "darkred", name = "") +
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8)) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8),
    legend.box.spacing = unit(0, "pt"),
    legend.spacing.y = unit(0, "pt"),
    legend.key.height = unit(8, "pt")
  )
ggsave("fig/error0_corr_dist_90CI_outliers.pdf", width = 6, height = 4)

## wind generation generic PC ####

sp_corr_est <- spatial_corr_by_distance_fast(
  gb_thin,
  value_col = "norm_power_est0",
  n_bins = 15,
  bin_type = "quantile",
  bands = c(0.05, 0.95)
)

corr_df_est <- sp_corr_est$summary

ggplot(corr_df_est, aes(x = dist_mean, y = corr_mean)) +
  geom_ribbon(
    aes(ymin = corr_lower, ymax = corr_upper, fill = "90% CI"),
    # alpha = 0.2,
    # fill = blues9[3]
  ) +
  geom_line(color = blues9[5]) +
  geom_point(color = blues9[7]) +
  # geom_line(aes(y = threshold_95), linetype = "dashed", color = "red") +
  # geom_line(aes(y = -threshold_95), linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Distance (km)",
    y = "Correlation",
    # title = "Spatial correlation vs distance with variability and significance"
  ) +
  theme_minimal() +
  scale_fill_manual(values = blues9[3], name = "") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8),
    legend.box.spacing = unit(0, "pt"),
    legend.spacing.y = unit(0, "pt"),
    legend.key.height = unit(8, "pt")
  )

ggsave("fig/power_est_corr_dist_90CI.pdf", width = 6, height = 4)

outliers_df_est <- sp_corr_est$pairs %>%
  filter(n_obs_per_pair >= 50) %>%
  group_by(bin) %>%
  mutate(
    Q1 = quantile(corr, 0.25, na.rm = TRUE),
    Q3 = quantile(corr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower = Q1 - 1.5 * IQR,
    upper = Q3 + 1.5 * IQR,
    is_outlier = corr < lower | corr > upper
  ) %>%
  filter(is_outlier)

bin_centers_est <- sp_corr_est$pairs %>%
  group_by(bin) %>%
  summarise(dist_mean = mean(dist_km), .groups = "drop")

outliers_df_est <- outliers_df_est %>%
  left_join(bin_centers_est, by = "bin")

ggplot(corr_df_est, aes(x = dist_mean, y = corr_mean)) +
  geom_ribbon(aes(ymin = corr_lower, ymax = corr_upper, fill = "90% CI")) +
  geom_line(color = blues9[5]) +
  geom_point(color = blues9[7]) +
  geom_point(
    data = outliers_df_est,
    aes(x = dist_mean, y = corr, col = "outliers"),
    # color = "red",
    alpha = 0.6,
    size = 1.5
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Distance (km)", y = "Correlation") +
  theme_minimal() +
  scale_fill_manual(values = blues9[3], name = "") +
  scale_color_manual(values = "darkred", name = "") +
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8)) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8),
    legend.box.spacing = unit(0, "pt"),
    legend.spacing.y = unit(0, "pt"),
    legend.key.height = unit(8, "pt")
  )
ggsave("fig/power_est_corr_dist_90CI_outliers.pdf", width = 6, height = 4)

## wind speed ####

sp_corr_ws <- spatial_corr_by_distance_fast(
  gb_thin,
  value_col = "ws_h_wmean",
  n_bins = 15,
  bin_type = "quantile",
  bands = c(0.05, 0.95)
)

corr_df_ws <- sp_corr_ws$summary

ggplot(corr_df_ws, aes(x = dist_mean, y = corr_mean)) +
  geom_ribbon(
    aes(ymin = corr_lower, ymax = corr_upper, fill = "90% CI"),
    # alpha = 0.2,
    # fill = blues9[3]
  ) +
  geom_line(color = blues9[5]) +
  geom_point(color = blues9[7]) +
  # geom_line(aes(y = threshold_95), linetype = "dashed", color = "red") +
  # geom_line(aes(y = -threshold_95), linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Distance (km)",
    y = "Correlation",
    # title = "Spatial correlation vs distance with variability and significance"
  ) +
  theme_minimal() +
  scale_fill_manual(values = blues9[3], name = "") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8),
    legend.box.spacing = unit(0, "pt"),
    legend.spacing.y = unit(0, "pt"),
    legend.key.height = unit(8, "pt")
  )

ggsave("fig/wind_speed_corr_dist_90CI.pdf", width = 6, height = 4)

outliers_df_ws <- sp_corr_ws$pairs %>%
  filter(n_obs_per_pair >= 50) %>%
  group_by(bin) %>%
  mutate(
    Q1 = quantile(corr, 0.25, na.rm = TRUE),
    Q3 = quantile(corr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower = Q1 - 1.5 * IQR,
    upper = Q3 + 1.5 * IQR,
    is_outlier = corr < lower | corr > upper
  ) %>%
  filter(is_outlier)

bin_centers_ws <- sp_corr_ws$pairs %>%
  group_by(bin) %>%
  summarise(dist_mean = mean(dist_km), .groups = "drop")

outliers_df_ws <- outliers_df_ws %>%
  left_join(bin_centers_ws, by = "bin")

ggplot(corr_df_ws, aes(x = dist_mean, y = corr_mean)) +
  geom_ribbon(aes(ymin = corr_lower, ymax = corr_upper, fill = "90% CI")) +
  geom_line(color = blues9[5]) +
  geom_point(color = blues9[7]) +
  geom_point(
    data = outliers_df_ws,
    aes(x = dist_mean, y = corr, col = "outliers"),
    # color = "red",
    alpha = 0.6,
    size = 1.5
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Distance (km)", y = "Correlation") +
  theme_minimal() +
  scale_fill_manual(values = blues9[3], name = "") +
  scale_color_manual(values = "darkred", name = "") +
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8)) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8),
    legend.box.spacing = unit(0, "pt"),
    legend.spacing.y = unit(0, "pt"),
    legend.key.height = unit(8, "pt")
  )
ggsave("fig/wind_speed_corr_dist_90CI_outliers.pdf", width = 6, height = 4)

# circulation patterns and error #####
gb_daily_df %>%
  ggplot(aes(x = tech_typ, y = err)) +
  geom_boxplot()


gb_daily_df %>%
  ggplot(aes(x = nao_group, y = err, group = nao_group)) +
  geom_boxplot() +
  facet_wrap(~tech_typ) +
  theme_bw()


gb_daily_df %>%
  ggplot(aes(x = err, y = nao_group, group = nao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "NAO index") +
  theme(legend.position = "none")

gb_daily_df %>%
  ggplot(aes(x = norm_potential, y = nao_group, group = nao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "NAO index") +
  theme(legend.position = "none")


gb_daily_df %>%
  ggplot(aes(x = err, y = nao_group, group = nao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "NAO index") +
  theme(legend.position = "none")


## AO
gb_daily_df %>%
  ggplot(aes(x = err, y = ao_group, group = ao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "AO index") +
  theme(legend.position = "none")

gb_daily_df %>%
  ggplot(aes(x = norm_potential, y = ao_group, group = ao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "AO index") +
  theme(legend.position = "none")

gb_daily_df %>%
  ggplot(aes(x = err, y = ao_group, group = ao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "AO index") +
  theme(legend.position = "none")

### Wind, NAO, AO, EA, SCAN

library(GGally)


vars <- gb_monthly_df[, c("ws_h_wmean", "nao", "ao", "ea", "scan")]

ggpairs(
  vars,
  diag = list(continuous = wrap("densityDiag")),
  lower = list(continuous = wrap("points", alpha = 0.4, size = 0.8)),
  upper = list(continuous = wrap("cor", size = 3))
)

## circulation patterns and generation ------------------------------
# NAO
gb_monthly_df %>%
  ggplot(aes(x = norm_potential, y = nao_group, group = nao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "wind generation", y = "NAO index") +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "nao_generation.pdf"
  ),
  width = 6,
  height = 4
)

# AO
gb_monthly_df %>%
  ggplot(aes(x = norm_potential, y = ao_group, group = ao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "wind generation", y = "AO index") +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "ao_generation.pdf"
  ),
  width = 6,
  height = 4
)
# EA
gb_monthly_df %>%
  ggplot(aes(x = norm_potential, y = ea_group, group = ea_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "wind generation", y = "EA index") +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "ea_generation.pdf"
  ),
  width = 6,
  height = 4
)
# SCAN
gb_monthly_df %>%
  ggplot(aes(x = norm_potential, y = scan_group, group = scan_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "wind generation", y = "SCAN index") +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "scan_generation.pdf"
  ),
  width = 6,
  height = 4
)

## circulation patterns and error ------------------------------
# NAO
gb_monthly_df %>%
  ggplot(aes(x = err, y = nao_group, group = nao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "NAO index") +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "nao_error.pdf"
  ),
  width = 6,
  height = 4
)

# AO
gb_monthly_df %>%
  ggplot(aes(x = err, y = ao_group, group = ao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "AO index") +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "ao_error.pdf"
  ),
  width = 6,
  height = 4
)
# EA

gb_monthly_df %>%
  ggplot(aes(x = err, y = ea_group, group = ea_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "EA index") +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "ea_error.pdf"
  ),
  width = 6,
  height = 4
)
# SCAN
gb_monthly_df %>%
  ggplot(aes(x = err, y = scan_group, group = scan_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "SCAN index") +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "scan_error.pdf"
  ),
  width = 6,
  height = 4
)

# spatial features impact ------

spatial_feat <- read_parquet("data/spatial_features.gpkg")

spatial_feat %>%
  ggplot(aes(x = dist_coast)) +
  geom_histogram(fill = blues9[7])

ggsave(
  file.path(
    fig_path,
    "dist_coast_hist.pdf"
  ),
  width = 6,
  height = 4
)

spatial_feat %>%
  ggplot(aes(x = elevation)) +
  geom_histogram(fill = blues9[7])

ggsave(
  file.path(
    fig_path,
    "elevation_hist.pdf"
  ),
  width = 6,
  height = 4
)

spatial_feat %>%
  ggplot(aes(x = elevation)) +
  geom_histogram(fill = blues9[7]) +
  facet_wrap(~tech_typ)

ggsave(
  file.path(
    fig_path,
    "elevation_hist_by_tech.pdf"
  ),
  width = 6,
  height = 4
)

fig_df <- pwr_curv_df %>%
  mutate(
    dist_coast_g5 = cut(
      dist_coast,
      breaks = quantile(dist_coast, probs = seq(0, 1, 0.2)),
      include.lowest = TRUE
    ),
    elevation_g5 = cut(
      elevation,
      breaks = quantile(elevation, probs = seq(0, 1, 0.2)),
      include.lowest = TRUE
    ),
    norm_potential = pmin(1, potential / capacity),
    norm_power_est0 = power_est0 / capacity,
    error0 = norm_potential - norm_power_est0
  )

## error by distance to coast -----
fig_df %>%
  ggplot() +
  geom_density_ridges(
    aes(x = error0, y = dist_coast_g5, fill = dist_coast_g5),
    alpha = 0.7
  ) +
  scale_fill_manual(
    values = blues9[seq(1, 9, length.out = 5)],
    name = "Distance to Coast (km)"
  ) +
  labs(x = "Normalised Error", y = "Distance to Coast Quintiles") +
  theme_ridges() +
  theme(legend.position = "none")

ggsave(
  file.path(
    fig_path,
    "error_by_dist_coast_rdens.pdf"
  ),
  width = 6,
  height = 4
)

fig_df %>%
  ggplot() +
  geom_density(
    aes(x = error0, fill = dist_coast_g5),
    alpha = 0.5
  ) +
  scale_fill_manual(
    values = blues9[seq(1, 9, length.out = 5)],
    name = "Distance to Coast (km)"
  ) +
  labs(x = "Normalised Error", y = "density") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(
  file.path(
    fig_path,
    "error_by_dist_coast_dens.pdf"
  ),
  width = 6,
  height = 4
)


## error by distance to coast and tech type -----
fig_df %>%
  ggplot() +
  geom_density_ridges(
    aes(x = error0, y = dist_coast_g5, fill = dist_coast_g5),
    alpha = 0.7
  ) +
  facet_wrap(~tech_typ) +
  scale_fill_manual(
    values = blues9[seq(1, 9, length.out = 5)],
    name = "Distance to Coast (km)"
  ) +
  labs(x = "Normalised Error", y = "Distance to Coast Quintiles") +
  theme_ridges() +
  theme(legend.position = "none")

ggsave(
  file.path(
    fig_path,
    "error_by_dist_coast_tech_rdens.pdf"
  ),
  width = 6,
  height = 4
)

fig_df %>%
  ggplot() +
  geom_density(
    aes(x = error0, fill = dist_coast_g5),
    alpha = 0.5
  ) +
  facet_wrap(~tech_typ) +
  scale_fill_manual(
    values = blues9[seq(1, 9, length.out = 5)],
    name = "Distance to Coast (km)"
  ) +
  labs(x = "Normalised Error", y = "density") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(
  file.path(
    fig_path,
    "error_by_dist_coast_tech_dens.pdf"
  ),
  width = 6,
  height = 4
)

## error by elevation -----
fig_df %>%
  ggplot() +
  geom_density_ridges(
    aes(x = error0, y = elevation_g5, fill = elevation_g5),
    alpha = 0.7
  ) +
  scale_fill_manual(
    values = blues9[seq(1, 9, length.out = 5)],
    name = "Elevation (m)"
  ) +
  labs(x = "Normalised Error", y = "Elevation Quintiles") +
  theme_ridges() +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "error_by_elevation_rdens.pdf"
  ),
  width = 6,
  height = 4
)

fig_df %>%
  ggplot() +
  geom_density(
    aes(x = error0, fill = elevation_g5),
    alpha = 0.5
  ) +
  scale_fill_manual(
    values = blues9[seq(1, 9, length.out = 5)],
    name = "Elevation (m)"
  ) +
  labs(x = "Normalised Error", y = "density") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(
  file.path(
    fig_path,
    "error_by_elevation_dens.pdf"
  ),
  width = 6,
  height = 4
)

# pwr_curv_df %>%
#   ggplot(aes(x = dist_coast, y = ws_h)) +
#   geom_hex()
