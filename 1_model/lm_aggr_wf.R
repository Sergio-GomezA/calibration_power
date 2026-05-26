# 1. Load libraries and data ------------------------------------------------

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE


require(parallel)

if (local_run) {
  data_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
} else {
  data_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  gen_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  model_path <- "/exports/eddie/scratch/s2441782/calibration/model_objects"
  temp_lib <- "/exports/eddie3_homes_local/s2441782/lib"
  .libPaths(temp_lib)
}


# install.packages(
#   c("qmap"),
#   temp_lib,
#   dependencies = TRUE
# )

require(arrow)
require(dplyr)
require(tidyr)
# require(rnaturalearth)
# require(rnaturalearthdata)
require(sf)
require(ggplot2)
require(ggthemes)
require(ggsci)
# require(FNN)
require(data.table)
# require(parallel)
library(purrr)
# require(brms)
require(INLA)
require(inlabru)
require(fmesher)
require(qmap)
require(ggridges)
require(ggspatial)
require(elevatr)

require(ModelMetrics)


source("aux_funct.R")


# 1.1 Load data ----------------------------------------------------------------
source("read_data.R")

pwr_curv_df <- read_parquet(file.path(
  gen_path,
  "power_curve_all_enriched.parquet"
))

aggr_cat <- read_parquet(file.path(
  "data",
  "aggregated_catalogue.parquet"
))

spatial_feat <- st_read("data/spatial_features.gpkg")

# 1.2 choose days to run ----------------------------------------------------------------
## by regime

gb_day_df <- GB_df %>%
  group_by(date, tech_typ) %>%
  summarise(
    across(
      c(norm_power_est0, norm_potential),
      ~ sum(. * capacity) / sum(capacity)
    ),
    across(c(ws_h_wmean), ~ sum(. * capacity) / sum(capacity)),
    across(c(capacity), mean)
  ) %>%
  summarise(
    across(
      c(norm_power_est0, norm_potential),
      ~ sum(. * capacity) / sum(capacity)
    ),
    across(c(ws_h_wmean), ~ sum(. * capacity) / sum(capacity)),
    across(c(capacity), sum),
    .groups = "drop"
  )

cutprobs3 <- c(0.25, 0.75)
p_quant3 <- quantile(gb_day_df$norm_potential, probs = cutprobs3)
cutprobs7 <- c(0.1, 0.2, 0.25, 0.75, 0.8, 0.9)
p_quant7 <- quantile(gb_day_df$norm_potential, probs = cutprobs7)

gb_day_df <- gb_day_df %>%
  mutate(
    p_group3 = cut(
      norm_potential,
      breaks = c(-Inf, p_quant3, Inf),
      labels = c("low", "mid", "high")
    ),
    p_group7 = cut(
      norm_potential,
      breaks = c(-Inf, p_quant7, Inf)
    )
  )

gb_day_df %>%
  ggplot() +
  geom_histogram(aes(norm_potential), bins = 30) +
  geom_vline(xintercept = p_quant3, col = "red") +
  geom_vline(xintercept = p_quant7, col = "blue") +
  theme_minimal() +
  labs(
    title = "Distribution of daily generation",
    x = "Wind generation (% of capacity)",
    y = "Count"
  ) +
  theme(legend.position = "bottom")

gb_day_df %>%
  ggplot() +
  geom_line(aes(date, norm_potential), col = "blue") +
  theme_minimal() +
  labs(
    title = "Daily generation over time",
    x = "Date",
    y = "Wind generation (% of capacity)"
  ) +
  theme(legend.position = "bottom")


regime_palette <- pal_lancet()(3)[c(2, 3, 1)]
gb_day_df %>%
  ggplot(aes(norm_potential, p_group3, fill = p_group3)) +
  geom_density_ridges(alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Daily generation by regime",
    x = "Regime",
    y = "Wind generation (% of capacity)"
  ) +
  theme(legend.position = "bottom") +
  labs(fill = "Regime") +
  scale_fill_manual(values = regime_palette)

set.seed(1)
# sample 1 day per regime
sampled_days <- gb_day_df %>%
  group_by(p_group3) %>%
  slice_sample(n = 1) %>%
  pull(date)

sample_df <- GB_df %>%
  filter(date %in% sampled_days) %>%
  group_by(time) %>%
  summarise(
    across(
      c(norm_power_est0, norm_potential),
      ~ sum(. * capacity) / sum(capacity)
    ),
    across(c(ws_h_wmean), ~ sum(. * capacity) / sum(capacity)),
    across(c(capacity), sum),
    date = first(date),
    .groups = "drop"
  ) %>%
  left_join(gb_day_df %>% dplyr::select(date, p_group3), by = "date") %>%
  mutate(
    hour = format(time, "%H:%M"),
    legend_label = factor(
      paste(p_group3, format(date, "%y-%m-%d")),
      levels = paste(
        c("low", "mid", "high"),
        format(sampled_days, "%y-%m-%d")
      )
    )
  )

ggplot() +
  geom_line(
    data = sample_df,
    aes(
      hour,
      norm_potential,
      col = legend_label,
      group = legend_label
    )
  ) +
  labs(
    title = "Daily generation time series for sampled days",
    x = "Hour",
    y = "Wind generation (% of capacity)",
    col = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 1)
  ) +
  scale_color_manual(values = regime_palette)

#
ggsave(
  "fig/daily_generation_time_series_sampled_days.pdf",
  width = 6,
  height = 4
)


# 2.1 Models for WF ####

d0 <- sampled_days
d0_tag <- "samp_days"
# n.days <- 1
wf_df_frag <- pwr_curv_df %>%
  rename(time = halfHourEndTime) %>%
  mutate(
    date = as.Date(time),
    site_name = site_name %>%
      gsub("\\b(wind\\s*farm|wf)\\b", "", ., ignore.case = TRUE) %>%
      trimws()
  ) %>%
  filter(date %in% sampled_days) %>%
  arrange(site_name) %>%
  mutate(
    site_id = as.integer(factor(site_name)),
    coord_id = as.integer(factor(paste(lon, lat)))
  ) %>%
  group_by(lon, lat, time) %>%
  summarise(
    site_name = first(site_name),
    coord_id = first(coord_id),
    elevation = first(elevation),
    dist_coast = first(dist_coast),
    tech_typ = first(tech_typ),
    across(c(ws_h, wd10, wd100), mean),
    ws_h_wmean = sum(ws_h * capacity) / sum(capacity),
    across(c(potential, power_est0, capacity, curtailment), sum),
    .groups = "drop"
  ) %>%
  mutate(
    norm_potential = pmin(1, potential / capacity),
    norm_power_est0 = power_est0 / capacity,
    error0 = norm_potential - norm_power_est0
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  mutate(lon = st_coordinates(.)[, 1], lat = st_coordinates(.)[, 2]) %>%
  st_transform(crs = 27700) %>%
  mutate(
    x = st_coordinates(.)[, 1] / 1000,
    y = st_coordinates(.)[, 2] / 1000,
  ) %>%
  # st_drop_geometry() %>%
  mutate(
    # site_id = as.integer(factor(site_name)),
    ws_group = inla.group(ws_h, n = 20, method = "quantile"),
    pow_group = inla.group(norm_power_est0, n = 20, method = "quantile"),
    time_id = as.integer(factor(time)),
    # loc = cbind(x, y)
  )

wf_df_frag <- wf_df_frag %>%
  st_geometry() %>%
  (\(g) g / 1000)() %>%
  st_set_geometry(wf_df_frag, .)


## 2.1.1 Linear model ####

base_model <- lm(
  norm_potential ~ norm_power_est0,
  data = wf_df_frag
)

full_model0 <- lm(
  norm_potential ~
    tech_typ *
    norm_power_est0 +
    # norm_power_est0 * month +
    # hour +
    # dist_coast * tech_typ +
    # elevation * tech_typ +
    # dist_coast:tech_typ +
    # elevation:tech_typ +
    tech_typ * poly(dist_coast, 2) +
    tech_typ * poly(elevation, 3) +
    tech_typ * poly(ws_h_wmean, 3),
  data = wf_df_frag
)
summary(full_model0)

model_AIC0 <- step(
  base_model,
  scope = list(lower = base_model, upper = full_model0),
  # steps = 5,
  k = 2
)

summary(model_AIC0)

full_model1 <- lm(
  norm_potential ~
    tech_typ *
    norm_power_est0 +
    # norm_power_est0 * month +
    # hour +
    dist_coast * tech_typ +
    elevation * tech_typ +
    # dist_coast:tech_typ +
    # elevation:tech_typ +
    # tech_typ * poly(dist_coast, 2) +
    # tech_typ * poly(elevation, 3) +
    tech_typ * poly(ws_h_wmean, 3),
  data = wf_df_frag
)
summary(full_model1)

model_AIC1 <- step(
  base_model,
  scope = list(lower = base_model, upper = full_model1),
  # steps = 5,
  k = 2
)

summary(model_AIC1)

## 2.1.2 Inlabru ####

## 2.1.3 Quantile mapping ####
qqmod <- fitQmap(
  obs = wf_df_frag %>% pull(norm_potential),
  mod = wf_df_frag %>% pull(norm_power_est0),
  method = "QUANT"
)

wgen_qm <- with(
  wf_df_frag,
  doQmapQUANT(norm_power_est0, qqmod, type = "linear")
)

## 2.1.1 Model predictions ####

# 2.2 Model for aggregate ####

samp_gb <- GB_df
# %>%
# filter(date %in% sampled_days) %>%

base_model_agg <- lm(
  norm_potential ~ norm_power_est0,
  data = samp_gb
)

full_model0_agg <- lm(
  norm_potential ~
    tech_typ *
    norm_power_est0 +
    # norm_power_est0 * month +
    # hour +
    # dist_coast * tech_typ +
    # elevation * tech_typ +
    # dist_coast:tech_typ +
    # elevation:tech_typ +
    # tech_typ * poly(dist_coast, 2) +
    # tech_typ * poly(elevation, 3) +
    tech_typ * poly(ws_h_wmean, 3),
  data = samp_gb
)
summary(full_model0_agg)

model_AIC0_agg <- step(
  base_model_agg,
  scope = list(lower = base_model_agg, upper = full_model0_agg),
  # steps = 5,
  k = 2
)

summary(model_AIC0_agg)

## 2.2.1 Model predictions ####

mod_labels <- c(
  "Generic PC",
  "Linear model",
  "Spatio-temporal model",
  "QM",
  "GB LM"
)
est_cols <- c(
  "norm_power_est0",
  "lm",
  "st",
  "qm",
  "agg_lm"
)
n <- nrow(wf_df_frag)
names(mod_labels) <- est_cols


model_df0 <- wf_df_frag %>%
  mutate(
    date = as.Date(time),
    lm = predict(full_model1, newdata = .),
    qm = wgen_qm,
    agg_lm = predict(model_AIC0_agg, newdata = wf_df_frag)
  ) %>%
  left_join(
    gb_day_df %>% dplyr::select(date, p_group3),
    by = c("date" = "date")
  )

st_write(
  model_df0,
  sprintf("data/calibration_df_%s.gpkg", d0_tag),
  driver = "GPKG",
  append = FALSE,
)

## Reading fitted values ####

model_df0 <- st_read(sprintf("data/calibration_df_%s.gpkg", d0_tag))

# 3.1 Comparison

df_long0 <- model_df0 %>%
  dplyr::select(
    date,
    time,
    site_name,
    coord_id,
    tech_typ,
    p_group3,
    norm_potential,
    any_of(est_cols)
  ) %>%
  mutate(hour = hour(time)) %>%
  pivot_longer(
    cols = any_of(est_cols),
    names_to = "model",
    values_to = "estimate"
  ) %>%
  mutate(,
    err = estimate - norm_potential,
    p_group3 = forcats::fct_rev(p_group3),
    model = factor(model, levels = est_cols, labels = mod_labels)
  )

df_long0 %>%
  ggplot() +
  geom_density(
    aes(x = err, fill = model),
    alpha = 0.5
  ) +
  facet_wrap(~model) +
  theme_minimal() +
  labs(
    title = "Error distribution by model",
    x = "Error (model estimate - observed)",
    y = "Density",
    fill = "Model"
  ) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = pal_lancet()(5))

ggsave(
  sprintf("fig/error_distribution_by_model_%s.pdf", d0_tag),
  width = 6,
  height = 4
)

metrics_table <- df_long0 %>%
  st_drop_geometry() %>%
  group_by(model) %>%
  summarise(
    RMSE = rmse(actual = norm_potential, predicted = estimate),
    MAE = mae(actual = norm_potential, predicted = estimate),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    model = mod_labels[model]
  ) %>%
  arrange(desc(RMSE))

write.csv(
  metrics_table,
  sprintf("summaries/calib_metrics_%s.csv", d0_tag),
  row.names = FALSE
)

# by tech type
df_long0 %>%
  st_drop_geometry() %>%
  group_by(tech_typ, model) %>%
  summarise(
    RMSE = rmse(actual = norm_potential, predicted = estimate),
    MAE = mae(actual = norm_potential, predicted = estimate),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # mutate(model = mod_labels[model]) %>%
  arrange(tech_typ, desc(RMSE))

df_long0 %>%
  st_drop_geometry() %>%
  ggplot() +
  geom_density_ridges(
    aes(err, tech_typ, fill = model),
    alpha = 0.5,
    scale = 1
  ) +
  theme_ridges() +
  labs(
    title = "Error distribution by technology type",
    x = "Error (model estimate - observed)",
    y = "Technology Type",
    fill = "Model"
  ) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = pal_lancet()(4))

ggsave(
  sprintf("fig/error_distribution_by_tech_type_%s.pdf", d0_tag),
  width = 6,
  height = 4
)

# by regime
df_long0 %>%
  st_drop_geometry() %>%
  group_by(p_group3, model) %>%
  summarise(
    RMSE = rmse(actual = norm_potential, predicted = estimate),
    MAE = mae(actual = norm_potential, predicted = estimate),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(p_group3, desc(RMSE))


df_long0 %>%
  st_drop_geometry() %>%
  ggplot() +
  geom_density_ridges(aes(err, p_group3, fill = model), alpha = 0.5) +
  theme_ridges() +
  labs(
    title = "Error distribution by regime",
    x = "Error (model estimate - observed)",
    y = "Regime",
    fill = "Model"
  ) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = pal_lancet()(4))
ggsave(
  sprintf("fig/error_distribution_by_regime_%s.pdf", d0_tag),
  width = 6,
  height = 4
)


# by hour of day
df_long0 %>%
  mutate(hour = factor(hour, levels = 0:23)) %>%
  st_drop_geometry() %>%
  ggplot() +
  geom_density_ridges(
    aes(x = err, y = hour, fill = model),
    alpha = 0.5
  ) +
  facet_wrap(~model) +
  theme_ridges() +
  labs(
    title = "Error distribution by hour of day",
    x = "Error (model estimate - observed)",
    y = "Hour of day",
    fill = "Model"
  ) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = pal_lancet()(4))

ggsave(
  sprintf("fig/error_distribution_by_hourA_%s.pdf", d0_tag),
  width = 6,
  height = 4
)

df_long0 %>%
  st_drop_geometry() %>%
  group_by(hour, model) %>%
  summarise(
    RMSE = rmse(actual = norm_potential, predicted = estimate),
    MAE = mae(actual = norm_potential, predicted = estimate),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(model = mod_labels[model]) %>%
  arrange(hour, desc(RMSE)) %>%
  ggplot(aes(hour, RMSE, col = model, group = model)) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Model performance by hour of day",
    x = "Hour of day",
    y = "RMSE",
    col = "Model"
  ) +
  scale_color_manual(values = pal_lancet()(4))

ggsave(
  sprintf("fig/model_performance_by_hourS_%s.pdf", d0_tag),
  width = 6,
  height = 4
)

#how performance changes by onshore/offshore, season, hour, and regime and how well farm level predictions aggregate

# aggregator GB, Onshore, Offshore, Scotland, rest of GB

# validation
# time GB Onshore, Offshore,
# dense regions, sparse regions,
# low, mid, high generation regimes
# Scotland, rest of GB
