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

ggsave(
  "fig/daily_generation_by_regime.pdf",
  width = 6,
  height = 4
)

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
# modelsummary::modelsummary(
#   list(base_model, model_AIC0, model_AIC0_agg),
#   output = "latex",
#   stars = TRUE,
#   # model_names = "LM",
#   escape = FALSE,
#   shape = term ~ model + statistic,
#   # align = "l",
#   fmt = 3
# )
saveRDS(
  model_AIC0,
  file.path(model_path, sprintf("lm_model_aic0_%s.rds", d0_tag))
)
saveRDS(
  model_AIC0_agg,
  file.path(model_path, sprintf("lm_model_aic0_agg_%s.rds", d0_tag))
)
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

### 2.1.2.1 Model fitting ####

#### mesh building #####
cat("Building spatial mesh\n")

loc_unique <- wf_df_frag %>%
  distinct(x, y) %>%
  as.matrix()

# bnd <- fm_extensions(loc_unique, convex = c(-.1, -.15))
# bnd <- fm_extensions(loc_unique, convex = c(-.08, -.3))
bnd <- fm_extensions(loc_unique, convex = c(-.1, -.35))
# bnd <- fm_extensions(loc_unique, convex = c(-.1, -.15))
# ggplot() + geom_sf(data = bnd[[2]])
bndin <- bnd[[1]]
bndout <- bnd[[2]]
# uk_map <- rnaturalearth::ne_countries(
#   scale = "medium",
#   country = "United Kingdom",
#   returnclass = "sf"
# )
uk_map <- rnaturalearth::ne_countries(
  scale = "medium",
  country = "United Kingdom",
  returnclass = "sf"
)
uk_map <- uk_map %>%
  st_transform(crs = 27700) %>%
  st_geometry() %>%
  (\(g) g / 1000)() %>%
  st_set_geometry(uk_map, .)


hex_0 <- fm_hexagon_lattice(bnd[[1]], edge_len = 60)
# ggplot() +
#   geom_sf(data = uk_map, fill = NA, color = "black") +
#   geom_sf(data = hex_0) +
#   geom_point(data = loc_unique, aes(x, y), color = "darkred") +
#   theme_void()

# wf.mesh <- fm_mesh_2d(
#   # loc = loc_unique,
#   loc = fm_hexagon_lattice(bnd[[1]], edge_len = 30),
#   boundary = bnd,
#   max.edge = c(60, 120), # km
#   # offset = -0.2,
#   cutoff = 25
# )

wf.mesh <- fm_mesh_2d(
  # loc = loc_unique,
  loc = hex_0,
  boundary = bnd,
  max.edge = c(100, 150), # km
  min.angle = 30,
  # offset = -0.2,
  cutoff = 30,
  max.n.strict = c(400, 100)
)
saveRDS(
  wf.mesh,
  file.path(model_path, sprintf("spatial_mesh_coarse_%s.rds", d0_tag))
)
ggplot() +
  geom_sf(data = uk_map, fill = NA, color = "black") +
  gg(wf.mesh) +
  geom_point(data = loc_unique, aes(x, y), color = "darkred") +
  annotation_scale(location = "bl", width_hint = 0.25, plot_unit = "km") +
  theme_void()
ggsave(
  sprintf("fig/spatial_mesh_coarse_%s.pdf", d0_tag),
  width = 4,
  height = 6
)

#### mesh assessment #####
mesh_assessment <- fm_assess(mesh = wf.mesh, spatial.range = 60) %>%
  st_filter(., bndin)

ggplot() +
  geom_sf(data = mesh_assessment, aes(col = edge.len)) +
  geom_point(data = loc_unique, aes(x, y), color = "darkred") +
  geom_sf(data = uk_map, fill = NA, color = "white") +
  annotation_scale(location = "bl", width_hint = 0.25, plot_unit = "km") +
  theme_void() +
  scale_color_viridis_c(option = "D")
ggsave(
  sprintf("fig/spatial_mesh_coarse_assessment_edgelen_%s.pdf", d0_tag),
  width = 4,
  height = 6
)
# sd.dev should be close to 1
ggplot() +
  geom_sf(data = mesh_assessment, aes(col = sd.dev)) +
  geom_point(data = loc_unique, aes(x, y), color = "darkred") +
  geom_sf(data = uk_map, fill = NA, color = "white") +
  annotation_scale(location = "bl", width_hint = 0.25, plot_unit = "km") +
  theme_void() +
  scale_color_viridis_c(option = "D")
ggsave(
  sprintf("fig/spatial_mesh_coarse_assessment_sddev_%s.pdf", d0_tag),
  width = 4,
  height = 6
)

#### SPDE model ####
spde <- INLA::inla.spde2.pcmatern(
  mesh = wf.mesh,
  prior.range = c(50, 0.5), # P(range < 100 km)=0.5
  prior.sigma = c(0.2, 0.5) # P(sd > 0.2)=0.5
)

# components0 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
#   power_correction(norm_power_est0, model = "linear") + # fixed slope
#   tech_typ(tech_typ, model = "iid") + # random intercept by tech_typ
#   tech_power(tech_typ, model = "iid", weights = norm_power_est0) + # random slope
#   wind(ws_group, model = "rw2") +
#   st_field(
#     geometry,
#     model = spde,
#     group = time_id,
#     control.group = list(model = "ar1")
#   )
components0 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  # tech_typ(tech_typ, model = "iid") + # random intercept by tech_typ
  power_correction(
    pow_group,
    model = "rw2",
    replicate = tech_typ,
    constr = TRUE
  ) + # smooth correction power
  wind(ws_group, model = "rw2", constr = TRUE) + # smooth correction wind
  st_field(
    geometry,
    model = spde,
    group = time_id,
    control.group = list(model = "ar1")
  ) +
  dist_coast +
  elevation

bru0 <- bru(
  components = components0,
  formula = norm_potential ~ Intercept +
    # tech_typ +
    power_correction +
    wind +
    st_field +
    dist_coast +
    elevation,
  family = "gaussian",
  data = wf_df_frag
)

# saveRDS(
#   bru0,
#   file = file.path(model_path, sprintf("st_bru0_coarse_mesh_%s.rds", d0_tag))
# )
bru0 <- readRDS(file.path(
  model_path,
  sprintf("st_bru0_coarse_mesh_%s.rds", d0_tag)
))

### summary and effect plots ####
summary(bru0)

### 2.1.2.2 Model effects #####

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


## Looping through days #####

for (d0 in sampled_days) {
  cat(sprintf("Running models for day %s", format(d0, "%Y-%m-%d")))
}

for (d in sampled_days) {
  message(sprintf(
    "Running models for day %s",
    base::format(as.Date(d), "%Y-%m-%d")
  ))
  d0 <- d
  source("1_model/inlabru_1d.R")
}
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
n_models <- length(est_cols)
n <- nrow(wf_df_frag)
names(mod_labels) <- est_cols


model_df0 <- wf_df_frag %>%
  mutate(
    date = as.Date(time),
    lm = predict(model_AIC0, newdata = .),
    qm = wgen_qm,
    st = bru0$summary.fitted.values[1:n, "mean"],
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
    elevation,
    dist_coast,
    capacity,
    tech_typ,
    p_group3,
    norm_potential,
    any_of(est_cols)
  ) %>%
  mutate(hour = hour(time)) %>%
  mutate(
    p_group3 = factor(p_group3, levels = c("low", "mid", "high")),
    dist_coast_g4 = cut(
      dist_coast,
      breaks = quantile(dist_coast, probs = seq(0, 1, 0.25)),
      include.lowest = TRUE
    ),
    elevation_g4 = cut(
      elevation,
      breaks = quantile(elevation, probs = seq(0, 1, 0.25)),
      include.lowest = TRUE
    )
  ) %>%
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
  coord_cartesian(xlim = c(-0.25, 0.25)) +
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
    MAPE = mape(actual = norm_potential, predicted = estimate, pos_only = TRUE),
    MDAPE = mdape(
      actual = norm_potential,
      predicted = estimate,
      pos_only = TRUE
    ),
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

# aggregated time series version
metrics_table <- df_long0 %>%
  st_drop_geometry() %>%
  group_by(time, model) %>%
  summarise(
    across(
      c(norm_potential, estimate),
      ~ sum(. * capacity) / sum(capacity)
    )
  ) %>%
  group_by(model) %>%
  summarise(
    RMSE = rmse(actual = norm_potential, predicted = estimate),
    MAE = mae(actual = norm_potential, predicted = estimate),
    MAPE = mape(actual = norm_potential, predicted = estimate, pos_only = TRUE),
    MDAPE = mdape(
      actual = norm_potential,
      predicted = estimate,
      pos_only = TRUE
    ),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    model = mod_labels[model]
  ) %>%
  arrange(desc(RMSE))
metrics_table
write.csv(
  metrics_table,
  sprintf("summaries/calib_metrics_%s_gb.csv", d0_tag),
  row.names = FALSE
)

df_long0 %>%
  group_by(time, model) %>%
  summarise(
    across(
      c(norm_potential, estimate),
      ~ sum(. * capacity) / sum(capacity)
    )
  ) %>%
  mutate(
    err = estimate - norm_potential,
    # model = factor(model, levels = est_cols, labels = mod_labels)
  ) %>%
  group_by(model) %>%
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
  sprintf("fig/error_distribution_by_model_%s_gb.pdf", d0_tag),
  width = 6,
  height = 4
)


# by tech type
metrics_table <- df_long0 %>%
  st_drop_geometry() %>%
  group_by(tech_typ, model) %>%
  summarise(
    RMSE = rmse(actual = norm_potential, predicted = estimate),
    MAE = mae(actual = norm_potential, predicted = estimate),
    MAPE = mape(actual = norm_potential, predicted = estimate, pos_only = TRUE),
    MDAPE = mdape(
      actual = norm_potential,
      predicted = estimate,
      pos_only = TRUE
    ),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # mutate(model = mod_labels[model]) %>%
  arrange(tech_typ, desc(RMSE))
metrics_table
write.csv(
  metrics_table,
  sprintf("summaries/calib_metrics_%s_tech.csv", d0_tag),
  row.names = FALSE
)

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
  scale_fill_manual(values = pal_lancet()(n_models))

ggsave(
  sprintf("fig/error_distribution_by_tech_type_%s.pdf", d0_tag),
  width = 6,
  height = 4
)

# by regime
metrics_table <- df_long0 %>%
  st_drop_geometry() %>%
  group_by(p_group3, model) %>%
  summarise(
    RMSE = rmse(actual = norm_potential, predicted = estimate),
    MAE = mae(actual = norm_potential, predicted = estimate),
    MAPE = mape(actual = norm_potential, predicted = estimate, pos_only = TRUE),
    MDAPE = mdape(
      actual = norm_potential,
      predicted = estimate,
      pos_only = TRUE
    ),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(p_group3, desc(RMSE))
metrics_table

write.csv(
  metrics_table,
  sprintf("summaries/calib_metrics_%s_regime.csv", d0_tag),
  row.names = FALSE
)

df_long0 %>%
  st_drop_geometry() %>%
  ggplot() +
  geom_density_ridges(
    aes(err, p_group3, fill = model),
    alpha = 0.5,
    scale = 1
  ) +
  theme_ridges() +
  labs(
    title = "Error distribution by regime",
    x = "Error (model estimate - observed)",
    y = "Regime",
    fill = "Model"
  ) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = pal_lancet()(n_models))
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
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(-0.5, 0.5)) +
  scale_fill_manual(values = pal_lancet()(n_models))

ggsave(
  sprintf("fig/error_distribution_by_hourA_%s.pdf", d0_tag),
  width = 12,
  height = 8
)

df_long0 %>%
  st_drop_geometry() %>%
  group_by(hour, model) %>%
  summarise(
    RMSE = rmse(actual = norm_potential, predicted = estimate),
    MAE = mae(actual = norm_potential, predicted = estimate),
    MAPE = mape(actual = norm_potential, predicted = estimate, pos_only = TRUE),
    MDAPE = mdape(
      actual = norm_potential,
      predicted = estimate,
      pos_only = TRUE
    ),
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
  scale_color_manual(values = pal_lancet()(n_models))

ggsave(
  sprintf("fig/model_performance_by_hourS_%s.pdf", d0_tag),
  width = 6,
  height = 4
)


# by distance to coast
metrics_table <- df_long0 %>%
  st_drop_geometry() %>%
  group_by(dist_coast_g4, model) %>%
  summarise(
    RMSE = rmse(actual = norm_potential, predicted = estimate),
    MAE = mae(actual = norm_potential, predicted = estimate),
    MAPE = mape(actual = norm_potential, predicted = estimate, pos_only = TRUE),
    MDAPE = mdape(
      actual = norm_potential,
      predicted = estimate,
      pos_only = TRUE
    ),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(dist_coast_g4, desc(RMSE))
metrics_table

write.csv(
  metrics_table,
  sprintf("summaries/calib_metrics_%s_dist_coast.csv", d0_tag),
  row.names = FALSE
)

df_long0 %>%
  st_drop_geometry() %>%
  ggplot() +
  geom_density_ridges(
    aes(err, dist_coast_g4, fill = model),
    alpha = 0.5,
    scale = 1
  ) +
  theme_ridges() +
  labs(
    title = "Error distribution by distance to coast",
    x = "Error (model estimate - observed)",
    y = "Distance to Coast",
    fill = "Model"
  ) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = pal_lancet()(n_models))
ggsave(
  sprintf("fig/error_distribution_by_dist_coast_%s.pdf", d0_tag),
  width = 6,
  height = 4
)


# by elevation
metrics_table <- df_long0 %>%
  st_drop_geometry() %>%
  group_by(elevation_g4, model) %>%
  summarise(
    RMSE = rmse(actual = norm_potential, predicted = estimate),
    MAE = mae(actual = norm_potential, predicted = estimate),
    MAPE = mape(actual = norm_potential, predicted = estimate, pos_only = TRUE),
    MDAPE = mdape(
      actual = norm_potential,
      predicted = estimate,
      pos_only = TRUE
    ),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(elevation_g4, desc(RMSE))
metrics_table

write.csv(
  metrics_table,
  sprintf("summaries/calib_metrics_%s_elevation.csv", d0_tag),
  row.names = FALSE
)

df_long0 %>%
  st_drop_geometry() %>%
  ggplot() +
  geom_density_ridges(
    aes(err, elevation_g4, fill = model),
    alpha = 0.5,
    scale = 1
  ) +
  theme_ridges() +
  labs(
    title = "Error distribution by elevation",
    x = "Error (model estimate - observed)",
    y = "Elevation",
    fill = "Model"
  ) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = pal_lancet()(n_models))
ggsave(
  sprintf("fig/error_distribution_by_elevation_%s.pdf", d0_tag),
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

# 2. Aggregated model ------------------------------------------------------------------
# Updates:
# - Adding entire dataset
# - Adding more covariates NAO, AO, EA, SCAN

## 2.1 LM step AIC selection --------------------------------------------------------------

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


# model_AIC <- step(
#   base_model,
#   scope = list(lower = base_model, upper = full_model0),
#   # steps = 5,
#   k = 2
# )
model_AIC0 <- step(
  base_model,
  scope = list(lower = base_model, upper = full_model0),
  # steps = 5,
  k = 2
)
model_AIC0 %>%
  saveRDS(
    file.path(
      model_path,
      sprintf("lm0_cir%s.rds", d0_tag)
    )
  )

# ModelMetrics::rmse(GB_df$norm_potential, model_AIC$fitted.values)
# ModelMetrics::rmse(wf_df_frag$norm_potential, model_AIC0$fitted.values)

# 3. Spatially disaggregated model -----------------------------------------------------
# Updates:
# -- Add more locations (previously only 25)
# -- Add anomaly detection to handle outliers separately
# -- Add more covariates terrain, NAO, AO, EA, SCAN
# -- Add non-stationary covariance

## 3.1 Single run attempt --------------------------------------------------------------

## time series figures ####
source("ts_figures.R")
