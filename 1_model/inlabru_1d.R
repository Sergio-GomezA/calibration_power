require(tidyverse)
require(sf)
require(inlabru)
require(fmesher)
require(ggspatial)
require(ModelMetrics)
require(qmap)
require(ggridges)

# inlabru model 1 day

# inlabru models by day ------

# sampled_days <- c("2020-08-14", "2024-04-17", "2024-04-12")
# d0 <- sampled_days[2] %>% as.Date()

d0_tag <- format(d0, "%y%m%d")
n.days <- 1
wf_df_frag <- pwr_curv_df %>%
  rename(time = halfHourEndTime) %>%
  mutate(
    date = as.Date(time),
    site_name = site_name %>%
      gsub("\\b(wind\\s*farm|wf)\\b", "", ., ignore.case = TRUE) %>%
      trimws()
  ) %>%
  # filter(date %in% sampled_days) %>%
  filter(date >= d0, date <= d0 + n.days - 1) %>%
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


## mesh building #####
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
coastline <- uk_map %>%
  st_transform(crs = 27700) %>%
  st_boundary()
uk_map <- uk_map %>%
  st_transform(crs = 27700) %>%
  st_geometry() %>%
  (\(g) g / 1000)() %>%
  st_set_geometry(uk_map, .)

edge_target <- 20 # km
hex_0 <- fm_hexagon_lattice(bnd[[1]], edge_len = edge_target)

wf.mesh <- fm_mesh_2d(
  # loc = loc_unique,
  loc = hex_0,
  boundary = bnd,
  max.edge = c(100, 150), # km
  min.angle = 25,
  # offset = -0.2,
  cutoff = edge_target,
  max.n.strict = c(900, 150)
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

### mesh assessment #####
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

## SPDE model ####
spde <- INLA::inla.spde2.pcmatern(
  mesh = wf.mesh,
  prior.range = c(50, 0.5), # P(range < 100 km)=0.5
  prior.sigma = c(0.2, 0.5) # P(sd > 0.2)=0.5
)


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
  )

bru0 <- bru(
  components = components0,
  formula = norm_potential ~ Intercept +
    power_correction +
    wind +
    st_field,
  family = "gaussian",
  data = wf_df_frag,
  options = bru_options(
    bru_verbose = 3,
    control.inla = list(verbose = TRUE)
  )
)

saveRDS(
  bru0,
  file = file.path(model_path, sprintf("st_bru0_coarse_mesh_%s.rds", d0_tag))
)
# bru0 <- readRDS(file.path(
#   model_path,
#   sprintf("st_bru0_coarse_mesh_%s.rds", d0_tag)
# ))

## summary and effect plots ####
summary(bru0)
# bru0$summary.fixed[, 1:6]
# bru0$summary.random$tech_typ[, 1:6]
# bru0$summary.random$tech_power[, 1:6]

# source("aux_funct.R")
plot.effects(bru0, "wind", show.plot = TRUE)
ggsave(sprintf("fig/wind_effect_coarse_%s.pdf", d0_tag), width = 6, height = 4)
plot.effects(
  bru0,
  "power_correction",
  show.plot = TRUE,
  n.replicate = 2,
  replicate_names = c("Offshore", "Onshore")
)
ggsave(
  sprintf("fig/power_correction_effect_coarse_%s.pdf", d0_tag),
  width = 6,
  height = 4
)
plot.hyper.dens(bru0)
ggsave(
  sprintf("fig/hyperparameters_coarse_%s.pdf", d0_tag),
  width = 6,
  height = 4
)
# wf_df_frag %>%
#   pull(pow_group) %>%
#   range()

# wf_df_frag %>%
#   ggplot()+ geom_density(aes(pow_group, fill = tech_typ), alpha = 0.5)+theme_minimal()

### plot intensity of spatial field ####

ppxl <- fm_pixels(wf.mesh, mask = bnd[[2]], format = "sf")
ppxl_all <- fm_cprod(
  ppxl,
  data.frame(
    # time_id = unique(wf_df_frag$time_id)
    time_id = c(9, 12, 18)
  )
)

set.seed(1)
pow_est_st <- predict(
  bru0,
  ppxl_all,
  ~ data.frame(
    time_id = time_id,
    norm_potential_est = st_field
  ),
  n.samples = 100
)

# ppxl <- fm_pixels(wf.mesh, mask = bnd[[2]], format = "sf", dims = c(50, 50))

# pow_est_st <- predict(
#   bru0,
#   ppxl,
#   ~ data.frame(
#     time_id = time_id,
#     norm_potential_est = st_field
#   ),
#   n.samples = 20
# )
# predict(bru0, ppxl, n.samples = 20)

p_median <- ggplot() +
  gg(pow_est_st, geom = "tile", aes(fill = q0.5)) +
  geom_sf(data = uk_map, fill = NA, color = "black", alpha = 0.5) +
  # gg(wf.mesh, alpha = 0.5) +
  geom_point(data = loc_unique, aes(x, y), color = "darkred", size = 0.5) +
  facet_wrap(
    ~time_id,
    labeller = as_labeller(c("9" = "9:00", "12" = "12:00", "18" = "18:00"))
  ) +
  coord_sf() +
  scale_fill_viridis_c() +
  theme_void()
p_median
ggsave(
  sprintf("fig/coarse_spatial_field_median_%s.pdf", d0_tag),
  width = 10,
  height = 6
)

p_sd <- ggplot() +
  gg(pow_est_st, geom = "tile", aes(fill = sd)) +
  geom_sf(data = uk_map, fill = NA, color = "white", alpha = 0.5) +
  # gg(wf.mesh, alpha = 0.5) +
  geom_point(data = loc_unique, aes(x, y), color = "darkred", size = 0.5) +
  facet_wrap(
    ~time_id,
    labeller = as_labeller(c("9" = "9:00", "12" = "12:00", "18" = "18:00"))
  ) +
  coord_sf() +
  scale_fill_viridis_c(option = "inferno") +
  theme_void()
p_sd
ggsave(
  sprintf("fig/coarse_spatial_field_sd_%s.pdf", d0_tag),
  width = 10,
  height = 6
)
# plot ts
# mesh triangle size
# different days
# covariate for terrain if necessary
#

# ppxl_all <- fm_cprod(
#   ppxl,
#   data.frame(
#     time_id = seq(9, 12, 18),
#     tech_typ = unique(wf_df_frag$tech_typ),
#     norm_power_est0 = seq(0, 1, length.out = 10),
#     ws_group = unique(wf_df_frag$ws_group) %>% sort()
#   )
# )

# pow_est_st <- predict(
#   bru0,
#   ppxl_all,
#   ~ data.frame(
#     time_id = time_id,
#     norm_potential_est = (Intercept +
#       power_correction +
#       tech_typ +
#       tech_power +
#       wind +
#       st_field)
#   )
# )

# lm wf version ####

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


model_AIC0 <- step(
  base_model,
  scope = list(lower = base_model, upper = full_model0),
  # steps = 5,
  k = 2
)

saveRDS(
  model_AIC0,
  file.path(model_path, sprintf("lm_model_aic0_%s.rds", d0_tag))
)


# GB lm version #####

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
saveRDS(
  model_AIC0_agg,
  file.path(model_path, sprintf("lm_model_aic0_agg_%s.rds", d0_tag))
)

# model comparison ####
model_AIC0 <- readRDS(file.path(
  model_path,
  sprintf("lm_model_aic0_%s.rds", d0_tag)
))
qqmod <- fitQmap(
  obs = wf_df_frag %>% pull(norm_potential),
  mod = wf_df_frag %>% pull(norm_power_est0),
  method = "QUANT"
)

wgen_qm <- with(
  wf_df_frag,
  doQmapQUANT(norm_power_est0, qqmod, type = "linear")
)

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
# length(wgen_qm)
# length(wf_df_frag$norm_potential)
# length(bru0$summary.fitted.values[1:n, "mean"])
# length(model_AIC0$fitted.values)
n <- nrow(wf_df_frag)
names(mod_labels) <- est_cols

model_df0 <- wf_df_frag %>%
  mutate(
    date = as.Date(time),
    lm = predict(model_AIC0, newdata = .),
    st = bru0$summary.fitted.values[1:n, "mean"],
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

## Exploring fitted values ####

model_df0 <- st_read(sprintf("data/calibration_df_%s.gpkg", d0_tag))


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


## ts plots #####

mod_labels2 <- c(mod_labels, "norm_potential" = "Observed")
model_df_ts <- model_df0 %>%
  dplyr::select(
    time,
    site_name,
    tech_typ,
    norm_potential,
    capacity,
    ws_h,
    all_of(est_cols)
  ) %>%
  pivot_longer(
    cols = all_of(c("norm_potential", est_cols)),
    names_to = "model",
    values_to = "estimate"
  )

model_df_ts %>%
  ggplot() +
  geom_line(
    aes(time, estimate, group = site_name),
    alpha = 0.5,
    col = "gray50"
  ) +
  geom_line(
    data = model_df_ts %>%
      group_by(time, model) %>%
      summarise(
        power = sum(estimate * capacity) / sum(capacity),
        .groups = "drop"
      ),
    aes(time, power, col = "capacity weighted avg."),
    lwd = 1
  ) +
  geom_line(
    data = model_df_ts %>%
      group_by(time, model) %>%
      summarise(power = mean(estimate), .groups = "drop"),
    aes(time, power, col = "simple avg."),
    lwd = 1
  ) +
  theme_minimal() +
  facet_wrap(~model, ncol = 3, labeller = as_labeller(mod_labels2)) +
  scale_x_datetime(date_labels = "%H:%M") +
  labs(
    title = sprintf("Power estimates Time Series %s", d0),
    x = "Time",
    y = "Generation (% of capacity)",
    col = ""
  ) +
  scale_color_manual(
    values = c("capacity weighted avg." = "blue", "simple avg." = "red")
  ) +
  theme(legend.position = "bottom")
ggsave(
  sprintf("fig/power_estimates_time_series_%s.pdf", d0_tag),
  width = 10,
  height = 6
)

model_df_ts2 <- model_df0 %>%
  dplyr::select(
    time,
    site_name,
    tech_typ,
    norm_potential,
    capacity,
    ws_h,
    all_of(est_cols)
  ) %>%
  mutate(across(
    all_of(est_cols),
    ~ . - norm_potential
  )) %>%
  pivot_longer(
    cols = all_of(est_cols),
    names_to = "model",
    values_to = "error"
  )
model_df_ts2 %>%
  ggplot() +
  geom_line(
    aes(time, error, group = site_name),
    alpha = 0.5,
    col = "gray50"
  ) +
  geom_line(
    data = model_df_ts2 %>%
      group_by(time, model) %>%
      summarise(
        error = sum(error * capacity) / sum(capacity),
        .groups = "drop"
      ),
    aes(time, error, col = "capacity weighted avg."),
    lwd = 1
  ) +
  geom_line(
    data = model_df_ts2 %>%
      group_by(time, model) %>%
      summarise(error = mean(error), .groups = "drop"),
    aes(time, error, col = "simple avg."),
    lwd = 1
  ) +
  theme_minimal() +
  facet_wrap(~model, ncol = 3, labeller = as_labeller(mod_labels)) +
  scale_x_datetime(date_labels = "%H:%M") +
  labs(
    title = sprintf("Power estimates Time Series %s", d0),
    x = "Time",
    y = "Error (% of capacity)",
    col = ""
  ) +
  scale_color_manual(
    values = c("capacity weighted avg." = "blue", "simple avg." = "red")
  ) +
  theme(legend.position = "bottom")
ggsave(
  sprintf("fig/power_estimates_error_time_series_%s.pdf", d0_tag),
  width = 10,
  height = 6
)
