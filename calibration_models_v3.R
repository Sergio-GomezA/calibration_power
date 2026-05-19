#  Calibration models revised

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

source("aux_funct.R")


# 1.1 Load data ----------------------------------------------------------------
source("read_data.R")


# 1.2 EDA ------------------------------------------------------------------------

source("eda_figures.R")

# 1.3 1 day data frame -------------------------------------------------
pwr_curv_df <- read_parquet(file.path(gen_path, "power_curve_all.parquet"))

d0 <- as.Date("2024-08-10")
d0 <- as.Date("2025-08-10")

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
  filter(date >= d0, date <= d0 + n.days - 1) %>%
  group_by(lon, lat, time) %>%
  summarise(
    site_name = first(site_name),
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
    site_id = as.integer(factor(site_name)),
    ws_group = inla.group(ws_h, n = 20, method = "quantile"),
    pow_group = inla.group(norm_power_est0, n = 20, method = "quantile"),
    time_id = as.integer(factor(time)),
    # loc = cbind(x, y)
  )
wf_df_frag <- wf_df_frag %>%
  st_geometry() %>%
  (\(g) g / 1000)() %>%
  st_set_geometry(wf_df_frag, .)

# 2. Aggregated model ------------------------------------------------------------------
# Updates:
# - Adding entire dataset
# - Adding more covariates NAO, AO, EA, SCAN

## 2.1 LM step AIC selection --------------------------------------------------------------

base_model <- lm(
  norm_potential ~ norm_power_est0,
  data = wf_df_frag
)

# full_model <- lm(
#   norm_potential ~ tech_typ *
#     norm_power_est0 +
#     # norm_power_est0 * month +
#     # hour +
#     poly(ws_h_wmean, 3) +
#     nao * tech_typ +
#     ao * tech_typ +
#     ea * tech_typ +
#     scan * tech_typ,
#   data = wf_df_frag
# )

full_model0 <- lm(
  norm_potential ~ tech_typ *
    norm_power_est0 +
    # norm_power_est0 * month +
    # hour +
    poly(ws_h_wmean, 3),
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
model_AIC %>%
  saveRDS(
    file.path(
      model_path,
      sprintf("lm0_cir%s.rds", d0_tag)
    )
  )

# ModelMetrics::rmse(GB_df$norm_potential, model_AIC$fitted.values)
ModelMetrics::rmse(wf_df_frag$norm_potential, model_AIC0$fitted.values)

# 3. Spatially disaggregated model -----------------------------------------------------
# Updates:
# -- Add more locations (previously only 25)
# -- Add anomaly detection to handle outliers separately
# -- Add more covariates terrain, NAO, AO, EA, SCAN
# -- Add non-stationary covariance

## 3.1 Single run attempt --------------------------------------------------------------

## time series figures ####
source("ts_figures.R")

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
  )

bru0 <- bru(
  components = components0,
  formula = norm_potential ~ Intercept +
    # tech_typ +
    power_correction +
    wind +
    st_field,
  family = "gaussian",
  data = wf_df_frag
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
bru0$summary.fixed[, 1:6]
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
  n.samples = 20
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
  gg(wf.mesh, alpha = 0.5) +
  geom_point(data = loc_unique, aes(x, y), color = "darkred", size = 0.5) +
  facet_wrap(~time_id) +
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
  geom_sf(data = uk_map, fill = NA, color = "black", alpha = 0.5) +
  gg(wf.mesh, alpha = 0.5) +
  geom_point(data = loc_unique, aes(x, y), color = "darkred", size = 0.5) +
  facet_wrap(~time_id) +
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

# model comparison ####

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
  "QM"
)
est_cols <- c(
  "norm_power_est0",
  "lm",
  "st",
  "qm"
)
n <- nrow(wf_df_frag)
names(mod_labels) <- est_cols
model_df <- wf_df_frag %>%
  mutate(
    lm = model_AIC0$fitted.values,
    st = bru0$summary.fitted.values[1:n, "mean"],
    # ar = full_model_ar1$fitted,
    # lm_bru = bru0$summary.fitted.values[1:n, "mean"],
    # ar_bru = bru_ar$summary.fitted.values[1:n, "mean"],
    # ar_bru2 = bru_ar$summary.fitted.values[1:n, "mean"] -
    # bru_ar$summary.random$u[1:n, "mean"],
    qm = wgen_qm
  )

# write_parquet(model_df, sprintf("data/calibration_df_%s.parquet", d0_tag))
st_write(
  model_df,
  sprintf("data/calibration_df_%s.gpkg", d0_tag),
  driver = "GPKG",
)

## Reading fitted values ####

model_df <- st_read(sprintf("data/calibration_df_%s.gpkg", d0_tag))

# model_df %>% head()

df_long <- model_df %>%
  dplyr::select(tech_typ, norm_potential, all_of(est_cols)) %>%
  pivot_longer(
    cols = all_of(est_cols),
    names_to = "model",
    values_to = "estimate"
  )
require(ModelMetrics)
metrics_table <- df_long %>%
  st_drop_geometry() %>%
  group_by(model) %>%
  summarise(
    RMSE = rmse(actual = norm_potential, predicted = estimate),
    MAE = mae(actual = norm_potential, predicted = estimate),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(model = mod_labels[model]) %>%
  arrange(desc(RMSE))

write.csv(
  metrics_table,
  sprintf("summaries/calib_metrics_%s.csv", d0_tag),
  row.names = FALSE
)

df_long %>%
  st_drop_geometry() %>%
  group_by(tech_typ, model) %>%
  summarise(
    RMSE = rmse(actual = norm_potential, predicted = estimate),
    MAE = mae(actual = norm_potential, predicted = estimate),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(model = mod_labels[model]) %>%
  arrange(tech_typ, desc(RMSE))

## ts plots #####

# wf_df_frag %>% pull(site_name) %>% unique() %>% length()

## 3.2 Partitioned run attempt ---------------------------------------------------------

# 4. Spatially disaggregated model with non-stationary parameters ----------------------
