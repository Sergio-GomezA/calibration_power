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

source("aux_funct.R")


# 1.1 Load data ----------------------------------------------------------------
source("read_data.R")

# 1.2 EDA ------------------------------------------------------------------------

source("eda_figures.R")

# 2. Aggregated model ------------------------------------------------------------------
# Updates:
# - Adding entire dataset
# - Adding more covariates NAO, AO, EA, SCAN

## 2.1 LM step AIC selection --------------------------------------------------------------

base_model <- lm(
  norm_potential ~ norm_power_est0,
  data = GB_df
)

full_model <- lm(
  norm_potential ~ tech_typ *
    norm_power_est0 +
    norm_power_est0 * month +
    hour +
    poly(ws_h_wmean, 3) +
    nao * tech_typ +
    ao * tech_typ +
    ea * tech_typ +
    scan * tech_typ,
  data = GB_df
)

full_model0 <- lm(
  norm_potential ~ tech_typ *
    norm_power_est0 +
    norm_power_est0 * month +
    hour +
    poly(ws_h_wmean, 3),
  data = GB_df
)
summary(full_model)


model_AIC <- step(
  base_model,
  scope = list(lower = base_model, upper = full_model),
  # steps = 5,
  k = 2
)
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
      "lm0_cir.rds"
    )
  )

ModelMetrics::rmse(GB_df$norm_potential, model_AIC$fitted.values)
ModelMetrics::rmse(GB_df$norm_potential, model_AIC0$fitted.values)

# 3. Spatially disaggregated model -----------------------------------------------------
# Updates:
# -- Add more locations (previously only 25)
# -- Add anomaly detection to handle outliers separately
# -- Add more covariates terrain, NAO, AO, EA, SCAN
# -- Add non-stationary covariance

## 3.1 Single run attempt --------------------------------------------------------------
pwr_curv_df <- read_parquet(file.path(gen_path, "power_curve.parquet"))

d0 <- as.Date("2024-08-10")
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
    time_id = as.integer(factor(time)),
    # loc = cbind(x, y)
  )
wf_df_frag <- wf_df_frag %>%
  st_geometry() %>%
  (\(g) g / 1000)() %>%
  st_set_geometry(wf_df_frag, .)


cat("Building spatial mesh\n")

loc_unique <- wf_df_frag %>%
  distinct(x, y) %>%
  as.matrix()
bnd <- fm_extensions(loc_unique, convex = c(-.1, -.15))
# bnd <- fm_extensions(loc_unique, convex = c(-.1, -.15))
# ggplot() + geom_sf(data = bnd[[2]])

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

wf.mesh <- fm_mesh_2d(
  # loc = loc_unique,
  loc = fm_hexagon_lattice(bnd[[1]], edge_len = 30),
  boundary = bnd,
  max.edge = c(60, 120), # km
  # offset = -0.2,
  cutoff = 25
)

fm_assess()
# close to 1

# plot(wf.mesh)
# points(
#   loc_unique,
#   col = 2,
#   pch = 16,
#   cex = 1
# )

ggplot() +
  geom_sf(data = uk_map, fill = NA, color = "black") +
  gg(wf.mesh) +
  geom_point(data = loc_unique, aes(x, y), color = "darkred") +
  theme_void()
wf.mesh$n

# SPDE model
spde <- INLA::inla.spde2.pcmatern(
  mesh = wf.mesh,
  prior.range = c(100, 0.5), # P(range < 100 km)=0.5
  prior.sigma = c(0.2, 0.5) # P(sd > 0.2)=0.5
)


components0 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  power_correction(norm_power_est0, model = "linear") + # fixed slope
  tech_typ(tech_typ, model = "iid") + # random intercept by tech_typ
  tech_power(tech_typ, model = "iid", weights = norm_power_est0) + # random slope
  wind(ws_group, model = "rw2") +
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
    tech_typ +
    tech_power +
    wind +
    st_field,
  family = "gaussian",
  data = wf_df_frag
)
summary(bru0)

# saveRDS(bru0, file = file.path(model_path, "st_bru0.rds"))
bru0 <- readRDS(file.path(model_path, "st_bru0.rds"))

bru0$summary.fixed[, 1:6]
bru0$summary.random$tech_typ[, 1:6]
bru0$summary.random$tech_power[, 1:6]
# bru0$summary.random$wind
plot.effects(bru0, "wind", show.plot = TRUE)

# plot intensity of spatial field

ppxl <- fm_pixels(wf.mesh, mask = bnd[[2]], format = "sf")
ppxl_all <- fm_cprod(
  ppxl,
  data.frame(
    time_id = c(9, 12, 18)
  )
)

pow_est_st <- predict(
  bru0,
  ppxl_all,
  ~ data.frame(
    time_id = time_id,
    norm_potential_est = st_field
  )
)

p_median <- ggplot() +
  gg(pow_est_st, geom = "tile", aes(fill = q0.5)) +
  geom_sf(data = uk_map, fill = NA, color = "black", alpha = 0.5) +
  geom_point(data = loc_unique, aes(x, y), color = "darkred") +
  facet_wrap(~time_id) +
  coord_sf() +
  scale_fill_viridis_c() +
  theme_void()
p_median


p_sd <- ggplot() +
  gg(pow_est_st, geom = "tile", aes(fill = sd)) +
  geom_sf(data = uk_map, fill = NA, color = "black", alpha = 0.5) +
  geom_point(data = loc_unique, aes(x, y), color = "darkred") +
  facet_wrap(~time_id) +
  coord_sf() +
  scale_fill_viridis_c(option = "inferno") +
  theme_void()
p_sd

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

# wf_df_frag %>% pull(site_name) %>% unique() %>% length()

## 3.2 Partitioned run attempt ---------------------------------------------------------

# 4. Spatially disaggregated model with non-stationary parameters ----------------------
