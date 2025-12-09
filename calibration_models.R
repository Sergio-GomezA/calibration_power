# Modelling calibration

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
library(purrr)
require(brms)
require(INLA)

require(inlabru)

source("aux_funct.R")
# local({
#   r <- c(
#     INLA = "https://inla.r-inla-download.org/R/testing",
#     CRAN = "https://cloud.r-project.org/",
#     inlabru_universe = "https://inlabru-org.r-universe.dev"
#   )
#   options(repos = r)
# })

# install.packages(c("INLA", "inlabru"), dependencies = TRUE)
# install.packages(
#   c("INLA", "inlabru"),
#   "/exports/eddie/scratch/s2441782/calibration/lib",
#   dependencies = TRUE
# )
# inla.binary.install()
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install(c("graph", "Rgraphviz"), dep = TRUE)

# read Data ####
data_path <- "~/Documents/ERA5_at_wf/"
gen_path <- "~/Documents/elexon/"
model_path <- "~/Documents/elexon/model_objects"
# era_df <- read_parquet(
#   file.path(data_path, "era5_combined.parquet")
# ) %>%
#   mutate(
#     ws100 = sqrt(u100^2 + v100^2),
#     ws10 = sqrt(u10^2 + v10^2),
#     wd100 = (atan2(-u100, -v100) * 180 / pi) %% 360,
#     wd10 = (atan2(-u10, -v10) * 180 / pi) %% 360
#   )

# # historical generation with curtailment and outages list
# gen_adj <- read_parquet(
#   file.path(gen_path, "gen_adj_v2.parquet")
# )
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

## aggregated to GB level ####
GB_df <- pwr_curv_df %>%
  filter(halfHourEndTime >= "2023-01-01", halfHourEndTime < "2025-01-01") %>%
  group_by(tech_typ, halfHourEndTime) %>%
  summarise(
    # ws_h_wmean = sum(ws_h * capacity),
    across(c(ws_h, wd10, wd100), ~ sum(. * capacity), .names = "{.col}_wmean"),
    across(
      c(power_est0, potential, capacity),
      sum
    ),
    across(c(ws_h, wd10, wd100), mean, .names = "{.col}_mean")
  ) %>%
  mutate(
    across(
      c(power_est0, potential),
      ~ . / capacity,
      .names = "norm_{.col}"
    ),
    # ws_h_wmean = ws_h_wmean / capacity,
    across(matches("_wmean"), ~ . / capacity),
    month = factor(month(halfHourEndTime)),
    hour = factor(hour(halfHourEndTime))
  ) %>%
  mutate(ws_group = inla.group(ws_h_wmean, n = 20, method = "quantile")) %>%
  group_by(tech_typ) %>%
  arrange(halfHourEndTime, .by_group = TRUE) %>%
  mutate(t = row_number()) %>%
  ungroup()

# write_parquet(GB_df, file.path(gen_path, "GB_aggr.parquet"))
GB_df <- read_parquet(file.path(gen_path, "GB_aggr.parquet"))

# write_parquet(GB_df, file.path(gen_path, "GB_aggr_frag.parquet"))
GB_df <- read_parquet(file.path(gen_path, "GB_aggr_frag.parquet"))

# Linear ####

base_model <- lm(
  potential ~ power_est0,
  data = GB_df
)

full_model <- lm(
  potential ~ tech_typ *
    power_est0 +
    power_est0 * month +
    hour +
    poly(ws_h_wmean, 3),
  data = GB_df
)
summary(full_model)

model_uBIC <- step(
  base_model,
  scope = list(lower = base_model, upper = full_model),
  # steps = 5,
  k = log(nrow(GB_df))
)

anova(model_uBIC)

model_uAIC <- step(
  base_model,
  scope = list(lower = base_model, upper = full_model),
  # steps = 5,
  k = 2
)

anova(model_uAIC)

ModelMetrics::rmse(model_uAIC$fitted.values, GB_df$potential)
ModelMetrics::rmse(model_uBIC$fitted.values, GB_df$potential)
GB_df$lin.fit <- model_uAIC$fitted.values
## Normalised linear ####
base_model <- lm(
  norm_potential ~ norm_power_est0,
  data = GB_df
)

full_model <- lm(
  norm_potential ~ tech_typ *
    norm_power_est0 +
    norm_power_est0 * month +
    hour +
    poly(ws_h_wmean, 3),
  data = GB_df
)
summary(full_model)

model_BIC <- step(
  base_model,
  scope = list(lower = base_model, upper = full_model),
  # steps = 5,
  k = log(nrow(GB_df))
)

anova(model_BIC)

model_AIC <- step(
  base_model,
  scope = list(lower = base_model, upper = full_model),
  # steps = 5,
  k = 2
)

anova(model_AIC)

plot(model_BIC)
summary(model_AIC)

GB_df$nlin.fit <- pmax(0, model_AIC$fitted.values)

GB_df %>%
  ggplot(aes(norm_potential, nlin.fit)) +
  geom_hex() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency"
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Elexon", y = "Calibrated ERA5 CF %")
ggsave("fig/gb_calib_lm_hexbin_all.pdf", width = 6, height = 4)


GB_df %>%
  ggplot(aes(potential, lin.fit)) +
  geom_hex() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  coord_fixed()


acf_df <- acf(
  model_AIC$residuals[GB_df$tech_typ == "Wind Offshore"],
  plot = TRUE,
  main = "Wind Offshore"
)
acf_df <- acf(
  model_AIC$residuals[GB_df$tech_typ == "Wind Onshore"],
  plot = TRUE,
  main = "Wind Onshore"
)

pacf_df <- pacf(
  model_AIC$residuals[GB_df$tech_typ == "Wind Offshore"],
  plot = TRUE,
  main = "Wind Offshore"
)
pacf_df <- pacf(
  model_AIC$residuals[GB_df$tech_typ == "Wind Onshore"],
  plot = TRUE,
  main = "Wind Onshore"
)


## Inlabru ####
GB_df$tech_typ <- as.factor(GB_df$tech_typ)
GB_df$month <- as.factor(GB_df$month)
GB_df$hour <- as.factor(GB_df$hour)
components <- ~ Intercept(1) +
  tech_typ(tech_typ) +
  power_est0(power_est0) +
  power_tech(power_est0, tech_typ, model = "iid") +
  month(month) +
  power_month(power_est0, month, model = "iid") +
  hour(hour) +
  ws_smooth(ws_h_wmean, model = "rw2")


full_model_bru <- bru(
  components = components,
  formula = potential ~ Intercept +
    tech_typ +
    power_est0 +
    power_tech +
    month +
    power_month +
    hour +
    ws_smooth,
  family = "gaussian",
  data = GB_df
)


# components0 <- ~ Intercept(1, prec.linear = exp(-7)) + # only declare the intercept (and any latent fields like smooths)
#   norm_power_est0(norm_power_est0, model = "linear") +
#   tech_typ(tech_typ, model = "iid")

# bru0 <- bru(
#   components = components0,
#   formula = norm_potential ~ Intercept + norm_power_est0 + tech_typ,
#   family = "gaussian",
#   data = GB_df
# )
# summary(bru0)
# brufitted <- predict(bru0)

with(GB_df, ModelMetrics::rmse(nlin.fit, norm_potential))
ModelMetrics::rmse(
  bru0$summary.fitted.values$mean[1:nrow(GB_df)],
  GB_df$norm_potential
)

components0 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  norm_power_est0(norm_power_est0, model = "linear") + # fixed slope
  tech_typ(tech_typ, model = "iid") + # random intercept by tech_typ
  tech_power(tech_typ, model = "iid", replicate = norm_power_est0) + # random slope
  month(month, model = "rw2", cyclic = TRUE) +
  # month_power(month, model = "iid", replicate = norm_power_est0) + # random slope
  hour(hour, model = "rw2", cyclic = TRUE) +
  wind(ws_group, model = "rw2")

bru0 <- bru(
  components = components0,
  formula = norm_potential ~ Intercept +
    norm_power_est0 +
    tech_typ +
    tech_power +
    month +
    hour +
    wind,
  family = "gaussian",
  data = GB_df
)
summary(bru0)
saveRDS(bru0, file = file.path(model_path, "lm_bru0.rds"))


componentsar1 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  norm_power_est0(norm_power_est0, model = "linear") + # fixed slope
  # tech_typ(tech_typ, model = "iid") + # random intercept by tech_typ
  # tech_power(tech_typ, model = "iid", replicate = norm_power_est0) + # random slope
  month(month, model = "rw2", cyclic = TRUE) +
  # hour(hour, model = "rw2", cyclic = TRUE) +
  wind(ws_group, model = "rw2") +
  u(t, model = "ar1", replicate = tech_typ)

bruar1 <- bru(
  components = componentsar1,
  formula = norm_potential ~ Intercept +
    norm_power_est0 +
    # tech_typ +
    # tech_power +
    month +
    # hour +
    wind +
    u,
  family = "gaussian",
  data = GB_df,
  options = list(
    control.predictor = list(
      compute = FALSE # don't store posterior for all latent effects
    ),
    control.compute = list(
      dic = TRUE, # you can keep DIC or WAIC if needed
      waic = TRUE,
      cpo = FALSE # skip CPO to save memory
    ),
    verbose = TRUE
  )
)
saveRDS(bruar1, file.path(model_path, "bru_ar1.rds"))
# arimax option ####
library(nlme)

full_model_ar1 <- gls(
  norm_potential ~ tech_typ *
    norm_power_est0 +
    norm_power_est0 * month +
    month +
    hour +
    poly(ws_h_wmean, 3),
  data = GB_df,
  correlation = corAR1(form = ~ t | tech_typ), # AR1 per tech type
  verbose = TRUE
)
saveRDS(full_model_ar1, file = file.path(model_path, "gls_ar1.rds"))
full_model_ar1 <- readRDS(file = file.path(model_path, "gls_ar1.rds"))
summary(full_model_ar1)
anova(full_model_ar1)
coef(full_model_ar1$modelStruct$corStruct)
phi <- coef(full_model_ar1$modelStruct$corStruct)
phi_est <- 2 / (1 + exp(-phi)) - 1 # inverse logit to [-1, 1]
phi_est


# saveRDS(full_model_ar1, file = file.path(model_path, "lm_ar1.rds"))

ModelMetrics::rmse(
  full_model_ar1$fitted,
  GB_df$norm_potential
)


full_model <- lm(
  norm_potential ~ tech_typ *
    norm_power_est0 +
    # norm_power_est0 * month +
    month +
    hour +
    poly(ws_h_wmean, 3),
  data = GB_df
)
summary(full_model)
ModelMetrics::rmse(
  full_model$fitted.values,
  GB_df$norm_potential
)

GB_df$lin.fit <- pmax(full_model$fitted.values, 0)
GB_df$lar1.fit <- pmax(full_model_ar1$fitted, 0)


GB_df %>%
  ggplot(aes(norm_potential, lin.fit)) +
  geom_hex() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  coord_fixed()

GB_df %>%
  ggplot(aes(norm_potential, lar1.fit)) +
  geom_hex() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Elexon CF", y = "ERA5 calibrated CF with AR1")
ggsave("fig/gb_calib_lar1_hexbin_all.pdf", width = 6, height = 4)


GB_df %>%
  filter(between(halfHourEndTime, "2024-01-01", "2024-01-31")) %>%
  ggplot() +
  geom_line(aes(t, norm_potential, col = "obs")) +
  geom_line(aes(t, lin.fit, col = "linear fit")) +
  geom_line(aes(t, lar1.fit, col = "ar1 fit")) +
  facet_wrap(~tech_typ, nrow = 2) +
  scale_color_aaas()


# Spatial correlation ####

pwr_curv_df <- read_parquet(file.path(gen_path, "power_curve_all.parquet"))
head(pwr_curv_df)


library(spdep)
library(sf)
library(gstat)
library(sp)
set.seed(0)
sample_hour <- sample(unique(pwr_curv_df$halfHourEndTime), 1)
hour_df <- pwr_curv_df %>%
  filter(halfHourEndTime == sample_hour) %>%
  mutate(
    norm_potential = potential / capacity,
    norm_power_est0 = power_est0 / capacity,
    error0 = norm_potential - norm_power_est0
  ) %>%
  group_by(lon, lat) %>%
  summarise(
    across(c(norm_potential, norm_power_est0, error0), sum),
    .groups = "drop"
  )

hour_df <- pwr_curv_df %>%
  filter(between(halfHourEndTime, "2024-01-01", "2024-12-31")) %>%
  group_by(lon, lat, )
mutate(
  norm_potential = potential / capacity,
  norm_power_est0 = power_est0 / capacity,
  error0 = norm_potential - norm_power_est0
) %>%
  group_by(lon, lat) %>%
  summarise(
    across(c(norm_potential, norm_power_est0, error0), sum),
    .groups = "drop"
  )

# Convert to sf object
hour_sf <- st_as_sf(hour_df, coords = c("lon", "lat"), crs = 4326)

# Create neighbors (k-nearest, k=4)
coords <- st_coordinates(hour_sf)
nb <- knn2nb(knearneigh(coords, k = 4))
lw <- nb2listw(nb, style = "W")

# Moran's I
moran.test(hour_sf$norm_potential, lw)
moran.test(hour_sf$error0, lw)


coordinates(hour_df) <- ~ lon + lat
vgram <- variogram(norm_potential ~ 1, data = hour_df)
plot(vgram)

hour_sf_proj <- st_transform(hour_sf, 27700)
# gstat variogram works with sf objects directly
vgram <- variogram(norm_potential ~ 1, data = hour_sf_proj)
vgram$dist_km <- vgram$dist / 1000
pdf("fig/variogram_norm_power.pdf", width = 5, height = 4)
plot(
  vgram$dist_km,
  vgram$gamma,
  type = "b",
  xlab = "Distance (km)",
  ylab = "Semivariance"
)
dev.off()


vgram <- variogram(norm_power_est0 ~ 1, data = hour_sf_proj)
vgram$dist_km <- vgram$dist / 1000
pdf("fig/variogram_norm_power_est.pdf", width = 5, height = 4)
plot(
  vgram$dist_km,
  vgram$gamma,
  type = "b",
  xlab = "Distance (km)",
  ylab = "Semivariance"
)
dev.off()

# Spatial correlation more times ####
n.days <- 7
time_skips <- 29
# slice data to 12 hour time steps
gb_thin <- pwr_curv_df %>%
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
    across(c(norm_potential, norm_power_est0, error0), sum),
    .groups = "drop"
  )


max_lag <- 20

pacf_df <- gb_thin %>%
  group_by(site_name) %>%
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
ggplot(pacf_df, aes(lag, site_name, fill = significant)) +
  geom_tile() +
  scale_fill_manual(values = c("grey90", "red")) +
  theme_minimal()


ggplot(pacf_df, aes(lag, site_name, fill = pacf)) +
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


library(spdep)
library(sf)
library(gstat)
library(sp)
gb_thin_sf <- st_as_sf(gb_thin, coords = c("lon", "lat"), crs = 4326)
gb_thin_proj <- st_transform(gb_thin_sf, 27700)
# gstat variogram works with sf objects directly
vgram_obs <- variogram(norm_potential ~ 1, data = gb_thin_proj, cutoff = 600e3)
vgram_obs$dist_km <- vgram_obs$dist / 1000
pdf("fig/variogram_norm_power_thin.pdf", width = 5, height = 4)
plot(
  vgram_obs$dist_km,
  vgram_obs$gamma,
  type = "b",
  xlab = "Distance (km)",
  ylab = "Semivariance"
)
dev.off()


vgram_est <- variogram(norm_power_est0 ~ 1, data = gb_thin_proj, cutoff = 600e3)
vgram_est$dist_km <- vgram_est$dist / 1000
pdf("fig/variogram_norm_power_est_thin.pdf", width = 5, height = 4)
plot(
  vgram_est$dist_km,
  vgram_est$gamma,
  type = "b",
  xlab = "Distance (km)",
  ylab = "Semivariance"
)
dev.off()

source("aux_funct.R")

corr_df <- spatial_corr_by_distance_fast(
  gb_thin,
  value_col = "norm_potential",
  n_bins = 15,
  bin_type = "quantile"
)

ggplot(corr_df, aes(x = dist_mean, y = corr_mean)) +
  geom_ribbon(
    aes(ymin = corr_lower, ymax = corr_upper),
    # alpha = 0.2,
    fill = blues9[3]
  ) +
  geom_line(color = blues9[5]) +
  geom_point(color = blues9[7]) +
  # geom_line(aes(y = threshold_95), linetype = "dashed", color = "red") +
  # geom_line(aes(y = -threshold_95), linetype = "dashed", color = "red") +
  labs(
    x = "Distance (km)",
    y = "Correlation",
    # title = "Spatial correlation vs distance with variability and significance"
  ) +
  theme_minimal()
ggsave("fig/norm_pot_corr_dist.pdf", width = 6, height = 4)

corr_df_e <- spatial_corr_by_distance_fast(
  gb_thin,
  value_col = "error0",
  n_bins = 20,
  bin_type = "quantile"
)
ggplot(corr_df_e, aes(x = dist_mean, y = corr_mean)) +
  geom_ribbon(
    aes(ymin = corr_lower, ymax = corr_upper),
    # alpha = 0.2,
    fill = blues9[3]
  ) +
  geom_line(color = blues9[5]) +
  geom_point(color = blues9[7]) +
  geom_line(aes(y = threshold_95), linetype = "dashed", color = "red") +
  # geom_line(aes(y = -threshold_95), linetype = "dashed", color = "red") +
  labs(
    x = "Distance (km)",
    y = "Correlation",
    # title = "Spatial correlation vs distance with variability and significance"
  ) +
  theme_minimal()
ggsave("fig/error0_corr_dist.pdf", width = 6, height = 4)


spatial_corr_by_distance_fast(
  gb_thin,
  value_col = "norm_power_est0",
  n_bins = 20,
  bin_type = "quantile"
) %>%
  ggplot(aes(x = dist_mean, y = corr_mean)) +
  geom_ribbon(
    aes(ymin = corr_lower, ymax = corr_upper),
    # alpha = 0.2,
    fill = blues9[3]
  ) +
  geom_line(color = blues9[5]) +
  geom_point(color = blues9[7]) +
  geom_line(aes(y = threshold_95), linetype = "dashed", color = "red") +
  # geom_line(aes(y = -threshold_95), linetype = "dashed", color = "red") +
  labs(
    x = "Distance (km)",
    y = "Correlation",
    # title = "Spatial correlation vs distance with variability and significance"
  ) +
  theme_minimal()
ggsave("fig/power_est_corr_dist.pdf", width = 6, height = 4)


# regression of spatial prop ###
lm_spatial <- lm(
  norm_potential ~ lon + lat + I(lon^2) + I(lat^2) + lon:lat,
  data = gb_thin
)
gb_thin$resid <- resid(lm_spatial)

spatial_corr_by_distance_fast(
  gb_thin,
  value_col = "resid",
  n_bins = 20,
  bin_type = "quantile"
) %>%
  ggplot(aes(x = dist_mean, y = corr_mean)) +
  geom_ribbon(
    aes(ymin = corr_lower, ymax = corr_upper),
    # alpha = 0.2,
    fill = blues9[3]
  ) +
  geom_line(color = blues9[5]) +
  geom_point(color = blues9[7]) +
  geom_line(aes(y = threshold_95), linetype = "dashed", color = "red") +
  # geom_line(aes(y = -threshold_95), linetype = "dashed", color = "red") +
  labs(
    x = "Distance (km)",
    y = "Correlation",
    # title = "Spatial correlation vs distance with variability and significance"
  ) +
  theme_minimal()

gb_thin_proj <- gb_thin %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(27700)

set.seed(1)

vgram_res <- variogram(
  resid ~ 1,
  data = gb_thin_proj,
  cutoff = 900e3,
  # width = 20e3,
  # maxlag = 15e3,
)
vgram_res$dist_km <- vgram_res$dist / 1000
pdf("fig/variogram_resid_thin.pdf", width = 5, height = 4)
plot(
  vgram_res$dist_km,
  vgram_res$gamma,
  type = "b",
  xlab = "Distance (km)",
  ylab = "Semivariance"
)
dev.off()


gb_sites <- gb_thin_proj %>%
  group_by(site_name) %>%
  summarise(resid = mean(resid, na.rm = TRUE), .groups = "drop")

vg_sites <- variogram(resid ~ 1, gb_sites, cutoff = 1e6)
vg_sites$dist_km <- vg_sites$dist / 1000
plot(
  vg_sites$dist_km,
  vg_sites$gamma,
  type = "b",
  xlab = "Distance (km)",
  ylab = "Semivariance"
)
