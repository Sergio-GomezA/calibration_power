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

require(brms)
require(INLA)

require(inlabru)

source("aux_funct.R")
local({
  r <- c(
    INLA = "https://inla.r-inla-download.org/R/testing",
    CRAN = "https://cloud.r-project.org/",
    inlabru_universe = "https://inlabru-org.r-universe.dev"
  )
  options(repos = r)
})

install.packages(c("INLA", "inlabru"), dependencies = TRUE)
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
  filter(halfHourEndTime > "2023-01-01") %>%
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
  ) %>%
  mutate(ws_group = inla.group(ws_h_wmean, n = 20, method = "quantile")) %>%
  group_by(tech_typ) %>%
  arrange(halfHourEndTime, .by_group = TRUE) %>%
  mutate(t = row_number()) %>%
  ungroup()

write_parquet(GB_df, file.path(gen_path, "GB_aggr.parquet"))
GB_df <- read_parquet(file.path(gen_path, "GB_aggr.parquet"))

write_parquet(GB_df, file.path(gen_path, "GB_aggr_frag.parquet"))
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

GB_df$nlin.fit <- pmax(0, model_AIC$fitted.values)

GB_df %>%
  ggplot(aes(norm_potential, nlin.fit)) +
  geom_hex() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency"
  ) +
  coord_fixed()


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
    # norm_power_est0 * month +
    month +
    hour +
    poly(ws_h_wmean, 3),
  data = GB_df,
  correlation = corAR1(form = ~ t | tech_typ) # AR1 per tech type
)
summary(full_model_ar1)
anova(full_model_ar1)
coef(full_model_ar1$modelStruct$corStruct)
phi <- coef(full_model_ar1$modelStruct$corStruct)
phi_est <- 2 / (1 + exp(-phi)) - 1 # inverse logit to [-1, 1]
phi_est


saveRDS(full_model_ar1, file = file.path(model_path, "lm_ar1.rds"))

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
  coord_fixed()


GB_df %>%
  filter(between(halfHourEndTime, "2024-01-01", "2024-01-31")) %>%
  ggplot() +
  geom_line(aes(t, norm_potential, col = "obs")) +
  geom_line(aes(t, lin.fit, col = "linear fit")) +
  geom_line(aes(t, lar1.fit, col = "ar1 fit")) +
  facet_wrap(~tech_typ, nrow = 2) +
  scale_color_aaas()
