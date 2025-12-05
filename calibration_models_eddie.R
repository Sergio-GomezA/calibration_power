# packages

temp_lib <- "/exports/eddie/scratch/s2441782/calibration/lib"
.libPaths(temp_lib)
if (!require(inlabru)) {
  local({
    r <- c(
      INLA = "https://inla.r-inla-download.org/R/testing",
      CRAN = "https://cloud.r-project.org/",
      inlabru_universe = "https://inlabru-org.r-universe.dev"
    )
    options(repos = r)
  })
  install.packages(
    c("INLA", "inlabru"),
    temp_lib,
    dependencies = TRUE
  )
  require(INLA)
  inla.binary.install()
}
install.packages(
  c("ggsci", "terra", "fst", "arrow"),
  temp_lib,
  dependencies = TRUE
)
install.packages('R.utils', temp_lib, dependencies = TRUE)


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
inla.binary.install()

require(inlabru)

setwd("/exports/eddie/scratch/s2441782/calibration")

source("aux_funct.R")

# inla.setOption(num.threads = parallel::detectCores())
inla.setOption(num.threads = paste0("0:", parallel::detectCores()))
gen_path <- "data"
model_path <- "model_objects"

# GB_df <- read_parquet(file.path(gen_path, "GB_aggr_frag.parquet"))
GB_df <- read_parquet(file.path(gen_path, "GB_aggr.parquet"))

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
saveRDS(bruar1, file.path(model_path, "bru_ar1_v2.rds"))
ModelMetrics::rmse(
  bruar1$summary.fitted.values$mean[1:nrow(GB_df)],
  GB_df$norm_potential
)

componentsar1 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  norm_power_est0(norm_power_est0, model = "linear") + # fixed slope
  tech_typ(tech_typ, model = "iid") + # random intercept by tech_typ
  tech_power(tech_typ, model = "iid", replicate = norm_power_est0) + # random slope
  month(month, model = "rw2", cyclic = TRUE) +
  hour(hour, model = "rw2", cyclic = TRUE) +
  wind(ws_group, model = "rw2") +
  u(t, model = "ar1", replicate = tech_typ)

bruar1 <- bru(
  components = componentsar1,
  formula = norm_potential ~ Intercept +
    norm_power_est0 +
    tech_typ +
    tech_power +
    month +
    hour +
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
saveRDS(bruar1, file.path(model_path, "bru_ar1_full_v2.rds"))


componentsar2 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  norm_power_est0(norm_power_est0, model = "linear") + # fixed slope
  tech_typ(tech_typ, model = "iid") + # random intercept by tech_typ
  tech_power(tech_typ, model = "iid", replicate = norm_power_est0) + # random slope
  month(month, model = "rw2", cyclic = TRUE) +
  hour(hour, model = "rw2", cyclic = TRUE) +
  wind(ws_group, model = "rw2") +
  u(t, model = "ar", order = 2, replicate = tech_typ)

bruar2 <- bru(
  components = componentsar2,
  formula = norm_potential ~ Intercept +
    norm_power_est0 +
    tech_typ +
    tech_power +
    month +
    hour +
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
saveRDS(bruar2, file.path(model_path, "bru_ar2_full_v2.rds"))
ModelMetrics::rmse(
  bruar2$summary.fitted.values$mean[1:nrow(GB_df)],
  GB_df$norm_potential
)


bruar1 <- readRDS(file.path(model_path, "bru_ar1_full_v2.rds"))
bru_ar1_fit = bru_fitted_exclude(bruar1, GB_df, exclude = "u")
ModelMetrics::rmse(
  bru_ar1_fit,
  GB_df$norm_potential
)

bru_ar1_fit <- bruar1$summary.fitted.values$mean[1:nrow(GB_df)] -
  bruar1$summary.random$u$mean
GB_df$bar1_fit <- bru_ar1_fit
GB_df$bar1_u <- bruar1$summary.random$u$mean
GB_df %>%
  ggplot(aes(norm_potential, bar1_fit)) +
  geom_hex() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  coord_fixed() +
  labs(x = "Elexon", y = "Calibration with ar1")
ggsave("fig/gb_calib_bar1_hexbin_all_test.png", width = 6, height = 4)

GB_df %>%
  filter(between(
    halfHourEndTime,
    as.Date("2024-01-01"),
    as.Date("2024-01-31")
  )) %>%
  ggplot() +
  geom_line(aes(halfHourEndTime, norm_potential, col = "obs")) +
  geom_line(aes(halfHourEndTime, bar1_fit, col = "bru no ar1")) +
  geom_line(aes(halfHourEndTime, bar1_u, col = "ar1")) +
  facet_wrap(~tech_typ, nrow = 2) +
  scale_color_aaas() +
  labs(x = "date", y = "Capacity factor %", col = "") +
  theme(legend.position = "bottom")
ggsave("fig/gb_calib_bar1_ts_2401_test.png", width = 6, height = 4)

GB_df$bar1_fit <- bru_ar1_fit
GB_df %>%
  ggplot(aes(
    norm_potential,
    bruar1$summary.fitted.values$mean[1:nrow(GB_df)]
  )) +
  geom_hex() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  coord_fixed()
ggsave("fig/gb_calib_bar1_hexbin_all_test.png", width = 6, height = 4)


# penalised u term ####
componentsar1 <- ~ Intercept(1, prec.linear = exp(-7)) +
  norm_power_est0(norm_power_est0, model = "linear") +
  tech_typ(tech_typ, model = "iid") +
  tech_power(tech_typ, model = "iid", replicate = norm_power_est0) +
  month(month, model = "rw2", cyclic = TRUE) +
  hour(hour, model = "rw2", cyclic = TRUE) +
  wind(ws_group, model = "rw2") +
  u(
    t,
    model = "ar1",
    replicate = tech_typ,
    hyper = list(
      rho = list(prior = "pc.cor1", param = c(0.8, 0.35)), # Pr(rho > 0.8) = 0.35
      prec = list(prior = "pc.prec", param = c(0.2, 0.01)) # P(sd > 0.2) = 0.01
    )
  )
bruar1_shrunk <- bru(
  components = componentsar1, # with hyper as above
  formula = norm_potential ~ Intercept +
    norm_power_est0 +
    tech_typ +
    tech_power +
    month +
    hour +
    wind +
    u,
  family = "gaussian",
  data = GB_df,
  options = list(
    control.predictor = list(compute = FALSE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = FALSE),
    control.inla = list(strategy = "simplified.laplace"),
    verbose = TRUE
  )
)
saveRDS(bruar1_shrunk, file.path(model_path, "bru_ar1s_full_v3.rds"))
bru_ar1s_fit0 = bru_fitted_exclude(bruar1_shrunk, GB_df, exclude = "u")

bru_ar1s_fit <- bruar1_shrunk$summary.fitted.values$mean[1:nrow(GB_df)] -
  bruar1_shrunk$summary.random$u$mean
ModelMetrics::rmse(
  bru_ar1s_fit,
  GB_df$norm_potential
)
GB_df$bar1s_fit0 <- bru_ar1s_fit0
GB_df$bar1s_fit <- bru_ar1s_fit
GB_df$bar1s_u <- bruar1_shrunk$summary.random$u$mean
GB_df %>%
  ggplot(aes(norm_potential, bar1s_fit)) +
  geom_hex() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  coord_fixed() +
  labs(x = "Elexon", y = "Calibration with ar1")
ggsave("fig/gb_calib_bar1s_hexbin_all_test.png", width = 6, height = 4)

GB_df %>%
  filter(between(
    halfHourEndTime,
    as.Date("2024-01-01"),
    as.Date("2024-01-31")
  )) %>%
  ggplot() +
  geom_line(aes(halfHourEndTime, norm_potential, col = "obs")) +
  geom_line(aes(halfHourEndTime, bar1s_fit, col = "bru no ar1")) +
  geom_line(aes(halfHourEndTime, bar1s_u, col = "ar1")) +
  facet_wrap(~tech_typ, nrow = 2) +
  scale_color_aaas() +
  labs(x = "date", y = "Capacity factor %", col = "") +
  theme(legend.position = "bottom")
ggsave("fig/gb_calib_bar1s_ts_2401_test.png", width = 6, height = 4)
