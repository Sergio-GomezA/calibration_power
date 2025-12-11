# Power curve model ####

# read Data ####
if (grepl("exports", getwd())) {
  # running in cluster
  cluster_run <- TRUE
  data_path <- gen_path <- "data"
  model_path <- "../calibration/model_objects"
  gen_path <- "../calibration/data"
  temp_lib <- "/exports/eddie/scratch/s2441782/calibration/lib"
  .libPaths(temp_lib)
  setwd("/exports/eddie/scratch/s2441782/calibration_power")
} else {
  cluster_run <- FALSE
  data_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
}

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
require(patchwork)

if (cluster_run) {
  n.cores <- detectCores()
  pwr_curv_df <- read_parquet(file.path(gen_path, "power_curve_all.parquet"))
} else {
  n.cores <- detectCores() - 2
  pwr_curv_df <- read_parquet(file.path(gen_path, "power_curve.parquet"))
}
source("aux_funct.R")
inla.setOption(num.threads = paste0(n.cores, ":1"))
# pc model : norm_power ~ pc(ws_h, group = site_name)
# penalty :  replicate = pc

# Basic models ####
## Smooth RW beta model #####

pwr_curv_df <- pwr_curv_df %>%
  group_by(lon, lat, halfHourEndTime) %>%
  summarise(
    site_name = first(site_name),
    ws_h = mean(ws_h),
    across(c(potential, power_est0, capacity), sum),
    .groups = "drop"
  ) %>%
  mutate(
    norm_potential = pmin(1, potential / capacity),
    norm_power_est0 = power_est0 / capacity,
    error0 = norm_potential - norm_power_est0
  ) %>%
  mutate(
    site_id = as.integer(factor(site_name)),
    pos_val = ifelse(
      norm_potential > 0 & norm_potential < 1,
      pmin(pmax(norm_potential, 1e-6), 1 - 1e-6),
      NA_real_
    ),
    is_zero = as.integer(norm_potential == 0),
    generic_logit = qlogis(
      pmin(pmax(norm_power_est0, 1e-6), 1 - 1e-6)
    )
  )

sum(pwr_curv_df$norm_potential <= 0) / nrow(pwr_curv_df) * 100
sum(pwr_curv_df$norm_potential >= 1) / nrow(pwr_curv_df) * 100

set.seed(0)
sites_vec <- pwr_curv_df$site_name %>% unique()
n_sites <- 12
sites_samp <- sample(sites_vec, n_sites)
df <- pwr_curv_df %>%
  filter(site_name %in% sites_samp) #%>%
# filter(!is.na(pos_val)) %>% # keep only rows used in beta model
# slice_sample(n = 100000)
n_groups <- 20 # start here; increase if you want a finer curve
brks <- seq(
  min(df$ws_h, na.rm = TRUE),
  max(df$ws_h, na.rm = TRUE),
  length.out = n_groups + 1
)

df$ws_group <- cut(
  df$ws_h,
  breaks = brks,
  include.lowest = TRUE,
  labels = FALSE
)
ws_midpoints <- (brks[-1] + brks[-length(brks)]) / 2

model_name <- "RW2"
print(
  sprintf("%s model --- initialisation", model_name)
)

components <- ~ Intercept(1) +
  curve(
    ws_group,
    model = "rw2",
    group = site_name,
    constr = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
  )

like_beta <- bru_obs(
  formula = pos_val ~ Intercept + curve,
  family = "beta", # or "betar" - check your INLA version
  data = df,
  control.family = list(link = "logit")
)
fit_rw2 <- bru(
  components,
  like_beta,
  options = list(
    control.inla = list(int.strategy = "auto"),
    verbose = TRUE,
    control.compute = list(config = TRUE) # if you later want posterior samples
  )
)

print(
  sprintf("%s model --- estimation finished", model_name)
)
saveRDS(
  fit_rw2,
  file = file.path(model_path, sprintf("%s-24-wfsamp.rds", model_name))
)
summary(fit_rw2)

ws_grid <- seq(min(df$ws_h), max(df$ws_h), length.out = 25)
pred_df <- expand.grid(
  ws_h = ws_grid,
  site_name = sites_samp
)
pred_df$ws_group <- sapply(pred_df$ws_h, function(x) {
  which.min(abs(ws_midpoints - x))
})
# head(pred_df)
# class(df$ws_group)

print(
  sprintf("%s model --- sampling and plotting", model_name)
)
n_samp <- 1000 # number of posterior samples
# predict linear predictor (logit mean)
pred_lp <- predict(
  fit_rw2,
  pred_df,
  ~ plogis(Intercept + curve),
  num.threads = n.cores,
  n.samples <- n_samp
)

ggplot() +
  gg(pred_lp) +
  facet_wrap(~site_name)
ggsave("fig/rw2_24_wfsamp_curve.pdf")
plot(fit_rw2, "Intercept")
ggsave("fig/rw2_24_wfsamp_intercept.pdf")

## Smooth 1D SPDE beta model ####
model_name <- "1D SPDE"
print(
  sprintf("%s model --- initialisation", model_name)
)

x <- seq(-5, 35, length.out = 20) # this sets mesh points - try others if you like
(mesh1D <- fm_mesh_1d(x, degree = 2, boundary = "dirichlet"))
ggplot() +
  geom_fm(data = mesh1D)

the_spde <- inla.spde2.pcmatern(
  mesh1D,
  prior.range = c(1, 0.01),
  prior.sigma = c(1, 0.01)
)

components_spde <- ~ Intercept(1) +
  curve(
    ws_group,
    model = the_spde,
    group = site_name
  )
like_beta <- bru_obs(
  formula = pos_val ~ Intercept + curve,
  family = "beta", # or "betar" - check your INLA version
  data = df,
  control.family = list(link = "logit")
)
fit_spde <- bru(
  components_spde,
  like_beta,
  options = list(
    control.inla = list(int.strategy = "auto"),
    verbose = TRUE,
    control.compute = list(config = TRUE) # if you later want posterior samples
  )
)
print(
  sprintf("%s model --- estimation finished", model_name)
)
saveRDS(
  fit_spde,
  file = file.path(model_path, sprintf("%s-24-wfsamp.rds", model_name))
)
summary(fit_spde)

print(
  sprintf("%s model --- sampling and plotting", model_name)
)
pred_spde <- predict(
  fit_spde,
  pred_df,
  ~ plogis(Intercept + curve),
  num.threads = n.cores,
  n.samples <- n_samp
)

ggplot() +
  gg(pred_spde) +
  facet_wrap(~site_name)
ggsave("fig/spde_24_wfsamp_curve.pdf")
plot(fit_spde, "Intercept")
ggsave("fig/spde_24_wfsamp_intercept.pdf")

spde.range <- spde.posterior(fit_spde, "curve", what = "range")
spde.logvar <- spde.posterior(fit_spde, "curve", what = "variance")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)
(range.plot | var.plot)
ggsave("fig/spde_24_wfsamp_pars.pdf")
# Models with zero inflated properties ####
## Smooth RW beta model #####

## Smooth 1D SPDE beta model ####
model_name <- "ZIB SPDE"
print(
  sprintf("%s model --- initialisation", model_name)
)
spde_beta <- inla.spde2.pcmatern(
  mesh1D,
  prior.range = c(1, 0.01),
  prior.sigma = c(1, 0.01)
)
spde_bern <- inla.spde2.pcmatern(
  mesh1D,
  prior.range = c(1, 0.01),
  prior.sigma = c(1, 0.01)
)

components_zero <- ~ Intercept_pc(1) +
  power_curve(
    ws_group,
    model = spde_beta,
    group = site_name
  ) +
  Intercept_bern(1) +
  bern_curve(
    ws_group,
    model = spde_bern,
    group = site_name
  )
like_beta <- bru_obs(
  formula = pos_val ~ Intercept_pc + power_curve,
  family = "beta",
  data = df %>% filter(!is_zero),
  control.family = list(link = "logit")
)

like_zero <- bru_obs(
  formula = is_zero ~ Intercept_bern + bern_curve,
  family = "binomial",
  data = df,
  control.family = list(link = "logit")
)

fit_zib <- bru(
  components_zero,
  like_beta,
  like_zero,
  options = list(
    control.inla = list(int.strategy = "auto"),
    verbose = TRUE,
    control.compute = list(config = TRUE)
  )
)
print(
  sprintf("%s model --- estimation finished", model_name)
)
saveRDS(
  fit_zib,
  file = file.path(model_path, sprintf("%s-24-wfsamp.rds", model_name))
)
summary(fit_zib)


print(
  sprintf("%s model --- sampling and plotting", model_name)
)
plot(fit_zib, "Intercept_bern")
ggsave("fig/ZIBspde_24_wfsamp_interZero.pdf")
plot(fit_zib, "Intercept_pc")
ggsave("fig/ZIBspde_24_wfsamp_interCurve.pdf")
# plot(fit_zib, "bern_curve")
# plot(fit_zib, "power_curve")
# ?plot.bru

pp_zero <- predict(
  fit_zib,
  newdata = pred_df,
  formula = ~ plogis(Intercept_bern + bern_curve),
  n.samples = n_samp,
  num.threads = n.cores
)
ggplot() +
  gg(pp_zero) +
  facet_wrap(~site_name)
ggsave("fig/ZIBspde_24_wfsamp_zeroProb.pdf")
pp_beta <- predict(
  fit_zib,
  newdata = pred_df,
  formula = ~ plogis(Intercept_pc + power_curve),
  n.samples = n_samp,
  num.threads = n.cores
)
ggplot() +
  gg(pp_beta) +
  facet_wrap(~site_name)
ggsave("fig/ZIBspde_24_wfsamp_curve.pdf")
# Models with penalty to generic curves ####
## Smooth RW beta model #####

## Smooth 1D SPDE beta model ####
model_name <- "PenalisedF ZIB SPDE"
print(
  sprintf("%s model --- initialisation", model_name)
)
x <- seq(-5, 35, length.out = 20) # this sets mesh points - try others if you like
(mesh1D <- fm_mesh_1d(x, degree = 2, boundary = "dirichlet"))
spde_beta <- inla.spde2.pcmatern(
  mesh1D,
  prior.range = c(1, 0.01),
  prior.sigma = c(1, 0.01)
)
spde_bern <- inla.spde2.pcmatern(
  mesh1D,
  prior.range = c(1, 0.01),
  prior.sigma = c(1, 0.01)
)

components_zero <- ~ Intercept_pc(1) +
  power_curve(
    ws_group,
    model = spde_beta,
    group = site_name
  ) +
  Intercept_bern(1) +
  bern_curve(
    ws_group,
    model = spde_bern,
    group = site_name
  )
like_beta <- bru_obs(
  formula = pos_val ~ Intercept_pc + power_curve,
  family = "beta",
  data = df %>% filter(!is_zero),
  control.family = list(link = "logit")
)

like_zero <- bru_obs(
  formula = is_zero ~ Intercept_bern + bern_curve,
  family = "binomial",
  data = df,
  control.family = list(link = "logit")
)
pseudo_precision <- 100
like_pseudo <- bru_obs(
  generic_logit ~ Intercept_pc + power_curve,
  family = "gaussian",
  data = df %>% filter(!is_zero),
  control.family = list(
    hyper = list(prec = list(initial = log(pseudo_precision), fixed = TRUE))
  )
)

model_fname <- file.path(model_path, sprintf("%s-24-wfsamp.rds", model_name))
if (!file.exists(model_fname)) {
  fit_with_penalty <- bru(
    components_zero,
    like_beta,
    like_zero,
    like_pseudo,
    options = list(
      control.inla = list(int.strategy = "auto"),
      verbose = TRUE,
      control.compute = list(config = TRUE)
    )
  )

  print(
    sprintf("%s model --- estimation finished", model_name)
  )
  saveRDS(
    fit_with_penalty,
    file = model_fname
  )
} else {
  print(
    sprintf("%s model file found --- loading file", model_name)
  )
  fit_with_penalty <- readRDS(model_fname)
}
summary(fit_with_penalty)
print(
  sprintf("%s model --- sampling and plotting", model_name)
)
plot(fit_with_penalty, "Intercept_bern")
ggsave("fig/ZIBspdePenF_24_wfsamp_interZero.pdf")
plot(fit_with_penalty, "Intercept_pc")
ggsave("fig/ZIBspdePenF_24_wfsamp_interCurve.pdf")

plot.hyper.dens(fit_with_penalty)
ggsave("fig/ZIBspdePenF_24_wfsamp_hyper.pdf")
# plot(fit_zib, "bern_curve")
# plot(fit_zib, "power_curve")
# ?plot.bru

pp_zero <- predict(
  fit_with_penalty,
  newdata = pred_df,
  formula = ~ plogis(Intercept_bern + bern_curve),
  n.samples = n_samp,
  num.threads = n.cores
)
ggplot() +
  gg(pp_zero) +
  facet_wrap(~site_name)
ggsave("fig/ZIBspdePenF_24_wfsamp_zeroProb.pdf")

pp_beta <- predict(
  fit_with_penalty,
  newdata = pred_df,
  formula = ~ plogis(Intercept_pc + power_curve),
  n.samples = n_samp,
  num.threads = n.cores
)
ggplot() +
  gg(pp_beta) +
  facet_wrap(~site_name)
ggsave("fig/ZIBspdePenF_24_wfsamp_curve.pdf")
print("Process finished")
