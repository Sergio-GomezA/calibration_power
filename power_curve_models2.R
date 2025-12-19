# Power curve model ####

# Global settings ####
n_groups <- 40
n_sites <- 25
force_zero_at_endpoints <- TRUE
end_points <- c(0, 30)
outer_bounds <- c(-5, 35)
force_prob1_at_endpoints <- FALSE
train_prop <- 0.7

x <- sort(unique(c(
  end_points,
  seq(outer_bounds[1], outer_bounds[2], length.out = n_groups)
)))
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
require(ModelMetrics)
if (cluster_run) {
  n.cores <- detectCores()
  pwr_curv_df <- read_parquet(file.path(gen_path, "power_curve.parquet"))
} else {
  n.cores <- detectCores() - 2
  pwr_curv_df <- read_parquet(file.path(gen_path, "power_curve.parquet"))
}
source("aux_funct.R")
inla.setOption(num.threads = paste0(n.cores, ":1"))
# pc model : norm_power ~ pc(ws_h, group = site_name)
# penalty :  replicate = pc

# Basic models ####
model_df_fname <- file.path(
  model_path,
  "power_curve_model_data_wfsamp.parquet"
)
test_df_fname <- file.path(
  model_path,
  "power_curve_test_data_wfsamp.parquet"
)
if (!file.exists(model_df_fname)) {
  pwr_curv_df <- pwr_curv_df %>%
    group_by(lon, lat, halfHourEndTime) %>%
    summarise(
      site_name = first(site_name),
      tech_typ = first(tech_typ),
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
    ) #%>%
  # mutate(
  #   site_name = ifelse(grepl("Dogger", site_name), "Dogger Bank", site_name)
  # )

  sum(pwr_curv_df$norm_potential <= 0) / nrow(pwr_curv_df) * 100
  sum(pwr_curv_df$norm_potential >= 1) / nrow(pwr_curv_df) * 100

  set.seed(1)
  sites_vec <- pwr_curv_df$site_name %>% unique()

  sites_samp <- sample(sites_vec, n_sites)
  sites_labels <- gsub(
    "Dogger Bank A & B (was Creyke Beck A & B)",
    "Dogger Bank",
    sites_samp,
    fixed = TRUE
  )
  names(sites_labels) <- sites_samp

  set.seed(123) # for reproducibility
  df <- pwr_curv_df %>%
    filter(site_name %in% sites_samp) %>%
    group_by(site_name) %>%
    slice_sample(prop = train_prop) %>%
    ungroup()

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
  write_parquet(
    df,
    model_df_fname
  )
  test_df <- pwr_curv_df %>%
    filter(site_name %in% sites_samp) %>%
    anti_join(df, by = colnames(pwr_curv_df))
  test_df$ws_group <- cut(
    test_df$ws_h,
    breaks = brks,
    include.lowest = TRUE,
    labels = FALSE
  )
  write_parquet(
    test_df,
    test_df_fname
  )
} else {
  df <- read_parquet(model_df_fname)
  brks <- seq(
    min(df$ws_h, na.rm = TRUE),
    max(df$ws_h, na.rm = TRUE),
    length.out = n_groups + 1
  )
  sites_samp <- df$site_name %>% unique()
  sites_labels <- gsub(
    "Dogger Bank A & B (was Creyke Beck A & B)",
    "Dogger Bank",
    sites_samp,
    fixed = TRUE
  )
  names(sites_labels) <- sites_samp

  test_df <- read_parquet(test_df_fname)
}

ws_midpoints <- (brks[-1] + brks[-length(brks)]) / 2

pred_df <- expand.grid(
  ws_h = x,
  site_name = sites_samp
)
pred_df$ws_group <- sapply(pred_df$ws_h, function(x) {
  which.min(abs(ws_midpoints - x))
})
n_samp <- 1000 # number of posterior samples
## Smooth RW beta model #####

model_name <- "RW2"
model_code <- ("rw2")
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

model_fname <- file.path(model_path, sprintf("%s-24-wfsamp.rds", model_name))
if (!file.exists(model_fname)) {
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
} else {
  print(
    sprintf("%s model file found --- loading file", model_name)
  )
  fit_rw2 <- readRDS(model_fname)
}
summary(fit_rw2)


# head(pred_df)
# class(df$ws_group)

print(
  sprintf("%s model --- sampling and plotting", model_name)
)

# predict linear predictor (logit mean)
pred_lp <- predict(
  fit_rw2,
  pred_df,
  ~ plogis(Intercept + curve),
  num.threads = n.cores,
  n.samples <- n_samp
)

ggplot() +
  # geom_point(data = df, aes(x = ws_h, pos_val), alpha = 0.2) +
  geom_hex(data = df, aes(x = ws_h, pos_val)) +
  # gg(pred_lp) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency"
  ) +
  facet_wrap(~site_name, labeller = as_labeller(sites_labels)) +
  labs(x = "ERA 5 wind speed", y = "Normalised power output (P>0)")
ggsave(sprintf("fig/%s_24_wfsamp_curve.pdf", model_code))
plot(fit_rw2, "Intercept")
ggsave(sprintf("fig/%s_24_wfsamp_intercept.pdf", model_code))
plot.hyper.dens(fit_rw2)
ggsave(sprintf("fig/%s_24_wfsamp_hyper.pdf", model_code))

## Smooth 1D SPDE beta model ####
model_name <- "1D SPDE"
model_code <- ("spde")
print(
  sprintf("%s model --- initialisation", model_name)
)

# x <- sort(unique(c(
#   end_points,
#   seq(outer_bounds[1], outer_bounds[2], length.out = n_groups)
# )))
x <- 1:40
mesh1D <- fm_mesh_1d(x, degree = 2, boundary = c("dirichlet"))

ggplot() +
  geom_fm(data = mesh1D)

if (force_zero_at_endpoints) {
  A <- matrix(0, nrow = 2, ncol = mesh1D$m)
  zero_idx <- which.min(abs(brks - end_points[1]))
  A[1, zero_idx] <- 1
  zero_idx <- min(mesh1D$m, which.min(abs(brks - end_points[2])))
  A[2, zero_idx] <- 1
  # undebug(inla.spde2.pcmatern)
  the_spde <- inla.spde2.pcmatern(
    mesh1D,
    prior.range = c(1, 0.01),
    prior.sigma = c(1, 0.01),
    extraconstr = list(
      A = A,
      e = rep(1e-6, nrow(A))
    )
  )
} else {
  the_spde <- inla.spde2.pcmatern(
    mesh1D,
    prior.range = c(1, 0.01),
    prior.sigma = c(1, 0.01)
  )
}
# ?inla.spde2.pcmatern
# inla.doc("^rw2$")
# inla.doc("^inla.spde2.pcmatern$")
components_spde <- ~ Intercept(1) +
  curve(
    ws_group,
    model = the_spde,
    group = site_name
  )
like_beta <- bru_obs(
  formula = pos_val ~ Intercept + curve,
  # formula = pos_val ~ curve,
  family = "beta", # or "betar" - check your INLA version
  data = df,
  control.family = list(link = "logit")
)
model_fname <- file.path(model_path, sprintf("%s-24-wfsamp.rds", model_name))
if (!file.exists(model_fname)) {
  fit_spde <- bru(
    components_spde,
    like_beta,
    options = list(
      control.inla = list(int.strategy = "auto"),
      verbose = FALSE,
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
} else {
  print(
    sprintf("%s model file found --- loading file", model_name)
  )
  fit_spde <- readRDS(model_fname)
}
summary(fit_spde)

print(
  sprintf("%s model --- sampling and plotting", model_name)
)
pred_spde <- predict(
  fit_spde,
  pred_df,
  ~ plogis(Intercept + curve),
  # ~ plogis(curve),
  num.threads = n.cores,
  n.samples <- n_samp
)

ggplot() +
  geom_hex(data = df, aes(x = ws_h, norm_potential)) +
  gg(pred_spde) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency"
  ) +
  facet_wrap(~site_name, labeller = as_labeller(sites_labels)) +
  labs(x = "ERA 5 wind speed", y = "Normalised power output (P>0)")
ggsave(sprintf("fig/%s_24_wfsamp_curve.pdf", model_code))
plot(fit_spde, "Intercept")
ggsave(sprintf("fig/%s_24_wfsamp_intercept.pdf", model_code))

spde.range <- spde.posterior(fit_spde, "curve", what = "range")
spde.logvar <- spde.posterior(fit_spde, "curve", what = "variance")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)
(range.plot | var.plot)
ggsave(sprintf("fig/%s_24_wfsamp_pars.pdf", model_code))
plot.hyper.dens(fit_spde)
ggsave(sprintf("fig/%s_24_wfsamp_hyper.pdf", model_code))
# Models with zero inflated properties ####
## Smooth RW beta model #####

## Smooth 1D SPDE beta model ####
model_name <- "ZIB SPDE"
model_code <- ("ZIBspde")
print(
  sprintf("%s model --- initialisation", model_name)
)
spde_beta <- the_spde

if (force_prob1_at_endpoints) {
  A <- matrix(0, nrow = 2, ncol = mesh1D$m)
  zero_idx <- which.min(abs(brks - end_points[1]))
  A[1, zero_idx] <- 1
  zero_idx <- min(mesh1D$m, which.min(abs(brks - end_points[2])))
  A[2, zero_idx] <- 1
  # undebug(inla.spde2.pcmatern)
  spde_bern <- inla.spde2.pcmatern(
    mesh1D,
    prior.range = c(1, 0.01),
    prior.sigma = c(1, 0.01),
    extraconstr = list(
      A = A,
      e = rep(1 - 1e-6, nrow(A))
    )
  )
} else {
  spde_bern <- inla.spde2.pcmatern(
    mesh1D,
    prior.range = c(1, 0.01),
    prior.sigma = c(1, 0.01)
  )
}

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

model_fname <- file.path(model_path, sprintf("%s-24-wfsamp.rds", model_name))
if (!file.exists(model_fname)) {
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
} else {
  print(
    sprintf("%s model file found --- loading file", model_name)
  )
  fit_zib <- readRDS(model_fname)
}
summary(fit_zib)

print(
  sprintf("%s model --- sampling and plotting", model_name)
)
plot(fit_zib, "Intercept_bern")
ggsave(sprintf("fig/%s_24_wfsamp_interZero.pdf", model_code))
plot(fit_zib, "Intercept_pc")
ggsave(sprintf("fig/%s_24_wfsamp_interCurve.pdf", model_code))
# spde.range <- spde.posterior(fit_zib, "curve", what = "range")
# spde.logvar <- spde.posterior(fit_zib, "curve", what = "variance")
# range.plot <- plot(spde.range)
# var.plot <- plot(spde.logvar)
# (range.plot | var.plot)
# ggsave(sprintf("fig/%s_24_wfsamp_pars.pdf", model_code))
plot.hyper.dens(fit_zib)
ggsave(sprintf("fig/%s_24_wfsamp_hyper.pdf", model_code))
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
  geom_col(
    data = df %>%
      group_by(site_name, ws_group) %>%
      summarise(
        prop_zero = mean(is_zero, na.rm = TRUE),
        ws_hmean = mean(ws_h, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ),
    aes(x = ws_hmean, prop_zero, fill = n)
  ) +
  gg(pp_zero) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency"
  ) +
  facet_wrap(~site_name, labeller = as_labeller(sites_labels)) +
  labs(x = "ERA 5 wind speed", y = "Probability of zero generation")
ggsave(
  sprintf("fig/%s_24_wfsamp_zeroProb.pdf", model_code),
  width = 8,
  height = 6
)
pp_beta <- predict(
  fit_zib,
  newdata = pred_df,
  formula = ~ plogis(Intercept_pc + power_curve),
  n.samples = n_samp,
  num.threads = n.cores
)
ggplot() +
  geom_hex(data = df, aes(x = ws_h, pos_val)) +
  gg(pp_beta) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency"
  ) +
  facet_wrap(~site_name, labeller = as_labeller(sites_labels)) +
  labs(x = "ERA 5 wind speed", y = "Normalised power output (P>0)")
ggsave(sprintf("fig/%s_24_wfsamp_curve.pdf", model_code), width = 8, height = 6)

pp_EPC <- predict(
  fit_zib,
  newdata = pred_df,
  formula = ~ (1 - plogis(Intercept_bern + bern_curve)) *
    plogis(Intercept_pc + power_curve),
  n.samples = n_samp,
  num.threads = n.cores
)
ggplot() +
  geom_hex(data = df, aes(x = ws_h, norm_potential)) +
  gg(pp_EPC) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency"
  ) +
  facet_wrap(~site_name, labeller = as_labeller(sites_labels)) +
  labs(x = "ERA 5 wind speed", y = "Expected normalised power %")
ggsave(sprintf("fig/%s_24_wfsamp_EPC.pdf", model_code), width = 8, height = 6)

# Models with penalty to generic curves ####
## Smooth RW beta model #####

## Smooth 1D SPDE beta model ####
pseudo_precision <- 0.1
fixed_penalty <- FALSE
model_name <- sprintf(
  "Penalised%0.2f%s ZIB SPDE",
  pseudo_precision,
  ifelse(fixed_penalty, "", "free")
)
model_code <- sprintf(
  "ZIBspdePen%0.2f%s",
  pseudo_precision,
  ifelse(fixed_penalty, "", "free")
)

print(
  sprintf("%s model --- initialisation", model_name)
)
x <- seq(-5, 35, length.out = 20) # this sets mesh points - try others if you like
(mesh1D <- fm_mesh_1d(x, degree = 2, boundary = "dirichlet"))
# spde_beta <- inla.spde2.pcmatern(
#   mesh1D,
#   prior.range = c(1, 0.01),
#   prior.sigma = c(1, 0.01)
# )
# spde_bern <- inla.spde2.pcmatern(
#   mesh1D,
#   prior.range = c(1, 0.01),
#   prior.sigma = c(1, 0.01)
# )

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

like_pseudo <- bru_obs(
  generic_logit ~ Intercept_pc + power_curve,
  family = "gaussian",
  data = df %>% filter(!is_zero),
  control.family = list(
    hyper = list(
      prec = list(initial = log(pseudo_precision), fixed = fixed_penalty)
    )
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
ggsave(sprintf("fig/%s_24_wfsamp_interZero.pdf", model_code))
plot(fit_with_penalty, "Intercept_pc")
ggsave(sprintf("fig/%s_24_wfsamp_interCurve.pdf", model_code))
# spde.range <- spde.posterior(fit_with_penalty, "curve", what = "range")
# spde.logvar <- spde.posterior(fit_with_penalty, "curve", what = "variance")
# range.plot <- plot(spde.range)
# var.plot <- plot(spde.logvar)
# (range.plot | var.plot)
# ggsave(sprintf("fig/%s_24_wfsamp_pars.pdf", model_code))
plot.hyper.dens(fit_with_penalty)
ggsave(sprintf("fig/%s_24_wfsamp_hyper.pdf", model_code))
# plot(fit_zib, "bern_curve")
# plot(fit_zib, "power_curve")
# ?plot.bru

pp_zero_pen <- predict(
  fit_with_penalty,
  newdata = pred_df,
  formula = ~ plogis(Intercept_bern + bern_curve),
  n.samples = n_samp,
  num.threads = n.cores
)
ggplot() +
  geom_col(
    data = df %>%
      group_by(site_name, ws_group) %>%
      summarise(
        prop_zero = mean(is_zero, na.rm = TRUE),
        ws_hmean = mean(ws_h, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ),
    aes(x = ws_hmean, prop_zero, fill = n)
  ) +
  gg(pp_zero_pen) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency"
  ) +
  facet_wrap(~site_name, labeller = as_labeller(sites_labels)) +
  labs(x = "ERA 5 wind speed", y = "Probability of zero generation")
ggsave(
  sprintf("fig/%s_24_wfsamp_zeroProb.pdf", model_code),
  width = 8,
  height = 6
)
pp_beta_pen <- predict(
  fit_with_penalty,
  newdata = pred_df,
  formula = ~ plogis(Intercept_pc + power_curve),
  n.samples = n_samp,
  num.threads = n.cores
)
ggplot() +
  geom_hex(data = df, aes(x = ws_h, pos_val)) +
  gg(pp_beta_pen) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency"
  ) +
  facet_wrap(~site_name, labeller = as_labeller(sites_labels)) +
  labs(x = "ERA 5 wind speed", y = "Normalised power output (P>0)")
ggsave(sprintf("fig/%s_24_wfsamp_curve.pdf", model_code), width = 8, height = 6)

pp_EPC_pen <- predict(
  fit_with_penalty,
  newdata = pred_df,
  formula = ~ (1 - plogis(Intercept_bern + bern_curve)) *
    plogis(Intercept_pc + power_curve),
  n.samples = n_samp,
  num.threads = n.cores
)
ggplot() +
  geom_hex(data = df, aes(x = ws_h, norm_potential)) +
  gg(pp_EPC_pen) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency"
  ) +
  facet_wrap(~site_name) +
  labs(x = "ERA 5 wind speed", y = "Expected normalised power %")
ggsave(sprintf("fig/%s_24_wfsamp_EPC.pdf", model_code), width = 8, height = 6)
print("Modelling process finished")


# Comparative figures ####

## Power Curves ####
mod_names <- c("B-RW2", "B-SPDE", rep("ZIB-SPDE", 3), rep("ZIB-SPDE-P", 3))
components <- c(
  rep("power curve", 2),
  rep(c("zero prob", "power curve", "EPC"), 2)
)
mod_names_short <- mod_names %>% unique()
pc_pred_fname <- file.path(model_path, "PC_model_wfsamp_comparison.parquet")

if (!file.exists(pc_pred_fname)) {
  pred_df_list <- list(
    pred_lp,
    pred_spde,
    pp_zero,
    pp_beta,
    pp_EPC,
    pp_zero_pen,
    pp_beta_pen,
    pp_EPC_pen
  )

  df_pc_pred <- lapply(
    seq_along(mod_names),
    \(i) {
      pred_df_list[[i]] %>%
        mutate(
          model = mod_names[i],
          component = components[i]
        )
    }
  ) %>%
    bind_rows()
  arrow::write_parquet(
    df_pc_pred,
    pc_pred_fname
  )
} else {
  df_pc_pred <- read_parquet(
    pc_pred_fname
  )
}


sites_subsamp <- "West of Duddon Sands" #sites_samp[23]
df_pc_pred_1wf <- df_pc_pred %>%
  filter(site_name %in% sites_subsamp, component == "power curve")


ggplot() +
  geom_hex(
    data = df %>% filter(site_name %in% sites_subsamp),
    aes(x = ws_h, pos_val)
  ) +
  gg(df_pc_pred_1wf) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency"
  ) +
  facet_wrap(~model) +
  labs(x = "ERA 5 wind speed", y = "Normalised power output (P>0)")
ggsave(
  sprintf("fig/%s_24_1wf_curve_raw.pdf", "comparison"),
  width = 8,
  height = 6
)

df_pc_pred_1wf <- df_pc_pred %>%
  filter(
    site_name %in% sites_subsamp,
    component == "EPC" | (model %in% c(mod_names[1:2]))
  )
ggplot() +
  geom_hex(
    data = df %>% filter(site_name %in% sites_subsamp),
    aes(x = ws_h, pos_val)
  ) +
  gg(df_pc_pred_1wf) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency"
  ) +
  facet_wrap(~model) +
  labs(x = "ERA 5 wind speed", y = "Normalised power output (P>0)")
ggsave(
  sprintf("fig/%s_24_1wf_curve.pdf", "comparison"),
  width = 8,
  height = 6
)

## Calibration log frequency plots ####
n <- nrow(df)

# fitted_EPC <- pp_EPC %>%
#   select(site_name, ws_h, mean) %>%
#   rename(fit_zib = mean)

# fitted_EPC <- pp_EPC_pen %>%
#   select(site_name, ws_h, mean) %>%
#   rename(fit_zibpen = mean)

pp_EPC <- df_pc_pred %>%
  filter(component == "EPC", model == "ZIB-SPDE")

pp_EPC_pen <- df_pc_pred %>%
  filter(component == "EPC", model == "ZIB-SPDE-P")

pcmod_obs_df <- df %>%
  mutate(
    fit_rw2 = fit_rw2$summary.fitted.values$mean[1:n],
    fit_spde = fit_spde$summary.fitted.values$mean[1:n]
  )
pcmod_obs_df$fit_zib <- pcmod_obs_df$fit_zibp <- 0
for (site in sites_samp) {
  site_ind <- which(df$site_name == site)
  pred_site_ind <- which(pp_EPC$site_name == site)

  pcmod_obs_df$fit_zib[site_ind] <- approx(
    x = pp_EPC$ws_h[pred_site_ind],
    y = pp_EPC$mean[pred_site_ind],
    xout = pcmod_obs_df$ws_h[site_ind],
    rule = 2
  )$y

  pcmod_obs_df$fit_zibp[site_ind] <- approx(
    x = pp_EPC_pen$ws_h[pred_site_ind],
    y = pp_EPC_pen$mean[pred_site_ind],
    xout = pcmod_obs_df$ws_h[site_ind],
    rule = 2
  )$y
}


ws_cols <- c("ws_h", "wd10", "wd100")
ws_cols <- c("ws_h")
aggr_pcmod_df <- pcmod_obs_df %>%
  # filter(halfHourEndTime == as.POSIXct("2024-01-01 01:00:00", tz = "UTC")) %>%
  group_by(tech_typ, halfHourEndTime) %>%
  # mutate(across(matches("fit_"), ~ . * capacity)) %>%
  summarise(
    # ws_h_wmean = sum(ws_h * capacity),
    # across(c(ws_h, wd10, wd100), ~ sum(. * capacity), .names = "{.col}_wmean"),
    across(any_of(ws_cols), ~ sum(. * capacity), .names = "{.col}_wmean"),
    across(matches("fit_"), ~ sum(. * capacity)),
    across(
      c(power_est0, potential, capacity),
      sum
    ),
    across(any_of(ws_cols), mean, .names = "{.col}_mean")
  ) %>%
  mutate(
    across(
      c(power_est0, potential),
      ~ . / capacity,
      .names = "norm_{.col}"
    ),
    # ws_h_wmean = ws_h_wmean / capacity,
    across(matches("_wmean|fit_"), ~ . / capacity),
    month = factor(month(halfHourEndTime)),
    hour = factor(hour(halfHourEndTime))
  ) %>%
  mutate(ws_group = inla.group(ws_h_wmean, n = 20, method = "quantile")) %>%
  group_by(tech_typ) %>%
  arrange(halfHourEndTime, .by_group = TRUE) %>%
  mutate(t = row_number()) %>%
  ungroup()

write_parquet(aggr_pcmod_df, "data/pcmodel_wfsamp_aggregated.parquet")
aggr_pcmod_df <- read_parquet("data/pcmodel_wfsamp_aggregated.parquet")
mod_labels <- c("Generic PC", mod_names_short)
est_cols <- c(
  "norm_power_est0",
  "fit_rw2",
  "fit_spde",
  "fit_zib",
  "fit_zibp"
)
names(mod_labels) <- est_cols
aggr_pcmod_df %>%
  pivot_longer(
    cols = c(norm_power_est0, starts_with("fit_")),
    names_to = "model",
    values_to = "fitted"
  ) %>%
  ggplot(aes(norm_potential, fitted)) +
  geom_hex() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency"
  ) +
  facet_wrap(~model, labeller = as_labeller(mod_labels)) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Elexon", y = "ERA5+PC model CF %")
ggsave("fig/gbwfsamp_pcmodel_hexbin_2024.pdf", width = 6, height = 4)

# Comparative scores ####
df_long <- aggr_pcmod_df %>%
  select(norm_potential, all_of(est_cols)) %>%
  pivot_longer(
    cols = all_of(est_cols),
    names_to = "model",
    values_to = "estimate"
  )
metrics_table <- df_long %>%
  group_by(model) %>%
  summarise(
    RMSE = rmse(actual = norm_potential, predicted = estimate),
    MAE = mae(actual = norm_potential, predicted = estimate),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(model = mod_labels[model])


# Out of sample evaluation ####

## predicted values ####

pred_rw2 <- df_pc_pred %>%
  filter(component == "power curve", model == "B-RW2")

pred_spde <- df_pc_pred %>%
  filter(component == "power curve", model == "B-SPDE")

pp_EPC <- df_pc_pred %>%
  filter(component == "EPC", model == "ZIB-SPDE")

pp_EPC_pen <- df_pc_pred %>%
  filter(component == "EPC", model == "ZIB-SPDE-P")

pcmod_obs_tdf <- test_df %>%
  mutate(
    pred_rw2 = 0,
    pred_spde = 0,
    pred_zib = 0,
    pred_zibp = 0
  )

for (site in sites_samp) {
  site_ind <- which(test_df$site_name == site)
  pred_site_ind <- which(pp_EPC$site_name == site)

  pcmod_obs_tdf$pred_rw2[site_ind] <- approx(
    x = pred_rw2$ws_h[pred_site_ind],
    y = pred_rw2$mean[pred_site_ind],
    xout = pcmod_obs_tdf$ws_h[site_ind],
    rule = 2
  )$y

  pcmod_obs_tdf$pred_spde[site_ind] <- approx(
    x = pred_spde$ws_h[pred_site_ind],
    y = pred_spde$mean[pred_site_ind],
    xout = pcmod_obs_tdf$ws_h[site_ind],
    rule = 2
  )$y

  pcmod_obs_tdf$pred_zib[site_ind] <- approx(
    x = pp_EPC$ws_h[pred_site_ind],
    y = pp_EPC$mean[pred_site_ind],
    xout = pcmod_obs_tdf$ws_h[site_ind],
    rule = 2
  )$y

  pcmod_obs_tdf$pred_zibp[site_ind] <- approx(
    x = pp_EPC_pen$ws_h[pred_site_ind],
    y = pp_EPC_pen$mean[pred_site_ind],
    xout = pcmod_obs_tdf$ws_h[site_ind],
    rule = 2
  )$y
}

ws_cols <- c("ws_h", "wd10", "wd100")
ws_cols <- c("ws_h")
aggr_pcmod_tdf <- pcmod_obs_tdf %>%
  # filter(halfHourEndTime == as.POSIXct("2024-01-01 01:00:00", tz = "UTC")) %>%
  group_by(tech_typ, halfHourEndTime) %>%
  # mutate(across(matches("fit_"), ~ . * capacity)) %>%
  summarise(
    # ws_h_wmean = sum(ws_h * capacity),
    # across(c(ws_h, wd10, wd100), ~ sum(. * capacity), .names = "{.col}_wmean"),
    across(any_of(ws_cols), ~ sum(. * capacity), .names = "{.col}_wmean"),
    across(matches("pred_"), ~ sum(. * capacity)),
    across(
      c(power_est0, potential, capacity),
      sum
    ),
    across(any_of(ws_cols), mean, .names = "{.col}_mean")
  ) %>%
  mutate(
    across(
      c(power_est0, potential),
      ~ . / capacity,
      .names = "norm_{.col}"
    ),
    # ws_h_wmean = ws_h_wmean / capacity,
    across(matches("_wmean|pred_"), ~ . / capacity),
    month = factor(month(halfHourEndTime)),
    hour = factor(hour(halfHourEndTime))
  ) %>%
  mutate(ws_group = inla.group(ws_h_wmean, n = 20, method = "quantile")) %>%
  group_by(tech_typ) %>%
  arrange(halfHourEndTime, .by_group = TRUE) %>%
  mutate(t = row_number()) %>%
  ungroup()

write_parquet(aggr_pcmod_tdf, "data/pcmodel_wfsamp_aggr_test.parquet")
aggr_pcmod_tdf <- read_parquet("data/pcmodel_wfsamp_aggr_test.parquet")

pred_cols <- gsub("fit_", "pred_", est_cols)
names(mod_labels) <- pred_cols
tdf_long <- aggr_pcmod_tdf %>%
  select(norm_potential, all_of(pred_cols)) %>%
  pivot_longer(
    cols = all_of(pred_cols),
    names_to = "model",
    values_to = "estimate"
  )
tmetrics_table <- tdf_long %>%
  group_by(model) %>%
  summarise(
    RMSE = rmse(actual = norm_potential, predicted = estimate),
    MAE = mae(actual = norm_potential, predicted = estimate),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(model = mod_labels[model])

tmetrics_table
