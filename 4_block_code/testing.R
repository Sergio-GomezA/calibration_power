test2 <- predict(
  model_list$st_bru0_very_coarse_250817.rds,
)


test <- predict(
  brulm,
  newdata = wf_df_frag,
  ~ data.frame(
    coord_id = coord_id,
    time = time,
    lin_pred = Intercept +
      techno +
      slope +
      d_coast +
      elev +
      wind
  ),
  n.samples = 1000
)

test %>% dplyr::select(mean, sd, matches("q0."), median)
brulm$summary.fitted.values %>% head()

alphas <- c(0.01, seq(0.05, 0.95, by = 0.05), 0.99)
alphas <- c(0.025, 0.975)
debug(bru_ci_plot)
source("aux_funct.R")
test2 <- bru_ci_plot(
  bru_model = brulm,
  newdata = wf_df_frag,
  n.samples = 1000,
  show.fig = FALSE,
  alphas = alphas,
  oos_type = "time",
  family = "gaussian",
  t_start = min(wf_df_frag$time)
)
test2$wf_summary %>% head()
wf_df_frag %>%
  left_join(
    test2$wf_summary,
    by = c("coord_id", "time")
  ) %>%
  dplyr::select(site_name.x, time, coord_id, lp, fit, lwr, upr) %>%
  head()
model_df0 %>% dplyr::select(site_name, time, coord_id, fit, lwr, upr) %>% head()
brulm$summary.fitted.values %>% head()


test <- model_df0 %>%
  head(1000) %>%
  left_join(pred_band_summary$st0_m2$wf_summary, by = c("coord_id", "time"))
test %>%
  filter(!is.na(fit)) %>%
  dplyr::select(time, site_name.x, norm_potential.x, st0_m2, lwr, upr, fit)

test %>%
  filter(!is.na(fit), coord_id != 130) %>%
  dplyr::select(time, site_name.x, norm_potential.x, st0_m2, lwr, upr, fit)

test <- model_df0 %>%
  head(1000) %>%
  left_join(pred_band_summary$lm_bru$wf_summary, by = c("coord_id", "time"))
test %>%
  filter(!is.na(fit)) %>%
  dplyr::select(time, coord_id, norm_potential.x, lm_bru, fit)


wf_fig_df$model %>% unique() %>% sort()
wf_fig_df %>%
  filter(
    model == "st0_m2",
  ) %>%
  dplyr::select(
    coord_id,
    time,
    site_name,
    norm_potential,
    lp,
    fit,
    lwr,
    upr
  ) %>%
  tail() %>%
  mutate(
    lptrans = plogis(lp)
  )

# fix day
which(sampled_days == "2025-01-22")
which(sampled_days == "2025-07-30")
day_id <- 12
d0 <- sampled_days[day_id]
d0_tag <- base::format(d0, "%y%m%d")

# read model df0
output_path <- "~/Documents/elexon/caloutput"
batch_name <- "batchY25d150_v3"
mesh_label <- "very_coarse"
extension <- "rds"
model_df_fname <- sprintf(
  "%s/%s/data/fit/calibration_df_%s_%s.%s",
  output_path,
  batch_name,
  mesh_label,
  d0_tag,
  extension
)
# file.exists(model_df_fname)
model_df0_r <- readRDS(model_df_fname)

# read pred band summary
task_prefix0 <- "time"
pred_summary_fname <- sprintf(
  "%s/%s/summaries/oos/pred_band_summary_%s_%s.rds",
  output_path,
  batch_name,
  task_prefix0,
  d0_tag
)

pred_summary <- readRDS(pred_summary_fname)
pred_summary$lm_bru$wf_summary %>%
  dplyr::select(coord_id, time, lp, fit, lwr, upr) %>%
  head()

df <- model_df0_r %>%
  dplyr::select(
    coord_id,
    site_name,
    time,
    norm_potential,
    lm_bru,
    st0_m2,
    lm_beta
  ) %>%
  inner_join(
    pred_summary$lm_bru$wf_summary %>%
      dplyr::select(coord_id, time, lp, fit, lwr, upr),
    by = c("coord_id", "time")
  ) %>%
  inner_join(
    pred_summary$st0_m2$wf_summary %>%
      dplyr::select(coord_id, time, lp, fit, lwr, upr),
    by = c("coord_id", "time"),
    suffix = c("_lm_bru", "_st0_m2")
  ) %>%
  inner_join(
    pred_summary$lm_beta$wf_summary %>%
      dplyr::select(
        coord_id,
        time,
        lp_lm_beta = lp,
        fit_lm_beta = fit,
        lwr_lm_beta = lwr,
        upr_lm_beta = upr
      ),
    by = c("coord_id", "time"),
    suffix = c("", "_lm_beta")
  )
# df$time %>% range()
# d0
df %>%
  ggplot() +
  geom_point(aes(norm_potential, lm_bru), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "darkred")

df %>%
  ggplot() +
  geom_point(aes(norm_potential, st0_m2), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "darkred")

df %>%
  ggplot() +
  geom_point(aes(norm_potential, lm_beta), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "darkred")


# rmse
ModelMetrics::rmse(df$norm_potential, df$lm_bru)
ModelMetrics::rmse(df$norm_potential, df$lm_beta)
ModelMetrics::rmse(df$norm_potential, df$st0_m2)

ModelMetrics::rmse(df$norm_potential, df$fit_lm_bru)
ModelMetrics::rmse(df$norm_potential, df$fit_lm_beta)
ModelMetrics::rmse(df$norm_potential, df$fit_st0_m2)

ModelMetrics::rmse(model_df0$norm_potential, model_df0$lm_bru)
ModelMetrics::rmse(model_df0$norm_potential, model_df0$lm_beta)
ModelMetrics::rmse(model_df0$norm_potential, model_df0$st0_m2)

with(
  model_df0 %>% filter(time %in% df$time),
  ModelMetrics::rmse(norm_potential, lm_bru)
)
with(
  df,
  ModelMetrics::rmse(norm_potential, fit_lm_bru)
)

with(
  model_df0 %>% filter(time %in% df$time),
  ModelMetrics::rmse(norm_potential, lm_beta)
)
with(
  df,
  ModelMetrics::rmse(norm_potential, fit_lm_beta)
)

with(
  model_df0,
  ModelMetrics::rmse(norm_potential, st0_m2)
)
with(
  model_df0 %>% filter(time %in% df$time),
  ModelMetrics::rmse(norm_potential, st0_m2)
)
with(
  df,
  ModelMetrics::rmse(norm_potential, fit_st0_m2)
)

df %>%
  ggplot() +
  geom_point(aes(lm_bru, fit_lm_bru), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "darkred")

df %>%
  ggplot() +
  geom_point(aes(lm_beta, fit_lm_beta), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "darkred")

df %>%
  ggplot() +
  geom_point(aes(st0_m2, lp_st0_m2), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "darkred")

df %>%
  ggplot() +
  geom_point(aes(norm_potential, fit_st0_m2), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "darkred")

df %>%
  ggplot() +
  geom_point(aes(norm_potential, st0_m2), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "darkred")

sub_df <- df %>%
  dplyr::select(
    coord_id,
    time,
    site_name,
    norm_potential,
    st0_m2,
    fit_st0_m2,
    lp_st0_m2
  )


###
# post runing model_fit

model_df0 %>%
  dplyr::select(
    coord_id,
    time,
    site_name,
    norm_potential,
    st0_m2
  ) %>%
  ggplot() +
  geom_point(aes(norm_potential, st0_m2), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "darkred")

ModelMetrics::rmse(model_df0$norm_potential, model_df0$st0_m2)
d0

source("aux_funct.R")

model_name <- "st_bru0_very_coarse_250730.rds"

get_bru_formula(model_list[[model_name]])
# summary(model_list$ts_bru0_ar1_250730.rds)
pred_df <- predict(
  object = model_list[[model_name]],
  newdata = model_df0,
  formula = as.formula(sprintf(
    "~ data.frame(
    coord_id = coord_id,
    time = time,
    norm_potential = norm_potential,
    lin_pred = %s
  )",
    get_bru_formula(
      model_list[[model_name]]
    )
  )),
  n.samples = 10
)
pred_df


pred_df %>%
  ggplot() +
  geom_point(aes(norm_potential, mean), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "darkred")

alphas <- c(0.025, 0.975)
source("aux_funct.R")
pred_band <- bru_ci_plot(
  bru_model = model_list[[model_name]],
  newdata = model_df0,
  n.samples = n_samp,
  show.fig = TRUE,
  alphas = alphas,
  oos_type = "time",
  family = model_list[[model_name]]$.args$family,
  t_start = min(model_df0$time)
)

pred_band$wf_summary %>%
  ggplot() +
  geom_point(aes(norm_potential, fit), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "darkred")

ModelMetrics::rmse(
  pred_band$wf_summary$norm_potential,
  pred_band$wf_summary$fit
)


### after building prediction df
pred_df2 <- predict(
  object = model_list[[model_name]],
  newdata = wf_df_pred,
  formula = as.formula(sprintf(
    "~ data.frame(
    coord_id = coord_id,
    time = time,
    norm_potential = norm_potential,
    lin_pred = %s
  )",
    get_bru_formula(
      model_list[[model_name]]
    )
  )),
  n.samples = 10
)

pred_df2 %>%
  ggplot() +
  geom_point(aes(norm_potential, mean), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "darkred")


wf_df_pred %>%
  dplyr::select(coord_id, time, site_name, norm_potential) %>%
  head()
model_df0 %>%
  filter(time %in% wf_df_pred$time) %>%
  dplyr::select(coord_id, time, site_name, norm_potential) %>%
  head()

fit_loc <- model_df0 %>%
  distinct(coord_id, site_name)
pred_loc <- wf_df_pred %>%
  distinct(coord_id, site_name)

# identical(fit_loc %>% arrange(coord_id), pred_loc %>% arrange(coord_id))

test <- fit_loc %>%
  arrange(coord_id) %>%
  bind_cols(
    pred_loc %>%
      arrange(coord_id) %>%
      rename(pred_coord_id = coord_id, pred_site_name = site_name)
  )

# alphas <- c(0.025, 0.975)
# source("aux_funct.R")
# pred_band <- bru_ci_plot(
#   bru_model = model_list[[model_name]],
#   newdata = model_df0,
#   n.samples = n_samp,
#   show.fig = TRUE,
#   alphas = alphas,
#   oos_type = "time",
#   family = model_list[[model_name]]$.args$family,
#   t_start = min(model_df0$time)
# )

sub_df_pred <- wf_df_pred %>%
  filter(time %in% model_df0$time)
pred_df3 <- predict(
  object = model_list[[model_name]],
  newdata = sub_df_pred,
  formula = as.formula(sprintf(
    "~ data.frame(
    coord_id = coord_id,
    time = time,
    norm_potential = norm_potential,
    lin_pred = %s
  )",
    get_bru_formula(
      model_list[[model_name]]
    )
  )),
  n.samples = 10
)

pred_df3 %>%
  ggplot() +
  geom_point(aes(norm_potential, mean), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "darkred")


dfa <- sub_df_pred %>%
  dplyr::select(
    coord_id,
    time,
    time_id,
    site_name,
    norm_power_est0,
    ws_group,
    tech_typ,
    dist_coast,
    d_coast_group,
    elev_group,
    norm_potential
  )

dfb <- model_df0 %>%
  filter(time %in% sub_df_pred$time) %>%
  dplyr::select(
    coord_id,
    time,
    time_id,
    site_name,
    norm_power_est0,
    ws_group,
    tech_typ,
    dist_coast,
    d_coast_group,
    elev_group,
    norm_potential
  )
pred_dfa <- predict(
  object = model_list[[model_name]],
  newdata = dfa,
  formula = as.formula(sprintf(
    "~ data.frame(
    coord_id = coord_id,
    time = time,
    norm_potential = norm_potential,
    lin_pred = %s
  )",
    get_bru_formula(
      model_list[[model_name]]
    )
  )),
  n.samples = 10
)

pred_dfa %>%
  ggplot() +
  geom_point(aes(norm_potential, mean), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "darkred")

pred_dfb <- predict(
  object = model_list[[model_name]],
  newdata = dfb,
  formula = as.formula(sprintf(
    "~ data.frame(
    coord_id = coord_id,
    time = time,
    norm_potential = norm_potential,
    lin_pred = %s
  )",
    get_bru_formula(
      model_list[[model_name]]
    )
  )),
  n.samples = 10
)

pred_dfb %>%
  ggplot() +
  geom_point(aes(norm_potential, mean), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "darkred")

dfa
dfb
