# read model df0

# read pred band summary

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
