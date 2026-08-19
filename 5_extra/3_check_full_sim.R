# true data frame
prefixfull <- "elexon_flpsim"
time_sel <- wf_df_frag$time %>% max()

extension <- ifelse(local_run, local_ext, cluster_ext)
df_pattern <- sprintf("^cal_aggr_test_df_.*_%s\\.%s$", d0_tag, extension)
files_found <- list.files("data", pattern = df_pattern, full.names = TRUE)

if (!override_objects && length(files_found) > 0) {
  cat(
    "Calibration data file already exists for this day. Loading existing data.\n"
  )
  if (extension != "rds") {
    true_df <- st_read(files_found[1])
  } else {
    true_df <- readRDS(files_found[1])
  }
} else {
  cat(
    "No existing calibration data file found for this day. Preparing new data.\n"
  )

  pwr_curv_df <- read_parquet(file.path(
    gen_path,
    "power_curve_all_enriched.parquet"
  ))

  coord_list_fname <- "data/coord_list.csv"
  if (!file.exists(coord_list_fname)) {
    cat("Creating new coordinate list\n")
    set.seed(1)
    coord_list <- pwr_curv_df %>%
      distinct(coord_id, tech_typ) %>%
      arrange(coord_id)

    coord_samp <- coord_list %>%
      group_by(tech_typ) %>%
      slice_sample(prop = 0.8)

    coord_list <- coord_list %>%
      mutate(
        sampled = ifelse(coord_id %in% coord_samp$coord_id, TRUE, FALSE)
      )
    write.csv(
      coord_list,
      file = coord_list_fname,
      row.names = FALSE
    )
  } else {
    cat("Loading existing coordinate list\n")
    coord_list <- read.csv(coord_list_fname)
  }

  n.days <- 0
  # n.days.before <- 7

  true_df <- pwr_curv_df %>%
    rename(time = halfHourEndTime) %>%
    filter(time %in% time_sel) %>%
    mutate(
      date = as.Date(time),
      elevation = pmax(0, elevation),
      site_name = site_name %>%
        gsub("\\b(wind\\s*farm|wf)\\b", "", ., ignore.case = TRUE) %>%
        trimws()
    ) %>%
    # filter(date %in% sampled_days) %>%
    # filter(date >= d0 - n.days.before, date <= d0 + n.days - 1) %>%
    # filter(coord_id %in% coord_list$coord_id[coord_list$sampled]) %>%
    arrange(site_name) %>%
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
    mutate(t = difftime(time, min(time), units = "hours") %>% as.numeric()) %>%
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
      d_coast_group = inla.group(dist_coast, n = 10, method = "quantile"),
      elev_group = inla.group(elevation, n = 10, method = "quantile"),
      time_id = as.integer(factor(time)),
      date = as.Date(time),
      oos = coord_id %in% coord_list$coord_id[!coord_list$sampled],
      observed = norm_potential
    ) %>%
    left_join(
      gb_day_df %>% dplyr::select(date, p_group3),
      by = c("date" = "date")
    )

  # plot candidate anomalies, norm_potential == 0 and p_group3 == "mid"
  cutprobs3 <- c(0.25, 0.75)
  # p_quant3 <- quantile(gb_day_df$norm_power_est0, probs = cutprobs3)
  p_quant3 <- quantile(gb_day_df$norm_potential, probs = cutprobs3)

  tol <- 0.01
  norm_dist_tol <- 0.50
  true_df <- true_df %>%
    mutate(
      anomaly = case_when(
        norm_potential <= tol & norm_power_est0 >= p_quant3[1] ~ TRUE,
        norm_power_est0 >= 1 - tol & norm_potential <= p_quant3[2] ~ TRUE,
        abs(norm_power_est0 - norm_potential) >= norm_dist_tol ~ TRUE,
        TRUE ~ FALSE
      )
    )
  true_df %>%
    ggplot() +
    geom_point(
      aes(norm_potential, norm_power_est0, color = anomaly),
      size = 0.5
    ) +
    scale_color_manual(values = c("grey50", "darkred")) +
    theme(legend.position = "bottom") +
    labs(
      x = "Elexon CF (%)",
      y = "Generic PC(ERA5) CF (%)",
      color = "Anomaly"
    ) +
    coord_fixed(ratio = 1)

  ggsave(
    sprintf("fig/%s/anomalies_%s.pdf", batch_name, d0_tag),
    width = 4,
    height = 4.5
  )

  # mask anomalies for model fitting

  true_df <- true_df %>%
    mutate(
      norm_potential_orig = norm_potential,
      norm_potential = ifelse(anomaly, NA, norm_potential)
    )

  x <- true_df$pow_group %>% unique() %>% sort()
  min_jump <- min(diff(sort(x))) / diff(range(x))
  if (min_jump <= 1e-4) {
    true_df <- true_df %>%
      mutate(pow_group = inla.group(norm_power_est0, n = 20, method = "cut"))
  }

  cat("Converting coordinates to km\n")
  true_df <- true_df %>%
    st_geometry() %>%
    (\(g) g / 1000)() %>%
    st_set_geometry(true_df, .)

  true_df_fname <- sprintf(
    "data/cal_aggr_test_df_%s_%s.%s",
    "base",
    d0_tag,
    extension
  )
  if (extension == "geojson" & file.exists(true_df_fname)) {
    file.remove(true_df_fname)
  }
  if (extension != "rds") {
    st_write(
      true_df,
      true_df_fname,
      driver = driver,
      append = FALSE,
      quiet = TRUE
    )
  } else {
    saveRDS(
      true_df,
      true_df_fname
    )
  }
}
cat("Number of unique locations:", nrow(true_df %>% distinct(x, y)), "\n")
n <- nrow(true_df)
cat("Number of records in the dataset:", n, "\n")

# model fit only one hour

## 2.4 ST SPDE model ####
spde <- INLA::inla.spde2.pcmatern(
  mesh = wf.mesh,
  prior.range = c(50, 0.5), # P(range < 100 km)=0.5
  prior.sigma = c(0.2, 0.5) # P(sd > 0.2)=0.5
)

components0 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  # tech_typ(tech_typ, model = "iid") + # random intercept by tech_typ
  norm_power_est0 +
  # power_correction(
  #   pow_group,
  #   model = "rw2",
  #   # replicate = tech_typ,
  #   constr = TRUE
  # ) + # smooth correction power
  d_coast(
    d_coast_group,
    model = "rw2",
    constr = TRUE
  ) + # smooth correction distance to coast
  elev(
    elev_group,
    model = "rw2",
    constr = TRUE
  ) + # smooth correction elevation
  wind(ws_group, model = "rw2", replicate = tech_typ, constr = TRUE) + # smooth correction wind
  st_field(
    geometry,
    model = spde
    # group = time_id,
    # control.group = list(model = "ar1")
  )

model_code <- sprintf("st_bru1H_%s_%s.rds", mesh_label, d0_tag)
model_fname <- file.path(
  model_path,
  model_code
)
true_intercept <- bru0$summary.fixed["Intercept", "mean"]
if (!file.exists(model_fname) || re_run_st) {
  cat(
    "-------------------------------------------------------------------------------------------------\n"
  )
  cat("Fitting spatiotemporal model\n")
  cat(
    "-------------------------------------------------------------------------------------------------\n"
  )
  bru1H <- bru(
    components = components0,
    formula = norm_potential ~ Intercept +
      norm_power_est0 +
      # power_correction +
      d_coast +
      elev +
      wind +
      st_field,
    family = "gaussian",
    data = true_df %>%
      filter(oos == FALSE),
    options = base_bru_options
  )

  # scores_df[[model_code]] <- extract_score_model(bru1H)
  # pit_list[[model_code]] <- extract_pit_model(bru1H)

  if (save_models) {
    saveRDS(
      bru1H,
      file = model_fname
    )
  } else {
    model_list[[model_code]] <- bru1H
  }
} else {
  cat("Loading existing spatiotemporal model\n")
  bru1H <- readRDS(model_fname)
}

### summary and effect plots ####
summary(bru1H)


# residuals
fitted <- bru1H$summary.fitted.values$mean[1:nrow(true_df)]
residuals <- fitted - true_df$observed
# residuals diagnostics

ggplot() +
  geom_point(aes(x = fitted, y = residuals), color = blues9[7]) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs Fitted", x = "Fitted", y = "Residuals") +
  theme_minimal()
ggsave(
  sprintf(
    "fig/%s/%s_residuals_vs_fitted_range%d_sd%s.pdf",
    batch_name,
    prefixfull,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 4,
  height = 4
)
ggplot() +
  geom_histogram(aes(residuals), bins = 30, fill = blues9[7]) +
  labs(title = "Histogram of Residuals", x = "Residuals", y = "Frequency") +
  theme_minimal()
ggsave(
  sprintf(
    "fig/%s/%s_histogram_residuals_range%d_sd%s.pdf",
    batch_name,
    prefixfull,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 4,
  height = 4
)
ggplot() +
  geom_qq(aes(sample = residuals), color = blues9[7]) +
  geom_qq_line(aes(sample = residuals), color = "red") +
  labs(
    title = "Q-Q Plot of Residuals",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_minimal()
ggsave(
  sprintf(
    "fig/%s/%s_qqplot_residuals_range%d_sd%s.pdf",
    batch_name,
    prefixfull,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 4,
  height = 4
)
# normality test on residuals
# plot(residuals, main = "Residuals", ylab = "Residuals", xlab = "Index")
norm_test <- shapiro.test(residuals)
cat(
  "Shapiro-Wilk normality test p-value: ",
  norm_test$p.value,
  "\n"
)
cat(ifelse(
  norm_test$p.value < 0.05,
  "Residuals are not normally distributed.\n",
  "Residuals could be normally distributed.\n"
))


# validate we recover estimates similar to original model fit
int.plot <- plot(bru1H, "Intercept") +
  geom_vline(
    xintercept = true_intercept,
    col = "red",
    linetype = "dashed"
  ) +
  labs(title = "Intercept posterior", x = "Intercept", y = "Density") +
  theme_bw()
spde.range <- spde.posterior(bru1H, "st_field", what = "range")
spde.var <- spde.posterior(bru1H, "st_field", what = "variance")
range.plot <- plot(spde.range) +
  geom_vline(
    xintercept = true_range,
    col = "red",
    linetype = "dashed"
  ) +
  labs(title = "Range posterior", x = "Range", y = "Density") +
  theme_bw()
var.plot <- plot(spde.var) +
  geom_vline(
    xintercept = (true_sigma^2 + extra_noise^2),
    col = "red",
    linetype = "dashed"
  ) +
  labs(title = "Variance posterior", x = "Variance", y = "Density") +
  theme_bw()

(range.plot / var.plot / int.plot)
ggsave(
  sprintf(
    "fig/%s/%s_est_hyper_range%d_sd%s.pdf",
    batch_name,
    prefixfull,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 4,
  height = 6
)
# prediction on full data set
the_formula_str <- get_bru_formula(bru1H)
the_formula <- as.formula(paste0("~ ", the_formula_str))
pred_on_fulldf <- predict(
  bru1H,
  true_df,
  the_formula
)
true_df$fit <- pred_on_fulldf$mean
true_df$lwr <- pred_on_fulldf$`q0.025`
true_df$upr <- pred_on_fulldf$`q0.975`

true_df <- true_df %>%
  mutate(
    above_lwr = observed >= lwr,
    below_upr = observed <= upr,
    in_ci = above_lwr & below_upr
  )

# pix <- fm_pixels(mesh, dims = c(200, 200))
# pred <- predict(
#   bru1H,
#   pix,
#   the_formula
# )
# samp <- generate(
#   bru1H,
#   pix,
#   the_formula,
#   n.samples = 1
# )
# pred$sample <- samp[, 1]

# pl_posterior_mean <- ggplot() +
#   gg(pred, geom = "tile") +
#   gg(bndsim, alpha = 0) +
#   geom_sf(data = true_df, aes(col = in_ci), inherit.aes = FALSE, size = 0.5) +
#   ggtitle("Post. mean")
# pl_posterior_sample <- ggplot() +
#   gg(pred, aes(fill = sample), geom = "tile") +
#   gg(bndsim, alpha = 0) +
#   geom_sf(data = true_df, aes(col = in_ci), inherit.aes = FALSE, size = 0.5) +
#   ggtitle("Post. sample")

# # Common colour scale for the truth and estimate:
# csc <- colsc(truth$field, pred$mean, pred$sample)
# samples to aggregate

# aggregation

# samples for aggregation ####

grab_prec_name <- bru1H$.args$control.family[[1]]$hyper[[
  "theta1"
]]$output.name %>%
  gsub("[ -]+", "_", .) %>%
  as.character()

# samp_loc <- generate(
#   fit,
#   mydata,
#   ~ field + Intercept,
#   n.samples = 1000
# )
samp_loc <- generate(
  bru1H,
  true_df,
  as.formula(sprintf(
    "~ data.frame(
    geometry = geometry,
    coord_id = coord_id,
    time_id = time_id,
    lin_pred = (%s),
    prec = %s,
    lin_pred_noise = rnorm(n, mean = (%s), sd = 1/sqrt((%s))),
    lwr = qnorm(0.025, mean = (%s), sd = 1/sqrt((%s))),
    upr = qnorm(0.975, mean = (%s), sd = 1/sqrt((%s)))
  )",
    the_formula_str,
    grab_prec_name,
    the_formula_str,
    grab_prec_name,
    the_formula_str,
    grab_prec_name,
    the_formula_str,
    grab_prec_name
  )),
  n.samples = 1000
)


## get quantiles
samp_df <- lapply(
  seq_along(samp_loc),
  function(s) {
    data.frame(
      geometry = samp_loc[[s]]$geometry,
      coord_id = samp_loc[[s]]$coord_id,
      time_id = samp_loc[[s]]$time_id,
      fit = samp_loc[[s]]$lin_pred_noise,
      lwr = samp_loc[[s]]$lwr,
      upr = samp_loc[[s]]$upr
    ) %>%
      mutate(sim = s)
  }
) %>%
  bind_rows()

# variance of samples vs variance of observations check
var_samples <- samp_df %>%
  # group_by(geometry) %>%
  summarise(var_samples = var(fit)) %>%
  # ungroup() %>%
  # summarise(mean_var_samples = mean(var_samples)) %>%
  pull(var_samples)
var_observed <- true_df %>%
  summarise(var_observed = var(observed)) %>%
  pull(var_observed)
cat("Mean variance of samples:", var_samples, "\n")
cat("Variance of observed data:", var_observed, "\n")

true_df2 <- samp_df %>%
  group_by(geometry) %>%
  summarise(
    q0.025 = quantile(fit, probs = 0.025),
    median = quantile(fit, probs = 0.5),
    q0.975 = quantile(fit, probs = 0.975)
  ) %>%
  bind_cols(
    true_df %>%
      st_drop_geometry()
  ) %>%
  mutate(
    above_lwr_noise = observed >= q0.025,
    below_upr_noise = observed <= q0.975,
    in_ci_noise = above_lwr_noise & below_upr_noise
  )

# plot field again with updated in_ci variable

coverage_noise <- mean(true_df2$in_ci_noise)
coverage_noise_is <- mean(true_df2$in_ci_noise[!true_df2$oos])
coverage_noise_oos <- mean(true_df2$in_ci_noise[true_df2$oos])

## scatter plots ####
true_df2 %>%
  ggplot() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  geom_errorbar(
    aes(x = fit, ymin = lwr, ymax = upr),
    col = "blue",
    alpha = 0.5
  ) +
  geom_point(aes(fit, observed, color = in_ci)) +
  labs(x = "Posterior mean", y = "Observed data")
ggsave(
  sprintf(
    "fig/%s/%s_scatter+error_range%d_sd%s_obs.pdf",
    batch_name,
    prefixfull,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 4,
  height = 4
)

true_df2 %>%
  ggplot() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  geom_errorbar(
    aes(x = fit, ymin = q0.025, ymax = q0.975),
    col = "blue",
    alpha = 0.5
  ) +
  geom_point(aes(fit, observed, color = in_ci_noise)) +
  labs(x = "Posterior mean", y = "Observed data")
ggsave(
  sprintf(
    "fig/%s/%s_scatter+obserror_range%d_sd%s_obs_noise.pdf",
    batch_name,
    prefixfull,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 4,
  height = 4
)

## aggregation ####
aggr_samples_fsim <- lapply(
  seq_along(samp_loc),
  function(s) {
    data.frame(
      coord_id = samp_loc[[s]]$coord_id,
      time_id = samp_loc[[s]]$time_id,
      geometry = samp_loc[[s]]$geometry,
      fit = samp_loc[[s]]$lin_pred_noise
    ) %>%
      left_join(
        true_df %>%
          st_drop_geometry() %>%
          dplyr::select(coord_id, time_id, observed, capacity, norm_power_est0),
        by = c("coord_id", "time_id")
      ) %>%
      summarise(
        fit = sum(fit * capacity) / sum(capacity),
        observed = sum(observed * capacity) / sum(capacity),
        norm_power_est0 = sum(norm_power_est0 * capacity) / sum(capacity),
        fit_unweighted = mean(fit),
        observed_unweighted = mean(observed),
        norm_power_est0_unweighted = mean(norm_power_est0)
      ) %>%
      mutate(sim = s)
  }
) %>%
  bind_rows() %>%
  summarise(
    date = d0,
    variance = var(fit),
    variance_unweighted = var(fit_unweighted),
    fit = mean(fit),
    observed = mean(observed),
    norm_power_est0 = mean(norm_power_est0),
    q0.025 = quantile(fit, probs = 0.025),
    median = quantile(fit, probs = 0.5),
    q0.975 = quantile(fit, probs = 0.975),
    fit_unweighted = mean(fit_unweighted),
    observed_unweighted = mean(observed_unweighted),
    norm_power_est0_unweighted = mean(norm_power_est0_unweighted),
    q0.025_unweighted = quantile(fit_unweighted, probs = 0.025),
    median_unweighted = quantile(fit_unweighted, probs = 0.5),
    q0.975_unweighted = quantile(fit_unweighted, probs = 0.975)
  ) %>%
  bind_cols(
    data.frame(
      # observed = sum(true_df2$observed),
      # fit = mean(true_df2$fit),
      cov_nonoise = coverage,
      cov_nonoise_is = coverage_is,
      cov_nonoise_oos = coverage_oos,
      cov_noise = coverage_noise,
      cov_noise_is = coverage_noise_is,
      cov_noise_oos = coverage_noise_oos
    )
  ) %>%
  mutate(
    aggr_in_ci = ifelse(observed >= q0.025 & observed <= q0.975, 1, 0),
    aggr_in_ci_unweighted = ifelse(
      observed_unweighted >= q0.025_unweighted &
        observed_unweighted <= q0.975_unweighted,
      1,
      0
    )
  )
aggr_samples_fsim

write.csv(
  aggr_samples_fsim,
  file = sprintf(
    "summaries/%s/%s_aggr_samples_day%s_range%d_sd%s.csv",
    batch_name,
    prefixfull,
    d0_tag,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  row.names = FALSE
)

cat("Variance of aggregated samples:", aggr_samples_fsim$variance, "\n")
cat(
  "Variance of unweighted aggregated samples:",
  aggr_samples_fsim$variance_unweighted,
  "\n"
)
