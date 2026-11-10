local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

pow_threshold <- 0.05
pow_threshold_label <- gsub("\\.", "_", as.character(pow_threshold))
tables_path <- "~/ownCloud-s2441782@pwr_curv_dfasync.ed.ac.uk/projects/calibration/calibration_power_main_doc/tables"
dir.exists(tables_path) || dir.create(tables_path, recursive = TRUE)

override_objects <- FALSE
batch_name <- "bias_eda_test"
norm_dist_tol <- 1

if (local_run) {
  cat("Running in local mode\n")
} else {
  cat("Running in cluster mode\n")
}

# 0.2 libraries and paths ####
require(parallel)

if (local_run) {
  pwr_curv_dfa_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
  extra_path <- "~/Documents/elexon/extra"
  output_path <- "~/Documents/elexon/caloutput"
  sample_path <- "~/Documents/elexon/samples"
  n_samp <- 100
  pixel_dims <- c(150, 150)
} else {
  pwr_curv_dfa_path <- "/exports/eddie/scratch/s2441782/calibration/pwr_curv_dfa"
  gen_path <- "/exports/eddie/scratch/s2441782/calibration/pwr_curv_dfa"
  model_path <- "/exports/eddie/scratch/s2441782/calibration/model_objects"
  extra_path <- "/exports/eddie/scratch/s2441782/calibration/extra"
  output_path <- "/exports/eddie/scratch/s2441782/calibration/caloutput"
  sample_path <- "/exports/eddie/scratch/s2441782/calibration/samples"
  temp_lib <- "/exports/eddie3_homes_local/s2441782/lib"
  pixel_dims <- c(300, 300)
  n_samp <- 1000
  .libPaths(temp_lib)
}

require(tidyverse)
require(sf)
require(INLA)
require(inlabru)
require(fmesher)
require(ggspatial)
require(ModelMetrics)
require(qmap)
require(ggridges)
require(ggthemes)
require(ggsci)
require(arrow)
require(kableExtra)
# require(ggspatial)

source("aux_funct.R")
mc <- ifelse(local_run, 1, available_cores())

lwe_output_figpath <- file.path("fig", "lwe_eda")
dir.exists(lwe_output_figpath) ||
  dir.create(lwe_output_figpath, recursive = TRUE)
lwe_output_tabpath <- file.path("tables", "lwe_eda")
dir.exists(lwe_output_tabpath) ||
  dir.create(lwe_output_tabpath, recursive = TRUE)

# pwr_curv_dfa

extension <- "rds"
coord_list_fname <- "pwr_curv_dfa/coord_list.csv"

cat("Loading existing coordinate list\n")
coord_list <- read.csv(coord_list_fname)

ref_catalog_2025 <- fread(
  file.path("pwr_curv_dfa/ref_catalog_wind_2025.csv.gz")
)


national <- read_csv(
  'tables/lwe_eda/national_hourly_generation.csv',
  show_col_types = FALSE
) %>%
  mutate(time = as.POSIXct(time, tz = 'UTC')) %>%
  arrange(time)


national <- national %>%
  rename(
    capacity = capacity_MW,
    potential = potential_MW,
    power_est0 = power_est0_MW
  ) %>%
  mutate(
    norm_potential = if_else(capacity > 0, potential / capacity, NA_real_),
    norm_power_est0 = if_else(capacity > 0, power_est0 / capacity, NA_real_),
    error_mw = potential - power_est0,
    error0 = norm_potential - norm_power_est0,
    year = year(time),
    month = factor(month(time), levels = 1:12, labels = month.abb),
    month_num = month(time),
    hour = factor(hour(time), levels = 0:23),
    season = factor(
      case_when(
        month_num %in% c(12, 1, 2) ~ 'Winter',
        month_num %in% 3:5 ~ 'Spring',
        month_num %in% 6:8 ~ 'Summer',
        TRUE ~ 'Autumn'
      ),
      levels = c('Winter', 'Spring', 'Summer', 'Autumn')
    ),
    gap1 = as.numeric(difftime(time, lag(time), units = 'hours')),
    gap2 = as.numeric(difftime(time, lag(time, 2), units = 'hours')),
    gap3 = as.numeric(difftime(time, lag(time, 3), units = 'hours')),
    dest = if_else(gap1 == 1, norm_power_est0 - lag(norm_power_est0), NA_real_),
    dest_l1 = if_else(gap2 == 2, lag(dest), NA_real_),
    dest_l2 = if_else(gap3 == 3, lag(dest, 2), NA_real_),
    err_l1 = if_else(gap1 == 1, lag(error0), NA_real_),
    err_l2 = if_else(gap2 == 2, lag(error0, 2), NA_real_)
  )


# ---- 3. Descriptive summaries ----
error_summary <- function(data) {
  data %>%
    summarise(
      n_hours = sum(!is.na(error0)),
      bias_cf = mean(error0, na.rm = TRUE),
      mae_cf = mean(abs(error0), na.rm = TRUE),
      rmse_cf = sqrt(mean(error0^2, na.rm = TRUE)),
      sd_cf = sd(error0, na.rm = TRUE),
      q01_cf = as.numeric(quantile(error0, .01, na.rm = TRUE)),
      q50_cf = median(error0, na.rm = TRUE),
      q99_cf = as.numeric(quantile(error0, .99, na.rm = TRUE)),
      .groups = 'drop'
    )
}

annual_bias <- national %>% group_by(year) %>% error_summary()
monthly_bias <- national %>% group_by(month_num, month) %>% error_summary()
seasonal_bias <- national %>% group_by(season) %>% error_summary()
year_month_bias <- national %>% group_by(year, month_num) %>% error_summary()
hourly_bias <- national %>% group_by(hour) %>% error_summary()


print(annual_bias)
print(seasonal_bias)

# scatter plots norm_potential vs norm_power_est0
p_scatter <- ggplot(national, aes(norm_power_est0, norm_potential)) +
  geom_point(alpha = 0.5) +
  labs(
    title = 'Scatter plot: normalised potential vs normalised initial estimate',
    x = 'Normalised initial estimate',
    y = 'Normalised potential'
  ) +
  theme_minimal()

print(p_scatter)

# ---- 4. Initial-error plots ----
p_error_hist <- ggplot(national, aes(error0)) +
  geom_histogram(bins = 80, fill = blues9[6]) +
  # geom_vline(xintercept = 0, linetype = 'dashed') +
  labs(
    title = 'Initial national estimate: normalised error',
    x = 'Potential minus initial estimate / national capacity',
    y = 'Hours'
  ) +
  theme_minimal()

p_annual_bias <- ggplot(annual_bias, aes(year, bias_cf)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(
    title = 'Annual mean initial-estimate bias',
    x = 'Year',
    y = 'Mean normalised error'
  ) +
  theme_minimal()

p_monthly_bias <- ggplot(monthly_bias, aes(month_num, bias_cf)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    title = 'Mean initial-estimate bias by month',
    x = 'Month',
    y = 'Mean normalised error'
  ) +
  theme_minimal()

p_year_month <- ggplot(
  year_month_bias,
  aes(month_num, factor(year), fill = bias_cf)
) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    title = 'Monthly bias by year',
    x = 'Month',
    y = 'Year',
    fill = 'Mean error'
  ) +
  theme_minimal()

p_season_box <- ggplot(national, aes(season, error0)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = quantile(national$error0, c(.01, .99), na.rm = TRUE)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(
    title = 'Initial-estimate error by season (1st–99th percentile view)',
    x = 'Season',
    y = 'Normalised error'
  ) +
  theme_minimal()

p_hourly_bias <- ggplot(
  hourly_bias,
  aes(as.integer(as.character(hour)), bias_cf)
) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(
    title = 'Bias by UTC hour',
    x = 'UTC hour',
    y = 'Mean normalised error'
  ) +
  theme_minimal()


# Bin by initial estimated CF; avoids a very dense 7-year scatter plot.
by_estimate <- national %>%
  filter(is.finite(norm_power_est0), is.finite(error0)) %>%
  mutate(
    est_bin = cut(
      norm_power_est0,
      breaks = unique(quantile(norm_power_est0, seq(0, 1, .05))),
      include.lowest = TRUE
    )
  ) %>%
  group_by(est_bin) %>%
  summarise(
    est_cf = mean(norm_power_est0),
    mean_error = mean(error0),
    q10 = quantile(error0, .1),
    q90 = quantile(error0, .9),
    .groups = 'drop'
  )

p_bias_vs_est <- ggplot(by_estimate, aes(est_cf, mean_error)) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = .2) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(
    title = 'Bias conditional on initial estimated capacity factor',
    x = 'Initial estimated capacity factor',
    y = 'Normalised error'
  ) +
  theme_minimal()


# ---- 5. Fit national LM without technology type ----
# In-sample model, then chronological test to assess generalisation.
if (!'ws_h_wmean' %in% names(national)) {
  stop(
    'ws_h_wmean is missing. Add it to the CSV or use the optional farm-level join above.'
  )
}

lm_data <- national %>%
  filter(if_all(
    all_of(c(
      'norm_potential',
      'norm_power_est0',
      'dest',
      'dest_l1',
      'dest_l2',
      'ws_h_wmean'
    )),
    is.finite
  )) %>%
  mutate(
    month = factor(month, levels = levels(national$month)),
    hour = factor(hour, levels = levels(national$hour))
  )

# The model has no tech_typ, so interactions with tech_typ are removed.
# Polynomial basis is learned on the training set and used for prediction.
lm_formula <- norm_potential ~ norm_power_est0 +
  month +
  hour +
  dest +
  dest_l1 +
  dest_l2 +
  poly(ws_h_wmean, 3)

# Last 20% of time as holdout (not a random split).
split_time <- quantile(as.numeric(lm_data$time), .8)
train <- lm_data %>% filter(as.numeric(time) <= split_time)
test <- lm_data %>% filter(as.numeric(time) > split_time)
if (nrow(train) < 100 || nrow(test) < 24) {
  stop('Insufficient train/test rows.')
}

lm_fit <- lm(lm_formula, data = train)
print(summary(lm_fit))

train <- train %>%
  mutate(
    fitted_cf = as.numeric(predict(lm_fit, newdata = train)),
    residual_cf = norm_potential - fitted_cf,
    subset = 'Train'
  )
test <- test %>%
  mutate(
    fitted_cf = as.numeric(predict(lm_fit, newdata = test)),
    residual_cf = norm_potential - fitted_cf,
    subset = 'Test'
  )
model_diagnostics <- bind_rows(train, test)

# Note: these are prediction errors on the test set, not in-sample LM residuals.
score_comparison <- model_diagnostics %>%
  group_by(subset) %>%
  summarise(
    n = n(),
    initial_bias = mean(error0),
    corrected_bias = mean(residual_cf),
    initial_mae = mean(abs(error0)),
    corrected_mae = mean(abs(residual_cf)),
    initial_rmse = sqrt(mean(error0^2)),
    corrected_rmse = sqrt(mean(residual_cf^2)),
    .groups = 'drop'
  )
print(score_comparison)


# ---- 6. Residual and predictive-error plots ----
p_residual_hist <- ggplot(model_diagnostics, aes(residual_cf, fill = subset)) +
  geom_histogram(bins = 70, alpha = .55, position = 'identity') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  facet_wrap(~subset, scales = 'free_y') +
  labs(
    title = 'LM calibration error distribution',
    x = 'Observed minus fitted normalised potential',
    y = 'Hours'
  ) +
  theme_minimal() +
  guides(fill = 'none')

p_residual_vs_fitted <- ggplot(model_diagnostics, aes(fitted_cf, residual_cf)) +
  geom_point(alpha = .08, size = .4) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  facet_wrap(~subset) +
  labs(
    title = 'Calibration errors versus fitted capacity factor',
    x = 'Fitted capacity factor',
    y = 'Observed minus fitted'
  ) +
  theme_minimal()

residual_by_month <- model_diagnostics %>%
  group_by(subset, month_num) %>%
  summarise(
    bias = mean(residual_cf),
    rmse = sqrt(mean(residual_cf^2)),
    .groups = 'drop'
  )
p_residual_month <- ggplot(
  residual_by_month,
  aes(month_num, bias, colour = subset)
) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    title = 'Residual bias by month',
    x = 'Month',
    y = 'Mean calibration error'
  ) +
  theme_minimal()

p_qq <- ggplot(model_diagnostics, aes(sample = residual_cf)) +
  stat_qq(alpha = .15) +
  stat_qq_line() +
  facet_wrap(~subset) +
  labs(title = 'Normal Q–Q plots of calibration errors') +
  theme_minimal()

# ACF: use the longest uninterrupted hourly run in each subset.
# This avoids interpreting missing hours as adjacent lags.
plot_longest_run_acf <- function(dat, variable, main, max_lag = 168) {
  dat <- dat %>%
    arrange(time) %>%
    mutate(
      run_id = cumsum(
        is.na(lag(time)) |
          as.numeric(difftime(time, lag(time), units = 'hours')) != 1
      )
    )
  run_lengths <- dat %>% count(run_id, sort = TRUE)
  selected <- dat %>% filter(run_id == run_lengths$run_id[1])
  if (nrow(selected) > 3) {
    stats::acf(
      selected[[variable]],
      lag.max = min(max_lag, nrow(selected) - 1),
      na.action = na.omit,
      main = main
    )
  }
}

# PACF: use the longest uninterrupted hourly run in each subset.
plot_longest_run_pacf <- function(dat, variable, main, max_lag = 168) {
  dat <- dat %>%
    arrange(time) %>%
    mutate(
      run_id = cumsum(
        is.na(lag(time)) |
          as.numeric(difftime(time, lag(time), units = 'hours')) != 1
      )
    )
  run_lengths <- dat %>% count(run_id, sort = TRUE)
  selected <- dat %>% filter(run_id == run_lengths$run_id[1])
  if (nrow(selected) > 3) {
    stats::pacf(
      selected[[variable]],
      lag.max = min(max_lag, nrow(selected) - 1),
      na.action = na.omit,
      main = main
    )
  }
}

# ---- 7. Display ----
plots <- list(
  p_error_hist,
  p_annual_bias,
  p_monthly_bias,
  p_year_month,
  p_season_box,
  p_hourly_bias,
  p_bias_vs_est,
  p_residual_hist,
  p_residual_vs_fitted,
  p_residual_month,
  p_qq
)
names(plots) <- c(
  "Initial error histogram",
  "Annual bias",
  "Monthly bias",
  "Year-month bias",
  "Seasonal boxplot",
  "Hourly bias",
  "Bias vs initial estimate",
  "Residual histogram",
  "Residual vs fitted",
  "Residual by month",
  "Q-Q plot"
)
bias_output_figpath <- gsub("lwe_eda", "bias_eda", lwe_output_figpath)
for (i in seq_along(plots)) {
  p <- plots[[i]]
  cat("Saving plot:", names(plots)[i], "\n")
  print(p)
  # browser()
  ggsave(
    filename = file.path(
      bias_output_figpath,
      paste0(names(plots)[i], ".png")
    ),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300,
    create.dir = TRUE
  )
}
plot_longest_run_acf(national, 'error0', 'Initial normalised error ACF')
plot_longest_run_acf(train, 'residual_cf', 'Training LM residual ACF')
plot_longest_run_acf(test, 'residual_cf', 'Test LM prediction-error ACF')
plot_longest_run_pacf(national, 'error0', 'Initial normalised error PACF')
plot_longest_run_pacf(train, 'residual_cf', 'Training LM residual PACF')
plot_longest_run_pacf(test, 'residual_cf', 'Test LM prediction-error PACF')
