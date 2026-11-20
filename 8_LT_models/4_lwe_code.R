local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

pow_threshold <- 0.1
pow_threshold_label <- gsub("\\.", "_", as.character(pow_threshold))
tol <- 0.01
norm_dist_tol <- 1

override_objects <- TRUE
# rerun_samples <- FALSE
# prec_init <- log(200)
# batch_name <- "batch2025"
# options
n.days.before.fit <- 365 * 2
n.days.after.fit <- 92
n_samp <- 100
d0 <- as.Date("2025-06-01")
d0_tag <- base::format(d0, "%y%m%d")

batch_name <- "batchLT"
batch_name <- sprintf(
  "%s_%d%s_%d%s_%ds",
  batch_name,
  round(
    if (n.days.before.fit >= 365) n.days.before.fit / 365 else n.days.before.fit
  ),
  if (n.days.before.fit >= 365) "Y" else "D",
  round(
    if (n.days.after.fit >= 365) n.days.after.fit / 365 else n.days.after.fit
  ),
  if (n.days.after.fit >= 365) "y" else "d",
  n_samp
)
# batch_name <- "batchLT_test"

if (local_run) {
  cat("Running in local mode\n")
} else {
  cat("Running in cluster mode\n")
}

# 0.2 libraries and paths ####
require(parallel)

if (local_run) {
  data_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
  extra_path <- "~/Documents/elexon/extra"
  output_path <- "~/Documents/elexon/caloutput"
  sample_path <- "~/Documents/elexon/samples"
  # n_samp <- 100
  pixel_dims <- c(150, 150)
} else {
  data_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  gen_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  model_path <- "/exports/eddie/scratch/s2441782/calibration/model_objects"
  extra_path <- "/exports/eddie/scratch/s2441782/calibration/extra"
  output_path <- "/exports/eddie/scratch/s2441782/calibration/caloutput"
  sample_path <- "/exports/eddie/scratch/s2441782/calibration/samples"
  temp_lib <- "/exports/eddie3_homes_local/s2441782/lib"
  pixel_dims <- c(300, 300)
  # n_samp <- 1000
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
# sampled_days_df <- read.csv("data/sample_days_df_25_n150.csv") %>%
#   mutate(date = as.Date(date))

# sampled_days_df <- read.csv("data/sample_days_df.csv") %>%
#   mutate(date = as.Date(date))

# sampled_days <- sampled_days_df %>%
#   pull(date) %>%
#   sort()
sampled_days <- c(d0)
coord_list_fname <- "data/coord_list.csv"

cat("Loading existing coordinate list\n")
coord_list <- read.csv(coord_list_fname)
coord_list <- coord_list %>%
  mutate(
    sampled_old = sampled,
    sampled = TRUE
  )
# Read models ####
model_catalog <- read.csv("data/model_catalog.csv") %>%
  na.omit()

# if (local_run) {
#   model_catalog <- model_catalog %>%
#     # filter(!grepl("fine", mod_labels)) %>%
#     filter(!grepl("st0", est_cols))
# }
mod_labels <- model_catalog$mod_labels
est_cols <- model_catalog$est_cols
n_models <- length(est_cols)
names(mod_labels) <- est_cols
# excluded_models0 <- c("lm_bru")
# excluded_models <- c("lm_bru", "qm")

excluded_models0 <- c("lm", "", "spde1d", "lm_t")
excluded_models <- c("lm", "", "spde1d", "lm_t", "qm")
# mod_labels["lm_bru"] <- "Linear Model"

model_df <- model_catalog %>%
  rename(code = est_cols, label = mod_labels) %>%
  arrange(desc(nchar(mode_code_prefix))) %>%
  mutate(
    type = case_when(
      grepl("lm_", code) ~ "bru",
      grepl("lm", code) ~ "lm",
      grepl("qm", code) ~ "qm",
      TRUE ~ "bru"
    )
  ) %>%
  filter(!code %in% c("st0_m0", "st0_m1"))
bru_df <- model_df %>% filter(type == "bru")

# LWE diagram ####

cat(
  "--------------------------------------------------------------------\n"
)
cat(
  " Low wind events analysis for days:",
  paste(format(sampled_days, "%b-%d"), collapse = ", "),
  ")\n"
)

cat("override_objects:", override_objects, "\n")

cat(
  "--------------------------------------------------------------------\n"
)


extension <- "rds"
model_df0 <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days[i] %>% as.Date()
    # print(d0)
    d0_tag <- base::format(d0, "%y%m%d")
    file_name <- sprintf(
      "%s/%s/summaries/oos/GB_fig_band_summary_time_%s.rds",
      output_path,
      batch_name,
      d0_tag
    )
    if (!file.exists(file_name)) {
      cat("File not found:", file_name, "\n")
      return(NULL)
    }
    # coarse
    readRDS(file_name) %>%
      st_drop_geometry()
  }
) %>%
  bind_rows() %>%
  filter(!model %in% excluded_models0)

# pos_breaks <- with(
#   model_df0,
#   quantile(elevation[elevation > 0], probs = seq(0, 1, 1 / 3))
# )
# pos_levels <- levels(cut(
#   model_df0$elevation[model_df0$elevation > 0],
#   breaks = pos_breaks,
#   include.lowest = TRUE
# ))

# GB aggregate model df
#
# LWE stats after aggregation
# LWE stats after aggregation
## 2.0 Low wind events in aggregated observed data ####
library(data.table)

cat("--------------------------------------------------------------------\n")
cat("Low wind events in observed GB data\n")
cat("--------------------------------------------------------------------\n")

max_h_duration <- 24 * 10
rolling_window <- 12

hist_ref_date <- as.Date(sampled_days)
stopifnot(length(hist_ref_date) == 1L, !is.na(hist_ref_date))

hist_start_date <- lubridate::`%m-%`(hist_ref_date, lubridate::years(7))
hist_end_date <- hist_ref_date + n.days.after.fit # exclusive

hist_start <- as.POSIXct(hist_start_date, tz = "UTC")
hist_end <- as.POSIXct(hist_end_date, tz = "UTC")

lwe_obs_hist_fname <- file.path(
  "summaries",
  sprintf(
    "lwe_obs_hist_GB_%s_%s_roll%dh_t%s_tol%s_cap%sh.%s",
    format(hist_start_date, "%Y%m%d"),
    format(hist_end_date, "%Y%m%d"),
    rolling_window,
    pow_threshold_label,
    norm_dist_tol,
    max_h_duration,
    extension
  )
)

if (!file.exists(lwe_obs_hist_fname) || override_objects) {
  sampled_coords <- coord_list$coord_id[which(coord_list$sampled)]
  cols <- c(
    "halfHourEndTime",
    "coord_id",
    "lon",
    "lat",
    "potential",
    "power_est0",
    "capacity"
  )

  cat("Reading and filtering observed wind-farm data\n")
  obs <- arrow::read_parquet(
    file.path(gen_path, "power_curve_all_enriched.parquet"),
    col_select = tidyselect::all_of(cols)
  )
  setDT(obs)
  setnames(obs, "halfHourEndTime", "time")

  obs <- obs[
    time >= hist_start & time < hist_end & coord_id %in% sampled_coords,
    .(lon, lat, time, potential, power_est0, capacity)
  ]
  if (!nrow(obs)) {
    stop(
      "No observations found for the historical window and sampled locations."
    )
  }

  # Combine co-located units before applying clipping and anomaly detection.
  site_hist <- obs[,
    .(
      potential = sum(potential),
      power_est0 = sum(power_est0),
      capacity = sum(capacity)
    ),
    by = .(lon, lat, time)
  ]
  rm(obs)

  site_hist[, `:=`(
    norm_potential = pmin(1, potential / capacity),
    norm_power_est0 = power_est0 / capacity
  )]
  site_hist[, anomaly := abs(norm_power_est0 - norm_potential) >= norm_dist_tol]
  site_hist[is.na(anomaly), anomaly := FALSE]

  cat(sprintf(
    "Percentage of location-hours flagged at %.1f%% tol: %.2f%%\n",
    norm_dist_tol * 100,
    100 * mean(site_hist$anomaly)
  ))

  # Use the same eligible locations in the numerator and denominator.
  # Keep hours with no eligible locations so they can break rolling windows.
  site_hist[,
    valid := !anomaly &
      is.finite(norm_potential) &
      is.finite(capacity) &
      capacity > 0
  ]
  site_hist[, `:=`(
    potential_used = fifelse(valid, norm_potential * capacity, 0),
    capacity_used = fifelse(valid, capacity, 0)
  )]

  gb_hist <- site_hist[,
    .(
      potential = sum(potential_used),
      capacity = sum(capacity_used)
    ),
    by = time
  ]
  rm(site_hist)

  gb_hist[, gb_cf := NA_real_]
  gb_hist[capacity > 0, gb_cf := potential / capacity]
  setorder(gb_hist, time)

  # Start a new block whenever there is a gap in the hourly series.
  gb_hist[,
    block_id := cumsum(
      c(TRUE, diff(as.numeric(time)) != 3600)
    )
  ]

  gb_hist[,
    gb_cf_roll := {
      if (.N < rolling_window) {
        NA_real_
      } else {
        frollmean(gb_cf, n = rolling_window, align = "right", na.rm = FALSE)
      }
    },
    by = block_id
  ]

  gb_hist[, below := !is.na(gb_cf_roll) & gb_cf_roll < pow_threshold]
  gb_hist[, run_id := rleid(below), by = block_id]

  lwe_obs_hist <- gb_hist[
    below == TRUE,
    .(
      start_time = first(time),
      duration_h = min(max_h_duration, .N)
    ),
    by = .(block_id, run_id)
  ]
  lwe_obs_hist[, c("block_id", "run_id") := NULL]
  lwe_obs_hist[, model := "observed"]

  dir.create(
    dirname(lwe_obs_hist_fname),
    recursive = TRUE,
    showWarnings = FALSE
  )
  saveRDS(lwe_obs_hist, lwe_obs_hist_fname)
} else {
  cat("Loading existing low wind events in observed GB data\n")
  lwe_obs_hist <- as.data.table(readRDS(lwe_obs_hist_fname))
}


## 2.1 Low wind events in observed WF data ####
cat("--------------------------------------------------------------------\n")
cat("Low wind events in observed data\n")
cat("--------------------------------------------------------------------\n")
max_h_duration <- 72
lwe_obs_pred_fname <- file.path(
  "summaries",
  sprintf(
    "lwe_obs_pred_%s_t%s.%s",
    "sampled_days",
    pow_threshold_label,
    extension
  )
)
lwe_obs_newloc_fname <- file.path(
  "summaries",
  sprintf(
    "lwe_obs_newloc_%s_t%s.%s",
    "sampled_days",
    pow_threshold_label,
    extension
  )
)

library(data.table)

cache_exists <- file.exists(lwe_obs_pred_fname)

if (!cache_exists || override_objects) {
  if (cache_exists) {
    cat(
      "LWE processed file found, but override_objects is TRUE. Recreating summary\n"
    )
  } else {
    cat("LWE processed file not found, creating new summary\n")
  }

  days_before <- 7L
  days_after <- n.days.after.fit
  pred_days <- unique(as.vector(outer(
    as.integer(as.Date(sampled_days)),
    seq.int(-days_before, days_after - 1L),
    "+"
  )))
  sampled_coords <- coord_list$coord_id[which(coord_list$sampled)]

  cols <- c(
    "halfHourEndTime",
    "coord_id",
    "lon",
    "lat",
    "site_name",
    "potential",
    "power_est0",
    "capacity"
  )

  cat("Reading and filtering observed wind-farm data\n")
  obs <- arrow::read_parquet(
    file.path(gen_path, "power_curve_all_enriched.parquet"),
    col_select = tidyselect::all_of(cols)
  )
  setDT(obs)
  setnames(obs, "halfHourEndTime", "time")

  # Filter coordinates first so date conversion runs on fewer rows.
  obs <- obs[coord_id %in% sampled_coords]
  obs <- obs[as.integer(as.IDate(time, tz = "UTC")) %in% pred_days]

  obs[,
    site_name := trimws(gsub(
      "\\b(wind\\s*farm|wf)\\b",
      "",
      site_name,
      ignore.case = TRUE
    ))
  ]

  # Reproduce the original coord_id reassignment after filtering.
  obs[, coord_id := as.integer(factor(paste(lon, lat)))]

  # Keep the first cleaned site name in alphabetical order within each group.
  setorderv(obs, "site_name", na.last = TRUE)
  hourly <- obs[,
    .(
      site_name = first(site_name),
      coord_id = first(coord_id),
      potential = sum(potential),
      power_est0 = sum(power_est0),
      capacity = sum(capacity)
    ),
    by = .(lon, lat, time)
  ]
  rm(obs)

  hourly[, `:=`(
    norm_potential = pmin(1, potential / capacity),
    norm_power_est0 = power_est0 / capacity
  )]

  cat("Identifying anomalies in the dataset\n")
  hourly[, anomaly := abs(norm_power_est0 - norm_potential) >= norm_dist_tol]
  hourly[is.na(anomaly), anomaly := FALSE]

  cat(sprintf(
    "Percentage of anomalies at %.1f%% tol in the dataset: %.2f%%\n",
    norm_dist_tol * 100,
    100 * mean(hourly$anomaly)
  ))

  # Masked anomalies and missing observations both break a low-wind run.
  hourly[,
    below := !anomaly & !is.na(norm_potential) & norm_potential < pow_threshold
  ]

  hourly <- hourly[, .(coord_id, site_name, time, below)]
  setorderv(hourly, c("coord_id", "site_name", "time"), na.last = TRUE)
  hourly[, run_id := rleid(below), by = .(coord_id, site_name)]

  # Keep all rows until run IDs are assigned so intervening FALSE rows split events.
  lwe_obs_pred <- hourly[
    below == TRUE,
    .(
      start_time = first(time),
      end_time = last(time)
    ),
    by = .(coord_id, site_name, run_id)
  ]

  lwe_obs_pred[, `:=`(
    duration_h = pmin(
      max_h_duration,
      as.numeric(difftime(end_time, start_time, units = "hours")) + 1
    ),
    model = "observed"
  )]
  lwe_obs_pred[, c("site_name", "run_id", "end_time") := NULL]

  saveRDS(lwe_obs_pred, lwe_obs_pred_fname)
  rm(hourly)
  gc()
} else {
  cat("Loading existing low wind events observed data\n")
  lwe_obs_pred <- as.data.table(readRDS(lwe_obs_pred_fname))
}

## 3. Low wind events in model predictions ####
cat("--------------------------------------------------------------------\n")
cat("Low wind events in model mean predictions\n")
cat("--------------------------------------------------------------------\n")

low_events_model <- model_df0 %>%
  filter(model %in% est_cols) %>%
  st_drop_geometry() %>%
  # filter(model == "Spatio-temporal model") %>%
  # filter(coord_id %in% coord_list$coord_id[coord_list$sampled]) %>%
  arrange(model, time) %>%
  group_by(model) %>%
  mutate(
    below = mean < pow_threshold,
    run_id = cumsum(below != lag(below, default = first(below)))
  ) %>%
  filter(below) %>%
  group_by(model, run_id) %>%
  summarise(
    start_time = first(time),
    end_time = last(time),
    duration_h = as.numeric(difftime(end_time, start_time, units = "hours")) +
      1,
    .groups = "drop"
  ) %>%
  dplyr::select(model, start_time, duration_h) %>%
  bind_rows(lwe_obs_pred) %>%
  mutate(
    model = factor(
      model,
      levels = c(est_cols, "observed"),
      labels = c(mod_labels, "observed" = "Observed")
    )
  )

rm(model_df0)
gc()

# low_events_model$model %>% unique() %>% sort() %>% print()
low_events_model %>%
  filter(!model %in% mod_labels[excluded_models0]) %>%
  filter(duration_h < max_h_duration) %>%
  # filter(
  #   model %in%
  #     c("observed", "Linear model", "AR1 model", "QM", "Spatio-temporal model")
  # ) %>%
  ggplot(aes(x = duration_h)) +
  # geom_density(aes(fill = model), alpha = 0.5) +
  geom_histogram(
    aes(fill = model),
    position = "identity",
    alpha = 0.5,
    binwidth = 1
  ) +
  # scale_x_log10(n.breaks = 8) +
  labs(
    title = "Distribution of Low Wind Events Duration",
    x = "Duration (hours)",
    y = "Frequency"
  ) +
  facet_wrap(
    ~model,
    scales = "free",
    # labeller = as_labeller(mod_labels)
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_npg()

ggsave(
  sprintf(
    "fig/%s/low_wind_duration_dist_t%s.pdf",
    batch_name,
    pow_threshold_label
  ),
  width = 10,
  height = 6
)

mods <- unique(low_events_model$model)
# Default hue palette
cols <- scales::hue_pal()(length(mods))
names(cols) <- mods

cols["Observed"] <- "darkred"
# excluded_models0
# mod_labels
# density plot of low wind event durations
# low_events_model %>%
#   filter(!model %in% mod_labels[excluded_models0]) %>%
#   filter(duration_h < max_h_duration) %>%
#   ggplot(aes(x = duration_h)) +
#   geom_density(
#     # data = ~ dplyr::filter(.x, model != "observed"),
#     aes(colour = model),
#     alpha = 0.5,
#     key_glyph = "path",
#     lwd = 1
#   ) +
#   # geom_density(
#   #   data = ~ dplyr::filter(.x, model == "observed"),
#   #   aes(colour = model),
#   #   # colour = "darkred",
#   #   # linewidth = 1.2,
#   #   fill = NA
#   # ) +
#   labs(
#     title = "Distribution of Low Wind Events Duration",
#     x = "Duration (hours)",
#     y = "Frequency"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "inside",
#     legend.position.inside = c(0.7, 0.7)
#   ) +
#   guides(
#     colour = guide_legend(ncol = 2)
#   ) +
#   scale_color_manual(values = cols)

# ggsave(
#   sprintf(
#     "fig/%s/low_wind_duration_dens_t%s.pdf",
#     batch_name,
#     pow_threshold_label
#   ),
#   width = 6,
#   height = 4
# )

# CDF plot of low wind event durations
low_events_model %>%
  filter(!model %in% mod_labels[excluded_models0]) %>%
  # filter(model %in% c("Observed", "QM")) %>%
  filter(duration_h < max_h_duration) %>%
  ggplot(aes(x = duration_h)) +
  stat_ecdf(
    aes(colour = model),
    alpha = 0.5,
    key_glyph = "path",
    lwd = 1
  ) +
  labs(
    title = "Cumulative Distribution of Low Wind Events Duration",
    x = "Duration (hours)",
    y = "Cumulative Probability"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "inside",
    legend.position.inside = c(0.7, 0.3)
  ) +
  guides(
    colour = guide_legend(ncol = 2)
  ) +
  scale_color_manual(values = cols)
ggsave(
  sprintf(
    "fig/%s/low_wind_duration_cdf_t%s.pdf",
    batch_name,
    pow_threshold_label
  ),
  width = 6,
  height = 4
)
obs <- low_events_model %>%
  filter(
    model == "Observed",
    duration_h < max_h_duration
  ) %>%
  pull(duration_h)

probs <- seq(0, 1, length.out = 101)

qq_df <- low_events_model %>%
  filter(model != "Observed") %>%
  filter(!model %in% mod_labels[excluded_models0]) %>%
  filter(duration_h < max_h_duration) %>%
  group_by(model) %>%
  summarise(
    obs_q = list(quantile(obs, probs = probs, na.rm = TRUE)),
    model_q = list(quantile(duration_h, probs = probs, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  unnest(c(obs_q, model_q))

ggplot(qq_df, aes(x = obs_q, y = model_q)) +
  geom_point(size = 1) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  facet_wrap(
    ~model,
    scales = "free",
    #  labeller = as_labeller(mod_labels)
  ) +
  labs(
    title = "Q-Q Plot of Low Wind Event Durations",
    x = "Observed quantiles",
    y = "Model quantiles"
  ) +
  # coord_equal() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  sprintf(
    "fig/%s/low_wind_duration_qq_t%s.pdf",
    batch_name,
    pow_threshold_label
  ),
  width = 10,
  height = 6
)

ggplot(qq_df, aes(x = obs_q, y = model_q)) +
  geom_line(aes(col = model)) +
  geom_abline(slope = 1, intercept = 0, colour = "darkred") +
  # facet_wrap(
  #   ~model,
  #   scales = "free",
  #   #  labeller = as_labeller(mod_labels)
  # ) +
  labs(
    title = "Q-Q Plot of Low Wind Event Durations",
    x = "Observed quantiles",
    y = "Model quantiles"
  ) +
  # coord_equal() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_manual(
    values = cols
  ) +
  coord_equal(ratio = 1)

ggsave(
  sprintf(
    "fig/%s/low_wind_duration_qq_line_t%s.pdf",
    batch_name,
    pow_threshold_label
  ),
  width = 6,
  height = 5
)
