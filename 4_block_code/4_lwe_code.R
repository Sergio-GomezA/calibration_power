local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

pow_threshold <- 0.05
pow_threshold_label <- gsub("\\.", "_", as.character(pow_threshold))


override_objects <- FALSE
# rerun_samples <- FALSE
# prec_init <- log(200)
# batch_name <- "batch2025"
batch_name <- "batchY25d150"


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
  n_samp <- 100
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
sampled_days_df <- read.csv("data/sample_days_df_25_n150.csv") %>%
  mutate(date = as.Date(date))

# sampled_days_df <- read.csv("data/sample_days_df.csv") %>%
#   mutate(date = as.Date(date))

sampled_days <- sampled_days_df %>%
  pull(date) %>%
  sort()

coord_list_fname <- "data/coord_list.csv"

cat("Loading existing coordinate list\n")
coord_list <- read.csv(coord_list_fname)

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

excluded_models0 <- c("")
excluded_models <- c("", "qm")

model_catalog <- read.csv("data/model_catalog.csv") %>%
  na.omit()
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
  paste(format(sampled_days, "%Y-%m-%d"), collapse = ", "),
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
      "%s/%s/summaries/oos/WF_fig_band_summary_time_%s.rds",
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
  bind_rows()

# pos_breaks <- with(
#   model_df0,
#   quantile(elevation[elevation > 0], probs = seq(0, 1, 1 / 3))
# )
# pos_levels <- levels(cut(
#   model_df0$elevation[model_df0$elevation > 0],
#   breaks = pos_breaks,
#   include.lowest = TRUE
# ))

## 2.1 Low wind events in observed data ####
cat("--------------------------------------------------------------------\n")
cat("Low wind events in observed data\n")
cat("--------------------------------------------------------------------\n")
max_h_duration <- 100
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

if (!file.exists(lwe_obs_pred_fname) | override_objects) {
  if (!file.exists(lwe_obs_pred_fname)) {
    cat("LWE processed file not found, creating new summary\n")
  } else {
    cat(
      "LWE processed file found, but override_objects is TRUE. Recreating summary\n"
    )
  }

  pwr_curv_df <- read_parquet(file.path(
    gen_path,
    "power_curve_all_enriched.parquet"
  ))

  cat("This may take a while...\n")

  # d0 <- as.Date("2024-01-01")
  # n.days <- 365

  days_before <- 0
  days_after <- 3
  pred_days <- lapply(
    sampled_days %>% as.Date(),
    \(d0) {
      seq(d0 - days_before, d0 + days_after - 1, by = "days") %>%
        as.Date()
    }
  ) %>%
    unlist() %>%
    unique() %>%
    sort() %>%
    as.Date()

  lwe_obs_pred <- pwr_curv_df %>%
    rename(time = halfHourEndTime) %>%
    mutate(
      date = as.Date(time)
    ) %>%
    # filter(date >= d0, date <= d0 + n.days) %>%
    filter(date %in% pred_days) %>%
    filter(coord_id %in% coord_list$coord_id[coord_list$sampled]) %>%
    mutate(
      elevation = pmax(0, elevation),
      site_name = site_name %>%
        gsub("\\b(wind\\s*farm|wf)\\b", "", ., ignore.case = TRUE) %>%
        trimws()
    ) %>%
    # filter(date %in% sampled_days) %>%
    # filter(date >= d0, date <= d0 + n.days - 1) %>%
    arrange(site_name) %>%
    mutate(
      site_id = as.integer(factor(site_name)),
      coord_id = as.integer(factor(paste(lon, lat)))
    ) %>%
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
    arrange(coord_id, time) %>%
    group_by(coord_id, site_name) %>%
    mutate(
      below = norm_potential < pow_threshold,
      run_id = cumsum(below != lag(below, default = first(below)))
    ) %>%
    group_by(coord_id, site_name, run_id, below) %>%
    summarise(
      start_time = first(time),
      end_time = last(time),
      duration_h = pmin(
        max_h_duration,
        as.numeric(difftime(end_time, start_time, units = "hours")) +
          1
      ),
      .groups = "drop"
    ) %>%
    filter(below) %>%
    dplyr::select(coord_id, start_time, duration_h) %>%
    mutate(
      model = "observed"
    )

  # low_events %>% pull(duration_h) %>% summary()

  saveRDS(lwe_obs_pred, lwe_obs_pred_fname)

  n.days <- 0
  n.days.before <- 7
  n.hours <- 1

  days_newloc <- lapply(
    sampled_days %>% as.Date(),
    \(d0) {
      seq(d0 - n.days.before, d0 + n.days - 1, by = "days") %>%
        as.Date()
    }
  ) %>%
    unlist() %>%
    unique() %>%
    sort() %>%
    as.Date()

  lwe_obs_newloc <- pwr_curv_df %>%
    rename(time = halfHourEndTime) %>%
    mutate(
      date = as.Date(time)
    ) %>%
    # filter(date >= d0, date <= d0 + n.days) %>%
    filter(date %in% days_newloc) %>%
    filter(coord_id %in% coord_list$coord_id[!coord_list$sampled]) %>%
    mutate(
      elevation = pmax(0, elevation),
      site_name = site_name %>%
        gsub("\\b(wind\\s*farm|wf)\\b", "", ., ignore.case = TRUE) %>%
        trimws()
    ) %>%
    # filter(date %in% sampled_days) %>%
    # filter(date >= d0, date <= d0 + n.days - 1) %>%
    arrange(site_name) %>%
    mutate(
      site_id = as.integer(factor(site_name)),
      coord_id = as.integer(factor(paste(lon, lat)))
    ) %>%
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
    arrange(coord_id, time) %>%
    group_by(coord_id, site_name) %>%
    mutate(
      below = norm_potential < pow_threshold,
      run_id = cumsum(below != lag(below, default = first(below)))
    ) %>%
    group_by(coord_id, site_name, run_id, below) %>%
    summarise(
      start_time = first(time),
      end_time = last(time),
      duration_h = pmin(
        max_h_duration,
        as.numeric(difftime(end_time, start_time, units = "hours")) +
          1
      ),
      .groups = "drop"
    ) %>%
    filter(below) %>%
    dplyr::select(coord_id, start_time, duration_h) %>%
    mutate(
      model = "observed"
    )

  saveRDS(lwe_obs_newloc, lwe_obs_newloc_fname)
} else {
  cat("Loading existing low wind events observed data\n")
  lwe_obs_pred <- readRDS(lwe_obs_pred_fname)
  lwe_obs_newloc <- readRDS(lwe_obs_newloc_fname)
}

## 3. Low wind events in model predictions ####
cat("--------------------------------------------------------------------\n")
cat("Low wind events in model mean predictions\n")
cat("--------------------------------------------------------------------\n")

low_events_model <- model_df0 %>%
  filter(model %in% est_cols) %>%
  st_drop_geometry() %>%
  # filter(model == "Spatio-temporal model") %>%
  filter(coord_id %in% coord_list$coord_id[coord_list$sampled]) %>%
  arrange(model, coord_id, time) %>%
  group_by(model, coord_id) %>%
  mutate(
    below = fit < pow_threshold,
    run_id = cumsum(below != lag(below, default = first(below)))
  ) %>%
  filter(below) %>%
  group_by(model, coord_id, run_id) %>%
  summarise(
    start_time = first(time),
    end_time = last(time),
    duration_h = as.numeric(difftime(end_time, start_time, units = "hours")) +
      1,
    .groups = "drop"
  ) %>%
  dplyr::select(model, coord_id, start_time, duration_h) %>%
  bind_rows(lwe_obs_pred) %>%
  mutate(
    model = factor(
      model,
      levels = c(est_cols, "observed"),
      labels = c(mod_labels, "observed" = "Observed")
    )
  )

# low_events_model$model %>% unique() %>% sort() %>% print()
low_events_model %>%
  filter(!model %in% excluded_models0) %>%
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
low_events_model %>%
  filter(!model %in% mod_labels[excluded_models0]) %>%
  filter(duration_h < max_h_duration) %>%
  ggplot(aes(x = duration_h)) +
  geom_density(
    # data = ~ dplyr::filter(.x, model != "observed"),
    aes(colour = model),
    alpha = 0.5,
    key_glyph = "path",
    lwd = 1
  ) +
  # geom_density(
  #   data = ~ dplyr::filter(.x, model == "observed"),
  #   aes(colour = model),
  #   # colour = "darkred",
  #   # linewidth = 1.2,
  #   fill = NA
  # ) +
  labs(
    title = "Distribution of Low Wind Events Duration",
    x = "Duration (hours)",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "inside",
    legend.position.inside = c(0.7, 0.7)
  ) +
  guides(
    colour = guide_legend(ncol = 2)
  ) +
  scale_color_manual(values = cols)

ggsave(
  sprintf(
    "fig/%s/low_wind_duration_dens_t%s.pdf",
    batch_name,
    pow_threshold_label
  ),
  width = 6,
  height = 4
)

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

probs <- seq(0, 1, length.out = 51)

qq_df <- low_events_model %>%
  filter(model != "Observed") %>%
  filter(!model %in% mod_labels[c("agg_lm", "qm", "lm_bru")]) %>%
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
