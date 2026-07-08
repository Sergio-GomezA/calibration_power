# 0. Setup ####

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

# 0.1 global parameter #####
day_id <- 1
# mesh_edge_par <- 20 # km, target edge length for the spatial mesh. 10 is fine, 20 is coarse but faster
override_objects <- FALSE
rerun_samples <- FALSE
prec_init <- log(200)


if (local_run) {
  cat("Running in local mode\n")
} else {
  cat("Running in cluster mode\n")
}

# Get task ID and others from command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Override defaults only if arguments are provided
if (length(args) > 0) {
  day_id <- as.numeric(args[1])
}
if (length(args) > 1) {
  override_objects <- as.logical(args[2])
}
if (length(args) > 2) {
  rerun_samples <- as.logical(args[3])
}
# 0.2 libraries and paths ####
require(parallel)

if (local_run) {
  data_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
  n_samp <- 100
  pixel_dims <- c(150, 150)
} else {
  data_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  gen_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  model_path <- "/exports/eddie/scratch/s2441782/calibration/model_objects"
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

sampled_days <- c("2020-08-14", "2024-04-17", "2024-04-12")
d0 <- sampled_days[day_id] %>% as.Date()
d0_tag <- base::format(d0, "%y%m%d")

n.days <- 1
n.hours <- 12
t1 <- d0 + n.days
th <- t1 + hours(n.hours)


# Read models ####
mod_labels <- c(
  # "Generic PC",
  "Linear model",
  "GB LM",
  "QM",
  "Spatio-temporal fine",
  "Spatio-temporal coarse",
  "1D SPDE model",
  "AR1 model",
  "AR2 model"
)
est_cols <- c(
  # "norm_power_est0",
  "lm",
  "agg_lm",
  "qm",
  "st0_m1",
  "st0_m2",
  "spde1d",
  "ar1",
  "ar2"
)
n_models <- length(est_cols)
names(mod_labels) <- est_cols

mod_vec <- list.files(model_path, pattern = d0_tag, full.names = TRUE)
mod_vec <- mod_vec[!grepl("spatial", mod_vec)] %>% sort() # exclude meshes from st model

model_df <- tibble(
  label = mod_labels,
  code = est_cols,
  fname = mod_vec
) %>%
  mutate(
    type = case_when(
      grepl("lm", code) ~ "lm",
      grepl("qm", code) ~ "qm",
      TRUE ~ "bru"
    )
  )

# Predictions for next hours ####

## prediction df ####
gb_day_df_fname <- sprintf("data/GB_daily_summary_%s.parquet", d0_tag)

if (!file.exists(gb_day_df_fname) || override_objects) {
  if (!file.exists(gb_day_df_fname)) {
    cat("GB daily summary file not found, creating new summary\n")
  } else {
    cat(
      "GB daily summary file found, but override_objects is TRUE. Recreating summary\n"
    )
  }
  GB_df <- read_parquet(file.path(gen_path, "GB_aggr.parquet")) %>%
    rename(time = halfHourEndTime) %>%
    mutate(
      err = norm_power_est0 - norm_potential,
      error0 = norm_potential - norm_power_est0,
      date = as.Date(time)
    )

  gb_day_df <- GB_df %>%
    group_by(date, tech_typ) %>%
    summarise(
      across(
        c(norm_power_est0, norm_potential),
        ~ sum(. * capacity) / sum(capacity)
      ),
      across(c(ws_h_wmean), ~ sum(. * capacity) / sum(capacity)),
      across(c(capacity), mean)
    ) %>%
    summarise(
      across(
        c(norm_power_est0, norm_potential),
        ~ sum(. * capacity) / sum(capacity)
      ),
      across(c(ws_h_wmean), ~ sum(. * capacity) / sum(capacity)),
      across(c(capacity), sum),
      .groups = "drop"
    )

  cutprobs3 <- c(0.25, 0.75)
  p_quant3 <- quantile(gb_day_df$norm_potential, probs = cutprobs3)
  cutprobs7 <- c(0.1, 0.2, 0.25, 0.75, 0.8, 0.9)
  p_quant7 <- quantile(gb_day_df$norm_potential, probs = cutprobs7)

  gb_day_df <- gb_day_df %>%
    mutate(
      p_group3 = cut(
        norm_potential,
        breaks = c(-Inf, p_quant3, Inf),
        labels = c("low", "mid", "high")
      ),
      p_group7 = cut(
        norm_potential,
        breaks = c(-Inf, p_quant7, Inf)
      )
    )

  write_parquet(gb_day_df, gb_day_df_fname)
} else {
  cat("Loading existing GB daily summary\n")
  gb_day_df <- read_parquet(gb_day_df_fname)
}

# extension <- ifelse(local_run, "gpkg", "geojson")
extension <- "rds"
df_pattern <- sprintf("^calibration_preddf_.*_%s\\.%s$", d0_tag, extension)
files_found <- list.files("data", pattern = df_pattern, full.names = TRUE)

coord_list_fname <- "data/coord_list.csv"

cat("Loading existing coordinate list\n")
coord_list <- read.csv(coord_list_fname)


if (!override_objects && length(files_found) > 0) {
  cat(
    "Calibration data file already exists for this day. Loading existing data.\n"
  )
  wf_df_pred <- readRDS(files_found[1])
} else {
  if (override_objects) {
    cat(
      "Override_objects is TRUE. Preparing new calibration data for this day.\n"
    )
  } else {
    cat(
      "No existing calibration data file found for this day. Preparing new data.\n"
    )
  }

  pwr_curv_df <- read_parquet(file.path(
    gen_path,
    "power_curve_all_enriched.parquet"
  ))

  wf_df_pred <- pwr_curv_df %>%
    rename(time = halfHourEndTime) %>%
    mutate(
      date = as.Date(time),
      elevation = pmax(0, elevation),
      site_name = site_name %>%
        gsub("\\b(wind\\s*farm|wf)\\b", "", ., ignore.case = TRUE) %>%
        trimws()
    ) %>%
    # filter(date %in% sampled_days) %>%
    # filter(date >= d0, date <= d0 + n.days - 1) %>%
    filter(time >= d0, time < th) %>%
    filter(coord_id %in% coord_list$coord_id[coord_list$sampled]) %>%
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
      date = as.Date(time)
    ) %>%
    left_join(
      gb_day_df %>% dplyr::select(date, p_group3),
      by = c("date" = "date")
    )

  x <- wf_df_pred$pow_group %>% unique() %>% sort()
  min_jump <- min(diff(sort(x))) / diff(range(x))
  if (min_jump <= 1e-4) {
    wf_df_pred <- wf_df_pred %>%
      mutate(pow_group = inla.group(norm_power_est0, n = 20, method = "cut"))
  }

  cat("Converting coordinates to km\n")
  wf_df_pred <- wf_df_pred %>%
    st_geometry() %>%
    (\(g) g / 1000)() %>%
    st_set_geometry(wf_df_pred, .)

  wf_df_fname <- sprintf(
    "data/calibration_preddf_%s_%s.%s",
    "base",
    d0_tag,
    extension
  )
  saveRDS(
    wf_df_pred,
    wf_df_fname
  )
}
n_loc <- nrow(wf_df_pred %>% distinct(x, y))
cat("Number of unique locations:", n_loc, "\n")
n <- nrow(wf_df_pred)

# read summary tables of predictions ######

gb_fig_df <- lapply(
  1:3,
  function(i) {
    d0 <- sampled_days[i] %>% as.Date()
    # print(i)
    d0_tag <- base::format(d0, "%y%m%d")
    readRDS(sprintf("summaries/GB_fig_band_summary_%s.rds", d0_tag)) %>%
      mutate(
        pgroup3 = case_when(
          i == 1 ~ "low",
          i == 2 ~ "mid",
          i == 3 ~ "high"
        ) %>%
          factor(levels = c("low", "mid", "high"))
      )
  }
) %>%
  bind_rows()


wf_fig_df <- lapply(
  1:3,
  function(i) {
    d0 <- sampled_days[i] %>% as.Date()
    d0_tag <- base::format(d0, "%y%m%d")
    readRDS(sprintf("summaries/WF_fig_band_summary_%s.rds", d0_tag))
  }
) %>%
  bind_rows()

# calibration fit scatter####

gb_fig_df %>%
  filter(!oos) %>%
  ggplot(aes(x = norm_potential, y = mean, col = pgroup3)) +
  geom_point() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkred"
  ) +
  facet_wrap(~model, labeller = as_labeller(mod_labels)) +
  theme_bw() +
  # scale_color_lancet() +
  scale_color_manual(values = regime_palette) +
  labs(
    x = "Normalised potential power",
    y = "Normalised predicted power",
    col = "regime"
  )

ggsave("fig/GB_fit_scatter.pdf", width = 10, height = 6, dpi = 300)

wf_fig_df %>%
  filter(!oos) %>%
  ggplot(aes(x = norm_potential, y = fit)) +
  geom_hex() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkred"
  ) +
  facet_wrap(~model, labeller = as_labeller(mod_labels)) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  theme_bw() +
  labs(
    x = "Normalised potential power",
    y = "Normalised predicted power",
    fill = "Count"
  )

ggsave("fig/WF_fit_scatter.pdf", width = 10, height = 6, dpi = 300)

## Coverage bands #####

### wf level
cov_bands_wf <- wf_fig_df %>%
  filter(oos) %>%
  group_by(model, coord_id) %>%
  summarise(
    coverage = mean(norm_potential >= lwr & norm_potential <= upr),
    .groups = "drop"
  ) %>%
  group_by(model) %>%
  summarise(
    mean_coverage = mean(coverage),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_coverage)) %>%
  mutate(
    model = factor(model, levels = model)
  )
cov_bands_wf %>%
  ggplot(aes(x = model, y = mean_coverage)) +
  geom_col(fill = blues9[7]) +
  geom_text(aes(label = round(mean_coverage, 3)), vjust = -0.5) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Mean coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf("fig/WF_pred_band_coverage.pdf"),
  width = 10,
  height = 6,
  # dpi = 300
)
### aggregated #####
cov_bands <- gb_fig_df %>%
  filter(oos) %>%
  group_by(model) %>%
  summarise(
    coverage = mean(norm_potential >= lwr & norm_potential <= upr),
    .groups = "drop"
  ) %>%
  arrange(desc(coverage)) %>%
  mutate(
    model = factor(model, levels = model)
  )

cov_bands %>%
  ggplot(aes(x = model, y = coverage)) +
  geom_col(fill = blues9[7]) +
  geom_text(aes(label = round(coverage, 3)), vjust = -0.5) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf("fig/GB_pred_band_coverage.pdf"),
  width = 10,
  height = 6,
  # dpi = 300
)

# error metrics ####
mod_labels <- c(
  # "Generic PC",
  "Linear model",
  "GB LM",
  "QM",
  "Spatio-temporal fine",
  "Spatio-temporal coarse",
  "1D SPDE model",
  "AR1 model",
  "AR2 model"
)
est_cols <- c(
  # "norm_power_est0",
  "lm",
  "agg_lm",
  "qm",
  "st0_m1",
  "st0_m2",
  "spde1d",
  "ar1",
  "ar2"
)
n_models <- length(est_cols)
names(mod_labels) <- est_cols
metrics_table <- wf_fig_df %>%
  group_by(oos, model) %>%
  summarise(
    RMSE = ModelMetrics::rmse(actual = norm_potential, predicted = fit),
    MAE = ModelMetrics::mae(actual = norm_potential, predicted = fit),
    MDAPE = mdape(
      actual = norm_potential,
      predicted = fit,
      pos_only = TRUE
    ),
    Bias = mean(fit - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    model = mod_labels[model],
    oos = ifelse(oos, "OOS", "IS")
  ) %>%
  pivot_wider(
    names_from = oos,
    values_from = c(RMSE, MAE, MDAPE, Bias)
  )
metrics_table
metrics_table %>%
  mutate(
    across(
      c(RMSE_IS, RMSE_OOS, MAE_IS, MAE_OOS, Bias_IS, Bias_OOS),
      ~ round(., 3)
    ),
    across(c(MDAPE_IS, MDAPE_OOS), ~ round(., 1))
  ) %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    align = "lcccccccc",
    col.names = c(
      "Model",
      "IS",
      "OOS",
      "IS",
      "OOS",
      "IS",
      "OOS",
      "IS",
      "OOS"
    ),
    caption = "Performance metrics for in-sample (IS) and out-of-sample (OOS) predictions."
  ) %>%
  add_header_above(c(
    " " = 1,
    "RMSE" = 2,
    "MAE" = 2,
    "MDAPE (%)" = 2,
    "Bias" = 2
  )) %>%
  kable_styling(latex_options = "hold_position")


# update currently used figures in overleaf ####

path1 <- "~/ownCloud-s2441782@datasync.ed.ac.uk/projects/calibration/calibration_power_main_doc/spfig"
path2 <- "~/ownCloud-s2441782@datasync.ed.ac.uk/projects/calibration/calibration_power/fig"

# Files to update in overleaf
files1 <- list.files(path1)

# Full paths to matching files in path2
source_files <- file.path(path2, files1)

# Keep only files that actually exist in path2
source_files <- source_files[file.exists(source_files)]

# Destination paths in path1
dest_files <- file.path(path1, basename(source_files))

# Copy and overwrite
file.copy(source_files, dest_files, overwrite = TRUE)
