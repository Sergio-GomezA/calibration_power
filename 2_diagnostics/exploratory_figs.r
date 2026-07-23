# 0. Setup ####

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

# 0.1 global parameter #####
day_id <- 1
# mesh_edge_par <- 20 # km, target edge length for the spatial mesh. 10 is fine, 20 is coarse but faster
override_objects <- FALSE
# rerun_samples <- FALSE
# prec_init <- log(200)
batch_name <- "batch2025"

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

# sampled_days <- c("2020-08-14", "2024-04-17", "2024-04-12")

sampled_days_df <- read.csv("data/sample_days_df.csv") %>%
  mutate(date = as.Date(date))

sampled_days <- sampled_days_df %>%
  pull(date) %>%
  sort()
d0 <- sampled_days[day_id] %>% as.Date()
d0_tag <- base::format(d0, "%y%m%d")

# Read models ####
mod_labels <- c(
  # "Generic PC",
  "Linear model",
  "GB LM",
  "QM",
  "Spatio-temporal fine",
  "Spatio-temporal coarse",
  "LM+hour model",
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

# model_df <- tibble(
#   label = mod_labels,
#   code = est_cols,
#   fname = mod_vec
# ) %>%
#   mutate(
#     type = case_when(
#       grepl("lm", code) ~ "lm",
#       grepl("qm", code) ~ "qm",
#       TRUE ~ "bru"
#     )
#   )

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


# sampled locations ####
coord_list_fname <- "data/coord_list.csv"

cat("Loading existing coordinate list\n")
coord_list <- read.csv(coord_list_fname)


sample_loc_fname <- "data/coord_list_wloc.csv"

if (file.exists(sample_loc_fname)) {
  cat("Loading existing coordinate list with location names\n")
  loc_cat <- read.csv(sample_loc_fname)
} else {
  cat("Creating coordinate list with location names\n")
  pwr_curv_df <- read_parquet(file.path(
    gen_path,
    "power_curve_all_enriched.parquet"
  ))
  # names(pwr_curv_df)
  loc_cat <- pwr_curv_df %>%
    distinct(coord_id, bmUnitName, site_name, lon, lat, tech_typ) %>%
    mutate(
      bmUnitName = trimws(bmUnitName),
      short_name = gsub(
        "WF|Wind|farm|Windfarm|Offshore|BMU",
        "",
        bmUnitName
      ) %>%
        gsub("\\s*\\d+$", "", .) %>% # Remove trailing numbers (and preceding spaces)
        trimws(),
      tech_typ2 = gsub("Wind", "", tech_typ) %>% trimws() %>% tolower()
    ) %>%
    group_by(coord_id) %>%
    summarise(
      short_name = first(short_name),
      site_name = first(site_name),
      tech_typ = first(tech_typ2),
      lon = first(lon),
      lat = first(lat)
    ) %>%
    mutate(
      tooltip_text = paste0(
        "Site: ",
        site_name,
        "<br>Type: ",
        tech_typ,
        "<br>coord id: ",
        coord_id
      )
    ) %>%
    arrange(coord_id) %>%
    left_join(
      coord_list %>% dplyr::select(coord_id, sampled),
      by = "coord_id"
    ) %>%
    mutate(
      sampled_label = factor(
        sampled,
        levels = c(FALSE, TRUE),
        labels = c("validation", "training")
      ),
      tech_label = factor(
        tech_typ,
        levels = c("offshore", "onshore"),
        labels = c("Offshore", "Onshore")
      )
    )

  write.csv(loc_cat, sample_loc_fname, row.names = FALSE)
}

uk_map <- rnaturalearth::ne_countries(
  scale = "medium",
  country = c("United Kingdom", "Ireland"),
  returnclass = "sf"
)

tsplit_p <- ggplot() +
  geom_sf(data = uk_map, fill = "lightgrey", color = "black") +
  geom_point(
    data = loc_cat,
    aes(
      x = lon,
      y = lat,
      color = sampled_label,
      shape = tech_label,
      text = tooltip_text
    ),
    size = 3,
    alpha = 0.7
  ) +
  scale_color_manual(
    name = "Sampled",
    values = c(
      "training" = blues9[7],
      "validation" = "darkred"
    )
  ) +
  scale_shape_manual(
    name = "Type",
    values = c(17, 16),
    labels = c("Offshore", "Onshore")
  ) +
  coord_sf(xlim = c(-10, 3), ylim = c(49.5, 61)) +
  theme_map() +
  annotation_scale(location = "br", width_hint = 0.25) +
  theme(
    legend.background = element_rect(fill = NA, color = NA),
    legend.key = element_rect(fill = NA, color = NA),
    legend.position = "inside",
    legend.position.inside = c(0.0, 0.6)
  )
tsplit_p


ggsave(
  filename = sprintf(
    "fig/%s/trainsplit_wind_farm_locations_v0.pdf",
    batch_name
  ),
  width = 4,
  height = 6
)

tsplit_p %>%
  plotly::ggplotly(tooltip = "text")


# figures for sampled days ####

## ts #####
GB_df <- read_parquet(file.path(gen_path, "GB_aggr.parquet")) %>%
  rename(time = halfHourEndTime) %>%
  mutate(
    err = norm_power_est0 - norm_potential,
    error0 = norm_potential - norm_power_est0,
    date = as.Date(time)
  )
sample_df <- GB_df %>%
  filter(date %in% sampled_days) %>%
  group_by(time) %>%
  summarise(
    across(
      c(norm_power_est0, norm_potential),
      ~ sum(. * capacity) / sum(capacity)
    ),
    across(c(ws_h_wmean), ~ sum(. * capacity) / sum(capacity)),
    across(c(capacity), sum),
    date = first(date),
    .groups = "drop"
  ) %>%
  left_join(gb_day_df %>% dplyr::select(date, p_group3), by = "date") %>%
  mutate(
    hour = format(time, "%H:%M"),
    legend_label = factor(
      paste(p_group3, format(date, "%y-%m-%d")),
      levels = paste(
        rep(c("low", "mid", "high"), each = 5),
        format(
          sampled_days_df %>%
            pull(date),
          "%y-%m-%d"
        )
      )
    )
  )

ggplot() +
  geom_line(
    data = sample_df,
    aes(
      hour,
      norm_potential,
      col = legend_label,
      group = legend_label
    )
  ) +
  labs(
    title = "Daily generation time series for sampled days",
    x = "Hour",
    y = "Wind generation (% of capacity)",
    col = ""
  ) +
  scale_x_discrete(breaks = c("00:00", "06:00", "12:00", "18:00", "23:30")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    # axis.text.x = element_text(angle = 90, vjust = 1)
  ) +
  scale_color_manual(values = rep(regime_palette, each = 5))

#
ggsave(
  sprintf("fig/%s/daily_generation_time_series_sampled_days.pdf", batch_name),
  width = 6,
  height = 4
)


# Error metrics summary ####

sampled_days_df <- read.csv("data/sample_days_df.csv") %>%
  mutate(
    date = as.Date(date),
    day_code = format(date, "%y%m%d"),
    summ_fname = paste0("summaries/calib_metrics_coarse_", day_code, ".csv"),
    calib_fname = paste0("data/calibration_df_coarse_", day_code, ".rds"),
    calib_fname2 = paste0("data/calibration_df_very_coarse_", day_code, ".rds")
  )

calib_tbl <- lapply(
  seq_along(sampled_days_df$calib_fname),
  function(x) {
    # browser()
    cutoff <- as.POSIXct(sampled_days_df$date[x])
    df1 <- readRDS(sampled_days_df$calib_fname[x]) %>%
      st_drop_geometry() %>%
      mutate(
        oos = as.POSIXct(time) >= cutoff
      ) %>%
      rename(st0_m1 = st)
    df2 <- readRDS(sampled_days_df$calib_fname2[x]) %>%
      st_drop_geometry() %>%
      mutate(
        oos = as.POSIXct(time) >= cutoff
      ) %>%
      rename(st0_m2 = st)
    df1 %>%
      left_join(
        df2 %>% dplyr::select(coord_id, time, st0_m2),
        by = c("coord_id", "time")
      )
  }
) %>%
  bind_rows()

df_long0 <- calib_tbl %>%
  dplyr::select(
    date,
    time,
    site_name,
    coord_id,
    elevation,
    dist_coast,
    capacity,
    tech_typ,
    p_group3,
    oos,
    norm_potential,
    any_of(est_cols)
  ) %>%
  mutate(
    hour = hour(time),
    elevation = pmax(0, elevation),
    p_group3 = factor(p_group3, levels = c("low", "mid", "high")),
    dist_coast_g4 = cut(
      dist_coast,
      breaks = quantile(dist_coast, probs = seq(0, 1, 0.25)),
      include.lowest = TRUE
    )
  ) %>%
  pivot_longer(
    cols = any_of(est_cols),
    names_to = "model",
    values_to = "estimate"
  ) %>%
  mutate(
    estimate = pmin(1, pmax(0, estimate)), # clipping estimates to [0, 1]
    err = estimate - norm_potential,
    p_group3 = forcats::fct_rev(p_group3),
    model = factor(model, levels = est_cols, labels = mod_labels)
  )

metrics_table <- df_long0 %>%
  st_drop_geometry() %>%
  group_by(model) %>%
  summarise(
    RMSE = ModelMetrics::rmse(actual = norm_potential, predicted = estimate),
    MAE = ModelMetrics::mae(actual = norm_potential, predicted = estimate),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    model = mod_labels[model]
  ) %>%
  arrange(desc(RMSE))

metrics_table
write.csv(
  metrics_table,
  sprintf("summaries/calib_metrics_%s.csv", batch_name),
  row.names = FALSE
)
## regime ####
metrics_table <- df_long0 %>%
  st_drop_geometry() %>%
  group_by(p_group3, model) %>%
  summarise(
    RMSE = ModelMetrics::rmse(norm_potential, estimate),
    MAE = ModelMetrics::mae(norm_potential, estimate),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    p_group3 = factor(p_group3, levels = c("low", "mid", "high"))
  ) %>%
  pivot_longer(
    cols = c(RMSE, MAE, Bias),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = c(p_group3, Metric),
    values_from = Value,
    names_vary = "fastest"
  ) %>%
  dplyr::select(
    model,
    low_RMSE,
    low_MAE,
    low_Bias,
    mid_RMSE,
    mid_MAE,
    mid_Bias,
    high_RMSE,
    high_MAE,
    high_Bias
  ) %>%
  arrange(desc(low_RMSE))
colnames(metrics_table) <- c(
  "Model",
  "RMSE",
  "MAE",
  "Bias",
  "RMSE",
  "MAE",
  "Bias",
  "RMSE",
  "MAE",
  "Bias"
)
metrics_table
kbl(
  metrics_table,
  # format = "latex",
  digits = 3,
  booktabs = TRUE,
  align = c("l", rep("c", 9))
) %>%
  add_header_above(c(
    " " = 1,
    "Low" = 3,
    "Mid" = 3,
    "High" = 3
  ))
write.csv(
  metrics_table,
  sprintf("summaries/calib_metrics_%s_regime.csv", batch_name),
  row.names = FALSE
)
## by tech type ####
metrics_table <- df_long0 %>%
  st_drop_geometry() %>%
  mutate(tech_typ = gsub("Wind", "", tech_typ) %>% trimws()) %>%
  group_by(tech_typ, model) %>%
  summarise(
    RMSE = ModelMetrics::rmse(norm_potential, estimate),
    MAE = ModelMetrics::mae(norm_potential, estimate),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(RMSE, MAE, Bias),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = c(tech_typ, Metric),
    values_from = Value,
    names_vary = "fastest"
  ) %>%
  dplyr::select(
    model,
    Offshore_RMSE,
    Offshore_MAE,
    Offshore_Bias,
    Onshore_RMSE,
    Onshore_MAE,
    Onshore_Bias
  ) %>%
  arrange(desc(Offshore_RMSE))
colnames(metrics_table) <- c(
  "Model",
  "RMSE",
  "MAE",
  "Bias",
  "RMSE",
  "MAE",
  "Bias"
)
metrics_table
kbl(
  metrics_table,
  # format = "latex",
  digits = 3,
  booktabs = TRUE,
  align = c("l", rep("c", 9))
) %>%
  add_header_above(c(
    " " = 1,
    "Offshore" = 3,
    "Onshore" = 3
  ))
write.csv(
  metrics_table,
  sprintf("summaries/calib_metrics_%s_tech.csv", batch_name),
  row.names = FALSE
)


# read summary tables of prediction bands ######
gb_fig_df <- lapply(
  seq_along(sampled_days[-15]),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    # print(i)
    d0_tag <- base::format(d0, "%y%m%d")
    readRDS(sprintf("summaries/GB_fig_band_summary_%s.rds", d0_tag)) %>%
      mutate(
        pgroup3 = sampled_days_df$p_group3[i] %>%
          factor(levels = c("low", "mid", "high"))
      )
  }
) %>%
  bind_rows()


wf_fig_df <- lapply(
  seq_along(sampled_days[-15]),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    d0_tag <- base::format(d0, "%y%m%d")
    readRDS(sprintf("summaries/WF_fig_band_summary_%s.rds", d0_tag)) %>%
      mutate(
        pgroup3 = sampled_days_df$p_group3[i] %>%
          factor(levels = c("low", "mid", "high"))
      )
  }
) %>%
  bind_rows()

## calibration fit scatter####
# #
gb_fig_df %>%
  filter(oos) %>%
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
    x = "Observed power",
    y = "Predicted power",
    col = "regime"
  )

ggsave(
  sprintf("fig/%s/GB_fit_scatter_oos.pdf", batch_name),
  width = 10,
  height = 6,
  dpi = 300
)


wf_fig_df %>%
  mutate(hour = hour(time)) %>%
  filter(oos, hour <= 1) %>%
  # filter(oos, between(hour, 10, 14)) %>%
  # filter(oos) %>%
  ggplot(aes(x = norm_potential, y = fit)) +
  geom_hex() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkred"
  ) +
  facet_wrap(~model, labeller = as_labeller(mod_labels)) +
  # scale_fill_viridis_c(
  #   trans = "log10",
  #   name = "frequency",
  #   limits = c(1, NA)
  # ) +
  scale_fill_gradient(
    trans = "log10",
    low = "grey90",
    high = blues9[5],
    # limits = c(10, NA),
  ) +
  theme_bw() +
  labs(
    x = "Normalised potential power",
    y = "Normalised predicted power",
    fill = "Count"
  )
# scales::show_col(blues9)
ggsave(
  sprintf("fig/%s/WF_fit_scatter_oos.pdf", batch_name),
  width = 10,
  height = 6,
  dpi = 300
)

## Coverage bands #####

### wf level ####
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
  filename = sprintf("fig/%s/WF_pred_band_coverage.pdf", batch_name),
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
  filename = sprintf("fig/%s/GB_pred_band_coverage.pdf", batch_name),
  width = 10,
  height = 6,
  # dpi = 300
)

cov_bands <- wf_fig_df %>%
  left_join(loc_cat %>% dplyr::select(coord_id, tech_typ), by = "coord_id") %>%
  filter(oos) %>%
  group_by(model, tech_typ) %>%
  summarise(
    coverage = mean(norm_potential >= lwr & norm_potential <= upr),
    .groups = "drop"
  ) %>%
  arrange(desc(coverage)) %>%
  mutate(
    model = factor(model, levels = model %>% unique())
  )

cov_bands %>%
  ggplot(aes(x = model, y = coverage)) +
  geom_col(fill = blues9[7]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkred") +
  geom_text(aes(label = round(coverage, 3)), vjust = -0.5) +
  facet_wrap(~tech_typ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf("fig/%s/GB_pred_band_coverage_tech.pdf", batch_name),
  width = 10,
  height = 6,
  # dpi = 300
)

cov_bands <- wf_fig_df %>%
  # left_join(loc_cat %>% dplyr::select(coord_id, tech_typ), by = "coord_id") %>%
  filter(oos) %>%
  group_by(model, pgroup3) %>%
  summarise(
    coverage = mean(norm_potential >= lwr & norm_potential <= upr),
    .groups = "drop"
  ) %>%
  arrange(desc(coverage)) %>%
  mutate(
    model = factor(model, levels = model %>% unique())
  )

cov_bands %>%
  ggplot(aes(x = model, y = coverage)) +
  geom_col(fill = blues9[7]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkred") +
  geom_text(aes(label = round(coverage, 3)), vjust = -0.5) +
  facet_wrap(~pgroup3) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf("fig/%s/GB_pred_band_coverage_regime.pdf", batch_name),
  width = 14,
  height = 6,
  # dpi = 300
)

## error metrics ####
mod_labels <- c(
  # "Generic PC",
  "Linear model",
  "GB LM",
  "QM",
  "Spatio-temporal fine",
  "Spatio-temporal coarse",
  "LM+hour model",
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
    # MDAPE = mdape(
    #   actual = norm_potential,
    #   predicted = fit,
    #   pos_only = TRUE
    # ),
    Bias = mean(fit - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    model = mod_labels[model],
    oos = ifelse(oos, "OOS", "IS")
  ) %>%
  pivot_wider(
    names_from = oos,
    values_from = c(
      RMSE,
      MAE,
      # MDAPE,
      Bias
    )
  )
metrics_table
metrics_table %>%
  mutate(
    across(
      c(RMSE_IS, RMSE_OOS, MAE_IS, MAE_OOS, Bias_IS, Bias_OOS),
      ~ round(., 3)
    ),
    # across(c(MDAPE_IS, MDAPE_OOS), ~ round(., 1))
  ) %>%
  kbl(
    # format = "latex",
    booktabs = TRUE,
    align = "lcccccccc",
    col.names = c(
      "Model",
      "IS",
      "OOS",
      "IS",
      "OOS",
      # "IS",
      # "OOS",
      "IS",
      "OOS"
    ),
    caption = "Performance metrics for in-sample (IS) and out-of-sample (OOS) predictions."
  ) %>%
  add_header_above(c(
    " " = 1,
    "RMSE" = 2,
    "MAE" = 2,
    # "MDAPE (%)" = 2,
    "Bias" = 2
  )) %>%
  kable_styling(latex_options = "hold_position")

# read summary tables of prediction bands spaceoos ######
gb_fig_df <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    # print(i)
    d0_tag <- base::format(d0, "%y%m%d")
    readRDS(sprintf(
      "summaries/GB_fig_band_summary_spaceoos_%s.rds",
      d0_tag
    )) %>%
      # mutate(
      #   pgroup3 = sampled_days_df$p_group3[i] %>%
      #     factor(levels = c("low", "mid", "high"))
      # ) %>%
      mutate(oos = time < as.POSIXct(d0))
  }
) %>%
  bind_rows() %>%
  mutate(date = as.Date(time)) %>%
  left_join(
    gb_day_df %>% dplyr::select(date, p_group3) %>% rename(pgroup3 = p_group3),
    by = "date"
  )

wf_fig_df <- lapply(
  seq_along(sampled_days[-15]),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    d0_tag <- base::format(d0, "%y%m%d")
    readRDS(sprintf(
      "summaries/WF_fig_band_summary_spaceoos_%s.rds",
      d0_tag
    )) %>%
      # mutate(
      #   pgroup3 = sampled_days_df$p_group3[i] %>%
      #     factor(levels = c("low", "mid", "high"))
      # ) %>%
      mutate(oos = time < as.POSIXct(d0))
  }
) %>%
  bind_rows() %>%
  mutate(date = as.Date(time)) %>%
  left_join(
    gb_day_df %>% dplyr::select(date, p_group3) %>% rename(pgroup3 = p_group3),
    by = "date"
  )

## calibration fit scatter####
# #
gb_fig_df %>%
  filter(oos) %>%
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
    x = "Observed power",
    y = "Predicted power",
    col = "regime"
  )

ggsave(
  sprintf("fig/%s/GB_fit_scatter_oos_space.pdf", batch_name),
  width = 10,
  height = 6,
  dpi = 300
)


wf_fig_df %>%
  filter(oos) %>%
  ggplot(aes(x = norm_potential, y = fit)) +
  geom_hex() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkred"
  ) +
  facet_wrap(~model, labeller = as_labeller(mod_labels)) +
  # scale_fill_viridis_c(
  #   trans = "log10",
  #   name = "frequency",
  #   limits = c(1, NA)
  # ) +
  scale_fill_gradient(
    trans = "log10",
    low = "grey90",
    high = blues9[8]
  ) +
  theme_bw() +
  labs(
    x = "Normalised potential power",
    y = "Normalised predicted power",
    fill = "Count"
  )

ggsave(
  sprintf("fig/%s/WF_fit_scatter_oos_space.pdf", batch_name),
  width = 10,
  height = 6,
  dpi = 300
)

## Coverage bands #####

### wf level ####
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
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkred") +
  geom_text(aes(label = round(mean_coverage, 3)), vjust = -0.5) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Mean coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf("fig/%s/WF_pred_band_coverage_spaceoos.pdf", batch_name),
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
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkred") +
  geom_text(aes(label = round(coverage, 3)), vjust = -0.5) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf("fig/%s/GB_pred_band_coverage_spaceoos.pdf", batch_name),
  width = 10,
  height = 6,
  # dpi = 300
)

cov_bands <- wf_fig_df %>%
  left_join(loc_cat %>% dplyr::select(coord_id, tech_typ), by = "coord_id") %>%
  filter(oos) %>%
  group_by(model, tech_typ) %>%
  summarise(
    coverage = mean(norm_potential >= lwr & norm_potential <= upr),
    .groups = "drop"
  ) %>%
  arrange(desc(coverage)) %>%
  mutate(
    model = factor(model, levels = model %>% unique())
  )

cov_bands %>%
  ggplot(aes(x = model, y = coverage)) +
  geom_col(fill = blues9[7]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkred") +
  geom_text(aes(label = round(coverage, 3)), vjust = -0.5) +
  facet_wrap(~tech_typ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf(
    "fig/%s/GB_pred_band_coverage_tech_spaceoos.pdf",
    batch_name
  ),
  width = 10,
  height = 6,
  # dpi = 300
)

cov_bands <- wf_fig_df %>%
  # left_join(loc_cat %>% dplyr::select(coord_id, tech_typ), by = "coord_id") %>%
  filter(oos) %>%
  group_by(model, pgroup3) %>%
  summarise(
    coverage = mean(norm_potential >= lwr & norm_potential <= upr),
    .groups = "drop"
  ) %>%
  arrange(desc(coverage)) %>%
  mutate(
    model = factor(model, levels = model %>% unique())
  )

cov_bands %>%
  ggplot(aes(x = model, y = coverage)) +
  geom_col(fill = blues9[7]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkred") +
  geom_text(aes(label = round(coverage, 3)), vjust = -0.5) +
  facet_wrap(~pgroup3) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf(
    "fig/%s/GB_pred_band_coverage_regime_spaceoos.pdf",
    batch_name
  ),
  width = 14,
  height = 6,
  # dpi = 300
)

## reliability diagram ####

# reading cov_summaries for time
cov_gbl <- lapply(
  seq_along(sampled_days[-15]),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    d0_tag <- base::format(d0, "%y%m%d")

    cov_obj <- readRDS(sprintf(
      "summaries/pred_band_coverage_summary_%s.rds",
      d0_tag
    ))
    cov_gbl <- lapply(
      seq_along(cov_obj),
      \(x) {
        cov_obj[[x]]$cov_gbl %>%
          mutate(
            model = cov_obj %>% names() %>% .[x],
            date = d0
          )
      }
    )
  }
) %>%
  bind_rows()
list.files("summaries/pred_band_coverage_summary_*.rds") %>%
  length()
model_palette <- c(
  "Observed" = "darkred",
  "Generic PC" = "#E69F00",
  "Linear model" = "#56B4E9",
  "AR1 model" = "#009E73",
  "AR2 model" = "#F0E442",
  "LM+hour model" = "#0072B2",
  "ST model coarse" = "#D55E00",
  "ST model coarser" = "#CC79A7",
  "QM" = "#999999",
  "GB LM" = "#000000"
)
model_catalog <- read.csv("data/model_catalog.csv") %>%
  na.omit()
mod_labels <- model_catalog$mod_labels
est_cols <- model_catalog$est_cols
n_models <- length(est_cols)
names(mod_labels) <- est_cols
# diabrams for time
rel_df <- cov_gbl %>%
  pivot_longer(
    cols = matches("coverage"),
    names_to = "level",
    values_to = "empirical"
  ) %>%
  mutate(
    nominal = as.numeric(gsub("coverage_", "", level)) / 100,
    model = factor(model, levels = names(mod_labels), labels = mod_labels)
  )
rel_df$model %>% unique()
rel_df %>%
  ggplot(aes(x = nominal, y = empirical, col = model)) +
  geom_line() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkred"
  ) +
  scale_color_manual(values = model_palette)

## error metrics ####
mod_labels <- c(
  # "Generic PC",
  "Linear model",
  "GB LM",
  "QM",
  "Spatio-temporal fine",
  "Spatio-temporal coarse",
  "LM+hour model",
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
    # MDAPE = mdape(
    #   actual = norm_potential,
    #   predicted = fit,
    #   pos_only = TRUE
    # ),
    Bias = mean(fit - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    model = mod_labels[model],
    oos = ifelse(oos, "OOS", "IS")
  ) %>%
  pivot_wider(
    names_from = oos,
    values_from = c(
      RMSE,
      MAE,
      # MDAPE,
      Bias
    )
  )
metrics_table
metrics_table %>%
  mutate(
    across(
      c(RMSE_IS, RMSE_OOS, MAE_IS, MAE_OOS, Bias_IS, Bias_OOS),
      ~ round(., 3)
    ),
    # across(c(MDAPE_IS, MDAPE_OOS), ~ round(., 1))
  ) %>%
  kbl(
    # format = "latex",
    booktabs = TRUE,
    align = "lcccccccc",
    col.names = c(
      "Model",
      "IS",
      "OOS",
      "IS",
      "OOS",
      # "IS",
      # "OOS",
      "IS",
      "OOS"
    ),
    caption = "Performance metrics for in-sample (IS) and out-of-sample (OOS) predictions."
  ) %>%
  add_header_above(c(
    " " = 1,
    "RMSE" = 2,
    "MAE" = 2,
    # "MDAPE (%)" = 2,
    "Bias" = 2
  )) %>%
  kable_styling(latex_options = "hold_position")

# update currently used figures in overleaf ####

path1 <- "~/ownCloud-s2441782@datasync.ed.ac.uk/projects/calibration/calibration_power_main_doc/spfig"
path2 <- "~/ownCloud-s2441782@datasync.ed.ac.uk/projects/calibration/calibration_power/fig/batch2025"

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
