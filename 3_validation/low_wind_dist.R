# low wind events distribution

# 0. Setup ####

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

# 0.1 global parameter #####
day_id <- 2
mesh_edge_par <- 20 # km, target edge length for the spatial mesh. 10 is fine, 20 is coarse but faster
override_objects <- TRUE
# prec_init <- log(200)

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
  mesh_edge_par <- as.numeric(args[2])
}
if (length(args) > 2) {
  override_objects <- as.logical(args[3])
}

# 0.2 libraries and paths ####
require(parallel)

if (local_run) {
  data_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
  pixel_dims <- c(150, 150)
} else {
  data_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  gen_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  model_path <- "/exports/eddie/scratch/s2441782/calibration/model_objects"
  temp_lib <- "/exports/eddie3_homes_local/s2441782/lib"
  pixel_dims <- c(300, 300)
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
# require(ggspatial)

source("aux_funct.R")


# 1. data preparation ####

cat("Preparing data for low wind events analysis\n")

sampled_days <- c("2020-08-14", "2024-04-17", "2024-04-12")

d0 <- sampled_days[day_id] %>% as.Date()
d0_tag <- base::format(d0, "%y%m%d")
extension <- ifelse(local_run, "gpkg", "geojson")


model_df0 <- lapply(
  1:3,
  function(i) {
    d0 <- sampled_days[i] %>% as.Date()
    print(d0)
    d0_tag <- base::format(d0, "%y%m%d")

    st_read(sprintf(
      "data/calibration_df_%s_%s.%s",
      # mesh_label,
      "coarse",
      d0_tag,
      extension
    ))
  }
) %>%
  bind_rows()

pos_breaks <- with(
  model_df0,
  quantile(elevation[elevation > 0], probs = seq(0, 1, 1 / 3))
)
pos_levels <- levels(cut(
  model_df0$elevation[model_df0$elevation > 0],
  breaks = pos_breaks,
  include.lowest = TRUE
))
mod_labels <- c(
  "Generic PC",
  "Linear model",
  "AR1 model",
  "AR2 model",
  # "1D SPDE model",
  "Spatio-temporal model",
  "QM",
  "GB LM"
)
est_cols <- c(
  "norm_power_est0",
  "lm",
  "ar1",
  "ar2",
  # "spde1d",
  "st",
  "qm",
  "agg_lm"
)
n_models <- length(est_cols)
# n <- nrow(wf_df_frag)
names(mod_labels) <- est_cols
df_long0 <- model_df0 %>%
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
    ),
    elevation_g4 = ifelse(
      elevation == 0,
      "0",
      as.character(
        cut(
          elevation,
          breaks = pos_breaks,
          include.lowest = TRUE,
          # labels = c("Low", "Mid", "High")
        )
      )
    ),
    elevation_g4 = factor(
      elevation_g4,
      levels = c("0", pos_levels)
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


pow_threshold <- 0.1

pwr_curv_df <- read_parquet(file.path(
  gen_path,
  "power_curve_all_enriched.parquet"
))
d0 <- as.Date("2024-01-01")
n.days <- 365
pwr_coord_df <- pwr_curv_df %>%
  rename(time = halfHourEndTime) %>%
  mutate(
    date = as.Date(time)
  ) %>%
  filter(date >= d0, date <= d0 + n.days) %>%
  # filter(date %in% sampled_days) %>%
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
  )

low_events <- pwr_coord_df %>%
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
      100,
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

low_events %>% pull(duration_h) %>% summary()


low_events %>%
  filter(duration_h < 1000) %>%
  ggplot(aes(x = duration_h)) +
  geom_density(fill = blues9[5], alpha = 0.5, bw = 0.14) +
  scale_x_log10(n.breaks = 8) +
  labs(
    title = "Distribution of Low Wind Events Duration",
    x = "Duration (hours)",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# test <- df_long0 %>%
#   st_drop_geometry() %>%
#   # filter(model == "Spatio-temporal model") %>%
#   arrange(model, coord_id, time) %>%
#   mutate(
#     below = estimate < pow_threshold,
#     run_id = cumsum(below != lag(below, default = first(below)))
#   )

low_events_model <- df_long0 %>%
  st_drop_geometry() %>%
  # filter(model == "Spatio-temporal model") %>%
  arrange(model, coord_id, time) %>%
  group_by(model, coord_id) %>%
  mutate(
    below = estimate < pow_threshold,
    run_id = cumsum(below != lag(below, default = first(below)))
  ) %>%
  group_by(model, coord_id, run_id, below) %>%
  summarise(
    start_time = first(time),
    end_time = last(time),
    duration_h = as.numeric(difftime(end_time, start_time, units = "hours")) +
      1,
    .groups = "drop"
  ) %>%
  filter(below) %>%
  dplyr::select(model, coord_id, start_time, duration_h) %>%
  bind_rows(low_events)


low_events_model %>%
  filter(duration_h < 1000) %>%
  filter(
    model %in%
      c("observed", "Linear model", "AR1 model", "QM", "Spatio-temporal model")
  ) %>%
  ggplot(aes(x = duration_h)) +
  geom_density(aes(fill = model), alpha = 0.5) +
  # geom_histogram(
  #   aes(fill = model),
  #   position = "identity",
  #   alpha = 0.5,
  #   binwidth = 1
  # ) +
  # scale_x_log10(n.breaks = 14) +
  labs(
    title = "Distribution of Low Wind Events Duration",
    x = "Duration (hours)",
    y = "Frequency"
  ) +
  facet_wrap(~model, scales = "free") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave("fig/low_wind_duration_dist.pdf", width = 10, height = 6)
