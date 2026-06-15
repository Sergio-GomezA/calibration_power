# 0. Setup ####

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

# 0.1 global parameter #####
day_id <- 1
mesh_edge_par <- 20 # km, target edge length for the spatial mesh. 10 is fine, 20 is coarse but faster
override_objects <- TRUE
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

sampled_days <- c("2020-08-14", "2024-04-17", "2024-04-12")
d0 <- sampled_days[day_id] %>% as.Date()
d0_tag <- base::format(d0, "%y%m%d")

# Read models ####
mod_vec <- list.files(model_path, pattern = d0_tag, full.names = TRUE)
mod_vec <- mod_vec[!grepl("mesh", mod_vec)]

# Predictions for next hours ####

## prediction df ####

extension <- ifelse(local_run, "gpkg", "geojson")
df_pattern <- sprintf("^calibration_df_.*_%s\\.%s$", d0_tag, extension)
files_found <- list.files("data", pattern = df_pattern, full.names = TRUE)

if (!override_objects && length(files_found) > 0) {
  cat(
    "Calibration data file already exists for this day. Loading existing data.\n"
  )
  wf_df_frag <- st_read(files_found[1])
} else {
  cat(
    "No existing calibration data file found for this day. Preparing new data.\n"
  )

  pwr_curv_df <- read_parquet(file.path(
    gen_path,
    "power_curve_all_enriched.parquet"
  ))

  n.days <- 1

  wf_df_frag <- pwr_curv_df %>%
    rename(time = halfHourEndTime) %>%
    mutate(
      date = as.Date(time),
      elevation = pmax(0, elevation),
      site_name = site_name %>%
        gsub("\\b(wind\\s*farm|wf)\\b", "", ., ignore.case = TRUE) %>%
        trimws()
    ) %>%
    # filter(date %in% sampled_days) %>%
    filter(date >= d0, date <= d0 + n.days - 1) %>%
    arrange(site_name) %>%
    # mutate(
    #   site_id = as.integer(factor(site_name)),
    #   coord_id = as.integer(factor(paste(lon, lat)))
    # ) %>%
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

  x <- wf_df_frag$pow_group %>% unique() %>% sort()
  min_jump <- min(diff(sort(x))) / diff(range(x))
  if (min_jump <= 1e-4) {
    wf_df_frag <- wf_df_frag %>%
      mutate(pow_group = inla.group(norm_power_est0, n = 20, method = "cut"))
  }

  cat("Converting coordinates to km\n")
  wf_df_frag <- wf_df_frag %>%
    st_geometry() %>%
    (\(g) g / 1000)() %>%
    st_set_geometry(wf_df_frag, .)

  wf_df_fname <- sprintf(
    "data/calibration_df_%s_%s.%s",
    "base",
    d0_tag,
    extension
  )
  if (extension == "geojson" & file.exists(wf_df_fname)) {
    file.remove(wf_df_fname)
  }
  st_write(
    wf_df_frag,
    wf_df_fname,
    driver = ifelse(local_run, "GPKG", "GeoJSON"),
    append = FALSE,
    quiet = TRUE
  )
}
cat("Number of unique locations:", nrow(wf_df_frag %>% distinct(x, y)), "\n")
n <- nrow(wf_df_frag)
## linear models ####

## quantile mapping ####

## bru models ####

# Plot uncertainty bands ####
