# global_vars

# initialising workflow

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

## 0.1 global parameter #####
day_id <- 102
mesh_edge_par <- 50 # km, target edge length for the spatial mesh. 10 is fine, 20 is coarse but faster
override_objects <- FALSE
save_daily_files <- FALSE
perform_mesh_assess <- FALSE
generate_old_figs <- FALSE
re_run_st <- FALSE
prec_init <- log(200) # for u
prec_init_gau <- log(30) # for gaussian family 1DSPDE
fixed_ucomp <- FALSE
fixed_gaus_1DSPE <- FALSE
n.days.before <- 3
n.days.before.heavy <- 3
batch_name <- "batchY25d150"

task_prefix <- "spaceoos"
save_models <- FALSE
save_samples <- FALSE

rerun_samples <- TRUE

# time validation
n.days.time <- 0
n.days.before.time <- 1
n.hours.time <- 72

# space validation
n.days.space <- 0
n.days.before.space <- 10
n.hours.space <- 72

cluster_ext <- "rds" # "geojson" previously

pow_threshold <- 0.05

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
if (length(args) > 3) {
  re_run_st <- as.logical(args[4])
}
if (length(args) > 4) {
  n.days.before <- as.numeric(args[5])
}
if (length(args) > 5) {
  batch_name <- as.character(args[6])
}
if (length(args) > 6) {
  save_models <- as.logical(args[7])
}
if (length(args) > 7) {
  pow_threshold <- as.numeric(args[8])
}

## 0.2 libraries and paths ####
require(parallel)

if (local_run) {
  data_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
  extra_path <- "~/Documents/elexon/extra"
  output_path <- "~/Documents/elexon/caloutput"
  sample_path <- "~/Documents/elexon/samples"
  pixel_dims <- c(150, 150)
  n_samp <- 10
  local_ext <- "rds" # previously "gpkg"
  driver <- "GPKG"
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
  driver <- "GeoJSON"
  .libPaths(temp_lib)
}

if (!dir.exists(file.path(output_path, batch_name, "fig"))) {
  dir.create(file.path(output_path, batch_name, "fig/fit"), recursive = TRUE)
  dir.create(file.path(output_path, batch_name, "fig/oos"), recursive = TRUE)
  dir.create(file.path(output_path, batch_name, "fig/checks"), recursive = TRUE)
}
if (!dir.exists(file.path(output_path, batch_name, "summaries"))) {
  dir.create(
    file.path(output_path, batch_name, "summaries/fit"),
    recursive = TRUE
  )
  dir.create(
    file.path(output_path, batch_name, "summaries/oos"),
    recursive = TRUE
  )
  dir.create(
    file.path(output_path, batch_name, "summaries/checks"),
    recursive = TRUE
  )
}
if (!dir.exists(file.path(output_path, batch_name, "data"))) {
  dir.create(file.path(output_path, batch_name, "data/fit"), recursive = TRUE)
  dir.create(file.path(output_path, batch_name, "data/oos"), recursive = TRUE)
}
if (!dir.exists(file.path(extra_path, batch_name, "meshes"))) {
  dir.create(file.path(extra_path, batch_name, "meshes"), recursive = TRUE)
  dir.create(file.path(extra_path, batch_name, "mesh_figs"), recursive = TRUE)
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
mc <- available_cores() - ifelse(local_run, 2, 0)
inla_core_option <- "%d:1"
cat("Setting INLA to use", mc, "cores\n")
inla.setOption(num.threads = sprintf(inla_core_option, mc))
