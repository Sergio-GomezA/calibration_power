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

cat("Preparing data for model fitting\n")

sampled_days <- c("2020-08-14", "2024-04-17", "2024-04-12")

d0 <- sampled_days[day_id] %>% as.Date()
d0_tag <- base::format(d0, "%y%m%d")
extension <- ifelse(local_run, "gpkg", "geojson")

model_df0 <- st_read(sprintf(
  "data/calibration_df_%s_%s.%s",
  mesh_label,
  d0_tag,
  extension
))


threshold <- 0.1
