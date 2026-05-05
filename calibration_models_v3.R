#  Calibration models revised

# 1. Load libraries and data ------------------------------------------------

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE


require(parallel)

if (local_run) {
  data_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
} else {
  data_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  gen_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  model_path <- "/exports/eddie/scratch/s2441782/calibration/model_objects"
  temp_lib <- "/exports/eddie3_homes_local/s2441782/lib"
  .libPaths(temp_lib)
}


# install.packages(
#   c("qmap"),
#   temp_lib,
#   dependencies = TRUE
# )

require(arrow)
require(dplyr)
require(tidyr)
# require(rnaturalearth)
# require(rnaturalearthdata)
require(sf)
require(ggplot2)
require(ggthemes)
require(ggsci)
# require(FNN)
require(data.table)
# require(parallel)
library(purrr)
# require(brms)
require(INLA)
require(inlabru)
require(qmap)
require(ggridges)

source("aux_funct.R")


# 1.1 Load data ----------------------------------------------------------------
source("read_data.R")

# 1.2 EDA ------------------------------------------------------------------------

source("eda_figures.R")

# 2. Aggregated model ------------------------------------------------------------------
# Updates:
# - Adding entire dataset
# - Adding more covariates NAO, AO, EA, SCAN

## 2.1 LM step AIC selection --------------------------------------------------------------

base_model <- lm(
  norm_potential ~ norm_power_est0,
  data = GB_df
)

full_model <- lm(
  norm_potential ~ tech_typ *
    norm_power_est0 +
    norm_power_est0 * month +
    hour +
    poly(ws_h_wmean, 3) +
    nao * tech_typ +
    ao * tech_typ +
    ea * tech_typ +
    scan * tech_typ,
  data = GB_df
)

full_model0 <- lm(
  norm_potential ~ tech_typ *
    norm_power_est0 +
    norm_power_est0 * month +
    hour +
    poly(ws_h_wmean, 3),
  data = GB_df
)
summary(full_model)


model_AIC <- step(
  base_model,
  scope = list(lower = base_model, upper = full_model),
  # steps = 5,
  k = 2
)
model_AIC0 <- step(
  base_model,
  scope = list(lower = base_model, upper = full_model0),
  # steps = 5,
  k = 2
)
model_AIC %>%
  saveRDS(
    file.path(
      model_path,
      "lm0_cir.rds"
    )
  )

ModelMetrics::rmse(GB_df$norm_potential, model_AIC$fitted.values)
ModelMetrics::rmse(GB_df$norm_potential, model_AIC0$fitted.values)

# 3. Spatially disaggregated model -----------------------------------------------------
# Updates:
# -- Add more locations (previously only 25)
# -- Add anomaly detection to handle outliers separately
# -- Add more covariates terrain, NAO, AO, EA, SCAN
# -- Add non-stationary covariance

## 3.1 Single run attempt --------------------------------------------------------------
## 3.2 Partitioned run attempt ---------------------------------------------------------

# 4. Spatially disaggregated model with non-stationary parameters ----------------------
