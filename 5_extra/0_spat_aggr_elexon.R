# based on new tutorial for inlabru ####

gbl_start_time <- Sys.time()

# 0. Setup ####

source('5_extra/0_setup.R')
require(dplyr)
library(INLA)
library(inlabru)
library(fmesher)
library(ggplot2)

require(sf)
require(patchwork)

source('5_extra/1_fit_forAggTest.R')

source('5_extra/2_check_field_sim.R')

source('5_extra/3_check_full_sim.R')


gbl_end_time <- Sys.time()
cat(
  "Total time taken: ",
  round(difftime(gbl_end_time, gbl_start_time, units = "mins"), 2),
  " minutes\n"
)
