# master code
gbl_start_time <- Sys.time()

# 0. Setup ####
# Get task ID and others from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
myecho <- FALSE
source('8_LT_models/0_setup.R')

# data preparation ####
source('8_LT_models/1_dataprep.R', echo = myecho)

# model fitting ####
source('8_LT_models/2_fitting.R', echo = myecho)

## predictions ####
source('8_LT_models/3_prediction.R', echo = myecho)

gbl_end_time <- Sys.time()
cat(
  "Total time taken: ",
  round(difftime(gbl_end_time, gbl_start_time, units = "mins"), 2),
  " minutes\n"
)
