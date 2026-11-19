# master code
gbl_start_time <- Sys.time()

# 0. Setup ####
# Get task ID and others from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
myecho <- FALSE
source('8_LT_models/0_setup.R')

# data preparation ####
report_memory("Before training data prep")
source('8_LT_models/1_dataprep.R', echo = myecho)

# model fitting ####
report_memory("Before fitting")
source('8_LT_models/2_fitting.R', echo = myecho)
report_memory("After fitting")

## predictions ####
source('8_LT_models/3_prediction.R', echo = myecho)
report_memory("After prediction")
gbl_end_time <- Sys.time()
cat(
  "Total time taken: ",
  round(difftime(gbl_end_time, gbl_start_time, units = "mins"), 2),
  " minutes\n"
)
