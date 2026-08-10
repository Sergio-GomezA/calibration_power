# master code
gbl_start_time <- Sys.time()

# 0. Setup ####

source('4_block_code/0_setup.R')


# model fitting ####
source('4_block_code/1_models.R')

## Lite models #####

## Heavy models #####

# Sample extraction ####

## OOS time ####
source('4_block_code/2_1_oos_time.R')

## OOS space ####
source('4_block_code/2_2_oos_space.R')

## OOS space-time ####
# source('4_block_code/2_3_oos_space_time.R')

# Diagnostics and Validation #####

gbl_end_time <- Sys.time()
cat(
  "Total time taken: ",
  round(difftime(gbl_end_time, gbl_start_time, units = "mins"), 2),
  " minutes\n"
)
