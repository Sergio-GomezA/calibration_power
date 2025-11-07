# REMIT data compilation ####

require(arrow)
require(dplyr)
require(rnaturalearth)
require(sf)
require(ggplot2)
require(ggsci)
require(ggthemes)
require(FNN)
require(data.table)
require(geosphere)

source("aux_funct.R")

## read data ####
data_path <- "~/Documents/elexon/data_by_year"

source("aux_funct.R")

# test <- get_remit(year = 2024, path = data_path, end_time = "2024-01-31")

test <- get_remit(year = 2025, path = data_path, end_time = "2025-09-30")

years <- 2019:2024

remit_full <- lapply(
  years,
  \(y) get_remit(year = y, path = data_path)
)
