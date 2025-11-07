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
data_path <- "~/Documents/elexon/"
