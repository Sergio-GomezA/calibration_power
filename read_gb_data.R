require(arrow)
require(dplyr)

gen_path <- "~/Documents/elexon/"


GB_df <- read_parquet("GB_aggr_wind_gen.parquet")

GB_df %>%
  head()


GB_df %>% names()

# GB_df <- read_parquet(file.path(gen_path, "GB_aggr.parquet")) %>%
#   rename(time = halfHourEndTime)

# write_parquet(GB_df, "GB_aggr_wind_gen.parquet")
