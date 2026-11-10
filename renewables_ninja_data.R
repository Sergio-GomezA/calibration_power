# renewables ninja data

rn_fnames <- c(
  "~/Downloads/ninja-wind-country-GB-current_offshore-merra2.csv",
  "~/Downloads/ninja-wind-country-GB-current_onshore-merra2.csv"
)

rn_on <- data.table::fread(
  rn_fnames[1],
  skip = 3
)
