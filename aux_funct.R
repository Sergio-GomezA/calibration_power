## functions

## get GWA data

get_gwa_data <- function(
  country,
  variable,
  height,
  url0 = "https://globalwindatlas.info/api/gis/country/GBR/wind-speed/200",
  filename = NULL,
  path = NULL,
  overwrite = FALSE
) {
  if (is.null(filename)) {
    filename <- sprintf("%s_%s_%sm.tif", country, variable, height)
  }
  if (is.null(path)) {
    dest <- filename
  } else {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
    }
    dest <- file.path(path, filename)
  }

  # construct URL
  url <- sprintf("%s/%s/%s/%s", url0, country, variable, height)

  if (file.exists(dest) && !overwrite) {
    r <- terra::rast(dest)
    return(r)
  }

  resp <- httr::GET(url, httr::write_disk(dest, overwrite = TRUE))
  httr::stop_for_status(resp)
  ct <- httr::headers(resp)[["content-type"]]

  if (is.null(ct)) {
    ct <- ""
  }
  if (!grepl("tiff|geotiff|image/tiff", ct, ignore.case = TRUE)) {
    unlink(dest)
    stop("Download did not return a GeoTIFF (content-type: ", ct, ").")
  }
  r <- terra::rast(dest)
  r
}
