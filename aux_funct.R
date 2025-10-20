require(sf)
require(dplyr)
require(rnaturalearth)
require(ggplot2)
require(terra)
require(ggthemes)
require(fst)

## functions

## get GWA data

get_gwa_data <- function(
  country,
  variable,
  height,
  url0 = "https://globalwindatlas.info/api/gis/country",
  filename = NULL,
  path = NULL,
  overwrite = FALSE
) {
  # browser()
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


process_gwa_tif <- function(
  data,
  agg_fact = 5L,
  ...
) {
  if (agg_fact > 1) {
    df_agg <- data |>
      terra::aggregate(fact = agg_fact, fun = mean, na.rm = TRUE) |>
      terra::as.data.frame(., xy = TRUE, na.rm = TRUE)
  } else {
    df_agg <- terra::as.data.frame(data, xy = TRUE, na.rm = TRUE)
  }

  var_name <- names(df_agg)[3]
  meta_info <- strsplit(var_name, "_|m")[[1]]
  df_agg <- df_agg |>
    mutate(
      country = meta_info[1],
      variable = meta_info[2],
      height = meta_info[3]
    ) |>
    rename(value = !!var_name)
  df_agg
}


combined_gwa_data <- function(
  countries = c("GBR", "IRL", "IMN"),
  variable = "wind-speed",
  height = 100,
  agg_fact = 20L,
  path = "~/Documents/GWA/",
  overwrite = TRUE,
  show_fig = TRUE,
  save_df = FALSE,
  mcores = 8,
  ...
) {
  d_countries <- lapply(
    c("GBR", "IRL", "IMN"),
    \(x) {
      get_gwa_data(
        country = x,
        var = variable,
        height = height,
        path = "~/Documents/GWA/"
      )
    }
  )

  d_combined <- Reduce(merge, d_countries)

  out_path <- file.path(path, sprintf("merged_%s_%sm.tif", variable, height))

  writeRaster(d_combined, out_path, overwrite = TRUE)

  df_combined <-
    lapply(
      d_countries,
      function(x) {
        process_gwa_tif(
          data = resample(x, d_combined, threads = mcores),
          agg_fact = agg_fact
        )
      }
    ) |>
    bind_rows()

  if (save_df) {
    write_fst(
      df_combined,
      file.path(path, sprintf("df_%dm_agg%d.fst", height, agg_fact)),
      compress = 100
    )
  }

  plot <- df_combined |>
    ggplot(aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(
      option = "mako",
      na.value = "transparent",
      name = "m/s"
    ) +
    coord_quickmap() +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = "lon", y = "lat", title = sprintf("Wind speed at %dm", height))

  if (show_fig) {
    print(plot)
  }
  fig_name <- sprintf("fig/ws_%dm_merged_agg%d.png", height, agg_fact)
  ggsave(plot = plot, filename = fig_name, width = 5, height = 7)

  invisible(list(data = df_combined, plot = plot))
}
