require(data.table)
require(dplyr)
require(tidyr)
require(ggplot2)
require(ggsci)
require(terra)
require(fst)
require(stringr)
require(httr)
require(jsonlite)
require(arrow)
require(geosphere)

theme_set(theme_bw())

wind.bmus <- fread(file.path("data", "wind_bmu_2.csv.gz"))
wind.bmus.alt <- read.csv("data/wind_bmu_alt.csv")

mypalette <- pal_aaas()(3)[c(1, 3)]
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
      file.path(
        path,
        sprintf("df_%s_%dm_agg%d.fst", variable, height, agg_fact)
      ),
      compress = 100
    )
  }

  plot <- df_combined |>
    ggplot(aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(
      option = "mako",
      na.value = "transparent",
      name = case_when(
        grepl("Weibull", variable) ~ "",
        TRUE ~ "m/s"
      ),
      ...
    ) +
    coord_quickmap() +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
      x = "lon",
      y = "lat",
      title = sprintf("%s at %dm", gsub("-", " ", variable), height)
    )

  if (show_fig) {
    print(plot)
  }
  fig_name <- sprintf("fig/%s_%dm_merged_agg%d.png", variable, height, agg_fact)
  ggsave(plot = plot, filename = fig_name, width = 5, height = 7)

  invisible(list(data = df_combined, plot = plot))
}


clean_names <- function(name) {
  name %>%
    iconv(from = "", to = "UTF-8", sub = "") %>% # Fix encoding
    tolower() %>% # Lowercase
    str_replace_all("[^a-z0-9 ]", " ") %>% # Remove punctuation
    str_replace_all("\\s+", " ") %>% # Remove extra spaces
    str_trim()
}


get_agpt <- function(
  year = 2024,
  path,
  file_name = sprintf("aggregated_gen_type_%d.csv.gz", year),
  end_time = NULL
) {
  if (is.null(end_time)) {
    time_seq <- dates <- seq.Date(
      from = as.Date(sprintf("%d-01-01", year)),
      to = as.Date(sprintf("%d-12-31", year)),
      by = "4 days"
    )
    time_seq <- c(time_seq, as.Date(sprintf("%d-01-01", year + 1)))
  } else {
    time_seq <- dates <- seq.Date(
      from = as.Date(sprintf("%d-01-01", year)),
      to = as.Date(end_time),
      by = "4 days"
    )
    if (time_seq[length(time_seq)] < end_time) time_seq <- c(time_seq, end_time)
  }

  # browser()
  extraction <- lapply(
    seq_along(time_seq[-1]),

    \(i) {
      t0 <- as.character(time_seq[i])
      t1 <- as.character(time_seq[i + 1])
      dt0 <- paste(t0, "01:00")
      dt1 <- paste(t1, "01:00")
      dt0_enc <- URLencode(dt0, reserved = TRUE)
      dt1_enc <- URLencode(dt1, reserved = TRUE)

      base_url <- "https://data.elexon.co.uk/bmrs/api/v1/datasets/AGPT/stream?"

      url <- sprintf(
        "%spublishDateTimeFrom=%s&publishDateTimeTo=%s&format=json",
        base_url,
        dt0_enc,
        dt1_enc
      )

      response <- GET(url, accept("text/plain"))
      json_data <- content(response, "text", encoding = "UTF-8")
      # browser()
      agpt0 <- json_data %>%
        fromJSON(flatten = T) %>%
        mutate(across(
          matches("Time"),
          ~ as.POSIXct(., format = "%Y-%m-%dT%H:%M:%OS")
        ))
      agpt0
    }
  )

  agptyear <- extraction %>% bind_rows()

  agptyear %>%
    write.csv(
      gzfile(file.path(path, file_name))
    )

  invisible(agptyear)
}


get_dwgs <- function(
  year = 2024,
  path,
  file_name = sprintf("dayahead_generation_windsolar_%d.csv.gz", year)
) {
  time_seq <- dates <- seq.Date(
    from = as.Date(sprintf("%d-01-01", year)),
    to = as.Date(sprintf("%d-12-31", year)),
    by = "4 days"
  )

  # if(day(time_seq[length(time_seq)])!=31)

  time_seq <- c(time_seq, as.Date(sprintf("%d-01-01", year + 1)))

  # browser()
  extraction <- lapply(
    seq_along(time_seq[-1]),

    \(i) {
      t0 <- as.character(time_seq[i])
      t1 <- as.character(time_seq[i + 1])
      dt0 <- paste(t0, "01:00")
      dt1 <- paste(t1, "01:00")
      dt0_enc <- URLencode(dt0, reserved = TRUE)
      dt1_enc <- URLencode(dt1, reserved = TRUE)

      base_url <- "https://data.elexon.co.uk/bmrs/api/v1/datasets/DGWS/stream?"

      url <- sprintf(
        "%spublishDateTimeFrom=%s&publishDateTimeTo=%s&format=json",
        base_url,
        dt0_enc,
        dt1_enc
      )

      response <- GET(url, accept("text/plain"))
      json_data <- content(response, "text", encoding = "UTF-8")
      # browser()
      dwgs0 <- tryCatch(
        json_data %>%
          fromJSON(flatten = T), #%>%
        # mutate(across(matches("Time"), ~ as.POSIXct(., format = "%Y-%m-%dT%H:%M:%OS")))
        error = function(e) {
          warning(sprintf(
            "Data extraction for day %s failed. Response: %s\n",
            t0,
            response$status_code
          ))
          return(NULL)
        }
      )
      dwgs0
    }
  )

  dwgsyear <- extraction %>% bind_rows()

  dwgsyear %>%
    write.csv(
      gzfile(file.path(path, file_name))
    )

  invisible(dwgsyear)
}

get_bmugen <- function(
  year = 2024,
  path,
  file_name = sprintf("wind_gen_bmu_%d.csv.gz", year),
  end_time = NULL
) {
  if (is.null(end_time)) {
    time_seq <- dates <- seq.Date(
      from = as.Date(sprintf("%d-01-01", year)),
      to = as.Date(sprintf("%d-12-31", year)),
      by = "4 days"
    )
    time_seq <- c(time_seq, as.Date(sprintf("%d-01-01", year + 1)))
  } else {
    time_seq <- dates <- seq.Date(
      from = as.Date(sprintf("%d-01-01", year)),
      to = as.Date(end_time),
      by = "4 days"
    )
    if (time_seq[length(time_seq)] < end_time) time_seq <- c(time_seq, end_time)
  }

  # browser()
  extraction <- lapply(
    seq_along(time_seq[-1]),

    \(i) {
      t0 <- as.character(time_seq[i])
      t1 <- as.character(time_seq[i + 1])
      dt0 <- paste(t0, "01:00")
      dt1 <- paste(t1, "01:00")
      dt0_enc <- URLencode(dt0, reserved = TRUE)
      dt1_enc <- URLencode(dt1, reserved = TRUE)

      base_url <- "https://data.elexon.co.uk/bmrs/api/v1/datasets/B1610/stream?"

      url <- sprintf(
        "%sfrom=%s&to=%s&bmUnit=",
        base_url,
        dt0_enc,
        dt1_enc
      )

      fragments <- 3
      n_bmu <- nrow(wind.bmus)
      parts_vec <- seq(1, n_bmu, length.out = fragments + 1) %>% trunc()

      # list strings in fragments
      bmu_strings <- mapply(
        \(left, right) {
          bmu.string <- wind.bmus %>%
            slice(left:right) %>%
            pull(elexonBmUnit) %>%
            paste(., collapse = "&bmUnit=")

          # print(bmu.string)
        },
        left = parts_vec[1:fragments] + c(0, rep(1, fragments - 1)),
        right = parts_vec[-1]
      )

      # query data by fragment
      bmu_df <- lapply(
        bmu_strings,
        \(string) {
          url <- paste0(url, string)
          response <- GET(url, accept("text/plain"))
          json_data <- content(response, "text", encoding = "UTF-8")
          # browser() sf
          generation.bmu <- tryCatch(
            fromJSON(json_data) %>%
              mutate(across(
                matches("Time"),
                ~ as.POSIXct(., format = "%Y-%m-%dT%H:%M:%OS")
              )),
            error = function(e) {
              warning(sprintf(
                "Data extraction for day %s failed. Response: %s\n",
                t0,
                response$status_code
              ))
              return(NULL)
            }
          )
        }
      ) %>%
        bind_rows()

      bmu_df
    }
  )

  bmugenyear <- extraction %>% bind_rows()

  bmugenyear %>%
    write.csv(
      gzfile(file.path(path, file_name))
    )

  invisible(bmugenyear)
}

get_bod <- function(
  year = 2024,
  path,
  file_name = sprintf("bid_offer_data_%d.csv.gz", year),
  end_time = NULL
) {
  if (is.null(end_time)) {
    time_seq <- dates <- seq.Date(
      from = as.Date(sprintf("%d-01-01", year)),
      to = as.Date(sprintf("%d-12-31", year)),
      by = "4 days"
    )
    time_seq <- c(time_seq, as.Date(sprintf("%d-01-01", year + 1)))
  } else {
    time_seq <- dates <- seq.Date(
      from = as.Date(sprintf("%d-01-01", year)),
      to = as.Date(end_time),
      by = "4 days"
    )
    if (time_seq[length(time_seq)] < end_time) time_seq <- c(time_seq, end_time)
  }

  # browser()
  extraction <- lapply(
    seq_along(time_seq[-1]),

    \(i) {
      t0 <- as.character(time_seq[i])
      t1 <- as.character(time_seq[i + 1])
      dt0 <- paste(t0, "01:00")
      dt1 <- paste(t1, "01:00")
      dt0_enc <- URLencode(dt0, reserved = TRUE)
      dt1_enc <- URLencode(dt1, reserved = TRUE)

      base_url <- "https://data.elexon.co.uk/bmrs/api/v1/datasets/BOD/stream?"

      url <- sprintf(
        "%sfrom=%s&to=%s&bmUnit=",
        base_url,
        dt0_enc,
        dt1_enc
      )

      fragments <- 3
      n_bmu <- nrow(wind.bmus)
      parts_vec <- seq(1, n_bmu, length.out = fragments + 1) %>% trunc()

      # list strings in fragments
      bmu_strings <- mapply(
        \(left, right) {
          bmu.string <- wind.bmus %>%
            slice(left:right) %>%
            pull(elexonBmUnit) %>%
            paste(., collapse = "&bmUnit=")

          # print(bmu.string)
        },
        left = parts_vec[1:fragments] + c(0, rep(1, fragments - 1)),
        right = parts_vec[-1]
      )

      # query data by fragment
      bmu_df <- lapply(
        bmu_strings,
        \(string) {
          url <- paste0(url, string)
          response <- GET(url, accept("text/plain"))
          json_data <- content(response, "text", encoding = "UTF-8")
          generation.bmu <- tryCatch(
            fromJSON(json_data) %>%
              mutate(across(
                matches("Time"),
                ~ as.POSIXct(., format = "%Y-%m-%dT%H:%M:%OS")
              )),
            error = function(e) {
              warning(sprintf(
                "Data extraction for day %s failed. Response: %s\n",
                t0,
                response$status_code
              ))
              return(NULL)
            }
          )
        }
      ) %>%
        bind_rows()

      bmu_df
    }
  )

  bmugenyear <- extraction %>% bind_rows()

  bmugenyear %>%
    write.csv(
      gzfile(file.path(path, file_name))
    )

  invisible(bmugenyear)
}

get_boa <- function(
  year = 2024,
  path,
  file_name = sprintf("bid_offer_accept_%d.csv.gz", year),
  end_time = NULL
) {
  if (is.null(end_time)) {
    time_seq <- dates <- seq.Date(
      from = as.Date(sprintf("%d-01-01", year)),
      to = as.Date(sprintf("%d-12-31", year)),
      by = "4 days"
    )
    time_seq <- c(time_seq, as.Date(sprintf("%d-01-01", year + 1)))
  } else {
    time_seq <- dates <- seq.Date(
      from = as.Date(sprintf("%d-01-01", year)),
      to = as.Date(end_time),
      by = "4 days"
    )
    if (time_seq[length(time_seq)] < end_time) time_seq <- c(time_seq, end_time)
  }

  # browser()
  extraction <- lapply(
    seq_along(time_seq[-1]),

    \(i) {
      t0 <- as.character(time_seq[i])
      t1 <- as.character(time_seq[i + 1])
      dt0 <- paste(t0, "01:00")
      dt1 <- paste(t1, "01:00")
      dt0_enc <- URLencode(dt0, reserved = TRUE)
      dt1_enc <- URLencode(dt1, reserved = TRUE)

      base_url <- "https://data.elexon.co.uk/bmrs/api/v1/datasets/BOALF/stream?"

      url <- sprintf(
        "%sfrom=%s&to=%s&bmUnit=",
        base_url,
        dt0_enc,
        dt1_enc
      )

      fragments <- 3
      n_bmu <- nrow(wind.bmus)
      parts_vec <- seq(1, n_bmu, length.out = fragments + 1) %>% trunc()

      # list strings in fragments
      bmu_strings <- mapply(
        \(left, right) {
          bmu.string <- wind.bmus %>%
            slice(left:right) %>%
            pull(elexonBmUnit) %>%
            paste(., collapse = "&bmUnit=")

          # print(bmu.string)
        },
        left = parts_vec[1:fragments] + c(0, rep(1, fragments - 1)),
        right = parts_vec[-1]
      )

      # query data by fragment
      bmu_df <- lapply(
        bmu_strings,
        \(string) {
          url <- paste0(url, string)
          response <- GET(url, accept("text/plain"))
          json_data <- content(response, "text", encoding = "UTF-8")
          generation.bmu <- tryCatch(
            fromJSON(json_data) %>%
              mutate(across(
                matches("Time"),
                ~ as.POSIXct(., format = "%Y-%m-%dT%H:%M:%OS")
              )),
            error = function(e) {
              warning(sprintf(
                "Data extraction for day %s failed. Response: %s\n",
                t0,
                response$status_code
              ))
              return(NULL)
            }
          )
        }
      ) %>%
        bind_rows()

      bmu_df
    }
  )

  bmugenyear <- extraction %>% bind_rows()

  bmugenyear %>%
    write.csv(
      gzfile(file.path(path, file_name))
    )

  invisible(bmugenyear)
}

get_pn <- function(
  year = 2024,
  path,
  file_name = sprintf("phys_notif_%d.csv.gz", year),
  end_time = NULL
) {
  if (is.null(end_time)) {
    time_seq <- dates <- seq.Date(
      from = as.Date(sprintf("%d-01-01", year)),
      to = as.Date(sprintf("%d-12-31", year)),
      by = "4 days"
    )
    time_seq <- c(time_seq, as.Date(sprintf("%d-01-01", year + 1)))
  } else {
    time_seq <- dates <- seq.Date(
      from = as.Date(sprintf("%d-01-01", year)),
      to = as.Date(end_time),
      by = "4 days"
    )
    if (time_seq[length(time_seq)] < end_time) time_seq <- c(time_seq, end_time)
  }

  # browser()
  extraction <- lapply(
    seq_along(time_seq[-1]),

    \(i) {
      t0 <- as.character(time_seq[i])
      t1 <- as.character(time_seq[i + 1])
      dt0 <- paste(t0, "01:00")
      dt1 <- paste(t1, "01:00")
      dt0_enc <- URLencode(dt0, reserved = TRUE)
      dt1_enc <- URLencode(dt1, reserved = TRUE)

      base_url <- "https://data.elexon.co.uk/bmrs/api/v1/datasets/PN/stream?"

      url <- sprintf(
        "%sfrom=%s&to=%s&bmUnit=",
        base_url,
        dt0_enc,
        dt1_enc
      )

      fragments <- 3
      n_bmu <- nrow(wind.bmus)
      parts_vec <- seq(1, n_bmu, length.out = fragments + 1) %>% trunc()

      # list strings in fragments
      bmu_strings <- mapply(
        \(left, right) {
          bmu.string <- wind.bmus %>%
            slice(left:right) %>%
            pull(elexonBmUnit) %>%
            paste(., collapse = "&bmUnit=")

          # print(bmu.string)
        },
        left = parts_vec[1:fragments] + c(0, rep(1, fragments - 1)),
        right = parts_vec[-1]
      )

      # query data by fragment
      bmu_df <- lapply(
        bmu_strings,
        \(string) {
          url <- paste0(url, string)
          response <- GET(url, accept("text/plain"))
          json_data <- content(response, "text", encoding = "UTF-8")
          generation.bmu <- tryCatch(
            fromJSON(json_data) %>%
              mutate(across(
                matches("Time"),
                ~ as.POSIXct(., format = "%Y-%m-%dT%H:%M:%OS")
              )),
            error = function(e) {
              warning(sprintf(
                "Data extraction for day %s failed. Response: %s\n",
                t0,
                response$status_code
              ))
              return(NULL)
            }
          )
        }
      ) %>%
        bind_rows()

      bmu_df
    }
  )

  bmugenyear <- extraction %>% bind_rows()

  bmugenyear %>%
    write.csv(
      gzfile(file.path(path, file_name))
    )

  invisible(bmugenyear)
}

get_outturn <- function(
  year = 2024,
  path,
  file_name = sprintf("outturn_gen_type_%d.csv.gz", year),
  end_time = NULL
) {
  if (is.null(end_time)) {
    time_seq <- dates <- seq.Date(
      from = as.Date(sprintf("%d-01-01", year)),
      to = as.Date(sprintf("%d-12-31", year)),
      by = "4 days"
    )
    time_seq <- c(time_seq, as.Date(sprintf("%d-01-01", year + 1)))
  } else {
    time_seq <- dates <- seq.Date(
      from = as.Date(sprintf("%d-01-01", year)),
      to = as.Date(end_time),
      by = "4 days"
    )
    if (time_seq[length(time_seq)] < end_time) time_seq <- c(time_seq, end_time)
  }

  # browser()
  extraction <- lapply(
    seq_along(time_seq),

    \(i) {
      t0 <- as.character(time_seq[i])
      # t1 <- as.character(time_seq[i + 1])
      dt0 <- paste(t0, "01:00")
      # dt1 <- paste(t1, "01:00")
      dt0_enc <- URLencode(dt0, reserved = TRUE)
      # dt1_enc <- URLencode(dt1, reserved = TRUE)

      base_url <- "https://data.elexon.co.uk/bmrs/api/v1/balancing/settlement/indicative/volumes/all/bid/"

      url <- sprintf(
        "%s%s?&bmUnit=",
        base_url,
        dt0_enc
      )

      fragments <- 3
      n_bmu <- nrow(wind.bmus)
      parts_vec <- seq(1, n_bmu, length.out = fragments + 1) %>% trunc()

      # list strings in fragments
      bmu_strings <- mapply(
        \(left, right) {
          bmu.string <- wind.bmus %>%
            slice(left:right) %>%
            pull(elexonBmUnit) %>%
            paste(., collapse = "&bmUnit=")

          # print(bmu.string)
        },
        left = parts_vec[1:fragments] + c(0, rep(1, fragments - 1)),
        right = parts_vec[-1]
      )

      # query data by fragment
      bmu_df <- lapply(
        bmu_strings,
        \(string) {
          url <- paste0(url, string)
          response <- GET(url, accept("text/plain"))
          json_data <- content(response, "text", encoding = "UTF-8")
          generation.bmu <- tryCatch(
            fromJSON(json_data) %>%
              mutate(across(
                matches("Time"),
                ~ as.POSIXct(., format = "%Y-%m-%dT%H:%M:%OS")
              )),
            error = function(e) {
              warning(sprintf(
                "Data extraction for day %s failed. Response: %s\n",
                t0,
                response$status_code
              ))
              return(NULL)
            }
          )
        }
      ) %>%
        bind_rows()
    }
  )

  df <- extraction %>% bind_rows()

  df %>%
    write.csv(
      gzfile(file.path(path, file_name))
    )

  invisible(df)
}


get_curtailment <- function(
  year = 2024,
  path,
  file_name = sprintf("wind_curt_bmu_%d.csv.gz", year),
  end_time = NULL
) {
  # browser()
  if (is.null(end_time)) {
    time_seq <- dates <- seq.Date(
      from = as.Date(sprintf("%d-01-01", year)),
      to = as.Date(sprintf("%d-12-31", year)),
      by = "1 day"
    )
    time_seq <- c(time_seq, as.Date(sprintf("%d-01-01", year + 1)))
  } else {
    time_seq <- dates <- seq.Date(
      from = as.Date(sprintf("%d-01-01", year)),
      to = as.Date(end_time),
      by = "1 day"
    )
    if (time_seq[length(time_seq)] < end_time) time_seq <- c(time_seq, end_time)
  }

  extraction <- lapply(
    seq_along(time_seq),

    \(i) {
      # if (i == length(time_seq)) {
      #   browser()
      # }
      t0 <- as.character(time_seq[i])
      # t1 <- as.character(time_seq[i + 1])
      # dt0 <- paste(t0, "01:00")
      # dt1 <- paste(t1, "01:00")
      # dt0_enc <- URLencode(dt0, reserved = TRUE)
      # dt1_enc <- URLencode(dt1, reserved = TRUE)

      base_url <- "https://data.elexon.co.uk/bmrs/api/v1/balancing/settlement/indicative/volumes/all/bid/"

      url <- sprintf(
        "%s%s?&bmUnit=",
        base_url,
        t0
      )

      fragments <- 3
      n_bmu <- nrow(wind.bmus)
      parts_vec <- seq(1, n_bmu, length.out = fragments + 1) %>% trunc()

      # list strings in fragments
      bmu_strings <- mapply(
        \(left, right) {
          bmu.string <- wind.bmus %>%
            slice(left:right) %>%
            pull(elexonBmUnit) %>%
            paste(., collapse = "&bmUnit=")

          # print(bmu.string)
        },
        left = parts_vec[1:fragments] + c(0, rep(1, fragments - 1)),
        right = parts_vec[-1]
      )

      # query data by fragment
      bmu_df <- lapply(
        bmu_strings,
        \(string) {
          # browser()
          url <- paste0(url, string)
          response <- GET(url, accept("text/plain"))
          json_data <- content(response, "text", encoding = "UTF-8")
          generation.bmu <- tryCatch(
            fromJSON(json_data)$data,
            error = function(e) {
              warning(sprintf(
                "Data extraction for day %s failed. Response: %s\n",
                t0,
                response$status_code
              ))
              return(NULL)
            }
          )
        }
      ) %>%
        bind_rows() %>%
        mutate(across(
          matches("Time"),
          ~ as.POSIXct(., format = "%Y-%m-%dT%H:%M:%OS")
        ))
    }
  )

  df <- extraction %>% bind_rows()

  df %>%
    write.csv(
      gzfile(file.path(path, file_name))
    )

  invisible(df)
}

get_remit <- function(
  year = 2024,
  path,
  file_name = sprintf("remit_all_%d.parquet", year),
  end_time = NULL
) {
  if (is.null(end_time)) {
    time_seq <- dates <- seq.Date(
      from = as.Date(sprintf("%d-01-01", year)),
      to = as.Date(sprintf("%d-12-31", year)),
      by = "7 days"
    )
    time_seq <- c(time_seq, as.Date(sprintf("%d-01-01", year + 1)))
  } else {
    time_seq <- dates <- seq.Date(
      from = as.Date(sprintf("%d-01-01", year)),
      to = as.Date(end_time),
      by = "7 days"
    )
    if (time_seq[length(time_seq)] < end_time) time_seq <- c(time_seq, end_time)
  }

  # browser()
  extraction <- lapply(
    seq_along(time_seq[-1]),

    \(i) {
      t0 <- as.character(time_seq[i])
      t1 <- as.character(time_seq[i + 1])
      dt0 <- paste(t0, "01:00")
      dt1 <- paste(t1, "01:00")
      dt0_enc <- URLencode(dt0, reserved = TRUE)
      dt1_enc <- URLencode(dt1, reserved = TRUE)

      base_url <- "https://data.elexon.co.uk/bmrs/api/v1/remit/list/by-event/stream?"

      url <- sprintf(
        "%sfrom=%s&to=%s&latestRevisionOnly=true&profileOnly=false",
        base_url,
        dt0_enc,
        dt1_enc
      )
      # browser()
      response <- GET(url, accept("text/plain"))
      json_data <- content(response, "text", encoding = "UTF-8")
      msg_tbl <- fromJSON(json_data)

      msg_tbl <- tryCatch(
        msg_tbl %>%
          mutate(
            json = purrr::map(url, function(u) {
              resp <- request(u) |> req_perform()
              resp_body_json(resp, simplifyVector = TRUE)
            })
          ), #
        error = function(e) {
          warning(sprintf(
            "Data extraction for day %s failed. message: %s\n",
            t0,
            msg_tbl
          ))
          browser()
          return(NULL)
        }
      )
      purrr::map(msg_tbl$json, "data") %>% bind_rows()
    }
  )

  df <- extraction %>% bind_rows()

  write_parquet(
    df,
    file.path(path, file_name)
  )
  invisible(df)
}

interp_log_ws <- function(h, z1, z2, u1, u2) {
  u1 + (log(h / z1) / log(z2 / z1)) * (u2 - u1)
}


# pick generic pc and rescale

generic_pow_conv <- function(
  wind_speed,
  turb_class = "offshore",
  turb_capacity = 4
) {
  # browser()
  class_curve <- fread("data/generic_powerCurves.csv.gz") %>%
    filter(class == turb_class)

  # max_rated_power = max(class_curve$ratedPower)

  class_curve <- class_curve %>%
    mutate(power_scaled = power_kw * turb_capacity / ratedPower)

  power_est <- approx(
    x = class_curve$wind_speed,
    y = class_curve$power_scaled,
    xout = wind_speed,
    rule = 2 # flat extrapolation
  )$y

  power_est
}


power_curve_data1f <- function(
  bmu_code,
  t0 = "2024-01-01",
  t1 = "2024-12-31",
  col_palette = c(
    "observed" = blues9[7],
    "generic power curve" = "darkred",
    "outages" = "darkorange"
  )
) {
  # browser()
  turb_class <- ref_catalog_2025 %>%
    filter(grepl(bmu_code, bmUnit)) %>%
    pull(turb_class) %>%
    unique()

  pwr_curv_1wf <- gen_adj %>%
    filter(grepl(bmu_code, bmUnit)) %>%
    mutate(
      quantity = quantity + lag(quantity),
      curtailment = curtailment + lag(curtailment),
      potential = potential + lag(potential)
    ) %>%
    filter(minute(halfHourEndTime) == 0) %>%
    left_join(
      ref_catalog_2025 %>%
        select(
          bmUnit,
          matches("lon|lat"),
          site_name,
          tech_typ,
          turb_class,
          height_turb_imp
        ),
      by = c("bmUnit")
    ) %>%
    left_join(
      era_df %>% select(time, longitude, latitude, ws100, wd100, ws10, wd10),
      by = c(
        "halfHourEndTime" = "time",
        "era5lon" = "longitude",
        "era5lat" = "latitude"
      )
    ) %>%
    # wind speed vertical interpolation
    mutate(
      ws_log = log(height_turb_imp / 10) /
        log(100 / 10) *
        (ws100 - ws10) +
        ws10,
      ws_h = ws100 * (height_turb_imp / 100)^(1 / 7)
    )

  scaled_pc <- generic_pc %>%
    filter(class %in% turb_class) %>%
    mutate(power_scaled = power_kw * pwr_curv_1wf$capacity[1] / ratedPower)

  p_quant <- pwr_curv_1wf %>%
    filter(between(halfHourEndTime, t0, t1)) %>%
    ggplot(aes(ws_h, quantity, col = "observed")) +
    geom_point(alpha = 0.2) +
    geom_line(
      data = scaled_pc,
      aes(wind_speed, power_scaled, col = "generic power curve")
    ) +
    scale_color_manual(
      values = col_palette
    ) +
    theme(
      legend.position = "bottom"
    ) +
    labs(col = "", x = "wind speed", y = "generation (MW)")

  p_pot <- pwr_curv_1wf %>%
    filter(between(halfHourEndTime, t0, t1)) %>%
    ggplot(aes(ws_h, potential, col = "observed")) +
    geom_point(alpha = 0.2) +
    geom_line(
      data = scaled_pc,
      aes(wind_speed, power_scaled, col = "generic power curve")
    ) +
    scale_color_manual(
      values = col_palette
    ) +
    theme(
      legend.position = "bottom"
    ) +
    labs(col = "", x = "wind speed", y = "potential generation (MW)")

  p_out <- pwr_curv_1wf %>%
    filter(between(halfHourEndTime, t0, t1)) %>%
    ggplot(aes(ws_h, potential, col = "observed")) +
    geom_point(alpha = 0.2) +
    geom_point(
      aes(ws_h, potential, col = "outages"),
      data = pwr_curv_1wf %>%
        filter(between(halfHourEndTime, t0, t1), !is.na(outageCapacity)),
      alpha = 0.2
    ) +
    geom_line(
      data = scaled_pc,
      aes(wind_speed, power_scaled, col = "generic power curve")
    ) +
    scale_color_manual(
      values = col_palette
    ) +
    theme(
      legend.position = "bottom"
    ) +
    labs(col = "", x = "wind speed", y = "potential generation (MW)")

  p_nout <- pwr_curv_1wf %>%
    filter(between(halfHourEndTime, t0, t1), is.na(outageCapacity)) %>%
    ggplot(aes(ws_h, potential, col = "observed")) +
    geom_point(alpha = 0.2) +
    geom_line(
      data = scaled_pc,
      aes(wind_speed, power_scaled, col = "generic power curve")
    ) +
    scale_color_manual(
      values = col_palette
    ) +
    theme(
      legend.position = "bottom"
    ) +
    labs(col = "", x = "wind speed", y = "potential generation (MW)")
  invisible(list(
    data = pwr_curv_1wf,
    p_quant = p_quant,
    p_pot = p_pot,
    p_out = p_out,
    p_nout = p_nout
  ))
}


est_pwr_curv <- function(
  df,
  ws_col = "ws_h_wmean",
  power_col = "norm_potential",
  avg_fun = mean,
  n_bins = 20,
  quantile_bins = FALSE,
  plot = FALSE
) {
  # Check inputs
  # browser()
  stopifnot(is.data.frame(df))
  stopifnot(all(c(ws_col, power_col) %in% names(df)))

  ws <- df[[ws_col]]
  power <- df[[power_col]]

  # Create bins
  if (quantile_bins) {
    # quantile-based bins (equal number of points)
    bins <- cut(
      ws,
      breaks = quantile(
        ws,
        probs = seq(0, 1, length.out = n_bins + 1),
        na.rm = TRUE
      ),
      include.lowest = TRUE
    )
  } else {
    # equal-width bins (default)
    bins <- cut(ws, n_bins, include.lowest = TRUE)
  }

  # Compute mean wind speed & mean power in each bin
  curve_df <- aggregate(
    cbind(ws, power),
    by = list(bin = bins),
    FUN = avg_fun,
    na.rm = TRUE
  )

  # Rename nicely
  names(curve_df) <- c("bin", "ws_mean", "power_mean")

  # Sort by wind speed
  curve_df <- curve_df[order(curve_df$ws_mean), ]

  # Optionally plot
  if (plot) {
    plot(
      curve_df$ws_mean,
      curve_df$power_mean,
      type = "b",
      pch = 19,
      xlab = "Wind Speed (m/s)",
      ylab = "Expected Power",
      main = "Binned Power Curve Estimate"
    )
  }

  return(curve_df)
}


bru_fitted_exclude <- function(bru_fit, data, exclude = NULL) {
  # bru_fit : inlabru model object
  # data    : dataset used to fit the model
  # exclude : vector of component names to remove (e.g. "u")

  # 1. Start with fixed effects ---------------------------------------------
  lp <- rep(0, nrow(data))

  if (!is.null(bru_fit$summary.fixed)) {
    for (nm in rownames(bru_fit$summary.fixed)) {
      if (nm %in% names(data)) {
        # classic fixed slope
        lp <- lp + bru_fit$summary.fixed[nm, "mean"] * data[[nm]]
      } else if (nm == "Intercept") {
        lp <- lp + bru_fit$summary.fixed[nm, "mean"]
      }
    }
  }

  # 2. Add random/smooth components -----------------------------------------
  comp_list <- bru_fit$summary.random

  for (comp_name in names(comp_list)) {
    if (comp_name %in% exclude) {
      next
    } # skip excluded components

    comp_summary <- comp_list[[comp_name]]

    # Component ID column may be name "ID" or the variable name
    id_col <- if ("ID" %in% names(comp_summary)) {
      "ID"
    } else {
      names(comp_summary)[1]
    }

    # Try to match to data column
    if (comp_name %in% names(data)) {
      idx <- match(data[[comp_name]], comp_summary[[id_col]])
      lp <- lp + comp_summary$mean[idx]
    } else {
      # Try to guess replicated/random slope structure
      # e.g., tech_power replicated on norm_power_est0
      if (!is.null(data[[comp_name]])) {
        idx <- match(data[[comp_name]], comp_summary[[id_col]])
        lp <- lp + comp_summary$mean[idx]
      }
    }
  }

  return(lp)
}


# Function for multiple observations per site
spatial_corr_by_distance <- function(
  df,
  value_col,
  lon_col = "lon",
  lat_col = "lat",
  time_col = "halfHourEndTime",
  n_bins = 10,
  bin_type = c("equal", "quantile")
) {
  bin_type <- match.arg(bin_type)

  # Unique sites
  sites <- df %>%
    select(site = site_name, lon = {{ lon_col }}, lat = {{ lat_col }}) %>%
    distinct()

  # All pairwise site combinations
  site_pairs <- expand.grid(i = 1:nrow(sites), j = 1:nrow(sites)) %>%
    filter(i < j) %>%
    mutate(
      lon1 = sites$lon[i],
      lat1 = sites$lat[i],
      lon2 = sites$lon[j],
      lat2 = sites$lat[j],
      site1 = sites$site[i],
      site2 = sites$site[j],
      dist_km = geosphere::distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)) /
        1000
    )

  # Compute correlation per pair, matching by time
  cor_list <- site_pairs %>%
    rowwise() %>%
    mutate(
      corr = {
        x <- df %>%
          filter(site_name == site1) %>%
          select({{ time_col }}, {{ value_col }})
        y <- df %>%
          filter(site_name == site2) %>%
          select({{ time_col }}, {{ value_col }})
        paired <- inner_join(x, y, by = as.character(time_col))
        # Use complete cases to avoid missing values
        if (nrow(paired) < 2) NA else cor(paired[[2]], paired[[3]])
      }
    ) %>%
    ungroup()

  # Bin distances
  if (bin_type == "equal") {
    cor_list <- cor_list %>% mutate(bin = cut(dist_km, breaks = n_bins))
  } else {
    cor_list <- cor_list %>%
      mutate(
        bin = cut(
          dist_km,
          breaks = quantile(
            dist_km,
            probs = seq(0, 1, length.out = n_bins + 1),
            na.rm = TRUE
          ),
          include.lowest = TRUE
        )
      )
  }

  # Aggregate correlations by distance bin
  df_corr <- cor_list %>%
    group_by(bin) %>%
    summarise(
      dist_mean = mean(dist_km),
      corr_mean = mean(corr, na.rm = TRUE),
      n_pairs = n(),
      .groups = "drop"
    )

  return(df_corr)
}


spatial_corr_by_distance_time <- function(
  df,
  value_col,
  lon_col = "lon",
  lat_col = "lat",
  time_col = "halfHourEndTime",
  n_bins = 10,
  bin_type = c("equal", "quantile")
) {
  bin_type <- match.arg(bin_type)

  # Unique sites
  sites <- df %>%
    select(site = site_name, lon = {{ lon_col }}, lat = {{ lat_col }}) %>%
    distinct()

  # All pairwise site combinations
  site_pairs <- expand.grid(i = 1:nrow(sites), j = 1:nrow(sites)) %>%
    filter(i < j) %>%
    mutate(
      lon1 = sites$lon[i],
      lat1 = sites$lat[i],
      lon2 = sites$lon[j],
      lat2 = sites$lat[j],
      site1 = sites$site[i],
      site2 = sites$site[j],
      dist_km = distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)) / 1000
    )

  # Compute correlation per pair, matching by time
  cor_list <- site_pairs %>%
    rowwise() %>%
    mutate(
      corr = {
        x <- df %>%
          filter(site_name == site1) %>%
          select({{ time_col }}, {{ value_col }})
        y <- df %>%
          filter(site_name == site2) %>%
          select({{ time_col }}, {{ value_col }})
        paired <- inner_join(x, y, by = as.character(time_col))
        # Use complete cases to avoid missing values
        if (nrow(paired) < 2) NA else cor(paired[[2]], paired[[3]])
      }
    ) %>%
    ungroup()

  # Bin distances
  if (bin_type == "equal") {
    cor_list <- cor_list %>% mutate(bin = cut(dist_km, breaks = n_bins))
  } else {
    cor_list <- cor_list %>%
      mutate(
        bin = cut(
          dist_km,
          breaks = quantile(
            dist_km,
            probs = seq(0, 1, length.out = n_bins + 1),
            na.rm = TRUE
          ),
          include.lowest = TRUE
        )
      )
  }

  # Aggregate correlations by distance bin
  df_corr <- cor_list %>%
    group_by(bin) %>%
    summarise(
      dist_mean = mean(dist_km),
      corr_mean = mean(corr, na.rm = TRUE),
      n_pairs = n(),
      .groups = "drop"
    )

  return(df_corr)
}

spatial_corr_by_distance_fast <- function(
  df,
  value_col,
  lon_col = "lon",
  lat_col = "lat",
  time_col = "halfHourEndTime",
  n_bins = 10,
  bin_type = c("equal", "quantile")
) {
  bin_type <- match.arg(bin_type)

  # 1. Pivot data to wide format: each site is a column, rows are timestamps
  df_wide <- df %>%
    select(site_name, {{ time_col }}, {{ value_col }}) %>%
    pivot_wider(names_from = site_name, values_from = {{ value_col }})

  # 2. Compute correlation matrix (pairwise correlations ignoring NAs)
  cor_mat <- cor(df_wide[-1], use = "pairwise.complete.obs") # remove time_col column

  # 3. Get site coordinates
  sites <- df %>%
    select(site = site_name, lon = {{ lon_col }}, lat = {{ lat_col }}) %>%
    distinct()

  # 4. Get pairwise distances and correlations
  site_names <- sites$site
  pair_idx <- which(upper.tri(cor_mat), arr.ind = TRUE)

  cor_list <- data.frame(
    site1 = site_names[pair_idx[, 1]],
    site2 = site_names[pair_idx[, 2]],
    corr = cor_mat[upper.tri(cor_mat)]
  ) %>%
    left_join(sites, by = c("site1" = "site")) %>%
    rename(lon1 = lon, lat1 = lat) %>%
    left_join(sites, by = c("site2" = "site")) %>%
    rename(lon2 = lon, lat2 = lat) %>%
    mutate(dist_km = distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)) / 1000)

  # 5. Bin distances
  if (bin_type == "equal") {
    cor_list <- cor_list %>% mutate(bin = cut(dist_km, breaks = n_bins))
  } else {
    cor_list <- cor_list %>%
      mutate(
        bin = cut(
          dist_km,
          breaks = quantile(
            dist_km,
            probs = seq(0, 1, length.out = n_bins + 1),
            na.rm = TRUE
          ),
          include.lowest = TRUE
        )
      )
  }

  # 6. Aggregate correlations by distance bin
  df_corr <- cor_list %>%
    group_by(bin) %>%
    summarise(
      dist_mean = mean(dist_km),
      corr_mean = mean(corr, na.rm = TRUE),
      corr_lower = quantile(corr, 0.25, na.rm = TRUE),
      corr_upper = quantile(corr, 0.75, na.rm = TRUE),
      n_pairs = n(),
      .groups = "drop"
    )

  return(df_corr)
}

spatial_corr_by_distance_fast <- function(
  df,
  value_col,
  lon_col = "lon",
  lat_col = "lat",
  time_col = "halfHourEndTime",
  n_bins = 10,
  bin_type = c("equal", "quantile")
) {
  bin_type <- match.arg(bin_type)

  # 1. Pivot data to wide format: each site is a column, rows are timestamps
  df_wide <- df %>%
    select(site_name, {{ time_col }}, {{ value_col }}) %>%
    pivot_wider(names_from = site_name, values_from = {{ value_col }})

  site_names <- colnames(df_wide)[-1] # all site columns

  # 2. Correlation matrix
  cor_mat <- cor(df_wide[-1], use = "pairwise.complete.obs")

  # 3. n_obs_per_pair: number of overlapping non-missing observations
  n_obs_mat <- outer(
    site_names,
    site_names,
    Vectorize(function(x, y) {
      sum(!is.na(df_wide[[x]]) & !is.na(df_wide[[y]]))
    })
  )

  # 4. Get site coordinates
  sites <- df %>%
    select(site = site_name, lon = {{ lon_col }}, lat = {{ lat_col }}) %>%
    distinct()

  # 5. Build cor_list with distances and n_obs_per_pair
  pair_idx <- which(upper.tri(cor_mat), arr.ind = TRUE)

  cor_list <- data.frame(
    site1 = site_names[pair_idx[, 1]],
    site2 = site_names[pair_idx[, 2]],
    corr = cor_mat[upper.tri(cor_mat)],
    n_obs_per_pair = n_obs_mat[upper.tri(n_obs_mat)]
  ) %>%
    left_join(sites, by = c("site1" = "site")) %>%
    rename(lon1 = lon, lat1 = lat) %>%
    left_join(sites, by = c("site2" = "site")) %>%
    rename(lon2 = lon, lat2 = lat) %>%
    mutate(dist_km = distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)) / 1000)

  # 6. Bin distances
  if (bin_type == "equal") {
    cor_list <- cor_list %>% mutate(bin = cut(dist_km, breaks = n_bins))
  } else {
    cor_list <- cor_list %>%
      mutate(
        bin = cut(
          dist_km,
          breaks = quantile(
            dist_km,
            probs = seq(0, 1, length.out = n_bins + 1),
            na.rm = TRUE
          ),
          include.lowest = TRUE
        )
      )
  }

  # 7. Aggregate correlations by distance bin
  df_corr <- cor_list %>%
    group_by(bin) %>%
    summarise(
      dist_mean = mean(dist_km),
      corr_mean = mean(corr, na.rm = TRUE),
      corr_lower = quantile(corr, 0.25, na.rm = TRUE),
      corr_upper = quantile(corr, 0.75, na.rm = TRUE),
      n_pairs = n(),
      n_obs_mean = mean(n_obs_per_pair, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # 8. Compute approximate significance thresholds under null hypothesis
    mutate(
      se_null = 1 / sqrt(n_obs_mean - 3),
      threshold_95 = 1.96 * se_null
    )

  return(df_corr)
}


plot.hyper.dens <- function(
  inla_model,
  facet.cols = 2,
  logx = FALSE,
  ...
) {
  # Extract the hyperparameter densities into a data frame for plotting
  # browser()

  # if(inherits(inla_model,"bru")){
  #   marginals <- inla_model$internal.margina
  # }

  hyperpar_densities <- do.call(
    rbind,
    lapply(names(inla_model$marginals.hyperpar), function(param) {
      data.frame(
        x = inla_model$marginals.hyperpar[[param]][, "x"],
        y = inla_model$marginals.hyperpar[[param]][, "y"],
        parameter = param
      )
    })
    # mapply(function(param, logtrans) {
    #   data.frame(
    #     x = inla_model$marginals.hyperpar[[param]][, "x"],
    #     y = inla_model$marginals.hyperpar[[param]][, "y"],
    #     parameter = param
    #   ) #%>%
    #     # mutate(logx = ifelse(logtrans, log(x), x))
    # },
    # names(inla_model$marginals.hyperpar),
    # logx,
    # SIMPLIFY = FALSE
    # )
  )
  # browser()
  # hyperpar_densities <- hyperpar_densities %>%
  #   mutate(
  #     x = ifelse(logx, log(x), x)
  #   )

  # Plot the densities using ggplot2
  p.dens <- ggplot(hyperpar_densities, aes(x = x, y = y)) +
    geom_line() +
    facet_wrap(~parameter, scales = "free", ncol = facet.cols) +
    {
      if (logx) {
        scale_x_log10() # Apply log scale to the x-axis
      } else {
        NULL
      }
    } +
    labs(
      title = "Density of Hyperparameters",
      x = "Value",
      y = "Density"
    ) +
    theme_minimal()

  print(p.dens)
}
