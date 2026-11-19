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

bmu_list_fname <- "data/wind_bmu_2.csv.gz"
if (!file.exists(bmu_list_fname)) {
  cat("BMU list file not found: ", bmu_list_fname, "\n")
  cat("Recreating...\n")
  # wind.bmus <- bmus %>%
  #   filter(grepl("WIND",fuelType), fpnFlag) %>%
  #   filter(productionOrConsumptionFlag %in% c("P"))
  url <- "https://data.elexon.co.uk/bmrs/api/v1/reference/bmunits/all"
  response <- GET(url, accept("text/plain"))
  json_data <- content(response, "text", encoding = "UTF-8")
  bmus <- fromJSON(json_data)
  wind.bmus <- bmus %>%
    filter(
      grepl("WIND", fuelType) |
        grepl("wind|offshore", leadPartyName, ignore.case = TRUE) |
        grepl("wind farm|windfarm", bmUnitName, ignore.case = TRUE)
    ) %>%
    filter(!is.na(elexonBmUnit)) %>%
    mutate(generationCapacity = as.numeric(generationCapacity))
  # wind.bmus

  wind.bmus %>%
    write.csv(
      .,
      gzfile(bmu_list_fname),
    )
} else {
  wind.bmus <- fread(bmu_list_fname)
}

wind_bmus_alt_fname <- "data/wind_bmu_alt.csv"
if (!file.exists(wind_bmus_alt_fname)) {
  power_park <- readxl::read_xlsx(
    "data/power-park-modules.xlsx",
    sheet = 1,
    n_max = 252
  ) %>%
    rename(
      elexonBmUnit = `Settlement BMU name`,
      capacity_alt = `BMU or Large Reg Cap only`
    )
  names(power_park) <- tolower(gsub(" ", "_", names(power_park)))
  wind.bmus.alt <- wind.bmus %>%
    left_join(
      power_park %>% dplyr::select(elexonbmunit, capacity_alt),
      by = c("elexonBmUnit" = "elexonbmunit")
    ) %>%
    mutate(
      diff_cap = abs(generationCapacity - capacity_alt),
      capacity = case_when(
        is.na(capacity_alt) ~ generationCapacity,
        capacity_alt <= 0 ~ generationCapacity,
        diff_cap > generationCapacity * 0.1 ~ capacity_alt,
        TRUE ~ generationCapacity
      )
    )
} else {
  wind.bmus.alt <- read.csv(wind_bmus_alt_fname)
}

mypalette <- pal_aaas()(3)[c(1, 3)] # offshore / onshore
regime_palette <- pal_lancet()(3)[c(2, 3, 1)]
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
        dplyr::select(
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

# spatial_corr_by_distance_fast <- function(
#   df,
#   value_col,
#   lon_col = "lon",
#   lat_col = "lat",
#   time_col = "halfHourEndTime",
#   n_bins = 10,
#   bin_type = c("equal", "quantile")
# ) {
#   bin_type <- match.arg(bin_type)

#   # 1. Pivot data to wide format: each site is a column, rows are timestamps
#   df_wide <- df %>%
#     select(site_name, {{ time_col }}, {{ value_col }}) %>%
#     pivot_wider(names_from = site_name, values_from = {{ value_col }})

#   # 2. Compute correlation matrix (pairwise correlations ignoring NAs)
#   cor_mat <- cor(df_wide[-1], use = "pairwise.complete.obs") # remove time_col column

#   # 3. Get site coordinates
#   sites <- df %>%
#     select(site = site_name, lon = {{ lon_col }}, lat = {{ lat_col }}) %>%
#     distinct()

#   # 4. Get pairwise distances and correlations
#   site_names <- sites$site
#   pair_idx <- which(upper.tri(cor_mat), arr.ind = TRUE)

#   cor_list <- data.frame(
#     site1 = site_names[pair_idx[, 1]],
#     site2 = site_names[pair_idx[, 2]],
#     corr = cor_mat[upper.tri(cor_mat)]
#   ) %>%
#     left_join(sites, by = c("site1" = "site")) %>%
#     rename(lon1 = lon, lat1 = lat) %>%
#     left_join(sites, by = c("site2" = "site")) %>%
#     rename(lon2 = lon, lat2 = lat) %>%
#     mutate(dist_km = distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)) / 1000)

#   # 5. Bin distances
#   if (bin_type == "equal") {
#     cor_list <- cor_list %>% mutate(bin = cut(dist_km, breaks = n_bins))
#   } else {
#     cor_list <- cor_list %>%
#       mutate(
#         bin = cut(
#           dist_km,
#           breaks = quantile(
#             dist_km,
#             probs = seq(0, 1, length.out = n_bins + 1),
#             na.rm = TRUE
#           ),
#           include.lowest = TRUE
#         )
#       )
#   }

#   # 6. Aggregate correlations by distance bin
#   df_corr <- cor_list %>%
#     group_by(bin) %>%
#     summarise(
#       dist_mean = mean(dist_km),
#       corr_mean = mean(corr, na.rm = TRUE),
#       corr_lower = quantile(corr, 0.25, na.rm = TRUE),
#       corr_upper = quantile(corr, 0.75, na.rm = TRUE),
#       n_pairs = n(),
#       .groups = "drop"
#     )

#   return(df_corr)
# }

spatial_corr_by_distance_fast <- function(
  df,
  value_col,
  lon_col = "lon",
  lat_col = "lat",
  time_col = "halfHourEndTime",
  n_bins = 10,
  bin_type = c("equal", "quantile"),
  bands = c(0.25, 0.75)
) {
  bin_type <- match.arg(bin_type)

  # 1. Pivot data to wide format: each site is a column, rows are timestamps
  df_wide <- df %>%
    dplyr::select(site_name, {{ time_col }}, {{ value_col }}) %>%
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
    dplyr::select(
      site = site_name,
      lon = {{ lon_col }},
      lat = {{ lat_col }}
    ) %>%
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
      corr_lower = quantile(corr, bands[1], na.rm = TRUE),
      corr_upper = quantile(corr, bands[2], na.rm = TRUE),
      n_pairs = n(),
      n_obs_mean = mean(n_obs_per_pair, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # 8. Compute approximate significance thresholds under null hypothesis
    mutate(
      se_null = 1 / sqrt(n_obs_mean - 3),
      threshold_95 = 1.96 * se_null
    )

  return(list(
    summary = df_corr,
    pairs = cor_list
  ))
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


plot.effects <- function(
  inla.model,
  rand.effect,
  trans = \(x) x,
  show.fig = TRUE,
  window = "all",
  n.replicate = 1,
  replicate_names = NULL,
  id_override = NULL,
  ...
) {
  if (rand.effect == "etaderiv" | rand.effect == "eta+deriv") {
    lagged_eta <- inla.model$summary.random$eta.2[
      -nrow(inla.model$summary.random$eta.2),
      c(1, 2, 3, 6, 5, 4, 7, 8)
    ]
    spline_summary <- inla.model$summary.random$eta.1[-1, ] -
      lagged_eta
    # ID fix
    spline_summary$ID <- inla.model$summary.random$eta.1$ID[-1]

    if (rand.effect == "eta+deriv") {
      spline_summary[, 2:ncol(spline_summary)] <- spline_summary[,
        2:ncol(spline_summary)
      ] +
        inla.model$summary.random$eta[-1, 2:ncol(spline_summary)]
    }
  } else {
    # Extract the summary for the specified random effect
    spline_summary <- inla.model$summary.random[[rand.effect]]
  }

  if (window != "all" & is.numeric(window)) {
    spline_summary <- spline_summary %>% tail(trunc(window))
  }
  n_row <- nrow(spline_summary)
  # browser()
  # Create a data frame for ggplot
  plot_data <- data.frame(
    group = if (!is.null(id_override)) {
      id_override
    } else {
      spline_summary$ID
    },
    mean = spline_summary$mean, # Posterior mean
    lower = spline_summary$`0.025quant`, # 2.5% quantile (lower credible interval)
    upper = spline_summary$`0.975quant` # 97.5% quantile (upper credible interval)
  ) %>%
    {
      if (rand.effect %in% c("hour", "month")) {
        mutate(., across(group, as.numeric))
      } else {
        .
      }
    } %>%
    {
      # browser()
      if (n.replicate > 1) {
        mutate(
          .,
          type = rep(
            if (is.null(replicate_names)) {
              paste0("rep", 1:n.replicate)
            } else {
              replicate_names
            },
            each = n_row / n.replicate
          )
        )
      } else {
        .
      }
    } %>%
    mutate(across(mean:upper, trans))

  # Plot using ggplot2
  layers <-
    if (rand.effect %in% c("techno", "slope")) {
      list(
        geom_col(aes(y = mean), fill = blues9[7], width = 0.7),
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)
      )
    } else {
      list(
        geom_line(aes(y = mean), color = "blue", linewidth = 1),
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2)
      )
    }
  p1 <- ggplot(plot_data, aes(x = group)) +
    layers +
    labs(
      title = paste("Estimated effect for", rand.effect),
      x = paste(rand.effect),
      y = "Estimated Effect"
    ) +
    {
      if (n.replicate > 1) facet_wrap(~type, scales = "free_y")
    } +
    theme_minimal()
  if (show.fig) {
    print(p1)
  }
  invisible(list(data = plot_data, fig = p1))
}


mape <- function(actual, predicted, pos_only = FALSE, ...) {
  if (pos_only) {
    idx <- actual > 0 & !is.na(actual) & !is.na(predicted)
    return(mean(abs((actual[idx] - predicted[idx]) / actual[idx]), ...) * 100)
  } else {
    return(
      mean(
        abs((actual - predicted) / actual) * 100,
        ...
      )
    )
  }
}
mdape <- function(actual, predicted, pos_only = FALSE, ...) {
  if (pos_only) {
    idx <- actual > 0 & !is.na(actual) & !is.na(predicted)
    return(median(abs((actual[idx] - predicted[idx]) / actual[idx]), ...) * 100)
  } else {
    return(
      median(
        abs((actual - predicted) / actual) * 100,
        ...
      )
    )
  }
}


get_bru_formula <- function(mod) {
  # get linear predictor formula
  form_raw <- mod$.args$formula
  lin_pred <- str_extract_all(as.character(form_raw)[3], "(?<=f\\()[^,]+")[[
    1
  ]] %>%
    paste(collapse = " + ")
  return(lin_pred)
}


bru_ci_plot <- function(
  bru_model,
  newdata,
  n.samples = 100,
  show.fig = TRUE,
  alphas = c(0.025, 0.975),
  oos_type = "time",
  family = "gaussian",
  t_start = NULL
) {
  family <- match.arg(family, c("gaussian", "beta", "t"))
  n <- nrow(newdata)
  stopifnot(
    n > 0L,
    length(n.samples) == 1L,
    is.finite(n.samples),
    n.samples >= 1L,
    n.samples == floor(n.samples),
    all(is.finite(alphas) & alphas > 0 & alphas < 1)
  )

  if (is.null(t_start)) {
    t_start <- t1
  }
  end_24 <- t_start + 24 * 3600
  end_12 <- t_start + 12 * 3600

  lin_pred <- get_bru_formula(bru_model)
  hyper <- bru_model$.args$control.family[[1L]]$hyper
  hyper_name <- function(parameter) {
    gsub("[ -]+", "_", parameter$output.name)
  }
  # Keep the original parameter convention for comparisons with existing runs.
  prec_expr <- if (isTRUE(hyper[[1L]]$fixed)) {
    as.character(hyper[[1L]]$initial)
  } else {
    hyper_name(hyper[[1L]])
  }
  dof_expr <- if (family == "t") {
    paste0(", dof = ", hyper_name(hyper[[2L]]))
  } else {
    ""
  }
  formula_draw <- stats::as.formula(
    sprintf(
      "~ data.frame(lp = (%s), prec = (%s)%s)",
      lin_pred,
      prec_expr,
      dof_expr
    )
  )

  keys <- c("coord_id", "site_name", "time")
  meta <- tibble::tibble(
    coord_id = newdata$coord_id,
    site_name = newdata$site_name,
    time = newdata$time,
    anomaly = newdata$anomaly,
    norm_potential = newdata$norm_potential,
    norm_potential_orig = newdata$norm_potential_orig
  )
  observed <- meta$norm_potential
  valid <- which(!is.na(observed))
  capacity <- newdata$capacity[valid]
  gb_time <- sort(unique(meta$time[valid]), na.last = TRUE)
  time_id <- match(meta$time[valid], gb_time)
  sum_by_time <- function(x) {
    if (!length(x)) {
      return(numeric())
    }
    as.numeric(rowsum(
      matrix(x, ncol = 1L),
      time_id,
      reorder = TRUE,
      na.rm = TRUE
    ))
  }
  gb_capacity <- sum_by_time(capacity)
  gb_denominator <- gb_capacity
  gb_denominator[sum_by_time(as.integer(is.na(capacity))) > 0L] <- NA_real_
  potential <- sum_by_time(observed[valid] * capacity)
  potential_orig <- sum_by_time(meta$norm_potential_orig[valid] * capacity)

  probs <- c(0.025, 0.975)
  alphas <- rev(alphas)
  coverage_names <- paste0("coverage_", round(100 * (1 - alphas), 1))
  if (anyDuplicated(coverage_names)) {
    stop("Supply alphas with distinct coverage levels after rounding to 0.1%.")
  }
  hits <- matrix(0L, nrow = n, ncol = length(alphas))
  cdf_count <- integer(n)
  lp_sum <- numeric(n)
  fits <- matrix(NA_real_, nrow = n, ncol = n.samples)
  gb_mw <- matrix(NA_real_, nrow = length(gb_time), ncol = n.samples)

  # One joint call keeps the posterior state aligned across sites and hours.
  draws <- inlabru::generate(
    bru_model,
    newdata = newdata,
    formula = formula_draw,
    n.samples = n.samples
  )
  for (s in seq_len(n.samples)) {
    draw <- draws[[s]]
    draws[s] <- list(NULL)
    lp <- draw$lp
    prec <- draw$prec
    lp_sum <- lp_sum + lp

    if (family == "beta") {
      shape1 <- stats::plogis(lp) * prec
      shape2 <- prec - shape1
      fit <- stats::rbeta(n, shape1, shape2)
      if (length(alphas)) {
        cdf <- stats::pbeta(observed, shape1, shape2)
      }
    } else if (family == "t") {
      sd <- 1 / sqrt(prec)
      fit <- lp + stats::rt(n, draw$dof) * sd
      if (length(alphas)) {
        cdf <- stats::pt((observed - lp) / sd, draw$dof)
      }
    } else {
      sd <- 1 / sqrt(prec)
      fit <- stats::rnorm(n, lp, sd)
      if (length(alphas)) {
        cdf <- stats::pnorm(observed, lp, sd)
      }
    }

    fits[, s] <- fit
    gb_mw[, s] <- sum_by_time(fit[valid] * capacity)
    if (length(alphas)) {
      usable <- !is.na(cdf)
      cdf_count <- cdf_count + usable
      for (j in seq_along(alphas)) {
        inside <- usable & cdf >= alphas[j] / 2 & cdf <= 1 - alphas[j] / 2
        hits[, j] <- hits[, j] + inside
      }
    }
  }
  rm(draws, draw)

  clip <- function(x) pmin(1, pmax(0, x))
  # Leave predictions on their original scale until reporting the summaries.
  gb_fits <- gb_mw / gb_denominator
  gb_quant <- matrixStats::rowQuantiles(
    gb_fits,
    probs = probs,
    useNames = FALSE
  )
  GB_summary <- tibble::tibble(
    time = gb_time,
    mean = clip(rowMeans(gb_fits)),
    lwr = clip(gb_quant[, 1L]),
    upr = clip(gb_quant[, 2L]),
    width = gb_quant[, 2L] - gb_quant[, 1L],
    norm_potential_orig = potential_orig / gb_denominator,
    norm_potential = potential / gb_denominator,
    capacity = gb_capacity,
    potential = potential,
    potential_orig = potential_orig,
    est_mw = rowMeans(gb_mw)
  )

  grouped <- dplyr::group_by(meta, coord_id, site_name, time)
  groups <- dplyr::group_rows(grouped)
  wf_summary <- dplyr::group_keys(grouped)
  duplicate_groups <- which(lengths(groups) > 1L)
  row_index <- if (length(duplicate_groups)) {
    vapply(groups, function(rows) rows[1L], integer(1L))
  } else {
    unlist(groups, use.names = FALSE)
  }
  wf_quant <- matrixStats::rowQuantiles(
    fits,
    rows = row_index,
    probs = probs,
    useNames = FALSE
  )
  wf_summary$lwr <- wf_quant[, 1L]
  wf_summary$upr <- wf_quant[, 2L]
  wf_summary$fit <- rowMeans(fits)[row_index]
  wf_summary$lp <- lp_sum[row_index] / n.samples
  wf_summary$anomaly <- meta$anomaly[row_index]
  wf_summary$norm_potential <- observed[row_index]
  wf_summary$norm_potential_orig <- meta$norm_potential_orig[row_index]

  # Support inputs with multiple rows for the same site and time.
  for (g in duplicate_groups) {
    rows <- groups[[g]]
    pooled <- fits[rows, , drop = FALSE]
    limits <- stats::quantile(pooled, probs, names = FALSE)
    wf_summary$lwr[g] <- limits[1L]
    wf_summary$upr[g] <- limits[2L]
    wf_summary$fit[g] <- mean(pooled)
    wf_summary$lp[g] <- mean(lp_sum[rows]) / n.samples
    wf_summary$anomaly[g] <- dplyr::first(meta$anomaly[rows], na_rm = TRUE)
    wf_summary$norm_potential[g] <- mean(observed[rows], na.rm = TRUE)
    wf_summary$norm_potential_orig[g] <- mean(
      meta$norm_potential_orig[rows],
      na.rm = TRUE
    )
  }
  wf_summary$width <- wf_summary$upr - wf_summary$lwr
  wf_summary <- dplyr::mutate(
    wf_summary,
    dplyr::across(c(lwr, upr, fit, lp), clip)
  )
  wf_summary <- wf_summary[, c(
    keys,
    "lwr",
    "upr",
    "fit",
    "lp",
    "width",
    "anomaly",
    "norm_potential",
    "norm_potential_orig"
  )]

  hours_24 <- meta$time >= t_start & meta$time < end_24
  hours_12 <- meta$time >= t_start & meta$time < end_12
  score_rows <- function(rows) {
    sample_mat <- fits[rows, , drop = FALSE]
    colnames(sample_mat) <- paste0("samp_", seq_len(n.samples))
    compute_scores(observed = observed[rows], sample_mat = sample_mat)
  }
  score_oos <- if (oos_type == "space") {
    meta$time < t_start
  } else {
    meta$time >= t_start
  }
  scores_24 <- scores_12 <- NULL
  score_indices <- which(score_oos & !is.na(observed))
  if (identical(score_indices, seq_len(n))) {
    colnames(fits) <- paste0("samp_", seq_len(n.samples))
    scores <- compute_scores(observed = observed, sample_mat = fits)
  } else {
    scores <- score_rows(score_indices)
  }
  if (oos_type == "time") {
    scores_24 <- score_rows(which(hours_24 & !is.na(observed)))
    scores_12 <- score_rows(which(hours_12 & !is.na(observed)))
  }
  rm(fits, score_rows)

  cov_24h <- cov_12h <- cov_gbl <- cov_time <- cov_loc <- NULL
  if (length(alphas)) {
    coverage <- hits[row_index, , drop = FALSE] / cdf_count[row_index]
    for (g in duplicate_groups) {
      rows <- groups[[g]]
      coverage[g, ] <- colSums(hits[rows, , drop = FALSE]) /
        sum(cdf_count[rows])
    }
    colnames(coverage) <- coverage_names
    coverage_df <- dplyr::bind_cols(
      wf_summary[, keys],
      tibble::as_tibble(coverage)
    )
    summarize_coverage <- function(x) {
      dplyr::summarise(
        x,
        dplyr::across(dplyr::all_of(coverage_names), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      )
    }
    if (oos_type == "time") {
      oos <- dplyr::filter(coverage_df, time >= t_start)
      cov_24h <- summarize_coverage(dplyr::filter(oos, time < end_24))
      cov_12h <- summarize_coverage(dplyr::filter(oos, time < end_12))
    } else {
      oos <- dplyr::filter(coverage_df, time <= t_start)
    }
    cov_time <- summarize_coverage(dplyr::group_by(oos, time))
    cov_loc <- summarize_coverage(dplyr::group_by(oos, coord_id, site_name))
    cov_gbl <- summarize_coverage(oos)
  }

  mean_width <- function(x) {
    rows <- which(x$time >= t_start & x$time < end_24)
    mean(x$width[rows], na.rm = TRUE)
  }
  band_width <- tibble::tibble(
    gbl = mean_width(GB_summary),
    wf = mean_width(wf_summary)
  )

  if (show.fig) {
    print(
      ggplot2::ggplot(GB_summary, ggplot2::aes(x = time)) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lwr, ymax = upr),
          fill = blues9[5L],
          alpha = 0.5
        ) +
        ggplot2::geom_line(
          ggplot2::aes(y = mean),
          colour = blues9[9L],
          linewidth = 1
        ) +
        ggplot2::geom_line(
          ggplot2::aes(y = norm_potential_orig),
          colour = "darkred",
          linewidth = 1
        )
    )
  }
  invisible(list(
    wf_summary = wf_summary,
    GB_summary = GB_summary,
    formula = lin_pred,
    cov_gbl = cov_gbl,
    cov_time = cov_time,
    cov_loc = cov_loc,
    cov_24h = cov_24h,
    cov_12h = cov_12h,
    scores = scores,
    scores_24 = scores_24,
    scores_12 = scores_12,
    band_width = band_width
  ))
}


available_cores <- function() {
  x <- system("taskset -pc $$", intern = TRUE)

  # extract CPU list after colon
  cpus <- sub(".*: ", "", x)

  # split ranges like 19-22 or 1,19,20,21
  parts <- strsplit(cpus, ",")[[1]]

  total <- 0
  for (p in parts) {
    if (grepl("-", p)) {
      r <- as.integer(strsplit(p, "-")[[1]])
      total <- total + (r[2] - r[1] + 1)
    } else {
      total <- total + 1
    }
  }

  total
}

extract_score_model <- function(mod.obj) {
  # browser()
  # check if model has two likelihoods
  # n <- if (length(mod.obj$.args$family) > 1) length(mod.obj$.args$data$intercept) / 2 else NULL

  if (length(mod.obj$.args$family) > 1) {
    n <- length(mod.obj$.args$data$intercept) / 2
    result <- data.frame(
      mlik = mod.obj$mlik[2],
      dic = if (!is.null(mod.obj$dic)) {
        sum(mod.obj$dic$dic[-c(1:n)], na.rm = TRUE)
      } else {
        NA
      },
      waic = if (!is.null(mod.obj$waic)) {
        sum(mod.obj$waic$local.waic[-c(1:n)], na.rm = TRUE)
      } else {
        NA
      },
      mean_log_cpo = if (!is.null(mod.obj$cpo$cpo)) {
        mean(log(mod.obj$cpo$cpo[-c(1:n)]), na.rm = TRUE)
      } else {
        NA
      },
      mean_log_gcpo = if (!is.null(mod.obj$gcpo$gcpo)) {
        mean(log(mod.obj$gcpo$gcpo[-c(1:n)]), na.rm = TRUE)
      } else {
        NA
      }
    )
  } else {
    result <- data.frame(
      mlik = mod.obj$mlik[2],
      dic = if (!is.null(mod.obj$dic)) mod.obj$dic$dic else NA,
      waic = if (!is.null(mod.obj$waic)) mod.obj$waic$waic else NA,
      mean_log_cpo = if (!is.null(mod.obj$cpo)) {
        mean(log(mod.obj$cpo$cpo), na.rm = TRUE)
      } else {
        NA
      },
      mean_log_gcpo = if (!is.null(mod.obj$gcpo$gcpo)) {
        mean(log(mod.obj$gcpo$gcpo), na.rm = TRUE)
      } else {
        NA
      }
    )
  }
  return(result)
}


extract_pit_model <- function(mod.obj, data = NULL) {
  # browser()
  # check if model has two likelihoods
  # n <- if (length(mod.obj$.args$family) > 1) length(mod.obj$.args$data$intercept) / 2 else NULL

  # if (length(mod.obj$.args$family) > 1) {
  #   n <- length(mod.obj$.args$data$intercept) / 2
  #   pit_vals <- mod.obj$cpo$pit[-c(1:n)]
  # } else {
  #   n <- length(mod.obj$cpo$pit)
  #   pit_vals <- mod.obj$cpo$pit
  # }
  n <- length(mod.obj$cpo$pit)
  pit_vals <- mod.obj$cpo$pit

  if (!is.null(data)) {
    time <- data$time
    site_name <- data$site_name
    coord_id <- data$coord_id
  } else {
    time <- rep(NA, n)
    site_name <- rep(NA, n)
    coord_id <- rep(NA, n)
  }

  result <- data.frame(
    id = 1:n,
    time = time,
    site_name = site_name,
    coord_id = coord_id,
    pit = pit_vals
  )
  return(result)
}
base_bru_options <- bru_options(
  bru_verbose = 3,
  # verbose = TRUE,
  control.compute = list(
    dic = TRUE,
    cpo = TRUE,
    waic = TRUE,
    control.gcpo = list(enable = TRUE, num.level.sets = 3)
  )
)
fit_bru_att <- function(n_attempts = 3, options = base_bru_options, ...) {
  for (i in 1:n_attempts) {
    message("Attempt ", i)

    current_options <- options

    if (i >= 2) {
      current_options$verbose <- TRUE
    }

    fit <- tryCatch(
      bru(..., options = current_options),
      error = function(e) NULL
    )

    valid_fit <- !is.null(fit) &&
      !is.null(fit$bru_info) &&
      !is.null(fit$summary.fixed)

    if (valid_fit) {
      return(fit)
    } else {
      message("bru failed on attempt ", i, ". Retrying...")
    }
  }

  stop("bru failed after ", n_attempts, " attempts")
}


compute_scores <- function(
  observed,
  sample_mat,
  bs_thresholds = c(0.01, 0.05, 0.1, 0.2)
) {
  scores <- list()
  # CRPS
  crps <- scoringRules::crps_sample(
    y = observed, # day-ahead response
    dat = sample_mat
  )
  scores$crps <- mean(crps)
  # Energy
  energy.s <- scoringRules::es_sample(
    y = observed, # day-ahead response
    dat = sample_mat
  )
  scores$energy <- energy.s
  # log
  log.s <- scoringRules::logs_sample(
    y = observed, # day-ahead response
    dat = sample_mat
  )
  scores$log <- mean(log.s[is.finite(log.s)])

  # Brier score

  bs_vec <- lapply(
    bs_thresholds,
    function(thresh) {
      prob_low <- rowMeans(sample_mat <= thresh)
      event_low <- observed <= thresh
      brier <- mean((prob_low - event_low)^2)
    }
  ) %>%
    unlist() %>%
    setNames(paste0("bs_", gsub("\\.", "_", bs_thresholds)))
  scores$brier <- bs_vec

  # scores$log_obs <- log.s
  scores
}


# make_groups <- function(
#   df,
#   ws_breaks,
#   pow_breaks,
#   d_coast_breaks,
#   elev_breaks
# ) {
#   group_left <- function(x, breaks) {
#     i <- findInterval(x, breaks)
#     i <- pmax(1, pmin(i, length(breaks)))
#     breaks[i]
#   }

#   df %>%
#     mutate(
#       ws_group = group_left(
#         ws_h,
#         ws_breaks
#       ),

#       pow_group = group_left(
#         norm_power_est0,
#         pow_breaks
#       ),

#       d_coast_group = group_left(
#         dist_coast,
#         d_coast_breaks
#       ),

#       elev_group = group_left(
#         elevation,
#         elev_breaks
#       )
#     )
# }

make_group_breaks <- function(x, n = 20) {
  # Same quantile probabilities as inla.group()
  breaks <- unique(
    quantile(
      x,
      probs = c(0, ppoints(n - 1), 1),
      na.rm = TRUE
    )
  )

  # Assign training observations to intervals
  idx <- cut(
    x,
    breaks = breaks,
    include.lowest = TRUE
  )

  # Same group representative as inla.group():
  # median of the observations within each interval
  group_values <- tapply(
    x,
    idx,
    median,
    na.rm = TRUE
  )

  list(
    breaks = as.numeric(breaks),
    values = as.numeric(group_values)
  )
}

make_groups <- function(
  df,
  ws_groups,
  pow_groups,
  d_coast_groups,
  elev_groups
) {
  group_value <- function(x, groups) {
    idx <- findInterval(
      x,
      groups$breaks,
      left.open = FALSE
    )

    # Keep observations outside the training range
    # in the first/last group
    idx <- pmax(
      1,
      pmin(idx, length(groups$values))
    )

    groups$values[idx]
  }

  df %>%
    mutate(
      ws_group = group_value(
        ws_h,
        ws_groups
      ),

      pow_group = group_value(
        norm_power_est0,
        pow_groups
      ),

      d_coast_group = group_value(
        dist_coast,
        d_coast_groups
      ),

      elev_group = group_value(
        elevation,
        elev_groups
      )
    )
}


cov_fig <- function(
  data,
  h_max = NULL,
  show.fig = TRUE,
  name_prefix = "WF"
) {
  # browser()
  data <- data %>%
    filter(oos) %>%
    group_by(model) %>%
    mutate(
      t1 = first(time),
      h = difftime(time, t1, units = "hours") %>% as.numeric()
    )

  if (is.null(h_max)) {
    h_max <- max(data$h)
  }

  for (j in seq_along(h_max)) {
    cov_bands_hmax <- data %>%
      filter(h <= h_max[j]) %>%
      filter(!model %in% excluded_models) %>%
      group_by(model) %>%
      summarise(
        coverage = mean(
          norm_potential >= lwr & norm_potential <= upr,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      arrange(desc(coverage)) %>%
      mutate(
        model = factor(model, levels = model)
      )

    cov_bar <- cov_bands_hmax %>%
      ggplot(aes(x = model, y = coverage)) +
      geom_col(fill = blues9[7]) +
      geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkred") +
      geom_text(aes(label = round(coverage, 3)), vjust = -0.5) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Model", y = "Coverage") +
      scale_x_discrete(labels = mod_labels) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    if (show.fig) {
      print(cov_bar)
    }
    ggsave(
      filename = sprintf(
        "fig/%s/%s_pred_band_coverage_%dh.pdf",
        batch_name,
        name_prefix,
        h_max[j]
      ),
      plot = cov_bar,
      width = 10,
      height = 6,
      # dpi = 300
    )
  }
}

sum_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
}


report_memory <- function(label = "") {
  cat("\n--- Memory report:", label, "---\n")

  print(gc())

  status <- readLines("/proc/self/status")

  cat(
    grep(
      "^(VmRSS|VmHWM|VmSize|VmPeak):",
      status,
      value = TRUE
    ),
    sep = "\n"
  )

  cat("\n")
}

sizes <- sapply(
  ls(envir = .GlobalEnv),
  function(x) {
    as.numeric(
      object.size(get(x, envir = .GlobalEnv))
    )
  }
)
