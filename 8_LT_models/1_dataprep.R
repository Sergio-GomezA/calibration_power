starttime <- Sys.time()

require(tidyverse)
require(sf)
require(INLA)
require(inlabru)
require(fmesher)
require(ggspatial)
require(ModelMetrics)
require(qmap)
require(ggridges)
require(ggthemes)
require(ggsci)
require(arrow)

spde1d_bru_opt <- as.bru_options(
  list(
    control.family = list(
      hyper = list(
        prec = list(
          initial = prec_init_gau,
          fixed = fixed_gaus_1DSPE
        )
      )
    )
  )
)
scores_df <- list()
pit_list <- list()
model_list <- list()

cat(
  "-------------------------------------------------------------------------------------------------\n"
)
cat("Preparing data for model fitting LT series\n")
cat(
  "-------------------------------------------------------------------------------------------------\n"
)
# 1. data preparation ####

dat_starttime <- Sys.time()
cat("Preparing data for model fitting\n")

cat("Creating coordinate list\n")
coord_list_fname <- "data/coord_list.csv"
if (!file.exists(coord_list_fname)) {
  cat("Creating new coordinate list\n")
  set.seed(1)
  coord_list <- pwr_curv_df %>%
    distinct(coord_id, tech_typ) %>%
    arrange(coord_id)

  coord_samp <- coord_list %>%
    group_by(tech_typ) %>%
    slice_sample(prop = 0.8)

  coord_list <- coord_list %>%
    mutate(
      sampled = ifelse(coord_id %in% coord_samp$coord_id, TRUE, FALSE)
    )
  write.csv(
    coord_list,
    file = coord_list_fname,
    row.names = FALSE
  )
} else {
  cat("Loading existing coordinate list\n")
  coord_list <- read.csv(coord_list_fname)
}

coord_list <- coord_list %>%
  mutate(
    sampled_old = sampled,
    sampled = TRUE
  )
## 1.0.1 GB daily summary ####

gb_day_df_fname <- sprintf("data/GB_daily_summary.parquet")
gb_df_fname <- file.path(
  gen_path,
  sprintf("GB_aggr_An-%0.1f.parquet", norm_dist_tol)
)

if (
  !file.exists(gb_day_df_fname) ||
    !file.exists(gb_df_fname) ||
    recalculate_gb
) {
  cat("GB daily summary file not found, creating new summary\n")
  cat("rebuilding aggr df for GB level models\n")

  pwr_curv_df <- read_parquet(file.path(
    gen_path,
    "power_curve_all_enriched.parquet"
  ))

  GB_df <- pwr_curv_df %>%
    rename(time = halfHourEndTime) %>%
    filter(coord_id %in% coord_list$coord_id[coord_list$sampled]) %>%
    arrange(site_name) %>%
    group_by(lon, lat, time) %>%
    summarise(
      site_name = first(site_name),
      coord_id = first(coord_id),
      elevation = first(elevation),
      dist_coast = first(dist_coast),
      tech_typ = first(tech_typ),
      across(c(ws_h, wd10, wd100), mean),
      ws_h_wmean = sum(ws_h * capacity) / sum(capacity),
      across(c(potential, power_est0, capacity, curtailment), sum),
      .groups = "drop"
    ) %>%
    # filter(time >= "2023-01-01", time < "2025-01-01") %>%
    mutate(
      norm_potential = pmin(1, potential / capacity),
      norm_power_est0 = power_est0 / capacity,
      anomaly = case_when(
        # norm_potential <= tol & norm_power_est0 >= 0.05 ~ TRUE,
        # norm_power_est0 >= 1 - tol & norm_potential <= 0.5 ~ TRUE,
        abs(norm_power_est0 - norm_potential) >= norm_dist_tol ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    # filter(!anomaly) %>%
    group_by(tech_typ, time) %>%
    summarise(
      # ws_h_wmean = sum(ws_h * capacity),
      across(
        c(ws_h, wd10, wd100),
        ~ sum(. * capacity),
        .names = "{.col}_wmean"
      ),
      across(
        c(power_est0, potential, capacity),
        sum
      ),
      across(c(ws_h, wd10, wd100), mean, .names = "{.col}_mean")
    ) %>%
    mutate(
      across(
        c(power_est0, potential),
        ~ . / capacity,
        .names = "norm_{.col}"
      ),
      # ws_h_wmean = ws_h_wmean / capacity,
      across(matches("_wmean"), ~ . / capacity),
      month = factor(month(time)),
      hour = factor(hour(time))
    ) %>%
    mutate(ws_group = inla.group(ws_h_wmean, n = 20, method = "quantile")) %>%
    group_by(tech_typ) %>%
    arrange(time, .by_group = TRUE) %>%
    mutate(t = row_number()) %>%
    ungroup() %>%
    # rename(time = halfHourEndTime) %>%
    mutate(
      err = norm_power_est0 - norm_potential,
      error0 = norm_potential - norm_power_est0,
      date = as.Date(time)
    )

  arrow::write_parquet(
    GB_df,
    gb_df_fname
  )

  gb_day_df <- GB_df %>%
    group_by(date, tech_typ) %>%
    summarise(
      across(
        c(norm_power_est0, norm_potential),
        ~ sum(. * capacity) / sum(capacity)
      ),
      across(c(ws_h_wmean), ~ sum(. * capacity) / sum(capacity)),
      across(c(capacity), mean)
    ) %>%
    summarise(
      across(
        c(norm_power_est0, norm_potential),
        ~ sum(. * capacity) / sum(capacity)
      ),
      across(c(ws_h_wmean), ~ sum(. * capacity) / sum(capacity)),
      across(c(capacity), sum),
      .groups = "drop"
    )

  cutprobs3 <- c(0.25, 0.75)
  p_quant3 <- quantile(gb_day_df$norm_potential, probs = cutprobs3)
  cutprobs7 <- c(0.1, 0.2, 0.25, 0.75, 0.8, 0.9)
  p_quant7 <- quantile(gb_day_df$norm_potential, probs = cutprobs7)

  gb_day_df <- gb_day_df %>%
    mutate(
      p_group3 = cut(
        norm_potential,
        breaks = c(-Inf, p_quant3, Inf),
        labels = c("low", "mid", "high")
      ),
      p_group7 = cut(
        norm_potential,
        breaks = c(-Inf, p_quant7, Inf)
      )
    )

  write_parquet(gb_day_df, gb_day_df_fname)
  rm(pwr_curv_df)
  gc()
} else {
  cat("Loading existing GB daily summary\n")
  GB_df <- read_parquet(gb_df_fname)
  gb_day_df <- read_parquet(gb_day_df_fname)
}

d0 <- as.Date("2025-01-01")
d0_tag <- base::format(d0, "%y%m%d")

extension <- ifelse(local_run, local_ext, cluster_ext)
df_pattern <- sprintf("^calibration_df_.*_%s\\.%s$", d0_tag, extension)
files_found <- list.files(
  file.path(output_path, batch_name, "data/fit"),
  pattern = df_pattern,
  full.names = TRUE
)

if (!override_objects && length(files_found) > 0) {
  cat(
    "Calibration data file already exists for this day. Loading existing data.\n"
  )
  if (extension != "rds") {
    wf_dt <- st_read(files_found[1])
  } else {
    wf_dt <- readRDS(files_found[1])
  }
  setDT(wf_dt)
  wf_dt <- wf_dt[date >= d0 - n.days.before.fit]
} else {
  if (override_objects) {
    cat("File found but override option is active. Rebuilding.\n")
  } else {
    cat(
      "No existing calibration data file found for this day. Preparing new data.\n"
    )
  }

  pwr_curv_dt <- as.data.table(
    read_parquet(file.path(
      gen_path,
      "power_curve_all_enriched.parquet"
    ))
  )
  setnames(pwr_curv_dt, "halfHourEndTime", "time")
  n.days <- 0
  wf_dt <- pwr_curv_dt[
    time >= as.POSIXct(d0 - max(n.days.before.fit, n.days.before.light)) &
      time < as.POSIXct(d0 + n.days) &
      coord_id %in% coord_list$coord_id[coord_list$sampled]
  ]
  wf_dt[, elevation := pmax(0, elevation)]
  wf_dt[,
    site_name := trimws(
      gsub("\\b(wind\\s*farm|wf)\\b", "", site_name, ignore.case = TRUE)
    )
  ]
  wf_dt <- wf_dt[,
    .(
      lon = first(lon),
      lat = first(lat),
      site_name = first(site_name),
      elevation = first(elevation),
      dist_coast = first(dist_coast),
      tech_typ = first(tech_typ),

      ws_h = mean(ws_h, na.rm = TRUE),

      ws_h_wmean = weighted.mean(
        ws_h,
        capacity,
        na.rm = TRUE
      ),

      potential = sum_or_na(potential),
      power_est0 = sum_or_na(power_est0),
      capacity = sum_or_na(capacity),
      curtailment = sum_or_na(curtailment)
    ),
    by = .(coord_id, time)
  ]

  # Extract unique locations
  coords_dt <- unique(
    wf_dt[, .(coord_id, lon, lat)]
  )
  coords_sf <- st_as_sf(
    coords_dt,
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  )
  coords_sf <- st_transform(coords_sf, 27700)
  xy <- st_coordinates(coords_sf)
  coords_dt[, `:=`(
    x = xy[, 1] / 1000,
    y = xy[, 2] / 1000
  )]
  wf_dt[
    coords_dt,
    on = "coord_id",
    `:=`(x = i.x, y = i.y)
  ]

  t0 <- min(wf_dt$time, na.rm = TRUE)
  wf_dt[, `:=`(
    t = as.numeric(difftime(time, t0, units = "hours")),
    norm_potential = pmin(1, potential / capacity),
    norm_power_est0 = power_est0 / capacity
  )]
  wf_dt[, error0 := norm_potential - norm_power_est0]
  wf_dt[, `:=`(
    ws_group = inla.group(
      ws_h,
      n = 20,
      method = "quantile"
    ),
    pow_group = inla.group(
      norm_power_est0,
      n = 20,
      method = "quantile"
    ),
    d_coast_group = inla.group(
      dist_coast,
      n = 10,
      method = "quantile"
    ),
    elev_group = inla.group(
      elevation,
      n = 10,
      method = "quantile"
    )
  )]
  wf_dt[, `:=`(
    time_id = as.integer(
      difftime(time, min(time), units = "hours")
    ) +
      1L,
    date = as.Date(time)
  )]
  gb_day_dt <- unique(
    as.data.table(gb_day_df)[,
      .(date, p_group3)
    ]
  )
  wf_dt[
    gb_day_dt,
    on = "date",
    p_group3 := i.p_group3
  ]

  wf_dt[,
    `:=`(
      err_l1 = data.table::shift(error0, 1L, fill = 0),
      err_l2 = data.table::shift(error0, 2L, fill = 0),

      dest = norm_power_est0 -
        data.table::shift(norm_power_est0, 1L, fill = 0)
    ),
    by = coord_id
  ]

  wf_dt[,
    `:=`(
      dest_l1 = data.table::shift(dest, 1L, fill = 0),
      dest_l2 = data.table::shift(dest, 2L, fill = 0)
    ),
    by = coord_id
  ]

  cat("Identifying anomalies in the dataset\n")
  cutprobs3 <- c(0.25, 0.75)
  p_quant3 <- quantile(
    gb_day_df$norm_potential,
    probs = cutprobs3,
    na.rm = TRUE
  )

  wf_dt[,
    anomaly := (norm_potential <= tol &
      norm_power_est0 >= p_quant3[1]) |
      (norm_power_est0 >= 1 - tol &
        norm_potential <= p_quant3[2]) |
      (abs(norm_power_est0 - norm_potential) >= norm_dist_tol)
  ]
  anomaly_perc <- 100 *
    wf_dt[,
      mean(anomaly, na.rm = TRUE)
    ]
  cat(sprintf(
    "Percentage of anomalies at %.1f%% tol in the dataset: %.2f%%\n",
    norm_dist_tol * 100,
    anomaly_perc
  ))

  set.seed(123)
  plot_dt <- wf_dt[
    sample.int(.N, min(.N, 100000L)),
    .(norm_potential, norm_power_est0, anomaly)
  ]

  p <- ggplot(
    plot_dt,
    aes(norm_potential, norm_power_est0, color = anomaly)
  ) +
    geom_point(size = 0.5, alpha = 0.4) +
    scale_color_manual(values = c("grey50", "darkred")) +
    theme(legend.position = "bottom") +
    labs(
      x = "Elexon CF",
      y = "Generic PC (ERA5) CF",
      color = "Anomaly"
    ) +
    coord_fixed(ratio = 1)

  ggsave(
    sprintf(
      "%s/%s/fig/fit/anomalies_%s.png",
      output_path,
      batch_name,
      d0_tag
    ),
    plot = p,
    width = 4,
    height = 4.5,
    dpi = 100
  )
  cat("Masking anomalies for model fitting\n")

  # mask anomalies for model fitting

  # Preserve original values
  wf_dt[, norm_potential_orig := norm_potential]

  # Mask only anomalous observations
  wf_dt[anomaly == TRUE, norm_potential := NA_real_]

  # lwe events
  lwe_output_tabpath <- file.path("tables", "lwe_eda")
  lwe_fname <- file.path(lwe_output_tabpath, "low_wind_event_stats.csv")
  event_dt <- fread(lwe_fname)

  time_dt <- unique(wf_dt[, .(time)])

  time_dt[
    event_dt,
    on = .(time >= start_time, time <= end_time),
    event_id := i.event_id
  ]
  # event_ids <- unique(time_dt$event_id)

  # missing_events <- event_dt[!event_id %in% event_ids]
  # event_id_missing <- which(event_dt$event_id %in% event_ids)

  wf_dt[
    time_dt,
    on = "time",
    event_id := i.event_id
  ]

  if (save_daily_files) {
    wf_df_fname <- sprintf(
      "%s/%s/data/fit/calibration_df_%s_%s.%s",
      output_path,
      batch_name,
      "base",
      d0_tag,
      extension
    )
    if (extension == "geojson" & file.exists(wf_df_fname)) {
      file.remove(wf_df_fname)
    }
    if (extension != "rds") {
      st_write(
        wf_dt,
        wf_df_fname,
        driver = driver,
        append = FALSE,
        quiet = TRUE
      )
    } else {
      saveRDS(
        wf_dt,
        wf_df_fname
      )
    }
  }
}
samp_gb <- wf_dt[
  !is.na(norm_potential),
  .(
    ws_h_wmean = sum(ws_h * capacity) / sum(capacity),

    norm_potential = sum(norm_potential * capacity, na.rm = TRUE) /
      sum(capacity, na.rm = TRUE),

    norm_potential_orig = sum(norm_potential_orig * capacity, na.rm = TRUE) /
      sum(capacity, na.rm = TRUE),

    norm_power_est0 = sum(norm_power_est0 * capacity, na.rm = TRUE) /
      sum(capacity, na.rm = TRUE),

    potential = sum(potential),
    power_est0 = sum(power_est0),
    capacity = sum(capacity)
  ),
  by = .(tech_typ, time)
]

setorder(samp_gb, tech_typ, time)

samp_gb[, `:=`(
  date = as.Date(time),
  error0 = norm_potential - norm_power_est0
)]

samp_gb[,
  `:=`(
    err_l1 = data.table::shift(error0, 1L, fill = 0),
    err_l2 = data.table::shift(error0, 2L, fill = 0),
    dest = norm_power_est0 -
      data.table::shift(norm_power_est0, 1L, fill = 0)
  ),
  by = tech_typ
]

samp_gb[,
  `:=`(
    dest_l1 = data.table::shift(dest, 1L, fill = 0),
    dest_l2 = data.table::shift(dest, 2L, fill = 0)
  ),
  by = tech_typ
]

cat("Number of unique locations:", nrow(wf_dt %>% distinct(x, y)), "\n")
n <- nrow(wf_dt)
cat("Number of records in the dataset:", n, "\n")

dat_endtime <- Sys.time()
dat_time_taken <- difftime(dat_endtime, dat_starttime, units = "auto")
cat(
  "Data preparation took: ",
  round(dat_time_taken, 2),
  " ",
  attr(dat_time_taken, "units"),
  "\n"
)
