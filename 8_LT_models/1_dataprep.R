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

# sample_days_fname <- "data/sample_days_df_25_n150.csv"

# if (!file.exists(sample_days_fname)) {
#   cat("List of days for validation not found, recreating...\n")
#   # generating sample days for 3nd round of modelling
#   set.seed(2)
#   # sample 1 day per regime
#   sampled_days_df <- gb_day_df %>%
#     filter(date >= "2025-01-01") %>%
#     group_by(p_group3) %>%
#     slice_sample(n = 50)

#   write.csv(
#     sampled_days_df,
#     file = "data/sample_days_df_25_n150.csv",
#     row.names = FALSE
#   )

#   cat("List of days for validation saved.\n")
# } else {
#   cat("Loading list of days for validation\n")
#   sampled_days_df <- read.csv(sample_days_fname)
# }

# sampled_days <- sampled_days_df %>%
#   pull(date)

# d0 <- sampled_days[day_id] %>% as.Date()
d0 <- as.Date("2025-01-01")
d0_tag <- base::format(d0, "%y%m%d")

extension <- ifelse(local_run, local_ext, cluster_ext)
df_pattern <- sprintf("^calibration_df_.*_%s\\.%s$", d0_tag, extension)
files_found <- list.files("data", pattern = df_pattern, full.names = TRUE)

if (!override_objects && length(files_found) > 0) {
  cat(
    "Calibration data file already exists for this day. Loading existing data.\n"
  )
  if (extension != "rds") {
    wf_df_frag <- st_read(files_found[1])
  } else {
    wf_df_frag <- readRDS(files_found[1])
  }
} else {
  cat(
    "No existing calibration data file found for this day. Preparing new data.\n"
  )

  pwr_curv_dt <- as.data.table(
    read_parquet(file.path(
      gen_path,
      "power_curve_all_enriched.parquet"
    ))
  )
  setnames(pwr_curv_dt, "halfHourEndTime", "time")
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
    time_id := as.integer(
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

  cat("Filtering data for the selected day and nearby days\n")
  # wf_df_frag <- pwr_curv_df %>%
  #   rename(time = halfHourEndTime) %>%
  #   mutate(
  #     date = as.Date(time),
  #     elevation = pmax(0, elevation),
  #     site_name = site_name %>%
  #       gsub("\\b(wind\\s*farm|wf)\\b", "", ., ignore.case = TRUE) %>%
  #       trimws()
  #   ) %>%
  #   # filter(date %in% sampled_days) %>%
  #   filter(
  #     date >= d0 - max(n.days.before.fit, n.days.before.light),
  #     date <= d0 + n.days - 1
  #   ) %>%
  #   filter(coord_id %in% coord_list$coord_id[coord_list$sampled]) %>%
  #   arrange(site_name) %>%
  #   group_by(lon, lat, time) %>%
  #   summarise(
  #     site_name = first(site_name),
  #     coord_id = first(coord_id),
  #     elevation = first(elevation),
  #     dist_coast = first(dist_coast),
  #     tech_typ = first(tech_typ),
  #     across(c(ws_h, wd10, wd100), mean),
  #     ws_h_wmean = sum(ws_h * capacity) / sum(capacity),
  #     across(c(potential, power_est0, capacity, curtailment), sum),
  #     .groups = "drop"
  #   ) %>%
  #   mutate(t = difftime(time, min(time), units = "hours") %>% as.numeric()) %>%
  #   mutate(
  #     norm_potential = pmin(1, potential / capacity),
  #     norm_power_est0 = power_est0 / capacity,
  #     error0 = norm_potential - norm_power_est0
  #   ) %>%
  #   st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  #   mutate(lon = st_coordinates(.)[, 1], lat = st_coordinates(.)[, 2]) %>%
  #   st_transform(crs = 27700) %>%
  #   mutate(
  #     x = st_coordinates(.)[, 1] / 1000,
  #     y = st_coordinates(.)[, 2] / 1000,
  #   ) %>%
  #   # st_drop_geometry() %>%
  #   mutate(
  #     # site_id = as.integer(factor(site_name)),
  #     ws_group = inla.group(ws_h, n = 20, method = "quantile"),
  #     pow_group = inla.group(norm_power_est0, n = 20, method = "quantile"),
  #     d_coast_group = inla.group(dist_coast, n = 10, method = "quantile"),
  #     elev_group = inla.group(elevation, n = 10, method = "quantile"),
  #     time_id = as.integer(factor(time)),
  #     date = as.Date(time)
  #   ) %>%
  #   left_join(
  #     gb_day_df %>% dplyr::select(date, p_group3),
  #     by = c("date" = "date")
  #   )
  # plot candidate anomalies, norm_potential == 0 and p_group3 == "mid"
  cutprobs3 <- c(0.25, 0.75)
  p_quant3 <- quantile(gb_day_df$norm_potential, probs = cutprobs3)

  cat("Identifying anomalies in the dataset\n")
  wf_df_frag <- wf_df_frag %>%
    mutate(
      anomaly = case_when(
        norm_potential <= tol & norm_power_est0 >= p_quant3[1] ~ TRUE,
        norm_power_est0 >= 1 - tol & norm_potential <= p_quant3[2] ~ TRUE,
        abs(norm_power_est0 - norm_potential) >= norm_dist_tol ~ TRUE,
        TRUE ~ FALSE
      )
    )
  wf_df_frag %>%
    ggplot() +
    geom_point(
      aes(norm_potential, norm_power_est0, color = anomaly),
      size = 0.5
    ) +
    scale_color_manual(values = c("grey50", "darkred")) +
    theme(legend.position = "bottom") +
    labs(
      x = "Elexon CF (%)",
      y = "Generic PC(ERA5) CF (%)",
      color = "Anomaly"
    ) +
    coord_fixed(ratio = 1)

  ggsave(
    sprintf("%s/%s/fig/fit/anomalies_%s.png", output_path, batch_name, d0_tag),
    width = 4,
    height = 4.5,
    dpi = 100
  )
  cat("Masking anomalies for model fitting\n")

  # mask anomalies for model fitting

  wf_df_frag <- wf_df_frag %>%
    mutate(
      norm_potential_orig = norm_potential,
      norm_potential = ifelse(anomaly, NA, norm_potential)
    )

  anomaly_perc <- mean(wf_df_frag$anomaly, na.rm = TRUE) * 100
  cat(sprintf(
    "Percentage of anomalies at %.1f%% tol in the dataset: %.2f%%\n",
    norm_dist_tol * 100,
    anomaly_perc
  ))

  # wf_df_frag <- wf_df_frag %>%
  #   filter(!anomaly)

  # x <- wf_df_frag$pow_group %>% unique() %>% sort()
  # min_jump <- min(diff(sort(x))) / diff(range(x))
  # if (min_jump <= 1e-4) {
  #   wf_df_frag <- wf_df_frag %>%
  #     mutate(pow_group = inla.group(norm_power_est0, n = 20, method = "cut"))
  # }

  # ws_groups <- make_group_breaks(wf_df_frag$ws_h, n = 20)

  # pow_groups <- make_group_breaks(
  #   wf_df_frag$norm_power_est0,
  #   n = 20
  # )

  # d_coast_groups <- make_group_breaks(
  #   wf_df_frag$dist_coast,
  #   n = 10
  # )

  # elev_groups <- make_group_breaks(
  #   wf_df_frag$elevation,
  #   n = 10
  # )

  # wf_df_frag <- wf_df_frag %>%
  #   make_groups(
  #     ws_groups = ws_groups,
  #     pow_groups = pow_groups,
  #     d_coast_groups = d_coast_groups,
  #     elev_groups = elev_groups
  #   ) %>%
  #   group_by(coord_id) %>%
  #   mutate(
  #     err_l1 = lag(error0, default = 0),
  #     err_l2 = lag(error0, n = 2, default = 0),
  #     dest = (norm_power_est0) - lag(norm_power_est0, default = 0),
  #     dest_l1 = lag(dest, default = 0),
  #     dest_l2 = lag(dest, n = 2, default = 0)
  #   ) %>%
  #   ungroup()

  # cat("Converting coordinates to km\n")
  # wf_df_frag <- wf_df_frag %>%
  #   st_geometry() %>%
  #   (\(g) g / 1000)() %>%
  #   st_set_geometry(wf_df_frag, .)

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
        wf_df_frag,
        wf_df_fname,
        driver = driver,
        append = FALSE,
        quiet = TRUE
      )
    } else {
      saveRDS(
        wf_df_frag,
        wf_df_fname
      )
    }
  }
}
samp_gb <- wf_df_frag %>%
  filter(!is.na(norm_potential)) %>%
  group_by(tech_typ, time) %>%
  summarise(
    ws_h_wmean = sum(ws_h * capacity) / sum(capacity),
    across(
      c(norm_potential, norm_potential_orig, norm_power_est0),
      ~ sum(. * capacity, na.rm = TRUE) / sum(capacity, na.rm = TRUE)
    ),
    across(
      c(potential, power_est0, capacity),
      sum
    ),
    .groups = "drop"
  ) %>%
  mutate(date = as.Date(time), error0 = norm_potential - norm_power_est0) %>%
  mutate(
    err_l1 = lag(error0, default = 0),
    err_l2 = lag(error0, n = 2, default = 0),
    dest = (norm_power_est0) - lag(norm_power_est0, default = 0),
    dest_l1 = lag(dest, default = 0),
    dest_l2 = lag(dest, n = 2, default = 0)
  )
wf_df_frag <- wf_df_frag %>% filter(date >= d0 - n.days.before.fit)
cat("Number of unique locations:", nrow(wf_df_frag %>% distinct(x, y)), "\n")
n <- nrow(wf_df_frag)
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
