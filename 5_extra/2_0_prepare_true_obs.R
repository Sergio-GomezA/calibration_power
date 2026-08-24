cat(
  "-------------------------------------------------------------------------------------------------\n"
)
cat("Preparing data for simulation\n")
cat(
  "-------------------------------------------------------------------------------------------------\n"
)

start_time <- Sys.time()

# true data frame
prefixfull <- "elex1H"
# time_sel <- wf_df_frag$time %>% max() + hours(1)
time_sel <- seq(
  from = wf_df_frag$time %>% max() + hours(1),
  by = "1 hour",
  length.out = 24
)
extension <- ifelse(local_run, local_ext, cluster_ext)
df_pattern <- sprintf("^cal_aggr_test_df_.*_%s\\.%s$", d0_tag, extension)
files_found <- list.files("data", pattern = df_pattern, full.names = TRUE)

if (!override_objects && length(files_found) > 0) {
  cat(
    "Calibration data file already exists for this day. Loading existing data.\n"
  )
  if (extension != "rds") {
    true_df <- st_read(files_found[1])
  } else {
    true_df <- readRDS(files_found[1])
  }
} else {
  cat(
    "No existing calibration data file found for this day. Preparing new data.\n"
  )

  pwr_curv_df <- read_parquet(file.path(
    gen_path,
    "power_curve_all_enriched.parquet"
  ))

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

  n.days <- 0
  # n.days.before <- 7

  true_df <- pwr_curv_df %>%
    rename(time = halfHourEndTime) %>%
    filter(time %in% time_sel) %>%
    mutate(
      date = as.Date(time),
      elevation = pmax(0, elevation),
      site_name = site_name %>%
        gsub("\\b(wind\\s*farm|wf)\\b", "", ., ignore.case = TRUE) %>%
        trimws()
    ) %>%
    # filter(date %in% sampled_days) %>%
    # filter(date >= d0 - n.days.before, date <= d0 + n.days - 1) %>%
    # filter(coord_id %in% coord_list$coord_id[coord_list$sampled]) %>%
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
    mutate(t = difftime(time, min(time), units = "hours") %>% as.numeric()) %>%
    mutate(
      norm_potential = pmin(1, potential / capacity),
      norm_power_est0 = power_est0 / capacity,
      error0 = norm_potential - norm_power_est0
    ) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    mutate(lon = st_coordinates(.)[, 1], lat = st_coordinates(.)[, 2]) %>%
    st_transform(crs = 27700) %>%
    mutate(
      x = st_coordinates(.)[, 1] / 1000,
      y = st_coordinates(.)[, 2] / 1000,
    ) %>%
    # st_drop_geometry() %>%
    mutate(
      # site_id = as.integer(factor(site_name)),
      ws_group = inla.group(ws_h, n = 20, method = "quantile"),
      pow_group = inla.group(norm_power_est0, n = 20, method = "quantile"),
      d_coast_group = inla.group(dist_coast, n = 10, method = "quantile"),
      elev_group = inla.group(elevation, n = 10, method = "quantile"),
      time_id = as.integer(factor(time)),
      date = as.Date(time),
      oos = coord_id %in% coord_list$coord_id[!coord_list$sampled],
      observed = norm_potential
    ) %>%
    left_join(
      gb_day_df %>% dplyr::select(date, p_group3),
      by = c("date" = "date")
    )

  # plot candidate anomalies, norm_potential == 0 and p_group3 == "mid"
  cutprobs3 <- c(0.25, 0.75)
  # p_quant3 <- quantile(gb_day_df$norm_power_est0, probs = cutprobs3)
  p_quant3 <- quantile(gb_day_df$norm_potential, probs = cutprobs3)

  true_df <- true_df %>%
    mutate(
      anomaly = case_when(
        norm_potential <= tol & norm_power_est0 >= p_quant3[1] ~ TRUE,
        norm_power_est0 >= 1 - tol & norm_potential <= p_quant3[2] ~ TRUE,
        abs(norm_power_est0 - norm_potential) >= norm_dist_tol ~ TRUE,
        TRUE ~ FALSE
      )
    )
  true_df %>%
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
    sprintf("fig/%s/anomalies_%s.pdf", batch_name, d0_tag),
    width = 4,
    height = 4.5
  )

  # mask anomalies for model fitting

  true_df <- true_df %>%
    mutate(
      norm_potential_orig = norm_potential,
      norm_potential = ifelse(anomaly, NA, norm_potential)
    )

  x <- true_df$pow_group %>% unique() %>% sort()
  min_jump <- min(diff(sort(x))) / diff(range(x))
  if (min_jump <= 1e-4) {
    true_df <- true_df %>%
      mutate(pow_group = inla.group(norm_power_est0, n = 20, method = "cut"))
  }

  cat("Converting coordinates to km\n")
  true_df <- true_df %>%
    st_geometry() %>%
    (\(g) g / 1000)() %>%
    st_set_geometry(true_df, .)

  true_df_fname <- sprintf(
    "data/cal_aggr_test_df_%s_%s.%s",
    "base",
    d0_tag,
    extension
  )
  if (extension == "geojson" & file.exists(true_df_fname)) {
    file.remove(true_df_fname)
  }
  if (extension != "rds") {
    st_write(
      true_df,
      true_df_fname,
      driver = driver,
      append = FALSE,
      quiet = TRUE
    )
  } else {
    saveRDS(
      true_df,
      true_df_fname
    )
  }
}
cat("Number of unique locations:", nrow(true_df %>% distinct(x, y)), "\n")
n <- nrow(true_df)
cat("Number of records in the dataset:", n, "\n")


end_time <- Sys.time()
cat(
  "Time taken to prepare true observations: ",
  round(difftime(end_time, start_time, units = "mins"), 2),
  " minutes\n"
)
