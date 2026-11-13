starttime <- Sys.time()
# 0. Setup ####
cat(
  "-------------------------------------------------------------------------------------------------\n"
)
cat("Running validation for out-of-sample temporal prediction\n")
cat(
  "-------------------------------------------------------------------------------------------------\n"
)
# local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

# 0.1 global parameter #####
# day_id <- 2
# mesh_edge_par <- 20 # km, target edge length for the spatial mesh. 10 is fine, 20 is coarse but faster
# override_objects <- FALSE
# rerun_samples <- FALSE
# prec_init <- log(200)
# batch_name <- "batch2025"

if (local_run) {
  cat("Running in local mode\n")
} else {
  cat("Running in cluster mode\n")
}

# Get task ID and others from command-line arguments
# args <- commandArgs(trailingOnly = TRUE)

# Override defaults only if arguments are provided
# if (length(args) > 0) {
#   day_id <- as.numeric(args[1])
# }
# if (length(args) > 1) {
#   override_objects <- as.logical(args[2])
# }
# if (length(args) > 2) {
#   rerun_samples <- as.logical(args[3])
# }
# if (length(args) > 3) {
#   batch_name <- args[4]
# }
# 0.2 libraries and paths ####
require(parallel)

# if (local_run) {
#   data_path <- "~/Documents/ERA5_at_wf/"
#   gen_path <- "~/Documents/elexon/"
#   model_path <- "~/Documents/elexon/model_objects"
#   sample_path <- "~/Documents/elexon/samples"
#   n_samp <- 10
#   pixel_dims <- c(150, 150)
# } else {
#   data_path <- "/exports/eddie/scratch/s2441782/calibration/data"
#   gen_path <- "/exports/eddie/scratch/s2441782/calibration/data"
#   model_path <- "/exports/eddie/scratch/s2441782/calibration/model_objects"
#   sample_path <- "/exports/eddie/scratch/s2441782/calibration/samples"
#   temp_lib <- "/exports/eddie3_homes_local/s2441782/lib"
#   pixel_dims <- c(300, 300)
#   n_samp <- 1000
#   .libPaths(temp_lib)
# }

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
# require(ggspatial)

# source("aux_funct.R")

# sampled_days_df <- read.csv("data/sample_days_df.csv") %>%
#   mutate(date = as.Date(date))

# sampled_days <- sampled_days_df %>%
#   pull(date)

# sampled_days <- c("2020-08-14", "2024-04-17", "2024-04-12")
# d0 <- sampled_days[day_id] %>% as.Date()
# d0_tag <- base::format(d0, "%y%m%d")

alphas <- c(0.01, seq(0.05, 0.95, by = 0.05), 0.99)
task_prefix0 <- "time"
cat(
  "-------------------------------------------------------------------------------------------------\n"
)
cat(
  " Running validation for day:",
  format(d0, "%Y-%m-%d"),
  " (tag:",
  d0_tag,
  ")\n"
)

cat("override_objects:", override_objects, "\n")
cat("rerun_samples:", rerun_samples, "\n")

cat(
  "-------------------------------------------------------------------------------------------------\n"
)
n.days <- n.days.time
n.days.before <- n.days.before.time
n.hours <- n.hours.time
t0 <- d0 - n.days.before
t1 <- d0 + n.days
th <- t1 + hours(n.hours)
t_ind_adj <- 24 * (n.days.before.fit - n.days.before.time)

cat("Traning sample length for assessment and figures:")
cat(
  "",
  n.days.before,
  ifelse(n.days.before == 1, "day before", "days before"),
  format(d0, "%Y-%m-%d"),
  "\n"
)
cat("Prediction next", n.hours, "hours after", format(t1, "%Y-%m-%d"), "\n")

cat(
  "-------------------------------------------------------------------------------------------------\n"
)

# Read models ####
model_catalog <- read.csv("data/model_catalog.csv") %>%
  na.omit()

if (local_run) {
  model_catalog <- model_catalog #%>%
  # filter(!grepl("fine", mod_labels)) %>%
  # filter(!grepl("st0", est_cols))
}
mod_labels <- model_catalog$mod_labels
est_cols <- model_catalog$est_cols
n_models <- length(est_cols)
names(mod_labels) <- est_cols
mod_vec <- list.files(model_path, pattern = d0_tag, full.names = TRUE)
exclusions <- if (local_run) {
  c("spatial", "fine")
} else {
  c("spatial")
}
mod_vec <- mod_vec[!grepl(paste(exclusions, collapse = "|"), mod_vec)] %>%
  sort() # exclude meshes from st model

mlist_names <- model_list %>% names()

model_df <- model_catalog %>%
  rename(code = est_cols, label = mod_labels) %>%
  arrange(desc(nchar(mode_code_prefix))) %>%
  mutate(
    fname = map_chr(
      mode_code_prefix,
      ~ {
        # browser()
        matches <- mod_vec[str_detect(mod_vec, fixed(.x))]
        if (length(matches) == 1) matches else matches[1]
      }
    )
  ) %>%
  mutate(
    type = case_when(
      # grepl("agg_", code) ~ "agg",
      grepl("lm_", code) ~ "bru",
      grepl("lm", code) ~ "lm",
      grepl("qm", code) ~ "qm",
      TRUE ~ "bru"
    )
  ) %>%
  mutate(
    # finding which models are in the model list
    inlist = vapply(
      mode_code_prefix,
      \(prefix) any(startsWith(mlist_names, prefix)),
      logical(1)
    )
  )

model_df %>%
  filter(is.na(fname), inlist == FALSE) %>%
  pull(code) -> missing_models

if (length(missing_models) > 0) {
  cat(
    "Warning: The following models are missing from the model path:\n",
    paste(missing_models, collapse = ", "),
    "\n"
  )

  model_df <- model_df %>% filter(!is.na(fname) | inlist == TRUE)
  mod_labels <- model_df$label
  est_cols <- model_df$code
  names(mod_labels) <- est_cols
}

extension <- "rds"
df_pattern <- sprintf("^calibration_preddf_.*_%s\\.%s$", d0_tag, extension)
files_found <- list.files("data", pattern = df_pattern, full.names = TRUE)

coord_list_fname <- "data/coord_list.csv"

cat("Loading existing coordinate list\n")
coord_list <- read.csv(coord_list_fname)
coord_list <- coord_list %>%
  mutate(
    sampled_old = sampled,
    sampled = TRUE
  )

if (!override_objects && length(files_found) > 0) {
  cat(
    "Calibration data file already exists for this day. Loading existing data.\n"
  )
  wf_df_pred <- readRDS(files_found[1])
} else {
  if (override_objects) {
    cat(
      "Override_objects is TRUE. Preparing new calibration data for this day.\n"
    )
  } else {
    cat(
      "No existing calibration data file found for this day. Preparing new data.\n"
    )
  }

  cat("Reading and filtering prediction data\n")

  pwr_curv_dt <- as.data.table(
    read_parquet(file.path(
      gen_path,
      "power_curve_all_enriched.parquet"
    ))
  )

  setnames(pwr_curv_dt, "halfHourEndTime", "time")

  wf_df_pred <- pwr_curv_dt[
    time >= t0 &
      time < th &
      coord_id %in% coord_list$coord_id[coord_list$sampled]
  ]

  # Clean variables
  wf_df_pred[, `:=`(
    date = as.Date(time),
    elevation = pmax(0, elevation),
    site_name = trimws(gsub(
      "\\b(wind\\s*farm|wf)\\b",
      "",
      site_name,
      ignore.case = TRUE
    ))
  )]

  # Preserve original ordering for first() in aggregation
  setorder(wf_df_pred, site_name)

  # ============================================================
  # 2. Aggregate observations by location and time
  # ============================================================

  cat("Aggregating wind farm observations\n")

  wf_df_pred <- wf_df_pred[,
    .(
      site_name = first(site_name),
      coord_id = first(coord_id),
      elevation = first(elevation),
      dist_coast = first(dist_coast),
      tech_typ = first(tech_typ),

      ws_h = mean(ws_h),
      wd10 = mean(wd10),
      wd100 = mean(wd100),

      ws_h_wmean = sum(ws_h * capacity) / sum(capacity),

      potential = sum(potential),
      power_est0 = sum(power_est0),
      capacity = sum(capacity),
      curtailment = sum(curtailment)
    ),
    by = .(lon, lat, time)
  ]

  # ============================================================
  # 3. Time and normalised variables
  # ============================================================

  wf_df_pred[, `:=`(
    t = t_ind_adj +
      as.numeric(difftime(time, t0, units = "hours")),

    norm_potential = pmin(1, potential / capacity),

    norm_power_est0 = power_est0 / capacity
  )]

  wf_df_pred[,
    error0 := norm_potential - norm_power_est0
  ]

  # ============================================================
  # 4. Transform unique coordinates to British National Grid
  # ============================================================

  cat("Converting unique coordinates to km\n")

  coords_dt <- unique(
    wf_df_pred[, .(lon, lat)]
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

  # Join coordinates back to observations
  wf_df_pred[
    coords_dt,
    on = .(lon, lat),
    `:=`(
      x = i.x,
      y = i.y
    )
  ]

  # ============================================================
  # 5. INLA grouping variables and time index
  # ============================================================

  cat("Creating INLA grouping variables\n")

  wf_df_pred[, `:=`(
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
    ),

    time_id = match(time, sort(unique(time))) + t_ind_adj,

    date = as.Date(time)
  )]

  # ============================================================
  # 6. Join GB daily classification
  # ============================================================

  gb_day_dt <- unique(
    as.data.table(gb_day_df)[,
      .(date, p_group3)
    ]
  )

  wf_df_pred[
    gb_day_dt,
    on = "date",
    p_group3 := i.p_group3
  ]

  # ============================================================
  # 7. Identify anomalies
  # ============================================================

  cat("Identifying anomalies in the dataset\n")

  wf_df_pred[,
    anomaly := abs(norm_power_est0 - norm_potential) >= norm_dist_tol
  ]

  # Match original case_when() behaviour for NA
  wf_df_pred[is.na(anomaly), anomaly := FALSE]

  anomaly_perc <- 100 *
    wf_df_pred[,
      mean(anomaly)
    ]

  cat(sprintf(
    "Percentage of anomalies at %.1f%% tol in the dataset: %.2f%%\n",
    norm_dist_tol * 100,
    anomaly_perc
  ))

  # ============================================================
  # 8. Mask anomalies for model fitting
  # ============================================================

  cat("Masking anomalies for model fitting\n")

  wf_df_pred[,
    norm_potential_orig := norm_potential
  ]

  wf_df_pred[
    anomaly == TRUE,
    norm_potential := NA_real_
  ]

  # ============================================================
  # 9. Apply previously calculated grouping breaks
  # ============================================================

  cat("Applying model grouping breaks\n")

  # Keep your existing make_groups() function.
  # This assumes it accepts and returns a data.frame/data.table.

  # wf_df_pred <- as.data.table(
  #   make_groups(
  #     wf_df_pred,
  #     ws_groups = ws_groups,
  #     pow_groups = pow_groups,
  #     d_coast_groups = d_coast_groups,
  #     elev_groups = elev_groups
  #   )
  # )

  # ============================================================
  # 10. Temporal lags by location
  # ============================================================

  cat("Calculating temporal lag variables\n")

  setorder(wf_df_pred, coord_id, time)

  wf_df_pred[,
    `:=`(
      err_l1 = shift(error0, 1L, fill = 0),
      err_l2 = shift(error0, 2L, fill = 0),

      dest = norm_power_est0 -
        shift(norm_power_est0, 1L, fill = 0)
    ),
    by = coord_id
  ]

  wf_df_pred[,
    `:=`(
      dest_l1 = shift(dest, 1L, fill = 0),
      dest_l2 = shift(dest, 2L, fill = 0)
    ),
    by = coord_id
  ]

  cat("Prediction dataset ready\n")

  wf_df_fname <- sprintf(
    "%s/%s/data/oos/calibration_preddf_%s_%s.%s",
    output_path,
    batch_name,
    task_prefix0,
    d0_tag,
    extension
  )
  saveRDS(
    wf_df_pred,
    wf_df_fname
  )
  rm(pwr_curv_df)
  gc()
}
n_loc <- nrow(wf_df_pred %>% distinct(x, y))
cat("Number of unique locations:", n_loc, "\n")
n <- nrow(wf_df_pred)
cat("Number of records in the dataset:", n, "\n")
cat(
  "-------------------------------------------------------------------------------------------------\n"
)
# wf_df_pred %>%
#   filter(time %in% seq_hours) %>%
#   st_drop_geometry() %>%
#   group_by(time) %>%
#   summarise(
#     n = n(),
#     across(
#       c(norm_potential_orig, norm_potential, norm_power_est0),
#       ~ sum(. * capacity, na.rm = TRUE) / sum(capacity, na.rm = TRUE)
#     ),
#     across(
#       c(potential, power_est0, capacity),
#       sum
#     ),
#     .groups = "drop"
#   )

## aggr models #####
agg_df <- model_df %>% filter(grepl("agg", code))

if (save_models) {
  agg_df %>% pull(fname) %>% map(readRDS) -> mod_list
  names(mod_list) <- agg_df %>% pull(code)
} else {
  agg_mods <- paste0(agg_df %>% pull(mode_code_prefix), d0_tag, ".rds")
  mod_list <- model_list[agg_mods]
}

# create aggregated data for prediction
agg_df_pred <- wf_df_pred %>%
  filter(!is.na(norm_potential)) %>%
  st_drop_geometry() %>%
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
  mutate(date = as.Date(time)) %>%
  mutate(
    error0 = norm_potential - norm_power_est0,
    err_l1 = lag(error0, default = 0),
    err_l2 = lag(error0, n = 2, default = 0),
    dest = (norm_power_est0) - lag(norm_power_est0, default = 0),
    dest_l1 = lag(dest, default = 0),
    dest_l2 = lag(dest, n = 2, default = 0)
  )

agg_pred <- lapply(
  seq(mod_list),
  function(mod) {
    predict(
      mod_list[[mod]],
      newdata = agg_df_pred,
      interval = "prediction"
    ) %>%
      as.data.frame() %>%
      rename(
        estimate = fit,
        lwr = lwr,
        upr = upr
      ) %>%
      bind_cols(
        agg_df_pred %>%
          dplyr::select(
            time,
            date,
            norm_potential,
            norm_potential_orig,
            norm_power_est0,
            ws_h_wmean,
            capacity,
            tech_typ
          ),
        .
      ) %>%
      mutate(
        estimate = pmin(1, pmax(0, estimate)),
        lwr = pmin(1, pmax(0, lwr)),
        upr = pmin(1, pmax(0, upr)),
        std_error = summary(mod_list[[mod]])$sigma,
        model = agg_df$code[mod]
      )
  }
) %>%
  bind_rows()

agg_pred_fig_df <- agg_pred %>%
  group_by(time, model) %>%
  summarise(
    mean = sum(estimate * capacity) / sum(capacity),
    std_error = if (first(model) == "lm") {
      sqrt(sum((capacity / sum(capacity))^2 * std_error^2))
    } else {
      mean(std_error)
    },
    norm_potential = sum(norm_potential * capacity, na.rm = TRUE) /
      sum(capacity),
    norm_potential_orig = sum(norm_potential_orig * capacity, na.rm = TRUE) /
      sum(capacity),
    norm_power_est0 = sum(norm_power_est0 * capacity) / sum(capacity),
    .groups = "drop"
  ) %>%
  mutate(
    lwr = mean - 1.96 * std_error,
    upr = mean + 1.96 * std_error,
    lwr = pmin(1, pmax(0, lwr)),
    upr = pmin(1, pmax(0, upr))
  )

if (local_run) {
  agg_pred_fig_df %>%
    filter(time >= t1, time < t1 + hours(24)) %>%
    ggplot(aes(x = norm_potential, y = mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.01, alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "red")

  agg_pred_fig_df %>%
    # filter(time >= t1, time < t1 + hours(72)) %>%
    filter(time >= t0 - hours(24), time <= t1 + hours(72)) %>%
    ggplot(aes(x = time, y = mean, color = model)) +
    geom_line() +
    geom_line(aes(time, y = norm_power_est0, col = "PC")) +
    geom_ribbon(
      aes(ymin = lwr, ymax = upr, fill = model),
      alpha = 0.2,
      color = NA
    ) +
    geom_point(aes(y = norm_potential), color = "black") +
    theme_minimal() +
    scale_x_datetime(date_labels = "%H:%M")
}


## quantile mapping ####
qm_df <- model_df %>% filter(type == "qm")

if (save_models) {
  qm_df %>% pull(fname) %>% map(readRDS) -> mod_list
  names(mod_list) <- qm_df %>% pull(code)
} else {
  qm_mods <- paste0(qm_df %>% pull(mode_code_prefix), d0_tag, ".rds")
  mod_list <- model_list[qm_mods]
}

qm_pred_df <- lapply(
  seq(mod_list),
  function(mod) {
    wf_df_pred %>%
      dplyr::select(
        coord_id,
        site_name,
        time,
        date,
        norm_potential,
        norm_potential_orig,
        norm_power_est0,
        capacity,
        tech_typ,
        p_group3,
        anomaly
      ) %>%
      mutate(
        estimate = doQmapQUANT(
          norm_power_est0,
          mod_list[[mod]],
          type = "linear"
        )
      ) %>%
      mutate(
        estimate = pmin(1, pmax(0, estimate)),
        model = qm_df$code[mod]
      )
  }
) %>%
  bind_rows()

qm_pred_fig_df <- qm_pred_df %>%
  filter(!is.na(norm_potential)) %>%
  st_drop_geometry() %>%
  group_by(time, model) %>%
  summarise(
    mean = sum(estimate * capacity) / sum(capacity),
    norm_potential = sum(norm_potential * capacity, na.rm = TRUE) /
      sum(capacity),
    norm_potential_orig = sum(norm_potential_orig * capacity, na.rm = TRUE) /
      sum(capacity),
    norm_power_est0 = sum(norm_power_est0 * capacity) / sum(capacity),
    .groups = "drop"
  )


## bru models ####

bru_df <- model_df %>% filter(type == "bru")

# source("aux_funct.R")
# mod_temp <- bruar1

### loop throgh all models #####

# pred band summary
pred_summary_fname <- sprintf(
  "%s/%s/summaries/oos/pred_band_summary_%s_%s.rds",
  output_path,
  batch_name,
  task_prefix0,
  d0_tag
)
pred_samples_fname <- file.path(
  sample_path,
  sprintf(
    "pred_samples_%s.rds",
    d0_tag
  )
)
cov_summary_fname <- sprintf(
  "%s/%s/summaries/oos/pred_band_coverage_summary_%s_%s.rds",
  output_path,
  batch_name,
  task_prefix0,
  d0_tag
)
scores_summary_fname <- sprintf(
  "%s/%s/summaries/oos/model_scores_summary_%s_%s.csv",
  output_path,
  batch_name,
  task_prefix0,
  d0_tag
)
# source("aux_funct.R")
if (!file.exists(pred_summary_fname) || rerun_samples) {
  if (!file.exists(pred_summary_fname)) {
    cat("Prediction band summary file not found, creating new summary\n")
  } else {
    cat(
      "Prediction band summary file found, but rerun_samples is TRUE. Recreating summary\n"
    )
  }

  bru_mods <- paste0(bru_df %>% pull(mode_code_prefix), d0_tag, ".rds")
  pred_band_summary <- lapply(
    seq_along(bru_df$fname),
    function(i) {
      initime <- Sys.time()
      # debug(bru_ci_plot)
      cat(
        "-------------------------------------------------------------------------------------------------\n"
      )
      cat("Processing model:", bru_df$label[i], "\n")
      cat(
        "-------------------------------------------------------------------------------------------------\n"
      )
      test <- tryCatch(
        {
          # browser()
          bru_ci_plot(
            bru_model = model_list[[bru_mods[i]]],
            newdata = wf_df_pred,
            n.samples = n_samp,
            show.fig = FALSE,
            alphas = alphas,
            oos_type = "time",
            family = bru_df$family[i]
          )
        },
        error = function(e) {
          cat("Error in processing model:", bru_df$label[i], "\n")
          cat("Error message:", e$message, "\n")
          return(NULL)
        }
      )
      finaltime <- Sys.time()
      cat(
        "Time taken for sampling model",
        bru_df$label[i],
        ":",
        round(difftime(finaltime, initime, units = "mins"), 2),
        "minutes\n"
      )
      test
    }
  )
  names(pred_band_summary) <- bru_df$code
  summary_only <- lapply(
    pred_band_summary,
    function(x) {
      list(
        GB_summary = x$GB_summary,
        wf_summary = x$wf_summary,
        formula = x$formula
      )
    }
  )
  coverage_summary <- lapply(
    pred_band_summary,
    function(x) {
      list(
        cov_gbl = x$cov_gbl,
        cov_time = x$cov_time,
        cov_loc = x$cov_loc,
        cov_24h = x$cov_24h,
        cov_12h = x$cov_12h
      )
    }
  )
  scores_summary <- pred_band_summary %>%
    lapply(\(x) {
      tibble(
        crps = x$scores$crps,
        energy = x$scores$energy,
        log = x$scores$log,
        !!!setNames(as.list(x$scores$brier), names(x$scores$brier)),
        crps_24 = x$scores_24$crps,
        energy_24 = x$scores_24$energy,
        log_24 = x$scores_24$log,
        !!!setNames(
          as.list(x$scores_24$brier),
          paste0(names(x$scores_24$brier), "_24")
        ),
        crps_12 = x$scores_12$crps,
        energy_12 = x$scores_12$energy,
        log_12 = x$scores_12$log,
        !!!setNames(
          as.list(x$scores_12$brier),
          paste0(names(x$scores_12$brier), "_12")
        ),
        gbl_width = x$band_width$gbl,
        wf_width = x$band_width$wf
      )
    }) %>%
    bind_rows(.id = "model")

  names(summary_only) <- bru_df$code
  names(coverage_summary) <- bru_df$code

  saveRDS(summary_only, pred_summary_fname)
  saveRDS(coverage_summary, cov_summary_fname)
  write.csv(scores_summary, scores_summary_fname, row.names = FALSE)

  # if (save_samples) {
  #   samples_only <- lapply(
  #     pred_band_summary,
  #     function(x) {
  #       list(
  #         sample_df = x$sample_df
  #       )
  #     }
  #   )
  #   names(samples_only) <- bru_df$code
  #   saveRDS(samples_only, pred_samples_fname)
  # }
  # rm(pred_band_summary)
  # gc()
} else {
  cat("Loading existing prediction band summary\n")
  summary_only <- readRDS(pred_summary_fname)
}


gb_fig_df <- bind_rows(
  # lm_pred_fig_df %>% filter(model == "lm"),
  agg_pred_fig_df,
  qm_pred_fig_df,
  lapply(
    bru_df$code,
    function(code) {
      summary_only[[code]]$GB_summary %>%
        as.data.frame() %>%
        mutate(
          model = code
        )
    }
  ) %>%
    bind_rows() %>%
    left_join(
      agg_pred_fig_df %>% dplyr::select(time, norm_power_est0) %>% unique(),
      by = c("time")
    )
) %>%
  mutate(
    oos = ifelse(time >= t1, TRUE, FALSE)
  )
saveRDS(
  gb_fig_df,
  sprintf(
    "%s/%s/summaries/oos/GB_fig_band_summary_%s_%s.rds",
    output_path,
    batch_name,
    task_prefix0,
    d0_tag
  )
)
# gb_fig_df$time %>% range()
excluded_models0 <- c("lm")
gb_fig_df %>%
  filter(!model %in% excluded_models0) %>%
  filter(time >= t1 - hours(3), time <= t1 + hours(24)) %>%
  ggplot() +
  geom_ribbon(
    aes(
      x = time,
      ymin = lwr,
      ymax = upr,
      fill = "95% CI"
    ),
    # fill = blues9[5],
    alpha = 0.5
  ) +
  geom_line(
    aes(x = time, y = norm_power_est0, col = "PC(ERA5)"),
    # color = blues9[9],
    lwd = 1
  ) +
  geom_line(
    aes(x = time, y = mean, col = "fit"),
    # color = blues9[9],
    lwd = 1
  ) +
  geom_line(
    aes(x = time, y = norm_potential_orig, col = "observed"),
    # color = "darkred",
    lwd = 1
  ) +
  geom_vline(
    xintercept = t1,
    linetype = "dashed",
    color = "black"
  ) +
  # coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(~model, nrow = 2, labeller = as_labeller(mod_labels)) +
  # scale_x_datetime(date_labels = "%m/%d") +
  scale_x_datetime(date_labels = "%H:%M") +
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("fit" = blues9[9], "observed" = "darkred", "PC(ERA5)" = "gray70")
  ) +
  scale_fill_manual(values = c("95% CI" = blues9[5])) +
  labs(fill = "", color = "", y = "Normalised power output")

ggsave(
  filename = sprintf(
    "%s/%s/fig/oos/GB_pred_band_%s_%s.pdf",
    output_path,
    batch_name,
    task_prefix0,
    d0_tag
  ),
  width = 10,
  height = 6,
  # dpi = 300
)


## lwe events ####
lwe_output_tabpath <- file.path("tables", "lwe_eda")
lwe_fname <- file.path(lwe_output_tabpath, "low_wind_event_stats.csv")
event_dt <- fread(lwe_fname)

time_dt <- unique(wf_df_pred[, .(time)])

time_dt[
  event_dt,
  on = .(time >= start_time, time <= end_time),
  event_id := i.event_id
]


wf_df_pred[
  time_dt,
  on = "time",
  event_id := i.event_id
]

wf_df_pred[, lwe := !is.na(event_id)]

t1_lwe <- "2025-07-28 22:00:00" %>% as.POSIXct(tz = "UTC")
t1_lwe_tag <- format(t1_lwe, "%y%m%d%H")
task_prefix01 <- "tlwe"
gb_fig_df %>%
  filter(!model %in% excluded_models0) %>%
  filter(time >= t1_lwe - hours(3), time <= t1_lwe + hours(24)) %>%
  ggplot() +
  geom_ribbon(
    aes(
      x = time,
      ymin = lwr,
      ymax = upr,
      fill = "95% CI"
    ),
    # fill = blues9[5],
    alpha = 0.5
  ) +
  geom_line(
    aes(x = time, y = norm_power_est0, col = "PC(ERA5)"),
    # color = blues9[9],
    lwd = 1
  ) +
  geom_line(
    aes(x = time, y = mean, col = "fit"),
    # color = blues9[9],
    lwd = 1
  ) +
  geom_line(
    aes(x = time, y = norm_potential_orig, col = "observed"),
    # color = "darkred",
    lwd = 1
  ) +
  geom_vline(
    xintercept = t1_lwe,
    linetype = "dashed",
    color = "black"
  ) +
  # coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(~model, nrow = 2, labeller = as_labeller(mod_labels)) +
  # scale_x_datetime(date_labels = "%m/%d") +
  scale_x_datetime(date_labels = "%H:%M") +
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("fit" = blues9[9], "observed" = "darkred", "PC(ERA5)" = "gray70")
  ) +
  scale_fill_manual(values = c("95% CI" = blues9[5])) +
  labs(fill = "", color = "", y = "Normalised power output")

ggsave(
  filename = sprintf(
    "%s/%s/fig/oos/GB_pred_band_%s_%s.pdf",
    output_path,
    batch_name,
    task_prefix01,
    t1_lwe_tag
  ),
  width = 10,
  height = 6,
  # dpi = 300
)

t2_lwe <- "2025-08-23 01:00:00" %>% as.POSIXct(tz = "UTC")
t2_lwe_tag <- format(t2_lwe, "%y%m%d%H")
gb_fig_df %>%
  filter(!model %in% excluded_models0) %>%
  filter(time >= t2_lwe - hours(3), time <= t2_lwe + hours(45)) %>%
  ggplot() +
  geom_ribbon(
    aes(
      x = time,
      ymin = lwr,
      ymax = upr,
      fill = "95% CI"
    ),
    # fill = blues9[5],
    alpha = 0.5
  ) +
  geom_line(
    aes(x = time, y = norm_power_est0, col = "PC(ERA5)"),
    # color = blues9[9],
    lwd = 1
  ) +
  geom_line(
    aes(x = time, y = mean, col = "fit"),
    # color = blues9[9],
    lwd = 1
  ) +
  geom_line(
    aes(x = time, y = norm_potential_orig, col = "observed"),
    # color = "darkred",
    lwd = 1
  ) +
  geom_vline(
    xintercept = t2_lwe,
    linetype = "dashed",
    color = "black"
  ) +
  # coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(~model, nrow = 2, labeller = as_labeller(mod_labels)) +
  # scale_x_datetime(date_labels = "%m/%d") +
  scale_x_datetime(date_labels = "%H:%M") +
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("fit" = blues9[9], "observed" = "darkred", "PC(ERA5)" = "gray70")
  ) +
  scale_fill_manual(values = c("95% CI" = blues9[5])) +
  labs(fill = "", color = "", y = "Normalised power output")
ggsave(
  filename = sprintf(
    "%s/%s/fig/oos/GB_pred_band_%s_%s.pdf",
    output_path,
    batch_name,
    task_prefix01,
    t2_lwe_tag
  ),
  width = 10,
  height = 6,
  # dpi = 300
)
