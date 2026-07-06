starttime <- Sys.time()
# low wind events distribution

# 0. Setup ####

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

# 0.1 global parameter #####
day_id <- 2
mesh_edge_par <- 20 # km, target edge length for the spatial mesh. 10 is fine, 20 is coarse but faster
override_objects <- FALSE
batch_name <- "batch2025"
# prec_init <- log(200)

if (local_run) {
  cat("Running in local mode\n")
} else {
  cat("Running in cluster mode\n")
}

# Get task ID and others from command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Override defaults only if arguments are provided
if (length(args) > 0) {
  day_id <- as.numeric(args[1])
}
if (length(args) > 1) {
  override_objects <- as.logical(args[2])
}
if (length(args) > 2) {
  batch_name <- args[3]
}

# 0.2 libraries and paths ####
require(parallel)

if (local_run) {
  data_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
  sample_path <- "~/Documents/elexon/samples"
  pixel_dims <- c(150, 150)
} else {
  data_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  gen_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  model_path <- "/exports/eddie/scratch/s2441782/calibration/model_objects"
  sample_path <- "/exports/eddie/scratch/s2441782/calibration/samples"
  temp_lib <- "/exports/eddie3_homes_local/s2441782/lib"
  pixel_dims <- c(300, 300)
  .libPaths(temp_lib)
}

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
require(data.table)
# require(ggspatial)

source("aux_funct.R")

mc <- ifelse(local_run, 1, available_cores())
# inla.setOption(num.threads = sprintf("%d:1", mc))

pow_threshold <- 0.1
pow_threshold_label <- as.character(pow_threshold) %>% gsub("\\.", "_", .)
# 1. data preparation ####
cat("--------------------------------------------------------------------\n")
cat("Preparing data for low wind events analysis\n")
cat("--------------------------------------------------------------------\n")


sampled_days_df <- read.csv("data/sample_days_df.csv") %>%
  mutate(date = as.Date(date))

sampled_days <- sampled_days_df %>%
  pull(date)
# sampled_days <- c("2020-08-14", "2024-04-17", "2024-04-12")

d0 <- sampled_days[day_id] %>% as.Date()
d0_tag <- base::format(d0, "%y%m%d")
extension <- ifelse(!local_run, "rds", "rds")


coord_list_fname <- "data/coord_list.csv"

cat("Loading existing coordinate list\n")
coord_list <- read.csv(coord_list_fname)

cat(
  "--------------------------------------------------------------------\n"
)
cat(
  " Low wind events analysis for days:",
  paste(format(sampled_days, "%Y-%m-%d"), collapse = ", "),
  ")\n"
)

cat("override_objects:", override_objects, "\n")

cat(
  "--------------------------------------------------------------------\n"
)

model_df0 <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days[i] %>% as.Date()
    # print(d0)
    d0_tag <- base::format(d0, "%y%m%d")

    readRDS(sprintf(
      "data/calibration_df_%s_%s.%s",
      # mesh_label,
      "coarse",
      d0_tag,
      extension
    ))
  }
) %>%
  bind_rows()

pos_breaks <- with(
  model_df0,
  quantile(elevation[elevation > 0], probs = seq(0, 1, 1 / 3))
)
pos_levels <- levels(cut(
  model_df0$elevation[model_df0$elevation > 0],
  breaks = pos_breaks,
  include.lowest = TRUE
))
mod_labels <- c(
  "Generic PC",
  "Linear model",
  "AR1 model",
  "AR2 model",
  "1D SPDE model",
  "Spatio-temporal model",
  "QM",
  "GB LM"
)
est_cols <- c(
  "norm_power_est0",
  "lm",
  "ar1",
  "ar2",
  "spde1d",
  "st",
  "qm",
  "agg_lm"
)

if (local_run) {
  mod_labels <- mod_labels[!grepl("fine", mod_labels)]
  est_cols <- est_cols[!grepl("st0_m2", est_cols)]
}
n_models <- length(est_cols)
names(mod_labels) <- est_cols
samp_vec <- list.files(sample_path, pattern = d0_tag, full.names = TRUE)
# exclusions <- if (local_run) {
#   c("spatial", "fine")
# } else {
#   c("spatial")
# }
# mod_vec <- mod_vec[!grepl(paste(exclusions, collapse = "|"), mod_vec)] %>%
#   sort() # exclude meshes from st model

df_long0 <- model_df0 %>%
  dplyr::select(
    date,
    time,
    site_name,
    coord_id,
    elevation,
    dist_coast,
    capacity,
    tech_typ,
    p_group3,
    norm_potential,
    any_of(est_cols)
  ) %>%
  mutate(
    hour = hour(time),
    elevation = pmax(0, elevation),
    p_group3 = factor(p_group3, levels = c("low", "mid", "high")),
    dist_coast_g4 = cut(
      dist_coast,
      breaks = quantile(dist_coast, probs = seq(0, 1, 0.25)),
      include.lowest = TRUE
    ),
    elevation_g4 = ifelse(
      elevation == 0,
      "0",
      as.character(
        cut(
          elevation,
          breaks = pos_breaks,
          include.lowest = TRUE,
          # labels = c("Low", "Mid", "High")
        )
      )
    ),
    elevation_g4 = factor(
      elevation_g4,
      levels = c("0", pos_levels)
    )
  ) %>%
  pivot_longer(
    cols = any_of(est_cols),
    names_to = "model",
    values_to = "estimate"
  ) %>%
  mutate(
    estimate = pmin(1, pmax(0, estimate)), # clipping estimates to [0, 1]
    err = estimate - norm_potential,
    p_group3 = forcats::fct_rev(p_group3),
    model = factor(model, levels = est_cols, labels = mod_labels)
  )

# 2.1 Low wind events in observed data ####
cat("--------------------------------------------------------------------\n")
cat("Low wind events in observed data\n")
cat("--------------------------------------------------------------------\n")


lwe_obs_pred_fname <- file.path(
  sample_path,
  sprintf(
    "lwe_obs_pred_%s_t%s.%s",
    "sampled_days",
    pow_threshold_label,
    extension
  )
)
lwe_obs_newloc_fname <- file.path(
  sample_path,
  sprintf(
    "lwe_obs_newloc_%s_t%s.%s",
    "sampled_days",
    pow_threshold_label,
    extension
  )
)

if (!file.exists(lwe_obs_pred_fname) | override_objects) {
  if (!file.exists(lwe_obs_pred_fname)) {
    cat("LWE processed file not found, creating new summary\n")
  } else {
    cat(
      "LWE processed file found, but override_objects is TRUE. Recreating summary\n"
    )
  }

  pwr_curv_df <- read_parquet(file.path(
    gen_path,
    "power_curve_all_enriched.parquet"
  ))

  cat("This may take a while...\n")

  # d0 <- as.Date("2024-01-01")
  # n.days <- 365

  lwe_obs_pred <- pwr_curv_df %>%
    rename(time = halfHourEndTime) %>%
    mutate(
      date = as.Date(time)
    ) %>%
    # filter(date >= d0, date <= d0 + n.days) %>%
    filter(date %in% sampled_days) %>%
    filter(coord_id %in% coord_list$coord_id[coord_list$sampled]) %>%
    mutate(
      elevation = pmax(0, elevation),
      site_name = site_name %>%
        gsub("\\b(wind\\s*farm|wf)\\b", "", ., ignore.case = TRUE) %>%
        trimws()
    ) %>%
    # filter(date %in% sampled_days) %>%
    # filter(date >= d0, date <= d0 + n.days - 1) %>%
    arrange(site_name) %>%
    mutate(
      site_id = as.integer(factor(site_name)),
      coord_id = as.integer(factor(paste(lon, lat)))
    ) %>%
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
    arrange(coord_id, time) %>%
    group_by(coord_id, site_name) %>%
    mutate(
      below = norm_potential < pow_threshold,
      run_id = cumsum(below != lag(below, default = first(below)))
    ) %>%
    group_by(coord_id, site_name, run_id, below) %>%
    summarise(
      start_time = first(time),
      end_time = last(time),
      duration_h = pmin(
        100,
        as.numeric(difftime(end_time, start_time, units = "hours")) +
          1
      ),
      .groups = "drop"
    ) %>%
    filter(below) %>%
    dplyr::select(coord_id, start_time, duration_h) %>%
    mutate(
      model = "observed"
    )

  # low_events %>% pull(duration_h) %>% summary()

  saveRDS(lwe_obs_pred, lwe_obs_pred_fname)

  n.days <- 0
  n.days.before <- 7
  n.hours <- 1

  days_newloc <- lapply(
    sampled_days %>% as.Date(),
    \(d0) {
      seq(d0 - n.days.before, d0 + n.days - 1, by = "days") %>%
        as.Date()
    }
  ) %>%
    unlist() %>%
    unique() %>%
    sort() %>%
    as.Date()

  lwe_obs_newloc <- pwr_curv_df %>%
    rename(time = halfHourEndTime) %>%
    mutate(
      date = as.Date(time)
    ) %>%
    # filter(date >= d0, date <= d0 + n.days) %>%
    filter(date %in% days_newloc) %>%
    filter(coord_id %in% coord_list$coord_id[!coord_list$sampled]) %>%
    mutate(
      elevation = pmax(0, elevation),
      site_name = site_name %>%
        gsub("\\b(wind\\s*farm|wf)\\b", "", ., ignore.case = TRUE) %>%
        trimws()
    ) %>%
    # filter(date %in% sampled_days) %>%
    # filter(date >= d0, date <= d0 + n.days - 1) %>%
    arrange(site_name) %>%
    mutate(
      site_id = as.integer(factor(site_name)),
      coord_id = as.integer(factor(paste(lon, lat)))
    ) %>%
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
    arrange(coord_id, time) %>%
    group_by(coord_id, site_name) %>%
    mutate(
      below = norm_potential < pow_threshold,
      run_id = cumsum(below != lag(below, default = first(below)))
    ) %>%
    group_by(coord_id, site_name, run_id, below) %>%
    summarise(
      start_time = first(time),
      end_time = last(time),
      duration_h = pmin(
        100,
        as.numeric(difftime(end_time, start_time, units = "hours")) +
          1
      ),
      .groups = "drop"
    ) %>%
    filter(below) %>%
    dplyr::select(coord_id, start_time, duration_h) %>%
    mutate(
      model = "observed"
    )

  saveRDS(lwe_obs_newloc, lwe_obs_newloc_fname)
} else {
  cat("Loading existing low wind events observed data\n")
  lwe_obs_pred <- readRDS(lwe_obs_pred_fname)
  lwe_obs_newloc <- readRDS(lwe_obs_newloc_fname)
}


# lwe_obs_newloc %>%
#   filter(duration_h < 1000) %>%
#   ggplot(aes(x = duration_h)) +
#   geom_density(fill = blues9[5], alpha = 0.5, bw = 0.14) +
#   scale_x_log10(n.breaks = 8) +
#   labs(
#     title = "Distribution of Low Wind Events Duration",
#     x = "Duration (hours)",
#     y = "Frequency"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )

# 3. Low wind events in model predictions ####
cat("--------------------------------------------------------------------\n")
cat("Low wind events in model mean predictions\n")
cat("--------------------------------------------------------------------\n")

low_events_model <- df_long0 %>%
  st_drop_geometry() %>%
  # filter(model == "Spatio-temporal model") %>%
  filter(coord_id %in% coord_list$coord_id[coord_list$sampled]) %>%
  arrange(model, coord_id, time) %>%
  group_by(model, coord_id) %>%
  mutate(
    below = estimate < pow_threshold,
    run_id = cumsum(below != lag(below, default = first(below)))
  ) %>%
  filter(below) %>%
  group_by(model, coord_id, run_id) %>%
  summarise(
    start_time = first(time),
    end_time = last(time),
    duration_h = as.numeric(difftime(end_time, start_time, units = "hours")) +
      1,
    .groups = "drop"
  ) %>%
  dplyr::select(model, coord_id, start_time, duration_h) %>%
  bind_rows(lwe_obs_pred)


low_events_model %>%
  filter(duration_h < 25) %>%
  # filter(
  #   model %in%
  #     c("observed", "Linear model", "AR1 model", "QM", "Spatio-temporal model")
  # ) %>%
  ggplot(aes(x = duration_h)) +
  # geom_density(aes(fill = model), alpha = 0.5) +
  geom_histogram(
    aes(fill = model),
    position = "identity",
    alpha = 0.5,
    binwidth = 1
  ) +
  # scale_x_log10(n.breaks = 8) +
  labs(
    title = "Distribution of Low Wind Events Duration",
    x = "Duration (hours)",
    y = "Frequency"
  ) +
  facet_wrap(~model, scales = "free", labeller = as_labeller(mod_labels)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(
  sprintf(
    "fig/%s/low_wind_duration_dist_t%s.pdf",
    batch_name,
    pow_threshold_label
  ),
  width = 10,
  height = 6
)

obs <- low_events_model %>%
  filter(
    model == "observed",
    duration_h < 25
  ) %>%
  pull(duration_h)

probs <- seq(0, 1, length.out = 100)

qq_df <- low_events_model %>%
  filter(duration_h < 25) %>%
  group_by(model) %>%
  summarise(
    obs_q = list(quantile(obs, probs = probs, na.rm = TRUE)),
    model_q = list(quantile(duration_h, probs = probs, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  unnest(c(obs_q, model_q))

ggplot(qq_df, aes(x = obs_q, y = model_q)) +
  geom_point(size = 1) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  facet_wrap(~model, scales = "free", labeller = as_labeller(mod_labels)) +
  labs(
    title = "Q-Q Plot of Low Wind Event Durations",
    x = "Observed quantiles",
    y = "Model quantiles"
  ) +
  # coord_equal() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  sprintf(
    "fig/%s/low_wind_duration_qq_t%s.pdf",
    batch_name,
    pow_threshold_label
  ),
  width = 10,
  height = 6
)

# 4. Low wind events in model samples ####
## 4.1 forecast samples ####
cat("--------------------------------------------------------------------\n")
cat("Low wind events in model samples forecast\n")
cat("--------------------------------------------------------------------\n")

# read samples

lwe_pred_fname <- file.path(
  sample_path,
  sprintf(
    "lwe_pred_%s_t%s.%s",
    d0_tag,
    pow_threshold_label,
    extension
  )
)

# if (!file.exists(lwe_pred_fname) | override_objects) {
#   if (!file.exists(lwe_pred_fname)) {
#     cat("LWE processed file not found, creating new summary\n")
#   } else {
#     cat(
#       "LWE processed file found, but override_objects is TRUE. Recreating summary\n"
#     )
#   }
#   # install.packages("profvis")
#   # library(profvis)
#   sampled_days2 <- sampled_days[1]
#   cat("This may take a while...\n")
#   samp_start_time <- Sys.time()
#   # profvis(
#   lwe_pred_df <- lapply(
#     seq_along(sampled_days2),
#     function(i) {
#       d0 <- sampled_days2[i] %>% as.Date()
#       d0_tag <- base::format(d0, "%y%m%d")

#       # browser()
#       pred_samples_fname <- file.path(
#         sample_path,
#         sprintf(
#           "pred_samples_%s.%s",
#           d0_tag,
#           extension
#         )
#       )
#       cat(
#         "--------------------------------------------------------------------\n"
#       )
#       cat("Reading prediction samples from:", pred_samples_fname, "\n")
#       full_file <- readRDS(pred_samples_fname)

#       cat(
#         "Processing prediction samples for day:",
#         format(d0, "%Y-%m-%d"),
#         "\n"
#       )
#       mod_codes <- names(full_file)
#       lapply(
#         seq_along(mod_codes),
#         function(j) {
#           # browser()
#           cat(
#             "Processing model:",
#             mod_codes[j],
#             "\n"
#           )
#           full_file[[j]]$sample_df %>%
#             mutate(
#               model = mod_codes[j],
#               date = d0
#             ) %>%
#             arrange(sim, coord_id, time) %>%
#             group_by(model, sim, coord_id) %>%
#             # apply threshold to get low wind events
#             mutate(
#               below = fit < pow_threshold,
#               run_id = cumsum(below != lag(below, default = first(below)))
#             ) %>%
#             filter(below, time >= d0) %>% # only lwe from oos
#             group_by(sim, coord_id, run_id) %>%
#             summarise(
#               model = first(model),
#               start_time = first(time),
#               end_time = last(time),
#               duration_h = as.numeric(difftime(
#                 end_time,
#                 start_time,
#                 units = "hours"
#               )) +
#                 1,
#               .groups = "drop"
#             ) %>%
#             dplyr::select(sim, model, coord_id, start_time, duration_h)
#         }
#       ) %>%
#         bind_rows()
#     }
#   ) %>%
#     bind_rows()
#   # )

#   samp_end_time <- Sys.time()
#   cat(
#     "Time taken to process prediction samples:",
#     difftime(samp_end_time, samp_start_time, units = "secs"),
#     "seconds\n"
#   )
#   cat("Saving low wind events prediction samples to:", lwe_pred_fname, "\n")
#   # saveRDS(lwe_pred_df, lwe_pred_fname)
# } else {
#   cat("Loading existing low wind events prediction samples\n")
#   lwe_pred_df <- readRDS(lwe_pred_fname)
# }
## alternative version using data.tables #####
if (!file.exists(lwe_pred_fname) | override_objects) {
  if (!file.exists(lwe_pred_fname)) {
    cat("LWE processed file not found, creating new summary\n")
  } else {
    cat(
      "LWE processed file found, but override_objects is TRUE. Recreating summary\n"
    )
  }
  # install.packages("profvis")
  # library(profvis)
  # sampled_days2 <- sampled_days[1]
  cat("This may take a while...\n")
  samp_start_time <- Sys.time()
  lwe_pred_df <- lapply(
    seq_along(sampled_days),
    function(i) {
      d0 <- sampled_days[i] %>% as.Date()
      d0_tag <- base::format(d0, "%y%m%d")

      # browser()
      pred_samples_fname <- file.path(
        sample_path,
        sprintf(
          "pred_samples_%s.%s",
          d0_tag,
          extension
        )
      )
      cat(
        "--------------------------------------------------------------------\n"
      )
      cat("Reading prediction samples from:", pred_samples_fname, "\n")
      full_file <- readRDS(pred_samples_fname)

      cat(
        "Processing prediction samples for day:",
        format(d0, "%Y-%m-%d"),
        "\n"
      )
      mod_codes <- names(full_file)
      mclapply(
        X = seq_along(mod_codes),
        FUN = function(j) {
          # browser()
          cat(
            "Processing model:",
            mod_codes[j],
            "\n"
          )
          dt <- as.data.table(full_file[[j]]$sample_df)

          dt[, `:=`(
            model = mod_codes[j],
            date = d0
          )]

          setorder(dt, sim, coord_id, time)

          dt[, below := fit < pow_threshold]

          dt[, run_id := rleid(below), by = .(sim, coord_id)]

          res <- dt[
            below == TRUE & time >= d0, # only lwe from oos
            .(
              model = first(model),
              start_time = first(time),
              end_time = last(time),
              duration_h = as.numeric(difftime(
                last(time),
                first(time),
                units = "hours"
              )) +
                1
            ),
            by = .(sim, coord_id, run_id)
          ][, .(sim, model, coord_id, start_time, duration_h)]

          return(res)
        },
        mc.cores = mc
      ) %>%
        bind_rows()
    }
  ) %>%
    bind_rows()
  # )

  samp_end_time <- Sys.time()
  cat(
    "Time taken to process prediction samples:",
    difftime(samp_end_time, samp_start_time, units = "secs"),
    "seconds\n",
    "Average time taken per day:",
    difftime(samp_end_time, samp_start_time, units = "secs") /
      length(sampled_days),
    "seconds\n"
  )
  # ~ 3mins per day per model
  # 5 models, 15 days is ~ 3.75 hours
  cat("Saving low wind events prediction samples to:", lwe_pred_fname, "\n")
  saveRDS(lwe_pred_df, lwe_pred_fname)
} else {
  cat("Loading existing low wind events prediction samples\n")
  lwe_pred_df <- readRDS(lwe_pred_fname)
}
## 4.2 new locations samples ####
cat("--------------------------------------------------------------------\n")
cat("Low wind events in model samples for new locations \n")
cat("--------------------------------------------------------------------\n")

# read samples

lwe_newloc_fname <- file.path(
  sample_path,
  sprintf(
    "lwe_newloc_%s_t%s.%s",
    d0_tag,
    pow_threshold_label,
    extension
  )
)

if (!file.exists(lwe_newloc_fname) | override_objects) {
  if (!file.exists(lwe_newloc_fname)) {
    cat("LWE processed file not found, creating new summary\n")
  } else {
    cat(
      "LWE processed file found, but override_objects is TRUE. Recreating summary\n"
    )
  }
  cat("This may take a while...\n")
  samp_start_time <- Sys.time()
  lwe_newloc_df <- lapply(
    seq_along(sampled_days),
    function(i) {
      d0 <- sampled_days[i] %>% as.Date()
      d0_tag <- base::format(d0, "%y%m%d")

      # browser()
      pred_samples_fname <- file.path(
        sample_path,
        sprintf(
          "pred_samples_spaceoos_%s.%s",
          d0_tag,
          extension
        )
      )
      cat("Reading prediction samples from:", pred_samples_fname, "\n")
      full_file <- readRDS(pred_samples_fname)

      cat(
        "Processing prediction samples for day:",
        format(d0, "%Y-%m-%d"),
        "\n"
      )
      mod_codes <- names(full_file)
      mclapply(
        X = seq_along(mod_codes),
        FUN = function(j) {
          # browser()

          dt <- as.data.table(full_file[[j]]$sample_df)

          dt[, `:=`(
            model = mod_codes[j],
            date = d0
          )]

          setorder(dt, sim, coord_id, time)

          dt[, below := fit < pow_threshold]

          dt[, run_id := rleid(below), by = .(sim, coord_id)]

          res <- dt[
            below == TRUE, # only lwe from oos
            .(
              model = first(model),
              start_time = first(time),
              end_time = last(time),
              duration_h = as.numeric(difftime(
                last(time),
                first(time),
                units = "hours"
              )) +
                1
            ),
            by = .(sim, coord_id, run_id)
          ][, .(sim, model, coord_id, start_time, duration_h)]

          res
        },
        mc.cores = mc
      ) %>%
        bind_rows()
    }
  ) %>%
    bind_rows()

  samp_end_time <- Sys.time()
  cat(
    "Time taken to process new location samples:",
    difftime(samp_end_time, samp_start_time, units = "secs"),
    "seconds\n",
    "Average time taken per day:",
    difftime(samp_end_time, samp_start_time, units = "secs") /
      length(sampled_days),
    "seconds\n"
  )
  cat("Saving low wind events new location samples to:", lwe_newloc_fname, "\n")
  saveRDS(lwe_newloc_df, lwe_newloc_fname)
} else {
  cat("Loading existing low wind events new location samples\n")
  lwe_newloc_df <- readRDS(lwe_newloc_fname)
}

# 5. Figures ####
cat("--------------------------------------------------------------------\n")
cat("Figures\n")
cat("--------------------------------------------------------------------\n")


# distribution of LWE
dur_lim <- 25
plot_df <-
  lwe_pred_df %>%
  bind_rows(lwe_obs_pred) %>%
  filter(duration_h < dur_lim) %>%
  mutate(bin = floor(duration_h)) %>% # matches binwidth = 1
  count(model, bin) %>%
  group_by(model) %>%
  mutate(percent = 100 * n / sum(n))

plot_df %>%
  ggplot(aes(bin, percent)) +
  geom_col(fill = "steelblue", alpha = 0.6, width = 1) +
  facet_wrap(~model, scales = "free_y", labeller = as_labeller(mod_labels)) +
  labs(
    x = "Duration (hours)",
    y = "Percentage (%)"
  ) +
  theme_minimal()

ggsave(
  sprintf(
    "fig/%s/pred_samples_lwe_duration_dist_t%s.pdf",
    batch_name,
    pow_threshold_label
  ),
  width = 10,
  height = 6
)


obs <- lwe_obs_pred %>%
  filter(
    model == "observed",
    duration_h < dur_lim
  ) %>%
  pull(duration_h)

probs <- seq(0, 1, length.out = 100)

qq_df <- lwe_pred_df %>%
  filter(duration_h < dur_lim) %>%
  group_by(model) %>%
  summarise(
    obs_q = list(quantile(obs, probs = probs, na.rm = TRUE)),
    model_q = list(quantile(duration_h, probs = probs, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  unnest(c(obs_q, model_q))

ggplot(qq_df, aes(x = obs_q, y = model_q)) +
  geom_point(size = 1) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  facet_wrap(~model, scales = "free", labeller = as_labeller(mod_labels)) +
  labs(
    title = "Q-Q Plot of Low Wind Event Durations",
    x = "Observed quantiles",
    y = "Model quantiles"
  ) +
  # coord_equal() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  sprintf(
    "fig/%s/pred_samples_lwe_duration_qq_t%s.pdf",
    batch_name,
    pow_threshold_label
  ),
  width = 10,
  height = 6
)


# distribution of LWE
dur_lim <- 100
plot_df <-
  lwe_newloc_df %>%
  bind_rows(lwe_obs_newloc) %>%
  filter(duration_h < dur_lim) %>%
  mutate(bin = floor(duration_h)) %>% # matches binwidth = 1
  count(model, bin) %>%
  group_by(model) %>%
  mutate(percent = 100 * n / sum(n))

plot_df %>%
  ggplot(aes(bin, percent)) +
  geom_col(fill = "steelblue", alpha = 0.6, width = 1) +
  facet_wrap(~model, scales = "free_y", labeller = as_labeller(mod_labels)) +
  labs(
    x = "Duration (hours)",
    y = "Percentage (%)"
  ) +
  theme_minimal()

ggsave(
  sprintf(
    "fig/%s/newloc_samples_lwe_duration_dist_t%s.pdf",
    batch_name,
    pow_threshold_label
  ),
  width = 10,
  height = 6
)


obs <- lwe_obs_newloc %>%
  filter(
    model == "observed",
    duration_h < dur_lim
  ) %>%
  pull(duration_h)

probs <- seq(0, 1, length.out = 100)

qq_df <- lwe_newloc_df %>%
  filter(duration_h < dur_lim) %>%
  group_by(model) %>%
  summarise(
    obs_q = list(quantile(obs, probs = probs, na.rm = TRUE)),
    model_q = list(quantile(duration_h, probs = probs, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  unnest(c(obs_q, model_q))

ggplot(qq_df, aes(x = obs_q, y = model_q)) +
  geom_point(size = 1) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  facet_wrap(~model, scales = "free", labeller = as_labeller(mod_labels)) +
  labs(
    title = "Q-Q Plot of Low Wind Event Durations",
    x = "Observed quantiles",
    y = "Model quantiles"
  ) +
  # coord_equal() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  sprintf(
    "fig/%s/newloc_samples_lwe_duration_qq_t%s.pdf",
    batch_name,
    pow_threshold_label
  ),
  width = 10,
  height = 6
)


# 6. score metrics ####

endtime <- Sys.time()

cat("--------------------------------------------------------------------\n")
cat("Finished whole process\n")
cat("--------------------------------------------------------------------\n")

timediff <- difftime(endtime, starttime, units = "auto")

cat(sprintf(
  "Total time taken: %.2f %s\n",
  timediff,
  units(timediff)
))
