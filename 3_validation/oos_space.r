# 0. Setup ####
cat(
  "--------------------------------------------------------------------\n"
)
cat("Running validation for out-of-sample spatial prediction\n")
cat(
  "--------------------------------------------------------------------\n"
)
local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

# 0.1 global parameter #####
day_id <- 1
# mesh_edge_par <- 20 # km, target edge length for the spatial mesh. 10 is fine, 20 is coarse but faster
override_objects <- TRUE
rerun_samples <- TRUE
prec_init <- log(200)
task_prefix <- "spaceoos"

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
  rerun_samples <- as.logical(args[3])
}
# 0.2 libraries and paths ####
require(parallel)

if (local_run) {
  data_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
  sample_path <- "~/Documents/elexon/samples"
  n_samp <- 100
  pixel_dims <- c(150, 150)
} else {
  data_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  gen_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  model_path <- "/exports/eddie/scratch/s2441782/calibration/model_objects"
  sample_path <- "/exports/eddie/scratch/s2441782/calibration/samples"
  temp_lib <- "/exports/eddie3_homes_local/s2441782/lib"
  pixel_dims <- c(300, 300)
  n_samp <- 1000
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
# require(ggspatial)

source("aux_funct.R")

sampled_days <- c("2020-08-14", "2024-04-17", "2024-04-12")
d0 <- sampled_days[day_id] %>% as.Date()
d0_tag <- base::format(d0, "%y%m%d")

n.days <- 0
n.days.before <- 7
n.hours <- 1
t0 <- d0 - n.days.before
t1 <- d0 + n.days
th <- t1 + hours(n.hours)


# Read models ####
mod_labels <- c(
  # "Generic PC",
  "Linear model",
  "GB LM",
  "QM",
  "Spatio-temporal coarse",
  # "Spatio-temporal fine",
  "1D SPDE model",
  "AR1 model",
  "AR2 model"
)
est_cols <- c(
  # "norm_power_est0",
  "lm",
  "agg_lm",
  "qm",
  "st0_m1",
  # "st0_m2",
  "spde1d",
  "ar1",
  "ar2"
)

n_models <- length(est_cols)
names(mod_labels) <- est_cols

mod_vec <- list.files(model_path, pattern = d0_tag, full.names = TRUE)
mod_vec <- mod_vec[!grepl("spatial|fine", mod_vec)] %>% sort() # exclude meshes from st model

model_df <- tibble(
  label = mod_labels,
  code = est_cols,
  fname = mod_vec
) %>%
  mutate(
    type = case_when(
      grepl("lm", code) ~ "lm",
      grepl("qm", code) ~ "qm",
      TRUE ~ "bru"
    )
  )

# Predictions for next hours ####

## prediction df ####
gb_day_df_fname <- sprintf("data/GB_daily_summary_%s.parquet", d0_tag)

if (!file.exists(gb_day_df_fname) || override_objects) {
  if (!file.exists(gb_day_df_fname)) {
    cat("GB daily summary file not found, creating new summary\n")
  } else {
    cat(
      "GB daily summary file found, but override_objects is TRUE. Recreating summary\n"
    )
  }
  GB_df <- read_parquet(file.path(gen_path, "GB_aggr.parquet")) %>%
    rename(time = halfHourEndTime) %>%
    mutate(
      err = norm_power_est0 - norm_potential,
      error0 = norm_potential - norm_power_est0,
      date = as.Date(time)
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
} else {
  cat("Loading existing GB daily summary\n")
  gb_day_df <- read_parquet(gb_day_df_fname)
}

# extension <- ifelse(local_run, "gpkg", "geojson")
extension <- "rds"
df_pattern <- sprintf(
  "^calibration_preddf_%s_.*_%s\\.%s$",
  task_prefix,
  d0_tag,
  extension
)
files_found <- list.files("data", pattern = df_pattern, full.names = TRUE)

coord_list_fname <- "data/coord_list.csv"

cat("Loading existing coordinate list\n")
coord_list <- read.csv(coord_list_fname)


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

  pwr_curv_df <- read_parquet(file.path(
    gen_path,
    "power_curve_all_enriched.parquet"
  ))

  wf_df_pred <- pwr_curv_df %>%
    rename(time = halfHourEndTime) %>%
    mutate(
      date = as.Date(time),
      elevation = pmax(0, elevation),
      site_name = site_name %>%
        gsub("\\b(wind\\s*farm|wf)\\b", "", ., ignore.case = TRUE) %>%
        trimws()
    ) %>%
    # filter(date %in% sampled_days) %>%
    # filter(date >= d0, date <= d0 + n.days - 1) %>%
    filter(time >= t0, time < th) %>%
    filter(coord_id %in% coord_list$coord_id[!coord_list$sampled]) %>%
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
      date = as.Date(time)
    ) %>%
    left_join(
      gb_day_df %>% dplyr::select(date, p_group3),
      by = c("date" = "date")
    )

  x <- wf_df_pred$pow_group %>% unique() %>% sort()
  min_jump <- min(diff(sort(x))) / diff(range(x))
  if (min_jump <= 1e-4) {
    wf_df_pred <- wf_df_pred %>%
      mutate(pow_group = inla.group(norm_power_est0, n = 20, method = "cut"))
  }

  cat("Converting coordinates to km\n")
  wf_df_pred <- wf_df_pred %>%
    st_geometry() %>%
    (\(g) g / 1000)() %>%
    st_set_geometry(wf_df_pred, .)

  wf_df_fname <- sprintf(
    "data/calibration_preddf_%s_%s.%s",
    "base",
    d0_tag,
    extension
  )
  saveRDS(
    wf_df_pred,
    wf_df_fname
  )
}
n_loc <- nrow(wf_df_pred %>% distinct(x, y))
cat("Number of unique locations:", n_loc, "\n")
n <- nrow(wf_df_pred)
cat("Number of records in the dataset:", n, "\n")
## linear models ####

lm_df <- model_df %>% filter(type == "lm")

lm_df %>% pull(fname) %>% map(readRDS) -> mod_list
names(mod_list) <- lm_df %>% pull(code)

lm_pred <- lapply(
  names(mod_list),
  function(mod) {
    predict(mod_list[[mod]], newdata = wf_df_pred, interval = "prediction") %>%
      as.data.frame() %>%
      rename(
        estimate = fit,
        lwr = lwr,
        upr = upr
      ) %>%
      bind_cols(
        wf_df_pred %>%
          dplyr::select(
            coord_id,
            site_name,
            time,
            date,
            norm_potential,
            norm_power_est0,
            capacity,
            tech_typ,
            p_group3
          ),
        .
      ) %>%
      mutate(
        estimate = pmin(1, pmax(0, estimate)),
        lwr = pmin(1, pmax(0, lwr)),
        upr = pmin(1, pmax(0, upr)),
        std_error = summary(mod_list[[mod]])$sigma,
        model = mod
      )
  }
) %>%
  bind_rows()


lm_pred_fig_df <- lm_pred %>%
  st_drop_geometry() %>%
  group_by(time, model) %>%
  summarise(
    mean = sum(estimate * capacity) / sum(capacity),
    std_error = if (first(model) == "lm") {
      sqrt(sum((capacity / sum(capacity))^2 * std_error^2))
    } else {
      mean(std_error)
    },
    norm_potential = sum(norm_potential * capacity) / sum(capacity),
    .groups = "drop"
  ) %>%
  mutate(
    lwr = mean - 1.96 * std_error,
    upr = mean + 1.96 * std_error,
    lwr = pmin(1, pmax(0, lwr)),
    upr = pmin(1, pmax(0, upr))
  )

# lm_pred_fig_df %>%
#   ggplot() +
#   geom_ribbon(
#     aes(
#       x = time,
#       ymin = lwr,
#       ymax = upr
#     ),
#     fill = blues9[5],
#     alpha = 0.5
#   ) +
#   geom_line(
#     aes(x = time, y = mean),
#     color = blues9[9],
#     lwd = 1
#   ) +
#   geom_line(
#     aes(x = time, y = norm_potential),
#     color = "darkred",
#     lwd = 1
#   ) +
#   coord_cartesian(ylim = c(0, 1)) +
#   facet_wrap(~model, nrow = 2)

## quantile mapping ####

## bru models ####

bru_df <- model_df %>% filter(type == "bru")
# mod_temp <- bruar1
# mod_temp <- readRDS(bru_df$fname[1])
# test <- bru_ci_plot(
#   bru_model = mod_temp,
#   newdata = wf_df_pred,
#   n.samples = 500,
#   show.fig = TRUE
# )

# test$GB_summary %>%
#   # filter(time >= t1) %>%
#   ggplot() +
#   geom_ribbon(
#     aes(
#       x = time,
#       ymin = lwr,
#       ymax = upr
#     ),
#     fill = blues9[5],
#     alpha = 0.5
#   ) +
#   geom_line(
#     aes(x = time, y = mean),
#     color = blues9[9],
#     lwd = 1
#   ) +
#   geom_line(
#     aes(x = time, y = norm_potential),
#     color = "darkred",
#     lwd = 1
#   ) +
#   # coord_cartesian(ylim = c(0, 1))+
#   scale_x_datetime()

# test$wf_summary %>%
#   filter(coord_id %in% c(0 + 1:30)) %>%
#   ggplot() +
#   geom_ribbon(
#     aes(
#       x = time,
#       ymin = lwr,
#       ymax = upr
#     ),
#     fill = blues9[5],
#     alpha = 0.5
#   ) +
#   geom_line(
#     aes(x = time, y = fit),
#     color = blues9[9],
#     lwd = 1
#   ) +
#   geom_line(
#     aes(x = time, y = norm_potential),
#     color = "darkred",
#     lwd = 1
#   ) +
#   facet_wrap(~site_name, scales = "free_y") +
#   coord_cartesian(ylim = c(0, 1)) +
#   scale_x_datetime(date_labels = "%H:%M")

# ?predict.bru

# mod_temp %>% summary()
# mod_temp$.args$control.family[[1]]$hyper$theta1$fixed

# lin_pred <- get_bru_formula(mod_temp)

source("aux_funct.R")
# mod_temp <- bruar1

### loop throgh all models #####

# pred band summary
pred_summary_fname <- sprintf(
  "summaries/pred_band_%s_summary_%s.rds",
  task_prefix,
  d0_tag
)
pred_samples_fname <- file.path(
  sample_path,
  sprintf(
    "pred_samples_%s_%s.rds",
    task_prefix,
    d0_tag
  )
)
if (!file.exists(pred_summary_fname) || rerun_samples) {
  if (!file.exists(pred_summary_fname)) {
    cat("Prediction band summary file not found, creating new summary\n")
  } else {
    cat(
      "Prediction band summary file found, but rerun_samples is TRUE. Recreating summary\n"
    )
  }
  pred_band_summary <- lapply(
    seq_along(bru_df$fname),
    function(i) {
      cat(
        "--------------------------------------------------------------------\n"
      )
      cat("Processing model:", bru_df$label[i], "\n")
      cat(
        "--------------------------------------------------------------------\n"
      )
      mod_temp <- readRDS(bru_df$fname[i])
      test <- bru_ci_plot(
        bru_model = mod_temp,
        newdata = wf_df_pred,
        n.samples = n_samp,
        show.fig = FALSE
      )
      test
    }
  )
  summary_only <- lapply(
    pred_band_summary,
    function(x) {
      list(
        GB_summary = x$GB_summary,
        wf_summary = x$wf_summary,
        formla = x$formula
      )
    }
  )
  names(pred_band_summary) <- bru_df$code
  names(summary_only) <- bru_df$code
  saveRDS(summary_only, pred_summary_fname)
  saveRDS(pred_band_summary, pred_samples_fname)
} else {
  cat("Loading existing prediction band summary\n")
  pred_band_summary <- readRDS(pred_summary_fname)
}
# pred_band_summary %>% lapply(., \(z) z$formula)

## Consolidated figures #####
### GB aggregation summary ####
gb_fig_df <- bind_rows(
  lm_pred_fig_df,
  lapply(
    bru_df$code,
    function(code) {
      pred_band_summary[[code]]$GB_summary %>%
        mutate(
          model = code
        )
    }
  ) %>%
    bind_rows()
) %>%
  mutate(
    oos = ifelse(time >= t1, TRUE, FALSE)
  )
saveRDS(
  gb_fig_df,
  sprintf("summaries/spaceoos_GB_fig_band_summary_%s.rds", d0_tag)
)

gb_fig_df %>%
  filter(time >= t1) %>%
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
    aes(x = time, y = mean, col = "fit"),
    # color = blues9[9],
    lwd = 1
  ) +
  geom_line(
    aes(x = time, y = norm_potential, col = "observed"),
    # color = "darkred",
    lwd = 1
  ) +
  # coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(~model, nrow = 2, labeller = as_labeller(mod_labels)) +
  scale_x_datetime(date_labels = "%H:%M") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("fit" = blues9[9], "observed" = "darkred")) +
  scale_fill_manual(values = c("95% CI" = blues9[5])) +
  labs(fill = "", color = "")

ggsave(
  filename = sprintf("fig/GB_pred_band_%s.pdf", d0_tag),
  width = 10,
  height = 6,
  # dpi = 300
)

### WF level summary ####

wf_fig_df <- bind_rows(
  lm_pred %>%
    dplyr::select(
      coord_id,
      site_name,
      time,
      norm_potential,
      # norm_power_est0,
      # capacity,
      model,
      estimate,
      lwr,
      upr
    ) %>%
    st_drop_geometry() %>%
    rename(fit = estimate),
  lapply(
    bru_df$code,
    function(code) {
      pred_band_summary[[code]]$wf_summary %>%
        mutate(
          model = code
        )
    }
  ) %>%
    bind_rows()
) %>%
  mutate(
    oos = ifelse(time >= t1, TRUE, FALSE)
  )

saveRDS(
  wf_fig_df,
  sprintf("summaries/WF_fig_band_summary_%s.rds", d0_tag)
)

for (mod in est_cols[!grepl("qm", est_cols)]) {
  for (k in 0:2) {
    # print(k)
    wf_fig_df %>%
      filter(model == mod) %>%
      filter(coord_id %in% c(k * 40 + 1:40)) %>%
      filter(time >= t1) %>%
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
        aes(x = time, y = fit, col = "fit"),
        # color = blues9[9],
        lwd = 1
      ) +
      geom_line(
        aes(x = time, y = norm_potential, col = "observed"),
        # color = "darkred",
        lwd = 1
      ) +
      facet_wrap(~site_name, scales = "free_y") +
      coord_cartesian(ylim = c(0, 1)) +
      scale_x_datetime(date_labels = "%H:%M") +
      theme(legend.position = "bottom") +
      scale_color_manual(
        values = c("fit" = blues9[9], "observed" = "darkred")
      ) +
      scale_fill_manual(values = c("95% CI" = blues9[5])) +
      labs(fill = "", color = "")
    ggsave(
      filename = sprintf("fig/WF_pred_band_%s_%s_%d.pdf", mod, d0_tag, k + 1),
      width = 10,
      height = 6,
      # dpi = 300
    )
  }
}
# wf_fig_df %>%
#   filter(model == "st0_m1") %>%
#   filter(coord_id %in% c(80 + 1:40)) %>%
#   filter(time >= t1) %>%
#   ggplot() +
#   geom_ribbon(
#     aes(
#       x = time,
#       ymin = lwr,
#       ymax = upr
#     ),
#     fill = blues9[5],
#     alpha = 0.5
#   ) +
#   geom_line(
#     aes(x = time, y = fit),
#     color = blues9[9],
#     lwd = 1
#   ) +
#   geom_line(
#     aes(x = time, y = norm_potential),
#     color = "darkred",
#     lwd = 1
#   ) +
#   facet_wrap(~site_name, scales = "free_y") +
#   # coord_cartesian(ylim = c(0, 1)) +

#   scale_x_datetime(date_labels = "%H:%M")
# test$wf_summary %>%
#   filter(coord_id %in% c(0 + 1:30)) %>%
#   ggplot() +
#   geom_ribbon(
#     aes(
#       x = time,
#       ymin = lwr,
#       ymax = upr
#     ),
#     fill = blues9[5],
#     alpha = 0.5
#   ) +
#   geom_line(
#     aes(x = time, y = fit),
#     color = blues9[9],
#     lwd = 1
#   ) +
#   geom_line(
#     aes(x = time, y = norm_potential),
#     color = "darkred",
#     lwd = 1
#   ) +
#   facet_wrap(~site_name, scales = "free_y") +
#   coord_cartesian(ylim = c(0, 1)) +
#   scale_x_datetime(date_labels = "%H:%M")
## Bands coverage by model ####
### wf level
cov_bands_wf <- wf_fig_df %>%
  filter(time >= t1) %>%
  group_by(model, coord_id) %>%
  summarise(
    coverage = mean(norm_potential >= lwr & norm_potential <= upr),
    .groups = "drop"
  ) %>%
  group_by(model) %>%
  summarise(
    mean_coverage = mean(coverage),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_coverage)) %>%
  mutate(
    model = factor(model, levels = model)
  )
cov_bands_wf %>%
  ggplot(aes(x = model, y = mean_coverage)) +
  geom_col(fill = blues9[5]) +
  geom_text(aes(label = round(mean_coverage, 3)), vjust = -0.5) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Mean coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf("fig/WF_pred_band_coverage_%s.pdf", d0_tag),
  width = 10,
  height = 6,
  # dpi = 300
)
### aggregated #####
cov_bands <- gb_fig_df %>%
  filter(time >= t1) %>%
  group_by(model) %>%
  summarise(
    coverage = mean(norm_potential >= lwr & norm_potential <= upr),
    .groups = "drop"
  ) %>%
  arrange(desc(coverage)) %>%
  mutate(
    model = factor(model, levels = model)
  )

cov_bands %>%
  ggplot(aes(x = model, y = coverage)) +
  geom_col(fill = blues9[5]) +
  geom_text(aes(label = round(coverage, 3)), vjust = -0.5) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf("fig/GB_pred_band_coverage_%s.pdf", d0_tag),
  width = 10,
  height = 6,
  # dpi = 300
)
## Correlation exploration ####
# gb_fig_df %>%
#   filter(model %in% c("spde1d", "ar2", "ar1", "st0_m1", "st0_m2")) %>%
#   group_by(model) %>%
#   mutate(
#     error = norm_potential - mean,
#   ) %>%
#   summarise(
#     cor = cor(error, use = "complete.obs")
#   )
gb_fig_df <- gb_fig_df %>%
  mutate(residual = norm_potential - mean)

var_emp <- gb_fig_df %>%
  group_by(model) %>%
  summarise(var_res = var(residual), sd_res = sd(residual), .groups = "drop")


var_wf <- wf_fig_df %>%
  group_by(model) %>%
  summarise(
    var_res = var(norm_potential - fit),
    sd_res = sd(norm_potential - fit),
    .groups = "drop"
  )

var_wf %>%
  mutate(
    scaled_var = var_res / n_loc
  ) %>%
  left_join(var_emp, by = "model", suffix = c("_wf", "_emp")) %>%
  mutate(
    rho_est = (var_res_emp - var_res_wf / n_loc) /
      (var_res_wf * (1 - 1 / n_loc))
  )
