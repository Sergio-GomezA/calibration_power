local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

pow_threshold <- 0.05
pow_threshold_label <- gsub("\\.", "_", as.character(pow_threshold))


override_objects <- FALSE
# rerun_samples <- FALSE
# prec_init <- log(200)
# batch_name <- "batch2025"
batch_name <- "batchY25d150"


if (local_run) {
  cat("Running in local mode\n")
} else {
  cat("Running in cluster mode\n")
}


# 0.2 libraries and paths ####
require(parallel)

if (local_run) {
  data_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
  extra_path <- "~/Documents/elexon/extra"
  output_path <- "~/Documents/elexon/caloutput"
  sample_path <- "~/Documents/elexon/samples"
  n_samp <- 100
  pixel_dims <- c(150, 150)
} else {
  data_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  gen_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  model_path <- "/exports/eddie/scratch/s2441782/calibration/model_objects"
  extra_path <- "/exports/eddie/scratch/s2441782/calibration/extra"
  output_path <- "/exports/eddie/scratch/s2441782/calibration/caloutput"
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
require(kableExtra)
# require(ggspatial)

source("aux_funct.R")
mc <- ifelse(local_run, 1, available_cores())

# sampled_days_df <- read.csv("data/sample_days_df.csv") %>%
#   mutate(date = as.Date(date))

sampled_days_df <- read.csv("data/sample_days_df_25_n150.csv") %>%
  mutate(date = as.Date(date))

# sampled_days_df <- read.csv("data/sample_days_df.csv") %>%
#   mutate(date = as.Date(date))

sampled_days <- sampled_days_df %>%
  pull(date) %>%
  sort()

extension <- "rds"
coord_list_fname <- "data/coord_list.csv"

cat("Loading existing coordinate list\n")
coord_list <- read.csv(coord_list_fname)


gb_day_df_fname <- sprintf("data/GB_daily_summary.parquet")

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

GB_df <- read_parquet(file.path(gen_path, "GB_aggr.parquet")) %>%
  rename(time = halfHourEndTime) %>%
  mutate(
    err = norm_power_est0 - norm_potential,
    error0 = norm_potential - norm_power_est0,
    date = as.Date(time)
  )
pacf_df <- GB_df %>%
  arrange(tech_typ, time) %>%
  group_by(tech_typ) %>%
  summarise(
    pacf = list(
      pacf(error0, plot = FALSE, na.action = na.pass)
    ),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lag = map(pacf, ~ .x$lag),
    acf = map(pacf, ~ .x$acf),
    ci = 1.96 / sqrt(n)
  ) %>%
  dplyr::select(-pacf) %>%
  unnest(c(lag, acf))

head(pacf_df)
pacf_df %>%
  ggplot(aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_segment(aes(xend = lag, yend = 0)) +
  facet_wrap(~tech_typ, scales = "free_y") +
  labs(
    x = "Lag",
    y = "Partial autocorrelation"
  ) +
  theme_bw() +
  coord_cartesian(ylim = c(-0.4, 1)) +
  scale_y_continuous(breaks = seq(-0.3, 1, by = 0.2))
ggsave(
  sprintf("fig/%s/GB_pacf.pdf", batch_name),
  width = 10,
  height = 6,
  create.dir = TRUE
)


# Read models ####
model_catalog <- read.csv("data/model_catalog.csv") %>%
  na.omit()

# if (local_run) {
#   model_catalog <- model_catalog %>%
#     # filter(!grepl("fine", mod_labels)) %>%
#     filter(!grepl("st0", est_cols))
# }
mod_labels <- model_catalog$mod_labels
est_cols <- model_catalog$est_cols
n_models <- length(est_cols)
names(mod_labels) <- est_cols
# excluded_models0 <- c("lm_bru")
# excluded_models <- c("lm_bru", "qm")

excluded_models0 <- c("")
excluded_models <- c("", "qm")

model_catalog <- read.csv("data/model_catalog.csv") %>%
  na.omit()
model_df <- model_catalog %>%
  rename(code = est_cols, label = mod_labels) %>%
  arrange(desc(nchar(mode_code_prefix))) %>%
  mutate(
    type = case_when(
      grepl("lm_", code) ~ "bru",
      grepl("lm", code) ~ "lm",
      grepl("qm", code) ~ "qm",
      TRUE ~ "bru"
    )
  ) %>%
  filter(!code %in% c("st0_m0", "st0_m1"))
bru_df <- model_df %>% filter(type == "bru")

# Time oos ####
## read summary tables of prediction bands ######
cat("Reading summary tables for GB prediction bands\n")
gb_fig_df <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    # print(i)
    d0_tag <- base::format(d0, "%y%m%d")

    file_name <- sprintf(
      "%s/%s/summaries/oos/GB_fig_band_summary_time_%s.rds",
      output_path,
      batch_name,
      d0_tag
    )
    if (!file.exists(file_name)) {
      cat("File not found:", file_name, "\n")
      return(NULL)
    }
    readRDS(file_name) %>%
      mutate(
        date = as.Date(time)
      )
  }
) %>%
  bind_rows() %>%
  left_join(
    gb_day_df %>% dplyr::select(date, p_group3) %>% rename(pgroup3 = p_group3),
    by = "date"
  )

cat("Reading summary tables for WF level prediction bands\n")
wf_fig_df <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    d0_tag <- base::format(d0, "%y%m%d")
    file_name <- sprintf(
      "%s/%s/summaries/oos/WF_fig_band_summary_time_%s.rds",
      output_path,
      batch_name,
      d0_tag
    )
    if (!file.exists(file_name)) {
      cat("File not found:", file_name, "\n")
      return(NULL)
    }
    readRDS(file_name) %>%
      mutate(
        date = as.Date(time)
      )
  }
) %>%
  bind_rows() %>%
  left_join(
    gb_day_df %>% dplyr::select(date, p_group3) %>% rename(pgroup3 = p_group3),
    by = "date"
  )

## calibration fit scatter####
gb_fig_df %>%
  filter(oos) %>%
  filter(!model %in% excluded_models0) %>%
  ggplot(aes(x = norm_potential, y = mean)) +
  geom_hex() +
  # geom_point(alpha = 0.5) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkred"
  ) +
  facet_wrap(~model, labeller = as_labeller(mod_labels)) +
  theme_bw() +
  # scale_color_lancet() +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  labs(
    x = "Observed power",
    y = "Predicted power",
    # col = "regime"
  )

ggsave(
  sprintf("fig/%s/GB_fit_scatter_oos.png", batch_name),
  width = 10,
  height = 6,
  dpi = 100
)

wf_fig_df %>%
  slice_sample(n = 1e6) %>%
  filter(!model %in% excluded_models0) %>%
  mutate(hour = hour(time)) %>%
  group_by(date) %>%
  mutate(leadh = difftime(time, min(time), units = "hours")) %>%
  ungroup() %>%
  filter(oos, leadh <= 1) %>%
  # filter(oos, between(hour, 10, 14)) %>%
  # filter(oos) %>%
  ggplot(aes(x = norm_potential, y = fit)) +
  geom_hex() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkred"
  ) +
  facet_wrap(~model, labeller = as_labeller(mod_labels)) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  # scale_fill_gradient(
  #   trans = "log10",
  #   low = "grey90",
  #   high = blues9[5],
  #   # limits = c(5, NA),
  # ) +
  theme_bw() +
  labs(
    x = "Normalised potential power",
    y = "Normalised predicted power",
    fill = "Count"
  )
# scales::show_col(blues9)
ggsave(
  sprintf("fig/%s/WF_fit_scatter_oos.png", batch_name),
  width = 10,
  height = 6,
  dpi = 100
)


## Coverage bands #####

### wf level ####

cov_bands_wf <- wf_fig_df %>%
  filter(!model %in% excluded_models) %>%
  filter(oos) %>%
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
  geom_col(fill = blues9[7]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkred") +
  geom_text(aes(label = round(mean_coverage, 3)), vjust = -0.5) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Mean coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf("fig/%s/WF_pred_band_coverage.pdf", batch_name),
  width = 10,
  height = 6,
  # dpi = 300
)
### aggregated #####
cov_bands <- gb_fig_df %>%
  filter(oos) %>%
  filter(!model %in% excluded_models) %>%
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
  geom_col(fill = blues9[7]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkred") +
  geom_text(aes(label = round(coverage, 3)), vjust = -0.5) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf("fig/%s/GB_pred_band_coverage.pdf", batch_name),
  width = 10,
  height = 6,
  # dpi = 300
)

### by tech type #####
sample_loc_fname <- "data/coord_list_wloc.csv"
cat("Loading existing coordinate list with location names\n")
loc_cat <- read.csv(sample_loc_fname)
cov_bands <- wf_fig_df %>%
  filter(!model %in% excluded_models) %>%
  left_join(loc_cat %>% dplyr::select(coord_id, tech_typ), by = "coord_id") %>%
  filter(oos) %>%
  group_by(model, tech_typ) %>%
  summarise(
    coverage = mean(norm_potential >= lwr & norm_potential <= upr),
    .groups = "drop"
  ) %>%
  arrange(desc(coverage)) %>%
  mutate(
    model = factor(model, levels = model %>% unique())
  )

cov_bands %>%
  ggplot(aes(x = model, y = coverage)) +
  geom_col(fill = blues9[7]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkred") +
  geom_text(aes(label = round(coverage, 3)), vjust = -0.5) +
  facet_wrap(~tech_typ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf("fig/%s/GB_pred_band_coverage_tech.pdf", batch_name),
  width = 10,
  height = 6,
  # dpi = 300
)
### by regime #####
cov_bands <- wf_fig_df %>%
  # left_join(loc_cat %>% dplyr::select(coord_id, tech_typ), by = "coord_id") %>%
  filter(oos) %>%
  filter(!model %in% excluded_models) %>%
  group_by(model, pgroup3) %>%
  summarise(
    coverage = mean(norm_potential >= lwr & norm_potential <= upr),
    .groups = "drop"
  ) %>%
  arrange(desc(coverage)) %>%
  mutate(
    model = factor(model, levels = model %>% unique())
  )

cov_bands %>%
  ggplot(aes(x = model, y = coverage)) +
  geom_col(fill = blues9[7]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkred") +
  geom_text(aes(label = round(coverage, 3)), vjust = -0.5) +
  facet_wrap(~pgroup3) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf("fig/%s/GB_pred_band_coverage_regime.pdf", batch_name),
  width = 14,
  height = 6,
  # dpi = 300
)

## Error metrics #####
metrics_table_t <- wf_fig_df %>%
  filter(!model %in% excluded_models0) %>%
  group_by(oos, model) %>%
  summarise(
    RMSE = ModelMetrics::rmse(actual = norm_potential, predicted = fit),
    MAE = ModelMetrics::mae(actual = norm_potential, predicted = fit),
    # MDAPE = mdape(
    #   actual = norm_potential,
    #   predicted = fit,
    #   pos_only = TRUE
    # ),
    Bias = mean(fit - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    model = mod_labels[model],
    oos = ifelse(oos, "OOS", "IS")
  ) %>%
  pivot_wider(
    names_from = oos,
    values_from = c(
      RMSE,
      MAE,
      # MDAPE,
      Bias
    )
  ) %>%
  arrange(desc(RMSE_OOS))
metrics_table_t
metrics_table_t %>%
  mutate(
    across(
      c(RMSE_IS, RMSE_OOS, MAE_IS, MAE_OOS, Bias_IS, Bias_OOS),
      ~ round(., 3)
    ),
    # across(c(MDAPE_IS, MDAPE_OOS), ~ round(., 1))
  ) %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    align = "lcccccccc",
    col.names = c(
      "Model",
      "IS",
      "OOS",
      "IS",
      "OOS",
      # "IS",
      # "OOS",
      "IS",
      "OOS"
    ),
    caption = "Performance metrics for in-sample (IS) and out-of-sample (OOS) predictions."
  ) %>%
  add_header_above(c(
    " " = 1,
    "RMSE" = 2,
    "MAE" = 2,
    # "MDAPE (%)" = 2,
    "Bias" = 2
  )) %>%
  kable_styling(latex_options = "hold_position")

# space oos ####
## read summary tables of prediction bands ######
cat("Reading summary tables for GB prediction bands\n")
gb_fig_df <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    # print(i)
    # browser()
    d0_tag <- base::format(d0, "%y%m%d")
    file_name <- sprintf(
      "%s/%s/summaries/oos/GB_fig_band_summary_spaceoos_%s.rds",
      output_path,
      batch_name,
      d0_tag
    )
    if (!file.exists(file_name)) {
      cat("File not found:", file_name, "\n")
      return(NULL)
    }
    readRDS(file_name) %>%
      # mutate(
      #   pgroup3 = sampled_days_df$p_group3[i] %>%
      #     factor(levels = c("low", "mid", "high"))
      # ) %>%
      mutate(oos = time < as.POSIXct(d0))
  }
) %>%
  bind_rows() %>%
  mutate(date = as.Date(time)) %>%
  left_join(
    gb_day_df %>% dplyr::select(date, p_group3) %>% rename(pgroup3 = p_group3),
    by = "date"
  )
cat("Reading summary tables for WF level prediction bands\n")
wf_fig_df <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    d0_tag <- base::format(d0, "%y%m%d")
    file_name <- sprintf(
      "%s/%s/summaries/oos/WF_fig_band_summary_spaceoos_%s.rds",
      output_path,
      batch_name,
      d0_tag
    )
    if (!file.exists(file_name)) {
      cat("File not found:", file_name, "\n")
      return(NULL)
    }
    readRDS(file_name) %>%
      # mutate(
      #   pgroup3 = sampled_days_df$p_group3[i] %>%
      #     factor(levels = c("low", "mid", "high"))
      # ) %>%
      mutate(oos = time < as.POSIXct(d0))
  }
) %>%
  bind_rows() %>%
  mutate(date = as.Date(time)) %>%
  left_join(
    gb_day_df %>% dplyr::select(date, p_group3) %>% rename(pgroup3 = p_group3),
    by = "date"
  )

## scatter #####
gb_fig_df %>%
  filter(oos) %>%
  filter(!model %in% excluded_models0) %>%
  ggplot(
    # aes(x = norm_potential, y = mean, col = pgroup3)
  ) +
  # geom_point() +
  geom_hex(aes(x = norm_potential, y = mean)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkred"
  ) +
  facet_wrap(~model, labeller = as_labeller(mod_labels)) +
  theme_bw() +
  # scale_color_lancet() +
  # scale_color_manual(values = regime_palette) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  labs(
    x = "Observed power",
    y = "Predicted power",
    # col = "regime"
  )

ggsave(
  sprintf("fig/%s/GB_fit_scatter_oos_space.png", batch_name),
  width = 10,
  height = 6,
  dpi = 100
)

wf_fig_df %>%
  slice_sample(n = 2e6) %>%
  filter(oos) %>%
  filter(!model %in% excluded_models0) %>%
  ggplot(aes(x = norm_potential, y = fit)) +
  geom_hex() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkred"
  ) +
  facet_wrap(~model, labeller = as_labeller(mod_labels)) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  # scale_fill_gradient(
  #   trans = "log10",
  #   low = "grey90",
  #   high = blues9[8]
  # ) +
  theme_bw() +
  labs(
    x = "Normalised potential power",
    y = "Normalised predicted power",
    fill = "Count"
  )

ggsave(
  sprintf("fig/%s/WF_fit_scatter_oos_space.png", batch_name),
  width = 10,
  height = 6,
  dpi = 100
)


## Coverage bands #####
### wf level ####
cov_bands_wf <- wf_fig_df %>%
  filter(!model %in% excluded_models) %>%
  filter(oos) %>%
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
  geom_col(fill = blues9[7]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkred") +
  geom_text(aes(label = round(mean_coverage, 3)), vjust = -0.5) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Mean coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf("fig/%s/WF_pred_band_coverage_spaceoos.pdf", batch_name),
  width = 10,
  height = 6,
  # dpi = 300
)
### aggregated #####
cov_bands <- gb_fig_df %>%
  filter(oos) %>%
  filter(!model %in% excluded_models) %>%
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
  geom_col(fill = blues9[7]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkred") +
  geom_text(aes(label = round(coverage, 3)), vjust = -0.5) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf("fig/%s/GB_pred_band_coverage_spaceoos.pdf", batch_name),
  width = 10,
  height = 6,
  # dpi = 300
)

### techtype #####
cov_bands <- wf_fig_df %>%
  left_join(loc_cat %>% dplyr::select(coord_id, tech_typ), by = "coord_id") %>%
  filter(oos) %>%
  filter(!model %in% excluded_models) %>%
  group_by(model, tech_typ) %>%
  summarise(
    coverage = mean(norm_potential >= lwr & norm_potential <= upr),
    .groups = "drop"
  ) %>%
  arrange(desc(coverage)) %>%
  mutate(
    model = factor(model, levels = model %>% unique())
  )

cov_bands %>%
  ggplot(aes(x = model, y = coverage)) +
  geom_col(fill = blues9[7]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkred") +
  geom_text(aes(label = round(coverage, 3)), vjust = -0.5) +
  facet_wrap(~tech_typ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf(
    "fig/%s/GB_pred_band_coverage_tech_spaceoos.pdf",
    batch_name
  ),
  width = 10,
  height = 6,
  # dpi = 300
)
### by regime #####
cov_bands <- wf_fig_df %>%
  # left_join(loc_cat %>% dplyr::select(coord_id, tech_typ), by = "coord_id") %>%
  filter(oos) %>%
  filter(!model %in% excluded_models) %>%
  group_by(model, pgroup3) %>%
  summarise(
    coverage = mean(norm_potential >= lwr & norm_potential <= upr),
    .groups = "drop"
  ) %>%
  arrange(desc(coverage)) %>%
  mutate(
    model = factor(model, levels = model %>% unique())
  )

cov_bands %>%
  ggplot(aes(x = model, y = coverage)) +
  geom_col(fill = blues9[7]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkred") +
  geom_text(aes(label = round(coverage, 3)), vjust = -0.5) +
  facet_wrap(~pgroup3) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Model", y = "Coverage") +
  scale_x_discrete(labels = mod_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = sprintf(
    "fig/%s/GB_pred_band_coverage_regime_spaceoos.pdf",
    batch_name
  ),
  width = 14,
  height = 6,
  # dpi = 300
)


## error metrics ####

metrics_table <- wf_fig_df %>%
  filter(!model %in% excluded_models0) %>%
  group_by(oos, model) %>%
  summarise(
    RMSE = ModelMetrics::rmse(actual = norm_potential, predicted = fit),
    MAE = ModelMetrics::mae(actual = norm_potential, predicted = fit),
    # MDAPE = mdape(
    #   actual = norm_potential,
    #   predicted = fit,
    #   pos_only = TRUE
    # ),
    Bias = mean(fit - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    model = mod_labels[model],
    oos = ifelse(oos, "OOS", "IS")
  ) %>%
  pivot_wider(
    names_from = oos,
    values_from = c(
      RMSE,
      MAE,
      # MDAPE,
      Bias
    )
  ) %>%
  arrange(desc(RMSE_OOS))
metrics_table
metrics_table %>%
  mutate(
    across(
      c(RMSE_IS, RMSE_OOS, MAE_IS, MAE_OOS, Bias_IS, Bias_OOS),
      ~ round(., 3)
    ),
    # across(c(MDAPE_IS, MDAPE_OOS), ~ round(., 1))
  ) %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    align = "lcccccccc",
    col.names = c(
      "Model",
      "IS",
      "OOS",
      "IS",
      "OOS",
      # "IS",
      # "OOS",
      "IS",
      "OOS"
    ),
    caption = "Performance metrics for in-sample (IS) and out-of-sample (OOS) predictions."
  ) %>%
  add_header_above(c(
    " " = 1,
    "RMSE" = 2,
    "MAE" = 2,
    # "MDAPE (%)" = 2,
    "Bias" = 2
  )) %>%
  kable_styling(latex_options = "hold_position")


# reliability diagrams ####

## Time oos ####
cov_gbl <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    d0_tag <- base::format(d0, "%y%m%d")
    # browser()

    file_name <- sprintf(
      "%s/%s/summaries/oos/pred_band_coverage_summary_time_%s.rds",
      output_path,
      batch_name,
      d0_tag
    )
    if (!file.exists(file_name)) {
      cat("File not found:", file_name, "\n")
      return(NULL)
    }
    cov_obj <- readRDS(file_name)
    cov_gbl <- lapply(
      seq_along(cov_obj),
      \(x) {
        cov_obj[[x]]$cov_gbl %>%
          mutate(
            model = cov_obj %>% names() %>% .[x],
            date = d0
          )
      }
    )
  }
) %>%
  bind_rows()
# list.files("summaries/pred_band_coverage_summary_*.rds") %>%
#   length()
# model_palette <- c(
#   "Observed" = "darkred",
#   "Generic PC" = "#E69F00",
#   "LM bru model" = "#56B4E9",
#   "LM t model" = "#F0E442",
#   "LM+hour model" = "#0072B2",
#   "AR1 model" = "#009E73",
#   "LM beta model" = "#D55E00",
#   "ST model coarser" = "#CC79A7",
#   "QM" = "#999999",
#   "GB LM" = "#000000"
# )
# model_catalog <- read.csv("data/model_catalog.csv") %>%
#   na.omit()
# mod_labels <- model_catalog$mod_labels
# est_cols <- model_catalog$est_cols
# n_models <- length(est_cols)
# names(mod_labels) <- est_cols
# diabrams for time
rel_df <- cov_gbl %>%
  pivot_longer(
    cols = matches("coverage"),
    names_to = "level",
    values_to = "empirical"
  ) %>%
  mutate(
    nominal = as.numeric(gsub("coverage_", "", level)) / 100,
    model = factor(model, levels = names(mod_labels), labels = mod_labels)
  ) %>%
  group_by(model, nominal) %>%
  summarise(empirical = mean(empirical, na.rm = TRUE), .groups = "drop")
# rel_df$model %>% unique()

nrel_mods <- rel_df$model %>% unique() %>% length()
cols <- scales::hue_pal()(nrel_mods)
rel_df %>%
  ggplot(aes(x = nominal, y = empirical, col = model)) +
  geom_line() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkred"
  ) +
  scale_color_manual(values = cols) +
  # scale_color_aaas() +
  guides(colour = guide_legend(nrow = 2)) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
ggsave(
  filename = sprintf("fig/%s/pred_band_reliability_diagram.pdf", batch_name),
  width = 6,
  height = 4,
  # dpi = 300
)

## Space oos ####
cov_gbl <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    d0_tag <- base::format(d0, "%y%m%d")
    file_name <- sprintf(
      "%s/%s/summaries/oos/pred_band_coverage_summary_spaceoos_%s.rds",
      output_path,
      batch_name,
      d0_tag
    )
    if (!file.exists(file_name)) {
      cat("File not found:", file_name, "\n")
      return(NULL)
    }
    cov_obj <- readRDS(file_name)
    cov_gbl <- lapply(
      seq_along(cov_obj),
      \(x) {
        cov_obj[[x]]$cov_gbl %>%
          mutate(
            model = cov_obj %>% names() %>% .[x],
            date = d0
          )
      }
    )
  }
) %>%
  bind_rows()

rel_df <- cov_gbl %>%
  pivot_longer(
    cols = matches("coverage"),
    names_to = "level",
    values_to = "empirical"
  ) %>%
  mutate(
    nominal = as.numeric(gsub("coverage_", "", level)) / 100,
    model = factor(model, levels = names(mod_labels), labels = mod_labels)
  ) %>%
  group_by(model, nominal) %>%
  summarise(empirical = mean(empirical, na.rm = TRUE), .groups = "drop")
# rel_df$model %>% unique()
rel_df %>%
  ggplot(aes(x = nominal, y = empirical, col = model)) +
  geom_line() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkred"
  ) +
  scale_color_manual(values = cols) +
  # scale_color_aaas() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
ggsave(
  filename = sprintf(
    "fig/%s/pred_band_reliability_diagram_spaceoos.pdf",
    batch_name
  ),
  width = 6,
  height = 4,
  # dpi = 300
)

# named listtest

# PIT diagrams ####
pit_df <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    d0_tag <- base::format(d0, "%y%m%d")
    file_name <- sprintf(
      "%s/%s/summaries/fit/model_pit_%s.csv.gz",
      output_path,
      batch_name,
      d0_tag
    )
    if (!file.exists(file_name)) {
      cat("File not found:", file_name, "\n")
      return(NULL)
    }
    fread(file_name) %>%
      mutate(date = d0)
  }
) %>%
  bind_rows() %>%
  filter(!code %in% c("ts_bru0_ar2")) %>%
  mutate(code = paste0(code, "_")) %>%
  left_join(
    model_df %>% dplyr::select(mode_code_prefix, code),
    by = c("code" = "mode_code_prefix")
  )
# mutate(
#   code = model_df$code[model_df$mode_code_prefix == paste0(code, "_")],
#   model = factor(code, levels = est_cols, labels = mod_labels)
# )

p <- pit_df %>%
  slice_sample(n = 1e6) %>%
  mutate(model = factor(code.y, levels = est_cols, labels = mod_labels)) %>%
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0), col = "darkgray") +
  stat_ecdf(aes(pit, col = model), na.rm = TRUE, lwd = 0.8) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(.75, .25),
    plot.title = element_text(size = 10), # Adjust title font size
    axis.text = element_text(size = 8), # Adjust axis text font size
    axis.title = element_text(size = 9), # Adjust axis label font size
    legend.text = element_text(size = 8), # Adjust legend text font size
    legend.title = element_text(size = 8), # Adjust legend title font size
    legend.background = element_blank(), # Makes background completely transparent
    legend.box.background = element_rect(fill = NA, color = NA) # No border
  ) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "PIT", y = "ECDF", col = "") +
  scale_color_manual(values = cols)
p
ggsave(
  filename = sprintf("fig/%s/pit_diagram.png", batch_name),
  width = 6,
  height = 4,
  dpi = 100
)
# Scores table ####
# CRPS LogScore Energy Brie
## time ####

scores_tbl_t <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    d0_tag <- base::format(d0, "%y%m%d")
    file_name <- sprintf(
      "%s/%s/summaries/oos/model_scores_summary_time_%s.csv",
      output_path,
      batch_name,
      d0_tag
    )
    if (!file.exists(file_name)) {
      cat("File not found:", file_name, "\n")
      return(NULL)
    }
    read.csv(file_name) %>%
      mutate(
        date = d0,
        # prefix = sub("\\d{6}\\.rds$", "", model),
        # model = factor(
        #   prefix,
        #   levels = model_catalog$mode_code_prefix,
        #   labels = model_catalog$est_cols
        # )
        model = bru_df$label[which(bru_df$code == model)]
      ) %>%
      filter(!is.na(model))
  }
) %>%
  bind_rows() %>%
  group_by(model) %>%
  summarise(
    across(c(crps, energy, log, matches("bs")), ~ mean(., na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(crps) %>%
  mutate(
    model = factor(model, levels = model)
  )
# scores_tbl_t %>%
#   dplyr::select(-bs_0_2) %>%
#   kbl(
#     format = "latex",
#     booktabs = TRUE,
#     digits = 3,
#     align = "lcccccc",
#     col.names = c(
#       "Model",
#       "CRPS",
#       "Energy",
#       "Log",
#       "BS (CF < 1%)",
#       "BS (CF < 5%)",
#       "BS (CF < 10%)"
#     ),
#     caption = "Average scores for out-of-sample predictions across all sampled days."
#   ) %>%
#   kable_styling(latex_options = "hold_position")

scores_tbl_t %>%
  dplyr::select(-bs_0_2) %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    digits = 3,
    align = "lcccccc",
    col.names = c(
      "Model",
      "CRPS",
      "Energy",
      "Log",
      "1\\%",
      "5\\%",
      "10\\%"
    ),
    caption = "Average scores for out-of-sample predictions across all sampled days."
  ) %>%
  add_header_above(c(
    " " = 4,
    "Brier score (CF threshold)" = 3
  )) %>%
  kable_styling(latex_options = "hold_position")

## space ####
scores_tbl <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    d0_tag <- base::format(d0, "%y%m%d")
    file_name <- sprintf(
      "%s/%s/summaries/oos/model_scores_summary_spaceoos_%s.csv",
      output_path,
      batch_name,
      d0_tag
    )
    if (!file.exists(file_name)) {
      cat("File not found:", file_name, "\n")
      return(NULL)
    }
    read.csv(file_name) %>%
      mutate(date = d0, model = bru_df$label[which(bru_df$code == model)])
  }
) %>%
  bind_rows() %>%
  group_by(model) %>%
  summarise(
    across(c(crps, energy, log, matches("bs")), ~ mean(., na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(crps) %>%
  mutate(
    model = factor(model, levels = model)
  )
scores_tbl %>%
  dplyr::select(-bs_0_2) %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    digits = 3,
    align = "lcccccc",
    col.names = c(
      "Model",
      "CRPS",
      "Energy",
      "Log",
      "1\\%",
      "5\\%",
      "10\\%"
    ),
    caption = "Average scores for out-of-sample predictions across all sampled days."
  ) %>%
  add_header_above(c(
    " " = 4,
    "Brier score (CF threshold)" = 3
  )) %>%
  kable_styling(latex_options = "hold_position")

# rm(gb_fig_df, wf_fig_df, pit_df)
# gc()
