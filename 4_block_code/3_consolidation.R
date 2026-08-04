local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

mc <- ifelse(local_run, 1, available_cores())
pow_threshold <- 0.05
pow_threshold_label <- as.character(pow_threshold) %>% gsub("\\.", "_", .)


override_objects <- FALSE
# rerun_samples <- FALSE
# prec_init <- log(200)
batch_name <- "batch2025"


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
  n_samp <- 100
  pixel_dims <- c(150, 150)
} else {
  data_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  gen_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  model_path <- "/exports/eddie/scratch/s2441782/calibration/model_objects"
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

# sampled_days_df <- read.csv("data/sample_days_df.csv") %>%
#   mutate(date = as.Date(date))

# sampled_days_df <- read.csv("data/sample_days_df.csv") %>%
#   mutate(date = as.Date(date))

sampled_days <- sampled_days_df %>%
  pull(date) %>%
  sort()


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
excluded_models <- c("qm", "lm_bru")

# Time oos ####
## read summary tables of prediction bands ######
gb_fig_df <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    # print(i)
    d0_tag <- base::format(d0, "%y%m%d")
    readRDS(sprintf("summaries/GB_fig_band_summary_%s.rds", d0_tag)) %>%
      mutate(
        pgroup3 = sampled_days_df$p_group3[i] %>%
          factor(levels = c("low", "mid", "high"))
      )
  }
) %>%
  bind_rows()

wf_fig_df <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    d0_tag <- base::format(d0, "%y%m%d")
    readRDS(sprintf("summaries/WF_fig_band_summary_%s.rds", d0_tag)) %>%
      mutate(
        pgroup3 = sampled_days_df$p_group3[i] %>%
          factor(levels = c("low", "mid", "high"))
      )
  }
) %>%
  bind_rows()

## calibration fit scatter####
gb_fig_df %>%
  filter(oos) %>%
  ggplot(aes(x = norm_potential, y = mean, col = pgroup3)) +
  geom_point() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkred"
  ) +
  facet_wrap(~model, labeller = as_labeller(mod_labels)) +
  theme_bw() +
  # scale_color_lancet() +
  scale_color_manual(values = regime_palette) +
  labs(
    x = "Observed power",
    y = "Predicted power",
    col = "regime"
  )

ggsave(
  sprintf("fig/%s/GB_fit_scatter_oos.pdf", batch_name),
  width = 10,
  height = 6,
  dpi = 300
)

wf_fig_df %>%
  mutate(hour = hour(time)) %>%
  filter(oos, hour <= 1) %>%
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
  # scale_fill_viridis_c(
  #   trans = "log10",
  #   name = "frequency",
  #   limits = c(1, NA)
  # ) +
  scale_fill_gradient(
    trans = "log10",
    low = "grey90",
    high = blues9[5],
    # limits = c(5, NA),
  ) +
  theme_bw() +
  labs(
    x = "Normalised potential power",
    y = "Normalised predicted power",
    fill = "Count"
  )
# scales::show_col(blues9)
ggsave(
  sprintf("fig/%s/WF_fit_scatter_oos.pdf", batch_name),
  width = 10,
  height = 6,
  dpi = 300
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
metrics_table <- wf_fig_df %>%
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
  )
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
    # format = "latex",
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
gb_fig_df <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    # print(i)
    d0_tag <- base::format(d0, "%y%m%d")
    readRDS(sprintf(
      "summaries/GB_fig_band_summary_spaceoos_%s.rds",
      d0_tag
    )) %>%
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

wf_fig_df <- lapply(
  seq_along(sampled_days[-15]),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    d0_tag <- base::format(d0, "%y%m%d")
    readRDS(sprintf(
      "summaries/WF_fig_band_summary_spaceoos_%s.rds",
      d0_tag
    )) %>%
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
  ggplot(aes(x = norm_potential, y = mean, col = pgroup3)) +
  geom_point() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkred"
  ) +
  facet_wrap(~model, labeller = as_labeller(mod_labels)) +
  theme_bw() +
  # scale_color_lancet() +
  scale_color_manual(values = regime_palette) +
  labs(
    x = "Observed power",
    y = "Predicted power",
    col = "regime"
  )

ggsave(
  sprintf("fig/%s/GB_fit_scatter_oos_space.pdf", batch_name),
  width = 10,
  height = 6,
  dpi = 300
)

wf_fig_df %>%
  filter(oos) %>%
  ggplot(aes(x = norm_potential, y = fit)) +
  geom_hex() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkred"
  ) +
  facet_wrap(~model, labeller = as_labeller(mod_labels)) +
  # scale_fill_viridis_c(
  #   trans = "log10",
  #   name = "frequency",
  #   limits = c(1, NA)
  # ) +
  scale_fill_gradient(
    trans = "log10",
    low = "grey90",
    high = blues9[8]
  ) +
  theme_bw() +
  labs(
    x = "Normalised potential power",
    y = "Normalised predicted power",
    fill = "Count"
  )

ggsave(
  sprintf("fig/%s/WF_fit_scatter_oos_space.pdf", batch_name),
  width = 10,
  height = 6,
  dpi = 300
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
  )
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
    # format = "latex",
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
    cov_obj <- readRDS(sprintf(
      "summaries/pred_band_coverage_summary_%s.rds",
      d0_tag
    ))
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
model_palette <- c(
  "Observed" = "darkred",
  "Generic PC" = "#E69F00",
  "LM bru model" = "#56B4E9",
  "LM t model" = "#F0E442",
  "LM+hour model" = "#0072B2",
  "AR1 model" = "#009E73",
  "LM beta model" = "#D55E00",
  "ST model coarser" = "#CC79A7",
  "QM" = "#999999",
  "GB LM" = "#000000"
)
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
rel_df %>%
  ggplot(aes(x = nominal, y = empirical, col = model)) +
  geom_line() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "darkred"
  ) +
  # scale_color_manual(values = model_palette) +
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

    cov_obj <- readRDS(sprintf(
      "summaries/pred_band_coverage_summary_spaceoos_%s.rds",
      d0_tag
    ))
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
  # scale_color_manual(values = model_palette) +
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
a <- list(
  "GB" = list(tbl = gb_fig_df, extra = 1),
  "WF" = list(tbl = wf_fig_df, extra = 2)
)
test <- lapply(
  a,
  function(element) {
    element$tbl %>% head()
  }
) %>%
  bind_rows(.id = "source")


# PIT diagrams ####

pit_df <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    d0_tag <- base::format(d0, "%y%m%d")
    read.csv(sprintf("summaries/model_pit_%s.csv", d0_tag)) %>%
      mutate(date = d0)
  }
) %>%
  bind_rows() %>%
  mutate(
    model = factor(code, levels = est_cols, labels = mod_labels)
  )


# Scores table ####
# CRPS LogScore Energy Brie
## time ####
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
scores_tbl <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    d0_tag <- base::format(d0, "%y%m%d")

    read.csv(sprintf("summaries/model_scores_summary_%s.csv", d0_tag)) %>%
      mutate(date = d0, model = bru_df$label)
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
## space ####
scores_tbl <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days_df$date[i] %>% as.Date()
    d0_tag <- base::format(d0, "%y%m%d")
    read.csv(sprintf(
      "summaries/model_scores_summary_spaceoos_%s.csv",
      d0_tag
    )) %>%
      mutate(date = d0, model = bru_df$label)
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

# LWE diagram ####
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
extension <- "rds"
model_df0 <- lapply(
  seq_along(sampled_days),
  function(i) {
    d0 <- sampled_days[i] %>% as.Date()
    # print(d0)
    d0_tag <- base::format(d0, "%y%m%d")

    # coarse
    readRDS(sprintf(
      "data/calibration_df_%s_%s.%s",
      # mesh_label,
      "very_coarse",
      d0_tag,
      extension
    )) %>%
      st_drop_geometry()
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
    # model = factor(model, levels = est_cols, labels = mod_labels)
  )

## 2.1 Low wind events in observed data ####
cat("--------------------------------------------------------------------\n")
cat("Low wind events in observed data\n")
cat("--------------------------------------------------------------------\n")

lwe_obs_pred_fname <- file.path(
  "summaries",
  sprintf(
    "lwe_obs_pred_%s_t%s.%s",
    "sampled_days",
    pow_threshold_label,
    extension
  )
)
lwe_obs_newloc_fname <- file.path(
  "summaries",
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

  days_before <- 0
  days_after <- 3
  pred_days <- lapply(
    sampled_days %>% as.Date(),
    \(d0) {
      seq(d0 - days_before, d0 + days_after - 1, by = "days") %>%
        as.Date()
    }
  ) %>%
    unlist() %>%
    unique() %>%
    sort() %>%
    as.Date()

  lwe_obs_pred <- pwr_curv_df %>%
    rename(time = halfHourEndTime) %>%
    mutate(
      date = as.Date(time)
    ) %>%
    # filter(date >= d0, date <= d0 + n.days) %>%
    filter(date %in% pred_days) %>%
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

## 3. Low wind events in model predictions ####
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
  bind_rows(lwe_obs_pred) %>%
  mutate(
    model = factor(model, levels = est_cols, labels = mod_labels)
  )

# low_events_model$model %>% unique() %>% sort() %>% print()
low_events_model %>%
  filter(duration_h < 24 * 7) %>%
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
  facet_wrap(
    ~model,
    scales = "free",
    # labeller = as_labeller(mod_labels)
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_npg()

ggsave(
  sprintf(
    "fig/%s/low_wind_duration_dist_t%s.pdf",
    batch_name,
    pow_threshold_label
  ),
  width = 10,
  height = 6
)

low_events_model %>%
  filter(duration_h < 24 * 7) %>%
  ggplot(aes(x = duration_h)) +
  geom_density(
    # data = ~ dplyr::filter(.x, model != "observed"),
    aes(colour = model),
    alpha = 0.5
  ) +
  # geom_density(
  #   data = ~ dplyr::filter(.x, model == "observed"),
  #   aes(colour = model),
  #   # colour = "darkred",
  #   # linewidth = 1.2,
  #   fill = NA
  # ) +
  labs(
    title = "Distribution of Low Wind Events Duration",
    x = "Duration (hours)",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "inside",
    legend.position.inside = c(0.7, 0.8)
  ) +
  guides(
    colour = guide_legend(ncol = 2)
  ) +
  scale_color_manual(
    values = c(
      "Observed" = "darkred",
      "Generic PC" = "#E69F00",
      "Linear model" = "#56B4E9",
      "AR1 model" = "#009E73",
      "AR2 model" = "#F0E442",
      "LM+hour model" = "#0072B2",
      "ST model fine" = "#D55E00",
      "ST model coarse" = "#CC79A7",
      "QM" = "#999999",
      "GB LM" = "#000000"
    )
  )

ggsave(
  sprintf(
    "fig/%s/low_wind_duration_dens_t%s.pdf",
    batch_name,
    pow_threshold_label
  ),
  width = 6,
  height = 4
)

obs <- low_events_model %>%
  filter(
    model == "observed",
    duration_h < 24 * 7
  ) %>%
  pull(duration_h)

probs <- seq(0, 1, length.out = 101)

qq_df <- low_events_model %>%
  filter(model != "observed") %>%
  filter(duration_h < 24 * 7) %>%
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
  facet_wrap(
    ~model,
    scales = "free",
    #  labeller = as_labeller(mod_labels)
  ) +
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

ggplot(qq_df, aes(x = obs_q, y = model_q)) +
  geom_line(aes(col = model)) +
  geom_abline(slope = 1, intercept = 0, colour = "darkred") +
  # facet_wrap(
  #   ~model,
  #   scales = "free",
  #   #  labeller = as_labeller(mod_labels)
  # ) +
  labs(
    title = "Q-Q Plot of Low Wind Event Durations",
    x = "Observed quantiles",
    y = "Model quantiles"
  ) +
  # coord_equal() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_manual(
    values = c(
      "Generic PC" = "#E69F00",
      "Linear model" = "#56B4E9",
      "AR1 model" = "#009E73",
      "AR2 model" = "#F0E442",
      "LM+hour model" = "#0072B2",
      "ST model fine" = "#D55E00",
      "ST model coarse" = "#CC79A7",
      "QM" = "#999999",
      "GB LM" = "#000000"
    )
  ) +
  coord_equal(ratio = 1)

ggsave(
  sprintf(
    "fig/%s/low_wind_duration_qq_line_t%s.pdf",
    batch_name,
    pow_threshold_label
  ),
  width = 6,
  height = 5
)
