# 1. Load libraries and data ------------------------------------------------

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE


require(parallel)

if (local_run) {
  data_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
} else {
  data_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  gen_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  model_path <- "/exports/eddie/scratch/s2441782/calibration/model_objects"
  temp_lib <- "/exports/eddie3_homes_local/s2441782/lib"
  .libPaths(temp_lib)
}


# install.packages(
#   c("qmap"),
#   temp_lib,
#   dependencies = TRUE
# )

require(arrow)
require(dplyr)
require(tidyr)
# require(rnaturalearth)
# require(rnaturalearthdata)
require(sf)
require(ggplot2)
require(ggthemes)
require(ggsci)
# require(FNN)
require(data.table)
# require(parallel)
library(purrr)
# require(brms)
require(INLA)
require(inlabru)
require(fmesher)
require(qmap)
require(ggridges)
require(ggspatial)
require(elevatr)

source("aux_funct.R")


# 1.1 Load data ----------------------------------------------------------------
source("read_data.R")

pwr_curv_df <- read_parquet(file.path(gen_path, "power_curve_all.parquet")) %>%
  mutate(
    site_name = site_name %>%
      gsub("\\b(wind\\s*farm|wf)\\b", "", ., ignore.case = TRUE) %>%
      trimws()
  )


# 1.2 choose days to run ----------------------------------------------------------------
## by regime

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

gb_day_df %>%
  ggplot() +
  geom_histogram(aes(norm_potential), bins = 30) +
  geom_vline(xintercept = p_quant3, col = "red") +
  geom_vline(xintercept = p_quant7, col = "blue") +
  theme_minimal() +
  labs(
    title = "Distribution of daily generation",
    x = "Wind generation (% of capacity)",
    y = "Count"
  ) +
  theme(legend.position = "bottom")

gb_day_df %>%
  ggplot() +
  geom_line(aes(date, norm_potential), col = "blue") +
  theme_minimal() +
  labs(
    title = "Daily generation over time",
    x = "Date",
    y = "Wind generation (% of capacity)"
  ) +
  theme(legend.position = "bottom")


regime_palette <- pal_lancet()(3)[c(2, 3, 1)]
gb_day_df %>%
  ggplot(aes(norm_potential, p_group3, fill = p_group3)) +
  geom_density_ridges(alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Daily generation by regime",
    x = "Regime",
    y = "Wind generation (% of capacity)"
  ) +
  theme(legend.position = "bottom") +
  labs(fill = "Regime") +
  scale_fill_manual(values = regime_palette)

set.seed(1)
# sample 1 day per regime
sampled_days <- gb_day_df %>%
  group_by(p_group3) %>%
  slice_sample(n = 1) %>%
  pull(date)


sample_df <- GB_df %>%
  filter(date %in% sampled_days) %>%
  group_by(time) %>%
  summarise(
    across(
      c(norm_power_est0, norm_potential),
      ~ sum(. * capacity) / sum(capacity)
    ),
    across(c(ws_h_wmean), ~ sum(. * capacity) / sum(capacity)),
    across(c(capacity), sum),
    date = first(date),
    .groups = "drop"
  ) %>%
  left_join(gb_day_df %>% dplyr::select(date, p_group3), by = "date") %>%
  mutate(
    hour = format(time, "%H:%M"),
    legend_label = factor(
      paste(p_group3, format(date, "%y-%m-%d")),
      levels = paste(
        c("low", "mid", "high"),
        format(sampled_days, "%y-%m-%d")
      )
    )
  )

ggplot() +
  geom_line(
    data = sample_df,
    aes(
      hour,
      norm_potential,
      col = legend_label,
      group = legend_label
    )
  ) +
  labs(
    title = "Daily generation time series for sampled days",
    x = "Hour",
    y = "Wind generation (% of capacity)",
    col = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 1)
  ) +
  scale_color_manual(values = regime_palette)

#

# 2.1 Model for WF ####

## 2.1.1 Model predictions ####

# 2.2 Model for aggregate ####

## 2.2.1 Model predictions ####

# 3.1 Comparison

#how performance changes by onshore/offshore, season, hour, and regime and how well farm level predictions aggregate

# aggregator GB, Onshore, Offshore, Scotland, rest of GB

# validation
# time GB Onshore, Offshore,
# dense regions, sparse regions,
# low, mid, high generation regimes
# Scotland, rest of GB
