# ============================================================
# National wind generation: distribution and lower-tail analysis
# Required columns: bmUnit, time, potential, capacity
# ============================================================

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

pow_threshold <- 0.05
pow_threshold_label <- gsub("\\.", "_", as.character(pow_threshold))
tables_path <- "~/ownCloud-s2441782@pwr_curv_dfasync.ed.ac.uk/projects/calibration/calibration_power_main_doc/tables"
dir.exists(tables_path) || dir.create(tables_path, recursive = TRUE)

override_objects <- FALSE
batch_name <- "bias_eda_test"
norm_dist_tol <- 1

if (local_run) {
  cat("Running in local mode\n")
} else {
  cat("Running in cluster mode\n")
}

# 0.2 libraries and paths ####
require(parallel)

if (local_run) {
  pwr_curv_dfa_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
  extra_path <- "~/Documents/elexon/extra"
  output_path <- "~/Documents/elexon/caloutput"
  sample_path <- "~/Documents/elexon/samples"
  n_samp <- 100
  pixel_dims <- c(150, 150)
} else {
  pwr_curv_dfa_path <- "/exports/eddie/scratch/s2441782/calibration/pwr_curv_dfa"
  gen_path <- "/exports/eddie/scratch/s2441782/calibration/pwr_curv_dfa"
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

lwe_output_figpath <- file.path("fig", "lwe_eda")
dir.exists(lwe_output_figpath) ||
  dir.create(lwe_output_figpath, recursive = TRUE)
lwe_output_tabpath <- file.path("tables", "lwe_eda")
dir.exists(lwe_output_tabpath) ||
  dir.create(lwe_output_tabpath, recursive = TRUE)

# pwr_curv_dfa

extension <- "rds"
coord_list_fname <- "pwr_curv_dfa/coord_list.csv"

cat("Loading existing coordinate list\n")
coord_list <- read.csv(coord_list_fname)

ref_catalog_2025 <- fread(
  file.path("pwr_curv_dfa/ref_catalog_wind_2025.csv.gz")
)

pwr_curv_df <- read_parquet(file.path(
  gen_path,
  "power_curve_all_enriched.parquet"
)) %>%
  mutate(
    settlementDate = as.Date(settlementDate),

    time_utc_reconstructed = as.POSIXct(
      settlementDate,
      tz = "UTC"
    ) +
      (settlementPeriod - 2) * 30 * 60,

    time_utc = if_else(
      settlementPeriod == 50,
      as.POSIXct(halfHourEndTime, tz = "UTC"),
      time_utc_reconstructed
    ),

    time = with_tz(time_utc, "Europe/London")
  )
# mutate(
#   time_utc = with_tz(halfHourEndTime, "UTC"),
#   time = with_tz(time_utc, "Europe/London")
# )
# mutate(
#   settlementDate = as.Date(settlementDate),
#   time_utc = as.POSIXct(
#     settlementDate,
#     tz = "UTC"
#   ) +
#     (settlementPeriod - 2) * 30 * 60,
#   time = with_tz(time_utc, "Europe/London")
# )
head(pwr_curv_df)

# ------------------------------------------------------------
# 0. Settings
# ------------------------------------------------------------

# Lowest 1% of national capacity factor
q_low <- 0.01

# Set to FALSE only if you deliberately want to analyse
# hours with incomplete fleet coverage.
require_complete_hours <- FALSE

# ------------------------------------------------------------
# 1. Basic pwr_curv_dfa checks
# ------------------------------------------------------------

required_cols <- c("bmUnit", "time", "potential", "capacity")

# if (anyNA(pwr_curv_df$bmUnit) || anyNA(pwr_curv_df$time)) {
#   stop("Missing bmUnit or time values. Resolve these before aggregation.")
# }

# if (any(pwr_curv_df$potential < 0, na.rm = TRUE)) {
#   warning("Negative potential generation values found.")
# }

# if (any(pwr_curv_df$capacity <= 0, na.rm = TRUE)) {
#   warning("Non-positive capacity values found.")
# }

# Check that each farm has at most one record per hour.
pwr_curv_df %>%
  filter(
    bmUnit == "2__PEDGE003" &
      # time <= as.POSIXct("2024-03-31 02:00:00") + hours(1) #&
      time >= as.POSIXct("2024-10-28", tz = "UTC") + hours(-25)
  ) %>%
  View()
pwr_curv_df %>%
  filter(
    bmUnit == "2__PEDGE003" &
      time >= as.POSIXct("2024-03-31 03:00:00", tz = "UTC") + hours(-2)
  )
duplicates <- pwr_curv_df %>%
  count(bmUnit, time, name = "n") %>%
  filter(n > 1)

if (nrow(duplicates) > 0) {
  print(head(duplicates, 20))
  stop(
    "Duplicate farm-hour records found. Resolve these before ",
    "aggregating; do not automatically sum duplicates."
  )
}

# Check whether installed capacity is constant within each farm.
capacity_check <- pwr_curv_df %>%
  group_by(bmUnit) %>%
  summarise(
    n_capacity_values = n_distinct(capacity, na.rm = TRUE),
    min_capacity = min(capacity, na.rm = TRUE),
    max_capacity = max(capacity, na.rm = TRUE),
    .groups = "drop"
  )

print(
  capacity_check %>%
    filter(n_capacity_values > 1)
)

# If capacity changes over time, use the time-varying capacity
# in the aggregation below. Investigate unexpected changes.

# ------------------------------------------------------------
# 2. Aggregate to national hourly figures
# ------------------------------------------------------------

n_farms_total <- n_distinct(pwr_curv_df$bmUnit)

national <- pwr_curv_df %>%
  group_by(time) %>%
  summarise(
    n_farms = n_distinct(bmUnit),
    n_valid = sum(
      !is.na(potential) &
        !is.na(capacity) &
        capacity > 0
    ),
    potential_MW = sum(
      potential[
        !is.na(potential) &
          !is.na(capacity) &
          capacity > 0
      ]
    ),
    capacity_MW = sum(
      capacity[
        !is.na(potential) &
          !is.na(capacity) &
          capacity > 0
      ]
    ),
    .groups = "drop"
  ) %>%
  mutate(
    coverage_fraction = n_valid / n_farms_total,
    capacity_factor = if_else(
      capacity_MW > 0,
      potential_MW / capacity_MW,
      NA_real_
    )
  )

# Inspect hourly coverage before proceeding.
print(summary(national$coverage_fraction))
print(summary(national$capacity_factor))

# Retain hours with all farms represented and valid.
# This is conservative: it can discard substantial pwr_curv_dfa if the
# fleet changes over the seven-year period.
if (require_complete_hours) {
  national <- national %>%
    filter(n_valid == n_farms_total)
}

# Check capacity factor bounds. Values above 1 may indicate
# inconsistent units, capacity definitions or pwr_curv_dfa errors.
print(
  national %>%
    summarise(
      min_CF = min(capacity_factor, na.rm = TRUE),
      max_CF = max(capacity_factor, na.rm = TRUE),
      n_hours = n(),
      n_missing_CF = sum(is.na(capacity_factor))
    )
)

if (any(national$capacity_factor > 1, na.rm = TRUE)) {
  warning(
    "Capacity factor exceeds 1. Check units and capacity definitions."
  )
}

# Retain valid observations for distributional analysis.
national <- national %>%
  filter(
    is.finite(capacity_factor),
    capacity_factor >= 0
  ) %>%
  mutate(
    month = lubridate::month(time, label = TRUE, abbr = TRUE),
    month_num = month(time),
    year = year(time),
    season = case_when(
      month_num %in% c(12, 1, 2) ~ "Winter",
      month_num %in% 3:5 ~ "Spring",
      month_num %in% 6:8 ~ "Summer",
      month_num %in% 9:11 ~ "Autumn"
    ),
    season = factor(
      season,
      levels = c("Winter", "Spring", "Summer", "Autumn")
    )
  )

# ------------------------------------------------------------
# 3. National distribution and lower quantiles
# ------------------------------------------------------------

national_quantiles <- quantile(
  national$capacity_factor,
  probs = c(0.01, 0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90),
  na.rm = TRUE,
  names = TRUE
)

print(national_quantiles)

# The 1% threshold in capacity-factor units
q01 <- unname(
  quantile(national$capacity_factor, q_low, na.rm = TRUE)
)
q05 <- unname(
  quantile(national$capacity_factor, 0.05, na.rm = TRUE)
)

cat("National 1% capacity-factor threshold:", q01, "\n")
cat("National 5% capacity-factor threshold:", q05, "\n")
cat("Equivalent percentage of installed capacity:", 100 * q01, "%\n")

# Distribution: histogram with 1% threshold
p_distribution <- ggplot(
  national,
  aes(x = capacity_factor)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 60,
    fill = blues9[6],
    colour = NA
  ) +
  geom_density(linewidth = 0.8, na.rm = TRUE) +
  geom_vline(
    xintercept = q01,
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = q05,
    linetype = "dashed",
    linewidth = 1,
    colour = "red"
  ) +
  labs(
    title = "Distribution of national wind capacity factor",
    subtitle = sprintf(
      "Dashed line: empirical 1st percentile = %.3f\nDashed line (red): empirical 5th percentile = %.3f",
      q01,
      q05
    ),
    x = "National capacity factor",
    y = "Density"
  ) +
  theme_minimal()

print(p_distribution)
ggsave(
  p_distribution,
  filename = file.path(lwe_output_figpath, "national_CF_distribution.png"),
  width = 8,
  height = 6,
  dpi = 300
)
# Optional: empirical cumulative distribution function (ECDF)
p_ecdf <- ggplot(
  national,
  aes(x = capacity_factor)
) +
  stat_ecdf(linewidth = 0.8) +
  geom_hline(yintercept = q_low, linetype = "dashed") +
  geom_vline(xintercept = q01, linetype = "dashed") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Empirical distribution of national capacity factor",
    x = "National capacity factor",
    y = "Cumulative probability"
  ) +
  theme_minimal()

print(p_ecdf)
ggsave(
  p_ecdf,
  filename = file.path(lwe_output_figpath, "national_CF_ecdf.png"),
  width = 8,
  height = 6,
  dpi = 300
)
# ------------------------------------------------------------
# 4. Monthly distributions
# ------------------------------------------------------------

# Boxplots show median, interquartile range and outliers.
p_monthly_box <- ggplot(
  national,
  aes(x = month, y = capacity_factor)
) +
  geom_boxplot(
    outlier.alpha = 0.15,
    na.rm = TRUE
  ) +
  labs(
    title = "Monthly distribution of national wind capacity factor",
    x = "Calendar month",
    y = "National capacity factor"
  ) +
  theme_minimal()

print(p_monthly_box)

ggsave(
  p_monthly_box,
  filename = file.path(lwe_output_figpath, "national_CF_monthly_boxplot.png"),
  width = 8,
  height = 6,
  dpi = 300
)
# Monthly density distributions, faceted for comparison
p_monthly_density <- ggplot(
  national,
  aes(x = capacity_factor)
) +
  geom_density(na.rm = TRUE, fill = blues9[6]) +
  facet_wrap(~month, ncol = 4) +
  labs(
    title = "National capacity-factor distributions by month",
    x = "National capacity factor",
    y = "Density"
  ) +
  theme_minimal()

ggsave(
  p_monthly_density,
  filename = file.path(lwe_output_figpath, "national_CF_monthly_density.png"),
  width = 8,
  height = 6,
  dpi = 300
)

print(p_monthly_density)

# Monthly quantiles: estimated separately for each calendar month
monthly_quantiles <- national %>%
  group_by(month, month_num) %>%
  summarise(
    n_hours = n(),
    q01 = quantile(capacity_factor, 0.01, na.rm = TRUE),
    q05 = quantile(capacity_factor, 0.05, na.rm = TRUE),
    q10 = quantile(capacity_factor, 0.10, na.rm = TRUE),
    median = median(capacity_factor, na.rm = TRUE),
    mean = mean(capacity_factor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(month_num)

print(monthly_quantiles)

# Plot monthly lower quantiles
monthly_quantiles_long <- monthly_quantiles %>%
  dplyr::select(month, q01, q05, q10) %>%
  pivot_longer(
    cols = c(q01, q05, q10),
    names_to = "quantile",
    values_to = "capacity_factor"
  )

p_monthly_quantiles <- ggplot(
  monthly_quantiles_long,
  aes(x = month, y = capacity_factor, group = quantile, linetype = quantile)
) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  labs(
    title = "Lower quantiles by calendar month",
    x = "Calendar month",
    y = "National capacity factor",
    linetype = "Quantile"
  ) +
  theme_minimal()

print(p_monthly_quantiles)
ggsave(
  p_monthly_quantiles,
  filename = file.path(lwe_output_figpath, "national_CF_monthly_quantiles.png"),
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 5. Seasonal distributions and seasonal quantiles
# ------------------------------------------------------------

p_seasonal_box <- ggplot(
  national,
  aes(x = season, y = capacity_factor)
) +
  geom_boxplot(outlier.alpha = 0.15) +
  labs(
    title = "Seasonal distribution of national wind capacity factor",
    x = "Season",
    y = "National capacity factor"
  ) +
  theme_minimal()

ggsave(
  p_seasonal_box,
  filename = file.path(lwe_output_figpath, "national_CF_seasonal_boxplot.png"),
  width = 8,
  height = 6,
  dpi = 300
)

print(p_seasonal_box)

seasonal_quantiles <- national %>%
  group_by(season) %>%
  summarise(
    n_hours = n(),
    q01 = quantile(capacity_factor, 0.01, na.rm = TRUE),
    q05 = quantile(capacity_factor, 0.05, na.rm = TRUE),
    q10 = quantile(capacity_factor, 0.10, na.rm = TRUE),
    median = median(capacity_factor, na.rm = TRUE),
    mean = mean(capacity_factor, na.rm = TRUE),
    .groups = "drop"
  )

print(seasonal_quantiles)

# ------------------------------------------------------------
# 6. Annual distributions and annual quantiles
# ------------------------------------------------------------

# Check whether lower-tail behaviour differs across years.
annual_quantiles <- national %>%
  group_by(year) %>%
  summarise(
    n_hours = n(),
    q01 = quantile(capacity_factor, 0.01, na.rm = TRUE),
    q05 = quantile(capacity_factor, 0.05, na.rm = TRUE),
    median = median(capacity_factor, na.rm = TRUE),
    mean = mean(capacity_factor, na.rm = TRUE),
    .groups = "drop"
  )

print(annual_quantiles)

p_annual <- ggplot(
  annual_quantiles,
  aes(x = year, y = q01)
) +
  geom_line(col = blues9[7]) +
  geom_line(aes(y = q05), colour = blues9[4]) +
  # geom_point() +
  labs(
    title = "Annual variation in the 1st and 5th percentiles",
    x = "Year",
    y = "Annual 1st and 5th percentiles of capacity factor"
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 0.1))

print(p_annual)
ggsave(
  p_annual,
  filename = file.path(lwe_output_figpath, "national_CF_annual_quantiles.png"),
  width = 8,
  height = 6,
  dpi = 300
)
# ------------------------------------------------------------
# 6. Optional: month-by-year and season-by-year diagnostics
# ------------------------------------------------------------
annualquarter_quantiles <- national %>%
  mutate(
    quarter = lubridate::quarter(time, type = "date_last")
  ) %>%
  group_by(quarter) %>%
  summarise(
    n_hours = n(),
    q01 = quantile(capacity_factor, 0.01, na.rm = TRUE),
    q05 = quantile(capacity_factor, 0.05, na.rm = TRUE),
    median = median(capacity_factor, na.rm = TRUE),
    mean = mean(capacity_factor, na.rm = TRUE),
    .groups = "drop"
  )

print(annualquarter_quantiles)

p_annualq <- ggplot(
  annualquarter_quantiles,
  aes(x = quarter, y = q01)
) +
  geom_line(col = blues9[7]) +
  # geom_point() +
  geom_line(
    aes(y = q05),
    colour = blues9[4]
  ) +
  labs(
    title = "Quarterly variation in the 1st and 5th percentiles",
    x = "Quarter (year + quarter)",
    y = "Quarterly 1st and 5th percentiles of capacity factor"
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 0.15)) #+
# scale_x_date(
#   date_breaks = "3 months",
#   labels = function(x) {
#     paste0(
#       "Q",
#       lubridate::quarter(x),
#       " ",
#       lubridate::year(x)
#     )
#   }
# ) +
# theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(p_annualq)
ggsave(
  p_annualq,
  filename = file.path(
    lwe_output_figpath,
    "national_CF_quarterly_quantiles.png"
  ),
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 7. Save the main outputs
# ------------------------------------------------------------

write.csv(
  national,
  file.path(lwe_output_tabpath, "national_hourly_generation.csv"),
  row.names = FALSE
)
write.csv(
  monthly_quantiles,
  file.path(lwe_output_tabpath, "national_monthly_quantiles.csv"),
  row.names = FALSE
)
write.csv(
  seasonal_quantiles,
  file.path(lwe_output_tabpath, "national_seasonal_quantiles.csv"),
  row.names = FALSE
)
write.csv(
  annual_quantiles,
  file.path(lwe_output_tabpath, "national_annual_quantiles.csv"),
  row.names = FALSE
)

write.csv(
  annualquarter_quantiles,
  file.path(lwe_output_tabpath, "national_quarterly_quantiles.csv"),
  row.names = FALSE
)
# Objects available for further analysis:
# national
# national_quantiles
# q01
# monthly_quantiles
# seasonal_quantiles
# annual_quantiles
# p_distribution, p_ecdf, p_monthly_box, p_monthly_density,
# p_monthly_quantiles, p_seasonal_box, p_annual

# Wind events #######

# 1. Calculate the global 1st-percentile threshold
q01 <- quantile(
  national$capacity_factor,
  probs = 0.05,
  na.rm = TRUE
) |>
  as.numeric()

# 2. Identify contiguous low-wind events
low_wind_hours <- national %>%
  arrange(time) %>%
  mutate(
    low_wind = !is.na(capacity_factor) & capacity_factor < q01,
    previous_low = lag(low_wind, default = FALSE),
    time_gap_hours = as.numeric(
      difftime(time, lag(time), units = "hours")
    ),
    new_event = low_wind &
      (!previous_low |
        is.na(time_gap_hours) |
        time_gap_hours != 1),
    event_id = cumsum(new_event)
  ) %>%
  filter(low_wind)

# 3. Summarise each event
event_catalogue <- low_wind_hours |>
  group_by(event_id) |>
  summarise(
    start_time = min(time),
    end_time = max(time),
    duration_hours = n(),
    min_capacity_factor = min(capacity_factor, na.rm = TRUE),
    mean_capacity_factor = mean(capacity_factor, na.rm = TRUE),
    min_output_mw = min(potential_MW, na.rm = TRUE),
    mean_output_mw = mean(potential_MW, na.rm = TRUE),

    # Sum of deficits below q01, in capacity-factor-hours
    cumulative_shortfall_cf_hours = sum(
      pmax(q01 - capacity_factor, 0),
      na.rm = TRUE
    ),

    # Energy shortfall relative to q01 * installed capacity, in MWh.
    # Valid when potential and capacity are in MW and records are hourly.
    cumulative_shortfall_mwh = sum(
      pmax(q01 * capacity_MW - potential_MW, 0),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  arrange(start_time)

# 4. Inspect the catalogue
print(event_catalogue)

# 5. Optional: 20 events with the largest cumulative energy shortfall
worst_by_shortfall <- event_catalogue |>
  arrange(desc(cumulative_shortfall_mwh)) |>
  slice_head(n = 20)

print(worst_by_shortfall)

# 6. Optional: 20 longest events
longest_events <- event_catalogue |>
  arrange(desc(duration_hours)) |>
  slice_head(n = 20)

print(longest_events)

# Part 3: Low-wind event statistics and visualisation

# -------------------------------------------------------------------
# 0. Basic checks and date fields
# -------------------------------------------------------------------

required_event_cols <- c(
  "event_id",
  "start_time",
  "end_time",
  "duration_hours",
  "min_capacity_factor",
  "mean_capacity_factor",
  "min_output_mw",
  "mean_output_mw",
  "cumulative_shortfall_cf_hours",
  "cumulative_shortfall_mwh"
)

missing_cols <- setdiff(required_event_cols, names(event_catalogue))
if (length(missing_cols) > 0) {
  stop(
    "event_catalogue is missing columns: ",
    paste(missing_cols, collapse = ", ")
  )
}


event_stats <- event_catalogue %>%
  mutate(
    start_time = as.POSIXct(start_time, tz = "UTC"),
    end_time = as.POSIXct(end_time, tz = "UTC"),
    start_year = year(start_time),
    # Season based on event start month; meteorological seasons
    start_season = case_when(
      month(start_time) %in% c(12, 1, 2) ~ "Winter",
      month(start_time) %in% 3:5 ~ "Spring",
      month(start_time) %in% 6:8 ~ "Summer",
      month(start_time) %in% 9:11 ~ "Autumn",
      TRUE ~ NA_character_
    ),
    start_season = factor(
      start_season,
      levels = c("Winter", "Spring", "Summer", "Autumn")
    )
  ) %>%
  arrange(start_time) %>%
  mutate(
    # Time between the previous event's end and this event's start.
    # This is elapsed time, not time from one event start to the next.
    inter_event_hours = as.numeric(
      difftime(start_time, lag(end_time), units = "hours")
    )
  )


# -------------------------------------------------------------------
# 1. Event frequency by year and season
# -------------------------------------------------------------------

annual_event_stats <- event_stats %>%
  group_by(start_year) %>%
  summarise(
    n_events = n(),
    mean_duration_hours = mean(duration_hours, na.rm = TRUE),
    median_duration_hours = median(duration_hours, na.rm = TRUE),
    max_duration_hours = max(duration_hours, na.rm = TRUE),
    mean_min_capacity_factor = mean(min_capacity_factor, na.rm = TRUE),
    median_min_capacity_factor = median(min_capacity_factor, na.rm = TRUE),
    total_shortfall_mwh = sum(cumulative_shortfall_mwh, na.rm = TRUE),
    mean_shortfall_mwh = mean(cumulative_shortfall_mwh, na.rm = TRUE),
    .groups = "drop"
  )

seasonal_event_stats <- event_stats %>%
  group_by(start_season, .drop = FALSE) %>%
  summarise(
    n_events = n(),
    mean_duration_hours = mean(duration_hours, na.rm = TRUE),
    median_duration_hours = median(duration_hours, na.rm = TRUE),
    max_duration_hours = if (all(is.na(duration_hours))) {
      NA_real_
    } else {
      max(duration_hours, na.rm = TRUE)
    },
    mean_min_capacity_factor = mean(min_capacity_factor, na.rm = TRUE),
    median_min_capacity_factor = median(min_capacity_factor, na.rm = TRUE),
    total_shortfall_mwh = sum(cumulative_shortfall_mwh, na.rm = TRUE),
    mean_shortfall_mwh = mean(cumulative_shortfall_mwh, na.rm = TRUE),
    .groups = "drop"
  )

# Annual event frequency
p_annual_frequency <- ggplot(
  annual_event_stats,
  aes(x = start_year, y = n_events)
) +
  geom_col(fill = blues9[6]) +
  scale_x_continuous(breaks = pretty(annual_event_stats$start_year)) +
  labs(
    title = "Annual frequency of low-wind events",
    x = "Year of event start",
    y = "Number of events"
  ) +
  theme_minimal()
ggsave(
  p_annual_frequency,
  filename = file.path(lwe_output_figpath, "annual_event_frequency.png"),
  width = 8,
  height = 6,
  dpi = 300
)

# Seasonal event frequency
p_seasonal_frequency <- ggplot(
  seasonal_event_stats,
  aes(x = start_season, y = n_events)
) +
  geom_col(fill = blues9[6]) +
  labs(
    title = "Low-wind event frequency by season",
    x = "Meteorological season of event start",
    y = "Number of events"
  ) +
  theme_minimal()

ggsave(
  p_seasonal_frequency,
  filename = file.path(lwe_output_figpath, "seasonal_event_frequency.png"),
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------
# 2. Duration distributions
# -------------------------------------------------------------------

p_duration_hist <- ggplot(
  event_stats,
  aes(x = duration_hours)
) +
  geom_histogram(bins = 30, fill = blues9[6]) +
  labs(
    title = "Distribution of low-wind event durations",
    x = "Event duration (hours)",
    y = "Number of events"
  ) +
  theme_minimal()
ggsave(
  p_duration_hist,
  filename = file.path(lwe_output_figpath, "event_duration_histogram.png"),
  width = 8,
  height = 6,
  dpi = 300
)
p_duration_box_season <- ggplot(
  event_stats,
  aes(x = start_season, y = duration_hours)
) +
  geom_boxplot(na.rm = TRUE, fill = blues9[6]) +
  labs(
    title = "Event duration by season",
    x = "Meteorological season of event start",
    y = "Event duration (hours)"
  ) +
  theme_minimal()

ggsave(
  p_duration_box_season,
  filename = file.path(lwe_output_figpath, "event_duration_boxplot.png"),
  width = 8,
  height = 6,
  dpi = 300
)

# Log-scale version can help if event durations are strongly right-skewed.
p_duration_hist_log <- ggplot(
  event_stats %>% filter(duration_hours > 0),
  aes(x = duration_hours)
) +
  geom_histogram(bins = 6, fill = blues9[6]) +
  scale_x_log10(n.breaks = 6) +
  labs(
    title = "Distribution of low-wind event durations (log scale)",
    x = "Event duration (hours, log scale)",
    y = "Number of events"
  ) +
  theme_minimal()

ggsave(
  p_duration_hist_log,
  filename = file.path(lwe_output_figpath, "event_duration_histogram_log.png"),
  width = 8,
  height = 6,
  dpi = 300
)
# -------------------------------------------------------------------
# 3. Inter-event time distributions
# -------------------------------------------------------------------

# The first event has no preceding event, so its inter-event time is NA.
inter_event_stats <- event_stats %>%
  filter(!is.na(inter_event_hours), inter_event_hours >= 0)

p_inter_event_hist <- ggplot(
  inter_event_stats,
  aes(x = inter_event_hours)
) +
  geom_histogram(bins = 30, fill = blues9[6]) +
  labs(
    title = "Distribution of inter-event times",
    x = "Time from previous event end to next event start (hours)",
    y = "Number of inter-event intervals"
  ) +
  theme_minimal()
ggsave(
  p_inter_event_hist,
  filename = file.path(lwe_output_figpath, "inter_event_time_histogram.png"),
  width = 8,
  height = 6,
  dpi = 300
)
p_inter_event_hist_log <- ggplot(
  inter_event_stats %>% filter(inter_event_hours > 0),
  aes(x = inter_event_hours)
) +
  geom_histogram(bins = 30, fill = blues9[6]) +
  scale_x_log10() +
  labs(
    title = "Distribution of inter-event times (log scale)",
    x = "Inter-event time (hours, log scale)",
    y = "Number of inter-event intervals"
  ) +
  theme_minimal()
ggsave(
  p_inter_event_hist_log,
  filename = file.path(
    lwe_output_figpath,
    "inter_event_time_histogram_log.png"
  ),
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------
# 4. Severity distributions
# -------------------------------------------------------------------

p_min_cf_hist <- ggplot(
  event_stats,
  aes(x = min_capacity_factor)
) +
  geom_histogram(bins = 30, fill = blues9[6]) +
  labs(
    title = "Distribution of minimum capacity factor per event",
    x = "Minimum capacity factor",
    y = "Number of events"
  ) +
  theme_minimal()
ggsave(
  p_min_cf_hist,
  filename = file.path(lwe_output_figpath, "min_capacity_factor_histogram.png"),
  width = 8,
  height = 6,
  dpi = 300
)

p_shortfall_hist <- ggplot(
  event_stats %>% filter(cumulative_shortfall_mwh >= 0),
  aes(x = cumulative_shortfall_mwh)
) +
  geom_histogram(bins = 30, fill = blues9[6]) +
  labs(
    title = "Distribution of cumulative energy shortfall",
    x = "Cumulative shortfall (MWh)",
    y = "Number of events"
  ) +
  theme_minimal()
ggsave(
  p_shortfall_hist,
  filename = file.path(
    lwe_output_figpath,
    "cumulative_shortfall_histogram.png"
  ),
  width = 8,
  height = 6,
  dpi = 300
)

p_shortfall_hist_log <- ggplot(
  event_stats %>% filter(cumulative_shortfall_mwh > 0),
  aes(x = cumulative_shortfall_mwh)
) +
  geom_histogram(bins = 30, fill = blues9[6]) +
  scale_x_log10() +
  labs(
    title = "Distribution of cumulative energy shortfall (log scale)",
    x = "Cumulative shortfall (MWh, log scale)",
    y = "Number of events"
  ) +
  theme_minimal()
ggsave(
  p_shortfall_hist_log,
  filename = file.path(
    lwe_output_figpath,
    "cumulative_shortfall_histogram_log.png"
  ),
  width = 8,
  height = 6,
  dpi = 300
)

p_duration_vs_shortfall <- ggplot(
  event_stats,
  aes(
    x = duration_hours,
    y = cumulative_shortfall_mwh,
    colour = min_capacity_factor
  )
) +
  geom_point(alpha = 0.75) +
  scale_colour_viridis_c() +
  labs(
    title = "Event duration versus cumulative shortfall",
    x = "Event duration (hours)",
    y = "Cumulative shortfall (MWh)",
    colour = "Minimum\ncapacity factor"
  ) +
  theme_minimal()
ggsave(
  p_duration_vs_shortfall,
  filename = file.path(
    lwe_output_figpath,
    "duration_vs_cumulative_shortfall.png"
  ),
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------
# 5. Annual and seasonal severity comparisons
# -------------------------------------------------------------------

p_annual_shortfall <- ggplot(
  annual_event_stats,
  aes(x = start_year, y = total_shortfall_mwh)
) +
  geom_col(fill = blues9[6]) +
  labs(
    title = "Total low-wind event shortfall by year",
    x = "Year of event start",
    y = "Total cumulative shortfall (MWh)"
  ) +
  theme_minimal()
ggsave(
  p_annual_shortfall,
  filename = file.path(lwe_output_figpath, "annual_shortfall.png"),
  width = 8,
  height = 6,
  dpi = 300
)

p_seasonal_duration <- ggplot(
  event_stats,
  aes(x = start_season, y = duration_hours)
) +
  geom_boxplot(na.rm = TRUE, fill = blues9[6]) +
  labs(
    title = "Low-wind event duration by season",
    x = "Meteorological season of event start",
    y = "Event duration (hours)"
  ) +
  theme_minimal()
ggsave(
  p_seasonal_duration,
  filename = file.path(lwe_output_figpath, "seasonal_duration_boxplot.png"),
  width = 8,
  height = 6,
  dpi = 300
)

p_seasonal_shortfall <- ggplot(
  event_stats,
  aes(x = start_season, y = cumulative_shortfall_mwh)
) +
  geom_boxplot(na.rm = TRUE, fill = blues9[6]) +
  labs(
    title = "Cumulative shortfall by season",
    x = "Meteorological season of event start",
    y = "Cumulative shortfall (MWh)"
  ) +
  theme_minimal() +
  scale_y_log10(n.breaks = 6)

ggsave(
  p_seasonal_shortfall,
  filename = file.path(
    lwe_output_figpath,
    "seasonal_shortfall_boxplot_log.png"
  ),
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------
# 6. Optional: plot the worst events and longest events
# -------------------------------------------------------------------

# Use existing objects if available; otherwise derive them from event_stats.
if (!exists("worst_by_shortfall")) {
  worst_by_shortfall <- event_stats %>%
    arrange(desc(cumulative_shortfall_mwh)) %>%
    slice_head(n = 20)
}

if (!exists("longest_events")) {
  longest_events <- event_stats %>%
    arrange(desc(duration_hours)) %>%
    slice_head(n = 20)
}

# Add a rank for legible horizontal bar plots.
worst_plot_data <- worst_by_shortfall %>%
  arrange(cumulative_shortfall_mwh) %>%
  mutate(
    event_label = paste0("Event ", event_id),
    event_label = start_time %>% format("%Y-%m-%d")
  )

p_worst_events <- ggplot(
  worst_plot_data,
  aes(
    x = reorder(event_label, cumulative_shortfall_mwh),
    y = cumulative_shortfall_mwh
  )
) +
  geom_col(fill = blues9[6]) +
  coord_flip() +
  labs(
    title = "Events with the largest cumulative energy shortfall",
    x = NULL,
    y = "Cumulative shortfall (MWh)"
  ) +
  theme_minimal()

ggsave(
  p_worst_events,
  filename = file.path(lwe_output_figpath, "worst_events_shortfall.png"),
  width = 8,
  height = 6,
  dpi = 300
)


longest_plot_data <- longest_events %>%
  arrange(duration_hours) %>%
  mutate(
    event_label = paste0("Event ", event_id),
    event_label = start_time %>% format("%Y-%m-%d")
  )

p_longest_events <- ggplot(
  longest_plot_data,
  aes(
    x = reorder(event_label, duration_hours),
    y = duration_hours
  )
) +
  geom_col(fill = blues9[6]) +
  coord_flip() +
  labs(
    title = "Longest low-wind events",
    x = NULL,
    y = "Event duration (hours)"
  ) +
  theme_minimal()

ggsave(
  p_longest_events,
  filename = file.path(lwe_output_figpath, "longest_events_duration.png"),
  width = 8,
  height = 6,
  dpi = 300
)

# printing
cat("\nAnnual event statistics:\n")
print(annual_event_stats)

cat("\nSeasonal event statistics:\n")
print(seasonal_event_stats)

cat("\nOverall event statistics:\n")
print(
  event_stats %>%
    summarise(
      n_events = n(),
      median_duration_hours = median(duration_hours, na.rm = TRUE),
      mean_duration_hours = mean(duration_hours, na.rm = TRUE),
      p90_duration_hours = quantile(duration_hours, 0.90, na.rm = TRUE),
      max_duration_hours = max(duration_hours, na.rm = TRUE),
      median_inter_event_hours = median(inter_event_hours, na.rm = TRUE),
      median_shortfall_mwh = median(cumulative_shortfall_mwh, na.rm = TRUE),
      p90_shortfall_mwh = quantile(
        cumulative_shortfall_mwh,
        0.90,
        na.rm = TRUE
      ),
      max_shortfall_mwh = max(cumulative_shortfall_mwh, na.rm = TRUE)
    )
)


write.csv(
  event_stats,
  file.path(tables_path, "low_wind_event_stats.csv"),
  row.names = FALSE
)
write.csv(
  annual_event_stats,
  file.path(tables_path, "low_wind_annual_stats.csv"),
  row.names = FALSE
)
write.csv(
  seasonal_event_stats,
  file.path(tables_path, "low_wind_seasonal_stats.csv"),
  row.names = FALSE
)
