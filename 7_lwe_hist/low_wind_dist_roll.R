# ============================================================
# National wind generation: distribution and lower-tail analysis
# Required columns: bmUnit, time, potential, capacity
# ============================================================

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

rolling_window <- 12L
q_low <- 0.1
event_quantile <- 0.05 # Preserve the original event cutoff.
require_complete_hours <- FALSE

stopifnot(
  length(rolling_window) == 1L,
  is.finite(rolling_window),
  rolling_window >= 1,
  rolling_window == floor(rolling_window)
)
rolling_window <- as.integer(rolling_window)
cf_label <- sprintf(
  "National capacity factor (%d-hour rolling mean)",
  rolling_window
)

pow_threshold <- 0.1
pow_threshold_label <- gsub("\\.", "_", as.character(pow_threshold))
tables_path <- "~/ownCloud-s2441782@pwr_curv_dfasync.ed.ac.uk/projects/calibration/calibration_power_main_doc/tables"
dir.exists(tables_path) || dir.create(tables_path, recursive = TRUE)

override_objects <- TRUE
batch_name <- "bias_eda_roll"
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

library(data.table)
library(ggplot2)
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
coord_list_fname <- "data/coord_list.csv"

cat("Loading existing coordinate list\n")
coord_list <- read.csv(coord_list_fname)

ref_catalog_2025 <- fread(
  file.path("data/ref_catalog_wind_2025.csv.gz")
)

input_columns <- c(
  "bmUnit",
  "settlementDate",
  "settlementPeriod",
  "halfHourEndTime",
  "potential",
  "capacity",
  "power_est0",
  "ws_h"
)
pwr_curv_df <- arrow::read_parquet(
  file.path(gen_path, "power_curve_all_enriched.parquet"),
  col_select = tidyselect::all_of(input_columns)
)
setDT(pwr_curv_df)

# Preserve the settlement-period reconstruction and period-50 exception.
pwr_curv_df[,
  time := as.POSIXct(as.Date(settlementDate), tz = "UTC") +
    (settlementPeriod - 2L) * 1800
]
pwr_curv_df[
  settlementPeriod == 50L,
  time := as.POSIXct(halfHourEndTime, tz = "UTC")
]
pwr_curv_df[,
  c("settlementDate", "settlementPeriod", "halfHourEndTime") := NULL
]

if (anyNA(pwr_curv_df$bmUnit) || anyNA(pwr_curv_df$time)) {
  stop("Missing bmUnit or time values. Resolve these before aggregation.")
}

# if (anyDuplicated(pwr_curv_df, by = c("bmUnit", "time"))) {
#   duplicates <- pwr_curv_df[, .N, by = .(bmUnit, time)][N > 1L]
#   print(head(duplicates, 20L))
#   stop("Duplicate farm-hour records found. Resolve these before aggregating.")
# }

capacity_check <- pwr_curv_df[
  !is.na(capacity),
  .(
    n_capacity_values = uniqueN(capacity),
    min_capacity = min(capacity),
    max_capacity = max(capacity)
  ),
  by = bmUnit
]
print(capacity_check[n_capacity_values > 1L])

n_farms_total <- uniqueN(pwr_curv_df$bmUnit)
pwr_curv_df[,
  valid := is.finite(potential) & is.finite(capacity) & capacity > 0
]
pwr_curv_df[,
  valid_est := is.finite(power_est0) & is.finite(capacity) & capacity > 0
]

national <- pwr_curv_df[,
  .(
    # Farm-hour uniqueness was checked above.
    n_farms = .N,
    n_valid = sum(valid),
    potential_MW = sum(potential[valid]),
    capacity_MW = sum(capacity[valid]),
    power_est0_MW = sum(power_est0[valid_est]),
    ws_h_wmean = sum(ws_h * capacity) / sum(capacity)
  ),
  by = time
]
rm(pwr_curv_df)

national[, coverage_fraction := n_valid / n_farms_total]
national[, capacity_factor := NA_real_]
national[capacity_MW > 0, capacity_factor := potential_MW / capacity_MW]

print(summary(national$coverage_fraction))
print(summary(national$capacity_factor))
if (require_complete_hours) {
  national <- national[n_valid == n_farms_total]
}
if (any(national$capacity_factor > 1, na.rm = TRUE)) {
  warning("Capacity factor exceeds 1. Check units and capacity definitions.")
}

national <- national[is.finite(capacity_factor) & capacity_factor >= 0]
if (!nrow(national)) {
  stop("No valid national capacity-factor observations.")
}
setorder(national, time)

# Roll over elapsed hours in UTC and restart at any missing/invalid hour.
national[, block_id := cumsum(c(TRUE, diff(as.numeric(time)) != 3600))]
national[,
  capacity_factor_roll := {
    if (.N < rolling_window) {
      NA_real_
    } else {
      frollmean(capacity_factor, n = rolling_window, align = "right")
    }
  },
  by = block_id
]

# Use GB local time for calendar groupings, while retaining UTC timestamps.
season_levels <- c("Winter", "Spring", "Summer", "Autumn")
season_for_month <- c(
  rep("Winter", 2L),
  rep("Spring", 3L),
  rep("Summer", 3L),
  rep("Autumn", 3L),
  "Winter"
)
local_time <- lubridate::with_tz(national$time, "Europe/London")
national[, `:=`(
  month = lubridate::month(local_time, label = TRUE, abbr = TRUE),
  month_num = lubridate::month(local_time),
  year = lubridate::year(local_time),
  quarter = lubridate::quarter(local_time, type = "date_last")
)]
national[,
  season := factor(season_for_month[month_num], levels = season_levels)
]
rm(local_time)

# The hourly output keeps the raw series and rolling warm-up rows.
national_rolling <- national[is.finite(capacity_factor_roll)]
if (!nrow(national_rolling)) {
  stop("No complete ", rolling_window, "-hour rolling windows are available.")
}

# ------------------------------------------------------------
# 3. National distribution and lower quantiles
# ------------------------------------------------------------

national_quantiles <- quantile(
  national_rolling$capacity_factor_roll,
  probs = c(0.01, 0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90),
  na.rm = TRUE,
  names = TRUE
)

print(national_quantiles)

# The 1% threshold in capacity-factor units
q01 <- unname(
  quantile(national_rolling$capacity_factor_roll, q_low, na.rm = TRUE)
)
q05 <- unname(
  quantile(national_rolling$capacity_factor_roll, 0.05, na.rm = TRUE)
)

cat("Rolling 10% capacity-factor threshold:", q01, "\n")
cat("Rolling 5% capacity-factor threshold:", q05, "\n")
cat("Equivalent percentage of installed capacity:", 100 * q01, "%\n")

# Distribution: histogram with 1% threshold
p_distribution <- ggplot(
  national_rolling,
  aes(x = capacity_factor_roll)
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
      "Dashed line: empirical 1st decile = %.3f\nDashed line (red): empirical 5th percentile = %.3f",
      q01,
      q05
    ),
    x = cf_label,
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
  national_rolling,
  aes(x = capacity_factor_roll)
) +
  stat_ecdf(linewidth = 0.8) +
  geom_hline(yintercept = q_low, linetype = "dashed") +
  geom_vline(xintercept = q01, linetype = "dashed") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Empirical distribution of national capacity factor",
    x = cf_label,
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
  national_rolling,
  aes(x = month, y = capacity_factor_roll)
) +
  geom_boxplot(
    outlier.alpha = 0.15,
    na.rm = TRUE
  ) +
  labs(
    title = "Monthly distribution of national wind capacity factor",
    x = "Calendar month",
    y = cf_label
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
  national_rolling,
  aes(x = capacity_factor_roll)
) +
  geom_density(na.rm = TRUE, fill = blues9[6]) +
  facet_wrap(~month, ncol = 4) +
  labs(
    title = "National capacity-factor distributions by month",
    x = cf_label,
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
monthly_quantiles <- national_rolling[,
  .(
    n_hours = .N,
    q01 = quantile(capacity_factor_roll, 0.01),
    q05 = quantile(capacity_factor_roll, 0.05),
    q10 = quantile(capacity_factor_roll, 0.10),
    median = median(capacity_factor_roll),
    mean = mean(capacity_factor_roll)
  ),
  keyby = .(month, month_num)
]

print(monthly_quantiles)

# Plot monthly lower quantiles
monthly_quantiles_long <- melt(
  monthly_quantiles,
  id.vars = "month",
  measure.vars = c("q01", "q05", "q10"),
  variable.name = "quantile",
  value.name = "capacity_factor"
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
    y = cf_label,
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
  national_rolling,
  aes(x = season, y = capacity_factor_roll)
) +
  geom_boxplot(outlier.alpha = 0.15) +
  labs(
    title = "Seasonal distribution of national wind capacity factor",
    x = "Season",
    y = cf_label
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

seasonal_quantiles <- national_rolling[,
  .(
    n_hours = .N,
    q01 = quantile(capacity_factor_roll, 0.01),
    q05 = quantile(capacity_factor_roll, 0.05),
    q10 = quantile(capacity_factor_roll, 0.10),
    median = median(capacity_factor_roll),
    mean = mean(capacity_factor_roll)
  ),
  keyby = season
]

print(seasonal_quantiles)

# ------------------------------------------------------------
# 6. Annual distributions and annual quantiles
# ------------------------------------------------------------

# Check whether lower-tail behaviour differs across years.
annual_quantiles <- national_rolling[,
  .(
    n_hours = .N,
    q01 = quantile(capacity_factor_roll, 0.01),
    q05 = quantile(capacity_factor_roll, 0.05),
    q10 = quantile(capacity_factor_roll, 0.10),
    median = median(capacity_factor_roll),
    mean = mean(capacity_factor_roll)
  ),
  keyby = year
]

print(annual_quantiles)

p_annual <- ggplot(
  annual_quantiles,
  aes(x = year, y = q05)
) +
  geom_line(col = blues9[7]) +
  geom_line(aes(y = q10), colour = blues9[4]) +
  # geom_point() +
  labs(
    title = "Annual variation in the 1st decile and 5th percentile",
    x = "Year",
    y = cf_label
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 0.15))

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
annualquarter_quantiles <- national_rolling[,
  .(
    n_hours = .N,
    q01 = quantile(capacity_factor_roll, 0.01),
    q05 = quantile(capacity_factor_roll, 0.05),
    q10 = quantile(capacity_factor_roll, 0.10),
    median = median(capacity_factor_roll),
    mean = mean(capacity_factor_roll)
  ),
  keyby = quarter
]

print(annualquarter_quantiles)

p_annualq <- ggplot(
  annualquarter_quantiles,
  aes(x = quarter, y = q10)
) +
  geom_line(col = blues9[7]) +
  # geom_point() +
  geom_line(
    aes(y = q05),
    colour = blues9[4]
  ) +
  labs(
    title = "Quarterly variation in the 1st decile and 5th percentile",
    x = "Quarter (year + quarter)",
    y = cf_label
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 0.25)) #+
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

fwrite(
  national[, !"block_id"],
  file.path(lwe_output_tabpath, "national_hourly_generation.csv")
)
fwrite(
  monthly_quantiles,
  file.path(lwe_output_tabpath, "national_monthly_quantiles.csv")
)
fwrite(
  seasonal_quantiles,
  file.path(lwe_output_tabpath, "national_seasonal_quantiles.csv")
)
fwrite(
  annual_quantiles,
  file.path(lwe_output_tabpath, "national_annual_quantiles.csv")
)

fwrite(
  annualquarter_quantiles,
  file.path(lwe_output_tabpath, "national_quarterly_quantiles.csv")
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
event_threshold <- unname(quantile(
  national_rolling$capacity_factor_roll,
  probs = event_quantile
))
cat(sprintf(
  "Event threshold (quantile %.3f of %d-hour rolling CF) = %.4f\n",
  event_quantile,
  rolling_window,
  event_threshold
))

national_rolling[, low_wind := capacity_factor_roll < event_threshold]
national_rolling[,
  event_id := cumsum(
    low_wind &
      (!data.table::shift(low_wind, fill = FALSE) |
        (as.numeric(time) - data.table::shift(as.numeric(time), fill = -Inf) !=
          3600))
  )
]
low_wind_hours <- national_rolling[low_wind == TRUE]

# A window contributes one hour at its ending timestamp. Count physical
# hourly energy deficit once, even though the rolling windows overlap.
low_wind_hours[, `:=`(
  shortfall_cf_hours = pmax(event_threshold - capacity_factor_roll, 0),
  shortfall_mwh = pmax(event_threshold * capacity_MW - potential_MW, 0)
)]

event_catalogue <- low_wind_hours[,
  .(
    start_time = min(time),
    end_time = max(time),
    duration_hours = .N,
    min_capacity_factor = min(capacity_factor_roll),
    mean_capacity_factor = mean(capacity_factor_roll),
    min_output_mw = min(potential_MW),
    mean_output_mw = mean(potential_MW),
    cumulative_shortfall_cf_hours = sum(shortfall_cf_hours),
    cumulative_shortfall_mwh = sum(shortfall_mwh)
  ),
  keyby = event_id
]
print(event_catalogue)

worst_by_shortfall <- head(
  event_catalogue[order(-cumulative_shortfall_mwh)],
  20L
)
longest_events <- head(event_catalogue[order(-duration_hours)], 20L)
print(worst_by_shortfall)
print(longest_events)

# Part 3: Low-wind event statistics and visualisation
event_stats <- copy(event_catalogue)
setorder(event_stats, start_time)
event_start_local <- lubridate::with_tz(event_stats$start_time, "Europe/London")
event_stats[, `:=`(
  start_year = lubridate::year(event_start_local),
  start_season = factor(
    season_for_month[lubridate::month(event_start_local)],
    levels = season_levels
  ),
  inter_event_hours = as.numeric(difftime(
    start_time,
    data.table::shift(end_time),
    units = "hours"
  ))
)]

summarise_events <- function(events, group_column) {
  events[,
    .(
      n_events = .N,
      mean_duration_hours = mean(duration_hours),
      median_duration_hours = median(duration_hours),
      max_duration_hours = max(duration_hours),
      mean_min_capacity_factor = mean(min_capacity_factor),
      median_min_capacity_factor = median(min_capacity_factor),
      total_shortfall_mwh = sum(cumulative_shortfall_mwh),
      mean_shortfall_mwh = mean(cumulative_shortfall_mwh)
    ),
    keyby = group_column
  ]
}

annual_event_stats <- summarise_events(event_stats, "start_year")
seasonal_event_stats <- summarise_events(event_stats, "start_season")
seasonal_event_stats <- seasonal_event_stats[
  data.table(start_season = factor(season_levels, levels = season_levels)),
  on = "start_season"
]
seasonal_event_stats[
  is.na(n_events),
  `:=`(n_events = 0L, total_shortfall_mwh = 0)
]

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
  event_stats[duration_hours > 0],
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
inter_event_stats <- event_stats[
  !is.na(inter_event_hours) & inter_event_hours >= 0
]

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
  inter_event_stats[inter_event_hours > 0],
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
    x = sprintf("Minimum %d-hour rolling CF", rolling_window),
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
  event_stats[cumulative_shortfall_mwh >= 0],
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
  event_stats[cumulative_shortfall_mwh > 0],
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
    colour = sprintf("Minimum\n%d-hour rolling CF", rolling_window)
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

worst_plot_data <- worst_by_shortfall[order(cumulative_shortfall_mwh)]
worst_plot_data[, event_label := format(start_time, "%Y-%m-%d")]

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


longest_plot_data <- longest_events[order(duration_hours)]
longest_plot_data[, event_label := format(start_time, "%Y-%m-%d")]

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
print(event_stats[, .(
  n_events = .N,
  median_duration_hours = median(duration_hours),
  mean_duration_hours = mean(duration_hours),
  p90_duration_hours = quantile(duration_hours, 0.90),
  max_duration_hours = if (.N) max(duration_hours) else NA_real_,
  median_inter_event_hours = median(inter_event_hours, na.rm = TRUE),
  median_shortfall_mwh = median(cumulative_shortfall_mwh),
  p90_shortfall_mwh = quantile(cumulative_shortfall_mwh, 0.90),
  max_shortfall_mwh = if (.N) max(cumulative_shortfall_mwh) else NA_real_
)])


fwrite(
  event_stats,
  file.path(lwe_output_tabpath, "low_wind_event_stats.csv")
)
fwrite(
  annual_event_stats,
  file.path(lwe_output_tabpath, "low_wind_annual_stats.csv")
)
fwrite(
  seasonal_event_stats,
  file.path(lwe_output_tabpath, "low_wind_seasonal_stats.csv")
)

fwrite(
  event_catalogue,
  file.path(lwe_output_tabpath, "low_wind_event_catalogue.csv")
)
