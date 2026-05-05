# ELEXON + ERA5 Data ------------------------------

data_path <- "~/Documents/ERA5_at_wf/"
gen_path <- "~/Documents/elexon/"
model_path <- "~/Documents/elexon/model_objects"

# wind farm catalog based on 2025 data
ref_catalog_2025 <- fread(
  file.path("data/ref_catalog_wind_2025_era.csv.gz")
)

# era 5 coordinates list
coords_tb <- read.csv("data/era5_loc_mapping.csv")
# elexon wind bmus with updated capacity
wind.bmus.alt <- read.csv("data/wind_bmu_alt.csv")

# potential energy historical summary
pot_summary <- read.csv("data/hist_pot2024.csv")

generic_pc <- fread("data/generic_powerCurves.csv.gz")

class_curves <- fread("data/generic_powerCurves.csv.gz") %>%
  group_by(class) %>%
  mutate(power_scaled = power_kw / ratedPower) %>%
  ungroup()

# pwr_curv_df <- read_parquet(file.path(gen_path, "power_curve.parquet"))

# pwr_curv_df <- read_parquet(file.path(gen_path, "power_curve_all.parquet")) %>%
#   rename(time = halfHourEndTime)

# ## aggregated to GB level ####
# GB_df <- pwr_curv_df %>%
#   # filter(time >= "2023-01-01", time < "2025-01-01") %>%
#   group_by(tech_typ, time) %>%
#   summarise(
#     # ws_h_wmean = sum(ws_h * capacity),
#     across(c(ws_h, wd10, wd100), ~ sum(. * capacity), .names = "{.col}_wmean"),
#     across(
#       c(power_est0, potential, capacity),
#       sum
#     ),
#     across(c(ws_h, wd10, wd100), mean, .names = "{.col}_mean")
#   ) %>%
#   mutate(
#     across(
#       c(power_est0, potential),
#       ~ . / capacity,
#       .names = "norm_{.col}"
#     ),
#     # ws_h_wmean = ws_h_wmean / capacity,
#     across(matches("_wmean"), ~ . / capacity),
#     month = factor(month(time)),
#     hour = factor(hour(time))
#   ) %>%
#   mutate(ws_group = inla.group(ws_h_wmean, n = 20, method = "quantile")) %>%
#   group_by(tech_typ) %>%
#   arrange(time, .by_group = TRUE) %>%
#   mutate(t = row_number()) %>%
#   ungroup()

# write_parquet(GB_df, file.path(gen_path, "GB_aggr_all.parquet"))
GB_df <- read_parquet(file.path(gen_path, "GB_aggr.parquet")) %>%
  rename(time = halfHourEndTime) %>%
  mutate(
    err = norm_power_est0 - norm_potential,
    date = as.Date(time)
  )


gb_daily_df <- GB_df %>%
  group_by(tech_typ, date) %>%
  summarise(
    err = sum(err * capacity) / sum(capacity),
    norm_power_est0 = sum(norm_power_est0 * capacity) / sum(capacity),
    norm_potential = sum(norm_potential * capacity) / sum(capacity),
    ws_h_wmean = sum(ws_h_wmean * capacity) / sum(capacity),
    capacity = sum(capacity)
  )


# gb_daily_df %>%
#   ggplot(aes(date, err, color = tech_typ)) +
#   geom_line() +
#   theme_bw()

# gb_daily_df %>%
#   ggplot(aes(norm_potential, err, color = tech_typ)) +
#   geom_point() +
#   theme_bw() +
#   scale_color_manual(values = mypalette)

# pwr_curv_df$halfHourEndTime %>% range()

# elexon_adj_data.R
# bmu_df: elexon all years
# gen_adj_v2: -remit events + curtailment
# historical_power_curve.R
# pwr_curv_df: hourly aggregation and era5 variables + vertical interpolation + power estimate 0

# NAO, AO, EA, SCAN indices

nao_df <- read.csv("data/norm.daily.nao.cdas.z500.19500101_current.csv") %>%
  mutate(
    date = as.Date(
      sprintf("%04d-%02d-%02d", year, month, day)
    )
  ) %>%
  rename(nao = nao_index_cdas)


ao_df <- read.csv("data/norm.daily.ao.cdas.z1000.19500101_current.csv") %>%
  mutate(
    date = as.Date(
      sprintf("%04d-%02d-%02d", year, month, day)
    )
  ) %>%
  rename(ao = ao_index_cdas)


gb_daily_df <- gb_daily_df %>%
  left_join(
    nao_df %>% dplyr::select(date, nao),
    by = "date"
  ) %>%
  left_join(
    ao_df %>% dplyr::select(date, ao),
    by = "date"
  ) %>%
  mutate(
    nao_group = inla.group(nao, n = 5, method = "quantile"),
    ao_group = inla.group(ao, n = 5, method = "quantile")
  )


nao_monthly_df <- fread("data/norm.nao.monthly.b5001.current.ascii") %>%
  setNames(c("year", "month", "nao")) %>%
  mutate(date = as.Date(sprintf("%04d-%02d-01", year, month))) %>%
  dplyr::select(date, nao)

ao_monthly_df <- fread("data/monthly.ao.index.b50.current.ascii") %>%
  setNames(c("year", "month", "ao")) %>%
  mutate(date = as.Date(sprintf("%04d-%02d-01", year, month))) %>%
  dplyr::select(date, ao)
ea_df <- read.csv("data/ea.csv") %>%
  setNames(c("date", "ea")) %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
  tibble()
scan_df <- read.csv("data/scand.csv") %>%
  setNames(c("date", "scan")) %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
  tibble()


gb_monthly_df <- GB_df %>%
  mutate(date = lubridate::floor_date(date, "month")) %>%
  group_by(tech_typ, date) %>%
  summarise(
    err = sum(err * capacity) / sum(capacity),
    norm_power_est0 = sum(norm_power_est0 * capacity) / sum(capacity),
    norm_potential = sum(norm_potential * capacity) / sum(capacity),
    ws_h_wmean = sum(ws_h_wmean * capacity) / sum(capacity),
    capacity = sum(capacity)
  ) %>%
  left_join(
    nao_monthly_df,
    by = "date"
  ) %>%
  left_join(
    ao_monthly_df %>% dplyr::select(date, ao),
    by = "date"
  ) %>%
  left_join(
    ea_df %>% dplyr::select(date, ea),
    by = "date"
  ) %>%
  left_join(
    scan_df %>% dplyr::select(date, scan),
    by = "date"
  ) %>%
  mutate(
    nao_group = inla.group(nao, n = 5, method = "quantile"),
    ao_group = inla.group(ao, n = 5, method = "quantile"),
    ea_group = inla.group(ea, n = 5, method = "quantile"),
    scan_group = inla.group(scan, n = 5, method = "quantile")
  )


GB_df <- GB_df %>%
  mutate(month0 = lubridate::floor_date(date, "month")) %>%
  left_join(
    nao_df %>% dplyr::select(date, nao),
    by = "date"
  ) %>%
  left_join(
    ao_df %>% dplyr::select(date, ao),
    by = "date"
  ) %>%
  left_join(
    ea_df %>% dplyr::select(date, ea),
    by = c("month0" = "date")
  ) %>%
  left_join(
    scan_df %>% dplyr::select(date, scan),
    by = c("month0" = "date")
  ) %>%
  mutate(
    nao_group = inla.group(nao, n = 5, method = "quantile"),
    ao_group = inla.group(ao, n = 5, method = "quantile"),
    ea_group = inla.group(ea, n = 5, method = "quantile"),
    scan_group = inla.group(scan, n = 5, method = "quantile")
  )
