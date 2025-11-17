# Elexon Data adjusted for curtailment and outages ####

# libraries ####

# read data ####

gen_adj <- bmu_df %>%
  select(
    settlementDate,
    settlementPeriod,
    bmUnit,
    halfHourEndTime,
    quantity
  ) %>%
  unique() %>%
  mutate(settlementDate = as.Date(settlementDate)) %>%
  # head(30) %>%
  left_join(
    curt_df %>%
      select(settlementDate, settlementPeriod, bmUnit, totalVolumeAccepted) %>%
      unique() %>%
      mutate(settlementDate = as.Date(settlementDate)),
    by = c("settlementDate", "settlementPeriod", "bmUnit")
  ) %>%
  mutate(
    curtailment = replace_na(-totalVolumeAccepted, 0),
    potential = quantity + curtailment
  ) %>%
  left_join(
    wind.bmus.alt %>%
      select(elexonBmUnit, bmUnitName, capacity) %>%
      unique(),
    by = c("bmUnit" = "elexonBmUnit")
  ) %>%
  filter(
    capacity > 0,
    potential > 0
  ) %>%
  mutate(
    cap_factor = ifelse(capacity > 0, potential / capacity)
  )

arrow::write_parquet(
  gen_adj,
  file.path(path, "gen_adj.parquet")
)
