# REMIT data compilation ####

require(arrow)
require(dplyr)
require(tidyr)
require(lubridate)
require(rnaturalearth)
require(sf)
require(ggplot2)
require(ggsci)
require(ggthemes)
require(FNN)
require(data.table)
require(geosphere)

source("aux_funct.R")

## read data ####
data_path <- "~/Documents/elexon/data_by_year"


# test <- get_remit(year = 2024, path = data_path, end_time = "2024-01-31")

test <- get_remit(year = 2025, path = data_path, end_time = "2025-09-30")

# Data Download ####

years <- 2023:2024

# remit_full <- lapply(
#   years,
#   \(y) get_remit(year = y, path = data_path)
# )

## joining years####
# catalog
ref_catalog_2025 <- fread(
  file.path("data/ref_catalog_wind_2025.csv.gz")
)
# Data exploration
remit_df <- lapply(
  2019:2025,
  \(y) {
    read_parquet(
      file.path("~/Documents/elexon/", sprintf("remit_all_%d.parquet", y))
    )
  }
) %>%
  bind_rows()

# get BmUnit
remit_df <- remit_df %>%
  mutate(
    elexonBmUnit = case_when(
      # asset matches bmUnit
      assetId %in% wind.bmus.alt$elexonBmUnit ~ assetId,
      # asset matches with T_
      paste0("T_", assetId) %in% wind.bmus.alt$elexonBmUnit ~ paste0(
        "T_",
        assetId
      ),
      # affected unit matches
      toupper(affectedUnit) %in% wind.bmus.alt$elexonBmUnit ~ affectedUnit,
      # affected unit matches with T_
      paste0("T_", toupper(affectedUnit)) %in%
        wind.bmus.alt$elexonBmUnit ~ paste0("T_", toupper(affectedUnit)),
      # special cases
      grepl("LARYO", assetId) ~ paste0("T_", gsub("O", "W", assetId)),
      grepl("RAMP", assetId) ~ paste0("T_", gsub("RAMP", "RMPNO", assetId)),
      # else
      TRUE ~ NA
    ),
    inelexon = !is.na(elexonBmUnit)
  ) %>%
  left_join(
    ref_catalog_2025 %>% select(bmUnit, tech_typ, lon, lat),
    by = c("elexonBmUnit" = "bmUnit")
  )

remit_wf <- remit_df %>%
  filter(
    eventStatus == "Active",
    !is.na(elexonBmUnit),
    # grepl("Wind", fuelType)
    eventStartTime < eventEndTime,
    !is.na(normalCapacity)
  ) %>%
  mutate(
    mincapacity = ifelse(
      is.na(outageProfile),
      availableCapacity,
      purrr::map_dbl(outageProfile, ~ min(.x$capacity))
    ),
    capacity_impact = normalCapacity - mincapacity,
    capacity_imp_perc = capacity_impact / normalCapacity * 100
  ) %>%
  filter(capacity_impact > 0) %>%
  unique() %>%
  mutate(
    eventEndTime = as.POSIXct(
      eventEndTime,
      format = "%Y-%m-%dT%H:%M:%OSZ",
      tz = "UTC"
    ),
    eventStartTime = as.POSIXct(
      eventStartTime,
      format = "%Y-%m-%dT%H:%M:%OSZ",
      tz = "UTC"
    ),
    duration = as.numeric(difftime(
      eventEndTime,
      eventStartTime,
      units = "hours"
    ))
  )

# check missing
# remit_df %>%
#   filter(
#     eventStatus == "Active",
#     !inelexon,
#     # !is_asset_in_elexon,
#     grepl("Wind", fuelType)
#   ) %>%
#   pull(assetId) %>%
#   unique()

# remit_df %>%
#   filter(
#     eventStatus == "Active",
#     is.na(elexonBmUnit),
#     grepl("Wind", fuelType)
#   ) %>%
#   # pull(assetId) %>%
#   pull(affectedUnit) %>%
#   unique()

# c("RAMP", "T_RMPNO", "LARYO", "LARYW", "FALGW-1", "T_FALGW-1", )
# Figures ####
# event statust type
remit_df %>%
  ggplot() +
  geom_bar(
    aes(x = eventStatus),
    stat = "count",
    position = "dodge",
    fill = "darkblue"
  )

# active by year
remit_df %>%
  filter(eventStatus == "Active") %>%
  unique() %>%
  mutate(year = factor(year(eventStartTime))) %>%
  ggplot() +
  geom_bar(
    aes(x = year),
    stat = "count",
    position = "dodge",
    fill = "darkblue"
  )

# fuel type
remit_df %>%
  unique() %>%
  ggplot() +
  geom_bar(
    aes(x = forcats::fct_infreq(fuelType)),
    stat = "count",
    position = "dodge",
    fill = "darkblue"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18)
  )


## included in elexon power data ####
# Wind related events
remit_wf %>%
  mutate(year = factor(year(eventStartTime))) %>%
  ggplot() +
  geom_bar(
    aes(x = forcats::fct_infreq(tech_typ)),
    stat = "count",
    position = "dodge",
    fill = "darkblue"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18)
  )
# remit_df %>%
#   filter(
#     eventStatus == "Active",
#     !is.na(elexonBmUnit),
#     is.na(tech_typ)
#   ) %>%
#   pull(elexonBmUnit) %>%
#   unique()

# Wind related by year
remit_wf %>%
  mutate(year = factor(year(eventStartTime))) %>%
  ggplot() +
  geom_bar(
    aes(x = year, fill = tech_typ),
    stat = "count",
    position = "stack"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 18),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = mypalette) +
  labs(fill = "", x = "year", y = "count")

ggsave("fig/remit_events_by_year.pdf", width = 6, height = 4)

# event duration distribution
remit_wf %>%
  ggplot() +
  geom_histogram(aes(x = duration), bins = 30, fill = "darkblue") +
  scale_x_log10(
    breaks = c(
      5 / 60, # 5 minutes
      1, # 1 hour
      1 * 24, # 1 day
      1 * 24 * 7, # 1 week
      1 * 24 * 30, # ~1 month
      1 * 24 * 365 # 1 year (optional)
    ),
    labels = c(
      "5 min",
      "1 hour",
      "1 day",
      "1 week",
      "1 month",
      "1 year"
    )
  ) +
  labs(x = "Event duration (log scale)")
ggsave("fig/remit_event_duration.pdf", width = 6, height = 4)

# capacity impact####

remit_wf %>%
  ggplot() +
  geom_histogram(aes(x = capacity_imp_perc, fill = tech_typ), bins = 30) +
  labs(x = "Unavailable capacity (%)") +
  facet_wrap(~tech_typ) +
  scale_fill_manual(values = mypalette) +
  theme(legend.position = "none")
ggsave("fig/remit_capacity_impact_perc.eps", width = 6, height = 4)


# spatial distribution ####
map_dat <- remit_wf %>%
  group_by(lon, lat, elexonBmUnit, tech_typ) %>%
  summarise(
    n_events = n(),
    duration = mean(duration),
    capacity_imp_perc = mean(capacity_imp_perc)
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

map_dat %>%
  ggplot() +
  geom_sf(
    data = ne_countries(scale = "medium", returnclass = "sf") %>%
      filter(admin %in% c("United Kingdom", "Ireland")),
    fill = "gray90",
    color = "white"
  ) +
  geom_sf(
    aes(geometry = geometry, size = n_events, color = tech_typ)
  ) +
  scale_color_manual(values = mypalette) +
  labs(color = "", size = "Number of events") +
  ggthemes::theme_map() +
  theme(legend.position = "right")
ggsave("fig/remit_nevents_map.pdf", width = 4, height = 4)

map_dat %>%
  ggplot() +
  geom_sf(
    data = ne_countries(scale = "medium", returnclass = "sf") %>%
      filter(admin %in% c("United Kingdom", "Ireland")),
    fill = "gray90",
    color = "white"
  ) +
  geom_sf(
    aes(geometry = geometry, size = duration, color = tech_typ)
  ) +
  scale_color_manual(values = mypalette) +
  scale_size_continuous(
    trans = "log10",
    breaks = c(1, 24, 24 * 30, 24 * 365), # hours: 5 min, 1 h, 1 day, 1 month, 1 year
    labels = c("1 hour", "< 1 day", "< 1 month", "1 year +")
  ) +
  labs(color = "", size = "Duration") +
  ggthemes::theme_map() +
  theme(legend.position = "right")
ggsave("fig/remit_duration_map.pdf", width = 4, height = 4)

map_dat %>%
  ggplot() +
  geom_sf(
    data = ne_countries(scale = "medium", returnclass = "sf") %>%
      filter(admin %in% c("United Kingdom", "Ireland")),
    fill = "gray90",
    color = "white"
  ) +
  geom_sf(
    aes(geometry = geometry, size = capacity_imp_perc, color = tech_typ)
  ) +
  scale_color_manual(values = mypalette) +
  labs(color = "", size = "Capacity impact (%)") +
  ggthemes::theme_map() +
  theme(legend.position = "right")
ggsave("fig/remit_cap_imp_map.pdf", width = 4, height = 4)
# outgage profile extraction ####
# remit_wf <- remit_df %>%
#   filter(
#     eventStatus == "Active",
#     !is.na(elexonBmUnit),
#     # grepl("Wind", fuelType)
#   ) %>%
#   unique() %>%
#   unnest(outageProfile, keep_empty = TRUE) %>%
#   mutate(
#     startTime = if_else(is.na(startTime), eventStartTime, startTime) %>%
#       ymd_hms(., tz = "UTC"),
#     endTime = if_else(is.na(endTime), eventEndTime, endTime) %>%
#       ymd_hms(., tz = "UTC"),
#     capacity = if_else(is.na(capacity), availableCapacity, capacity)
#   ) %>%
#   filter(!is.na(capacity))

# joining with generation adjusted by curtailment ####
# done in elexon_adj_data.R
