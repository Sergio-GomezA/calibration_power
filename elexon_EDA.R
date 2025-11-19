require(tidyr)
require(dplyr)
require(arrow)
require(lubridate)
require(ggplot2)
require(ggsci)
require(rnaturalearth)
require(rnaturalearthdata)
require(sf)


# Data reading ####

path <- "~/Documents/elexon/"

gen_adj <- read_parquet(
  file.path(path, "gen_adj_v2.parquet")
)
ref_catalog_2025 <- fread(
  file.path("data/ref_catalog_wind_2025.csv.gz")
)

# Figures ####

## Curtailment ####
source("wind_farm_data.R")
### monthly ####
curt_mdata <- gen_adj %>%
  filter(
    bmUnit %in% ref_catalog_2025$bmUnit
  ) %>%
  mutate(
    curt_perc = curtailment / capacity,
    month = month(halfHourEndTime)
  ) %>%
  group_by(month, bmUnit) %>%
  summarise(curt_perc = mean(curt_perc, na.rm = TRUE)) %>%
  left_join(
    ref_catalog_2025 %>% select(bmUnit, tech_typ, lon, lat),
    by = "bmUnit"
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_join(uk["geonunit"]) %>%
  st_drop_geometry() %>%
  group_by(geonunit, month) %>%
  summarise(
    curt_perc = mean(curt_perc, na.rm = TRUE) * 100
  )
curt_mdata %>%
  ggplot() +
  geom_line(aes(x = month, y = curt_perc, color = geonunit)) +
  # scale_color_manual(values = mypalette) +
  scale_x_continuous(breaks = 1:6 * 2) +
  theme(legend.position = "bottom") +
  labs(x = "month", y = "curtailment % of capacity", col = "")

### intraday ####

curt_pdata <- gen_adj %>%
  filter(
    bmUnit %in% ref_catalog_2025$bmUnit
  ) %>%
  mutate(
    curt_perc = curtailment / capacity,
    hour = hour(halfHourEndTime)
  ) %>%
  group_by(hour, bmUnit) %>%
  summarise(curt_perc = mean(curt_perc, na.rm = TRUE)) %>%
  left_join(
    ref_catalog_2025 %>% select(bmUnit, tech_typ),
    by = "bmUnit"
  ) %>%
  group_by(tech_typ, hour) %>%
  summarise(
    curt_perc = mean(curt_perc, na.rm = TRUE) * 100
  )

curt_pdata %>%
  ggplot() +
  geom_line(aes(x = hour, y = curt_perc, color = tech_typ)) +
  scale_color_manual(values = mypalette) +
  scale_x_continuous(breaks = 0:6 * 4) +
  theme(legend.position = "bottom") +
  labs(x = "hour", y = "curtailment % of capacity", col = "")

### map ####
curt_mapdata <- gen_adj %>%
  filter(
    bmUnit %in% ref_catalog_2025$bmUnit
  ) %>%
  mutate(
    curt_perc = curtailment / capacity * 100,
    hour = hour(halfHourEndTime)
  ) %>%
  group_by(bmUnit) %>%
  summarise(curt_perc = mean(curt_perc, na.rm = TRUE)) %>%
  left_join(
    ref_catalog_2025 %>% select(bmUnit, tech_typ, lon, lat),
    by = "bmUnit"
  ) %>%
  group_by(lon, lat, tech_typ) %>%
  summarise(curt_perc = mean(curt_perc, na.rm = TRUE)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

uk_ir <- ne_countries(scale = "medium") %>%
  filter(admin %in% c("United Kingdom", "Ireland")) %>%
  st_as_sf()

ggplot() +
  geom_sf(data = uk_ir) +
  geom_sf(
    aes(geometry = geometry, col = tech_typ, size = curt_perc),
    alpha = 0.5,
    data = curt_mapdata
  ) +
  scale_color_manual(values = mypalette) +
  scale_size(breaks = c(1, 2, 5, 10)) +
  ggthemes::theme_map() +
  labs(col = "", size = "curtailment\n% of capacity") +
  theme(legend.position = "right")

ggsave(
  filename = "fig/curtailment_map.svg",
  width = 4,
  height = 4
)
