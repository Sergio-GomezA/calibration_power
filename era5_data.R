# install.packages("arrow")
require(arrow)
require(dplyr)
require(rnaturalearth)
require(sf)
require(ggplot2)
require(ggthemes)
require(FNN)
require(data.table)

source("aux_funct.R")

## read data ####
data_path <- "~/Documents/ERA5_at_wf/"

era_df <- read_parquet(
  file.path(data_path, "era5_combined.parquet")
) %>%
  mutate(
    ws100 = sqrt(u100^2 + v100^2),
    ws10 = sqrt(u10^2 + v10^2),
    wd100 = (atan2(-u100, -v100) * 180 / pi) %% 360,
    wd10 = (atan2(-u10, -v10) * 180 / pi) %% 360
  )

ref_catalog_2025 <- fread("data/ref_catalog_wind_2025.csv.gz") %>%
  filter(quantity > 0)

## list of coordinates extracted
coords_tb <- era_df %>%
  select(longitude, latitude) %>%
  unique() %>%
  rename(era5lon = longitude, era5lat = latitude)

# uk border
uk_ie <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(admin %in% c("United Kingdom", "Ireland")) %>%
  st_transform(crs = 4326)

## map 0
uk_ie %>%
  ggplot() +
  geom_sf(data = uk_ie, fill = "gray95", color = "gray70") + # map background
  geom_point(data = coords_tb, aes(era5lon, era5lat), color = "gray50") +
  theme_map()
ggsave(
  "fig/era5_at_windfarms.eps",
  width = 4,
  height = 6,
  units = "in",
  dpi = 300
)
## mapping era 5 gridpoints to elexon locations ####

# Extract coordinate matrices
ref_mat <- as.matrix(ref_catalog_2025[, .(lon, lat)])
era5_mat <- as.matrix(coords_tb[, c("era5lon", "era5lat")])

# Find nearest ERA5 point for each catalog entry
nn <- get.knnx(era5_mat, ref_mat, k = 1)

# Add columns to catalog
ref_catalog_2025[, era5lon := coords_tb$era5lon[nn$nn.index]]
ref_catalog_2025[, era5lat := coords_tb$era5lat[nn$nn.index]]
ref_catalog_2025[, dist_deg := nn$nn.dist] # distance in degrees

## Group era5 data by tech_typ

coords_tb <- coords_tb %>%
  left_join(
    ref_catalog_2025 %>% select(era5lon, era5lat, tech_typ) %>% unique(),
    by = c("era5lon", "era5lat")
  ) %>%
  na.omit()

## plot ERA 5 series by aggregate ####

era_1ser <- era_df %>%
  group_by(
    time
  ) %>%
  summarise(
    ws_l = quantile(ws100, probs = 0.1),
    ws_h = quantile(ws100, probs = 0.9),
    across(
      c(ws100, ws10),
      mean
    )
  ) %>%
  mutate(
    year = year(time),
    quarter = paste0(year, "Q", quarter(time)),
    month = month(time)
  )

era_1ser %>%
  group_by(
    year,
    month
  ) %>%
  summarise(
    ws_l = quantile(ws100, probs = 0.05),
    ws_h = quantile(ws100, probs = 0.95),
    # across(c(ws_l, ws_h, ws100), mean)
    ws100 = mean(ws100)
  ) %>%
  mutate(date = as.Date(paste(year, month, 15, sep = "-"))) %>%
  ggplot(aes(date, ws100)) +
  geom_ribbon(aes(ymin = ws_l, ymax = ws_h), fill = "blue", alpha = 0.2) +
  geom_line(col = mypalette[1]) +
  labs(
    x = "Date",
    y = "wind speed m/s"
  )

ggsave("fig/ws_uk_1s.pdf", width = 6, height = 4, units = "in", dpi = 300)


era_1ser %>%
  group_by(
    year,
    month
  ) %>%
  summarise(
    ws100 = mean(ws100),
    .groups = "drop"
  ) %>%
  ggplot() +
  geom_line(aes(month, ws100, col = factor(year))) +
  labs(
    x = "month",
    y = "wind speed 100m m/s",
    col = "year"
  ) +
  scale_x_continuous(breaks = seq(2, 12, 2)) +
  scale_color_d3() +
  guides(col = guide_legend(ncol = 2)) + # make legend 2 columns
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.5, 0.78),
    legend.background = element_rect(fill = "transparent", color = NA), # transparent background
    legend.key = element_rect(fill = "transparent", color = NA) # transparent key boxes)
  )

ggsave("fig/ws_uk_1s_seas.eps", width = 6, height = 4, units = "in", dpi = 300)

## by tech type ####
era_tech <- era_df %>%
  left_join(
    coords_tb,
    by = c("longitude" = "era5lon", "latitude" = "era5lat")
  ) %>%
  group_by(
    tech_typ,
    time
  ) %>%
  summarise(
    across(
      c(ws100, ws10),
      mean
    )
  ) %>%
  mutate(
    year = year(time),
    quarter = paste0(year, "Q", quarter(time)),
    month = month(time)
  ) %>%
  filter(!is.na(tech_typ))

era_tech %>%
  group_by(
    tech_typ,
    year,
    month
  ) %>%
  summarise(
    ws_l = quantile(ws100, probs = 0.05),
    ws_h = quantile(ws100, probs = 0.95),
    # across(c(ws_l, ws_h, ws100), mean)
    ws100 = mean(ws100)
  ) %>%
  mutate(date = as.Date(paste(year, month, 15, sep = "-"))) %>%
  ggplot(aes(date, ws100, col = tech_typ, group = tech_typ)) +
  # geom_ribbon(aes(ymin = ws_l, ymax = ws_h), fill = "blue", alpha = 0.2) +
  geom_line() +
  scale_color_manual(values = mypalette) +
  labs(
    col = "",
    x = "Date",
    y = "wind speed m/s"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.7, 0.88), # place inside plot
    legend.background = element_rect(fill = "transparent", color = NA), # transparent background
    legend.key = element_rect(fill = "transparent", color = NA) # transparent key boxes
  )
ggsave("fig/ws_uk_tech.eps", width = 6, height = 4, units = "in", dpi = 300)

## wind direction plot ####

era_df %>%
  mutate(
    wd100_bin = cut(wd100, breaks = seq(0, 360, by = 10), include.lowest = TRUE)
  ) %>%
  count(wd100_bin) %>%
  mutate(mid_angle = seq(5, 355, by = 10)) %>% # center of bins
  ggplot(aes(x = mid_angle, y = n)) +
  geom_col(width = 10, fill = "steelblue", color = "white") +
  coord_polar(start = 0, direction = 1) + # make it circular
  scale_x_continuous(limits = c(0, 360), breaks = seq(0, 330, 30)) +
  labs(
    x = "Wind direction (° from North)",
    y = "Frequency",
    title = "Radial Histogram of Wind Direction (100 m)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )
ggsave("fig/wd_uk_1s.eps", width = 6, height = 4, units = "in", dpi = 300)


era_df %>%
  mutate(
    wd100_bin = cut(
      wd100,
      breaks = seq(0, 360, by = 10),
      include.lowest = TRUE
    ),
    month = factor(month(time), levels = 1:12, labels = month.abb) # month names
  ) %>%
  group_by(month) %>%
  count(wd100_bin) %>%
  mutate(mid_angle = seq(5, 355, by = 10)) %>% # center of bins
  ggplot(aes(x = mid_angle, y = n)) +
  geom_col(width = 10, fill = "steelblue", color = "white") +
  facet_wrap(~month) +
  coord_polar(start = 0, direction = 1) + # make it circular
  scale_x_continuous(limits = c(0, 360), breaks = seq(0, 330, 30)) +
  labs(
    x = "Wind direction (° from North)",
    y = "Frequency",
    title = "Radial Histogram of Wind Direction (100 m)"
  ) +
  # theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )
ggsave("fig/wd_uk_1s_seas.eps", width = 10, height = 6, units = "in", dpi = 300)
