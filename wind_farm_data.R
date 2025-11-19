require(data.table)
require(dplyr)
require(tidyr)
require(sf)

# Data ###############################

path <- "~/Documents/elexon"
elexon_path <- "/home/s2441782/Documents/elexon/"

## Renewable Energy Power Database ####
repd0 <- fread(file.path(
  path,
  "renew_energ_plan_db",
  "repd-q2-jul-2025.csv"
)) %>%
  rename(
    tech_typ = `Technology Type`,
    devel_stat = `Development Status (short)`,
    capacity_repd = `Installed Capacity (MWelec)`,
    capacity_turb = `Turbine Capacity`,
    n_turb = `No. of Turbines`,
    height_turb = `Height of Turbines (m)`,
    x_coord = `X-coordinate`,
    y_coord = `Y-coordinate`
  )
names(repd0) <- tolower(gsub(" ", "_", names(repd0)))

# filter wind only
repd_wind <- repd0 %>%
  filter(
    tech_typ %in% c("Wind Onshore", "Wind Offshore"),
    # devel_stat %in% c("Operational")
  ) %>%
  #   mutate(
  #     site_name = stringi::stri_enc_toutf8(as.character(site_name)),
  #     tech_typ  = stringi::stri_enc_toutf8(as.character(tech_typ))
  #   )
  mutate(
    site_name = iconv(site_name, from = "latin1", to = "UTF-8"),
    site_name = gsub("[^[:print:]]", "", site_name),
    ref_id = as.character(ref_id)
  )
rm(repd0)
## catalog v0 from REPD to ELEXON ########
id_catalog <- readxl::read_xlsx(
  path = file.path(path, "osuked", "id_dictionary.xlsx"),
  na = c("NA", "n/a", "-", "")
) %>%
  rename(
    rep_id = `REPD ID (New)`,
    bmu_id = `Settlement BMU ID`,
    dict_id = `Dictionary ID`,
    name = `Common Name`
  ) %>%
  filter(!is.na(rep_id)) %>%
  mutate(dict_id = as.character(dict_id))
names(id_catalog) <- tolower(gsub(" ", "_", names(id_catalog)))
bmu_to_repd <- id_catalog %>%
  tidyr::separate_rows(bmu_id, sep = ",\\s*") %>%
  tidyr::separate_rows(rep_id, sep = ",\\s*") %>%
  select(dict_id, name, bmu_id, rep_id) %>%
  unique()

## fixing mappings and manual maps ####
bmu_to_repd <- id_catalog %>%
  # fix aikengall
  mutate(rep_id = ifelse(rep_id == "7015", "4332, 4591", rep_id)) %>%
  # fix clyde wind farm: combine central and south / north: extension
  mutate(
    rep_id = ifelse(rep_id == "4011, 4178", "4011, 4178, 4011", rep_id)
  ) %>%
  # fix london array, all bmus to 2511, 2507 was abandoned
  mutate(
    rep_id = ifelse(
      rep_id == "2511, 2507",
      paste0(rep("2511", 4), collapse = ", "),
      rep_id
    )
  ) %>%
  # fix thanet. only THNTO are active
  mutate(
    rep_id = ifelse(rep_id == "6299, 2499", "2499, 2499, 6299, 6299", rep_id)
  ) %>%
  # fix walney: last id covers 2 BMUs  2500, 2506, 2533, 2533
  mutate(
    rep_id = ifelse(
      rep_id == "2533, 2500, 2506",
      "2500, 2506, 2533, 2533",
      rep_id
    ),
    bmu_id = ifelse(
      bmu_id == "T_WLNYW-1, E_WLNYW-2, T_WLNYO-3, T_WLNYO-4, T_WLNYO-2",
      "T_WLNYW-1, T_WLNYO-2, T_WLNYO-3, T_WLNYO-4",
      bmu_id
    )
  ) %>%
  # fix east anglia: EA1 and EA2 share BMUs
  mutate(
    rep_id = ifelse(rep_id == "2484, 6483, 2525, 2470", "2524, 2524", rep_id)
  ) %>%
  # fix hornsea 1
  mutate(
    rep_id = ifelse(rep_id == "2502, 2472, 2525", "2525, 2525, 2525", rep_id)
  ) %>%
  # add hornsea 2
  bind_rows(
    data.frame(
      dict_id = "HORSEA2",
      name = "Hornsea 2",
      bmu_id = c("T_HOWBO-1", "T_HOWBO-2", "T_HOWBO-3"),
      rep_id = c("2502", "2502", "2502")
    )
  ) %>%
  # fix clashindarroch
  mutate(rep_id = ifelse(rep_id == "3104", "4623", rep_id)) %>%
  # fix moray east
  mutate(rep_id = ifelse(rep_id == "2537", "2534", rep_id)) %>%
  # fix dorenell
  mutate(rep_id = ifelse(rep_id == "4244, 5051", "4244, 4244", rep_id)) %>%
  # add moray west
  bind_rows(
    data.frame(
      dict_id = "MORAYWEST",
      name = "Moray West",
      bmu_id = c("T_MOWWO-1", "T_MOWWO-2", "T_MOWWO-3", "T_MOWWO-4"),
      rep_id = c("6158")
    )
  ) %>%
  # update sanquhar
  mutate(rep_id = ifelse(rep_id == "6629", "4441", rep_id)) %>%
  bind_rows(
    data.frame(
      dict_id = "Neart na Gaoithe",
      name = "Neart na Gaoithe",
      bmu_id = c("T_NNGAO-1", "T_NNGAO-2"),
      rep_id = c("2522")
    )
  ) %>%
  # fix burbo bank
  mutate(rep_id = ifelse(rep_id == "2539, 2487", "2487, 2539", rep_id)) %>%
  mutate(across(c(bmu_id, rep_id), ~ str_split(.x, ",\\s*"))) %>%
  unnest(c(bmu_id, rep_id)) %>%
  filter(
    rep_id != "6797", #fix berry burn
    rep_id != "3619" # fix beinneun
  )

# new wind farms
manual_map <- read.csv("data/manual_bmu_to_repd.csv") %>%
  mutate(
    bmu_id = as.character(bmu_id),
    rep_id = as.character(rep_id)
  ) %>%
  filter(!is.na(rep_id)) %>%
  select(-c(rep_id2, rep_id3))
# write.csv(manual_map, "data/manual_bmu_to_repd.csv", row.names = FALSE)
bmu_to_repd <- bmu_to_repd %>%
  bind_rows(manual_map)

## elexon data ####
bmu_df <- fread(
  file.path(path, "data_by_year", "wind_gen_bmu_2025.csv.gz")
)
bmu_hours <- bmu_df %>%
  group_by(bmUnit) %>%
  arrange(halfHourEndTime) %>%
  mutate(quantity = quantity + lag(quantity)) %>%
  filter(minute(halfHourEndTime) == 0) %>%
  group_by(bmUnit) %>%
  summarise(
    quantity_raw = sum(quantity, na.rm = TRUE),
    n_hours = n_distinct(halfHourEndTime)
  )

bmu_y <- bmu_df %>%
  group_by(bmUnit) %>%
  arrange(halfHourEndTime) %>%
  mutate(quantity = quantity + lag(quantity)) %>%
  filter(quantity >= 0) %>%
  filter(minute(halfHourEndTime) == 0) %>%
  group_by(bmUnit) %>%
  summarise(
    quantity = sum(quantity, na.rm = TRUE),
    qpos_hours = n_distinct(halfHourEndTime)
  ) %>%
  arrange(desc(quantity)) %>%
  left_join(
    bmu_hours,
    by = "bmUnit"
  )
rm(bmu_df)
## adding REPD info ######
repd_vars <- c(
  "ref_id",
  "site_name",
  "tech_typ",
  "capacity_repd",
  "capacity_turb",
  "n_turb",
  "height_turb",
  "devel_stat",
  "x_coord",
  "y_coord",
  "country",
  "operational"
)

ref_catalog_2025 <- bmu_y %>%
  left_join(
    # wind bmus from elexon
    wind.bmus %>%
      select(
        elexonBmUnit,
        matches("name|capacity")
      ) %>%
      unique(),
    by = c("bmUnit" = "elexonBmUnit")
  ) %>%
  left_join(
    # mapping elexon to repd
    bmu_to_repd %>%
      select(bmu_id, rep_id, name),
    by = c("bmUnit" = "bmu_id")
  ) %>%
  left_join(
    # repd
    repd_wind %>%
      select(any_of(repd_vars)),
    by = c("rep_id" = "ref_id")
  ) %>%
  mutate(operational_date = as.Date(operational, format = "%d/%m/%Y")) %>%
  filter(!is.na(x_coord))

write.csv(
  ref_catalog_2025,
  gzfile(file.path("data/ref_catalog_wind_2025.csv.gz")),
  row.names = FALSE
)

## generic power curves ####

generic_pc_fname <- list.files(
  "data",
  pattern = "^generic.*\\.json$",
  full.names = TRUE
)
generic_pc <- lapply(
  generic_pc_fname,
  \(file) {
    dat <- fromJSON(file)
    # val_names <- names(dat)
    result <- as.data.frame(dat) %>%
      rename(
        wind_speed = powerCurveData.1,
        power_kw = powerCurveData.2,
      ) %>%
      mutate(across(airDensity, as.numeric))
    result
  }
) %>%
  bind_rows() %>%
  mutate(
    class = gsub(".* - ", "", title) %>% tolower()
  )

last_row <- generic_pc %>%
  group_by(class) %>%
  summarise(wind_speed = max(wind_speed) + 0.5) %>%
  mutate(power_kw = 0)

generic_pc <- bind_rows(generic_pc, last_row) %>%
  arrange(class, wind_speed) %>%
  group_by(class) %>%
  fill(everything(), .direction = "down") %>%
  ungroup()

write.csv(
  generic_pc,
  gzfile("data/generic_powerCurves.csv.gz"),
  row.names = FALSE
)

## adding GWA values to catalog ####
wf_loc <- ref_catalog_2025 %>%
  st_as_sf(coords = c("x_coord", "y_coord"), crs = 27700) %>%
  st_transform(crs = 4326) %>%
  st_coordinates()

ref_catalog_2025 <- fread(
  file.path("data/ref_catalog_wind_2025.csv.gz")
)
heights <- c(50, 100, 150, 200)
path <- "~/Documents/GWA/"
gwa_ws <- lapply(
  heights,
  \(h) {
    raster_fname <- sprintf("merged_wind-speed_%dm.tif", h)
    r <- rast(file.path(path, raster_fname))
    wind_values <- terra::extract(r, wf_loc)
    wind_values
  }
) %>%
  bind_cols() %>%
  setNames(paste0("gwa", heights))

gwa_A <- lapply(
  heights,
  \(h) {
    raster_fname <- sprintf("merged_combined-Weibull-A_%dm.tif", h)
    r <- rast(file.path(path, raster_fname))
    wind_values <- terra::extract(r, wf_loc)
    wind_values
  }
) %>%
  bind_cols() %>%
  setNames(paste0("gwa_A", heights))

gwa_k <- lapply(
  heights,
  \(h) {
    raster_fname <- sprintf("merged_combined-Weibull-k_%dm.tif", h)
    r <- rast(file.path(path, raster_fname))
    wind_values <- terra::extract(r, wf_loc)
    wind_values
  }
) %>%
  bind_cols() %>%
  setNames(paste0("gwa_k", heights))

ref_catalog_2025 <- ref_catalog_2025 %>%
  bind_cols(
    gwa_ws,
    gwa_A,
    gwa_k
  )
write.csv(
  ref_catalog_2025,
  gzfile(file.path("data/ref_catalog_wind_2025.csv.gz")),
  row.names = FALSE
)
set.seed(0)
## Interpolation of GWA
# filling in turbine characteristics
ref_catalog_2025 <- ref_catalog_2025 %>%
  # completing capacity per turbine
  mutate(
    capacity_turb = ifelse(
      is.na(capacity_turb),
      capacity_repd / n_turb,
      capacity_turb
    )
  ) %>%
  # imputing mean turbine height
  group_by(tech_typ) %>%
  mutate(
    height_turb_imp = ifelse(
      is.na(height_turb),
      mean(height_turb, na.rm = T),
      height_turb
    )
  ) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    # log interpolate GWA mean wind speed
    ws_ht = case_when(
      height_turb_imp <= 75 ~ interp_log_ws(
        height_turb_imp,
        50,
        100,
        gwa50,
        gwa100
      ),
      height_turb_imp <= 125 ~ interp_log_ws(
        height_turb_imp,
        100,
        150,
        gwa100,
        gwa150
      ),
      height_turb_imp <= 175 ~ interp_log_ws(
        height_turb_imp,
        150,
        200,
        gwa150,
        gwa200
      ),
      TRUE ~ interp_log_ws(height_turb_imp, 150, 200, gwa150, gwa200) # extrapolate slightly
    ),
    k_ht = approx(
      x = c(50, 100, 150, 200),
      y = c(gwa_k50, gwa_k100, gwa_k150, gwa_k200),
      xout = height_turb_imp,
      rule = 2 # extrapolate using end values if outside [50,200]
    )$y,
    A_ht = ws_ht / gamma(1 + 1 / k_ht)
  ) %>%
  ungroup() %>%
  mutate(
    turb_class = case_when(
      grepl("offshore", tech_typ, ignore.case = TRUE) ~ "offshore",
      ws_ht > 9.5 ~ "iec class 1",
      ws_ht > 8 ~ "iec class 2",
      TRUE ~ "iec class 3"
    )
  ) %>%
  group_by(rep_id) %>%
  mutate(
    n_bmu = n(),
    wf_cap_sum = sum(generationCapacity),
    capacity_diff = abs(capacity_repd - wf_cap_sum) / capacity_repd,
    turb_bmu = n_turb * capacity_repd / sum(capacity_repd),
    capacity_bmu = capacity_repd / n_bmu
  ) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    # power curve on mean wind speed PC(E(ws))
    power_mws = generic_pow_conv(
      ws_ht,
      turb_class = turb_class,
      turb_capacity = capacity_turb
    ),
    # mean power using weibull distribution E(PC(ws))
    mean_power = generic_pow_conv(
      wind_speed = rweibull(10000, shape = k_ht, scale = A_ht),
      turb_class = turb_class,
      turb_capacity = capacity_turb
    ) %>%
      mean()
  ) %>%
  ungroup() %>%
  mutate(
    gen_est0 = power_mws * turb_bmu * n_hours,
    gen_est = mean_power * turb_bmu * n_hours
  )

ref_catalog_2025 <- ref_catalog_2025 %>%
  st_as_sf(coords = c("x_coord", "y_coord"), crs = 27700) %>%
  st_transform(4326) %>%
  mutate(
    lon = st_coordinates(.)[, 1],
    lat = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()

write.csv(
  ref_catalog_2025,
  gzfile(file.path("data/ref_catalog_wind_2025.csv.gz")),
  row.names = FALSE
)
