# 0. Setup ####

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

# 0.1 global parameter #####
day_id <- 3
mesh_edge_par <- 20 # km, target edge length for the spatial mesh. 10 is fine, 20 is coarse but faster
override_objects <- TRUE
prec_init <- log(200)

if (local_run) {
  cat("Running in local mode\n")
} else {
  cat("Running in cluster mode\n")
}

# Get task ID and others from command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Override defaults only if arguments are provided
if (length(args) > 0) {
  day_id <- as.numeric(args[1])
}
if (length(args) > 1) {
  mesh_edge_par <- as.numeric(args[2])
}
if (length(args) > 2) {
  override_objects <- as.logical(args[3])
}

# 0.2 libraries and paths ####
require(parallel)

if (local_run) {
  data_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
  pixel_dims <- c(150, 150)
} else {
  data_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  gen_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  model_path <- "/exports/eddie/scratch/s2441782/calibration/model_objects"
  temp_lib <- "/exports/eddie3_homes_local/s2441782/lib"
  pixel_dims <- c(300, 300)
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
# require(ggspatial)

source("aux_funct.R")


# 1. data preparation ####

cat("Preparing data for model fitting\n")

sampled_days <- c("2020-08-14", "2024-04-17", "2024-04-12")

d0 <- sampled_days[day_id] %>% as.Date()
d0_tag <- base::format(d0, "%y%m%d")
df_pattern <- sprintf("^calibration_df_.*_%s\\.gpkg$", d0_tag)
files_found <- list.files("data", pattern = df_pattern, full.names = TRUE)

if (!override_objects && length(files_found) > 0) {
  cat(
    "Calibration data file already exists for this day. Loading existing data.\n"
  )
  wf_df_frag <- st_read(files_found[1])
} else {
  cat(
    "No existing calibration data file found for this day. Preparing new data.\n"
  )

  pwr_curv_df <- read_parquet(file.path(
    gen_path,
    "power_curve_all_enriched.parquet"
  ))

  n.days <- 1

  wf_df_frag <- pwr_curv_df %>%
    rename(time = halfHourEndTime) %>%
    mutate(
      date = as.Date(time),
      elevation = pmax(0, elevation),
      site_name = site_name %>%
        gsub("\\b(wind\\s*farm|wf)\\b", "", ., ignore.case = TRUE) %>%
        trimws()
    ) %>%
    # filter(date %in% sampled_days) %>%
    filter(date >= d0, date <= d0 + n.days - 1) %>%
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
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    mutate(lon = st_coordinates(.)[, 1], lat = st_coordinates(.)[, 2]) %>%
    st_transform(crs = 27700) %>%
    mutate(
      x = st_coordinates(.)[, 1] / 1000,
      y = st_coordinates(.)[, 2] / 1000,
    ) %>%
    # st_drop_geometry() %>%
    mutate(
      # site_id = as.integer(factor(site_name)),
      ws_group = inla.group(ws_h, n = 20, method = "quantile"),
      pow_group = inla.group(norm_power_est0, n = 20, method = "quantile"),
      d_coast_group = inla.group(dist_coast, n = 10, method = "quantile"),
      elev_group = inla.group(elevation, n = 10, method = "quantile"),
      time_id = as.integer(factor(time)),
      # loc = cbind(x, y)
    )

  x <- wf_df_frag$pow_group %>% unique() %>% sort()
  min_jump <- min(diff(sort(x))) / diff(range(x))
  if (min_jump <= 1e-4) {
    wf_df_frag <- wf_df_frag %>%
      mutate(pow_group = inla.group(norm_power_est0, n = 20, method = "cut"))
  }

  cat("Converting coordinates to km\n")
  wf_df_frag <- wf_df_frag %>%
    st_geometry() %>%
    (\(g) g / 1000)() %>%
    st_set_geometry(wf_df_frag, .)

  st_write(
    wf_df_frag,
    sprintf("data/calibration_df_%s_%s.gpkg", "base", d0_tag),
    driver = "GPKG",
    append = FALSE,
    quiet = TRUE
  )
}
cat("Number of unique locations:", nrow(wf_df_frag %>% distinct(x, y)), "\n")
n <- nrow(wf_df_frag)
# wf_df_frag %>%
#   ggplot() +
#   geom_point(aes(pow_group,error0), bins = 50)
# wf_df_frag %>%
#   ggplot() +
#   geom_point(aes(pow_group,error0), bins = 50)+
#   facet_wrap(~tech_typ, scales = "free_x")

## 1.0.1 GB daily summary ####

gb_day_df_fname <- sprintf("data/GB_daily_summary_%s.parquet", d0_tag)

if (!file.exists(gb_day_df_fname) || override_objects) {
  cat("GB daily summary file not found, creating new summary\n")
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

## 1.1 mesh building #####
cat("Building spatial mesh\n")

edge_target <- mesh_edge_par # km
mesh_label <- ifelse(edge_target >= 20, "coarse", "fine")
mesh_fname <- file.path(
  model_path,
  sprintf("spatial_mesh_%s_%s.rds", mesh_label, d0_tag)
)

if (!file.exists(mesh_fname) || override_objects) {
  cat("Mesh file not found, building new mesh\n")

  loc_unique <- wf_df_frag %>%
    distinct(x, y) %>%
    as.matrix()

  # bnd <- fm_extensions(loc_unique, convex = c(-.1, -.15))
  # bnd <- fm_extensions(loc_unique, convex = c(-.08, -.3))
  if (mesh_edge_par <= 20) {
    conv_par <- c(-.05, -.35)
    max_n <- c(900, 300)
  } else {
    conv_par <- c(-.1, -.35)
    max_n <- c(900, 150)
  }
  bnd <- fm_extensions(loc_unique, convex = conv_par)
  # bnd <- fm_extensions(loc_unique, convex = c(-.1, -.15))
  # ggplot() + geom_sf(data = bnd[[1]])
  bndin <- bnd[[1]]
  bndout <- bnd[[2]]

  uk_map <- rnaturalearth::ne_countries(
    scale = "medium",
    country = "United Kingdom",
    returnclass = "sf"
  )
  coastline <- uk_map %>%
    st_transform(crs = 27700) %>%
    st_boundary()
  uk_map <- uk_map %>%
    st_transform(crs = 27700) %>%
    st_geometry() %>%
    (\(g) g / 1000)() %>%
    st_set_geometry(uk_map, .)

  hex_0 <- fm_hexagon_lattice(bnd[[1]], edge_len = edge_target * 2)

  wf.mesh <- fm_mesh_2d(
    # loc = loc_unique,
    loc = hex_0,
    boundary = bnd,
    # max.edge = c(100, 150), # km
    min.angle = 30,
    # offset = -0.2,
    cutoff = edge_target,
    max.n.strict = max_n
  )
  saveRDS(
    wf.mesh,
    mesh_fname
  )
  ggplot() +
    geom_sf(data = uk_map, fill = NA, color = "black") +
    gg(wf.mesh) +
    geom_point(data = loc_unique, aes(x, y), color = "darkred") +
    annotation_scale(location = "bl", width_hint = 0.25, plot_unit = "km") +
    theme_void()
  ggsave(
    sprintf("fig/spatial_mesh_%s_%s.pdf", mesh_label, d0_tag),
    width = 4,
    height = 6
  )
} else {
  cat("Loading existing mesh\n")
  wf.mesh <- readRDS(mesh_fname)
}

### 1.2 mesh assessment #####
mesh_assess_fname <- sprintf(
  "fig/spatial_mesh_%s_assessment2_sddev_%s.pdf",
  mesh_label,
  d0_tag
)
if (!file.exists(mesh_assess_fname) || override_objects) {
  cat("Assessing spatial mesh\n")
  mesh_assessment <- fm_assess(mesh = wf.mesh, spatial.range = 70)

  ggplot() +
    geom_sf(
      data = mesh_assessment %>%
        st_filter(., bndout),
      aes(col = edge.len)
    ) +
    geom_point(data = loc_unique, aes(x, y), color = "darkred") +
    geom_sf(data = uk_map, fill = NA, color = "white") +
    annotation_scale(location = "bl", width_hint = 0.25, plot_unit = "km") +
    theme_void() +
    scale_color_viridis_c(option = "D")
  ggsave(
    sprintf(
      "fig/spatial_mesh_%s_assessment_edgelen_%s.pdf",
      mesh_label,
      d0_tag
    ),
    width = 4,
    height = 6
  )
  # sd.dev should be close to 1
  ggplot() +
    gg(
      data = mesh_assessment %>%
        st_filter(., bndout),
      aes(col = sd.dev)
    ) +
    geom_point(data = loc_unique, aes(x, y), color = "darkred") +
    geom_sf(data = uk_map, fill = NA, color = "white") +
    annotation_scale(location = "bl", width_hint = 0.25, plot_unit = "km") +
    theme_void() #+
  # scale_color_viridis_c(option = "D")
  ggsave(
    sprintf("fig/spatial_mesh_%s_assessment_sddev_%s.pdf", mesh_label, d0_tag),
    width = 4,
    height = 6
  )

  ggplot() +
    geom_sf(
      data = mesh_assessment %>%
        st_filter(., bndin),
      aes(col = edge.len)
    ) +
    geom_point(data = loc_unique, aes(x, y), color = "darkred") +
    geom_sf(data = uk_map, fill = NA, color = "white") +
    annotation_scale(location = "bl", width_hint = 0.25, plot_unit = "km") +
    theme_void() +
    scale_color_viridis_c(option = "D")
  ggsave(
    sprintf(
      "fig/spatial_mesh_%s_assessment2_edgelen_%s.pdf",
      mesh_label,
      d0_tag
    ),
    width = 4,
    height = 6
  )
  # sd.dev should be close to 1
  ggplot() +
    gg(
      data = mesh_assessment %>%
        st_filter(., bndin),
      aes(col = sd.dev)
    ) +
    geom_point(data = loc_unique, aes(x, y), color = "darkred") +
    geom_sf(data = uk_map, fill = NA, color = "white") +
    annotation_scale(location = "bl", width_hint = 0.25, plot_unit = "km") +
    theme_void() #+
  # scale_color_viridis_c(option = "D")
  ggsave(
    sprintf("fig/spatial_mesh_%s_assessment2_sddev_%s.pdf", mesh_label, d0_tag),
    width = 4,
    height = 6
  )
}

# 2. Model fitting ####
## 2.1 AR1 temporal model ####
ar_tag <- "ar1"
components0 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  # tech_typ(tech_typ, model = "iid") + # random intercept by tech_typ
  power_correction(
    pow_group,
    model = "rw2",
    # replicate = tech_typ,
    constr = TRUE
  ) + # smooth correction power
  d_coast(
    d_coast_group,
    model = "rw2",
    constr = TRUE
  ) + # smooth correction distance to coast
  elev(
    elev_group,
    model = "rw2",
    constr = TRUE
  ) + # smooth correction elevation
  wind(ws_group, model = "rw2", replicate = tech_typ, constr = TRUE) + # smooth correction wind
  u(
    t,
    model = "ar1",
    replicate = coord_id,
    hyper = list(
      rho = list(
        prior = "pc.cor1",
        param = c(0.6, 0.5)
      ),
      # prec = list(
      #   prior = "pc.prec",
      #   param = c(50, 0.05)
      # ),
      prec = list(initial = prec_init, fixed = TRUE)
    )
  )

model_fname <- file.path(
  model_path,
  sprintf("ts_bru0_%s_%s.rds", ar_tag, d0_tag)
)

if (!file.exists(model_fname) || override_objects) {
  cat("Fitting ar1 model\n")
  bruar1 <- bru(
    components = components0,
    formula = norm_potential ~ Intercept +
      power_correction +
      d_coast +
      elev +
      wind +
      u,
    family = "gaussian",
    control.family = list(
      hyper = list(
        prec = list(
          prior = "pc.prec",
          param = c(50, 0.05)
        )
      )
    ),
    data = wf_df_frag,
    options = bru_options(
      bru_verbose = 3,
      control.inla = list(verbose = TRUE)
    )
  )

  saveRDS(
    bruar1,
    file = model_fname
  )
} else {
  cat("Loading existing ar1 model\n")
  bruar1 <- readRDS(model_fname)
}

summary(bruar1)
# bruar1$summary.fixed[, 1:6]
# bruar1$summary.random$tech_typ[, 1:6]
# bruar1$summary.random$tech_power[, 1:6]

source("aux_funct.R")
effect_names <- names(bruar1$summary.random)
excluded_effects <- c("u")
effect_names <- setdiff(effect_names, excluded_effects)
for (effect in effect_names) {
  if (effect == "wind") {
    n_repl <- 2
    repl_names <- c("Offshore", "Onshore")
  } else {
    n_repl <- 1
    repl_names <- NULL
  }
  plot.effects(
    bruar1,
    effect,
    n.replicate = n_repl,
    replicate_names = repl_names,
    show.plot = TRUE
  )
  ggsave(
    sprintf("fig/%s_effect_%s_%s.pdf", effect, ar_tag, d0_tag),
    width = 6,
    height = 4
  )
}

plot.hyper.dens(bruar1)
ggsave(
  sprintf("fig/hyperparameters_%s_%s.pdf", ar_tag, d0_tag),
  width = 6,
  height = 4
)

# plot(bruar1$summary.fitted.values$mean[1:n], wf_df_frag$norm_potential)
# plot(bruar1$summary.random$u$mean, wf_df_frag$norm_potential)
## 2.2 AR2 temporal model ####
ar_tag <- "ar2"
components0 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  # tech_typ(tech_typ, model = "iid") + # random intercept by tech_typ
  power_correction(
    pow_group,
    model = "rw2",
    # replicate = tech_typ,
    constr = TRUE
  ) + # smooth correction power
  d_coast(
    d_coast_group,
    model = "rw2",
    constr = TRUE
  ) + # smooth correction distance to coast
  elev(
    elev_group,
    model = "rw2",
    constr = TRUE
  ) + # smooth correction elevation
  wind(ws_group, model = "rw2", replicate = tech_typ, constr = TRUE) + # smooth correction wind
  u(
    t,
    model = "ar",
    order = 2,
    replicate = coord_id,
    hyper = list(
      # rho = list(
      #   prior = "pc.cor1",
      #   param = c(0.6, 0.5)
      # ),
      # prec = list(
      #   prior = "pc.prec",
      #   param = c(50, 0.05)
      # )
      prec = list(initial = prec_init, fixed = TRUE)
    )
  )

model_fname <- file.path(
  model_path,
  sprintf("ts_bru0_%s_%s.rds", ar_tag, d0_tag)
)

if (!file.exists(model_fname) || override_objects) {
  cat("Fitting ar2 model\n")
  bruar2 <- bru(
    components = components0,
    formula = norm_potential ~ Intercept +
      power_correction +
      d_coast +
      elev +
      wind +
      u,
    family = "gaussian",
    data = wf_df_frag,
    options = bru_options(
      bru_verbose = 3,
      control.inla = list(verbose = TRUE)
    )
  )

  saveRDS(
    bruar2,
    file = model_fname
  )
} else {
  cat("Loading existing ar2 model\n")
  bruar2 <- readRDS(model_fname)
}

summary(bruar2)
# bruar2$summary.fixed[, 1:6]
# bruar2$summary.random$tech_typ[, 1:6]
# bruar2$summary.random$tech_power[, 1:6]

# source("aux_funct.R")
effect_names <- names(bruar2$summary.random)
excluded_effects <- c("u")
effect_names <- setdiff(effect_names, excluded_effects)
for (effect in effect_names) {
  if (effect == "power_correction") {
    n_repl <- 2
    repl_names <- c("Offshore", "Onshore")
  } else {
    n_repl <- 1
    repl_names <- NULL
  }
  plot.effects(
    bruar2,
    effect,
    n.replicate = n_repl,
    replicate_names = repl_names,
    show.plot = TRUE
  )
  ggsave(
    sprintf("fig/%s_effect_%s_%s.pdf", effect, ar_tag, d0_tag),
    width = 6,
    height = 4
  )
}
plot.hyper.dens(bruar2)
ggsave(
  sprintf("fig/hyperparameters_%s_%s.pdf", ar_tag, d0_tag),
  width = 6,
  height = 4
)
# plot(bruar2$summary.fitted.values$mean[1:n], wf_df_frag$norm_potential)
## 2.3 1D SPDE temporal model ####
ar_tag <- "1DSPDE"
mint <- 0
maxt <- 23
buffert <- 0
x <- seq(mint - buffert, maxt + buffert, by = 0.5)
mesh1D <- fm_mesh_1d(
  loc = x,
  interval = c(mint, maxt),
  degree = 2,
  boundary = "cyclic"
)

spde1D <- inla.spde2.pcmatern(
  mesh = mesh1D,
  prior.range = c(3, 0.5), # P(range < 1 hour) = 0.5
  prior.sigma = c(0.03, 0.5) # P(sd > 0.2) = 0.5
)

components0 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  # tech_typ(tech_typ, model = "iid") + # random intercept by tech_typ
  power_correction(
    pow_group,
    model = "rw2",
    # replicate = tech_typ,
    constr = TRUE
  ) + # smooth correction power
  d_coast(
    d_coast_group,
    model = "rw2",
    constr = TRUE
  ) + # smooth correction distance to coast
  elev(
    elev_group,
    model = "rw2",
    constr = TRUE
  ) + # smooth correction elevation
  wind(ws_group, model = "rw2", replicate = tech_typ, constr = TRUE) + # smooth correction wind
  hour(
    t,
    model = spde1D,
    replicate = coord_id
    # hyper = list(prec = list(initial = log(1000), fixed = TRUE))
  )

model_fname <- file.path(
  model_path,
  sprintf("ts_bru0_%s_%s.rds", ar_tag, d0_tag)
)

if (!file.exists(model_fname) || override_objects) {
  cat("Fitting 1D SPDE model\n")
  bru1d <- bru(
    components = components0,
    formula = norm_potential ~ Intercept +
      # tech_typ +
      power_correction +
      d_coast +
      elev +
      wind +
      hour,
    family = "gaussian",
    data = wf_df_frag,
    options = bru_options(
      bru_verbose = 3,
      control.family = list(
        hyper = list(
          prec = list(
            initial = log(30),
            fixed = TRUE
          )
        )
      ),
      control.inla = list(verbose = TRUE)
    )
  )

  saveRDS(
    bru1d,
    file = model_fname
  )
} else {
  cat("Loading existing 1D SPDE model\n")
  bru1d <- readRDS(model_fname)
}

summary(bru1d)
# bru1d$summary.fixed[, 1:6]
# bru1d$summary.random$tech_typ[, 1:6]
# test <- bru1d$summary.random$hour
# plot(bru1d$summary.fitted.values$mean[1:n], wf_df_frag$norm_potential)
# source("aux_funct.R")
effect_names <- names(bru1d$summary.random)
excluded_effects <- c("u", "hour")
effect_names <- setdiff(effect_names, excluded_effects)
for (effect in effect_names) {
  if (effect == "power_correction") {
    n_repl <- 2
    repl_names <- c("Offshore", "Onshore")
  } else {
    n_repl <- 1
    repl_names <- NULL
  }
  plot.effects(
    bru1d,
    effect,
    n.replicate = n_repl,
    replicate_names = repl_names,
    show.plot = TRUE
  )
  ggsave(
    sprintf("fig/%s_effect_%s_%s.pdf", effect, ar_tag, d0_tag),
    width = 6,
    height = 4
  )
}

plot.hyper.dens(bru1d)
ggsave(
  sprintf("fig/hyperparameters_%s_%s.pdf", ar_tag, d0_tag),
  width = 6,
  height = 4
)


## 2.4 ST SPDE model ####
spde <- INLA::inla.spde2.pcmatern(
  mesh = wf.mesh,
  prior.range = c(50, 0.5), # P(range < 100 km)=0.5
  prior.sigma = c(0.2, 0.5) # P(sd > 0.2)=0.5
)

components0 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  # tech_typ(tech_typ, model = "iid") + # random intercept by tech_typ
  power_correction(
    pow_group,
    model = "rw2",
    # replicate = tech_typ,
    constr = TRUE
  ) + # smooth correction power
  d_coast(
    d_coast_group,
    model = "rw2",
    constr = TRUE
  ) + # smooth correction distance to coast
  elev(
    elev_group,
    model = "rw2",
    constr = TRUE
  ) + # smooth correction elevation
  wind(ws_group, model = "rw2", replicate = tech_typ, constr = TRUE) + # smooth correction wind
  st_field(
    geometry,
    model = spde,
    group = time_id,
    control.group = list(model = "ar1")
  )

model_fname <- file.path(
  model_path,
  sprintf("st_bru0_%s_mesh_%s.rds", mesh_label, d0_tag)
)

if (!file.exists(model_fname) || override_objects) {
  cat("Fitting spatiotemporal model\n")
  bru0 <- bru(
    components = components0,
    formula = norm_potential ~ Intercept +
      power_correction +
      d_coast +
      elev +
      wind +
      st_field,
    family = "gaussian",
    data = wf_df_frag,
    options = bru_options(
      bru_verbose = 3,
      control.inla = list(verbose = TRUE)
    )
  )

  saveRDS(
    bru0,
    file = model_fname
  )
} else {
  cat("Loading existing spatiotemporalmodel\n")
  bru0 <- readRDS(model_fname)
}

## summary and effect plots ####
summary(bru0)
# bru0$summary.fixed[, 1:6]
# bru0$summary.random$tech_typ[, 1:6]
# bru0$summary.random$tech_power[, 1:6]

# source("aux_funct.R")
effect_names <- names(bru0$summary.random)
excluded_effects <- c("u", "hour", "st_field")
effect_names <- setdiff(effect_names, excluded_effects)
for (effect in effect_names) {
  if (effect == "power_correction") {
    n_repl <- 2
    repl_names <- c("Offshore", "Onshore")
  } else {
    n_repl <- 1
    repl_names <- NULL
  }
  plot.effects(
    bru0,
    effect,
    n.replicate = n_repl,
    replicate_names = repl_names,
    show.plot = TRUE
  )
  ggsave(
    sprintf("fig/%s_effect_%s_%s.pdf", effect, mesh_label, d0_tag),
    width = 6,
    height = 4
  )
}
plot.hyper.dens(bru0)
ggsave(
  sprintf("fig/hyperparameters_%s_%s.pdf", mesh_label, d0_tag),
  width = 6,
  height = 4
)
# wf_df_frag %>%
#   pull(pow_group) %>%
#   range()

# wf_df_frag %>%
#   ggplot()+ geom_density(aes(pow_group, fill = tech_typ), alpha = 0.5)+theme_minimal()

### plot intensity of spatial field ####

ppxl <- fm_pixels(wf.mesh, mask = bnd[[2]], format = "sf", dims = pixel_dims)
ppxl_all <- fm_cprod(
  ppxl,
  data.frame(
    # time_id = unique(wf_df_frag$time_id)
    time_id = c(9, 12, 18)
  )
)

set.seed(1)
safe_predict <- function(model, newdata, fun, n1 = 100, n2 = 10) {
  tryCatch(
    {
      fun(model, newdata, n.samples = n1)
    },
    error = function(e) {
      message("First predict failed: ", conditionMessage(e))
      message("Retrying with n.samples = ", n2)

      tryCatch(
        {
          fun(model, newdata, n.samples = n2)
        },
        error = function(e2) {
          message("Second predict also failed: ", conditionMessage(e2))
          stop(e2)
        }
      )
    }
  )
}
pow_est_st <- safe_predict(
  model = bru0,
  newdata = ppxl_all,
  fun = function(model, newdata, n.samples) {
    predict(
      model,
      newdata,
      ~ data.frame(
        time_id = time_id,
        norm_potential_est = pmax(-1, pmin(1, st_field)) # should i cap this?
      ),
      n.samples = n.samples
    )
  },
  n1 = 100,
  n2 = 10
)


p_median <- ggplot() +
  gg(pow_est_st, geom = "tile", aes(fill = q0.5)) +
  geom_sf(data = uk_map, fill = NA, color = "black", alpha = 0.5) +
  # gg(wf.mesh, alpha = 0.5) +
  geom_point(data = loc_unique, aes(x, y), color = "darkred", size = 0.5) +
  facet_wrap(
    ~time_id,
    labeller = as_labeller(c("9" = "9:00", "12" = "12:00", "18" = "18:00"))
  ) +
  coord_sf() +
  scale_fill_viridis_c() +
  theme_void()
p_median
ggsave(
  sprintf("fig/%s_spatial_field_median_%s.pdf", mesh_label, d0_tag),
  width = 10,
  height = 6
)

p_sd <- ggplot() +
  gg(pow_est_st, geom = "tile", aes(fill = sd)) +
  geom_sf(data = uk_map, fill = NA, color = "white", alpha = 0.5) +
  # gg(wf.mesh, alpha = 0.5) +
  geom_point(data = loc_unique, aes(x, y), color = "darkred", size = 0.5) +
  facet_wrap(
    ~time_id,
    labeller = as_labeller(c("9" = "9:00", "12" = "12:00", "18" = "18:00"))
  ) +
  coord_sf() +
  scale_fill_viridis_c(option = "inferno") +
  theme_void()
p_sd
ggsave(
  sprintf("fig/%s_spatial_field_sd_%s.pdf", mesh_label, d0_tag),
  width = 10,
  height = 6
)
# plot ts
# mesh triangle size
# different days
# covariate for terrain if necessary
#

# ppxl_all <- fm_cprod(
#   ppxl,
#   data.frame(
#     time_id = seq(9, 12, 18),
#     tech_typ = unique(wf_df_frag$tech_typ),
#     norm_power_est0 = seq(0, 1, length.out = 10),
#     ws_group = unique(wf_df_frag$ws_group) %>% sort()
#   )
# )

# pow_est_st <- predict(
#   bru0,
#   ppxl_all,
#   ~ data.frame(
#     time_id = time_id,
#     norm_potential_est = (Intercept +
#       power_correction +
#       tech_typ +
#       tech_power +
#       wind +
#       st_field)
#   )
# )

## 2.3 lm wf version ####

lm_wf_modfname <- sprintf("lm_model_aic0_%s.rds", d0_tag)

if (!file.exists(file.path(model_path, lm_wf_modfname)) || override_objects) {
  cat("Fitting LM model\n")
  base_model <- lm(
    norm_potential ~ norm_power_est0,
    data = wf_df_frag
  )

  full_model0 <- lm(
    norm_potential ~
      tech_typ *
      norm_power_est0 +
      # norm_power_est0 * month +
      # hour +
      # dist_coast * tech_typ +
      # elevation * tech_typ +
      # dist_coast:tech_typ +
      # elevation:tech_typ +
      tech_typ * poly(dist_coast, 2) +
      tech_typ * poly(elevation, 3) +
      tech_typ * poly(ws_h_wmean, 3),
    data = wf_df_frag
  )

  model_AIC0 <- step(
    base_model,
    scope = list(lower = base_model, upper = full_model0),
    # steps = 5,
    k = 2
  )

  saveRDS(
    model_AIC0,
    file.path(model_path, lm_wf_modfname)
  )
} else {
  cat("Loading existing LM model\n")
  model_AIC0 <- readRDS(file.path(model_path, lm_wf_modfname))
}

## 2.4 GB lm version #####

lm_agg_modfname <- sprintf("lm_model_aic0_agg_%s.rds", d0_tag)

if (!file.exists(file.path(model_path, lm_agg_modfname)) || override_objects) {
  cat("Fitting GB aggregated LM\n")

  samp_gb <- GB_df
  # %>%
  # filter(date %in% sampled_days) %>%

  base_model_agg <- lm(
    norm_potential ~ norm_power_est0,
    data = samp_gb
  )

  full_model0_agg <- lm(
    norm_potential ~
      tech_typ *
      norm_power_est0 +
      # norm_power_est0 * month +
      # hour +
      # dist_coast * tech_typ +
      # elevation * tech_typ +
      # dist_coast:tech_typ +
      # elevation:tech_typ +
      # tech_typ * poly(dist_coast, 2) +
      # tech_typ * poly(elevation, 3) +
      tech_typ * poly(ws_h_wmean, 3),
    data = samp_gb
  )
  summary(full_model0_agg)

  model_AIC0_agg <- step(
    base_model_agg,
    scope = list(lower = base_model_agg, upper = full_model0_agg),
    # steps = 5,
    k = 2
  )
  saveRDS(
    model_AIC0_agg,
    file.path(model_path, lm_agg_modfname)
  )
} else {
  cat("Loading existing GB aggregated LM\n")
  model_AIC0_agg <- readRDS(file.path(
    model_path,
    lm_agg_modfname
  ))
}

## 2.5 QM version ####

qm_fname <- sprintf("qm_model_%s.rds", d0_tag)

if (!file.exists(file.path(model_path, qm_fname)) || override_objects) {
  cat("Fitting quantile mapping model\n")
  qqmod <- fitQmap(
    obs = wf_df_frag %>% pull(norm_potential),
    mod = wf_df_frag %>% pull(norm_power_est0),
    method = "QUANT"
  )
  saveRDS(
    qqmod,
    file.path(model_path, qm_fname)
  )
} else {
  cat("Loading existing quantile mapping model\n")
  qqmod <- readRDS(file.path(model_path, qm_fname))
}

# 3. model comparison ####

wgen_qm <- with(
  wf_df_frag,
  doQmapQUANT(norm_power_est0, qqmod, type = "linear")
)

mod_labels <- c(
  "Generic PC",
  "Linear model",
  "AR1 model",
  "AR2 model",
  # "1D SPDE model",
  "Spatio-temporal model",
  "QM",
  "GB LM"
)
est_cols <- c(
  "norm_power_est0",
  "lm",
  "ar1",
  "ar2",
  # "spde1d",
  "st",
  "qm",
  "agg_lm"
)
n_models <- length(est_cols)
n <- nrow(wf_df_frag)
names(mod_labels) <- est_cols
# length(wgen_qm)
# length(wf_df_frag$norm_potential)
# length(bru0$summary.fitted.values[1:n, "mean"])
# length(model_AIC0$fitted.values)
# n <- nrow(wf_df_frag)
# names(mod_labels) <- est_cols

# ar_test <- predict(
#   bruar1,
#   newdata = wf_df_frag,
#   # ~ data.frame(
#   #   time_id = time_id,
#   #   # latent = x,
#   #   norm_potential_est = (Intercept +
#   #     power_correction +
#   #     wind)
#   # ),
#   ~ Intercept +
#     power_correction +
#     wind +
#     u,
#   n.samples = 10
# )
# length(ar_test$mean)
# names(ar_test)
# plot(ar_test$mean, wf_df_frag$norm_potential)

# set.seed(1)
# spde1d_pred <- predict(
#   bru1d,
#   newdata = wf_df_frag,
#   ~ Intercept +
#     power_correction +
#     wind +
#     hour,
#   type = "response",
#   n.samples = 100,
#   verbose = TRUE,
#   used.improved.mean = FALSE
# )

## fitted values df ####
model_df0 <- wf_df_frag %>%
  mutate(
    date = as.Date(time),
    lm = predict(model_AIC0, newdata = .),
    ar1 = bruar1$summary.fitted.values[1:n, "mean"],
    ar2 = bruar2$summary.fitted.values[1:n, "mean"],
    spde1d = bru1d$summary.fitted.values[1:n, "mean"],
    st = bru0$summary.fitted.values[1:n, "mean"],
    qm = wgen_qm,
    agg_lm = predict(model_AIC0_agg, newdata = wf_df_frag)
  ) %>%
  mutate(
    st_low = bru0$summary.fitted.values[1:n, "0.025quant"],
    st_high = bru0$summary.fitted.values[1:n, "0.975quant"]
  ) %>%
  left_join(
    gb_day_df %>% dplyr::select(date, p_group3),
    by = c("date" = "date")
  )

st_write(
  model_df0,
  sprintf("data/calibration_df_%s_%s.gpkg", mesh_label, d0_tag),
  driver = "GPKG",
  append = FALSE,
)

## 3.1 Exploring fitted values ####

# mod_labels <- c(
#   "Generic PC",
#   "Linear model",
#   "Spatio-temporal model",
#   "QM",
#   "GB LM"
# )
# est_cols <- c(
#   "norm_power_est0",
#   "lm",
#   "st",
#   "qm",
#   "agg_lm"
# )
# n_models <- length(est_cols)
# n <- nrow(wf_df_frag)
# names(mod_labels) <- est_cols
model_df0 <- st_read(sprintf(
  "data/calibration_df_%s_%s.gpkg",
  mesh_label,
  d0_tag
))

pos_breaks <- with(
  model_df0,
  quantile(elevation[elevation > 0], probs = seq(0, 1, 1 / 3))
)
pos_levels <- levels(cut(
  model_df0$elevation[model_df0$elevation > 0],
  breaks = pos_breaks,
  include.lowest = TRUE
))

## summary table for figures and tables #####
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
    model = factor(model, levels = est_cols, labels = mod_labels)
  )

df_long0 %>%
  ggplot() +
  geom_density(
    aes(x = err, fill = model),
    alpha = 0.5
  ) +
  facet_wrap(~model, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Error distribution by model",
    x = "Error (model estimate - observed)",
    y = "Density",
    fill = "Model"
  ) +
  # coord_cartesian(xlim = c(-0.25, 0.25)) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = pal_lancet()(n_models))

ggsave(
  sprintf("fig/error_distribution_by_model_%s_%s.pdf", mesh_label, d0_tag),
  width = 6,
  height = 4
)

metrics_table <- df_long0 %>%
  st_drop_geometry() %>%
  group_by(model) %>%
  summarise(
    RMSE = ModelMetrics::rmse(actual = norm_potential, predicted = estimate),
    MAE = ModelMetrics::mae(actual = norm_potential, predicted = estimate),
    # MAPE = mape(actual = norm_potential, predicted = estimate, pos_only = TRUE),
    MDAPE = mdape(
      actual = norm_potential,
      predicted = estimate,
      pos_only = TRUE
    ),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    model = mod_labels[model]
  ) %>%
  arrange(desc(RMSE))
metrics_table
write.csv(
  metrics_table,
  sprintf("summaries/calib_metrics_%s_%s.csv", mesh_label, d0_tag),
  row.names = FALSE
)

## 3.2 aggregated time series version ####
metrics_table <- df_long0 %>%
  st_drop_geometry() %>%
  group_by(time, model) %>%
  summarise(
    across(
      c(norm_potential, estimate),
      ~ sum(. * capacity) / sum(capacity)
    )
  ) %>%
  group_by(model) %>%
  summarise(
    RMSE = ModelMetrics::rmse(actual = norm_potential, predicted = estimate),
    MAE = ModelMetrics::mae(actual = norm_potential, predicted = estimate),
    MAPE = mape(actual = norm_potential, predicted = estimate, pos_only = TRUE),
    MDAPE = mdape(
      actual = norm_potential,
      predicted = estimate,
      pos_only = TRUE
    ),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    model = mod_labels[model]
  ) %>%
  arrange(desc(RMSE))
metrics_table
write.csv(
  metrics_table,
  sprintf("summaries/calib_metrics_%s_%s_gb.csv", mesh_label, d0_tag),
  row.names = FALSE
)

df_long0 %>%
  group_by(time, model) %>%
  summarise(
    across(
      c(norm_potential, estimate),
      ~ sum(. * capacity) / sum(capacity)
    )
  ) %>%
  mutate(
    err = estimate - norm_potential,
    # model = factor(model, levels = est_cols, labels = mod_labels)
  ) %>%
  group_by(model) %>%
  ggplot() +
  geom_density(
    aes(x = err, fill = model),
    alpha = 0.5
  ) +
  facet_wrap(~model, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Error distribution by model",
    x = "Error (model estimate - observed)",
    y = "Density",
    fill = "Model"
  ) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = pal_lancet()(n_models))

ggsave(
  sprintf("fig/error_distribution_by_model_%s_%s_gb.pdf", mesh_label, d0_tag),
  width = 6,
  height = 4
)

## 3.3 error distribution by covariates ####
### 3.3.1 by tech type #####
metrics_table <- df_long0 %>%
  st_drop_geometry() %>%
  group_by(tech_typ, model) %>%
  summarise(
    RMSE = ModelMetrics::rmse(actual = norm_potential, predicted = estimate),
    MAE = ModelMetrics::mae(actual = norm_potential, predicted = estimate),
    MAPE = mape(actual = norm_potential, predicted = estimate, pos_only = TRUE),
    MDAPE = mdape(
      actual = norm_potential,
      predicted = estimate,
      pos_only = TRUE
    ),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # mutate(model = mod_labels[model]) %>%
  arrange(tech_typ, desc(RMSE))
metrics_table
write.csv(
  metrics_table,
  sprintf("summaries/calib_metrics_%s_%s_tech.csv", mesh_label, d0_tag),
  row.names = FALSE
)

df_long0 %>%
  st_drop_geometry() %>%
  ggplot() +
  geom_density_ridges(
    aes(err, tech_typ, fill = model),
    alpha = 0.5,
    scale = 1
  ) +
  theme_ridges() +
  labs(
    title = "Error distribution by technology type",
    x = "Error (model estimate - observed)",
    y = "Technology Type",
    fill = "Model"
  ) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = pal_lancet()(n_models))

ggsave(
  sprintf("fig/error_distribution_by_tech_type_%s_%s.pdf", mesh_label, d0_tag),
  width = 6,
  height = 4
)

### 3.3.2 by regime #####
metrics_table <- df_long0 %>%
  st_drop_geometry() %>%
  group_by(p_group3, model) %>%
  summarise(
    RMSE = ModelMetrics::rmse(actual = norm_potential, predicted = estimate),
    MAE = ModelMetrics::mae(actual = norm_potential, predicted = estimate),
    MAPE = mape(actual = norm_potential, predicted = estimate, pos_only = TRUE),
    MDAPE = mdape(
      actual = norm_potential,
      predicted = estimate,
      pos_only = TRUE
    ),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(p_group3, desc(RMSE))
metrics_table

write.csv(
  metrics_table,
  sprintf("summaries/calib_metrics_%s_%s_regime.csv", mesh_label, d0_tag),
  row.names = FALSE
)

df_long0 %>%
  st_drop_geometry() %>%
  ggplot() +
  geom_density_ridges(
    aes(err, p_group3, fill = model),
    alpha = 0.5,
    scale = 1
  ) +
  theme_ridges() +
  labs(
    title = "Error distribution by regime",
    x = "Error (model estimate - observed)",
    y = "Regime",
    fill = "Model"
  ) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = pal_lancet()(n_models))
ggsave(
  sprintf("fig/error_distribution_by_regime_%s_%s.pdf", mesh_label, d0_tag),
  width = 6,
  height = 4
)


#### 3.3.3 by hour of day #####
df_long0 %>%
  mutate(hour = factor(hour, levels = 0:23)) %>%
  st_drop_geometry() %>%
  ggplot() +
  geom_density_ridges(
    aes(x = err, y = hour, fill = model),
    alpha = 0.5
  ) +
  facet_wrap(~model) +
  theme_ridges() +
  labs(
    title = "Error distribution by hour of day",
    x = "Error (model estimate - observed)",
    y = "Hour of day",
    fill = "Model"
  ) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(-0.5, 0.5)) +
  scale_fill_manual(values = pal_lancet()(n_models))

ggsave(
  sprintf("fig/error_distribution_by_hourA_%s_%s.pdf", mesh_label, d0_tag),
  width = 12,
  height = 8
)

df_long0 %>%
  st_drop_geometry() %>%
  group_by(hour, model) %>%
  summarise(
    RMSE = ModelMetrics::rmse(actual = norm_potential, predicted = estimate),
    MAE = ModelMetrics::mae(actual = norm_potential, predicted = estimate),
    MAPE = mape(actual = norm_potential, predicted = estimate, pos_only = TRUE),
    MDAPE = mdape(
      actual = norm_potential,
      predicted = estimate,
      pos_only = TRUE
    ),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(model = mod_labels[model]) %>%
  arrange(hour, desc(RMSE)) %>%
  ggplot(aes(hour, RMSE, col = model, group = model)) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Model performance by hour of day",
    x = "Hour of day",
    y = "RMSE",
    col = "Model"
  ) +
  scale_color_manual(values = pal_lancet()(n_models))

ggsave(
  sprintf("fig/model_performance_by_hourS_%s_%s.pdf", mesh_label, d0_tag),
  width = 6,
  height = 4
)


#### 3.3.4 by distance to coast #####
metrics_table <- df_long0 %>%
  st_drop_geometry() %>%
  group_by(dist_coast_g4, model) %>%
  summarise(
    RMSE = ModelMetrics::rmse(actual = norm_potential, predicted = estimate),
    MAE = ModelMetrics::mae(actual = norm_potential, predicted = estimate),
    MAPE = mape(actual = norm_potential, predicted = estimate, pos_only = TRUE),
    MDAPE = mdape(
      actual = norm_potential,
      predicted = estimate,
      pos_only = TRUE
    ),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(dist_coast_g4, desc(RMSE))
metrics_table

write.csv(
  metrics_table,
  sprintf("summaries/calib_metrics_%s_%s_dist_coast.csv", mesh_label, d0_tag),
  row.names = FALSE
)

df_long0 %>%
  st_drop_geometry() %>%
  ggplot() +
  geom_density_ridges(
    aes(err, dist_coast_g4, fill = model),
    alpha = 0.5,
    scale = 1
  ) +
  theme_ridges() +
  labs(
    title = "Error distribution by distance to coast",
    x = "Error (model estimate - observed)",
    y = "Distance to Coast",
    fill = "Model"
  ) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = pal_lancet()(n_models))
ggsave(
  sprintf("fig/error_distribution_by_dist_coast_%s_%s.pdf", mesh_label, d0_tag),
  width = 6,
  height = 4
)


#### 3.3.5 by elevation #####
metrics_table <- df_long0 %>%
  st_drop_geometry() %>%
  group_by(elevation_g4, model) %>%
  summarise(
    RMSE = ModelMetrics::rmse(actual = norm_potential, predicted = estimate),
    MAE = ModelMetrics::mae(actual = norm_potential, predicted = estimate),
    MAPE = mape(actual = norm_potential, predicted = estimate, pos_only = TRUE),
    MDAPE = mdape(
      actual = norm_potential,
      predicted = estimate,
      pos_only = TRUE
    ),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(elevation_g4, desc(RMSE))
metrics_table

write.csv(
  metrics_table,
  sprintf("summaries/calib_metrics_%s_%s_elevation.csv", mesh_label, d0_tag),
  row.names = FALSE
)

df_long0 %>%
  st_drop_geometry() %>%
  ggplot() +
  geom_density_ridges(
    aes(err, elevation_g4, fill = model),
    alpha = 0.5,
    scale = 1
  ) +
  theme_ridges() +
  labs(
    title = "Error distribution by elevation",
    x = "Error (model estimate - observed)",
    y = "Elevation",
    fill = "Model"
  ) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = pal_lancet()(n_models))
ggsave(
  sprintf("fig/error_distribution_by_elevation_%s_%s.pdf", mesh_label, d0_tag),
  width = 6,
  height = 4
)


## 3.4 time series plots #####

### 3.4.1 aggregated time series version #####
mod_labels2 <- c(mod_labels, "norm_potential" = "Observed")
model_df_ts <- model_df0 %>%
  dplyr::select(
    time,
    site_name,
    tech_typ,
    norm_potential,
    capacity,
    ws_h,
    all_of(est_cols)
  ) %>%
  pivot_longer(
    cols = all_of(c("norm_potential", est_cols)),
    names_to = "model",
    values_to = "estimate"
  ) %>%
  mutate(
    estimate = pmin(1, pmax(0, estimate)), # clipping estimates to [0, 1]
  )

model_df_ts %>%
  ggplot() +
  geom_line(
    aes(time, estimate, group = site_name),
    alpha = 0.5,
    col = "gray50"
  ) +
  geom_line(
    data = model_df_ts %>%
      group_by(time, model) %>%
      summarise(
        power = sum(estimate * capacity) / sum(capacity),
        .groups = "drop"
      ),
    aes(time, power, col = "capacity weighted avg."),
    lwd = 1
  ) +
  geom_line(
    data = model_df_ts %>%
      group_by(time, model) %>%
      summarise(power = mean(estimate), .groups = "drop"),
    aes(time, power, col = "simple avg."),
    lwd = 1
  ) +
  theme_minimal() +
  facet_wrap(~model, ncol = 3, labeller = as_labeller(mod_labels2)) +
  scale_x_datetime(date_labels = "%H:%M") +
  labs(
    title = sprintf("Power estimates Time Series %s", d0),
    x = "Time",
    y = "Generation (% of capacity)",
    col = ""
  ) +
  scale_color_manual(
    values = c("capacity weighted avg." = "blue", "simple avg." = "red")
  ) +
  theme(legend.position = "bottom")
ggsave(
  sprintf("fig/power_estimates_time_series_%s_%s.pdf", mesh_label, d0_tag),
  width = 10,
  height = 6
)

### 3.4.2 some locations time series ####

set.seed(1)
sample_sites <- model_df0 %>%
  group_by(tech_typ) %>%
  distinct(site_name) %>%
  slice_sample(n = 2) %>%
  pull(site_name)
sample_sites

model_df_ts %>%
  filter(site_name %in% sample_sites) %>%
  mutate(model = mod_labels2[model]) %>%
  ggplot() +
  geom_line(
    aes(time, estimate, group = model, col = model),
    # alpha = 0.5
  ) +
  theme_minimal() +
  facet_wrap(
    ~ sprintf(
      "%s (%s)",
      site_name,
      gsub("Wind", "", tech_typ) %>% trimws()
    ),
    ncol = 2
  ) +
  scale_x_datetime(date_labels = "%H:%M") +
  labs(
    title = sprintf("Power estimates Time Series %s", d0),
    x = "Time",
    y = "Generation (% of capacity)",
    col = ""
  ) +
  # scale_fill_manual(values = pal_lancet()(n_models)) +
  theme(legend.position = "bottom")
ggsave(
  sprintf(
    "fig/power_estimates_time_series_%s_sampWF_%s.pdf",
    mesh_label,
    d0_tag
  ),
  width = 10,
  height = 6
)


### 3.4.3 aggregated error time series version ####
model_df_ts2 <- model_df0 %>%
  dplyr::select(
    time,
    site_name,
    tech_typ,
    norm_potential,
    capacity,
    ws_h,
    all_of(est_cols)
  ) %>%
  mutate(across(
    all_of(est_cols),
    ~ . - norm_potential
  )) %>%
  pivot_longer(
    cols = all_of(est_cols),
    names_to = "model",
    values_to = "error"
  )
model_df_ts2 %>%
  ggplot() +
  geom_line(
    aes(time, error, group = site_name),
    alpha = 0.5,
    col = "gray50"
  ) +
  geom_line(
    data = model_df_ts2 %>%
      group_by(time, model) %>%
      summarise(
        error = sum(error * capacity) / sum(capacity),
        .groups = "drop"
      ),
    aes(time, error, col = "capacity weighted avg."),
    lwd = 1
  ) +
  geom_line(
    data = model_df_ts2 %>%
      group_by(time, model) %>%
      summarise(error = mean(error), .groups = "drop"),
    aes(time, error, col = "simple avg."),
    lwd = 1
  ) +
  theme_minimal() +
  facet_wrap(~model, ncol = 3, labeller = as_labeller(mod_labels)) +
  scale_x_datetime(date_labels = "%H:%M") +
  labs(
    title = sprintf("Estimation Error Time Series %s", d0),
    x = "Time",
    y = "Error (% of capacity)",
    col = ""
  ) +
  scale_color_manual(
    values = c("capacity weighted avg." = "blue", "simple avg." = "red")
  ) +
  theme(legend.position = "bottom")
ggsave(
  sprintf(
    "fig/power_estimates_error_time_series_%s_%s.pdf",
    mesh_label,
    d0_tag
  ),
  width = 10,
  height = 6
)
