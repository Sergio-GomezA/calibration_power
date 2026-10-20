# set up

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

pow_threshold <- 0.05
pow_threshold_label <- gsub("\\.", "_", as.character(pow_threshold))
tables_path <- "~/ownCloud-s2441782@datasync.ed.ac.uk/projects/calibration/calibration_power_main_doc/tables"
dir.exists(tables_path) || dir.create(tables_path, recursive = TRUE)

override_objects <- FALSE
batch_name <- "bias_eda"
norm_dist_tol <- 1

if (local_run) {
  cat("Running in local mode\n")
} else {
  cat("Running in cluster mode\n")
}

# 0.2 libraries and paths ####
require(parallel)

if (local_run) {
  data_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
  extra_path <- "~/Documents/elexon/extra"
  output_path <- "~/Documents/elexon/caloutput"
  sample_path <- "~/Documents/elexon/samples"
  n_samp <- 100
  pixel_dims <- c(150, 150)
} else {
  data_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  gen_path <- "/exports/eddie/scratch/s2441782/calibration/data"
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

# Data

extension <- "rds"
coord_list_fname <- "data/coord_list.csv"

cat("Loading existing coordinate list\n")
coord_list <- read.csv(coord_list_fname)


gb_day_df_fname <- sprintf("data/GB_daily_summary.parquet")

if (!file.exists(gb_day_df_fname) || override_objects) {
  if (!file.exists(gb_day_df_fname)) {
    cat("GB daily summary file not found, creating new summary\n")
  } else {
    cat(
      "GB daily summary file found, but override_objects is TRUE. Recreating summary\n"
    )
  }
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

gb_df_fname <- sprintf("GB_aggr_An-%0.1f.parquet", norm_dist_tol)
GB_df <- read_parquet(file.path(gen_path, gb_df_fname)) %>%
  # rename(time = halfHourEndTime) %>%
  mutate(
    err = norm_power_est0 - norm_potential,
    error0 = norm_potential - norm_power_est0,
    date = as.Date(time)
  ) %>%
  mutate(
    err_l1 = lag(error0, default = 0),
    err_l2 = lag(error0, n = 2, default = 0),
    dest = (norm_power_est0) - lag(norm_power_est0, default = 0),
    dest_l1 = lag(dest, default = 0),
    dest_l2 = lag(dest, n = 2, default = 0)
  )
pacf_df <- GB_df %>%
  arrange(tech_typ, time) %>%
  group_by(tech_typ) %>%
  summarise(
    pacf = list(
      pacf(error0, plot = FALSE, na.action = na.pass)
    ),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lag = map(pacf, ~ .x$lag),
    acf = map(pacf, ~ .x$acf),
    ci = 1.96 / sqrt(n)
  ) %>%
  dplyr::select(-pacf) %>%
  unnest(c(lag, acf))

# head(pacf_df)
pacf_df %>%
  ggplot(aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_segment(aes(xend = lag, yend = 0)) +
  facet_wrap(~tech_typ, scales = "free_y") +
  labs(
    x = "Lag",
    y = "Partial autocorrelation"
  ) +
  theme_bw() +
  coord_cartesian(ylim = c(-0.4, 1)) +
  scale_y_continuous(breaks = seq(-0.3, 1, by = 0.2))
ggsave(
  sprintf("fig/%s/GB_pacf.pdf", batch_name),
  width = 10,
  height = 6,
  create.dir = TRUE
)

# LM Model

base_model_agg <- lm(
  norm_potential ~ norm_power_est0,
  data = GB_df
)

full_model0_agg <- lm(
  norm_potential ~
    tech_typ +
    tech_typ *
      norm_power_est0 +
    month * tech_typ +
    hour * tech_typ +
    # err_l1 +
    # err_l2 +
    dest +
    dest_l1 +
    dest_l2 +
    # norm_power_est0 * month +
    # hour +
    # dist_coast * tech_typ +
    # elevation * tech_typ +
    # dist_coast:tech_typ +
    # elevation:tech_typ +
    # tech_typ * poly(dist_coast, 2) +
    # tech_typ * poly(elevation, 3) +
    tech_typ * poly(ws_h_wmean, 3),
  data = GB_df
)
summary(full_model0_agg)

model_AIC0_agg <- step(
  base_model_agg,
  scope = list(lower = base_model_agg, upper = full_model0_agg),
  # steps = 5,
  k = 2,
  trace = 0
)
summary(model_AIC0_agg)
# Residuals
GB_df <- GB_df %>%
  mutate(
    res_AIC0_agg = resid(model_AIC0_agg),
    res_AIC0_agg_std = scale(res_AIC0_agg)
  )
# AR PACF
pacf_df <- GB_df %>%
  arrange(tech_typ, time) %>%
  group_by(tech_typ) %>%
  summarise(
    pacf = list(
      pacf(res_AIC0_agg, plot = FALSE, na.action = na.pass)
    ),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lag = map(pacf, ~ .x$lag),
    acf = map(pacf, ~ .x$acf),
    ci = 1.96 / sqrt(n)
  ) %>%
  dplyr::select(-pacf) %>%
  unnest(c(lag, acf))

pacf_df %>%
  ggplot(aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_segment(aes(xend = lag, yend = 0)) +
  facet_wrap(~tech_typ, scales = "free_y") +
  labs(
    x = "Lag",
    y = "Partial autocorrelation"
  ) +
  theme_bw() +
  coord_cartesian(ylim = c(-0.4, 1)) +
  scale_y_continuous(breaks = seq(-0.3, 1, by = 0.2))
ggsave(
  sprintf("fig/%s/GB_pacf.pdf", batch_name),
  width = 10,
  height = 6,
  create.dir = TRUE
)

# fitting AR
# time index creation
GB_df <- GB_df %>%
  arrange(tech_typ, time) %>%
  group_by(tech_typ) %>%
  mutate(
    time_index = as.numeric(
      difftime(
        time,
        min(time, na.rm = TRUE),
        units = "hours"
      )
    )
  ) %>%
  ungroup()

missing_hours <- GB_df %>%
  filter(tech_typ == "Wind Onshore") %>%
  distinct(time_index) %>%
  anti_join(
    tibble(
      time_index = seq(
        min(.$time_index, na.rm = TRUE),
        max(.$time_index, na.rm = TRUE),
        by = 1
      )
    ),
    by = "time_index"
  ) %>%
  arrange(time_index)

missing_hours

all_hours <- with(
  GB_df,
  seq(
    min(time_index, na.rm = TRUE),
    max(time_index, na.rm = TRUE),
    by = 1
  )
)
missing_time_index <- which(!all_hours %in% GB_df$time_index)

GB_df <- GB_df %>%
  complete(
    time_index = all_hours,
    tech_typ,
    fill = list(
      norm_power_est0 = NA,
      norm_potential = NA,
      ws_h_wmean = NA,
      capacity = NA,
      err = NA,
      error0 = NA,
      date = as.Date(NA),
      err_l1 = NA,
      err_l2 = NA,
      dest = NA,
      dest_l1 = NA,
      dest_l2 = NA
    )
  )


mf <- model.frame(
  formula(model_AIC0_agg),
  data = GB_df,
  na.action = na.omit
)
X_all <- model.matrix(
  terms(model_AIC0_agg),
  data = mf
)
y_all <- model.response(mf)
idx <- mf$tech_typ == "Wind Onshore"


dat_on <- GB_df[
  match(rownames(mf)[idx], rownames(GB_df)),
]
y <- y_all[idx]
X <- X_all[idx, , drop = FALSE]

# Sort chronologically
ord <- order(dat_on$time_index)
dat_on <- dat_on[ord, ]
y <- y[ord]
X <- X[ord, , drop = FALSE]
stopifnot(!anyDuplicated(dat_on$time_index))

table(diff(dat_on$time_index))

segment <- cumsum(c(
  TRUE,
  diff(dat_on$time_index) != 1
))

segment_sizes <- table(segment)
longest_segment <- as.integer(names(segment_sizes)[
  which.max(segment_sizes)
])

# inla model
mod_starttime <- Sys.time()
prec_init <- log(200) # for u
ar_tag <- "ar1"
components0 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  techno(tech_typ, model = "iid") + # random intercept by tech_typ
  slope(
    tech_typ,
    model = "iid",
    weights = norm_power_est0
  ) +
  wind(ws_group, model = "rw2", replicate = tech_typ, constr = TRUE) + # smooth correction wind
  u(
    t,
    model = "ar1",
    replicate = tech_typ,
    hyper = list(
      rho = list(
        prior = "pc.cor1",
        param = c(0.6, 0.5)
      ),
      # prec = list(
      #   prior = "pc.prec",
      #   param = c(50, 0.05)
      # ),
      prec = list(initial = prec_init, fixed = FALSE)
    )
  )


model_code <- sprintf("ts_agg0_%s_%s.rds", ar_tag, "full")
model_fname <- file.path(
  model_path,
  model_code
)


cat(
  "-------------------------------------------------------------------------------------------------\n"
)
cat("Fitting ar1 model\n")
cat(
  "-------------------------------------------------------------------------------------------------\n"
)
bruar1 <- bru(
  components = components0,
  formula = norm_potential ~ Intercept +
    techno +
    slope +
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
  data = GB_df %>% filter(time >= "2024-01-01"),
  options = base_bru_options
)

summary(bruar1)
# Feature engineering
# Events

# Cluster
