# Power curve model ####

require(arrow)
require(dplyr)
require(tidyr)
require(rnaturalearth)
require(rnaturalearthdata)
require(sf)
require(ggplot2)
require(ggthemes)
require(ggsci)
# require(FNN)
require(data.table)
require(parallel)
library(purrr)
require(brms)
require(INLA)

require(inlabru)

source("aux_funct.R")


# read Data ####
if (grepl("exports", getwd())) {
  # running in cluster
  cluster_run <- TRUE
  data_path <- gen_path <- "data"
  model_path <- "model_objects"
} else {
  cluster_run <- FALSE
  data_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
}

source("aux_funct.R")


if (cluster_run) {
  inla.setOption(num.threads = paste0("0:", parallel::detectCores()))
  pwr_curv_df <- read_parquet(file.path(gen_path, "power_curve_all.parquet"))
} else {
  pwr_curv_df <- read_parquet(file.path(gen_path, "power_curve.parquet"))
}

sum(pwr_curv_df$norm_potential <= 0) / nrow(pwr_curv_df) * 100
sum(pwr_curv_df$norm_potential >= 1) / nrow(pwr_curv_df) * 100

# pc model : norm_power ~ pc(ws_h, group = site_name)
# penalty :  replicate = pc

# Basic models ####
## Smooth RW beta model #####

pwr_curv_df <- pwr_curv_df %>%
  mutate(
    norm_potential = potential / capacity,
    norm_power_est0 = power_est0 / capacity,
    error0 = norm_potential - norm_power_est0
  ) %>%
  group_by(lon, lat, halfHourEndTime) %>%
  summarise(
    site_name = first(site_name),
    ws_h = mean(ws_h),
    across(c(norm_potential, norm_power_est0, error0), sum),
    .groups = "drop"
  ) %>%
  mutate(
    site_id = as.integer(factor(site_name)),
    pos_val = ifelse(
      norm_potential > 0 & norm_potential < 1,
      pmin(pmax(norm_potential, 1e-6), 1 - 1e-6),
      NA_real_
    )
  )

df <- pwr_curv_df %>%
  filter(!is.na(pos_val)) # keep only rows used in beta model

n_groups <- 20 # start here; increase if you want a finer curve
df$ws_group <- inla.group(df$ws_h, n = n_groups, method = "quantile")
ws_breaks <- df$ws_group %>% unique() %>% sort()

components <- ~ Intercept(1) +
  curve(
    ws_group,
    model = "rw2",
    constr = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
  )

like_beta <- bru_obs(
  formula = pos_val ~ Intercept + curve,
  family = "beta", # or "betar" - check your INLA version
  data = df,
  control.family = list(link = "logit")
)
fit_rw2 <- bru(
  components,
  like_beta,
  options = list(
    control.inla = list(int.strategy = "auto"),
    verbose = TRUE,
    control.compute = list(config = TRUE) # if you later want posterior samples
  )
)
summary(fit_rw2)

ws_grid <- seq(min(df$ws_h), max(df$ws_h), length.out = 400)
pred_df <- data.frame(ws_h = ws_grid)
pred_df$ws_group <- sapply(pred_df$ws_h, function(x) {
  ws_breaks[which.min(abs(ws_breaks - x))]
})
# head(pred_df)
# class(df$ws_group)
# predict linear predictor (logit mean)
pred_lp <- predict(
  fit_rw2,
  pred_df,
  ~ Intercept + curve
)
n_samp <- 1000 # number of posterior samples
samples <- predict(fit_rw2, pred_df, ~ Intercept + curve, n.samples = n_samp)
# predict() returns a list with $mean and $sd for the linear predictor
?predict.bru
eta_mat <- samples$latent # or samples$sample depending on bru version
p_mat <- plogis(eta_mat) # transform each sample to [0,1]
pred_df$mean_p <- rowMeans(p_mat)
pred_df$lower_p <- apply(p_mat, 1, quantile, probs = 0.025)
pred_df$upper_p <- apply(p_mat, 1, quantile, probs = 0.975)

ggplot() +
  geom_line(data = pred_df, aes(x = ws_h, y = mean_p), lwd = 1) +
  geom_ribbon(
    data = pred_df,
    aes(x = ws_h, ymin = lower_p, ymax = upper_p),
    alpha = 0.2
  ) +
  labs(
    x = "Wind speed (m/s)",
    y = "Normalized power",
    title = "RW2 Beta model fit with posterior samples"
  )


pred_df$eta_mean <- pred_lp$mean
pred_df$eta_sd <- pred_lp$sd
pred_df$mean_p <- plogis(pred_df$eta_mean)
# approximate 95% band on link scale -> transform endpoints (delta approx)
pred_df$lower_p <- plogis(pred_df$eta_mean - 1.96 * pred_df$eta_sd)
pred_df$upper_p <- plogis(pred_df$eta_mean + 1.96 * pred_df$eta_sd)

# quick ggplot
ggplot() +
  # geom_point(data = df, aes(x = ws_h, y = pos_val), alpha = 0.3, size = 0.8) +
  geom_line(data = pred_df, aes(x = ws_h, y = mean_p), lwd = 1) +
  geom_ribbon(
    data = pred_df,
    aes(x = ws_h, ymin = lower_p, ymax = upper_p),
    alpha = 0.2
  ) +
  labs(
    x = "Wind speed (m/s)",
    y = "Normalized power",
    title = "RW2 Beta model fit"
  )

## Smooth 1D SPDE beta model ####

# Models with zero inflated properties ####
## Smooth RW beta model #####

## Smooth 1D SPDE beta model ####

# Models with penalty to generic curves ####
## Smooth RW beta model #####

## Smooth 1D SPDE beta model ####
