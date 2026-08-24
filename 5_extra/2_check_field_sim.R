cat(
  "-------------------------------------------------------------------------------------------------\n"
)
cat("Running spatial aggregation and coverage check for a simulated field\n")
cat(
  "-------------------------------------------------------------------------------------------------\n"
)

start_time <- Sys.time()

prefix <- "elexon"
set.seed(2026)
n <- nrow(true_df)
oos_perc <- 0.2

colsc <- function(...) {
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, color_name)),
    limits = range(..., na.rm = TRUE)
  )
}
colscB <- function(...) {
  scale_color_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, color_name)),
    limits = range(..., na.rm = TRUE)
  )
}

bndsim <- bnd[[2]] %>% st_as_sf()
## For fmesher 0.3.0:
##   mesh_fine <- fm_mesh_2d_inla(boundary = bnd, max.edge = 0.2)
## For fmesher >= 0.4.0:
range_estimate <- bru0$summary.hyperpar["Range for st_field", "mean"]
sigma_estimate <- bru0$summary.hyperpar["Stdev for st_field", "mean"]
edgeA <- range_estimate / 5
edgeB <- edgeA * 1.25
mesh_fine <- fm_mesh_2d(
  loc = fm_hexagon_lattice(bndsim, edge_len = edgeA),
  boundary = bndsim,
  max.edge = edgeB
)
# ggplot() +
#   geom_fm(data = mesh_fine)

# Note: the priors here will not be used in estimation
matern_fine <-
  inla.spde2.pcmatern(
    mesh_fine,
    prior.sigma = c(sigma_estimate, 0.5),
    prior.range = c(range_estimate, 0.5)
  )
true_range <- round(range_estimate, 0)
true_sigma <- round(sigma_estimate, 0)
extra_noise <- 1 /
  sqrt(bru0$summary.hyperpar["Precision for the Gaussian observations", "mean"])
true_intercept <- bru0$summary.fixed["Intercept", "mean"]
true_Q <- inla.spde.precision(
  matern_fine,
  theta = log(c(range_estimate, sigma_estimate))
)
## plot sd field ####
true_sd <- diag(inla.qinv(true_Q))^0.5
# ggplot() +
#   gg(mesh_fine, color = true_sd) +
#   coord_equal()

## generating samples from model ####

myseed <- trunc(abs(rnorm(1)) * 10000)
true_field <- inla.qsample(1, true_Q, seed = myseed)[, 1] + true_intercept

bb <- sf::st_bbox(bndsim)

truth <- expand.grid(
  easting = seq(bb["xmin"], bb["xmax"], length.out = 100),
  northing = seq(bb["ymin"], bb["ymax"], length.out = 100)
)
# Convert grid points to sf
truth <- sf::st_as_sf(
  truth,
  coords = c("easting", "northing"),
  crs = sf::st_crs(bndsim)
)

# Keep only points inside the polygon
truth <- truth[sf::st_within(truth, bndsim, sparse = FALSE), ]

truth$field <- fm_evaluate(
  mesh_fine,
  loc = truth,
  field = true_field
)

color_name <- "RdGy"
csc <- colsc(truth$field)
# pl_truth <- ggplot() +
#   gg(truth, aes(fill = field), geom = "tile") +
#   ggtitle("True field")
# pl_truth + csc
# ggsave(
#   sprintf(
#     "fig/%s/%s_sim_field_range%d_sd%s.pdf",
#     batch_name,
#     prefix,
#     true_range,
#     sub("\\.", "_", sprintf("%.2f", true_sigma))
#   ),
#   width = 4,
#   height = 4
# )

coord_list_fname <- "data/coord_list.csv"

cat("Loading existing coordinate list\n")
coord_list <- read.csv(coord_list_fname)

sample_loc_fname <- "data/coord_list_wloc.csv"
cat("Loading existing coordinate list with location names\n")
loc_cat <- read.csv(sample_loc_fname) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = st_crs(27700))
loc_cat <- loc_cat %>%
  st_geometry() %>%
  (\(g) g / 1000)() %>%
  st_set_geometry(loc_cat, .)

mydata <-
  true_df %>%
  left_join(
    loc_cat %>%
      dplyr::select(coord_id, short_name) %>%
      st_drop_geometry(),
    by = c("coord_id" = "coord_id")
  )

# mutate(oos = !sampled)

n_fit <- sum(!mydata$oos)
# predict a random sample from bru0 model at mydata locations
the_formula_str <- get_bru_formula(bru0)
grab_prec_name <- bru0$.args$control.family[[1]]$hyper[[
  "theta1"
]]$output.name %>%
  gsub("[ -]+", "_", .) %>%
  as.character()

the_formula_wnoise_str <- sprintf(
  "~ data.frame(
    geometry = geometry,
    coord_id = coord_id,
    time_id = time_id,
    lin_pred = (%s),
    prec = %s,
    lin_pred_noise = rnorm(n, mean = (%s), sd = (%s))
  )",
  the_formula_str,
  grab_prec_name,
  the_formula_str,
  extra_noise
)
the_formula <- as.formula(the_formula_wnoise_str)
simulating_obs <- generate(
  bru0,
  true_df,
  the_formula,
  n.samples = 1
)
# simulating_obs[[1]] %>%
#   left_join(
#     mydata %>% dplyr::select(coord_id, norm_potential, anomaly),
#     by = "coord_id"
#   ) %>%
#   ggplot() +
#   geom_point(
#     aes(x = lin_pred, y = norm_potential, col = anomaly),
#     alpha = 0.1
#   ) +
#   scale_color_manual(values = c("TRUE" = "darkred", "FALSE" = blues9[7])) +
#   labs(
#     x = "Simulated observed data",
#     y = "True normalised potential",
#     col = "Anomaly"
#   )
simulating_obs[[1]] %>%
  left_join(
    mydata %>% dplyr::select(coord_id, norm_potential, anomaly),
    by = "coord_id"
  ) %>%
  ggplot() +
  geom_hex(
    aes(x = lin_pred, y = norm_potential),
    bins = 50,
    # color = "grey50",
    # fill = blues9[5],
    alpha = 0.5
  ) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  labs(
    x = "Simulated observed data",
    y = "True normalised potential",
    col = "Anomaly"
  )
ggsave(
  sprintf(
    "fig/%s/%s_sim_obs_vs_true_range%d_sd%s.pdf",
    batch_name,
    prefix,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 4,
  height = 4
)

mydata <- left_join(
  mydata %>% dplyr::select(-observed),
  simulating_obs[[1]] %>%
    st_drop_geometry() %>%
    dplyr::select(coord_id, lin_pred_noise) %>%
    rename(observed = lin_pred_noise),
  by = "coord_id"
)
# plot(mydata$observed, mydata$norm_potential)
# mydata$observed <- fm_evaluate(
#   mesh_fine,
#   loc = mydata,
#   field = true_field
# ) +
#   rnorm(n, sd = extra_noise)
cscB <- colscB(mydata$observed)

ggplot() +
  geom_sf(data = uk_map, fill = NA, color = "black") +
  gg(mydata, aes(col = observed)) +
  cscB +
  geom_sf(data = bndsim, alpha = 0)


## Estimating the field ####
## For fmesher 0.3.0:
##   mesh <- fm_mesh_2d_inla(
##     boundary = bnd,
##     max.edge = 0.5
##  )
## For fmesher >= 0.4.0
if (mesh_edge_par <= 20) {
  max_n <- c(900, 300)
} else {
  max_n <- c(900, 150)
}

hex_0 <- fm_hexagon_lattice(bnd[[1]], edge_len = edge_target)
mesh <- fm_mesh_2d(
  # loc = loc_unique,
  loc = hex_0,
  boundary = bnd,
  # max.edge = c(100, 150), # km
  min.angle = 30,
  # offset = -0.2,
  cutoff = edge_target / 2,
  max.n.strict = max_n
)
ggplot() +
  geom_fm(data = mesh)

matern <- inla.spde2.pcmatern(
  mesh,
  prior.sigma = c(1, 0.01),
  prior.range = c(1, 0.01)
)

spde <- INLA::inla.spde2.pcmatern(
  mesh = wf.mesh,
  prior.range = c(50, 0.5), # P(range < 100 km)=0.5
  prior.sigma = c(0.2, 0.5) # P(sd > 0.2)=0.5
)

components0 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  # tech_typ(tech_typ, model = "iid") + # random intercept by tech_typ
  norm_power_est0 +
  # power_correction(
  #   pow_group,
  #   model = "rw2",
  #   # replicate = tech_typ,
  #   constr = TRUE
  # ) + # smooth correction power
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
    model = spde
    # group = time_id,
    # control.group = list(model = "ar1")
  )

cat(
  "-------------------------------------------------------------------------------------------------\n"
)
cat("Fitting spatiotemporal model on simulated data\n")
cat(
  "-------------------------------------------------------------------------------------------------\n"
)

# cmp <- observed ~ field(geometry, model = matern) + Intercept(1)
fit <- bru(
  components = components0,
  formula = norm_potential ~ Intercept +
    norm_power_est0 +
    # power_correction +
    d_coast +
    elev +
    wind +
    st_field,
  family = "gaussian",
  data = true_df %>%
    filter(oos == FALSE),
  options = base_bru_options
)
summary(fit)

# residuals
fitted <- fit$summary.fitted.values$mean[1:n_fit]
residuals <- fitted - mydata$observed[mydata$oos == FALSE]
# residuals diagnostics

ggplot() +
  geom_point(aes(x = fitted, y = residuals), color = blues9[7]) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs Fitted", x = "Fitted", y = "Residuals") +
  theme_minimal()
ggsave(
  sprintf(
    "fig/%s/%s_residuals_vs_fitted_range%d_sd%s.pdf",
    batch_name,
    prefix,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 4,
  height = 4
)
ggplot() +
  geom_histogram(aes(residuals), bins = 30, fill = blues9[7]) +
  labs(title = "Histogram of Residuals", x = "Residuals", y = "Frequency") +
  theme_minimal()
ggsave(
  sprintf(
    "fig/%s/%s_histogram_residuals_range%d_sd%s.pdf",
    batch_name,
    prefix,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 4,
  height = 4
)
ggplot() +
  geom_qq(aes(sample = residuals), color = blues9[7]) +
  geom_qq_line(aes(sample = residuals), color = "red") +
  labs(
    title = "Q-Q Plot of Residuals",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_minimal()
ggsave(
  sprintf(
    "fig/%s/%s_qqplot_residuals_range%d_sd%s.pdf",
    batch_name,
    prefix,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 4,
  height = 4
)

# normality test on residuals

norm_test <- shapiro.test(residuals)
cat(
  "Shapiro-Wilk normality test p-value: ",
  norm_test$p.value,
  "\n"
)
cat(ifelse(
  norm_test$p.value < 0.05,
  "Residuals are not normally distributed.\n",
  "Residuals could be normally distributed.\n"
))
# summaries from model fit

pred_on_fulldf <- predict(
  fit,
  mydata,
  as.formula(sprintf("~ %s", the_formula_str)),
)
mydata$fit <- pred_on_fulldf$mean
mydata$lwr <- pred_on_fulldf$`q0.025`
mydata$upr <- pred_on_fulldf$`q0.975`

mydata <- mydata %>%
  mutate(
    above_lwr = observed >= lwr,
    below_upr = observed <= upr,
    in_ci = above_lwr & below_upr
  )

# pix <- fm_pixels(mesh, dims = c(200, 200))
# pred <- predict(
#   fit,
#   pix,
#   ~ field + Intercept
# )
# samp <- generate(
#   fit,
#   pix,
#   ~ field + Intercept,
#   n.samples = 1
# )
# pred$sample <- samp[, 1]

# pl_posterior_mean <- ggplot() +
#   gg(pred, geom = "tile") +
#   gg(bndsim, alpha = 0) +
#   geom_sf(data = mydata, aes(col = in_ci), inherit.aes = FALSE, size = 0.5) +
#   ggtitle("Post. mean")
# pl_posterior_sample <- ggplot() +
#   gg(pred, aes(fill = sample), geom = "tile") +
#   gg(bndsim, alpha = 0) +
#   geom_sf(data = mydata, aes(col = in_ci), inherit.aes = FALSE, size = 0.5) +
#   ggtitle("Post. sample")

# # Common colour scale for the truth and estimate:
# csc <- colsc(truth$field, pred$mean, pred$sample)
# # est_plot <- multiplot(
# #   pl_truth + csc + labels(fill = ""),
# #   pl_posterior_mean + csc,
# #   pl_posterior_sample + csc,
# #   cols = 3
# # )

# (pl_truth +
#   geom_sf(data = mydata, aes(col = in_ci), inherit.aes = FALSE, size = 0.5) +
#   csc +
#   labs(fill = "") |
#   pl_posterior_mean + csc + labs(fill = "") |
#   pl_posterior_sample + csc + labs(fill = "")) +
#   plot_layout(ncol = 3, guides = "collect") &
#   theme(legend.position = "right")

# ggsave(
#   sprintf(
#     "fig/%s/%s_est_field_range%d_sd%s.pdf",
#     batch_name,
#     prefix,
#     true_range,
#     sub("\\.", "_", sprintf("%.2f", true_sigma))
#   ),
#   width = 12,
#   height = 4
# )

int.plot <- plot(fit, "Intercept") +
  geom_vline(
    xintercept = true_intercept,
    col = "red",
    linetype = "dashed"
  ) +
  labs(title = "Intercept posterior", x = "Intercept", y = "Density") +
  theme_bw()
spde.range <- spde.posterior(fit, "st_field", what = "range")
spde.var <- spde.posterior(fit, "st_field", what = "variance")
range.plot <- plot(spde.range) +
  geom_vline(
    xintercept = range_estimate,
    col = "red",
    linetype = "dashed"
  ) +
  labs(title = "Range posterior", x = "Range", y = "Density") +
  theme_bw()
var.plot <- plot(spde.var) +
  geom_vline(
    xintercept = (sigma_estimate^2 + extra_noise^2),
    col = "red",
    linetype = "dashed"
  ) +
  labs(title = "Variance posterior", x = "Variance", y = "Density") +
  theme_bw()

(range.plot / var.plot / int.plot)
ggsave(
  sprintf(
    "fig/%s/%s_est_hyper_range%d_sd%s.pdf",
    batch_name,
    prefix,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 4,
  height = 6
)

csc <- colsc(
  pred[["median"]],
  pred[["q0.025"]],
  pred[["q0.975"]]
) ## Common colour scale from SpatialPixelsDataFrame

# calculate 95% CI coverage ####

coverage <- mean(mydata$in_ci, na.rm = TRUE)
coverage_is <- mean(mydata$in_ci[!mydata$oos], na.rm = TRUE)
coverage_oos <- mean(mydata$in_ci[mydata$oos], na.rm = TRUE)

# gmedian <- ggplot() +
#   gg(pred["median"], geom = "tile") +
#   csc +
#   # gg(mesh) +
#   geom_sf(data = mydata, aes(col = in_ci), inherit.aes = FALSE, size = 0.5) +
#   labs(title = "Posterior median", fill = "") +
#   annotate(
#     "text",
#     x = 5,
#     y = 10.3,
#     label = sprintf("Coverage: %.2f", coverage),
#     size = 4
#   )
# glower95 <- ggplot() +
#   gg(pred["q0.025"], geom = "tile") +
#   csc +
#   geom_sf(data = mydata, aes(col = in_ci), inherit.aes = FALSE, size = 0.5) +
#   labs(title = "Posterior 2.5%", fill = "") +
#   annotate(
#     "text",
#     x = 5,
#     y = 10.3,
#     label = sprintf("Coverage: %.2f", coverage),
#     size = 4
#   )
# gupper95 <- ggplot() +
#   gg(pred["q0.975"], geom = "tile") +
#   csc +
#   geom_sf(data = mydata, aes(col = in_ci), inherit.aes = FALSE, size = 0.5) +
#   labs(title = "Posterior 97.5%", fill = "") +
#   annotate(
#     "text",
#     x = 5,
#     y = 10.3,
#     label = sprintf("Coverage: %.2f", coverage),
#     size = 4
#   )

# (gmedian |
#   glower95 |
#   gupper95) +
#   plot_layout(ncol = 3, guides = "collect") &
#   theme(legend.position = "right")

# ggsave(
#   sprintf(
#     "fig/%s/%s_est_field_range%d_sd%s_ci.pdf",
#     batch_name,
#     prefix,
#     true_range,
#     sub("\\.", "_", sprintf("%.2f", true_sigma))
#   ),
#   width = 12,
#   height = 4
# )

# samples for aggregation ####

# grab_prec_name <- fit$.args$control.family[[1]]$hyper[[
#   "theta1"
# ]]$output.name %>%
#   gsub("[ -]+", "_", .) %>%
#   as.character()

the_formula_wnoise_str <- as.formula(sprintf(
  "~ data.frame(
    geometry = geometry,
    coord_id = coord_id,
    time_id = time_id,
    lin_pred = (%s),
    prec = %s,
    lin_pred_noise = rnorm(n, mean = (%s), sd = 1/sqrt((%s))),
    lwr = qnorm(0.025, mean = (%s), sd = 1/sqrt((%s))),
    upr = qnorm(0.975, mean = (%s), sd = 1/sqrt((%s)))
  )",
  the_formula_str,
  grab_prec_name,
  the_formula_str,
  grab_prec_name,
  the_formula_str,
  grab_prec_name,
  the_formula_str,
  grab_prec_name
))
samp_loc <- generate(
  fit,
  mydata,
  the_formula_wnoise_str,
  n.samples = 1000
)
# samp_loc[[1]]
## double check coverages with samples at locations

## get quantiles
samp_df <- lapply(
  seq_along(samp_loc),
  function(s) {
    data.frame(
      geometry = samp_loc[[s]]$geometry,
      fit0 = samp_loc[[s]]$lin_pred,
      fit = samp_loc[[s]]$lin_pred_noise,
      lwr = samp_loc[[s]]$lwr,
      upr = samp_loc[[s]]$upr
    ) %>%
      mutate(sim = s)
  }
) %>%
  bind_rows()

# variance of samples vs variance of observations check
var_samples <- samp_df %>%
  # group_by(geometry) %>%
  summarise(var_samples = var(fit, na.rm = TRUE)) %>%
  # ungroup() %>%
  # summarise(mean_var_samples = mean(var_samples)) %>%
  pull(var_samples)
var_observed <- mydata %>%
  summarise(var_observed = var(observed, na.rm = TRUE)) %>%
  pull(var_observed)
cat("Mean variance of samples:", var_samples, "\n")
cat("Variance of observed data:", var_observed, "\n")

mydata_2 <- samp_df %>%
  group_by(geometry) %>%
  summarise(
    q0.025 = quantile(fit, probs = 0.025),
    median = quantile(fit, probs = 0.5),
    q0.975 = quantile(fit, probs = 0.975)
  ) %>%
  bind_cols(
    mydata %>%
      st_drop_geometry()
  ) %>%
  mutate(
    above_lwr_noise = observed >= q0.025,
    below_upr_noise = observed <= q0.975,
    in_ci_noise = above_lwr_noise & below_upr_noise
  )

# plot field again with updated in_ci variable

coverage_noise <- mean(mydata_2$in_ci_noise)
coverage_noise_is <- mean(mydata_2$in_ci_noise[!mydata_2$oos])
coverage_noise_oos <- mean(mydata_2$in_ci_noise[mydata_2$oos])

# pl_truth <- ggplot() +
#   gg(truth, aes(fill = field), geom = "tile") +
#   ggtitle("True field")
# pl_posterior_mean <- ggplot() +
#   gg(pred, geom = "tile") +
#   gg(bnd[[2]] %>% st_as_sf(), alpha = 0) +
#   ggtitle("Post. mean")
# pl_posterior_sample <- ggplot() +
#   gg(pred, aes(fill = sample), geom = "tile") +
#   gg(bnd[[2]] %>% st_as_sf(), alpha = 0) +
#   ggtitle("Post. sample")

# # Common colour scale for the truth and estimate:
# csc <- colsc(truth$field, pred$mean, pred$sample)
# updated_points <- geom_sf(
#   data = mydata_2,
#   aes(geometry = geometry, col = in_ci_noise),
#   inherit.aes = FALSE,
#   size = 0.5
# )
# (pl_truth +
#   updated_points +
#   csc +
#   labs(fill = "") |
#   pl_posterior_mean +
#     updated_points +
#     csc +
#     labs(fill = "") |
#   pl_posterior_sample +
#     updated_points +
#     csc +
#     labs(fill = "")) +
#   plot_layout(ncol = 3, guides = "collect") &
#   theme(legend.position = "right")

# ggsave(
#   sprintf(
#     "fig/%s/%s_est_field_range%d_sd%s.pdf",
#     batch_name,
#     prefix,
#     true_range,
#     sub("\\.", "_", sprintf("%.2f", true_sigma))
#   ),
#   width = 12,
#   height = 4
# )

# mydata <-
#   bind_cols(
#     mydata,
#     lapply(
#       samp_loc,
#       function(x) {
#         quantile(x$fit, probs = c(0.025, 0.5, 0.975))
#       }
#     ) %>%
#       t() %>%
#       as.data.frame() %>%
#       setNames(c("q0.025", "median", "q0.975"))
#   )
## scatter plots ####
mydata_2 %>%
  ggplot() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  geom_errorbar(
    aes(x = fit, ymin = lwr, ymax = upr),
    col = "blue",
    alpha = 0.5
  ) +
  geom_point(aes(fit, observed, color = in_ci)) +
  labs(x = "Posterior mean", y = "Observed data")
ggsave(
  sprintf(
    "fig/%s/%s_scatter+error_range%d_sd%s_obs.pdf",
    batch_name,
    prefix,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 4,
  height = 4
)

mydata_2 %>%
  ggplot() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  geom_errorbar(
    aes(x = fit, ymin = q0.025, ymax = q0.975),
    col = "blue",
    alpha = 0.5
  ) +
  geom_point(aes(fit, observed, color = in_ci_noise)) +
  labs(x = "Posterior mean", y = "Observed data")
ggsave(
  sprintf(
    "fig/%s/%s_scatter+obserror_range%d_sd%s_obs_noise.pdf",
    batch_name,
    prefix,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 4,
  height = 4
)

## compare with observed data

aggr_samples <- lapply(
  seq_along(samp_loc),
  function(s) {
    data.frame(
      geometry = samp_loc[[s]]$geometry,
      fit = samp_loc[[s]]$lin_pred_noise
    ) %>%
      summarise(
        fit = mean(fit)
      ) %>%
      mutate(sim = s)
  }
) %>%
  bind_rows() %>%
  summarise(
    date = d0,
    q0.025 = quantile(fit, probs = 0.025),
    median = quantile(fit, probs = 0.5),
    q0.975 = quantile(fit, probs = 0.975),
    fit = mean(fit)
  ) %>%
  bind_cols(
    data.frame(
      observed = mean(mydata_2$observed),
      # fit = mean(mydata_2$fit),
      cov_nonoise = coverage,
      cov_nonoise_is = coverage_is,
      cov_nonoise_oos = coverage_oos,
      cov_noise = coverage_noise,
      cov_noise_is = coverage_noise_is,
      cov_noise_oos = coverage_noise_oos
    )
  ) %>%
  mutate(
    variance_samples = var_samples,
    variance_observed = var_observed,
    aggr_in_ci = ifelse(observed >= q0.025 & observed <= q0.975, 1, 0),
    shapiro_p = norm_test$p.value,
    normality = ifelse(
      norm_test$p.value < 0.05,
      "not normal",
      "maybe normal"
    )
  )
aggr_samples
write.csv(
  aggr_samples,
  sprintf(
    "summaries/%s/%s_aggr_samples_day%s_range%d_sd%s.csv",
    batch_name,
    prefix,
    d0_tag,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  row.names = FALSE
)


end_time <- Sys.time()
cat(
  "Time taken to check aggregation coverage of simulations: ",
  round(difftime(end_time, start_time, units = "mins"), 2),
  " minutes\n"
)
