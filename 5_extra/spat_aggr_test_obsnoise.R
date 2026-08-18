# new tutorial in inlabru ####
require(dplyr)
library(INLA)
library(inlabru)
library(fmesher)
library(mgcv)
library(ggplot2)
require(patchwork)
require(sf)

theme_set(theme_bw())
prefix <- "wnoise"
set.seed(2026)
n <- 150

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

bnd <- spoly(
  data.frame(
    easting = c(0, 10, 10, 0),
    northing = c(0, 0, 10, 10)
  ),
  format = "sf"
)
## For fmesher 0.3.0:
##   mesh_fine <- fm_mesh_2d_inla(boundary = bnd, max.edge = 0.2)
## For fmesher >= 0.4.0:
edgeA <- 0.25
edgeB <- edgeA * 1.25
mesh_fine <- fm_mesh_2d(
  loc = fm_hexagon_lattice(bnd, edge_len = edgeA),
  boundary = bnd,
  max.edge = edgeB
)
ggplot() +
  geom_fm(data = mesh_fine)


# Note: the priors here will not be used in estimation
matern_fine <-
  inla.spde2.pcmatern(
    mesh_fine,
    prior.sigma = c(1, 0.01),
    prior.range = c(1, 0.01)
  )
true_range <- 4
true_sigma <- 1
true_Q <- inla.spde.precision(
  matern_fine,
  theta = log(c(true_range, true_sigma))
)
## plot sd field ####
true_sd <- diag(inla.qinv(true_Q))^0.5
ggplot() +
  gg(mesh_fine, color = true_sd) +
  coord_equal()

## generating samples from model ####

myseed <- trunc(abs(rnorm(1)) * 10000)
true_field <- inla.qsample(1, true_Q, seed = myseed)[, 1]

truth <- expand.grid(
  easting = seq(0, 10, length = 100),
  northing = seq(0, 10, length = 100)
)
truth <- sf::st_as_sf(truth, coords = c("easting", "northing"))
truth$field <- fm_evaluate(
  mesh_fine,
  loc = truth,
  field = true_field
)

color_name <- "RdGy"
csc <- colsc(truth$field)
pl_truth <- ggplot() +
  gg(truth, aes(fill = field), geom = "tile") +
  ggtitle("True field")
pl_truth + csc
ggsave(
  sprintf(
    "fig/testing/%s_sim_field_range%d_sd%s.pdf",
    prefix,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 4,
  height = 4
)

## Or with another colour scale:

# multiplot(pl_truth, pl_truth + csc, cols = 2)

mydata <- sf::st_as_sf(
  data.frame(easting = runif(n, 0, 10), northing = runif(n, 0, 10)),
  coords = c("easting", "northing")
)
mydata$observed <-
  fm_evaluate(
    mesh_fine,
    loc = mydata,
    field = true_field
  ) +
  rnorm(n, sd = 0.4)
cscB <- colscB(truth$field)
ggplot() +
  gg(mydata, aes(col = observed)) +
  cscB


## Estimating the field ####
## For fmesher 0.3.0:
##   mesh <- fm_mesh_2d_inla(
##     boundary = bnd,
##     max.edge = 0.5
##  )
## For fmesher >= 0.4.0
mesh <- fm_mesh_2d(
  loc = fm_hexagon_lattice(bnd, edge_len = edgeA * 2),
  boundary = bnd,
  max.edge = edgeB * 2
)
ggplot() +
  geom_fm(data = mesh)


matern <-
  inla.spde2.pcmatern(mesh, prior.sigma = c(10, 0.01), prior.range = c(1, 0.01))

cmp <- observed ~ field(geometry, model = matern) + Intercept(1)
fit <- bru(cmp, mydata, family = "gaussian")
summary(fit)

# summaries from model fit
mydata$fit <- fit$summary.fitted.values$mean[1:n]
mydata$lwr <- fit$summary.fitted.values$`0.025quant`[1:n]
mydata$upr <- fit$summary.fitted.values$`0.975quant`[1:n]

mydata <- mydata %>%
  mutate(
    above_lwr = observed >= lwr,
    below_upr = observed <= upr,
    in_ci = above_lwr & below_upr
  )

pix <- fm_pixels(mesh, dims = c(200, 200))
pred <- predict(
  fit,
  pix,
  ~ field + Intercept
)
samp <- generate(
  fit,
  pix,
  ~ field + Intercept,
  n.samples = 1
)
pred$sample <- samp[, 1]

pl_posterior_mean <- ggplot() +
  gg(pred, geom = "tile") +
  gg(bnd, alpha = 0) +
  geom_sf(data = mydata, aes(col = in_ci), inherit.aes = FALSE, size = 0.5) +
  ggtitle("Post. mean")
pl_posterior_sample <- ggplot() +
  gg(pred, aes(fill = sample), geom = "tile") +
  gg(bnd, alpha = 0) +
  geom_sf(data = mydata, aes(col = in_ci), inherit.aes = FALSE, size = 0.5) +
  ggtitle("Post. sample")

# Common colour scale for the truth and estimate:
csc <- colsc(truth$field, pred$mean, pred$sample)
# est_plot <- multiplot(
#   pl_truth + csc + labels(fill = ""),
#   pl_posterior_mean + csc,
#   pl_posterior_sample + csc,
#   cols = 3
# )

(pl_truth +
  geom_sf(data = mydata, aes(col = in_ci), inherit.aes = FALSE, size = 0.5) +
  csc +
  labs(fill = "") |
  pl_posterior_mean + csc + labs(fill = "") |
  pl_posterior_sample + csc + labs(fill = "")) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "right")

ggsave(
  sprintf(
    "fig/testing/%s_est_field_range%d_sd%s.pdf",
    prefix,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 12,
  height = 4
)

int.plot <- plot(fit, "Intercept")
spde.range <- spde.posterior(fit, "field", what = "range")
spde.var <- spde.posterior(fit, "field", what = "variance")
range.plot <- plot(spde.range) +
  geom_vline(
    xintercept = true_range,
    col = "red",
    linetype = "dashed"
  ) +
  labs(title = "Range posterior", x = "Range", y = "Density") +
  theme_bw()
var.plot <- plot(spde.var) +
  geom_vline(
    xintercept = (true_sigma^2),
    col = "red",
    linetype = "dashed"
  ) +
  labs(title = "Variance posterior", x = "Log variance", y = "Density") +
  theme_bw()

(range.plot / var.plot / int.plot)
ggsave(
  sprintf(
    "fig/testing/%s_est_hyper_range%d_sd%s.pdf",
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

coverage <- sum(
  mydata$in_ci
) /
  n

gmedian <- ggplot() +
  gg(pred["median"], geom = "tile") +
  csc +
  # gg(mesh) +
  geom_sf(data = mydata, aes(col = in_ci), inherit.aes = FALSE, size = 0.5) +
  labs(title = "Posterior median", fill = "") +
  annotate(
    "text",
    x = 5,
    y = 10.3,
    label = sprintf("Coverage: %.2f", coverage),
    size = 4
  )
glower95 <- ggplot() +
  gg(pred["q0.025"], geom = "tile") +
  csc +
  geom_sf(data = mydata, aes(col = in_ci), inherit.aes = FALSE, size = 0.5) +
  labs(title = "Posterior 2.5%", fill = "") +
  annotate(
    "text",
    x = 5,
    y = 10.3,
    label = sprintf("Coverage: %.2f", coverage),
    size = 4
  )
gupper95 <- ggplot() +
  gg(pred["q0.975"], geom = "tile") +
  csc +
  geom_sf(data = mydata, aes(col = in_ci), inherit.aes = FALSE, size = 0.5) +
  labs(title = "Posterior 97.5%", fill = "") +
  annotate(
    "text",
    x = 5,
    y = 10.3,
    label = sprintf("Coverage: %.2f", coverage),
    size = 4
  )

(gmedian |
  glower95 |
  gupper95) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "right")

ggsave(
  sprintf(
    "fig/testing/%s_est_field_range%d_sd%s_ci.pdf",
    prefix,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 12,
  height = 4
)

# samples for aggregation ####

grab_prec_name <- fit$.args$control.family[[1]]$hyper[[
  "theta1"
]]$output.name %>%
  gsub("[ -]+", "_", .) %>%
  as.character()

# samp_loc <- generate(
#   fit,
#   mydata,
#   ~ field + Intercept,
#   n.samples = 1000
# )
samp_loc <- generate(
  fit,
  mydata,
  as.formula(sprintf(
    "~ data.frame(
    geometry = geometry,
    lin_pred = field + Intercept,
    prec = %s,
    lin_pred_noise = rnorm(n, mean = field + Intercept, sd = 1/sqrt((%s))),
    lwr = qnorm(0.025, mean = field + Intercept, sd = 1/sqrt((%s))),
    upr = qnorm(0.975, mean = field + Intercept, sd = 1/sqrt((%s)))
  )",
    grab_prec_name,
    grab_prec_name,
    grab_prec_name,
    grab_prec_name
  )),
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
      fit = samp_loc[[s]]$lin_pred_noise,
      lwr = samp_loc[[s]]$lwr,
      upr = samp_loc[[s]]$upr
    ) %>%
      mutate(sim = s)
  }
) %>%
  bind_rows()

mydata <- samp_df %>%
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

coverage_noise <- sum(
  mydata$in_ci_noise
) /
  n
pl_truth <- ggplot() +
  gg(truth, aes(fill = field), geom = "tile") +
  ggtitle("True field")
pl_posterior_mean <- ggplot() +
  gg(pred, geom = "tile") +
  gg(bnd, alpha = 0) +
  ggtitle("Post. mean")
pl_posterior_sample <- ggplot() +
  gg(pred, aes(fill = sample), geom = "tile") +
  gg(bnd, alpha = 0) +
  ggtitle("Post. sample")

# Common colour scale for the truth and estimate:
csc <- colsc(truth$field, pred$mean, pred$sample)
updated_points <- geom_sf(
  data = mydata,
  aes(geometry = geometry, col = in_ci_noise),
  inherit.aes = FALSE,
  size = 0.5
)
(pl_truth +
  updated_points +
  csc +
  labs(fill = "") |
  pl_posterior_mean +
    updated_points +
    csc +
    labs(fill = "") |
  pl_posterior_sample +
    updated_points +
    csc +
    labs(fill = "")) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "right")

ggsave(
  sprintf(
    "fig/testing/%s_est_field_range%d_sd%s.pdf",
    prefix,
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 12,
  height = 4
)

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

mydata %>%
  ggplot() +
  geom_point(aes(fit, observed), col = "blue") +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  geom_errorbar(
    aes(x = fit, ymin = q0.025, ymax = q0.975),
    col = "blue",
    alpha = 0.5
  ) +
  labs(x = "Posterior mean", y = "Observed data")

## compare with observed data

coverage_loc <- sum(
  mydata$observed >= mydata$q0.025 &
    mydata$observed <= mydata$q0.975
) /
  n

aggr_samples <- lapply(
  seq_along(samp_loc),
  function(s) {
    data.frame(
      geometry = samp_loc[[s]]$geometry,
      fit = samp_loc[[s]]$lin_pred_noise
    ) %>%
      summarise(
        fit = sum(fit)
      ) %>%
      mutate(sim = s)
  }
) %>%
  bind_rows() %>%
  summarise(
    q0.025 = quantile(fit, probs = 0.025),
    median = quantile(fit, probs = 0.5),
    q0.975 = quantile(fit, probs = 0.975)
  ) %>%
  bind_cols(
    data.frame(
      observed = sum(mydata$observed),
      fit = sum(mydata$fit),
      cov_loc = coverage_loc
    )
  )
