# new tutorial in inlabru ####
require(dplyr)
library(INLA)
library(inlabru)
library(fmesher)
library(mgcv)
library(ggplot2)
require(patchwork)

theme_set(theme_bw())

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
    "fig/testing/sim_field_range%d_sd%s.pdf",
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
    "fig/testing/est_field_range%d_sd%s.pdf",
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 12,
  height = 4
)

int.plot <- plot(fit, "Intercept")
spde.range <- spde.posterior(fit, "field", what = "range")
spde.logvar <- spde.posterior(fit, "field", what = "variance")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)

(range.plot / var.plot / int.plot)


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
    "fig/testing/est_field_range%d_sd%s_ci.pdf",
    true_range,
    sub("\\.", "_", sprintf("%.2f", true_sigma))
  ),
  width = 12,
  height = 4
)

# samples for aggregation ####

samp_loc <- generate(
  fit,
  mydata,
  ~ field + Intercept,
  n.samples = 1000
)

## double check coverages with samples at locations

## get quantiles

mydata <-
  bind_cols(
    mydata,
    apply(
      samp_loc,
      1,
      function(x) {
        quantile(x, probs = c(0.025, 0.5, 0.975))
      }
    ) %>%
      t() %>%
      as.data.frame() %>%
      setNames(c("q0.025", "median", "q0.975"))
  )

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

aggr_samples <- samp_loc %>%
  as.data.frame() %>%
  summarise_all(sum) %>%
  t()

aggr_samples %>%
  quantile(probs = c(0.025, 0.5, 0.975)) %>%
  t() %>%
  as.data.frame() %>%
  setNames(c("q0.025", "median", "q0.975")) %>%
  bind_cols(
    data.frame(
      observed = sum(mydata$observed)
    )
  )
