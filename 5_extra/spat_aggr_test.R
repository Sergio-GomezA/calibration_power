# new tutorial in inlabru ####

library(INLA)
library(inlabru)
library(fmesher)
library(mgcv)
library(ggplot2)

colsc <- function(...) {
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
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
set.seed(2024)
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

pl_truth <- ggplot() +
  gg(truth, aes(fill = field), geom = "tile") +
  ggtitle("True field")
pl_truth
ggsave(sprintf(
  "fig/testing/sim_field_range%d_sd%s.pdf",
  true_range,
  sub("\\.", "_", sprintf("%.2f", true_sigma))
))
