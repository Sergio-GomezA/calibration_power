require(inlabru)

# inla.doc("matern")
# install.packages("viridis")

nrow <- 50
ncol <- 70
n <- nrow * ncol
s.noise <- 1
zi.mat <- matrix(NA, nrow = nrow, ncol = ncol)
i <- 1:nrow
for (j in 1:ncol) {
  zi.mat[i, j] <- 3 * exp(-(i - j)^2 / 30)
}

plot(zi.mat, main = "True signal")
## iid noise
noise.mat <- matrix(rnorm(nrow * ncol, sd = s.noise), nrow, ncol)
## make simulated data with no spatial component
y.mat <- zi.mat + noise.mat
plot(y.mat, main = "Observed data")
## convert matrices to the internal representation in INLA
y <- inla.matrix2vector(y.mat)
node <- 1:n
formula <- y ~ 1 +
  f(
    node,
    model = "matern2d",
    nu = 1,
    nrow = nrow,
    ncol = ncol,
    hyper = list(
      range = list(param = c(1, 1), prior = "loggamma", initial = 1),
      prec = list(param = c(1, 1))
    )
  )
data <- data.frame(y = y, node = node)
## fit the model
result <- inla(
  formula,
  family = "gaussian",
  data = data,
  verbose = TRUE,
  control.predictor = list(compute = TRUE),
  control.family = list(
    hyper = list(theta = list(initial = log(1 / s.noise^2), fixed = FALSE))
  ),
  keep = F
)
summary(result)


INLA:::inla.display.matrix(zi.mat)

INLA:::inla.display.matrix(INLA:::inla.vector2matrix(
  result$summary.linear.predictor$mean,
  nrow = nrow,
  ncol = ncol
))


inla.spde2.pcmatern()


## Spatial interpolation
n <- 100
field.fcn <- function(loc) (10 * cos(2 * pi * 2 * (loc[, 1] + loc[, 2])))
loc <- matrix(runif(n * 2), n, 2)
## One field, 2 observations per location
idx.y <- rep(1:n, 2)
y <- field.fcn(loc[idx.y, ]) + rnorm(length(idx.y))

mesh <- fm_mesh_2d_inla(loc, max.edge = 0.05, cutoff = 0.01)
spde <- inla.spde2.pcmatern(
  mesh,
  prior.range = c(0.01, 0.1),
  prior.sigma = c(100, 0.1)
)
data <- list(y = y, field = mesh$idx$loc[idx.y])
formula <- y ~ -1 + f(field, model = spde)
result <- inla(formula, data = data, family = "normal")

## Plot the mesh structure:
plot(mesh)

if (require(rgl)) {
  col.pal <- colorRampPalette(c("blue", "cyan", "green", "yellow", "red"))
  ## Plot the posterior mean:
  plot_rgl(
    mesh,
    result$summary.random$field[, "mean"],
    color.palette = col.pal
  )
  ## Plot residual field:
  plot_rgl(
    mesh,
    result$summary.random$field[, "mean"] - field.fcn(mesh$loc),
    color.palette = col.pal
  )
}


result.field <- inla.spde.result(result, "field", spde)
par(mfrow = c(2, 1))
plot(
  result.field$marginals.range.nominal[[1]],
  type = "l",
  main = "Posterior density for range"
)
plot(
  inla.tmarginal(sqrt, result.field$marginals.variance.nominal[[1]]),
  type = "l",
  main = "Posterior density for std.dev."
)
par(mfrow = c(1, 1))

## Spatial model
set.seed(1234234)

## Generate spatial locations
nObs <- 200
loc <- matrix(runif(nObs * 2), nrow = nObs, ncol = 2)

## Generate observation of spatial field
nu <- 1.0
rhoT <- 0.2
kappaT <- sqrt(8 * nu) / rhoT
sigT <- 1.0
Sig <- sigT^2 *
  inla.matern.cov(
    nu = nu,
    kappa = kappaT,
    x = as.matrix(dist(loc)),
    d = 2,
    corr = TRUE
  )
L <- t(chol(Sig))
u <- L %*% rnorm(nObs)

## Construct observation with nugget
sigN <- 0.1
y <- u + sigN * rnorm(nObs)

## Create the mesh and spde object
mesh <- fm_mesh_2d_inla(loc, max.edge = 0.05, cutoff = 0.01)
spde <- inla.spde2.pcmatern(
  mesh,
  prior.range = c(0.01, 0.05),
  prior.sigma = c(10, 0.05)
)

## Create projection matrix for observations
A <- fm_basis(mesh = mesh, loc = loc)

## Run model without any covariates
idx <- 1:spde$n.spde
res <- inla(
  y ~ f(idx, model = spde) - 1,
  data = list(y = y, idx = idx, spde = spde),
  control.predictor = list(A = A)
)

## Re-run model with fixed range
spde.fixed <- inla.spde2.pcmatern(
  mesh,
  prior.range = c(0.2, NA),
  prior.sigma = c(10, 0.05)
)

res.fixed <- inla(
  y ~ f(idx, model = spde) - 1,
  data = list(y = y, idx = idx, spde = spde.fixed),
  control.predictor = list(A = A)
)
