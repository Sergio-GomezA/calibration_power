require(inlabru)

inla.doc("matern")
install.packages("viridis")

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
