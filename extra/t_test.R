require(INLA)
set.seed(123)
n <- 500
phi <- 0.4
sig_ar <- sqrt(1 / (1 - phi^2))
mu <- 0
eta <- rep(0, n)
for (i in 2:n) {
  eta[i] <- mu + phi * (eta[i - 1] - mu) + rnorm(1)
}

plot(eta, type = "l", main = "AR(1) Process", xlab = "Time", ylab = "Value")

nu <- 2.5
t <- rt(n, df = nu)
plot(
  t,
  type = "l",
  main = "t-distributed Random Variables",
  xlab = "Index",
  ylab = "Value"
)
y <- eta + t / (sqrt(nu / (nu - 2)))
plot(
  y,
  type = "l",
  main = "AR(1) Process with t-distributed Noise",
  xlab = "Time",
  ylab = "Value"
)
data <- list(y = y, z = seq(1:n))
#define the model and fit
formula <- y ~ f(z, model = "ar1")
result <- inla(formula, family = "T", data = data)
summary(result)

result$summary.hyperpar

samples <- inla.hyperpar.sample(
  n = 5000,
  result
)

sd_t <- sapply(samples[, 1], function(x) sqrt(1 / x))
summary(sd_t)

sd_z <- sapply(samples[, 3], function(x) sqrt(1 / x))
summary(sd_z)

summary(samples[, 2])
nu_test <- mapply(
  function(tau, s) {
    nu <- 2 * s * tau / (s * tau - 1)
  },
  samples[, 1],
  samples[, 2]
)
summary(nu_test)


sd_samples <- 1 / sqrt(samples[, 1])
summary(sd_samples)


hyper <- list(
  theta2 = list(
    initial = log(5 - 2),
    fixed = TRUE
  )
)

result_fixed <- inla(
  formula,
  family = "T",
  data = data,
  control.family = list(
    hyper = hyper
  )
)
summary(result_fixed)
samples <- inla.hyperpar.sample(
  n = 5000,
  result_fixed
)

tau_samples <- samples[, 1]

sd_samples <- 1 / sqrt(tau_samples)
summary(sd_samples)
