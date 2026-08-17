# loop of field simulation, aggregation and coverage calculation ####
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
N <- 100

bnd <- spoly(
  data.frame(
    easting = c(0, 10, 10, 0),
    northing = c(0, 0, 10, 10)
  ),
  format = "sf"
)
edgeA <- 0.25
edgeB <- edgeA * 1.25
mesh_fine <- fm_mesh_2d(
  loc = fm_hexagon_lattice(bnd, edge_len = edgeA),
  boundary = bnd,
  max.edge = edgeB
)
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

# true field
# myseed <- trunc(abs(rnorm(1)) * 10000)

check_agg_cov <- function(n) {
  true_field <- inla.qsample(1, true_Q)[, 1]

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

  # generate data
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
  # fit model
  mesh <- fm_mesh_2d(
    loc = fm_hexagon_lattice(bnd, edge_len = edgeA * 2),
    boundary = bnd,
    max.edge = edgeB * 2
  )
  # ggplot() +
  #   geom_fm(data = mesh)

  matern <-
    inla.spde2.pcmatern(
      mesh,
      prior.sigma = c(10, 0.01),
      prior.range = c(1, 0.01)
    )

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

  grab_prec_name <- fit$.args$control.family[[1]]$hyper[[
    "theta1"
  ]]$output.name %>%
    gsub("[ -]+", "_", .) %>%
    as.character()
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

  # calculate coverage
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
      above_lwr = observed >= q0.025,
      below_upr = observed <= q0.975,
      in_ci = above_lwr & below_upr
    )

  coverage_loc <- sum(
    mydata$observed >= mydata$q0.025 &
      mydata$observed <= mydata$q0.975
  ) /
    n
  # calculate aggregate coverage
  aggr_samples <- samp_loc %>%
    as.data.frame() %>%
    summarise_all(sum) %>%
    t()

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
}


cat("Progress: ")
n_failed <- 0

coverage_results <- lapply(1:N, function(i) {
  if (i %in% round(seq(1, N, length.out = 50))) {
    cat("-")
    flush.console()
  }

  tryCatch(
    check_agg_cov(n),
    error = function(e) {
      n_failed <<- n_failed + 1
      NULL
    }
  )
}) %>%
  bind_rows()

cat(sprintf(" done (%d failed)\n", n_failed))

coverage_results %>%
  summarise(
    cov_loc = mean(cov_loc),
    cov_aggr = mean(observed >= q0.025 & observed <= q0.975)
  )
