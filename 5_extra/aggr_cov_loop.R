# loop of field simulation, aggregation and coverage calculation ####

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

if (local_run) {
  cat("Running in local mode\n")
} else {
  cat("Running in cluster mode\n")
}

require(parallel)

if (!local_run) {
  temp_lib <- "/exports/eddie3_homes_local/s2441782/lib"
  .libPaths(temp_lib)
}

# new tutorial in inlabru ####
require(dplyr)
library(INLA)
library(inlabru)
library(fmesher)
library(mgcv)
library(ggplot2)
# require(patchwork)

require(sf)

source("aux_funct.R")
mc <- available_cores() - ifelse(local_run, 2, 0)
inla_core_option <- "%d:1"
cat("Setting INLA to use", mc, "cores\n")
inla.setOption(num.threads = sprintf(inla_core_option, mc))

theme_set(theme_bw())

set.seed(2026)
n <- 150
N <- 100
oos_perc <- 0.2

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
true_sigma <- 2
extra_noise <- 1
true_intercept <- 2
true_Q <- inla.spde.precision(
  matern_fine,
  theta = log(c(true_range, true_sigma))
)

# true field
# myseed <- trunc(abs(rnorm(1)) * 10000)

check_agg_cov <- function(n) {
  true_field <- inla.qsample(1, true_Q)[, 1] + true_intercept

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
  ) %>%
    mutate(
      oos = FALSE
    )
  mydata$oos[sample(nrow(mydata), round(oos_perc * nrow(mydata)))] <- TRUE
  n_fit <- sum(!mydata$oos)
  mydata$observed <-
    fm_evaluate(
      mesh_fine,
      loc = mydata,
      field = true_field
    ) +
    rnorm(n, sd = extra_noise)
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
  fit <- bru(cmp, mydata %>% filter(!oos), family = "gaussian")
  summary(fit)

  # summaries from model fit
  pred_on_fulldf <- predict(
    fit,
    mydata,
    ~ field + Intercept
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

  samp_loc <- generate(
    fit,
    mydata,
    ~ field + Intercept,
    n.samples = 1000
  )

  # calculate coverage
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

  aggr_samples %>%
    quantile(probs = c(0.025, 0.5, 0.975)) %>%
    t() %>%
    as.data.frame() %>%
    setNames(c("q0.025", "median", "q0.975")) %>%
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

write.csv(
  coverage_results,
  file = sprintf("5_extra/coverage_results_nonoise_nloc%d_sim%d.csv", n, N),
  row.names = FALSE
)
