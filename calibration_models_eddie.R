# packages

temp_lib <- "/exports/eddie/scratch/s2441782/calibration/lib"
if (!require(inlabru)) {
  local({
    r <- c(
      INLA = "https://inla.r-inla-download.org/R/testing",
      CRAN = "https://cloud.r-project.org/",
      inlabru_universe = "https://inlabru-org.r-universe.dev"
    )
    options(repos = r)
  })
  install.packages(
    c("INLA", "inlabru"),
    temp_lib,
    dependencies = TRUE
  )
  inla.binary.install()
}
install.packages(
  c("ggsci", "terra", "fst", "arrow"),
  temp_lib,
  dependencies = TRUE
)
install.packages('R.utils', temp_lib, dependencies = TRUE)


.libPaths(temp_lib)

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

require(brms)
require(INLA, lib.loc = temp_lib)

require(inlabru, lib.loc = temp_lib)

setwd("/exports/eddie/scratch/s2441782/calibration")

source("aux_funct.R")


gen_path <- "data"
model_path <- "model_objects"

GB_df <- read_parquet(file.path(gen_path, "GB_aggr_frag.parquet"))

componentsar1 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  norm_power_est0(norm_power_est0, model = "linear") + # fixed slope
  # tech_typ(tech_typ, model = "iid") + # random intercept by tech_typ
  # tech_power(tech_typ, model = "iid", replicate = norm_power_est0) + # random slope
  month(month, model = "rw2", cyclic = TRUE) +
  # hour(hour, model = "rw2", cyclic = TRUE) +
  wind(ws_group, model = "rw2") +
  u(t, model = "ar1", replicate = tech_typ)

bruar1 <- bru(
  components = componentsar1,
  formula = norm_potential ~ Intercept +
    norm_power_est0 +
    # tech_typ +
    # tech_power +
    month +
    # hour +
    wind +
    u,
  family = "gaussian",
  data = GB_df,
  options = list(
    control.predictor = list(
      compute = FALSE # don't store posterior for all latent effects
    ),
    control.compute = list(
      dic = TRUE, # you can keep DIC or WAIC if needed
      waic = TRUE,
      cpo = FALSE # skip CPO to save memory
    ),
    verbose = TRUE
  )
)
saveRDS(bruar1, file.path(model_path, "bru_ar1.rds"))
ModelMetrics::rmse(
  bruar1$summary.fitted.values$mean[1:nrow(GB_df)],
  GB_df$norm_potential
)


componentsar1 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  norm_power_est0(norm_power_est0, model = "linear") + # fixed slope
  tech_typ(tech_typ, model = "iid") + # random intercept by tech_typ
  tech_power(tech_typ, model = "iid", replicate = norm_power_est0) + # random slope
  month(month, model = "rw2", cyclic = TRUE) +
  hour(hour, model = "rw2", cyclic = TRUE) +
  wind(ws_group, model = "rw2") +
  u(t, model = "ar1", replicate = tech_typ)

bruar1 <- bru(
  components = componentsar1,
  formula = norm_potential ~ Intercept +
    norm_power_est0 +
    tech_typ +
    tech_power +
    month +
    hour +
    wind +
    u,
  family = "gaussian",
  data = GB_df,
  options = list(
    control.predictor = list(
      compute = FALSE # don't store posterior for all latent effects
    ),
    control.compute = list(
      dic = TRUE, # you can keep DIC or WAIC if needed
      waic = TRUE,
      cpo = FALSE # skip CPO to save memory
    ),
    verbose = TRUE
  )
)
saveRDS(bruar1, file.path(model_path, "bru_ar1_full.rds"))
