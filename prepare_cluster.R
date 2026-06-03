# Intall packages in cluster if needed

temp_lib <- "/exports/eddie3_homes_local/s2441782/lib"
if (!dir.exists(temp_lib)) {
  # create directory
  dir.create(temp_lib, recursive = TRUE)
}
.libPaths(temp_lib)
if (!require(INLA)) {
  local({
    r <- c(
      INLA = "https://inla.r-inla-download.org/R/testing",
      CRAN = "https://cloud.r-project.org/",
      inlabru_universe = "https://inlabru-org.r-universe.dev"
    )
    options(repos = r)
  })
  install.packages(
    c("INLA"), #, "inlabru"),
    temp_lib,
    dependencies = TRUE
  )
  require(INLA)
  inla.binary.install()
}
install.packages(
  c("ggsci", "terra", "fst", "arrow"),
  temp_lib,
  dependencies = TRUE
)
install.packages('R.utils', temp_lib, dependencies = TRUE)

dir.create("/export/eddie/scratch/s2441782/Rtmp", showWarnings = FALSE)
Sys.setenv(TMPDIR = "/export/eddie/scratch/s2441782/Rtmp")
install.packages(
  c("rnaturalearth", "rnaturalearthdata", "ggthemes", "brms", "geosphere"),
  temp_lib,
  dependencies = TRUE
)

install.packages(
  c("ggthemes"),
  temp_lib,
  dependencies = TRUE
)

install.packages("ggridges", temp_lib, dependencies = TRUE)
install.packages("ggspatial", temp_lib, dependencies = TRUE)
install.packages("elevatr", temp_lib, dependencies = TRUE)
install.packages("ModelMetrics", temp_lib, dependencies = TRUE)
install.packages("geosphere", temp_lib, dependencies = TRUE)
install.packages("terra", temp_lib, dependencies = TRUE)
