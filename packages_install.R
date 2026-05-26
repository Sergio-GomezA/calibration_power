# installing packages #####

install.packages(c(
  "tidyverse",
  "parallel",
  "ggsci",
  "kableExtra",
  "arrow",
  "sf"
))

install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("reshape2")
install.packages("forecast")
install.packages("ggthemes")


install.packages("fields")
install.packages("ggh4x")


install.packages("data.table")

install.packages("terra")
install.packages("fst")
install.packages("stringr")
install.packages("httr")
install.packages("jsonlite")
install.packages("geosphere")
install.packages("lubridate")
install.packages("mgcv")


local({
  r <- c(
    INLA = "https://inla.r-inla-download.org/R/testing",
    CRAN = "https://cloud.r-project.org/",
    inlabru_universe = "https://inlabru-org.r-universe.dev"
  )
  options(repos = r)
})

install.packages(c("INLA", "inlabru"), dependencies = TRUE)
inla.binary.install()

df <- data.frame(y = rnorm(100) + 10)
fit <- INLA::inla(
  y ~ 1,
  data = df
)
summary(fit)


# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install(c("graph", "Rgraphviz"), dep = TRUE)

install.packages("janitor")
install.packages("ModelMetrics")


install.packages("qmap")
install.packages("ggridges")
install.packages('R.utils')

install.packages("GGally")
install.packages("ggspatial")

install.packages("plotly")

remove.packages(c("sf", "lwgeom", "stars", "terra", "units"))
install.packages(
  c(
    "sf",
    "terra",
    "lwgeom",
    "units",
    "stars"
  ),
  type = "source"
)


install.packages("sf")
install.packages("sf", type = "source")
remove.packages("stringi")
install.packages("stringi", type = "source")
install.packages("ggthemes", type = "source")

install.packages("elevatr")
# Sys.getenv("PATH")
# Sys.which("gdal-config")

# Sys.setenv(PATH = paste(
#   "/usr/bin",
#   "/bin",
#   "/usr/sbin",
#   "/sbin",
#   "/usr/local/bin",
#   sep = ":"
# ))

# Sys.which("gdal-config")
