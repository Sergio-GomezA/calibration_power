# model fit 2.0

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE


require(parallel)

if (local_run) {
  data_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
} else {
  data_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  gen_path <- "/exports/eddie/scratch/s2441782/calibration/data"
  model_path <- "/exports/eddie/scratch/s2441782/calibration/model_objects"
  temp_lib <- "/exports/eddie3_homes_local/s2441782/lib"
  .libPaths(temp_lib)
}


# install.packages(
#   c("qmap"),
#   temp_lib,
#   dependencies = TRUE
# )

require(arrow)
require(dplyr)
require(tidyr)
# require(rnaturalearth)
# require(rnaturalearthdata)
require(sf)
require(ggplot2)
require(ggthemes)
require(ggsci)
# require(FNN)
require(data.table)
# require(parallel)
library(purrr)
# require(brms)
require(INLA)
require(inlabru)
require(qmap)

source("aux_funct.R")

GB_df <- read_parquet(file.path(gen_path, "GB_aggr.parquet"))
GB_df %>%
  pull(halfHourEndTime) %>%
  range()

GB_df %>%
  ggplot(aes(norm_potential, norm_power_est0)) +
  geom_hex() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Elexon CF", y = "Generic PC estimate")
ggsave("fig/gb_est0_hexbin_all.pdf", width = 6, height = 4)

GB_df <- read_parquet(file.path(gen_path, "GB_aggr_frag.parquet"))

base_model <- lm(
  norm_potential ~ norm_power_est0,
  data = GB_df
)

full_model <- lm(
  norm_potential ~ tech_typ *
    norm_power_est0 +
    norm_power_est0 * month +
    hour +
    poly(ws_h_wmean, 3),
  data = GB_df
)
summary(full_model)


model_AIC <- step(
  base_model,
  scope = list(lower = base_model, upper = full_model),
  # steps = 5,
  k = 2
)
model_AIC %>%
  saveRDS(
    file.path(
      model_path,
      "lm0.rds"
    )
  )


components0 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  norm_power_est0(norm_power_est0, model = "linear") + # fixed slope
  tech_typ(tech_typ, model = "iid") + # random intercept by tech_typ
  tech_power(tech_typ, model = "iid", replicate = norm_power_est0) + # random slope
  month(month, model = "rw2", cyclic = TRUE) +
  # month_power(month, model = "iid", replicate = norm_power_est0) + # random slope
  hour(hour, model = "rw2", cyclic = TRUE) +
  wind(ws_group, model = "rw2")

bru0 <- bru(
  components = components0,
  formula = norm_potential ~ Intercept +
    norm_power_est0 +
    tech_typ +
    tech_power +
    month +
    hour +
    wind,
  family = "gaussian",
  data = GB_df
)
summary(bru0)
saveRDS(bru0, file = file.path(model_path, "lm_bru0.rds"))


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
# arimax option ####
library(nlme)

full_model_ar1 <- gls(
  norm_potential ~ tech_typ *
    norm_power_est0 +
    norm_power_est0 * month +
    month +
    hour +
    poly(ws_h_wmean, 3),
  data = GB_df,
  correlation = corAR1(form = ~ t | tech_typ), # AR1 per tech type
  verbose = TRUE
)
saveRDS(full_model_ar1, file = file.path(model_path, "gls_ar1.rds"))

qqmod <- fitQmap(
  obs = GB_df %>% pull(norm_potential),
  mod = GB_df %>% pull(norm_power_est0),
  method = "QUANT"
)

wgen_qm <- with(
  GB_df,
  doQmapQUANT(norm_power_est0, qqmod, type = "linear")
)


# read models ####

# seasonal lm
model_AIC <- readRDS(file.path(model_path, "lm0.rds"))
# model_AIC$fitted.values %>% length()
# ar lm
full_model_ar1 <- readRDS(file = file.path(model_path, "gls_ar1.rds"))

bru0 <- readRDS(file.path(model_path, "lm_bru0.rds"))
bru0$.args$formula
summary(bru0)
bru_ar <- readRDS(file.path(model_path, "bru_ar1.rds"))


model_AIC %>% summary()
anova(model_AIC)

# wind direction effect..

# install.packages("stargazer")
# install.packages("modelsummary")
# require(stargazer)
require(modelsummary)

modelsummary(model_AIC, output = "latex")

# qm

mod_labels <- c(
  "Generic PC",
  "Linear",
  "AR1",
  "Linear inlabru",
  "AR inlabru",
  "QM"
)
est_cols <- c(
  "norm_power_est0",
  "lm",
  "ar",
  "lm_bru",
  # "ar_bru",
  "ar_bru2", # excludes ar term (interpolates)
  "qm"
)
n <- nrow(GB_df)
names(mod_labels) <- est_cols
model_df <- GB_df %>%
  mutate(
    lm = model_AIC$fitted.values,
    ar = full_model_ar1$fitted,
    lm_bru = bru0$summary.fitted.values[1:n, "mean"],
    ar_bru = bru_ar$summary.fitted.values[1:n, "mean"],
    ar_bru2 = bru_ar$summary.fitted.values[1:n, "mean"] -
      bru_ar$summary.random$u[1:n, "mean"],
    qm = wgen_qm
  )

write_parquet(model_df, "data/calibration_df.parquet")


## Reading fitted values ####

model_df <- read_parquet("data/calibration_df.parquet")

model_df %>% head()


ModelMetrics::rmse(
  full_model_ar1$fitted,
  GB_df$norm_potential
)


df_long <- model_df %>%
  dplyr::select(tech_typ, norm_potential, all_of(est_cols)) %>%
  pivot_longer(
    cols = all_of(est_cols),
    names_to = "model",
    values_to = "estimate"
  )
require(ModelMetrics)
metrics_table <- df_long %>%
  group_by(model) %>%
  summarise(
    RMSE = rmse(actual = norm_potential, predicted = estimate),
    MAE = mae(actual = norm_potential, predicted = estimate),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(model = mod_labels[model]) %>%
  arrange(desc(RMSE))

write.csv(
  metrics_table,
  "summaries/calib_metrics.csv",
  row.names = FALSE
)

df_long %>%
  group_by(tech_typ, model) %>%
  summarise(
    RMSE = rmse(actual = norm_potential, predicted = estimate),
    MAE = mae(actual = norm_potential, predicted = estimate),
    Bias = mean(estimate - norm_potential, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(model = mod_labels[model]) %>%
  arrange(tech_typ, desc(RMSE))


est_cols <- c(
  "norm_power_est0",
  "lm",
  "ar",
  "lm_bru",
  # "ar_bru",
  "ar_bru2",
  "qm"
)
model_df %>%
  ggplot(aes(norm_potential, lm)) +
  geom_hex() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Elexon CF", y = "Calibrated ERA5 CF")
ggsave("fig/gb_lmcalib_hexbin_all.pdf", width = 6, height = 4)

model_df %>%
  ggplot(aes(norm_potential, ar)) +
  geom_hex() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Elexon CF", y = "Calibrated ERA5 CF")
ggsave("fig/gb_ar1calib_hexbin_all.pdf", width = 6, height = 4)


model_df %>%
  ggplot(aes(norm_potential, lm_bru)) +
  geom_hex() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Elexon CF", y = "Calibrated ERA5 CF")
ggsave("fig/gb_lmbrucalib_hexbin_all.pdf", width = 6, height = 4)


model_df %>%
  ggplot(aes(norm_potential, ar_bru2)) +
  geom_hex() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Elexon CF", y = "Calibrated ERA5 CF")
ggsave("fig/gb_arbrucalib_hexbin_all.pdf", width = 6, height = 4)


model_df %>%
  ggplot(aes(norm_potential, qm)) +
  geom_hex() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Elexon CF", y = "Calibrated ERA5 CF")
ggsave("fig/gb_qmcalib_hexbin_all.pdf", width = 6, height = 4)


source("aux_funct.R")
bru0$summary.random %>% names()

bru0

wind_eff <- plot.effects(bru0, "wind", show.fig = F)
plot.effects(bru0, "hour")
plot.effects(bru0, "month")
test <- plot.effects(bru0, "tech_typ", show.fig = F)
test +
  bru0$summary.random$wind
wind_eff$data %>% head()
facet_wrap()

# residual exploratory analysis ####

## Insample #####

model_df %>% names()

est_cols <- c(
  "norm_power_est0",
  "lm",
  "ar",
  "lm_bru",
  # "ar_bru",
  "ar_bru2", # excludes ar term (interpolates)
  "qm"
)
model_df <- model_df %>%
  mutate(
    across(
      any_of(est_cols),
      ~ . - norm_potential,
      .names = "{.col}_resid"
    ),
    month = month(halfHourEndTime),
    season = case_when(
      month %in% c(12, 1, 2) ~ "DJF",
      month %in% c(3, 4, 5) ~ "MAM",
      month %in% c(6, 7, 8) ~ "JJA",
      month %in% c(9, 10, 11) ~ "SON",
      TRUE ~ NA_character_
    ) %>%
      factor(levels = c("MAM", "JJA", "SON", "DJF"))
  )

# turbine type
model_df %>%
  ggplot(aes(tech_typ, norm_power_est0_resid)) +
  geom_boxplot(aes(fill = tech_typ)) +
  scale_fill_manual(values = mypalette) +
  theme(legend.position = "none") +
  labs(x = "", y = "Residuals of Generic PC")

ggsave("fig/gb_tech_genPC_resid_boxplot.pdf", width = 5, height = 4)

model_df %>%
  ggplot(aes(norm_power_est0_resid)) +
  geom_density(aes(fill = tech_typ), alpha = 0.5) +
  scale_fill_manual(values = mypalette) +
  theme(legend.position = "inside", legend.position.inside = c(0.2, 0.8)) +
  labs(x = "", y = "Residuals of Generic PC", fill = "type")

ggsave("fig/gb_tech_genPC_resid_density.pdf", width = 5, height = 4)

model_df %>%
  ggplot(aes(x = norm_potential, y = norm_power_est0_resid)) +
  # This creates the filled contour "topography"
  geom_hex(alpha = 0.8) +
  # Add contour lines for extra precision
  # geom_density_2d(color = "white", size = 0.2, alpha = 0.3) +
  scale_fill_viridis_c(
    trans = "log10",
    name = "frequency",
    limits = c(1, NA)
  ) +
  facet_wrap(~tech_typ) +
  # theme_minimal() +
  theme(legend.position = "right") +
  labs(
    x = "Normalised wind generation",
    y = "Residuals",
    fill = "Density Level"
  )
ggsave("fig/gb_tech+_genPC_resid_2Ddensity.pdf", width = 7, height = 4)

# season
# 1: Red, 2: Blue, 3: Green, 8: Brown/Dark Grey
seasonal_colors <- c(
  "DJF" = "#4DBBD5FF", # Blue (Winter)
  "MAM" = "#00A087FF", # Teal/Green (Spring)
  "JJA" = "#E64B35FF", # Red/Coral (Summer)
  "SON" = "#7E6148FF" # Brown (Autumn)
)
model_df %>%
  ggplot(aes(season, norm_power_est0_resid)) +
  geom_boxplot(aes(fill = season)) +
  facet_wrap(~tech_typ) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = seasonal_colors) +
  theme(legend.position = "none") +
  labs(x = "", y = "Residuals of Generic PC")

ggsave("fig/gb_tech+seas_genPC_resid_boxplot.pdf", width = 3.5, height = 3)

require(ggridges)
model_df %>%
  mutate(
    season = factor(season, levels = rev(c("DJF", "MAM", "JJA", "SON")))
  ) %>%
  ggplot(aes(norm_power_est0_resid, season, fill = season)) +
  # geom_density(aes(fill = season), alpha = 0.5) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, color = "white") +
  facet_wrap(~tech_typ, nrow = 2) +
  theme_ridges() +
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(-0.2, 0.6)) +
  scale_fill_manual(values = seasonal_colors) +
  theme(legend.position = "none", legend.position.inside = c(0.1, 0.8)) +
  labs(x = "", y = "Residuals of Generic PC")
ggsave("fig/gb_tech+seas_genPC_resid_ridge.pdf", width = 5, height = 4)

# month
model_df %>%
  ggplot(aes(x = norm_power_est0_resid, y = as.factor(month))) +
  geom_density_ridges(
    aes(fill = season),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  coord_cartesian(xlim = c(-0.2, 0.6)) +
  scale_fill_manual(values = seasonal_colors) +
  theme(legend.position = "none") +
  labs(x = "Residuals of Generic PC", y = "month")
ggsave("fig/gb_tech+month_genPC_resid_ridge.pdf", width = 5, height = 5)


model_df %>%
  ggplot(aes(as.factor(month), norm_power_est0_resid)) +
  geom_boxplot(aes(fill = season)) +
  facet_wrap(~tech_typ) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = seasonal_colors) +
  theme(legend.position = "none") +
  labs(x = "month", y = "Residuals of Generic PC")
ggsave("fig/gb_tech+month_genPC_resid_boxplot.pdf", width = 5, height = 4)

# hour

# Dark Blue (00:00) -> Light Blue (12:00) -> Dark Blue (23:00)
simple_blues <- c("#08306b", "#c6dbef", "#08306b")
model_df %>% names()
model_df %>%
  ggplot(aes(x = as.factor(hour), y = ws_h_wmean)) +
  # Use as.numeric(hour) to ensure the gradient mapping works
  geom_boxplot(
    aes(fill = as.numeric(hour)),
    outlier.size = 0.5,
    alpha = 0.8,
    outliers = FALSE
  ) +
  facet_wrap(~tech_typ, nrow = 2, scales = "free") +
  # Applying the gradient
  scale_fill_gradientn(colors = simple_blues) +
  # theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    # Academic journals often prefer slightly larger axis titles
    # axis.title = element_text(face = "bold")
  ) +
  labs(
    x = "Hour of Day",
    y = "Residuals of Generic PC"
  )
ggsave("fig/gb_tech+hour_genPC_resid_boxplot.pdf", width = 5, height = 4)


model_df %>%
  ggplot(aes(y = as.factor(hour), x = norm_power_est0_resid)) +
  geom_density_ridges(
    aes(fill = as.numeric(hour)),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ, nrow = 2, scales = "free") +
  scale_fill_gradientn(colors = simple_blues) +
  theme_ridges() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    # Academic journals often prefer slightly larger axis titles
    # axis.title = element_text(face = "bold")
  ) +
  labs(
    y = "Hour of Day",
    x = "Residuals of Generic PC"
  )
ggsave("fig/gb_tech+hour_genPC_resid_ridge.pdf", width = 5, height = 7)

# best model coefficients ####

# Reading wind farm level PC models ####
