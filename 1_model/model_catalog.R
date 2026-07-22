# model labels catalog
mod_labels <- c(
  "Generic PC",
  "Linear model",
  "GB LM",
  "QM",
  "ST model fine",
  "ST model coarse",
  "ST model coarser",
  "LM+hour model",
  "AR1 model",
  "AR2 model"
)
est_cols <- c(
  "norm_power_est0",
  "lm",
  "agg_lm",
  "qm",
  "st0_m0",
  "st0_m1",
  "st0_m2",
  "spde1d",
  "ar1",
  "ar2"
)

mode_code_prefix <- c(
  NA,
  "lm_model_aic0_",
  "lm_model_aic0_agg_",
  "qm_model_",
  "st_bru0_fine_",
  "st_bru0_coarse_",
  "st_bru0_very_coarse_",
  "ts_bru0_1DSPDE_",
  "ts_bru0_ar1_",
  "ts_bru0_ar2_"
)

data.frame(
  mod_labels = mod_labels,
  est_cols = est_cols,
  mode_code_prefix = mode_code_prefix
) -> model_catalog

write.csv(model_catalog, "data/model_catalog.csv", row.names = FALSE)
