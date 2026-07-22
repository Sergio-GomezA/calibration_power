# model labels catalog
mod_labels <- c(
  "Generic PC",
  "Linear model",
  "GB LM",
  "QM",
  "ST model fine",
  "ST model coarse",
  "LM+hour model",
  "AR1 model",
  "AR2 model"
)
est_cols <- c(
  "norm_power_est0",
  "lm",
  "agg_lm",
  "qm",
  "st0_m1",
  "st0_m2",
  "spde1d",
  "ar1",
  "ar2"
)

mode_code_prefix <- c(
  NA,
  "lm_model_aic0_",
  "ts_bru0_ar1_",
  "ts_bru0_ar2_",
  "ts_bru0_1DSPDE_",
  "st_bru0_coarse_",
  "st_bru0_very_coarse_",
  "qm_model_",
  "lm_model_aic0_agg_"
)
