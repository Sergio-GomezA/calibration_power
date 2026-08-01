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
  "AR2 model",
  "LM bru model",
  "LM beta model",
  "LM t model"
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
  "ar2",
  "lm_bru",
  "lm_beta",
  "lm_t"
)

mode_code_prefix <- c(
  NA,
  "lm_model_aic0_",
  "gblm_model_aic0_",
  "qm_model_",
  "st_bru0_fine_",
  "st_bru0_coarse_",
  "st_bru0_very_coarse_",
  "ts_bru0_1DSPDE_",
  "ts_bru0_ar1_",
  "ts_bru0_ar2_",
  "ts_bru0_lm_",
  "ts_bru0_lmbeta_",
  "ts_bru0_lmt_"
)

dplyr::mutate(
  data.frame(
    mod_labels = mod_labels,
    est_cols = est_cols,
    mode_code_prefix = mode_code_prefix
  ),
  family = dplyr::case_when(
    est_cols == "lm_beta" ~ "beta",
    est_cols == "lm_t" ~ "t",
    TRUE ~ "gaussian"
  )
) -> model_catalog

write.csv(model_catalog, "data/model_catalog.csv", row.names = FALSE)
