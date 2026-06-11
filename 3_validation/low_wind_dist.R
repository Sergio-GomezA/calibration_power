# low wind events distribution

model_df0 <- st_read(sprintf(
  "data/calibration_df_%s_%s.%s",
  mesh_label,
  d0_tag,
  extension
))
