starttime <- Sys.time()


## 2.0 bru lm model ####

mod_starttime <- Sys.time()

mod_tag <- "lm"
components0 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  techno(tech_typ, model = "iid") + # random intercept by tech_typ
  # norm_power_est0 +
  slope(
    tech_typ,
    model = "iid",
    weights = norm_power_est0
  ) +
  d_coast(
    d_coast_group,
    model = "rw2",
    constr = TRUE
  ) + # smooth correction distance to coast
  elev(
    elev_group,
    model = "rw2",
    constr = TRUE
  ) + # smooth correction elevation
  wind(ws_group, model = "rw2", replicate = tech_typ, constr = TRUE) # smooth correction wind

model_code <- sprintf("ts_bru0_%s_%s.rds", mod_tag, d0_tag)
model_fname <- file.path(
  model_path,
  model_code
)

if (!file.exists(model_fname) || override_objects) {
  cat(
    "-------------------------------------------------------------------------------------------------\n"
  )
  cat("Fitting bru lm model\n")
  cat(
    "-------------------------------------------------------------------------------------------------\n"
  )
  brulm <- bru(
    components = components0,
    formula = norm_potential ~ Intercept +
      techno +
      slope +
      # power_correction +
      d_coast +
      elev +
      wind,
    family = "gaussian",
    control.family = list(
      hyper = list(
        prec = list(
          prior = "pc.prec",
          param = c(50, 0.05)
        )
      )
    ),
    data = wf_dt,
    options = base_bru_options
  )

  scores_df[[model_code]] <- extract_score_model(brulm)
  pit_list[[model_code]] <- extract_pit_model(brulm, wf_dt)

  if (save_models) {
    saveRDS(
      brulm,
      file = model_fname
    )
  } else {
    model_list[[model_code]] <- brulm
  }
} else {
  cat("Loading existing lm model\n")
  brulm <- readRDS(model_fname)
}

summary(brulm)

effect_names <- names(brulm$summary.random)
excluded_effects <- c("u")
effect_names <- setdiff(effect_names, excluded_effects)
for (effect in effect_names) {
  if (effect %in% c("wind")) {
    n_repl <- 2
    repl_names <- c("Offshore", "Onshore")
  } else {
    n_repl <- 1
    repl_names <- NULL
  }
  # browser()
  plot.effects(
    brulm,
    effect,
    n.replicate = n_repl,
    replicate_names = repl_names,
    show.plot = TRUE
  )
  ggsave(
    sprintf(
      "%s/%s/fig/fit/%s_effect_%s_%s.pdf",
      output_path,
      batch_name,
      effect,
      mod_tag,
      d0_tag
    ),
    width = 6,
    height = 4
  )
}

plot.hyper.dens(brulm)
ggsave(
  sprintf(
    "%s/%s/fig/fit/hyperparameters_%s_%s.pdf",
    output_path,
    batch_name,
    mod_tag,
    d0_tag
  ),
  width = 6,
  height = 4
)
mod_endtime <- Sys.time()
f_time <- difftime(mod_endtime, mod_starttime, units = "auto")
cat(
  "Model fitting took: ",
  round(f_time, 2),
  " ",
  attr(f_time, "units"),
  "\n"
)


## 2.1 AR1 temporal model ####
mod_starttime <- Sys.time()
ar_tag <- "ar1"
components0 <- ~ Intercept(1, prec.linear = exp(-7)) + # latent intercept
  techno(tech_typ, model = "iid") + # random intercept by tech_typ
  # norm_power_est0 +
  # err_l1 +
  # err_l2 +
  # dest +
  slope(
    tech_typ,
    model = "iid",
    weights = norm_power_est0
  ) +
  # power_correction(
  #   pow_group,
  #   model = "rw2",
  #   # replicate = tech_typ,
  #   constr = TRUE
  # ) + # smooth correction power
  d_coast(
    d_coast_group,
    model = "rw2",
    constr = TRUE
  ) + # smooth correction distance to coast
  elev(
    elev_group,
    model = "rw2",
    constr = TRUE
  ) + # smooth correction elevation
  wind(ws_group, model = "rw2", replicate = tech_typ, constr = TRUE) + # smooth correction wind
  u(
    t,
    model = "ar1",
    replicate = tech_typ,
    hyper = list(
      rho = list(
        prior = "pc.cor1",
        param = c(0.6, 0.5)
      ),
      # prec = list(
      #   prior = "pc.prec",
      #   param = c(50, 0.05)
      # ),
      prec = list(initial = prec_init, fixed = fixed_ucomp)
    )
  )


model_code <- sprintf("ts_bru0_%s_%s.rds", ar_tag, d0_tag)
model_fname <- file.path(
  model_path,
  model_code
)

if (!file.exists(model_fname) || override_objects) {
  cat(
    "-------------------------------------------------------------------------------------------------\n"
  )
  cat("Fitting ar1 model\n")
  cat(
    "-------------------------------------------------------------------------------------------------\n"
  )
  bruar1 <- bru(
    components = components0,
    formula = norm_potential ~ Intercept +
      techno +
      # err_l1 +
      # err_l2 +
      # dest +
      slope +
      # power_correction +
      d_coast +
      elev +
      wind +
      u,
    family = "gaussian",
    control.family = list(
      hyper = list(
        prec = list(
          prior = "pc.prec",
          param = c(50, 0.05)
        )
      )
    ),
    data = wf_dt,
    options = base_bru_options
  )

  scores_df[[model_code]] <- extract_score_model(bruar1)
  pit_list[[model_code]] <- extract_pit_model(bruar1)

  if (save_models) {
    saveRDS(
      bruar1,
      file = model_fname
    )
  } else {
    model_list[[model_code]] <- bruar1
  }
} else {
  cat("Loading existing ar1 model\n")
  bruar1 <- readRDS(model_fname)
}

summary(bruar1)

effect_names <- names(bruar1$summary.random)
excluded_effects <- c("u")
effect_names <- setdiff(effect_names, excluded_effects)
for (effect in effect_names) {
  if (effect == "wind") {
    n_repl <- 2
    repl_names <- c("Offshore", "Onshore")
  } else {
    n_repl <- 1
    repl_names <- NULL
  }
  plot.effects(
    bruar1,
    effect,
    n.replicate = n_repl,
    replicate_names = repl_names,
    show.plot = TRUE
  )
  ggsave(
    sprintf(
      "%s/%s/fig/fit/%s_effect_%s_%s.pdf",
      output_path,
      batch_name,
      effect,
      ar_tag,
      d0_tag
    ),
    width = 6,
    height = 4
  )
}

plot.hyper.dens(bruar1)
ggsave(
  sprintf(
    "%s/%s/fig/fit/hyperparameters_%s_%s.pdf",
    output_path,
    batch_name,
    ar_tag,
    d0_tag
  ),
  width = 6,
  height = 4
)
mod_endtime <- Sys.time()
f_time <- difftime(mod_endtime, mod_starttime, units = "auto")
cat(
  "Model fitting took: ",
  round(f_time, 2),
  " ",
  attr(f_time, "units"),
  "\n"
)


endtime <- Sys.time()

cat(
  "-------------------------------------------------------------------------------------------------\n"
)
cat("Finished model fitting process\n")
cat(
  "-------------------------------------------------------------------------------------------------\n"
)


timediff <- difftime(endtime, starttime, units = "auto")

cat(sprintf(
  "Time taken for model fitting: %.2f %s\n",
  timediff,
  units(timediff)
))
