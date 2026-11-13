starttime <- Sys.time()

if (model_id == 1) {
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
}

if (model_id == 2) {
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
    brulm <- bru(
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

    scores_df[[model_code]] <- extract_score_model(brulm)
    pit_list[[model_code]] <- extract_pit_model(brulm)

    if (save_models) {
      saveRDS(
        brulm,
        file = model_fname
      )
    } else {
      model_list[[model_code]] <- brulm
    }
  } else {
    cat("Loading existing ar1 model\n")
    brulm <- readRDS(model_fname)
  }

  summary(brulm)

  effect_names <- names(brulm$summary.random)
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
        ar_tag,
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
}

# if (model_id == 3) {
## 2.4 GB lm version #####
mod_starttime <- Sys.time()
model_code <- sprintf("gblm_model_aic0_%s.rds", d0_tag)

if (!file.exists(file.path(model_path, model_code)) || override_objects) {
  cat(
    "-------------------------------------------------------------------------------------------------\n"
  )
  cat("Fitting GB aggregated LM\n")
  cat(
    "-------------------------------------------------------------------------------------------------\n"
  )

  # samp_gb <- GB_df %>%
  #   filter(date >= d0 - n.days.before.light, date <= d0 + n.days - 1)

  base_model_agg <- lm(
    norm_potential ~ norm_power_est0,
    data = samp_gb
  )

  full_model0_agg <- lm(
    norm_potential ~
      tech_typ +
      tech_typ *
        norm_power_est0 +
      # err_l1 +
      # err_l2 +
      dest +
      dest_l1 +
      dest_l2 +
      # norm_power_est0 * month +
      # hour +
      # dist_coast * tech_typ +
      # elevation * tech_typ +
      # dist_coast:tech_typ +
      # elevation:tech_typ +
      # tech_typ * poly(dist_coast, 2) +
      # tech_typ * poly(elevation, 3) +
      tech_typ * poly(ws_h_wmean, 3),
    data = samp_gb
  )
  summary(full_model0_agg)

  model_AIC0_agg <- step(
    base_model_agg,
    scope = list(lower = base_model_agg, upper = full_model0_agg),
    # steps = 5,
    k = 2,
    trace = 0
  )
  summary(model_AIC0_agg)
  scores_df[[model_code]] <- data.frame(
    AIC = AIC(model_AIC0_agg),
    BIC = BIC(model_AIC0_agg),
    logLik = as.numeric(logLik(model_AIC0_agg)),
    deviance = deviance(model_AIC0_agg),
    R2 = summary(model_AIC0_agg)$r.squared
  )
  if (save_models) {
    saveRDS(
      model_AIC0_agg,
      file.path(model_path, model_code)
    )
  } else {
    model_list[[model_code]] <- model_AIC0_agg
  }
} else {
  cat("Loading existing GB aggregated LM\n")
  model_AIC0_agg <- readRDS(file.path(
    model_path,
    model_code
  ))
}
mod_endtime <- Sys.time()
f_time <- difftime(mod_endtime, mod_starttime, units = "auto")
cat(
  "Model fitting took: ",
  round(f_time, 2),
  " ",
  attr(f_time, "units"),
  "\n"
)
# }

## 2.5 QM version ####
# if (model_id == 4) {
mod_starttime <- Sys.time()
qm_fname <- sprintf("qm_model_%s.rds", d0_tag)

if (!file.exists(file.path(model_path, qm_fname)) || override_objects) {
  cat(
    "-------------------------------------------------------------------------------------------------\n"
  )
  cat("Fitting quantile mapping model\n")
  cat(
    "-------------------------------------------------------------------------------------------------\n"
  )
  qqmod <- fitQmap(
    obs = wf_dt %>% filter(!anomaly) %>% pull(norm_potential),
    mod = wf_dt %>% filter(!anomaly) %>% pull(norm_power_est0),
    method = "QUANT"
  )

  if (save_models) {
    saveRDS(
      qqmod,
      file.path(model_path, qm_fname)
    )
  } else {
    model_list[[qm_fname]] <- qqmod
  }
} else {
  cat("Loading existing quantile mapping model\n")
  qqmod <- readRDS(file.path(model_path, qm_fname))
}

wgen_qm <- with(
  wf_dt,
  doQmapQUANT(norm_power_est0, qqmod, type = "linear")
)
scores_df[[qm_fname]] <- data.frame(
  R2 = cor(wgen_qm, wf_dt$norm_potential, use = "complete.obs")^2
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
# }
cat(
  "-------------------------------------------------------------------------------------------------\n"
)
cat("Finished model fitting process\n")
cat(
  "-------------------------------------------------------------------------------------------------\n"
)


## 3.1 Model scores object saved ####

library(DBI)
library(duckdb)
library(dplyr)
library(purrr)

# Prepare data
scores_df <- list_rbind(scores_df, names_to = "model")

pit_df <- list_rbind(pit_list, names_to = "model_code") %>%
  mutate(
    code = sub("_[0-9]{6}\\.rds$", "", model_code)
  )

# Database location
db_path <- file.path(
  output_path,
  batch_name,
  "summaries/fit",
  sprintf("model_results_%s.duckdb", d0_tag)
)

dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)

# Connect
con <- dbConnect(duckdb(), dbdir = db_path)

# Append if table exists; create otherwise
dbWriteTable(
  con,
  "scores",
  scores_df,
  append = dbExistsTable(con, "scores")
)

dbWriteTable(
  con,
  "pit",
  pit_df,
  append = dbExistsTable(con, "pit")
)

dbDisconnect(con)

## fitted values df ####

model_df0 <- wf_df_frag %>%
  ungroup() %>%
  mutate(
    # lm = predict(model_AIC0, newdata = .),
    # ar1 = bruar1$summary.fitted.values[1:n, "mean"],
    # ar2 = bruar2$summary.fitted.values[1:n, "mean"],
    # spde1d = bru1d$summary.fitted.values[1:n, "mean"],
    # st0_m2 = bru0$summary.fitted.values[1:n, "mean"],
    qm = wgen_qm,
    agg_lm = predict(model_AIC0_agg, newdata = wf_df_frag),
    lm_bru = brulm$summary.fitted.values[1:n, "mean"],
    # lm2 = brulm2$summary.fitted.values[1:n, "mean"],
    # lm_beta = brulmbeta$summary.fitted.values[1:n, "mean"]
    # lm_t = brulmt$summary.fitted.values[1:n, "mean"]
  )

model_df_fname <- sprintf(
  "%s/%s/data/fit/calibration_df_%s.%s",
  output_path,
  batch_name,
  d0_tag,
  extension
)

if (extension == "geojson" & file.exists(model_df_fname)) {
  cat("GeoJSON file already exists. Overwriting...\n")
  file.remove(model_df_fname)
}

if (extension != "rds") {
  st_write(
    model_df0,
    model_df_fname,
    driver = driver,
    append = FALSE,
  )
} else {
  saveRDS(
    model_df0,
    model_df_fname
  )
}


timediff <- difftime(endtime, starttime, units = "auto")

cat(sprintf(
  "Time taken for model fitting: %.2f %s\n",
  timediff,
  units(timediff)
))
