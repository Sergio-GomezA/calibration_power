require(tidyverse)
source('aux_funct.R')
# read output data ####
prefix <- "elexon"
prefixfull <- "elex1H"
tols <- c(0.3, 0.4, 0.5)
batch_names <- paste0("batchsim_tol", gsub("\\.", "_", sprintf("%.2f", tols)))

elexsim_files <- lapply(
  batch_names,
  function(batch_name) {
    list.files(
      file.path("summaries", batch_name),
      pattern = sprintf("^%s_aggr_samples.*\\.%s$", prefix, "csv")
    )
  }
)

elexsim_df <- lapply(
  seq(elexsim_files),
  function(i) {
    {
      # browser()
      {
        lapply(
          elexsim_files[[i]],
          function(x) {
            read.csv(
              file.path("summaries", batch_names[i], x),
              stringsAsFactors = FALSE
            )
          }
        ) %>%
          bind_rows() %>%
          mutate(anomaly_thres = tols[i])
      }
    }
  }
) %>%
  bind_rows() # %>%
# filter(!is.na(var_samples)) # remove cases where no samples were generated

elex1H_files <- lapply(
  batch_names,
  function(batch_name) {
    list.files(
      file.path("summaries", batch_name),
      pattern = sprintf("^%s_.*\\.%s$", prefixfull, "csv")
    )
  }
)

elex1H_df <- lapply(
  seq(elex1H_files),
  function(i) {
    {
      {
        lapply(
          elex1H_files[[i]],
          function(x) {
            read.csv(
              file.path("summaries", batch_names[i], x),
              stringsAsFactors = FALSE
            )
          }
        ) %>%
          bind_rows() %>%
          mutate(anomaly_thres = tols[i])
      }
    }
  }
) %>%
  bind_rows()
# filter(!is.na(variance_samples)) # remove cases where no samples were generated

elexsim_df %>%
  group_by(anomaly_thres) %>%
  summarise(
    aggr_in_ci = mean(aggr_in_ci)
  )
elex1H_df %>%
  group_by(anomaly_thres) %>%
  summarise(
    aggr_in_ci = mean(aggr_in_ci)
  )


elexsim_df %>%
  group_by(anomaly_thres) %>%
  summarise(
    cov_noise_oos = mean(cov_noise_oos),
    aggr_in_ci = mean(aggr_in_ci)
  )
# plot figures ####

## scatter plot of posterior mean vs observed data with error bars

### elexsim

elexsim_df %>%
  ggplot(aes(x = fit, y = observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  # geom_point() +
  geom_errorbar(
    aes(ymin = q0.025, ymax = q0.975),
    color = blues9[5],
    # width = 0.1,
    alpha = 0.5
  ) +
  geom_point(aes(fit, observed, color = (as.logical(aggr_in_ci)))) +
  facet_wrap(~anomaly_thres) +
  labs(x = "Posterior mean", y = "Simulated from model fit", col = "In CI") +
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("TRUE" = blues9[7], "FALSE" = "darkred"),
  )
# ggsave(
#   file.path("fig", batch_name, sprintf("%s_aggr_cov_scatter.pdf", prefix)),
#   width = 6,
#   height = 6
# )

### elex1H model

elex1H_df %>%
  ggplot(aes(x = fit, y = observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_errorbar(
    aes(ymin = q0.025, ymax = q0.975),
    color = blues9[5],
    # width = 0.1,
    alpha = 0.5
  ) +
  geom_point(aes(fit, observed, color = as.logical(aggr_in_ci))) +
  labs(x = "Posterior mean", y = "Observed data", col = "In CI") +
  scale_color_manual(
    values = c("TRUE" = blues9[7], "FALSE" = "darkred"),
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  facet_wrap(~anomaly_thres)
# ggsave(
#   file.path("fig", batch_name, sprintf("%s_aggr_cov_scatter.pdf", prefixfull)),
#   width = 6,
#   height = 6
# )

## coverage simulations vs real data

# elexsim_df %>%
#   pivot_longer(
#     cols = c(cov_nonoise_oos, cov_noise_oos),
#     names_to = "coverage_type",
#     values_to = "coverage"
#   ) %>%
#   mutate(
#     coverage_type = factor(
#       coverage_type,
#       levels = c("cov_nonoise_oos", "cov_noise_oos"),
#       labels = c("latent", "latent + obs noise")
#     )
#   ) %>%
#   ggplot() +
#   geom_boxplot(aes(x = coverage_type, y = coverage)) +
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
#   labs(x = "Coverage type", y = "Coverage probability")

# elex1H_df %>%
#   pivot_longer(
#     cols = c(cov_nonoise_oos, cov_noise_oos),
#     names_to = "coverage_type",
#     values_to = "coverage"
#   ) %>%
#   mutate(
#     coverage_type = factor(
#       coverage_type,
#       levels = c("cov_nonoise_oos", "cov_noise_oos"),
#       labels = c("latent", "latent + obs noise")
#     )
#   ) %>%
#   ggplot() +
#   geom_boxplot(aes(x = coverage_type, y = coverage)) +
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
#   labs(x = "Coverage type", y = "Coverage probability")

elexsim_df %>%
  dplyr::select(aggr_in_ci, cov_noise_oos, anomaly_thres) %>%
  mutate(
    type = "simulation"
  ) %>%
  bind_rows(
    elex1H_df %>%
      dplyr::select(aggr_in_ci, cov_noise_oos, anomaly_thres) %>%
      mutate(
        type = "real data"
      )
  ) %>%
  mutate(
    aggr_in_ci = as.logical(aggr_in_ci),
    anomaly_thres = factor(anomaly_thres, levels = tols)
  ) %>%
  ggplot() +
  geom_boxplot(aes(x = anomaly_thres, y = cov_noise_oos)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red")


elexsim_df %>%
  dplyr::select(aggr_in_ci, cov_noise_oos, anomaly_thres) %>%
  mutate(
    type = "simulation"
  ) %>%
  bind_rows(
    elex1H_df %>%
      dplyr::select(aggr_in_ci, cov_noise_oos, anomaly_thres) %>%
      mutate(
        type = "real data"
      )
  ) %>%
  mutate(
    aggr_in_ci = as.logical(aggr_in_ci)
  ) %>%
  ggplot() +
  geom_boxplot(aes(x = aggr_in_ci, y = cov_noise_oos, fill = type)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  facet_wrap(~anomaly_thres) +
  labs(
    x = "Aggregated in CI",
    y = "Coverage probability (WF level)",
    fill = "Data type"
  ) +
  scale_fill_manual(
    values = c("simulation" = blues9[7], "real data" = blues9[5])
  ) +
  theme(legend.position = "bottom")

ggsave(
  file.path("fig", batch_name, sprintf("%s_WF_coverage_boxplot.pdf", prefix)),
  width = 6,
  height = 6
)


# coverage by normality test

# elexsim_df %>%
#   filter(!is.na(normality)) %>%
#   pivot_longer(
#     cols = c(cov_nonoise_oos, cov_noise_oos),
#     names_to = "coverage_type",
#     values_to = "coverage"
#   ) %>%
#   mutate(
#     coverage_type = factor(
#       coverage_type,
#       levels = c("cov_nonoise_oos", "cov_noise_oos"),
#       labels = c("latent", "latent + obs noise")
#     )
#   ) %>%
#   ggplot() +
#   geom_boxplot(aes(x = coverage_type, y = coverage)) +
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
#   facet_wrap(~normality) +
#   labs(x = "Coverage type", y = "Coverage probability")

# elex1H_df %>%
#   filter(!is.na(normality)) %>%
#   pivot_longer(
#     cols = c(cov_nonoise_oos, cov_noise_oos),
#     names_to = "coverage_type",
#     values_to = "coverage"
#   ) %>%
#   mutate(
#     coverage_type = factor(
#       coverage_type,
#       levels = c("cov_nonoise_oos", "cov_noise_oos"),
#       labels = c("latent", "latent + obs noise")
#     )
#   ) %>%
#   ggplot() +
#   geom_boxplot(aes(x = coverage_type, y = coverage)) +
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
#   facet_wrap(~normality) +
#   labs(x = "Coverage type", y = "Coverage probability")

## variance of samples vs variance of observed data

elexsim_df %>%
  ggplot(aes(x = variance_samples, y = variance_observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(col = blues9[7]) +
  facet_wrap(~anomaly_thres) +
  labs(x = "Variance of samples", y = "Variance of observed data")

elex1H_df %>%
  ggplot(aes(x = variance_samples, y = variance_observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(
    aes(col = as.logical(aggr_in_ci)),
    # col = blues9[7]
    alpha = 1
  ) +
  facet_wrap(~anomaly_thres) +
  labs(
    x = "Variance of samples",
    y = "Variance of observed data",
    col = "Aggregated in CI"
  ) +
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("TRUE" = blues9[7], "FALSE" = "darkred"),
  )
ggsave(
  file.path(
    "fig",
    batch_name,
    sprintf("%s_WF_variance_scatter.pdf", prefixfull)
  ),
  width = 6,
  height = 6
)
