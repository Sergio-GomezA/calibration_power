# read output data ####

elexsim_files <- list.files(
  file.path("summaries", batch_name),
  pattern = sprintf("^%s_aggr.*\\.%s$", prefix, "csv")
)

elexsim_df <- lapply(
  elexsim_files,
  function(x) {
    read.csv(
      file.path("summaries", batch_name, x),
      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows() %>%
  filter(!is.na(var_samples)) # remove cases where no samples were generated

elex1H_files <- list.files(
  file.path("summaries", batch_name),
  pattern = sprintf("^%s_.*\\.%s$", prefixfull, "csv")
)

elex1H_df <- lapply(
  elex1H_files,
  function(x) {
    read.csv(
      file.path("summaries", batch_name, x),
      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows() %>%
  filter(!is.na(variance_samples)) # remove cases where no samples were generated


elexsim_df %>%
  summarise(
    aggr_in_ci = mean(aggr_in_ci)
  )
elex1H_df %>%
  summarise(
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
  geom_point(aes(fit, observed, color = factor(aggr_in_ci))) +
  labs(x = "Posterior mean", y = "Simulated from model fit")

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
  geom_point(aes(fit, observed, color = factor(aggr_in_ci))) +
  labs(x = "Posterior mean", y = "Observed data")


## coverage simulations vs real data

elexsim_df %>%
  pivot_longer(
    cols = c(cov_nonoise_oos, cov_noise_oos),
    names_to = "coverage_type",
    values_to = "coverage"
  ) %>%
  mutate(
    coverage_type = factor(
      coverage_type,
      levels = c("cov_nonoise_oos", "cov_noise_oos"),
      labels = c("latent", "latent + obs noise")
    )
  ) %>%
  ggplot() +
  geom_boxplot(aes(x = coverage_type, y = coverage)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(x = "Coverage type", y = "Coverage probability")


elex1H_df %>%
  pivot_longer(
    cols = c(cov_nonoise_oos, cov_noise_oos),
    names_to = "coverage_type",
    values_to = "coverage"
  ) %>%
  mutate(
    coverage_type = factor(
      coverage_type,
      levels = c("cov_nonoise_oos", "cov_noise_oos"),
      labels = c("latent", "latent + obs noise")
    )
  ) %>%
  ggplot() +
  geom_boxplot(aes(x = coverage_type, y = coverage)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(x = "Coverage type", y = "Coverage probability")


# coverage by normality test

elexsim_df %>%
  filter(!is.na(normality)) %>%
  pivot_longer(
    cols = c(cov_nonoise_oos, cov_noise_oos),
    names_to = "coverage_type",
    values_to = "coverage"
  ) %>%
  mutate(
    coverage_type = factor(
      coverage_type,
      levels = c("cov_nonoise_oos", "cov_noise_oos"),
      labels = c("latent", "latent + obs noise")
    )
  ) %>%
  ggplot() +
  geom_boxplot(aes(x = coverage_type, y = coverage)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  facet_wrap(~normality) +
  labs(x = "Coverage type", y = "Coverage probability")

elex1H_df %>%
  filter(!is.na(normality)) %>%
  pivot_longer(
    cols = c(cov_nonoise_oos, cov_noise_oos),
    names_to = "coverage_type",
    values_to = "coverage"
  ) %>%
  mutate(
    coverage_type = factor(
      coverage_type,
      levels = c("cov_nonoise_oos", "cov_noise_oos"),
      labels = c("latent", "latent + obs noise")
    )
  ) %>%
  ggplot() +
  geom_boxplot(aes(x = coverage_type, y = coverage)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  facet_wrap(~normality) +
  labs(x = "Coverage type", y = "Coverage probability")
