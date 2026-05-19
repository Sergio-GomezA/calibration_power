### wind speed #####
wf_df_frag %>%
  ggplot() +
  geom_line(aes(time, ws_h, group = site_name), alpha = 0.5, col = "gray50") +
  geom_line(
    data = wf_df_frag %>%
      group_by(time) %>%
      summarise(ws_h = sum(ws_h * capacity) / sum(capacity), .groups = "drop"),
    aes(time, ws_h, col = "capacity weighted avg."),
    lwd = 1
  ) +
  geom_line(
    data = wf_df_frag %>%
      group_by(time) %>%
      summarise(ws_h = mean(ws_h), .groups = "drop"),
    aes(time, ws_h, col = "simple avg."),
    lwd = 1
  ) +
  theme_minimal() +
  scale_x_datetime(date_labels = "%H:%M") +
  labs(
    title = sprintf("Wind Speed Time Series %s", d0),
    x = "Time",
    y = "Wind Speed (m/s)",
    col = ""
  ) +
  scale_color_manual(
    values = c("capacity weighted avg." = "blue", "simple avg." = "red")
  ) +
  theme(legend.position = "bottom")
ggsave(
  sprintf("fig/wind_speed_time_series_%s.pdf", d0_tag),
  width = 6,
  height = 4
)

### power estimate and potential ####
wf_df_frag %>%
  ggplot() +
  geom_line(
    aes(time, norm_power_est0, group = site_name),
    alpha = 0.5,
    col = "gray50"
  ) +
  geom_line(
    data = wf_df_frag %>%
      group_by(time) %>%
      summarise(
        power = sum(norm_power_est0 * capacity) / sum(capacity),
        .groups = "drop"
      ),
    aes(time, power, col = "capacity weighted avg."),
    lwd = 1
  ) +
  geom_line(
    data = wf_df_frag %>%
      group_by(time) %>%
      summarise(power = mean(norm_power_est0), .groups = "drop"),
    aes(time, power, col = "simple avg."),
    lwd = 1
  ) +
  theme_minimal() +
  scale_x_datetime(date_labels = "%H:%M") +
  labs(
    title = sprintf("ERA5 + Generic PC Time Series %s", d0),
    x = "Time",
    y = "Generation (% of capacity)",
    col = ""
  ) +
  scale_color_manual(
    values = c("capacity weighted avg." = "blue", "simple avg." = "red")
  ) +
  theme(legend.position = "bottom")
ggsave(
  sprintf("fig/era5_generic_pc_time_series_%s.pdf", d0_tag),
  width = 6,
  height = 4
)

wf_df_frag %>%
  ggplot() +
  geom_line(
    aes(time, norm_potential, group = site_name),
    alpha = 0.5,
    col = "gray50"
  ) +
  geom_line(
    data = wf_df_frag %>%
      group_by(time) %>%
      summarise(
        power = sum(norm_potential * capacity) / sum(capacity),
        .groups = "drop"
      ),
    aes(time, power, col = "capacity weighted avg."),
    lwd = 1
  ) +
  geom_line(
    data = wf_df_frag %>%
      group_by(time) %>%
      summarise(power = mean(norm_potential), .groups = "drop"),
    aes(time, power, col = "simple avg."),
    lwd = 1
  ) +
  theme_minimal() +
  scale_x_datetime(date_labels = "%H:%M") +
  labs(
    title = sprintf("Observed power Time Series %s", d0),
    x = "Time",
    y = "Generation (% of capacity)",
    col = ""
  ) +
  scale_color_manual(
    values = c("capacity weighted avg." = "blue", "simple avg." = "red")
  ) +
  theme(legend.position = "bottom")
ggsave(
  sprintf("fig/observed_power_time_series_%s.pdf", d0_tag),
  width = 6,
  height = 4
)
### error ####
wf_df_frag %>%
  ggplot() +
  geom_line(
    aes(time, error0, group = site_name),
    alpha = 0.5,
    col = "gray50"
  ) +
  geom_line(
    data = wf_df_frag %>%
      group_by(time) %>%
      summarise(
        error0 = sum(error0 * capacity) / sum(capacity),
        .groups = "drop"
      ),
    aes(time, error0, col = "capacity weighted avg."),
    lwd = 1
  ) +
  geom_line(
    data = wf_df_frag %>%
      group_by(time) %>%
      summarise(error0 = mean(error0), .groups = "drop"),
    aes(time, error0, col = "simple avg."),
    lwd = 1
  ) +
  theme_minimal() +
  scale_x_datetime(date_labels = "%H:%M") +
  labs(
    title = sprintf("Error Time Series %s", d0),
    x = "Time",
    y = "Error",
    col = ""
  ) +
  scale_color_manual(
    values = c("capacity weighted avg." = "blue", "simple avg." = "red")
  ) +
  theme(legend.position = "bottom")
ggsave(
  sprintf("fig/error_time_series_%s.pdf", d0_tag),
  width = 6,
  height = 4
)

### tech type versions ####
wf_df_frag %>%
  ggplot() +
  geom_line(
    aes(time, norm_power_est0, group = site_name),
    alpha = 0.5,
    col = "gray50"
  ) +
  geom_line(
    data = wf_df_frag %>%
      group_by(time, tech_typ) %>%
      summarise(
        power = sum(norm_power_est0 * capacity) / sum(capacity),
        .groups = "drop"
      ),
    aes(time, power, col = "capacity weighted avg."),
    lwd = 1
  ) +
  geom_line(
    data = wf_df_frag %>%
      group_by(time, tech_typ) %>%
      summarise(power = mean(norm_power_est0), .groups = "drop"),
    aes(time, power, col = "simple avg."),
    lwd = 1
  ) +
  facet_wrap(~tech_typ) +
  theme_minimal() +
  scale_x_datetime(date_labels = "%H:%M") +
  labs(
    title = sprintf("ERA5 + Generic PC Time Series %s", d0),
    x = "Time",
    y = "Generation (% of capacity)",
    col = ""
  ) +
  scale_color_manual(
    values = c("capacity weighted avg." = "blue", "simple avg." = "red")
  ) +
  theme(legend.position = "bottom")
ggsave(
  sprintf("fig/era5_generic_pc_time_series_tech_typ_%s.pdf", d0_tag),
  width = 6,
  height = 4
)

wf_df_frag %>%
  ggplot() +
  geom_line(
    aes(time, norm_potential, group = site_name),
    alpha = 0.5,
    col = "gray50"
  ) +
  geom_line(
    data = wf_df_frag %>%
      group_by(time, tech_typ) %>%
      summarise(
        power = sum(norm_potential * capacity) / sum(capacity),
        .groups = "drop"
      ),
    aes(time, power, col = "capacity weighted avg."),
    lwd = 1
  ) +
  geom_line(
    data = wf_df_frag %>%
      group_by(time, tech_typ) %>%
      summarise(power = mean(norm_potential), .groups = "drop"),
    aes(time, power, col = "simple avg."),
    lwd = 1
  ) +
  facet_wrap(~tech_typ) +
  theme_minimal() +
  scale_x_datetime(date_labels = "%H:%M") +
  labs(
    title = sprintf("Observed power Time Series %s", d0),
    x = "Time",
    y = "Generation (% of capacity)",
    col = ""
  ) +
  scale_color_manual(
    values = c("capacity weighted avg." = "blue", "simple avg." = "red")
  ) +
  theme(legend.position = "bottom")
ggsave(
  sprintf("fig/observed_power_time_series_tech_typ_%s.pdf", d0_tag),
  width = 6,
  height = 4
)


wf_df_frag %>%
  ggplot() +
  geom_line(
    aes(time, error0, group = site_name),
    alpha = 0.5,
    col = "gray50"
  ) +
  geom_line(
    data = wf_df_frag %>%
      group_by(time, tech_typ) %>%
      summarise(
        error0 = sum(error0 * capacity) / sum(capacity),
        .groups = "drop"
      ),
    aes(time, error0, col = "capacity weighted avg."),
    lwd = 1
  ) +
  geom_line(
    data = wf_df_frag %>%
      group_by(time, tech_typ) %>%
      summarise(error0 = mean(error0), .groups = "drop"),
    aes(time, error0, col = "simple avg."),
    lwd = 1
  ) +
  facet_wrap(~tech_typ) +
  theme_minimal() +
  scale_x_datetime(date_labels = "%H:%M") +
  labs(
    title = sprintf("Error Time Series %s", d0),
    x = "Time",
    y = "Error",
    col = ""
  ) +
  scale_color_manual(
    values = c("capacity weighted avg." = "blue", "simple avg." = "red")
  ) +
  theme(legend.position = "bottom")
ggsave(
  sprintf("fig/error_time_series_tech_typ_%s.pdf", d0_tag),
  width = 6,
  height = 4
)
