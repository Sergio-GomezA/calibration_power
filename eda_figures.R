fig_path <- "fig"


gb_daily_df %>%
  ggplot(aes(x = tech_typ, y = err)) +
  geom_boxplot()


gb_daily_df %>%
  ggplot(aes(x = nao_group, y = err, group = nao_group)) +
  geom_boxplot() +
  facet_wrap(~tech_typ) +
  theme_bw()


gb_daily_df %>%
  ggplot(aes(x = err, y = nao_group, group = nao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "NAO index") +
  theme(legend.position = "none")

gb_daily_df %>%
  ggplot(aes(x = norm_potential, y = nao_group, group = nao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "NAO index") +
  theme(legend.position = "none")


gb_daily_df %>%
  ggplot(aes(x = err, y = nao_group, group = nao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "NAO index") +
  theme(legend.position = "none")


## AO
gb_daily_df %>%
  ggplot(aes(x = err, y = ao_group, group = ao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "AO index") +
  theme(legend.position = "none")

gb_daily_df %>%
  ggplot(aes(x = norm_potential, y = ao_group, group = ao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "AO index") +
  theme(legend.position = "none")


gb_daily_df %>%
  ggplot(aes(x = err, y = ao_group, group = ao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "AO index") +
  theme(legend.position = "none")

### Wind, NAO, AO, EA, SCAN

library(GGally)


vars <- gb_monthly_df[, c("ws_h_wmean", "nao", "ao", "ea", "scan")]

ggpairs(
  vars,
  diag = list(continuous = wrap("densityDiag")),
  lower = list(continuous = wrap("points", alpha = 0.4, size = 0.8)),
  upper = list(continuous = wrap("cor", size = 3))
)

## circulation patterns and generation ------------------------------
# NAO
gb_monthly_df %>%
  ggplot(aes(x = norm_potential, y = nao_group, group = nao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "wind generation", y = "NAO index") +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "nao_generation.pdf"
  ),
  width = 6,
  height = 4
)

# AO
gb_monthly_df %>%
  ggplot(aes(x = norm_potential, y = ao_group, group = ao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "wind generation", y = "AO index") +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "ao_generation.pdf"
  ),
  width = 6,
  height = 4
)
# EA
gb_monthly_df %>%
  ggplot(aes(x = norm_potential, y = ea_group, group = ea_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "wind generation", y = "EA index") +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "ea_generation.pdf"
  ),
  width = 6,
  height = 4
)
# SCAN
gb_monthly_df %>%
  ggplot(aes(x = norm_potential, y = scan_group, group = scan_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "wind generation", y = "SCAN index") +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "scan_generation.pdf"
  ),
  width = 6,
  height = 4
)

## circulation patterns and error ------------------------------
# NAO
gb_monthly_df %>%
  ggplot(aes(x = err, y = nao_group, group = nao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "NAO index") +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "nao_error.pdf"
  ),
  width = 6,
  height = 4
)

# AO
gb_monthly_df %>%
  ggplot(aes(x = err, y = ao_group, group = ao_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "AO index") +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "ao_error.pdf"
  ),
  width = 6,
  height = 4
)
# EA

gb_monthly_df %>%
  ggplot(aes(x = err, y = ea_group, group = ea_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "EA index") +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "ea_error.pdf"
  ),
  width = 6,
  height = 4
)
# SCAN
gb_monthly_df %>%
  ggplot(aes(x = err, y = scan_group, group = scan_group)) +
  geom_density_ridges(
    aes(fill = tech_typ),
    alpha = 0.7,
    scale = 1.2,
    color = "white"
  ) +
  facet_wrap(~tech_typ) +
  theme_ridges() +
  scale_fill_manual(values = mypalette) +
  labs(x = "normalised error", y = "SCAN index") +
  theme(legend.position = "none")
ggsave(
  file.path(
    fig_path,
    "scan_error.pdf"
  ),
  width = 6,
  height = 4
)
