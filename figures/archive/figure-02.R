# Effects of N-3 supplementation on fatty acid composition
# This file saves two figures:
#   Figure 2: Omega 3 Index % with group contrasts
#   Figure 3: Fatty acid compositionsc with group contrasts

# Source figure options, packages, and main data prep.
source("./R/fig-opt-fun.R")
source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")

# Settings
refit.models <- FALSE


fat_dat <- AOData::omega3 |>
  inner_join(idkeys) |>

  filter(time %in% c("T1", "T2", "T4"), participant %in% total_participants) |>
  select(participant, time, treat, C16:C226n3) |>
  pivot_longer(cols = C16:C226n3, names_to = "fa", values_to = "value") |>
  mutate(
    .by = fa,
    value_std = log(value),
    value_std = (value_std - mean(value_std)) / sd(value_std)
  ) |>
  mutate(FA = dense_rank(fa), TX = dense_rank(treat))


# Fatty acid summarise
fat_sum <- fat_dat |>
  summarise(.by = c(FA, fa), m = mean(log(value)), s = sd(log(value)))


# Check n
fat_n <- fat_dat |>
  distinct(time, treat, participant) |>
  summarise(.by = c(time, treat), n = n())
# Prepare data set with complete cases at baseline, imputation of missing values
# post-treatment.
fat_dat_compl <- AOData::omega3 |>
  inner_join(idkeys) |>
  filter(time %in% c("T1", "T2", "T4"), participant %in% total_participants) |>

  # The omega-3 index is modelled as an additional "fatty acid" (FA = 12) in
  # m4 (see 02-fa.qmd)
  select(participant, time, treat, C16:C226n3, omega3index) |>

  pivot_longer(
    cols = c(C16:C226n3, omega3index),
    names_to = "fa",
    values_to = "value"
  ) |>

  pivot_wider(names_from = time, values_from = value) |>

  # Two participants have no data at T1
  filter(!is.na(T1)) |>
  pivot_longer(names_to = "time", values_to = "value", cols = T1:T4) |>
  mutate(
    .by = fa,
    value_std = log(value),
    value_std = (value_std - mean(value_std, na.rm = TRUE)) /
      sd(value_std, na.rm = TRUE)
  ) |>
  mutate(FA = dense_rank(fa), TX = dense_rank(treat))


# Fatty acid summarise on complete data (same data as used for
# standardisation in fat_dat_compl, includes the omega-3 index)
fat_sum_compl <- fat_dat_compl |>
  summarise(
    .by = c(FA, fa),
    m = mean(log(value), na.rm = TRUE),
    s = sd(log(value), na.rm = TRUE)
  )


# Check n
fat_n_compl <- fat_dat |>
  distinct(time, participant) |>
  summarise(.by = c(time), n = n()) |>
  mutate(missing = 1 - n / max(n))


# Filter for complete case and prep data
fat_dat_compl2 <- fat_dat_compl |>
  select(participant, FA, TX, time, value_std) |>
  pivot_wider(names_from = time, values_from = value_std) |>

  pivot_longer(names_to = "TI", values_to = "y", cols = T2:T4) |>
  mutate(
    TI = dense_rank(TI),
    P = dense_rank(participant),
    # Add missing data and fill y
    obs = if_else(is.na(y), 0, 1),
    y = if_else(is.na(y), -999, y)
  )


fit_wrap(m4, DATA = d, fit = "fit.m4p", iter = 2000)
draws <- fit.m4p$draws(format = "df")


## Labels

# Labels
fa_labels <- data.frame(
  fa = c(
    "C16",
    "C18",
    "C181c9",
    "C182n6",
    "C183n3",
    "C183n6",
    "C203n6",
    "C204n6",
    "C205n3",
    "C225n3",
    "C226n3"
  ),
  fa_lab = c(
    "Palmitic Acid",
    "Stearic Acid",
    "Oleic Acid",
    "Linoleic Acid (n-6)",
    "Gamma-linolenic Acid (n-6)",
    "Alpha-linolenic Acid (n-3)",
    "Dihomo-gamma-linolenic Acid (n-6)",
    "Arachidonic Acid (n-6)",
    "Eicosapentaenoic Acid (EPA, n-3)",
    "Docosapentaenoic Acid (n-3)",
    "Docosahexaenoic Acid (DHA, n-3)"
  )
) |>
  mutate(
    fa_lab = factor(
      fa_lab,
      levels = c(
        "Docosahexaenoic Acid (DHA, n-3)",
        "Eicosapentaenoic Acid (EPA, n-3)",
        "Docosapentaenoic Acid (n-3)",
        "Alpha-linolenic Acid (n-3)",
        "Arachidonic Acid (n-6)",
        "Dihomo-gamma-linolenic Acid (n-6)",
        "Gamma-linolenic Acid (n-6)",
        "Linoleic Acid (n-6)",
        "Oleic Acid",
        "Stearic Acid",
        "Palmitic Acid"
      )
    )
  )


### Posterior predictions

pre_int_means <- fat_dat_compl2 |>
  filter(TI == 1) |>
  summarise(.by = FA, pre = mean(T1))


sumdat_all <- draws |>

  select(starts_with("alpha"), .draw) |>
  pivot_longer(cols = starts_with("alpha")) |>
  separate(name, into = c("param", "FA", "TI", "TX"), sep = "[,\\[]") |>
  mutate(TX = gsub("\\]", "", TX)) |>
  select(-param, .draw, FA:TX, alpha = value, ) |>

  inner_join(
    draws |>

      select(starts_with("sigma"), .draw) |>
      pivot_longer(cols = starts_with("sigma")) |>
      separate(name, into = c("param", "FA"), sep = "\\[") |>
      mutate(FA = gsub("\\]", "", FA)) |>
      separate(FA, into = c("FA", "TX")) |>
      select(-param, .draw, FA, TX, sigma = value)
  ) |>
  inner_join(
    draws |>

      select(starts_with("beta"), .draw) |>
      pivot_longer(cols = starts_with("beta")) |>
      separate(name, into = c("param", "FA"), sep = "\\[") |>
      mutate(FA = gsub("\\]", "", FA)) |>
      select(-param, .draw, FA, beta = value)
  ) |>
  mutate(FA = as.integer(FA), TX = as.integer(TX)) |>

  inner_join(pre_int_means) |>
  rowwise() |>

  mutate(
    mu = alpha + beta * pre,
    y_pred = rnorm(1, alpha + beta * pre, sigma)
  ) |>
  ungroup() |>
  summarise(
    .by = c(FA, TI, TX),
    Mu = mean(mu),
    lwr99 = quantile(mu, 0.005),
    upr99 = quantile(mu, 0.995),
    lwr95 = quantile(mu, 0.025),
    upr95 = quantile(mu, 0.975),
    lwr50 = quantile(mu, 0.25),
    upr50 = quantile(mu, 0.75),

    lwr_pred95 = quantile(y_pred, 0.025),
    upr_pred95 = quantile(y_pred, 0.975)
  ) |>

  # Back-transform all estimates
  pivot_longer(cols = Mu:upr_pred95, names_to = "est", values_to = "value") |>
  inner_join(fat_sum_compl) |>

  mutate(value_org = exp(value * s + m)) |>

  select(FA:est, fa, value_org) |>
  pivot_wider(names_from = est, values_from = value_org) |>

  mutate(
    time = if_else(TI == 1, "T2", "T4"),
    treat = if_else(TX == 1, "ctrl", "omega"),
    treat = factor(
      treat,
      levels = c("ctrl", "omega"),
      labels = c("Placebo", "N-3")
    )
  ) |>

  add_row(
    fat_dat_compl |>
      filter(time == "T1") |>
      summarise(.by = c(time, fa), Mu = mean(value, na.rm = TRUE)) |>
      mutate(TX = NA, TI = NA, FA = NA) |>
      expand_grid(treat = c("ctrl", "omega")) |>
      mutate(
        treat = factor(
          treat,
          levels = c("ctrl", "omega"),
          labels = c("Placebo", "N-3")
        )
      )
  )

# Fatty acids (labels in fa_labels) and the omega-3 index
sumdat <- sumdat_all |>
  inner_join(fa_labels, by = "fa")

o3i_sumdat <- sumdat_all |>
  filter(fa == "omega3index")


## Plotting (layout follows panel A of the subscores plot in
## 07-metsyn-inference.qmd: observed means + model-based adjusted means, and
## posterior densities of treatment contrasts)

time_labels <- c(T1 = "Baseline", T2 = "Supplementation", T4 = "Combined")

# Shared strip settings so that strips in the two panels have the same height
# and form one continuous gray band across the figure
strip_size <- 6.4
strip_bg <- element_rect(fill = "gray95", color = NA)

# Observed group means (all available observations). Geometric means, as the
# model is fitted on the log scale, so that they are comparable with the
# back-transformed model estimates.
obs_means_all <- fat_dat_compl |>
  filter(!is.na(value)) |>
  summarise(.by = c(fa, time, treat), Mu = exp(mean(log(value))), n = n()) |>
  mutate(
    treat = factor(
      treat,
      levels = c("ctrl", "omega"),
      labels = c("Placebo", "N-3")
    )
  )

obs_means <- obs_means_all |>
  inner_join(fa_labels, by = "fa")

o3i_obs_means <- obs_means_all |>
  filter(fa == "omega3index")

# Model-based adjusted means at T2 and T4
pred_plot_dat <- sumdat |>
  filter(!is.na(TI))

pred_plot <- function(dat, obs) {
  # Position
  pd <- position_dodge(width = 0.1)

  obs <- obs |>
    mutate(
      time = factor(time_labels[time], levels = time_labels)
    )

  dat |>
    mutate(
      time = factor(time_labels[time], levels = time_labels)
    ) |>

    ggplot(aes(time, Mu, color = treat, group = treat)) +

    # Observed means at all time-points (descriptive)
    geom_point(
      data = obs,
      aes(time, Mu, color = treat, group = treat),
      shape = 1,
      size = 1.5,
      position = position_dodge(width = 0.5)
    ) +
    geom_line(
      data = obs,
      aes(time, Mu, color = treat, group = treat),
      linetype = "dotted",
      position = position_dodge(width = 0.5)
    ) +

    # Model-based adjusted means at T2 and T4
    geom_errorbar(
      aes(ymin = lwr95, ymax = upr95),
      width = 0,
      linewidth = 0.4,
      position = pd
    ) +
    geom_errorbar(
      aes(ymin = lwr50, ymax = upr50),
      width = 0,
      linewidth = 1.1,
      position = pd
    ) +
    geom_point(
      aes(fill = treat),
      size = 1.5,
      color = "black",
      position = pd,
      shape = 21
    ) +

    facet_wrap(~fa_lab, scales = "free_y", ncol = 1) +
    scale_color_manual(values = colscale[c(4, 2)]) +
    scale_fill_manual(values = colscale[c(4, 2)]) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +

    alpha_theme() +
    theme(legend.position = "top") +

    theme(
      strip.text = element_text(size = strip_size),
      strip.background = strip_bg,
      # Let long fatty acid names extend beyond the panel
      strip.clip = "off",
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      # No right margin, strips continue into the contrast panel
      plot.margin = margin(5.5, 0, 5.5, 5.5),
      # Transparent so that it does not hide axis labels of the contrast panel
      plot.background = element_blank()
    ) +

    labs(
      x = "Time-point",
      y = "Fatty acid concentration (&mu;g&times;mL<sup>-1</sup>)",
      color = "Treatment",
      fill = "Treatment"
    )
}

fa_scores <- pred_plot(pred_plot_dat, obs_means)


### Posterior contrasts

contrast_dat_all <- draws |>

  select(starts_with("alpha"), .draw) |>
  pivot_longer(cols = starts_with("alpha")) |>
  separate(name, into = c("param", "FA", "TI", "TX"), sep = "[,\\[]") |>
  mutate(TX = gsub("\\]", "", TX)) |>
  select(.draw, FA:TX, alpha = value) |>
  mutate(FA = as.integer(FA), TI = as.integer(TI)) |>

  # Compute contrasts
  mutate(treat = if_else(TX == 1, "ctrl", "omega")) |>
  select(.draw:TI, treat, alpha) |>

  pivot_wider(names_from = treat, values_from = alpha) |>
  inner_join(fat_sum_compl, by = "FA") |>
  mutate(
    contrast = omega - ctrl,
    # Fold-difference N-3 / Placebo on the natural scale
    contrast_out = exp(s * contrast),
    time = factor(
      if_else(TI == 1, "Supplementation", "Combined"),
      levels = c("Supplementation", "Combined")
    )
  )

contrast_dat <- contrast_dat_all |>
  inner_join(fa_labels, by = "fa")

# Omega-3 index contrasts (fold-difference N-3 / Placebo)
o3i_contrast_dat <- contrast_dat_all |>
  filter(fa == "omega3index")


contrast_plot <- function(dat) {
  dat |>
    ggplot(aes(contrast_out, color = time)) +
    geom_line(stat = "density") +

    facet_wrap(~fa_lab, scales = "free", ncol = 1) +
    alpha_theme(yaxis = FALSE, text_base_size = 7) +
    theme(
      axis.title.y = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      axis.title.x = element_markdown(size = 7, margin = margin(t = -6)),
      # Strip background without text (transparent text keeps the height)
      strip.text = element_text(size = strip_size, color = "transparent"),
      strip.background = strip_bg,
      # No left margin, strips continue from the prediction panel
      plot.margin = margin(5.5, 5.5, 5.5, 0)
    ) +
    labs(x = "Fold-difference (N-3 / Placebo)", color = "Time") +

    geom_vline(xintercept = 1, lty = 2, color = "gray50") +
    scale_color_manual(values = c(colscale[3], colscale[6]))
}

fa_contrast <- contrast_plot(contrast_dat)


# Omega 3 index plot #######################################

o3i_sumdat <- sumdat_all |>
  filter(fa == "omega3index") |>
  mutate(
    time = factor(
      time,
      levels = c("T1", "T2", "T4"),
      labels = c("Baseline", "Supplementation", "Combined")
    )
  )

o3i_obs_means <- obs_means_all |>
  filter(fa == "omega3index") |>
  mutate(
    time = factor(
      time,
      levels = c("T1", "T2", "T4"),
      labels = c("Baseline", "Supplementation", "Combined")
    )
  )

pd <- position_dodge(width = 0.17)
pd2 <- position_dodge(width = 0.02)

omega_index_fig <- o3i_sumdat |>
  mutate(Mu = if_else(is.na(FA), NA, Mu)) |>

  ggplot(aes(time, Mu, color = treat, group = treat)) +

  # Observed means at all time-points (descriptive)
  geom_point(
    data = o3i_obs_means,
    aes(time, Mu, color = treat, group = treat),
    shape = 1,
    size = 1.5,
    position = pd2
  ) +
  geom_line(
    data = o3i_obs_means,
    aes(time, Mu, color = treat, group = treat),
    linetype = "dotted",
    position = pd2
  ) +

  # Model-based adjusted means at T2 and T4
  geom_errorbar(
    aes(ymin = lwr95, ymax = upr95),
    width = 0,
    linewidth = 0.4,
    position = pd
  ) +
  geom_errorbar(
    aes(ymin = lwr50, ymax = upr50),
    width = 0,
    linewidth = 1,
    position = pd
  ) +
  geom_point(
    aes(fill = treat),
    size = 1.5,
    color = "black",
    position = pd,
    shape = 21
  ) +

  scale_color_manual(values = colscale[c(4, 2)]) +
  scale_fill_manual(values = colscale[c(4, 2)]) +
  scale_x_discrete() +
  scale_y_continuous(
    limits = c(5.5, 9),
    breaks = scales::breaks_pretty(n = 3),
    expand = c(0, 0)
  ) +

  alpha_theme() +
  theme(legend.position = "top") +

  theme(
    strip.text = element_text(size = strip_size),
    strip.background = strip_bg,
    # Let long fatty acid names extend beyond the panel
    strip.clip = "off",
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.85),
    legend.key.height = unit(3, "mm"),
    # Transparent so that it does not hide axis labels of the contrast panel
    plot.background = element_blank()
  ) +

  labs(
    x = "Time-point",
    y = "Omega 3 Index (%)",
    color = "Treatment",
    fill = "Treatment"
  )


# Omega 3 index contrast #####################################

o3i_contrast <- o3i_contrast_dat |>

  ggplot(aes(contrast_out, color = time)) +
  geom_line(stat = "density") +

  alpha_theme(yaxis = FALSE, text_base_size = 7) +
  theme(
    axis.title.y = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.2, 0.80),
    legend.key.height = unit(3, "mm"),
    legend.title = element_blank(),
    axis.title.x = element_markdown(size = 7),
    # Strip background without text (transparent text keeps the height)
    strip.text = element_text(size = strip_size, color = "transparent"),
    strip.background = strip_bg,
  ) +
  labs(x = "Fold-difference (N-3 / Placebo)", color = "Time") +

  geom_vline(xintercept = 1, lty = 2, color = "gray50") +
  scale_color_manual(values = c(colscale[3], colscale[6]))


# Combining figures for output. ##############################

fig2 <- cowplot::plot_grid(
  omega_index_fig,
  o3i_contrast,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 0.6),
  ncol = 1
)

fig3 <- cowplot::plot_grid(
  fa_scores,
  fa_contrast,
  align = "h",
  rel_widths = c(1, 0.9)
)

# Save figures
saveRDS(fig2, file = "figures/figure-02.RDS")
saveRDS(fig3, file = "figures/figure-03.RDS")

ggsave(
  fig2,
  filename = "figures/figure-2.pdf",
  width = 90,
  height = 75,
  unit = "mm"
)

ggsave(
  fig3,
  filename = "figures/figure-3.pdf",
  width = 90,
  height = 220,
  unit = "mm"
)


ggsave(
  fig2,
  filename = "figures/figure-2.svg",
  width = 90,
  height = 75,
  unit = "mm"
)

ggsave(
  fig3,
  filename = "figures/figure-3.svg",
  width = 90,
  height = 220,
  unit = "mm"
)
