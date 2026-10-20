# Effect of N-3 supplementation on Glucose Biology
#
# This figure combines OGTT, fasting glucose and HOMA2-IR.
#
# OGTT panels are extracted from 08-ogtt.qmd (model mo3, chunk fig-mo3):
#   Top: predicted glucose curves per time x treatment
#   Middle: within-group changes from baseline at measured time points
#   Bottom: between-group contrasts of changes from baseline

#

# Load packages and functions #####################################
source("./R/fig-opt-fun.R")
source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")
source("./R/sim-fun.R")

# The spline package is needed for basis spline construction
library(splines)

# Settings
refit.models <- FALSE

# cowplot::plot_grid() lays out plots on a null device. The default pdf(NULL)
# gives wrong word spacing in ggtext::element_markdown() labels (e.g.
# "Fastingglucose"), the agg device (ragg) measures text correctly.
cowplot::set_null_device("agg")


# OGTT data preparation ###########################################

# Get the full index of participants
part_full <- criteria_dat |>
  distinct(participant) |>
  pull(participant)

# Glucose data
glucose_dat <- AOData::glucose |>
  select(participant:glucose) |>
  # Treatment
  inner_join(
    AOData::idkeys |>
      select(participant, treat)
  )

# Extract the full data sets.
full_dat <- glucose_dat |>
  inner_join(AOData::idkeys) |>
  filter(participant %in% part_full) |>

  dplyr::select(participant, time, min, glucose, treat) |>
  # Wider to find NA -> longer -> complete data set
  pivot_wider(values_from = glucose, names_from = time) |>
  pivot_longer(values_to = "glucose", names_to = "time", cols = T1:T4) |>
  arrange(participant, time, min) |>
  # Indexing
  mutate(
    y_obs = (glucose - mean(glucose, na.rm = TRUE)) /
      sd(glucose, na.rm = TRUE),
    y_obs = if_else(is.na(glucose), -999, y_obs),
    y_miss = if_else(is.na(glucose), 0, 1),
    tx = paste0(time, treat),
    tx = factor(
      tx,
      levels = c("T1ctrl", "T1omega", "T2ctrl", "T2omega", "T4ctrl", "T4omega")
    ),

    G = dense_rank(tx),
    P = dense_rank(participant),
    Tt = dense_rank(time)
  )

# For post-fitting analysis, mean and SD of glucose
mglucose <- mean(full_dat$glucose, na.rm = TRUE)
sdglucose <- sd(full_dat$glucose, na.rm = TRUE)

B <- bs(full_dat$min, df = 5)

dat <- list(
  N = nrow(full_dat),
  y_obs = full_dat$y_obs,
  y_miss = full_dat$y_miss,
  K = ncol(B),
  B = B,
  N_G = length(unique(full_dat$G)),
  G = full_dat$G,
  N_P = length(unique(full_dat$P)),
  P = full_dat$P,
  N_T = length(unique(full_dat$Tt)),
  Tt = full_dat$Tt,
  prior_only = 0
)


# Model fit (loads the saved fit unless refit.models is TRUE) #####
mo3 <- cmdstan_model(
  stan_file = "models/mo3.stan",
  exe_file = "models/mo3.exe"
)

fit_wrap(mo3, fit = "fit.mo3", DATA = dat)


# Posterior predictions ###########################################

draws <- fit.mo3$draws(variables = c("a_g", "b_g"), format = "df")

# Basis evaluated on a dense grid (B was built as bs(full_dat$min, df = 5))
min_seq <- seq(from = 0, to = 120, by = 1)

B_pred_long <- predict(B, min_seq) |>
  data.frame() |>
  cbind(min_seq) |>
  pivot_longer(starts_with("X"), values_to = "B")

# Map cell index g back to time / treatment labels
cell_key <- full_dat |>
  distinct(time, treat, G, Tt) |>
  rename(g = G)

ag_long <- draws |>
  dplyr::select(.draw, starts_with("a_g[")) |>
  pivot_longer(starts_with("a_g["), names_to = "par", values_to = "a_g") |>
  mutate(g = as.integer(gsub("a_g\\[(\\d+)\\]", "\\1", par))) |>
  dplyr::select(.draw, g, a_g)

bg_long <- draws |>
  dplyr::select(.draw, starts_with("b_g[")) |>
  pivot_longer(starts_with("b_g["), names_to = "par", values_to = "b_g") |>
  mutate(
    k = as.integer(gsub("b_g\\[(\\d+),(\\d+)\\]", "\\1", par)),
    g = as.integer(gsub("b_g\\[(\\d+),(\\d+)\\]", "\\2", par)),
    name = paste0("X", k)
  )

# Group-level predictions - These are implied group level glucose curves
# at the three time points
pred <- B_pred_long |>
  inner_join(bg_long, by = "name", relationship = "many-to-many") |>
  summarise(.by = c(.draw, g, min_seq), mu_spl = sum(B * b_g)) |>
  inner_join(ag_long, by = c(".draw", "g")) |>
  # back transform the standardization
  mutate(mu = (a_g + mu_spl) * sdglucose + mglucose) |>
  inner_join(cell_key) |>
  mutate(time_treat = paste0(time, "_", treat)) |>
  dplyr::select(.draw, min_seq, mu, time_treat)


# Within-group contrasts at measured time points ##################

within_contrasts <- pred |>
  pivot_wider(values_from = mu, names_from = time_treat) |>
  # Calculate contrasts
  mutate(
    ctrl_t2 = T2_ctrl - T1_ctrl,
    ctrl_t4 = T4_ctrl - T1_ctrl,
    omega_t2 = T2_omega - T1_omega,
    omega_t4 = T4_omega - T1_omega,
  ) |>
  select(.draw, min_seq, ctrl_t2:omega_t4) |>
  pivot_longer(cols = ctrl_t2:omega_t4) |>
  filter(min_seq %in% c(15, 30, 45, 60, 60, 90, 120)) |>
  summarise(
    .by = c(name, min_seq),
    m = median(value),
    lwr50 = quantile(value, 0.25),
    upr50 = quantile(value, 0.75),
    lwr95 = quantile(value, 0.025),
    upr95 = quantile(value, 0.975)
  ) |>
  separate(name, into = c("treat", "time"), sep = "_") |>
  mutate(
    time = case_when(
      time == "t2" ~ "Supplementation vs. Baseline",
      time == "t4" ~ "Combined vs. Baseline"
    ),
    time = factor(
      time,
      levels = c("Supplementation vs. Baseline", "Combined vs. Baseline")
    ),
    treat = factor(
      treat,
      levels = c("ctrl", "omega"),
      labels = c("Placebo change from baseline", "N-3 change from baseline")
    )
  )


# Plots ###########################################################

pd <- position_dodge(width = 3)

p_wcontr <- within_contrasts |>
  ggplot(aes(min_seq, m, color = time, fill = time)) +

  geom_hline(yintercept = 0, lty = 2, color = "gray") +

  geom_errorbar(
    aes(ymin = lwr50, ymax = upr50),
    linewidth = 1.2,
    width = 0,
    position = pd
  ) +
  geom_errorbar(
    aes(ymin = lwr95, ymax = upr95),
    linewidth = 0.6,
    width = 0,
    position = pd
  ) +
  geom_point(shape = 21, color = "black", position = pd) +

  scale_color_manual(values = colscale[c(3, 6)]) +
  scale_fill_manual(values = colscale[c(3, 6)]) +
  scale_x_continuous(limits = c(0, 125), expand = c(0, 0)) +

  labs(
    x = "Time (min)",
    y = " ΔBlood  glucose<br>(mmol L<sup>−1</sup>)"
  ) +

  facet_wrap(~treat) +
  alpha_theme() +
  theme(
    axis.title.y = ggtext::element_markdown(),
    legend.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text = ggtext::element_markdown(size = 7)
  )


p_pred <- pred |>
  summarise(
    .by = c(time_treat, min_seq),
    m = median(mu),
    lwr50 = quantile(mu, 0.25),
    upr50 = quantile(mu, 0.75),
    lwr95 = quantile(mu, 0.025),
    upr95 = quantile(mu, 0.975)
  ) |>
  separate(time_treat, into = c("time", "treat"), sep = "_") |>
  mutate(
    time = case_when(
      time == "T1" ~ "Baseline",
      time == "T2" ~ "Supplementation",
      time == "T4" ~ "Combined"
    ),
    time = factor(time, levels = c("Baseline", "Supplementation", "Combined")),
    treat = factor(
      treat,
      levels = c("ctrl", "omega"),
      labels = c("Placebo predictions", "N-3 predictions")
    )
  ) |>

  ggplot(aes(min_seq, m, lty = time, color = treat)) +
  geom_line() +
  facet_wrap(~treat, ncol = 2) +
  alpha_theme() +
  labs(x = "Time (min)", y = "Blood glucose<br>(mmol L<sup>-1</sup>)") +

  scale_x_continuous(limits = c(0, 125), expand = c(0, 0)) +

  scale_color_manual(values = colscale[c(4, 2)], guide = "none") +
  theme(
    axis.title.y = ggtext::element_markdown(),
    legend.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )


# Contrasts - These are differences in delta from baseline
delta_contrasts <- pred |>
  pivot_wider(values_from = mu, names_from = time_treat) |>
  # Calculate contrasts
  mutate(
    ctrl_t2 = T2_ctrl - T1_ctrl,
    ctrl_t4 = T4_ctrl - T1_ctrl,
    omega_t2 = T2_omega - T1_omega,
    omega_t4 = T4_omega - T1_omega,
    omega_ctrl_t2 = omega_t2 - ctrl_t2,
    omega_ctrl_t4 = omega_t4 - ctrl_t4,
  ) |>
  pivot_longer(
    names_to = "contrast",
    values_to = "value",
    cols = T1_ctrl:omega_ctrl_t4
  ) |>
  summarise(
    .by = c(min_seq, contrast),
    m = mean(value),
    lwr50 = quantile(value, 0.25),
    upr50 = quantile(value, 0.75),
    lwr95 = quantile(value, 0.025),
    upr95 = quantile(value, 0.975)
  ) |>
  filter(contrast %in% c("omega_ctrl_t2", "omega_ctrl_t4")) |>
  mutate(
    contrast = case_when(
      contrast == "omega_ctrl_t2" ~ "Placebo vs. N-3: Supplementation period",
      contrast == "omega_ctrl_t4" ~ "Placebo vs. N-3: Combined period"
    ),
    contrast = factor(
      contrast,
      levels = c(
        "Placebo vs. N-3: Supplementation period",
        "Placebo vs. N-3: Combined period"
      )
    )
  ) |>

  ggplot(aes(min_seq, m, fill = contrast)) +
  scale_x_continuous(
    limits = c(0, 125),
    expand = c(0, 0),
    breaks = c(25, 50, 75, 100)
  ) +
  geom_ribbon(aes(ymin = lwr50, ymax = upr50), alpha = 0.3) +
  geom_ribbon(aes(ymin = lwr95, ymax = upr95), alpha = 0.3) +
  geom_line() +
  facet_wrap(~contrast) +
  alpha_theme() +
  scale_fill_manual(values = colscale[c(3, 6)]) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.5) +
  theme(legend.position = "none") +
  labs(
    x = "Time (min)",
    y = "&Delta;N-3 − &Delta;Placebo<br>(mmol L<sup>-1</sup>)"
  )


# Combine the full figure #########################################

ogtt <- cowplot::plot_grid(
  p_pred,
  p_wcontr,
  delta_contrasts,
  ncol = 1,
  rel_heights = c(1, 1, 1.4),
  align = "v"
)


# Fasting glucose and HOMA2-IR ####################################
#
# Extracted from 07-metsyn-inference.qmd (chunks load-data and
# infer-fig-data-prep). Left: observed and model-based adjusted means (s13),
# right: posterior distributions of treatment contrasts (N-3 - Placebo).

# Labels, data and model. ms_ prefix avoids clashing with the OGTT objects.
ms_dat <- data_load_s7s10()

ms_labels <- ms_dat$scale |>
  select(type, analyte) |>
  mutate(
    label = case_when(
      analyte == "HDL" ~ "High density lipoprotein",
      analyte == "LDL" ~ "Low density lipoprotein",
      analyte == "glucose" ~ "Fasting glucose",
      analyte == "TG" ~ "Triglycerides",
      analyte == "waist" ~ "Waist circumference",
      analyte == "sys" ~ "Systolic blood pressure",
      analyte == "dias" ~ "Diastolic blood pressure",
      analyte == "fat.rel" ~ "Relative fat mass",
      analyte == "HOMA2IR" ~ "HOMA2-IR",
      analyte == "bmi" ~ "Body mass index",
      analyte == "Adiponectin" ~ "Adiponectin",
      .default = analyte
    )
  ) |>
  arrange(type)

# S13 is the main model, multivariate t distribution.
s13 <- cmdstan_model(
  stan_file = "models/s13.stan",
  exe_file = "models/s13.exe"
)

fit_wrap(s13, DATA = ms_dat$d, fit = "fit.s13.post")


# Units on the natural scale (used for observed and predicted means)
unit_nat <- function(analyte) {
  case_when(
    analyte %in% c("HDL", "LDL", "TG", "glucose") ~ "mmol&times;L<sup>-1</sup>",
    analyte %in% c("sys", "dias") ~ "mmHg",
    analyte == "waist" ~ "cm",
    analyte == "bmi" ~ "kg&times;m<sup>-2</sup>",
    analyte == "fat.rel" ~ "%",
    analyte == "Adiponectin" ~ "mg&times;L<sup>-1</sup>",
    analyte == "HOMA2IR" ~ "index (log-scale)",
    .default = NA_character_
  )
}

analyte_levels <- c(
  "waist",
  "dias",
  "sys",
  "glucose",
  "HDL",
  "TG",
  "bmi",
  "fat.rel",
  "LDL",
  "Adiponectin",
  "HOMA2IR"
)

strip_lab <- function(label, unit) paste0(label, " (", unit, ")")


# Observed data on the natural scale (all available observations)
obs_dat <- criteria_dat |>
  select(participant, time, treat, sex, analyte, value) |>
  filter(!is.na(value)) |>
  inner_join(ms_labels, by = "analyte") |>
  inner_join(ms_dat$scale |> select(type, transform), by = "type") |>
  mutate(
    treat = factor(
      treat,
      levels = c("ctrl", "omega"),
      labels = c("Placebo", "N-3")
    )
  )

# Observed group means. Geometric means for the log-transformed outcomes
# (TG, HOMA2IR) so that they are comparable with the back-transformed model
# estimates, arithmetic means otherwise.
obs_means <- obs_dat |>
  summarise(
    .by = c(type, analyte, label, time, treat),
    Mu = if_else(
      first(transform) == "log",
      exp(mean(log(value))),
      mean(value)
    ),
    n = n()
  ) |>
  mutate(unit = unit_nat(analyte))

# Pre-intervention covariate: overall mean of the standardised baseline value
# per type (the same value for both groups gives adjusted group means). Sex is
# coded -0.5/0.5, so sex = 0 gives the average over sexes.
pre_means <- ms_dat$dat |>
  filter(pre_miss == 1) |>
  distinct(P, type, pre_obs) |>
  summarise(.by = type, pre = mean(pre_obs))

# Model-based adjusted means, back-transformed to the natural scale
# alpha[time, type, treat] are time-specific intercepts, beta[type, time] the
# pre-value slopes; mu = alpha + beta * pre on the standardised scale.
predict_means <- function(fit, model) {
  alpha <- fit$draws("alpha", format = "df") |>
    select(.draw, starts_with("alpha")) |>
    pivot_longer(cols = -.draw, values_to = "alpha") |>
    mutate(name = gsub("alpha\\[|\\]", "", name)) |>
    separate(
      name,
      into = c("time", "type", "treat"),
      sep = ",",
      convert = TRUE
    )

  beta <- fit$draws("beta", format = "df") |>
    select(.draw, starts_with("beta")) |>
    pivot_longer(cols = -.draw, values_to = "beta") |>
    mutate(name = gsub("beta\\[|\\]", "", name)) |>
    separate(name, into = c("type", "time"), sep = ",", convert = TRUE)

  alpha |>
    inner_join(beta, by = c(".draw", "time", "type")) |>
    inner_join(pre_means, by = "type") |>
    mutate(mu = alpha + beta * pre) |>
    summarise(
      .by = c(time, type, treat),
      Mu = median(mu),
      lwr95 = quantile(mu, 0.025),
      upr95 = quantile(mu, 0.975),
      lwr50 = quantile(mu, 0.25),
      upr50 = quantile(mu, 0.75)
    ) |>
    # Back-transform (monotone, so quantiles map directly)
    inner_join(ms_dat$scale, by = "type") |>
    mutate(across(
      c(Mu, lwr95, upr95, lwr50, upr50),
      ~ if_else(transform == "log", exp(mean + sd * .x), mean + sd * .x)
    )) |>
    mutate(
      time = if_else(time == 1, "T2", "T4"),
      treat = factor(treat, levels = c(1, 2), labels = c("Placebo", "N-3")),
      model = model
    ) |>
    select(model, time, type, analyte, treat, Mu, lwr95, upr95, lwr50, upr50)
}

pred_plot_dat <- predict_means(fit.s13.post, "s13") |>
  inner_join(ms_labels, by = c("type", "analyte")) |>
  mutate(
    unit = unit_nat(analyte),
    analyte = factor(analyte, levels = analyte_levels)
  ) |>
  mutate(
    strip.lab = strip_lab(label, unit),
    # Facet strips follow the order of the analyte levels
    strip.lab = fct_reorder(strip.lab, as.integer(analyte), .fun = min)
  )


# Contrasts on the standardised scale (alpha is the intercept at pre = 0 and
# sex = 0, and there is no treatment x pre interaction, so omega - ctrl is
# the adjusted mean difference at any pre value). Back-transformation uses
# the per-type sd:
# contrast_nat = sd * contrast  -> difference in natural units, or on the
#                                   log scale for TG and HOMA2IR
# ratio = exp(sd * contrast) -> N-3 / placebo ratio (log outcomes only)
contrast_dat <- fit.s13.post$draws("alpha", format = "df") |>
  pivot_longer(cols = starts_with("alpha")) |>
  mutate(name = gsub("alpha\\[|\\]", "", name)) |>
  separate(name, into = c("time", "type", "treat")) |>
  mutate(treat = if_else(treat == 1, "ctrl", "omega")) |>
  pivot_wider(names_from = treat, values_from = value) |>
  mutate(contrast = omega - ctrl, type = as.integer(type)) |>
  inner_join(
    ms_dat$scale |> select(type, analyte, transform, sd),
    by = "type"
  ) |>
  mutate(
    contrast_nat = sd * contrast,
    ratio = if_else(transform == "log", exp(contrast_nat), NA_real_)
  ) |>
  inner_join(ms_labels, by = c("type", "analyte")) |>
  mutate(
    unit = if_else(
      analyte %in% c("TG", "HOMA2IR"),
      "Fold-difference",
      unit_nat(analyte)
    ),
    analyte = factor(analyte, levels = analyte_levels)
  ) |>
  mutate(
    strip.lab = strip_lab(label, unit),
    strip.lab = fct_reorder(strip.lab, as.integer(analyte), .fun = min),
    contrast_out = if_else(transform == "log", ratio, contrast_nat),
    # No difference line (0 for natural scale and 1 for fold-difference)
    no_diff = if_else(analyte %in% c("TG", "HOMA2IR"), 1, 0)
  )


# Plot functions ##################################################

pred_plot <- function(dat, obs, analytes) {
  # Position
  pd <- position_dodge(width = 0.1)

  obs_means <- obs |>
    filter(analyte %in% analytes) |>
    # Same factor levels as dat so both layers share the panel order
    mutate(analyte = factor(analyte, levels = levels(dat$analyte))) |>
    mutate(
      time = factor(time, labels = c("Baseline", "Supplementation", "Combined"))
    )

  # Strip labels keyed by analyte
  strip_labels <- dat |>
    distinct(analyte, strip.lab) |>
    with(setNames(as.character(strip.lab), analyte))

  dat |>
    mutate(time = factor(time, labels = c("Supplementation", "Combined"))) |>
    filter(analyte %in% analytes) |>

    ggplot(aes(time, Mu, color = treat, group = interaction(treat, model))) +

    # Observed means at all time-points (descriptive)
    geom_point(
      data = obs_means,
      aes(time, Mu, color = treat, group = treat),
      shape = 1,
      size = 1.5,
      position = position_dodge(width = 0.5)
    ) +
    geom_line(
      data = obs_means,
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

    facet_wrap(
      ~analyte,
      scales = "free_y",
      ncol = 1,
      labeller = as_labeller(strip_labels)
    ) +
    scale_color_manual(values = colscale[c(4, 2)]) +
    scale_fill_manual(values = colscale[c(4, 2)]) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
    alpha_theme() +
    theme(legend.position = "top") +
    theme(
      strip.text = ggtext::element_markdown(),
      axis.title.x = element_blank(),
      legend.title = element_blank()
    ) +
    labs(
      x = "Time-point",
      y = "",
      color = "Treatment",
      fill = "Treatment"
    )
}


contrast_plot <- function(dat, analytes) {
  dat |>
    filter(analyte %in% analytes) |>
    mutate(
      time = factor(
        if_else(time == 1, "Supplementation", "Combined"),
        levels = c("Supplementation", "Combined")
      )
    ) |>
    ggplot(aes(contrast_out, color = time)) +
    geom_line(stat = "density") +
    facet_wrap(~strip.lab, scales = "free", ncol = 1) +
    alpha_theme(yaxis = FALSE, text_base_size = 7) +
    theme(
      axis.title.y = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      axis.title.x = ggtext::element_markdown(
        size = 7,
        margin = margin(t = -6)
      ),
      strip.text = element_blank()
    ) +
    labs(x = "Contrast (N-3 - Placebo)", color = "Time") +
    geom_vline(aes(xintercept = no_diff)) +
    scale_color_manual(values = c(colscale[3], colscale[6]))
}


# Glucose and HOMA2-IR panels #####################################

glu_analytes <- c("glucose", "HOMA2IR")

glu_pred <- pred_plot(pred_plot_dat, obs_means, analytes = glu_analytes)

glu_contrast <- contrast_plot(contrast_dat, analytes = glu_analytes)

glu_homa <- cowplot::plot_grid(glu_pred, glu_contrast, align = "h")

figure_04 <- cowplot::plot_grid(
  ogtt,
  glu_homa,
  ncol = 1,
  rel_heights = c(1, 0.8)
)


ggsave(
  filename = "figures/figure-4.pdf",
  figure_04,
  device = cairo_pdf,
  width = 18,
  height = 15,
  units = "cm"
)

ggsave(
  filename = "figures/figure-4.svg",
  figure_04,
  width = 18,
  height = 15,
  units = "cm"
)
