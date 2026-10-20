# Effect of N-3 supplementation on related outcomes
#
# Extracted from 07-metsyn-inference.qmd (chunks load-data,
# infer-fig-data-prep and fig-metsyn-related): observed and model-based
# adjusted means (left) and treatment contrasts (right) for relative fat mass,
# adiponectin, body mass index and LDL. HOMA2-IR is shown in figure 4.

# Load packages and functions #####################################
source("./R/fig-opt-fun.R")
source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")
source("./R/sim-fun.R")

# Settings
refit.models <- FALSE

# cowplot::plot_grid() lays out plots on a null device. The default pdf(NULL)
# gives wrong word spacing in ggtext::element_markdown() labels (e.g.
# "Fastingglucose"), the agg device (ragg) measures text correctly.
cowplot::set_null_device("agg")


# Data and model ##################################################

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


# Sub-score data ##################################################

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


# Sub-score plot functions ########################################

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
    scale_y_continuous(breaks = scales::breaks_pretty(n = 4)) +
    alpha_theme() +
    theme(legend.position = "top") +
    theme(
      strip.text = ggtext::element_markdown(),
      axis.title.x = element_blank(),
      legend.title = element_blank()
    ) +
    labs(
      x = "Time-point",
      y = "Outcome (natural scale)",
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


related_analytes <- c("fat.rel", "Adiponectin", "bmi", "LDL")

rel_scores <- pred_plot(pred_plot_dat, obs_means, analytes = related_analytes)

rel_contrast <- contrast_plot(contrast_dat, analytes = related_analytes)


# Combine the full figure #########################################

figure_06 <- cowplot::plot_grid(rel_scores, rel_contrast, align = "h")


ggsave(
  filename = "figures/figure-6.pdf",
  figure_06,
  device = cairo_pdf,
  width = 9,
  height = 13,
  units = "cm"
)

ggsave(
  filename = "figures/figure-6.svg",
  figure_06,
  width = 9,
  height = 13,
  units = "cm"
)
