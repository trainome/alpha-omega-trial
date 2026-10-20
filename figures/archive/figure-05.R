# Effect of N-3 supplementation on metabolic syndrome
#
# Extracted from 07-metsyn-inference.qmd (chunks load-data,
# infer-fig-data-prep, joint-burden-score, diagnostic-criteria-analysis and
# fig-metsyn-inference):
#   A: observed and model-based adjusted means (left) and treatment contrasts
#      (right) for the metabolic syndrome sub-scores
#   B: average burden score contrast (N-3 - Placebo)
#   C: percentage of participants meeting the diagnostic criteria
#   D: treatment contrast for metabolic syndrome diagnosis (log OR)

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

# Text formatting of posterior probabilities
fmt_prob <- function(p) {
  case_when(
    p > 0.99 ~ ">99%",
    p < 0.01 ~ "<1%",
    .default = paste0(round(100 * p), "%")
  )
}


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


met_syn_analytes <- c("HDL", "TG", "waist", "sys", "dias")

met_sym_scores <- pred_plot(
  pred_plot_dat,
  obs_means,
  analytes = met_syn_analytes
)

met_syn_contrast <- contrast_plot(contrast_dat, analytes = met_syn_analytes)


# Joint burden score ##############################################

# Extract all alpha effects
a <- fit.s13.post$draws("alpha", format = "df")

# A function for getting the t:th time point and k:th effect from the matrix
# and contrasted between treatment groups (2 vs 1).
eff <- function(t, k) {
  a[[sprintf("alpha[%d,%d,2]", t, k)]] - a[[sprintf("alpha[%d,%d,1]", t, k)]]
}

# Combined posterior effect. The k indices are from the model formulation
# (see s13.stan and ms_dat$scale):
# Type     # Type index
# HDL         2
# TG          5
# dias        7
# glucose     9
# sys         10
# waist       11
#
# The score is
# (0.5 * (dias + sys) + (- HDL) + TG + glucose + waist) / 5
burden_eff <- sapply(1:2, function(t) {
  (-eff(t, 2) +
    eff(t, 5) +
    0.5 * (eff(t, 7) + eff(t, 10)) +
    eff(t, 9) +
    eff(t, 11)) /
    5
})

colnames(burden_eff) <- c("T2", "T4")

# Posterior summaries used for the P(benefit) labels
burden_contrast <- burden_eff |>
  data.frame() |>
  pivot_longer(everything()) |>
  summarise(
    .by = c(name),
    # Median posterior on the standardized scale (sd units)
    med = median(value),
    # Density height at the median (same defaults as stat_density)
    dens_med = with(density(value), approx(x, y, xout = med)$y),
    # Proportion of posterior that shows benefit in treatment
    benefit = mean(value < 0)
  ) |>
  mutate(
    prob = paste("P(benefit)", fmt_prob(benefit)),
    name = factor(
      name,
      levels = c("T2", "T4"),
      labels = c("Supplementation", "Combined")
    )
  )

burden_fig <- burden_eff |>
  data.frame() |>
  pivot_longer(everything()) |>
  mutate(
    name = factor(
      name,
      levels = c("T2", "T4"),
      labels = c("Supplementation", "Combined")
    )
  ) |>
  ggplot(aes(value, color = name)) +
  # Shade the part of each density indicating benefit (N-3 < Placebo)
  geom_area(
    aes(y = after_stat(ifelse(x < 0, density, 0)), fill = name),
    stat = "density",
    color = NA,
    alpha = 0.3,
    position = "identity"
  ) +
  geom_line(stat = "density", alpha = 0.8, show.legend = FALSE) +

  geom_vline(xintercept = 0, lty = 2, color = "gray50") +

  geom_text(
    data = burden_contrast,
    aes(x = med - 0.12, y = dens_med, label = prob),
    show.legend = FALSE,
    size = 6 / .pt
  ) +

  alpha_theme(yaxis = FALSE, text_base_size = 7) +
  theme(
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.direction = "horizontal",
    legend.key.size = unit(3, "mm")
  ) +
  labs(
    color = "Treatment-period",
    fill = "Treatment-period",
    x = "Average burden score<br>N-3 - Placebo (SD-units)"
  ) +
  scale_fill_manual(values = colscale[c(3, 6)]) +
  scale_color_manual(values = colscale[c(3, 6)])


# Diagnostic criteria #############################################

# Treatment allocations per participant and time-point
dobs <- ms_dat$dat |>
  distinct(P, treat) |>
  expand_grid(time = c(1, 2, 3))

draws_diag <- fit.s13.post$draws("diagnosed", format = "df")


# Percentage diagnosed per group, uncertainty from missing-value imputation
infer_diagnosed <- function(draws, dobs, ndraws = 1000) {
  pd <- position_dodge(width = 0.15)

  dat_out <- draws |>
    filter(.draw %in% sample(1:8000, ndraws)) |>
    pivot_longer(cols = starts_with("diagnosed")) |>
    mutate(name = gsub("diagnosed\\[|\\]", "", name)) |>
    separate(name, into = c("P", "time")) |>
    mutate(P = as.integer(P), time = as.integer(time)) |>

    inner_join(dobs) |>

    mutate(
      time = if_else(
        time == 1,
        "Baseline",
        if_else(time == 2, "Supplementation", "Combined")
      ),
      time = factor(time, levels = c("Baseline", "Supplementation", "Combined")),

      treat = if_else(treat == 1, "Placebo", "N-3"),
      treat = factor(treat, levels = c("Placebo", "N-3"))
    ) |>

    summarise(
      .by = c(treat, time, .draw),
      frac.diagnosed = 100 * (sum(value) / n())
    ) |>
    summarise(
      .by = c(treat, time),
      m = median(frac.diagnosed),
      lwr95 = quantile(frac.diagnosed, 0.025),
      lwr50 = quantile(frac.diagnosed, 0.25),
      upr50 = quantile(frac.diagnosed, 0.75),
      upr95 = quantile(frac.diagnosed, 0.975)
    )

  dat_out |>
    ggplot(aes(
      time,
      m,
      fill = as.factor(treat),
      color = as.factor(treat),
      group = as.factor(treat)
    )) +

    geom_errorbar(
      aes(ymin = lwr95, ymax = upr95),
      linewidth = 0.4,
      width = 0,
      position = pd
    ) +
    geom_errorbar(
      aes(ymin = lwr50, ymax = upr50),
      linewidth = 1.5,
      width = 0,
      position = pd
    ) +

    geom_line(position = pd) +

    geom_point(shape = 21, size = 0.9, color = "black", position = pd) +

    scale_fill_manual(values = colscale[c(4, 2)]) +
    scale_color_manual(values = colscale[c(4, 2)]) +
    labs(y = "Diagnosis (%)") +

    alpha_theme(text_base_size = 7) +

    theme(
      legend.position = "inside",
      legend.position.inside = c(0.5, 0.9),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.key.height = unit(2, "mm"),
      legend.background = element_blank()
    ) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2))
}

diagnosis_per_treat <- infer_diagnosed(draws_diag, dobs, ndraws = 100)


# Diagnosis differences between treatment groups.
# The diagnosis is a 0/1 variable built from the sub-scores and their
# thresholds, computed per MCMC iteration in the generated quantities block.
# A logistic model is fitted to each of a set of iterations (seqwrap) and
# estimates are pooled with Rubin's rules, so that the intervals include both
# imputation and sampling uncertainty.
infer_seqwrap <- function(draws, dobs, ndraws = 500) {
  # Create metadata
  metadata <- draws |>
    filter(.draw == 1) |>
    select(.draw, starts_with("diagnosed")) |>
    pivot_longer(cols = starts_with("diag")) |>
    mutate(
      name = gsub("\\[|\\]", "", name),
      name = gsub(",", "_", name),
      part_time = gsub("diagnosed", "", name)
    ) |>
    separate(part_time, into = c("P", "time")) |>
    mutate(P = as.integer(P)) |>
    inner_join(dobs |> distinct(P, treat)) |>
    select(seq_sample_id = name, P, time, treat) |>
    mutate(treat = paste0("treat", treat), time = paste0("time", time))

  # Data set
  dat <- draws |>
    slice_sample(n = ndraws) |>
    select(.draw, starts_with("diagnosed")) |>
    pivot_longer(cols = starts_with("diag")) |>
    mutate(name = gsub("\\[|\\]", "", name), name = gsub(",", "_", name)) |>
    pivot_wider(values_from = value, names_from = name)

  # Prepare and run seqwrap
  model <- seqwrap::seqwrap_compose(
    modelfun = glmmTMB::glmmTMB,
    arguments = alist(
      formula = y ~ time * treat + (1 | P),
      family = stats::binomial,
      # Weak priors on fixed effects parameters in each run.
      priors = data.frame(
        prior = rep("normal(0,3)", 5),
        class = rep("fixef", 5),
        coef = c(
          "timetime2",
          "timetime3",
          "treattreat2",
          "timetime2:treattreat2",
          "timetime3:treattreat2"
        )
      )
    ),
    data = dat,
    samplename = "seq_sample_id",
    metadata = metadata
  )

  results <- seqwrap::seqwrap(model, cores = 12)

  seqwrap::seqwrap_summarise(results)
}

# Retrieved from saved results if already there
if (!file.exists("data-derived/ch-07-diag-crit-rubin-s13.RDS")) {
  sw_results13 <- infer_seqwrap(draws_diag, dobs = dobs, ndraws = 500)

  saveRDS(sw_results13, "data-derived/ch-07-diag-crit-rubin-s13.RDS")
} else {
  sw_results13 <- readRDS("data-derived/ch-07-diag-crit-rubin-s13.RDS")
}


# Rubin's rules pooling of the time x treatment interaction terms
pool_rubin <- function(sums) {
  sums |>
    filter(term %in% c("timetime2:treattreat2", "timetime3:treattreat2")) |>
    summarise(
      .by = term,
      m = n(),
      qbar = mean(estimate), # Pooled estimate
      ubar = mean(std.error^2), # Within-impute variance
      b = var(estimate) # Between-impute variance
    ) |>
    mutate(
      # Total variance
      t = ubar + (1 + 1 / m) * b,
      # Fraction of missing information
      lambda = (1 + 1 / m) * b / t,
      df = (m - 1) / lambda^2,
      se = sqrt(t),
      # Confidence intervals (95%)
      lwr = qbar - qt(0.975, df) * se,
      upr = qbar + qt(0.975, df) * se
    )
}

est_13 <- pool_rubin(sw_results13$summaries)

diag_criteria_fig <- est_13 |>
  mutate(
    term = if_else(
      term == "timetime2:treattreat2",
      "N-3 - Control\nSupplementation",
      "N-3 - Control\nCombined treatment"
    ),
    term = factor(
      term,
      levels = c(
        "N-3 - Control\nCombined treatment",
        "N-3 - Control\nSupplementation"
      )
    )
  ) |>

  ggplot(aes(term, qbar, color = term)) +
  geom_errorbar(
    aes(ymin = lwr, ymax = upr),
    width = 0,
    show.legend = FALSE
  ) +
  geom_point(size = 1.5, show.legend = FALSE) +
  scale_color_manual(values = colscale[c(3, 6)]) +
  alpha_theme(text_base_size = 7) +
  labs(x = "", y = "Log(OR) (&pm; 95% CI)") +
  geom_hline(yintercept = 0, lty = 2, color = "gray50") +
  coord_flip() +
  theme(axis.title.x = ggtext::element_markdown())


# Combine the full figure #########################################

figure_05 <- cowplot::plot_grid(
  cowplot::plot_grid(met_sym_scores, met_syn_contrast, align = "h", 
                     ncol = 2, rel_widths = c(1, 0.8)),
  cowplot::plot_grid(
    burden_fig,
    cowplot::plot_grid(diagnosis_per_treat, diag_criteria_fig, ncol = 1),
    rel_widths = c(0.5, 0.5)
  ),
  rel_heights = c(1, 0.3),
  ncol = 1
) +
  annotate(
    "text",
    x = c(0.025, 0.025, 0.48, 0.55),
    y = c(0.95, 0.22, 0.22, .1),
    label = c("A", "B", "C", "D"),
    size = 12 / .pt
  )


ggsave(
  filename = "figures/figure-5.pdf",
  figure_05,
  device = cairo_pdf,
  width = 9,
  height = 18,
  units = "cm"
)

ggsave(
  filename = "figures/figure-5.svg",
  figure_05,
  width = 9,
  height = 18,
  units = "cm"
)
