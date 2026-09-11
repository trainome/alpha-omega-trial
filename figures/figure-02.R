# Effect of N-3 supplementation on Metabolic Syndrom and related measures.
#
#
#
#
#
#
#
#
#
#

# Load packages and functions
source("./R/fig-opt-fun.R")
source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")
source("./R/sim-fun.R")

# Extra packages
library(ggh4x)


# Extract relevant model posterior
fit_wrap(s10, DATA = d, fit = "fit.s10.post", refit = FALSE)
draws <- fit.s10.post$draws("alpha", format = "df")


# Create labels and extract standardizers
labels <- d_full |>
  distinct(analyte) |>
  mutate(
    type = as.character(dense_rank(analyte)),
    label = case_when(
      analyte == "HDL" ~ "High density lipoprotein<br>(mmol L<sup>-1</sup>)",
      analyte == "glucose" ~ "Fasting glucose<br>(mmol L<sup>-1</sup>)",
      analyte == "TG" ~ "Triglycerides<br>(%)",
      analyte == "waist" ~ "Waist circumference<br>(cm)",
      analyte == "dias" ~ "Diastolic blood pressure<br>(mmHg)",
      analyte == "sys" ~ "Systolic blood pressure<br>(mmHg)",
      analyte == "fat.rel" ~ "Relative fat mass<br>(%-points)",
      analyte == "HOMA2IR" ~ "Insulin resistance index<br>(AU)",
      analyte == "bmi" ~ "Body mass index<br>(kg m<sup>-2</sup>)",
      analyte == "Adiponectin" ~ "Adiponectin<br>(&mu;g L<sup>-1</sup>)",
      .default = analyte
    )
  ) |>
  arrange(type) |>

  mutate(
    label = factor(
      label,
      levels = c(
        "High density lipoprotein<br>(mmol L<sup>-1</sup>)",
        "Fasting glucose<br>(mmol L<sup>-1</sup>)",
        "Triglycerides<br>(%)",
        "Waist circumference<br>(cm)",
        "Diastolic blood pressure<br>(mmHg)",
        "Systolic blood pressure<br>(mmHg)",
        "Relative fat mass<br>(%-points)",
        "Body mass index<br>(kg m<sup>-2</sup>)",
        "Insulin resistance index<br>(AU)",
        "Adiponectin<br>(&mu;g L<sup>-1</sup>)"
      )
    )
  )


# Double check thresholds and type index

#   Subscore = c("HDL Cholesterol",
#                "Fasting glucose",
#                "Triglycerides",
#                "Waist circumference",
#                "Blood pressure",),
#   Threshold = c("< 40 mg/dL for men, < 50 mg/dL for women",
#                 "100 mg/dL or higher (or taking medication)",
#                 "150 mg/dL or higher (or taking medication)",
#                 "> 102 cm for men, > 88 cm for women",
#                 "130/85 mm Hg or higher (or taking medication)")) |>
#

# This collects the thresholds values for validation
crit_thresholds <- d_full |>

  # TG is on the log-scale
  mutate(value = if_else(analyte == "TG", log(value), value)) |>

  summarise(
    .by = c(analyte),
    m = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE)
  ) |>
  mutate(
    threshold_male = c(1.03, 5.6, log(1.7), 94, 130, 85, NA, NA, NA, NA),
    threshold_female = c(1.29, 5.6, log(1.7), 80, 130, 85, NA, NA, NA, NA),
    std_threshold_male = (threshold_male - m) / sd,
    std_threshold_female = (threshold_female - m) / sd,
    type = as.character(dense_rank(analyte))
  ) |>
  select(analyte, type, m, sd, std_threshold_male, std_threshold_female) |>
  inner_join(labels) |>

  mutate(
    label = factor(
      label,
      levels = c(
        "High density lipoprotein<br>(mmol L<sup>-1</sup>)",
        "Fasting glucose<br>(mmol L<sup>-1</sup>)",
        "Triglycerides<br>(%)",
        "Waist circumference<br>(cm)",
        "Diastolic blood pressure<br>(mmHg)",
        "Systolic blood pressure<br>(mmHg)",
        "Relative fat mass<br>(%-points)",
        "Body mass index<br>(kg m<sup>-2</sup>)",
        "Insulin resistance index<br>(AU)",
        "Adiponectin<br>(&mu;g L<sup>-1</sup>)"
      )
    )
  ) |>
  arrange(labels)


# Plot 1: Effects of treatment on sub-scores and related scores
#
#
# This plot gives the sub-score contrast estimates as densities and
# probabilities N-3 > Placebo.
#
#
fig02_p1 <- function(draws) {
  # Save estimates
  estimates <- draws |>
    pivot_longer(cols = starts_with("alpha")) |>
    mutate(name = gsub("alpha\\[|\\]", "", name)) |>

    separate(name, into = c("time", "type", "treat")) |>
    mutate(treat = if_else(treat == 1, "ctrl", "omega"))

  # Restructure to get contrasts
  contrasts <- estimates |>
    pivot_wider(names_from = treat, values_from = value) |>
    mutate(contrast = omega - ctrl) |>
    # Multiply by sd to convert to org scale
    inner_join(crit_thresholds) |>
    mutate(
      contrast_org = if_else(
        analyte %in% c("HOMA2IR", "TG"),
        1 - exp(contrast * sd),
        contrast * sd
      )
    )

  # Summarise the posterior
  sums <- contrasts |>

    summarise(
      .by = c(time, type),
      med = median(contrast),
      q05 = quantile(contrast, 0.05),
      q25 = quantile(contrast, 0.25),
      q75 = quantile(contrast, 0.75),
      q95 = quantile(contrast, 0.95),
      p = sum(contrast > 0) / n()
    ) |>
    # Add probability label
    mutate(prob_label = paste0("Prob(N-3 - Placebo > 0)=", round(p, 2))) |>
    inner_join(crit_thresholds)

  # Compute symmetric limits per facet
  x_limits <- contrasts |>
    group_by(label) |>
    summarise(lim = max(abs(contrast_org)), .groups = "drop")

  # Build a named list of scale overrides
  x_scales <- x_limits |>
    mutate(scale = map(lim, ~ scale_x_continuous(limits = c(-.x, .x)))) |>
    pull(scale)

  contrasts |>
    ggplot(aes(contrast_org, fill = time)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~label, scale = "free", ncol = 1) +
    alpha_theme(yaxis = FALSE) +

    geom_text(
      data = filter(sums, time == 1),
      mapping = aes(x = -Inf, y = -Inf, label = prob_label),
      hjust = -0.1, # Adjust horizontal justification
      vjust = -3, # Adjust vertical justification
      size = 3,
      color = colscale[1]
    ) +

    geom_text(
      data = filter(sums, time == 2),
      mapping = aes(x = -Inf, y = -Inf, label = prob_label),
      hjust = -0.1, # Adjust horizontal justification
      vjust = -1.5, # Adjust vertical justification
      size = 3,
      color = colscale[2]
    ) +

    geom_vline(xintercept = 0) +
    scale_fill_manual(values = colscale) +
    labs(fill = "Time") +
    theme(strip.text = element_markdown(), axis.title.y = element_blank()) +
    facetted_pos_scales(x = x_scales)
}


p1 <- fig02_p1(draws = draws)


draws <- fit.s12.post$draws("alpha", format = "df")
draws_sigma <- fit.s12.post$draws("sigma", format = "df")

fig02_p2 <- function(draws, draws_sigma) {
  # Save estimates
  estimates <- draws |>
    pivot_longer(cols = starts_with("alpha")) |>
    mutate(name = gsub("alpha\\[|\\]", "", name)) |>

    separate(name, into = c("time", "type", "treat")) |>
    mutate(treat = if_else(treat == 1, "ctrl", "omega"))

  sigmas <- draws_sigma |>
    pivot_longer(cols = starts_with("sigma")) |>
    mutate(name = gsub("sigma\\[|\\]", "", name)) |>

    separate(name, into = c("time", "type", "treat")) |>
    mutate(treat = if_else(treat == 1, "ctrl", "omega")) |>
    rename(sigma = value)

  comb_estimate <- estimates |>
    filter(.draw %in% 1:100) |>
    inner_join(sigmas) |>
    rowwise() |>
    mutate(preds = rnorm(1, value, sigma)) |>
    ggplot(aes())
}

#
