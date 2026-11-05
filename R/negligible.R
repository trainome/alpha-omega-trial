# Probability of benefit, negligible effect and harm ##########################
#
# Thresholds for negligible effects (delta) are derived in 01-data.qmd
# (section "Determining thresholds for negligible effects") and stored in
# data-derived/negligible-thresholds.RDS with one row per outcome:
#
#   outcome  name used in the analysis data (analyte), or "ogtt_point"
#            (single minutes of the OGTT), "ogtt_mean" (average across the
#            OGTT), "burden_pop" (main burden score, population SD units) and
#            "burden" (burden score on the model scale, sensitivity analysis)
#   delta    threshold for a negligible effect, on the scale of the contrast:
#            natural units (scale = "identity"), log scale (scale = "log"; TG,
#            HOMA2IR, HbA1c) or SD units of the burden score (scale = "sd";
#            population SD units for "burden_pop", model-scale SD units for
#            "burden")
#   better   direction of benefit, "lower" or "higher"
#
# Contrasts (krill oil - placebo) must be on the same scale as delta. For log
# outcomes this is the difference in logs, i.e. log(ratio), NOT the %
# difference.

neg_file <- "data-derived/negligible-thresholds.RDS"

load_neg <- function(file = neg_file) readRDS(file)


# Probabilities from posterior draws of a contrast (x). Returns a named vector
# c(benefit, negligible, harm) that sums to 1.
p_bnh <- function(x, delta, better = c("lower", "higher")) {
  better <- match.arg(better)
  # Orient so that positive values are beneficial
  s <- if (better == "lower") -x else x
  c(benefit = mean(s > delta),
    negligible = mean(abs(s) <= delta),
    harm = mean(s < -delta))
}


# As p_bnh() with delta and direction looked up by outcome name
p_bnh_outcome <- function(x, outcome, thresholds = load_neg()) {
  th <- thresholds[thresholds$outcome == outcome, ]
  if (nrow(th) != 1) stop("No unique threshold for outcome: ", outcome)
  p_bnh(x, th$delta, th$better)
}


# Formatting. p is the output of p_bnh().
#   "long"  -> "P(benefit) 62 %, P(negligible) 35 %, P(harm) 3 %"
#   "short" -> "62/35/3 %" (benefit/negligible/harm)
fmt_pct_bnh <- function(p) {
  ifelse(p > 0.99, ">99", ifelse(p < 0.01, "<1", sprintf("%.0f", 100 * p)))
}

p_bnh_txt <- function(p, style = c("long", "short")) {
  style <- match.arg(style)
  v <- fmt_pct_bnh(p)
  if (style == "long") {
    sprintf("P(benefit) %s %%, P(negligible) %s %%, P(harm) %s %%",
            v[["benefit"]], v[["negligible"]], v[["harm"]])
  } else {
    sprintf("%s/%s/%s %%", v[["benefit"]], v[["negligible"]], v[["harm"]])
  }
}


# OGTT ########################################################################
# pred: posterior predicted population-level glucose (mmol/L) with columns
# .draw, min_seq, mu and time_treat ("T1_ctrl", ..., "T4_omega"), as created
# in the inference sections of 09-ogtt.qmd.
#
# Treatment effects (krill oil - placebo) as differences in change:
#   T1-T2: supplementation period, T1-T4: from baseline to after the combined
#   period, T2-T4: combined period (primary outcome)
ogtt_effects <- function(pred) {
  pred |>
    tidyr::pivot_wider(values_from = mu, names_from = time_treat) |>
    dplyr::transmute(.draw, min_seq,
                     `T1-T2` = (T2_omega - T1_omega) - (T2_ctrl - T1_ctrl),
                     `T1-T4` = (T4_omega - T1_omega) - (T4_ctrl - T1_ctrl),
                     `T2-T4` = (T4_omega - T2_omega) - (T4_ctrl - T2_ctrl)) |>
    tidyr::pivot_longer(-c(.draw, min_seq), names_to = "contrast", values_to = "value")
}

# P(benefit/negligible/harm) at each minute of the curve (time-point-specific
# delta, "ogtt_point")
ogtt_bnh_curve <- function(eff, thresholds = load_neg()) {
  eff |>
    dplyr::summarise(.by = c(contrast, min_seq),
                     bnh_df(value, "ogtt_point", thresholds))
}

# Average effect across a window of the curve (difference in AUC / window
# length) per draw, with P(benefit/negligible/harm) using the delta for the
# average glucose across the curve ("ogtt_mean")
ogtt_bnh_mean <- function(eff, from = 0, to = 120, thresholds = load_neg()) {
  eff |>
    dplyr::filter(min_seq >= from, min_seq <= to) |>
    dplyr::summarise(.by = c(contrast, .draw), value = mean(value)) |>
    dplyr::summarise(.by = contrast,
                     m = median(value),
                     lwr95 = quantile(value, 0.025),
                     upr95 = quantile(value, 0.975),
                     bnh_df(value, "ogtt_mean", thresholds))
}


# Data frame version for grouped summaries, e.g.
#   contrasts |> summarise(.by = c(analyte, time), bnh_df(contrast_nat, analyte[1]))
bnh_df <- function(x, outcome, thresholds = load_neg()) {
  p <- p_bnh_outcome(x, outcome, thresholds)
  data.frame(benefit = p[["benefit"]],
             negligible = p[["negligible"]],
             harm = p[["harm"]])
}
