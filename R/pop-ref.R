# Population reference values for the metabolic syndrome criteria ############
#
# Shared by 01-data.qmd (thresholds for negligible effects),
# 07-metsyn-inference.qmd and _methods-results.qmd (main burden score) and
# 08-metsyn-sensitivity.qmd (cMetS and alternative scalings). @tbl-cmets-ref-tromso to -saminor2 in
# 08-metsyn-sensitivity.qmd are built from pop_ref.
#
# Requires criteria_dat (R/data-prep.R) for the glucose reference.
#
# Objects:
#   pop_ref        study- and sex-specific published means and SDs (TG as
#                  ln(mmol/L)); use = FALSE rows are kept for documentation
#   map_approx()   MAP mean and SD from SBP and DBP summary statistics
#   cmets_ref_pop  sex-specific references pooled over studies, incl. MAP.
#                  Fasting glucose uses the trial's own sex-specific baseline
#                  mean and SD, as published estimates were not comparable
#                  (participants with diabetes, non-fasting samples).
#   burden_w       burden score weights (sign-aligned)
#   burden_pop_factors()  per-component factors converting model-scale
#                  contrasts to the re-weighted (population SD) burden score

# MAP = DBP + (SBP - DBP) / 3 = (2 * DBP + SBP) / 3, so
#   mean(MAP) = (2 * mean(DBP) + mean(SBP)) / 3
#   var(MAP)  = (4 * var(DBP) + var(SBP) + 4 * r * sd(DBP) * sd(SBP)) / 9
# where r is the SBP-DBP correlation (not reported in published tables).
map_approx <- function(sbp, sbp_sd, dbp, dbp_sd, r) {
  data.frame(r = r,
             mean = (2 * dbp + sbp) / 3,
             sd = sqrt((4 * dbp_sd^2 + sbp_sd^2 + 4 * r * dbp_sd * sbp_sd) / 9))
}

# Published reference values. TG as ln(mmol/L). Glucose and TG for SAMINOR
# derived from median (IQR) (see 08-metsyn-sensitivity.qmd). Published glucose
# estimates are not used (use = FALSE), see cmets_ref_pop below.
pop_ref <- tibble::tribble(
  ~study,              ~use,   ~sex, ~analyte,  ~n,   ~mean,  ~sd,
  "Tromso7",           TRUE,   "m",  "HDL",     7670, 1.4,    0.4,
  "Tromso7",           TRUE,   "m",  "sys",     7670, 130.9,  17.1,
  "Tromso7",           TRUE,   "m",  "dias",    7670, 78.9,   9.6,
  "Tromso7",           TRUE,   "f",  "HDL",     8896, 1.6,    0.5,
  "Tromso7",           TRUE,   "f",  "sys",     8896, 126.8,  18.3,
  "Tromso7",           TRUE,   "f",  "dias",    8896, 75.5,   10.01,
  "FinHealth2017",     TRUE,   "m",  "waist",   2621, 99.1,   12.9,
  "FinHealth2017",     TRUE,   "m",  "HDL",     2621, 1.3,    0.3,
  "FinHealth2017",     TRUE,   "m",  "TG",      2621, 0.190,  0.656,
  "FinHealth2017",     FALSE,  "m",  "glucose", 2621, 6.0,    1.4,
  "FinHealth2017",     TRUE,   "m",  "sys",     2621, 135.9,  17.9,
  "FinHealth2017",     TRUE,   "m",  "dias",    2621, 81.2,   11.2,
  "FinHealth2017",     TRUE,   "f",  "waist",   3066, 89.6,   14.1,
  "FinHealth2017",     TRUE,   "f",  "HDL",     3066, 1.6,    0.4,
  "FinHealth2017",     TRUE,   "f",  "TG",      3066, 0.135,  0.505,
  "FinHealth2017",     FALSE,  "f",  "glucose", 3066, 5.6,    1.1,
  "FinHealth2017",     TRUE,   "f",  "sys",     3066, 132.8,  20.4,
  "FinHealth2017",     TRUE,   "f",  "dias",    3066, 77.8,   10.6,
  # SAMINOR 1 and 2, non-Sami participants [@michalsen2019]
  "SAMINOR1",          FALSE,  "m",  "waist",   1864, 93.9,   10.2,
  "SAMINOR1",          FALSE,  "m",  "HDL",     1864, 1.28,   0.34,
  "SAMINOR1",          FALSE,  "m",  "TG",      1864, 0.457,  0.524,
  "SAMINOR1",          FALSE,  "m",  "glucose", 1864, 5.41,   0.852,
  "SAMINOR1",          FALSE,  "m",  "sys",     1864, 136.1,  17.6,
  "SAMINOR1",          FALSE,  "m",  "dias",    1864, 78.2,   10.0,
  "SAMINOR1",          FALSE,  "f",  "waist",   2176, 85.6,   12.0,
  "SAMINOR1",          FALSE,  "f",  "HDL",     2176, 1.49,   0.40,
  "SAMINOR1",          FALSE,  "f",  "TG",      2176, 0.300,  0.496,
  "SAMINOR1",          FALSE,  "f",  "glucose", 2176, 5.29,   0.808,
  "SAMINOR1",          FALSE,  "f",  "sys",     2176, 133.0,  20.1,
  "SAMINOR1",          FALSE,  "f",  "dias",    2176, 73.0,   10.5,
  "SAMINOR2",          TRUE,   "m",  "waist",   1571, 100.2,  10.7,
  "SAMINOR2",          TRUE,   "m",  "HDL",     1571, 1.28,   0.38,
  "SAMINOR2",          TRUE,   "m",  "TG",      1571, 0.405,  0.532,
  "SAMINOR2",          FALSE,  "m",  "glucose", 1571, 5.4,    0.815,
  "SAMINOR2",          TRUE,   "m",  "sys",     1571, 135.1,  17.2,
  "SAMINOR2",          TRUE,   "m",  "dias",    1571, 77.8,   9.4,
  "SAMINOR2",          TRUE,   "f",  "waist",   1899, 92.9,   12.0,
  "SAMINOR2",          TRUE,   "f",  "HDL",     1899, 1.55,   0.45,
  "SAMINOR2",          TRUE,   "f",  "TG",      1899, 0.336,  0.469,
  "SAMINOR2",          FALSE,  "f",  "glucose", 1899, 5.2,    0.741,
  "SAMINOR2",          TRUE,   "f",  "sys",     1899, 131.1,  18.6,
  "SAMINOR2",          TRUE,   "f",  "dias",    1899, 72.3,   9.0
)

# MAP per study from SBP and DBP, with SBP-DBP correlation r = 0.8
pop_map <- pop_ref |>
  dplyr::filter(use, analyte %in% c("sys", "dias")) |>
  tidyr::pivot_wider(id_cols = c(study, sex, n), names_from = analyte,
                     values_from = c(mean, sd)) |>
  dplyr::mutate(map_approx(mean_sys, sd_sys, mean_dias, sd_dias, r = 0.8),
                analyte = "map") |>
  dplyr::select(study, sex, analyte, n, mean, sd)

# Fasting glucose: sex-specific baseline mean and SD in the trial sample
glucose_ref_own <- criteria_dat |>
  dplyr::filter(time == "T1", analyte == "glucose", !is.na(value)) |>
  dplyr::summarise(.by = sex, n = dplyr::n(), mean = mean(value), sd = sd(value)) |>
  dplyr::mutate(study = "Present study", analyte = "glucose")

# Sex-specific reference values pooled over studies, weighted by study size:
# mean = sum(n * mean) / sum(n), SD from the pooled within-study variance
# sum((n - 1) * sd^2) / sum(n - 1)
cmets_ref_pop <- dplyr::bind_rows(
  pop_ref |> dplyr::filter(use) |> dplyr::select(-use),
  pop_map,
  glucose_ref_own
) |>
  dplyr::summarise(.by = c(sex, analyte),
                   n_studies = dplyr::n(),
                   mean = sum(n * mean) / sum(n),
                   sd = sqrt(sum((n - 1) * sd^2) / sum(n - 1)),
                   n = sum(n))

# Burden score weights: mean of the five criteria, HDL reversed, blood
# pressure as the average of systolic and diastolic
burden_w <- c(HDL = -1, TG = 1, dias = 0.5, sys = 0.5, glucose = 1, waist = 1) / 5

# Re-weighted burden score contrast. Model contrasts are on the model scale
# (SD_marg, dat$scale). For participant p, the contrast in population SD units
# is sum_k w_k * SD_k^marg * dz_k / SD_k,sex(p)^ref; averaged over the sample
# this is sum_k f_k * dz_k with f_k = w_k * SD_k^marg * mean_p(1 / SD_k,sex(p)^ref)
# (identical to cmets_contrast(..., ref = cmets_ref_pop, w = burden_w), sex =
# "average").
#
# scale, dat$scale; sex, sex of each participant in the model (dat$part$sex)
burden_pop_factors <- function(scale, sex, ref = cmets_ref_pop, w = burden_w) {
  sapply(names(w), function(a) {
    sd_marg <- scale$sd[scale$analyte == a]
    sd_ref <- ref$sd[match(paste(sex, a), paste(ref$sex, ref$analyte))]
    if (anyNA(sd_ref)) stop("Missing reference SD for ", a)
    w[[a]] * sd_marg * mean(1 / sd_ref)
  })
}

# Baseline re-weighted burden score per participant (observed, complete
# criteria only), for the threshold of a negligible effect
burden_pop_baseline <- function(ref = cmets_ref_pop, w = burden_w) {
  criteria_dat |>
    dplyr::filter(time == "T1", analyte %in% names(w), !is.na(value)) |>
    dplyr::mutate(value = dplyr::if_else(analyte == "TG", log(value), value)) |>
    dplyr::inner_join(ref |> dplyr::select(sex, analyte, ref_mean = mean, ref_sd = sd),
                      by = c("sex", "analyte")) |>
    dplyr::mutate(z = (value - ref_mean) / ref_sd) |>
    dplyr::filter(dplyr::n() == length(w), .by = participant) |>
    dplyr::summarise(.by = c(participant, sex), burden = sum(w[analyte] * z))
}
