# Data prep. ####################################################
# This script wrangles the main outcome variables for the met-syn
# analysis. Sourcing the script will place the following objects in
# the env.
#
# - criteria_dat: data including outcomes
#   * Metabolic syndrome (with sub-scores) (primary)
#     * Triglycerides
#     * Waist circumference
#     * Blood pressure (diastolic and systolic)
#     * HDL Cholesterol
#     * Resting glucose
#   * Body mass index
#   * Relative fat mass
#   * OGTT
#   * HOMA2-IR
#   * Adipondektin
#   * HbA1c
#
# - participant_complete_T2: Vector of participants completed T2
# - participant_complete_T4: Vector of participants completed T4

# The data load function preps observed data for models s7 to s10->
#
#
#

source("./R/packages.R")

## Body composition
part <- AOData::idkeys

## Descriptives data
baseline_dat <- AOData::bodycomp |>
  select(participant, time, lean.whole, BMC.whole, fat.whole) |>
  mutate(weight = lean.whole + BMC.whole + fat.whole) |>
  inner_join(part) |>
  filter(time == "T1") |>
  ## Calculate height based on weight (DXA) and BMI
  mutate(height = 100 * sqrt((weight / 1000) / bmi)) |>
  ## Calculate relative lean and fat
  mutate(
    lean.rel = 100 * (lean.whole / weight),
    fat.rel = 100 * (fat.whole / weight),
    # Weight in kg
    weight = weight / 1000
  ) |>
  select(participant, treat, completed, sex, age, height, weight, bmi) |>
  pivot_longer(cols = age:bmi)


## Included participants
total_participants <- baseline_dat |>
  distinct(participant) |>
  pull(participant)


body_comp <- AOData::bodycomp |>
  select(participant, time, lean.whole, BMC.whole, fat.whole) |>
  mutate(weight = lean.whole + BMC.whole + fat.whole) |>
  inner_join(part) |>
  filter(time %in% c("T1", "T2", "T4")) |>
  ## Calculate height based on weight (DXA) and BMI
  mutate(
    height = if_else(time == "T1", 100 * sqrt((weight / 1000) / bmi), NA)
  ) |>
  mutate(.by = participant, height = mean(height, na.rm = TRUE)) |>

  ## Calculate relative lean and fat
  mutate(
    lean.rel = 100 * (lean.whole / weight),
    fat.rel = 100 * (fat.whole / weight),
    # Weight in kg
    weight = weight / 1000,
    bmi = weight * (height / 100)^-2
  ) |>

  select(participant, time, sex, fat.rel, bmi) |>
  pivot_longer(cols = fat.rel:bmi, names_to = "analyte", values_to = "value") |>
  mutate(threshold = NA)


### Glucose data
glucose_dat <- AOData::glucose |>
  inner_join(
    idkeys |>
      select(participant, sex)
  ) |>
  filter(min == 0) |>
  select(participant, sex, time, HOMA2IR, Adiponectin) |>
  pivot_longer(
    cols = c(HOMA2IR, Adiponectin),
    names_to = "analyte",
    values_to = "value"
  ) |>
  mutate(threshold = NA)


### HbA1c data 
hba1c_dat <- AOData::bloodvar |>
  inner_join(
    AOData::idkeys |>
      select(participant, sex)
  ) |>
  filter(analyte == "HbA1c") |>
  mutate(threshold = NA) |>
  select(participant, time, sex, analyte, value, threshold)





# Combine all data from metabolic syndrome criteria variables
criteria_dat <- bind_rows(
  # HDL Cholesterol: <40 mg/dL for men, <50 mg/dL for women, or taking medication for it.
  AOData::bloodvar |>
    inner_join(
      AOData::idkeys |>
        select(participant, sex)
    ) |>
    filter(analyte == "HDL") |>
    # Convert to mg/dl: 1 mmol/L = 38,6654 mg/dL
    # Keeping on mmol/L
    mutate(
      hdl = if_else(
        sex == "m",
        if_else(value < 1.03, TRUE, FALSE),
        if_else(value < 1.29, TRUE, FALSE)
      )
    ) |>
    select(participant, time, sex, analyte, value, threshold = hdl),

  # Fasting Glucose: 100 mg/dL or higher, or taking medication for it.

  AOData::glucose |>
    filter(min == 0) |>
    filter(!is.na(glucose)) |>
    select(participant, time, glucose) |>
    inner_join(
      AOData::idkeys |>
        select(participant, sex)
    ) |>
    # Convert to mg/dl: 1 mmol/L = 18,0156 mg/dL
    # Keep on mmol/L glucose
    mutate(
      threshold = if_else(glucose >= 5.6, TRUE, FALSE),
      analyte = "glucose"
    ) |>
    select(participant, time, sex, analyte, value = glucose, threshold),

  # Triglycerides: 150 mg/dL or higher, or taking medication for it.

  AOData::bloodvar |>
    filter(analyte == "TG") |>
    inner_join(
      AOData::idkeys |>
        select(participant, sex)
    ) |>
    # Convert to mg/dl: 1 mmol/L = 88.57 mg/dL

    mutate(tg = if_else(value > 1.7, TRUE, FALSE)) |>
    select(participant, time, sex, analyte, value, threshold = tg),

  # Waist circumference: >40 inches (94 cm) for men,
  # >35 inches (80 cm) for women.

  AOData::waist_bp |>
    filter(measure == "waist") |>
    mutate(
      waist = if_else(
        sex == "m",
        if_else(value > 94, TRUE, FALSE),
        if_else(value > 80, TRUE, FALSE)
      )
    ) |>
    select(
      participant,
      time,
      sex,
      analyte = measure,
      value,
      threshold = waist
    ) |>
    filter(!is.na(value)),

  # Blood Pressure: 130/85 mm Hg or higher, or taking medication for it.
  #
  #

  AOData::waist_bp |>
    filter(measure != "waist") |>
    pivot_wider(names_from = measure, values_from = value) |>

    pivot_longer(
      cols = c(sys, dias),
      names_to = "analyte",
      values_to = "value"
    ) |>
    mutate(
      threshold = if_else(
        analyte == "sys" & value > 130 | analyte == "dias" & value > 85,
        TRUE,
        FALSE
      )
    ) |>
    select(participant, time, sex, analyte, value, threshold) |>
    filter(!is.na(value)),
  # Add body composition
  body_comp,
  # Add homa2IR and adipondectin
  glucose_dat,
  # Add HbA1c
  hba1c_dat,
  # Add LDL (no threshold)
  AOData::bloodvar |>
    inner_join(
      AOData::idkeys |>
        select(participant, sex)
    ) |>
    filter(analyte == "LDL") |>
    # Convert to mg/dl: 1 mmol/L = 38,6654 mg/dL
    # Keeping on mmol/L
    mutate(
      threshold = NA
    ) |>
    select(participant, time, sex, analyte, value, threshold),
) |>
  # Create binomial variable for individual criteria
  mutate(crit = if_else(threshold == TRUE, 1, 0)) |>
  inner_join(
    idkeys |>
      select(participant, treat, age, sex)
  ) |>

  # Filter to correspond to total participants
  filter(participant %in% total_participants)


# Data set for complete-case analysis ##########################################
# These are participants with complete cases for the main outcome (met syndrome)
# from T1 to T4. No missing values are allowed in the complete case analysis.
# For the missingness analysis, participants completing T2 are also sorted out.
participant_complete_T4 <- criteria_dat |>
  filter(
    analyte %in% c("HDL", "TG", "waist", "glucose", "sys", "dias")
  ) |>
  select(participant, time, treat, sex, analyte, value) |>
  # Critera measures
  pivot_wider(names_from = analyte, values_from = value) |>
  filter(complete.cases(pick(HDL:dias))) |>
  pivot_longer(cols = HDL:dias) |>

  pivot_wider(names_from = time, values_from = value) |>
  filter(complete.cases(pick(T2, T4))) |>
  distinct(participant) |>
  pull(participant)


participant_complete_T2 <- criteria_dat |>
  filter(
    analyte %in% c("HDL", "TG", "waist", "glucose", "sys", "dias")
  ) |>
  select(participant, time, treat, sex, analyte, value) |>
  # Critera measures
  pivot_wider(names_from = analyte, values_from = value) |>
  filter(complete.cases(pick(HDL:dias))) |>
  pivot_longer(cols = HDL:dias) |>

  pivot_wider(names_from = time, values_from = value) |>
  filter(complete.cases(pick(T2))) |>
  distinct(participant) |>
  pull(participant)


# Data load function - Prepare data for models s7 -
#
#
#
#
data_load_s7s10 <- function() {
  d_full <- criteria_dat |>
    select(participant, time, treat, sex, age, analyte, value) |>
    # Add missing values for all missing obs.
    pivot_wider(names_from = time, values_from = value) |>

    pivot_longer(
      cols = c(T1, T2, T4),
      values_to = "value",
      names_to = "time"
    ) |>

    pivot_wider(names_from = analyte, values_from = value) |>

    pivot_longer(cols = HDL:LDL, values_to = "value", names_to = "analyte") |>

    mutate(
      value_std = if_else(analyte %in% c("TG", "HOMA2IR", "HbA1c"), log(value), value)
    )

  # Scale parameters used for standardisation (and for back-transforming
  # posteriors to the natural scale). Pooled over all time-points, treatments
  # and participants (na.rm). TG and HOMA2IR are standardised on the log scale:
  #   value_std = (f(value) - mean) / sd,  f = log for TG/HOMA2IR/HbA1c, otherwise
  # identity.
  # Back-transform: f(value) = mean + sd * value_std; value = exp(.) if log.
  # A treatment difference on the std scale, delta_std, is sd * delta_std on
  # the natural (or log) scale; exp(sd * delta_std) is the ratio for TG/HOMA2IR.
  # type matches the type index used in the Stan data (dense_rank(analyte)).
  scale_pars <- d_full |>
    summarise(
      .by = analyte,
      mean = mean(value_std, na.rm = TRUE),
      sd = sd(value_std, na.rm = TRUE)
    ) |>
    mutate(
      type = dense_rank(analyte),
      transform = if_else(analyte %in% c("TG", "HOMA2IR", "HbA1c"), "log", "identity")
    ) |>
    select(type, analyte, transform, mean, sd) |>
    arrange(type)

  d_full <- d_full |>
    mutate(
      .by = analyte,
      value_std = (value_std - mean(value_std, na.rm = TRUE)) /
        sd(value_std, na.rm = TRUE)
    ) |>

    mutate(participant = paste0("p", participant))

  # Including all data
  dat <- d_full |>
    filter(time %in% c("T1", "T2", "T4")) |>

    mutate(
      type = dense_rank(analyte),
      participant = dense_rank(participant),
      treat_seq = dense_rank(treat),
      value_std = if_else(is.na(value_std), -999, value_std)
    ) |>

    select(-value) |>

    pivot_wider(names_from = time, values_from = c(value_std)) |>

    pivot_longer(cols = T2:T4, names_to = "time", values_to = "y_obs") |>
    rename(pre_obs = T1) |>

    mutate(
      time = if_else(time == "T2", 1, 2),
      sex = if_else(sex == "m", -0.5, 0.5),
      y_miss = if_else(y_obs == -999, 0, 1),
      pre_miss = if_else(pre_obs == -999, 0, 1)
    ) |>
    # Complete cases for log likelihood calculations
    mutate(.by = participant, n_y = sum(y_miss), n_pre = sum(pre_miss)) |>
    mutate(cc = if_else(n_y == max(n_y) & n_pre == max(n_pre), 1, 0)) |>

    select(
      y_obs,
      pre_obs,
      y_miss,
      pre_miss,
      P = participant,
      sex,
      type,
      type,
      time,
      treat = treat_seq,
      cc
    ) |>

    arrange(P, time, type)

  # Prepare a data list for models s7-s10
  d <- list(
    N = nrow(dat),
    N_P = length(unique(dat$P)),
    N_TX = length(unique(dat$treat)),
    K = length(unique(dat$type)),
    y_obs = dat$y_obs,
    pre_obs = dat$pre_obs,
    y_miss = dat$y_miss,
    pre_miss = dat$pre_miss,
    type = dat$type,
    treat = dat$treat,
    P = dat$P,
    prior_only = 0,
    # Adding sex and time to data structure
    time = dat$time,
    sex = dat$sex,
    sex_idx = if_else(dat$sex < 0, 1L, 2L),
    # Adding nu as a fixed parameter for model s12
    nu = 4,
    # Adding number of latent factors for models s14-s17
    D = 3,
    # Using only complete cases for log-lik (1 when all data observed)
    cc = dat$cc
  )

  # Treatment index lookup (treat_seq = dense_rank(treat)) for interpretation
  treat_lookup <- d_full |>
    distinct(treat) |>
    mutate(treat_idx = dense_rank(treat)) |>
    arrange(treat_idx)

  # Participant lookup (P = dense_rank(participant), as in dat)
  part_lookup <- d_full |>
    distinct(participant, sex, age, treat) |>
    mutate(P = dense_rank(participant)) |>
    select(P, participant, sex, age, treat) |>
    arrange(P)

  # d: Stan data list; dat: long data; scale: per-type mean/sd for
  # back-transformation to the natural scale; treat: treatment index lookup;
  # part: participant lookup (P, sex, age, treat).
  return(list(d = d, dat = dat, scale = scale_pars, treat = treat_lookup,
              part = part_lookup))
}


# Data load functions for models fitted in model-scripts/ ######################
# Each function mirrors the data preparation in the corresponding qmd chunk
# and returns the Stan data list (prior_only = 1, set to 0 for the posterior).


# Drop-out model m1 (03-metsyn-intro.qmd, chunk missing-data-prep)
data_load_m1 <- function() {
  miss_dat <- criteria_dat |>
    select(-threshold, -crit) |>
    # Filter for time-point 1
    filter(time == "T1") |>
    pivot_wider(names_from = analyte, values_from = value) |>
    pivot_longer(cols = age:LDL, values_to = "value", names_to = "analyte") |>
    mutate(
      .by = analyte,
      value_std = if_else(analyte %in% c("TG", "HOMA2IR", "HbA1c"), log(value), value),
      value_std = (value_std - mean(value_std, na.rm = TRUE)) /
        sd(value_std, na.rm = TRUE)
    )

  miss_dat2 <- miss_dat |>
    filter(!is.na(value_std)) |>
    # Missing at follow-up
    mutate(
      miss_t4 = if_else(participant %in% participant_complete_T4, 0, 1),
      miss_t2 = if_else(participant %in% participant_complete_T2, 0, 1),
      P = dense_rank(participant),
      A = dense_rank(analyte),
      TX = if_else(treat == "ctrl", -0.5, 0.5),
      S = if_else(sex == "m", -0.5, 0.5),
    )

  list(
    N = nrow(miss_dat2),
    N_A = length(unique(miss_dat2$A)),
    N_P = length(unique(miss_dat2$P)),
    N_TX = length(unique(miss_dat2$TX)),
    N_S = length(unique(miss_dat2$S)),
    P = miss_dat2$P,
    A = miss_dat2$A,
    TX = miss_dat2$TX,
    S = miss_dat2$S,
    x = miss_dat2$value_std,
    do4 = miss_dat2$miss_t4,
    do2 = miss_dat2$miss_t2,
    prior_only = 1
  )
}


# Fatty acid model m2, complete cases (02-fa.qmd, chunk fig-m2-priorpred)
data_load_m2 <- function() {
  fat_dat <- AOData::omega3 |>
    inner_join(AOData::idkeys) |>
    filter(time %in% c("T1", "T2", "T4"), participant %in% total_participants) |>
    select(participant, time, treat, C16:C226n3) |>
    pivot_longer(cols = C16:C226n3, names_to = "fa", values_to = "value") |>
    mutate(
      .by = fa,
      value_std = log(value),
      value_std = (value_std - mean(value_std)) / sd(value_std)
    ) |>
    mutate(FA = dense_rank(fa), TX = dense_rank(treat))

  fat_dat2 <- fat_dat |>
    select(participant, FA, TX, time, value_std) |>
    pivot_wider(names_from = time, values_from = value_std) |>
    filter(complete.cases(pick(T1:T4))) |>
    pivot_longer(names_to = "TI", values_to = "y", cols = T2:T4) |>
    mutate(TI = dense_rank(TI), P = dense_rank(participant))

  list(
    N = nrow(fat_dat2),
    N_FA = length(unique(fat_dat2$FA)),
    N_P = length(unique(fat_dat2$P)),
    N_TX = length(unique(fat_dat2$TX)),
    N_TI = length(unique(fat_dat2$TI)),
    P = fat_dat2$P,
    FA = fat_dat2$FA,
    TX = fat_dat2$TX,
    TI = fat_dat2$TI,
    pre = fat_dat2$T1,
    y = fat_dat2$y,
    prior_only = 1
  )
}


# Fatty acid models m3-m6, imputation of post-treatment values (02-fa.qmd,
# chunks m3-data-prep and fig-m3-prior)
data_load_m3m6 <- function() {
  fat_dat_compl <- AOData::omega3 |>
    inner_join(AOData::idkeys) |>
    filter(time %in% c("T1", "T2", "T4"), participant %in% total_participants) |>
    # The omega-3 index is modelled as an additional "fatty acid" (FA = 12)
    select(participant, time, treat, C16:C226n3, omega3index) |>
    pivot_longer(cols = c(C16:C226n3, omega3index), names_to = "fa", values_to = "value") |>
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

  list(
    N = nrow(fat_dat_compl2),
    N_FA = length(unique(fat_dat_compl2$FA)),
    N_P = length(unique(fat_dat_compl2$P)),
    N_TX = length(unique(fat_dat_compl2$TX)),
    N_TI = length(unique(fat_dat_compl2$TI)),
    P = fat_dat_compl2$P,
    FA = fat_dat_compl2$FA,
    TX = fat_dat_compl2$TX,
    TI = fat_dat_compl2$TI,
    pre = fat_dat_compl2$T1,
    y = fat_dat_compl2$y,
    obs = fat_dat_compl2$obs,
    prior_only = 1
  )
}


# Glucose (OGTT) data with treatment, as defined in 08-ogtt.qmd. Named
# ogtt_glucose to not overwrite glucose_dat above.
ogtt_glucose <- function() {
  AOData::glucose |>
    select(participant:glucose) |>
    inner_join(AOData::idkeys |> select(participant, treat))
}


# OGTT model mo1, complete pre-intervention curves (08-ogtt.qmd, chunks
# baseline-dat-prep and fig-mo1-pp)
data_load_mo1 <- function() {
  pre_dat <- ogtt_glucose() |>
    filter(time == "T1") |>
    filter(!(is.na(glucose))) |>
    mutate(
      stglucose = (glucose - mean(glucose)) / sd(glucose),
      participant = as.factor(participant)
    ) |>
    # Dropping participants with less than 7 measurements
    mutate(.by = participant, n = n()) |>
    filter(n == 7) |>
    mutate(P = dense_rank(participant))

  B <- splines::bs(pre_dat$min, df = 5)

  list(
    y = pre_dat$stglucose,
    N = nrow(pre_dat),
    K = ncol(B),
    B = B,
    N_P = max(pre_dat$P),
    P = pre_dat$P,
    prior_only = 1
  )
}


# OGTT model mo2, complete curves at T1, T2 and T4 (08-ogtt.qmd, chunk
# fig-prior-predictive-mo2)
data_load_mo2 <- function() {
  ml_dat <- ogtt_glucose() |>
    inner_join(AOData::idkeys) |>
    dplyr::select(participant, time, min, glucose, treat) |>
    filter(time %in% c("T1", "T2", "T4")) |>
    filter(!is.na(glucose)) |>
    mutate(.by = c(participant, time), n = n()) |>
    filter(n == 7) |>
    mutate(
      stglucose = (glucose - mean(glucose)) / sd(glucose),
      participant = as.factor(participant),
      tx = paste0(time, treat),
      tx = factor(tx, levels = c("T1ctrl", "T1omega", "T2ctrl", "T2omega", "T4ctrl", "T4omega"))
    ) |>
    mutate(P = dense_rank(participant), G = dense_rank(tx), Tt = dense_rank(time))

  B2 <- splines::bs(ml_dat$min, df = 5)

  list(
    y = ml_dat$stglucose,
    N = nrow(ml_dat),
    K = ncol(B2),
    B = B2,
    N_G = max(ml_dat$G),
    G = ml_dat$G,
    N_P = max(ml_dat$P),
    P = ml_dat$P,
    N_T = max(ml_dat$Tt),
    Tt = ml_dat$Tt,
    prior_only = 1
  )
}


# OGTT models mo3 and mo4, all participants with imputation of missing values
# (08-ogtt.qmd, chunks fig-prior-mo3 and fig-prior-mo4)
data_load_mo3mo4 <- function() {
  part_full <- criteria_dat |>
    distinct(participant) |>
    pull(participant)

  full_dat <- ogtt_glucose() |>
    inner_join(AOData::idkeys) |>
    filter(participant %in% part_full) |>
    dplyr::select(participant, time, min, glucose, treat) |>
    # Wider to find NA -> longer -> complete data set
    pivot_wider(values_from = glucose, names_from = time) |>
    pivot_longer(values_to = "glucose", names_to = "time", cols = T1:T4) |>
    arrange(participant, time, min) |>
    mutate(
      y_obs = (glucose - mean(glucose, na.rm = TRUE)) / sd(glucose, na.rm = TRUE),
      y_obs = if_else(is.na(glucose), -999, y_obs),
      y_miss = if_else(is.na(glucose), 0, 1),
      tx = paste0(time, treat),
      tx = factor(tx, levels = c("T1ctrl", "T1omega", "T2ctrl", "T2omega", "T4ctrl", "T4omega")),
      G = dense_rank(tx),
      P = dense_rank(participant),
      Tt = dense_rank(time)
    )

  B <- splines::bs(full_dat$min, df = 5)

  list(
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
    prior_only = 1
  )
}
