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
#
# - participant_complete_T2: Vector of participants completed T2
# - participant_complete_T4: Vector of participants completed T4

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
  glucose_dat
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
  filter(analyte %in% c("HDL", "TG", "waist", "glucose", "sys", "dias")) |>
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
  filter(analyte %in% c("HDL", "TG", "waist", "glucose", "sys", "dias")) |>
  select(participant, time, treat, sex, analyte, value) |>
  # Critera measures
  pivot_wider(names_from = analyte, values_from = value) |>
  filter(complete.cases(pick(HDL:dias))) |>
  pivot_longer(cols = HDL:dias) |>

  pivot_wider(names_from = time, values_from = value) |>
  filter(complete.cases(pick(T2))) |>
  distinct(participant) |>
  pull(participant)
