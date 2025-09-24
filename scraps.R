



# Waist circumference: >40 inches (102 cm) for men, >35 inches (88 cm) for women.
# Triglycerides: 150 mg/dL or higher, or taking medication for it.
# HDL Cholesterol: <40 mg/dL for men, <50 mg/dL for women, or taking medication for it.
# Blood Pressure: 130/85 mm Hg or higher, or taking medication for it.
# Fasting Glucose: 100 mg/dL or higher, or taking medication for it.


AOData::bloodvar |> 
  filter(analyte == "HDL")


AOData::glucose |> 
  filter(min == 0) |> 
  filter(!is.na(glucose)) |> 
  select(participant, time, glucose) |>
  # Convert to mg/dl
  mutate()
  print()


AOData::bloodvar |> 
  filter(analyte == "TG") |> 
  # Convert to mg/dl
  mutate(TG = if_else(value * 88.57 > 150, TRUE, FALSE)) |> 
  print()


AOData::waist_bp |> 
  filter(measure == "waist") |> 
  mutate(waist = if_else(sex == "m", if_else(value > 102, TRUE, FALSE), 
                         if_else(value > 88, TRUE, FALSE)))

AOData::waist_bp |> 
  filter(measure != "waist") |> 
  pivot_wider(names_from = measure, values_from = value) |> 
  mutate(bp = if_else(sys > 130 | dias > 85, TRUE, FALSE)) |> 
  print()
