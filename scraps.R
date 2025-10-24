



# Waist circumference: >40 inches (102 cm) for men, >35 inches (88 cm) for women.
# Triglycerides: 150 mg/dL or higher, or taking medication for it.
# HDL Cholesterol: <40 mg/dL for men, <50 mg/dL for women, or taking medication for it.
# Blood Pressure: 130/85 mm Hg or higher, or taking medication for it.
# Fasting Glucose: 100 mg/dL or higher, or taking medication for it.

bind_rows(
  AOData::bloodvar |> 
  inner_join(AOData::idkeys |> 
               select(participant, sex)) |> 
  filter(analyte == "HDL") |> 
  # Convert to mg/dl: 1 mmol/L = 38,6654 mg/dL
  mutate(hdl = if_else(sex == "m",
                       if_else(value * 38.6654 < 40, TRUE, FALSE), 
                       if_else(value * 38.6654 < 50, TRUE, FALSE))) |> 
  select(participant, time, sex, analyte, value, threshold = hdl), 
  
  AOData::glucose |> 
  filter(min == 0) |> 
  filter(!is.na(glucose)) |> 
  select(participant, time, glucose) |>
  inner_join(AOData::idkeys |> 
               select(participant, sex)) |> 
  # Convert to mg/dl: 1 mmol/L = 18,0156 mg/dL
  mutate(threshold = if_else(glucose * 18.0156 > 100, TRUE, FALSE), 
         analyte = "glucose") |> 
  select(participant, time, sex, analyte, value = glucose, threshold), 
  
  AOData::bloodvar |> 
  filter(analyte == "TG") |> 
  inner_join(AOData::idkeys |> 
               select(participant, sex)) |> 
  # Convert to mg/dl: 1 mmol/L = 88.57 mg/dL
  mutate(tg = if_else(value * 88.57 > 150, TRUE, FALSE)) |> 
  select(participant, time, sex, analyte, value, threshold = tg), 
  
  AOData::waist_bp |> 
  filter(measure == "waist") |> 
  mutate(waist = if_else(sex == "m", if_else(value > 102, TRUE, FALSE), 
                         if_else(value > 88, TRUE, FALSE))) |> 
  select(participant, time, sex, analyte = measure, value, threshold = waist) |> 
  filter(!is.na(value)),
  
  AOData::waist_bp |> 
  filter(measure != "waist") |> 
  pivot_wider(names_from = measure, values_from = value) |> 
  mutate(bp = if_else(sys > 130 | dias > 85, TRUE, FALSE), 
         map = dias + ((sys-dias)/3), 
         analyte = "bp") |> 
  select(participant, time, sex, analyte, value = map, threshold = bp) |> 
  filter(!is.na(value))
) |> 
  

  summarise(.by = c(participant, time), 
            percent = sum(threshold) / n()) |> 
  filter(percent >=0.6) |> 
  print()
  
  ggplot(aes(n)) + geom_histogram() + 
  facet_wrap(~ time)
  
  

