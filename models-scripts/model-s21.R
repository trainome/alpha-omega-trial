source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")

refit.models <- TRUE


dat <- data_load_s7s10()

d <- dat$d

# Moderator: baseline relative fat mass. The model takes the baseline value of
# this measurement type from the (partly imputed) baseline matrix and
# standardises it within each draw.
d$mod_type <- dat$scale |>
  filter(analyte == "fat.rel") |>
  pull(type)

# Run model

# Compile model (with methods)
s21 <- cmdstan_model(
  stan_file = "models/s21.stan",
  exe_file = "models/s21.exe",
  compile_model_methods = TRUE,
  force_recompile = TRUE
)


# Fit and save the model with fit wrap. This will be possible to extract again
# in the qmd file.
fit_wrap(
  s21,
  DATA = d,
  fit = "fit.s21.post",
  init = make_inits(s21, DATA = d),
  refit = TRUE,
  ad = 0.95,
  iterwu = 2000
)
# Initialize methods for loo
fit.s21.post$init_model_methods()

# Use loo with moment matching
loo.s21 <- fit.s21.post$loo(
  variables = "log_lik",
  moment_match = TRUE,
  cores = 4
)

saveRDS(loo.s21, "fits/paper1/loo_s21.RDS")
