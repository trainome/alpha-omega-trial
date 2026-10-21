source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")

refit.models <- TRUE


dat <- data_load_s7s10()

d <- dat$d

# Run model

# Compile model (with methods)
s9 <- cmdstan_model(
  stan_file = "models/s9.stan",
  exe_file = "models/s9.exe",
  compile_model_methods = TRUE,
  force_recompile = TRUE
)


# Fit and save the model with fit wrap. This will be possible to extract again
# in the qmd file.
fit_wrap(
  s9,
  DATA = d,
  fit = "fit.s9.post",
  init = make_inits(s9, DATA = d),
  refit = TRUE
)
# Initialize methods for loo
fit.s9.post$init_model_methods()

# Use loo with moment matching
loo.s9 <- fit.s9.post$loo(
  variables = "log_lik",
  moment_match = TRUE,
  cores = 4
)

saveRDS(loo.s9, "fits/paper1/loo_s9.RDS")
