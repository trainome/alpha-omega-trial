source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")

refit.models <- TRUE


dat <- data_load_s7s10()

d <- dat$d

# Run model

# Compile model (with methods)
s13 <- cmdstan_model(
  stan_file = "models/s13.stan",
  exe_file = "models/s13.exe",
  compile_model_methods = TRUE,
  force_recompile = TRUE
)


# Fit and save the model with fit wrap. This will be possible to extract again
# in the qmd file.
fit_wrap(
  s13,
  DATA = d,
  fit = "fit.s13.post",
  refit = TRUE,
  ad = 0.95,
  iterwu = 2000
)
# Initialize methods for loo
fit.s13.post$init_model_methods()

# Use loo with moment matching
loo.s13 <- fit.s13.post$loo(
  variables = "log_lik",
  moment_match = TRUE,
  cores = 4
)

saveRDS(loo.s13, "fits/paper1/loo_s13.RDS")
