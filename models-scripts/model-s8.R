source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")

refit.models <- TRUE


dat <- data_load_s7s10()

d <- dat$d

# Run model

# Compile model (with methods)
s8 <- cmdstan_model(
  stan_file = "models/s8.stan",
  exe_file = "models/s8.exe",
  compile_model_methods = TRUE,
  force_recompile = TRUE
)


# Fit and save the model with fit wrap. This will be possible to extract again
# in the qmd file.
fit_wrap(s8, DATA = d, fit = "fit.s8.post", refit = TRUE)
# Initialize methods for loo
fit.s8.post$init_model_methods()

# Use loo with moment matching
loo.s8 <- fit.s8.post$loo(
  variables = "log_lik",
  moment_match = TRUE,
  cores = 4
)

saveRDS(loo.s8, "fits/paper1/loo_s8.RDS")
