source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")

refit.models <- TRUE


dat <- data_load_s7s10()

d <- dat$d

# Run model

# Compile model (with methods)
s18 <- cmdstan_model(
  stan_file = "models/s18.stan",
  exe_file = "models/s18.exe",
  compile_model_methods = TRUE,
  force_recompile = TRUE
)


# Fit and save the model with fit wrap. This will be possible to extract again
# in the qmd file.
fit_wrap(s18, DATA = d, fit = "fit.s18.post", refit = TRUE)
# Initialize methods for loo
fit.s18.post$init_model_methods()

# Use loo with moment matching
loo.s18 <- fit.s18.post$loo(
  variables = "log_lik",
  moment_match = TRUE,
  cores = 4
)

saveRDS(loo.s18, "fits/paper1/loo_s18.RDS")
