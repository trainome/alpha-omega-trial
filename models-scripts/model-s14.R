source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")

refit.models <- TRUE


dat <- data_load_s7s10()

d <- dat$d

# Number of latent factors (D = 4 is the established sweet spot)
d$D <- 4

# Run model

# Compile model (with methods)
s14 <- cmdstan_model(
  stan_file = "models/s14.stan",
  exe_file = "models/s14.exe",
  compile_model_methods = TRUE,
  force_recompile = TRUE
)


# Fit and save the model with fit wrap. This will be possible to extract again
# in the qmd file.
fit_wrap(
  s14,
  DATA = d,
  fit = "fit.s14.post",
  init = make_inits(s14, DATA = d),
  refit = TRUE,
  ad = 0.95,
  iterwu = 2000
)
# Initialize methods for loo
fit.s14.post$init_model_methods()

# Use loo with moment matching
loo.s14 <- fit.s14.post$loo(
  variables = "log_lik",
  moment_match = TRUE,
  cores = 4
)

saveRDS(loo.s14, "fits/paper1/loo_s14.RDS")
