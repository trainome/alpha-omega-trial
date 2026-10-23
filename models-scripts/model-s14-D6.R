source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")

refit.models <- TRUE


dat <- data_load_s7s10()

d <- dat$d

# Number of latent factors (D = 6, compared to D = 4 in model-s14.R)
d$D <- 6

# Run model

# Compile model (with methods). Must be recompiled, otherwise the model methods
# needed by loo (moment matching) are not available. Do not run at the same
# time as model-s14.R, as both write to the same executable.
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
  fit = "fit.s14_D6.post",
  init = make_inits(s14, DATA = d),
  refit = TRUE,
  ad = 0.95,
  iterwu = 2000
)
# Initialize methods for loo
fit.s14_D6.post$init_model_methods()

# Use loo with moment matching
loo.s14_D6 <- fit.s14_D6.post$loo(
  variables = "log_lik",
  moment_match = TRUE,
  cores = 4
)

saveRDS(loo.s14_D6, "fits/paper1/loo_s14_D6.RDS")
