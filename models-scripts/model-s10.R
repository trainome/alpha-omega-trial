source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")

refit.models <- TRUE


dat <- data_load_s7s10()

d <- dat$d

# Run model

# Compile model (with methods)
s10 <- cmdstan_model(
  stan_file = "models/s10.stan",
  exe_file = "models/s10.exe",
  compile_model_methods = TRUE,
  force_recompile = TRUE
)


# Initial values: start the correlation Cholesky factors near the identity so
# chains begin away from the singular (rank-deficient) boundary of L_Omega that
# trapped a chain previously. Other parameters keep cmdstan's random init.
near_identity_chol <- function(dim, jitter = 0.01) {
  L <- diag(dim)
  L[lower.tri(L)] <- rnorm(sum(lower.tri(L)), 0, jitter)
  L / sqrt(rowSums(L^2)) # renormalise rows -> valid corr Cholesky (unit diag)
}
init_s10 <- function() {
  list(
    L_Omega = near_identity_chol(2 * d$K),
    L_Omega_pre = near_identity_chol(d$K)
  )
}

# Fit and save the model with fit wrap. This will be possible to extract again
# in the qmd file.
fit_wrap(
  s10,
  DATA = d,
  fit = "fit.s10.post",
  ad = 0.90,
  refresh = 250,

  chains = 4,

  init = init_s10,
  refit = TRUE
)
# Initialize methods for loo
fit.s10.post$init_model_methods()


# Use loo with moment matching
loo.s10 <- fit.s10.post$loo(
  variables = "log_lik",
  moment_match = TRUE,
  cores = 4
)

saveRDS(loo.s10, "fits/paper1/loo_s10.RDS")
