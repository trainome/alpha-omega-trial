source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")

# OGTT model, all participants with imputation (08-ogtt.qmd): prior predictive
# (fit.mo3.pp) and posterior (fit.mo3).

# Smoke test (REFIT_SMOKE=1): few iterations, fits saved with a "smoke." prefix
smoke <- Sys.getenv("REFIT_SMOKE") == "1"
fit_name <- function(x) if (smoke) paste0("smoke.", x) else x
n_iter <- if (smoke) 20 else 2000
n_wu <- if (smoke) 20 else 1000

dat <- data_load_mo3mo4()

mo3 <- cmdstan_model(stan_file = "models/mo3.stan", exe_file = "models/mo3.exe")

# Prior predictive
fit_wrap(mo3, DATA = dat, fit = fit_name("fit.mo3.pp"), iter = n_iter, iterwu = n_wu, refit = TRUE)

# Posterior
dat$prior_only <- 0
fit_wrap(mo3, DATA = dat, fit = fit_name("fit.mo3"), iter = n_iter, iterwu = n_wu, refit = TRUE)
