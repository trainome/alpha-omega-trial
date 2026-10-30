source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")

# OGTT model, complete pre-intervention curves (08-ogtt.qmd): prior predictive
# (fit.mo1.pp) and posterior (fit.mo1).

# Smoke test (REFIT_SMOKE=1): few iterations, fits saved with a "smoke." prefix
smoke <- Sys.getenv("REFIT_SMOKE") == "1"
fit_name <- function(x) if (smoke) paste0("smoke.", x) else x
n_iter <- if (smoke) 20 else 2000
n_wu <- if (smoke) 20 else 1000

dat <- data_load_mo1()

mo1 <- cmdstan_model(stan_file = "models/mo1.stan", exe_file = "models/mo1.exe")

# Prior predictive
fit_wrap(mo1, DATA = dat, fit = fit_name("fit.mo1.pp"), iter = n_iter, iterwu = n_wu, refit = TRUE)

# Posterior
dat$prior_only <- 0
fit_wrap(mo1, DATA = dat, fit = fit_name("fit.mo1"), iter = n_iter, iterwu = n_wu, refit = TRUE)
