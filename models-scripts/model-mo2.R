source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")

# OGTT model, complete curves at T1, T2 and T4 (08-ogtt.qmd): prior predictive
# (fit.mo2.pp) and posterior (fit.mo2).

# Smoke test (REFIT_SMOKE=1): few iterations, fits saved with a "smoke." prefix
smoke <- Sys.getenv("REFIT_SMOKE") == "1"
fit_name <- function(x) if (smoke) paste0("smoke.", x) else x
n_iter <- if (smoke) 20 else 2000
n_wu <- if (smoke) 20 else 1000

dat2 <- data_load_mo2()

mo2 <- cmdstan_model(stan_file = "models/mo2.stan", exe_file = "models/mo2.exe")

# Prior predictive
fit_wrap(mo2, DATA = dat2, fit = fit_name("fit.mo2.pp"), iter = n_iter, iterwu = n_wu, refit = TRUE)

# Posterior
dat2$prior_only <- 0
fit_wrap(mo2, DATA = dat2, fit = fit_name("fit.mo2"), iter = n_iter, iterwu = n_wu, refit = TRUE)
