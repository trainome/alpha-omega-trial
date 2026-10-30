source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")

# OGTT model mo4, same data as mo3 (08-ogtt.qmd): prior predictive
# (fit.mo4.pp) and posterior (fit.mo4).

# Smoke test (REFIT_SMOKE=1): few iterations, fits saved with a "smoke." prefix
smoke <- Sys.getenv("REFIT_SMOKE") == "1"
fit_name <- function(x) if (smoke) paste0("smoke.", x) else x
n_iter <- if (smoke) 20 else 2000
n_wu <- if (smoke) 20 else 1000

dat4 <- data_load_mo3mo4()

mo4 <- cmdstan_model(stan_file = "models/mo4.stan", exe_file = "models/mo4.exe")

# Prior predictive
fit_wrap(mo4, DATA = dat4, fit = fit_name("fit.mo4.pp"), iter = n_iter, iterwu = n_wu, refit = TRUE)

# Posterior
dat4$prior_only <- 0
fit_wrap(mo4, DATA = dat4, fit = fit_name("fit.mo4"), iter = n_iter, iterwu = n_wu, refit = TRUE)
