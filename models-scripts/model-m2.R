source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")

# Fatty acid model, complete cases (02-fa.qmd): prior predictive (fit.m2) and
# posterior (fit.m2p).

# Smoke test (REFIT_SMOKE=1): few iterations, fits saved with a "smoke." prefix
smoke <- Sys.getenv("REFIT_SMOKE") == "1"
fit_name <- function(x) if (smoke) paste0("smoke.", x) else x
n_iter <- if (smoke) 20 else 2000
n_wu <- if (smoke) 20 else 1000

d <- data_load_m2()

m2 <- cmdstan_model(stan_file = "models/m2.stan", exe_file = "models/m2.exe")

# Prior predictive
fit_wrap(m2, DATA = d, fit = fit_name("fit.m2"), iter = n_iter, iterwu = n_wu, refit = TRUE)

# Posterior
d$prior_only <- 0
fit_wrap(m2, DATA = d, fit = fit_name("fit.m2p"), iter = n_iter, iterwu = n_wu, refit = TRUE)
