source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")

# Fatty acid model with FA x time x treatment-specific sigma (02-fa.qmd):
# posterior (fit.m6p).

# Smoke test (REFIT_SMOKE=1): few iterations, fits saved with a "smoke." prefix
smoke <- Sys.getenv("REFIT_SMOKE") == "1"
fit_name <- function(x) if (smoke) paste0("smoke.", x) else x
n_iter <- if (smoke) 20 else 2000
n_wu <- if (smoke) 20 else 1000

d <- data_load_m3m6()
d$prior_only <- 0

m6 <- cmdstan_model(
  stan_file = "models/m6.stan",
  exe_file = "models/m6.exe",
  compile_model_methods = TRUE
)

fit_wrap(m6, DATA = d, fit = fit_name("fit.m6p"), iter = n_iter, iterwu = n_wu, refit = TRUE)
