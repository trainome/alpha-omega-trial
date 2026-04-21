# Comparing model performance between simulated data sets
#
#
#
# Model s7 and s8 are compared against simulated data sets to assess
# precision and bias.
library(furrr)
library(progressr)
handlers(global = TRUE)
handlers("txtprogressbar")
source("./R/fig-opt-fun.R")
source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")
source("./R/sim-fun.R")


sim_comp_wrap <- function(i, p = NULL) {
  if (!is.null(p)) p(sprintf("sim %d", i))
  # Get models
  s7 <- cmdstanr::cmdstan_model(
    exe_file = "models/s7.exe",
    stan_file = "models/s7.stan"
  )
  s8 <- cmdstanr::cmdstan_model(
    exe_file = "models/s8.exe",
    stan_file = "models/s8.stan"
  )

  # Simulate data sets. This will attempt to create the data
  # max_tries times, if no data is relaibly created, stop.
  max_tries <- 20
  d <- NULL
  for (attempt in seq_len(max_tries)) {
    d <- tryCatch(
      sim_wrap2(
        group1 = runif(8, -0.2, 0.5),
        group2 = runif(8, -0.2, 0.5),
        loadings_beta = matrix(c(2, 10, 1, 10, 1, 10), ncol = 2, byrow = TRUE),
        c_corr = rbeta(8, 5, 10),
        N = 148,
        K = 8
      ),
      error = function(e) {
        message(sprintf(
          "sim %d attempt %d failed: %s",
          i,
          attempt,
          conditionMessage(e)
        ))
        NULL
      }
    )
    if (!is.null(d)) break
  }
  if (is.null(d))
    stop(sprintf("sim %d: no valid data after %d tries", i, max_tries))

  # Run fitting and comparison
  comp_out <- comp_fun(d, s7, s8)

  # Save output

  saveRDS(comp_out$coverage, file = paste0("data-sim/coverage/d", i, ".RDS"))
  saveRDS(
    comp_out$error_comp,
    file = paste0("data-sim/error_comp/d", i, ".RDS")
  )
  saveRDS(
    comp_out$est_error_dist,
    file = paste0("data-sim/est_error_dist/d", i, ".RDS")
  )
  saveRDS(
    comp_out$average_hmc_metrics,
    file = paste0("data-sim/average_hmc_metrics/d", i, ".RDS")
  )
}


# Check if folders exists
if (!dir.exists("data-sim/error_comp")) dir.create("data-sim/error_comp")
if (!dir.exists("data-sim/coverage")) dir.create("data-sim/coverage")
if (!dir.exists("data-sim/est_error_dist"))
  dir.create("data-sim/est_error_dist")
if (!dir.exists("data-sim/average_hmc_metrics"))
  dir.create("data-sim/average_hmc_metrics")

n_sims <- 500

plan(multisession, workers = 3)

with_progress({
  p <- progressor(steps = n_sims)
  future_walk(
    1:n_sims,
    sim_comp_wrap,
    p = p,
    .options = furrr_options(
      seed = TRUE,
      packages = c("dplyr", "tidyr", "cmdstanr")
    )
  )
})
