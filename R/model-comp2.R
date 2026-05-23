# Comparing model performance between simulated data sets (heteroskedastic DGP)
#
# Companion to model-comp.R. Fits the same four models (s7, s8, s9, s10)
# but generates data with sim_wrap3, which uses treatment-specific residual
# SDs (shared correlation structure). Under this DGP:
#   - s7, s8 are mis-specified (forced to estimate a shared sigma across
#     treatments that genuinely differ); coverage of treatment contrasts
#     should drop below nominal, especially in cells with sigma_ratio != 1.
#   - s9, s10 are correctly specified; coverage should hit nominal.
#
# Output is written to data-sim2/ so the two DGP scenarios can be
# compared side by side.

library(furrr)
source("./R/fig-opt-fun.R")
source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")
source("./R/sim-fun.R")


sim_comp_wrap2 <- function(i) {
  # Get models
  s7 <- cmdstanr::cmdstan_model(
    exe_file = "models/s7.exe",
    stan_file = "models/s7.stan"
  )
  s8 <- cmdstanr::cmdstan_model(
    exe_file = "models/s8.exe",
    stan_file = "models/s8.stan"
  )
  s9 <- cmdstanr::cmdstan_model(
    exe_file = "models/s9.exe",
    stan_file = "models/s9.stan"
  )
  s10 <- cmdstanr::cmdstan_model(
    exe_file = "models/s10.exe",
    stan_file = "models/s10.stan"
  )

  mods <- list(s7 = s7, s8 = s8, s9 = s9, s10 = s10)

  # Simulate data with treatment-specific SDs
  max_tries <- 20
  d <- NULL
  for (attempt in seq_len(max_tries)) {
    d <- tryCatch(
      sim_wrap3(
        group1 = ifelse(runif(8) < 0.3, 0, runif(8, -0.4, 0.4)),
        group2 = ifelse(runif(8) < 0.3, 0, runif(8, -0.4, 0.4)),
        loadings_beta = matrix(c(2, 10, 1, 10, 1, 10), ncol = 2, byrow = TRUE),
        c_corr = rbeta(8, 5, 10),
        # Mix of homoskedastic (ratio = 1) and heteroskedastic (ratio != 1)
        # cells so recovery can be tested for both regimes within one sim.
        sigma_ratio = c(1.5, 0.6, 1.0, 1.0, 1.3, 0.8, 1.0, 1.0),
        N = 148,
        K = 8,
        prop_noise = 0.025,
        sd_noise = 1.5
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

  comp_out <- tryCatch(
    comp_fun2(d, mods),
    error = function(e) {
      message(sprintf("sim %d: comp_fun2 failed: %s", i, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(comp_out)) return(invisible(NULL))

  saveRDS(comp_out$coverage, file = paste0("data-sim2/coverage/d", i, ".RDS"))
  saveRDS(
    comp_out$error_comp,
    file = paste0("data-sim2/error_comp/d", i, ".RDS")
  )
  saveRDS(
    comp_out$est_error_dist,
    file = paste0("data-sim2/est_error_dist/d", i, ".RDS")
  )
  saveRDS(
    comp_out$average_hmc_metrics,
    file = paste0("data-sim2/average_hmc_metrics/d", i, ".RDS")
  )
}


# Check if folders exist
if (!dir.exists("data-sim2")) dir.create("data-sim2")
if (!dir.exists("data-sim2/error_comp")) dir.create("data-sim2/error_comp")
if (!dir.exists("data-sim2/coverage")) dir.create("data-sim2/coverage")
if (!dir.exists("data-sim2/est_error_dist"))
  dir.create("data-sim2/est_error_dist")
if (!dir.exists("data-sim2/average_hmc_metrics"))
  dir.create("data-sim2/average_hmc_metrics")

n_sims <- 150

plan(multisession, workers = 12)

future_walk(
  1:n_sims,
  sim_comp_wrap2,
  .options = furrr_options(
    seed = TRUE,
    packages = c("dplyr", "tidyr", "cmdstanr")
  )
)
