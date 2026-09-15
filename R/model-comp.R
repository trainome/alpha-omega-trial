# Comparing model performance between simulated data sets (homoskedastic DGP)
#
# Models s7, s8, s9, s10 are compared against simulated data sets to assess
# precision, bias, and coverage of treatment contrasts. Data are generated
# with sim_wrap3 and sigma_ratio = 1 in every cell, so the residual SD is
# shared across treatment arms within each (type, time) cell. Under this DGP:
#   - s7, s8 are correctly specified.
#   - s9, s10 are over-specified (extra sigma cells that collapse under the
#     exponential(1) prior).
#
# The companion script model-comp2.R tests recovery under a
# heteroskedastic DGP (sim_wrap3 with sigma_ratio != 1) with the same four
# models. Both scripts use the same simulation function and the same added
# observation noise; only sigma_ratio differs.
#
# Output is written to data-sim/. The earlier s7-vs-s8 output has been
# moved to data-sim-legacy/.

library(furrr)
source("./R/fig-opt-fun.R")
source("./R/packages.R")
source("./R/data-prep.R")
source("./R/helpers.R")
source("./R/sim-fun.R")

for (m in c("s7", "s8", "s9", "s10")) {
  cmdstanr::cmdstan_model(
    stan_file = paste0("models/", m, ".stan"),
    exe_file = paste0("models/", m, ".exe"),
    force_recompile = FALSE
  )
}


sim_comp_wrap <- function(i) {
  k <- 11 # number of measurement types (analytes)
  # Get models (compiled executables assumed to exist; pre-compile once
  # outside future_walk to avoid parallel workers racing).
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

  # Simulate data sets. Homoskedastic: sigma_ratio = 1 for all types.
  # Attempt up to max_tries; stop if none succeed.
  max_tries <- 20
  d <- NULL
  for (attempt in seq_len(max_tries)) {
    d <- tryCatch(
      sim_wrap3(
        # Symmetric around 0 with a null spike so coverage is tested
        # for zero, negative, and positive effects.
        group1 = ifelse(runif(k) < 0.3, 0, runif(k, -0.4, 0.4)),
        group2 = ifelse(runif(k) < 0.3, 0, runif(k, -0.4, 0.4)),
        loadings_beta = matrix(c(2, 10, 1, 10, 1, 10), ncol = 2, byrow = TRUE),
        c_corr = rbeta(k, 5, 10),
        sigma_ratio = rep(1, k),
        N = 148,
        K = k,
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

  # Run fitting and comparison -- guard against fit-level failures
  # (e.g. a dead chain, divergent pathology) so one bad sim doesn't
  # take down the whole future_walk.
  comp_out <- tryCatch(
    comp_fun2(d, mods),
    error = function(e) {
      message(sprintf("sim %d: comp_fun2 failed: %s", i, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(comp_out)) return(invisible(NULL))

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


# Check if folders exist
if (!dir.exists("data-sim")) dir.create("data-sim")
if (!dir.exists("data-sim/error_comp")) dir.create("data-sim/error_comp")
if (!dir.exists("data-sim/coverage")) dir.create("data-sim/coverage")
if (!dir.exists("data-sim/est_error_dist"))
  dir.create("data-sim/est_error_dist")
if (!dir.exists("data-sim/average_hmc_metrics"))
  dir.create("data-sim/average_hmc_metrics")

n_sims <- 150

plan(multisession, workers = 12)

future_walk(
  1:n_sims,
  sim_comp_wrap,
  .options = furrr_options(
    seed = TRUE,
    packages = c("dplyr", "tidyr", "cmdstanr")
  )
)
