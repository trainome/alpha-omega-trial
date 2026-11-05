# Comparing accuracy of s13 and s18 in simulated data sets
#
# s13 (multivariate Student-t, full 2K x 2K residual correlation, the model
# used for inference) is compared with s18 (its univariate Student-t
# counterpart, no correlation between types or time-points). Both use
# time x type residual scales and nu = 4, so the comparison isolates what the
# correlation structure adds to recovery of the true treatment contrasts.
#
# Data are generated with sim_wrap3 (homoskedastic, sigma_ratio = 1, with
# added outlier noise) over a grid of
#   - missingness: participant-level dropout plus values missing completely
#     at random. Approximate proportion of follow-up values missing:
#     low ~10%, mid ~28% (as in the trial data), high ~40%.
#   - correlation strength: weak (near independence; s13 and s18 should
#     agree) and strong (sim_wrap3 defaults; mean |r| between types ~0.27,
#     between time-points ~0.34).
# Baseline (pre) values are ~8% missing in all cells, as in the trial data.
# Within-time and between-time correlation cannot both be made much stronger
# with sim_type3 (R_full stops being positive definite).
#
# Each data set is fit on both the full and the missing variant (comp_fun2),
# with make_inits() starting values. Accuracy is summarised from the saved
# output: error_comp$bias_* is the error of the posterior mean
# (RMSE = sqrt(mean(bias^2))), coverage holds the 90% intervals.
#
# Cells are interleaved over the simulation index so a partial run covers
# all cells evenly. Existing output files are skipped, so an interrupted run
# can be restarted. Output is written to data-sim3/ (est_error_dist/ is
# git-ignored because of its size; the other summaries can be committed).
#
# Usage, from the project root:
#   Rscript R/model-comp3.R            # full run
#   Rscript R/model-comp3.R --smoke    # one simulation in a separate worker,
#                                      # output to a temporary folder
# The number of parallel workers is set with the environment variable
# SIM_WORKERS (default: available cores - 1). Each worker needs ~1-2 GB of
# memory and cmdstan writes a few GB of temporary CSV files per worker, so
# reduce SIM_WORKERS on machines with less memory.
#
# Requirements: CmdStan (cmdstanr::install_cmdstan()) and the R packages
# cmdstanr, posterior, dplyr, tidyr, tidyselect, MASS, furrr. The trial data
# (AOData) is not needed. Model executables are compiled on first use.

if (!file.exists("R/model-comp3.R")) {
  stop("Run from the project root (the folder containing R/ and models/)")
}

library(furrr)
library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
source("./R/helpers.R")
source("./R/sim-fun.R")

smoke <- "--smoke" %in% commandArgs(trailingOnly = TRUE)
out_root <- if (smoke) file.path(tempdir(), "data-sim3-smoke") else "data-sim3"
n_workers <- as.integer(Sys.getenv(
  "SIM_WORKERS",
  max(1, parallelly::availableCores() - 1)
))

# Pre-compile once outside future_walk so parallel workers do not race.
for (m in c("s13", "s18")) {
  cmdstanr::cmdstan_model(
    stan_file = paste0("models/", m, ".stan"),
    exe_file = paste0("models/", m, ".exe"),
    force_recompile = FALSE
  )
}

# Simulation grid
miss_levels <- list(
  low = list(dropout_int = -3, miss_prop_y = 0.01),
  mid = list(dropout_int = -1.7, miss_prop_y = 0.05),
  high = list(dropout_int = -1.25, miss_prop_y = 0.10)
)
corr_levels <- list(
  weak = list(
    loadings_beta = matrix(c(2, 10, 1, 40, 1, 40), ncol = 2, byrow = TRUE),
    c_corr = c(2, 10) # beta(2, 10) shape parameters, mean 0.17
  ),
  strong = list(
    loadings_beta = matrix(c(2, 10, 1, 10, 1, 10), ncol = 2, byrow = TRUE),
    c_corr = c(5, 10) # mean 0.33
  )
)
cells <- expand.grid(
  miss = names(miss_levels),
  corr = names(corr_levels),
  stringsAsFactors = FALSE
)

n_per_cell <- 50
n_sims <- n_per_cell * nrow(cells)

out_dirs <- c(
  "coverage",
  "error_comp",
  "est_error_dist",
  "average_hmc_metrics",
  "sim_info"
)


# Realised correlation of follow-up residuals (after adjusting for pre and
# treatment) in the full data: mean |r| between types within time, and mean
# r between time-points within type.
realised_cor <- function(dl) {
  K <- dl$K
  df <- data.frame(
    y = dl$y,
    pre = dl$pre,
    type = dl$type,
    time = dl$time,
    P = dl$P,
    treat = dl$treat
  )
  df$r <- NA
  for (ty in 1:K) {
    for (ti in 1:2) {
      i <- df$type == ty & df$time == ti
      df$r[i] <- resid(lm(y ~ pre + factor(treat), data = df[i, ]))
    }
  }
  w <- df |>
    select(P, type, time, r) |>
    pivot_wider(names_from = c(time, type), values_from = r)
  w <- as.matrix(w[, paste(rep(1:2, each = K), rep(1:K, 2), sep = "_")])
  cm <- cor(w)
  within <- cm[1:K, 1:K]
  c(
    cor_within = mean(abs(within[lower.tri(within)])),
    cor_time = mean(diag(cm[1:K, (K + 1):(2 * K)]))
  )
}


sim_comp_wrap3 <- function(i) {
  if (file.exists(file.path(out_root, "coverage", paste0("d", i, ".RDS")))) {
    return(invisible(NULL))
  }

  k <- 12 # number of measurement types (analytes)
  cell <- cells[(i - 1) %% nrow(cells) + 1, ]
  ml <- miss_levels[[cell$miss]]
  cl <- corr_levels[[cell$corr]]

  # Get models (compiled executables assumed to exist)
  s13 <- cmdstanr::cmdstan_model(
    exe_file = "models/s13.exe",
    stan_file = "models/s13.stan"
  )
  s18 <- cmdstanr::cmdstan_model(
    exe_file = "models/s18.exe",
    stan_file = "models/s18.stan"
  )
  mods <- list(s13 = s13, s18 = s18)

  # Simulate data. Attempt up to max_tries; stop if none succeed.
  max_tries <- 20
  d <- NULL
  for (attempt in seq_len(max_tries)) {
    d <- tryCatch(
      sim_wrap3(
        # Symmetric around 0 with a null spike so coverage is tested
        # for zero, negative, and positive effects.
        group1 = ifelse(runif(k) < 0.3, 0, runif(k, -0.4, 0.4)),
        group2 = ifelse(runif(k) < 0.3, 0, runif(k, -0.4, 0.4)),
        loadings_beta = cl$loadings_beta,
        c_corr = rbeta(k, cl$c_corr[1], cl$c_corr[2]),
        sigma_ratio = rep(1, k),
        N = 148,
        K = k,
        prop_noise = 0.025,
        sd_noise = 1.5,
        dropout_int = ml$dropout_int,
        miss_prop_y = ml$miss_prop_y,
        miss_prop_pre = 0.04
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
  if (is.null(d)) {
    stop(sprintf("sim %d: no valid data after %d tries", i, max_tries))
  }

  comp_out <- tryCatch(
    # s13/s18 take the standardised MetS thresholds as data (values from
    # diag_thresholds() for the trial data); the diagnosis generated
    # quantities are not used in the simulation comparison.
    comp_fun2(d, mods, extra_data = list(
      nu = 4,
      threshold_male = c(-1.2037, 0.9539, 0.2868, 1.8785, -0.1225, -0.4402),
      threshold_female = c(-0.1962, 0.9539, 0.2868, 1.8785, -0.1225, -1.4842)
    ), use_inits = TRUE),
    error = function(e) {
      message(sprintf("sim %d: comp_fun2 failed: %s", i, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(comp_out)) return(invisible(NULL))

  dm <- d$miss$dl_itt
  comp_out$sim_info <- data.frame(
    sim = i,
    miss = cell$miss,
    corr = cell$corr,
    prop_miss_y = 1 - mean(dm$y_miss),
    prop_miss_pre = 1 - mean(dm$pre_miss),
    t(realised_cor(d$full$dl_itt))
  )

  # Tag every output table with the grid cell
  for (nm in out_dirs) {
    comp_out[[nm]] <- comp_out[[nm]] |>
      mutate(sim = i, miss = cell$miss, corr = cell$corr)
  }

  # Write coverage last: its presence marks the simulation as complete
  for (nm in c(setdiff(out_dirs, "coverage"), "coverage")) {
    saveRDS(
      comp_out[[nm]],
      file = file.path(out_root, nm, paste0("d", i, ".RDS"))
    )
  }
}


# Check if folders exist
for (nm in out_dirs) {
  dir.create(file.path(out_root, nm), recursive = TRUE, showWarnings = FALSE)
}

# Smoke test: one simulation (mid missingness, strong correlation) in a
# separate worker process, so parallel export of functions is also tested.
sims <- if (smoke) 5 else 1:n_sims
if (smoke) n_workers <- 1

set.seed(2026)
plan(multisession, workers = n_workers)
message(sprintf(
  "Running %d simulation(s) on %d worker(s), output to %s",
  length(sims),
  n_workers,
  out_root
))

future_walk(
  sims,
  sim_comp_wrap3,
  .options = furrr_options(
    seed = TRUE,
    packages = c("dplyr", "tidyr", "cmdstanr", "posterior")
  )
)

if (smoke) {
  f <- function(nm) readRDS(file.path(out_root, nm, "d5.RDS"))
  if (!file.exists(file.path(out_root, "coverage", "d5.RDS"))) {
    stop("Smoke test failed: no output written (see messages above)")
  }
  print(f("sim_info"))
  print(
    f("coverage") |>
      summarise(
        .by = c(model, variant),
        coverage = mean(coverage),
        ci_width = mean(ci_width)
      )
  )
  print(
    f("error_comp") |>
      summarise(
        .by = model,
        rmse_full = sqrt(mean(bias_full^2)),
        rmse_miss = sqrt(mean(bias_miss^2)),
        mean_bias_full = mean(bias_full),
        mean_bias_miss = mean(bias_miss)
      )
  )
  print(
    f("average_hmc_metrics") |>
      filter(name %in% c("rhat", "ess_bulk")) |>
      select(name, min, max, model, variant)
  )
  message("Smoke test passed")
}
