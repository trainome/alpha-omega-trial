# Helper functions ############################################################

# Wrapper function for sampling from Stan models.
#
# When models are run and fits are stored the function retrieves the fit
# from storage. This will make it easier to work with the qmd document.
#
#
# model, a compiled Stan model object
# DATA, the data set used to fit the model
# fit, name of the fitted model
# iter, number of iterations
# iterwu, number of iterations in war-up
# chains, number of chains
# pchains, number of parrallel chains,
# refresh, progress refresh
fit_wrap <- function(
  model,
  DATA = d,
  fit = "fit.m1",
  iter = 2000,
  iterwu = 1000,
  chains = 4,
  pchains = 4,
  refresh = 1000,
  refit = refit.models
) {
  fit.filepath <- paste0("fits/paper1/", fit, ".RDS")

  # If refit is set to true in the document, all model fits will be refitted
  if (refit) file.remove(fit.filepath)

  # Sample from the model if a saved version does not exixts.
  if (!file.exists(fit.filepath)) {
    fittedmodel <- model$sample(
      data = DATA, # The list of variables
      iter_sampling = iter, # Number of iterations for sampling
      iter_warmup = iterwu, # Number of warm-up iterations
      chains = chains, # Number of Markov chains
      parallel_chains = pchains, # Number of parallel chains
      refresh = refresh
    ) # Number of iterations before each print
    # Save the model for easy access
    fittedmodel$save_object(file = fit.filepath) # <-- use this
  } else {
    fittedmodel <- readRDS(fit.filepath)
  }
  # Assign the model object to the environment.
  assign(fit, fittedmodel, envir = parent.frame())
}

#
#
#
#
#
#
# Fit and evaluate models using simulated data
#
#
#
# The function assesses models s7 and s8 and compares them.
# The function expects s7 and s8 to be compiled models in the environment.

# d, a data set from sim_wrap2 function

comp_fun <- function(d, s7, s8) {
  # d is a list with $full and $miss (each a struct_sim output) and
  # $true_effects. Each model is fit on both variants so that marginal
  # group contrasts can be compared with vs. without missingness.

  # Compile true effects. These are treatment differences.
  true_effects <- data.frame(
    time = as.character(c(rep(1, 8), rep(2, 8))),
    type = as.character(rep(1:8, 2)),
    true = c(d$true_effects[1, ], d$true_effects[2, ])
  )

  fit_model <- function(mod, dat) {
    dat$prior_only <- 0
    mod$sample(
      data = dat,
      chains = 4,
      iter_warmup = 500,
      iter_sampling = 500,
      parallel_chains = 4,
      refresh = 0,
      show_messages = FALSE
    )
  }

  fit.s7_full <- fit_model(s7, d$full$dl_itt)
  fit.s7_miss <- fit_model(s7, d$miss$dl_itt)
  fit.s8_full <- fit_model(s8, d$full$dl_itt)
  fit.s8_miss <- fit_model(s8, d$miss$dl_itt)

  # HMC summaries per model × variant
  summarise_hmc <- function(fit, model, variant) {
    fit$summary() |>
      filter(variable != "lp__") |>
      select(rhat:ess_tail) |>
      pivot_longer(cols = everything()) |>
      summarise(
        .by = name,
        min = min(value, na.rm = TRUE),
        m = mean(value, na.rm = TRUE),
        max = max(value, na.rm = TRUE)
      ) |>
      mutate(model = model, variant = variant)
  }

  sums <- bind_rows(
    summarise_hmc(fit.s7_full, "s7", "full"),
    summarise_hmc(fit.s7_miss, "s7", "miss"),
    summarise_hmc(fit.s8_full, "s8", "full"),
    summarise_hmc(fit.s8_miss, "s8", "miss")
  )

  # Posterior differences between treatments and errors
  extract_treat <- function(fit) {
    fit$draws("alpha", format = "df") |>
      pivot_longer(cols = starts_with("alpha")) |>
      mutate(name = gsub("alpha\\[|\\]", "", name)) |>
      separate(name, into = c("time", "type", "treat")) |>
      mutate(treat = paste0("treat", treat)) |>
      pivot_wider(names_from = treat, values_from = value) |>
      mutate(contrast = treat2 - treat1) |>
      inner_join(true_effects) |>
      mutate(
        error = contrast - true,
        abs_dev = abs(contrast - true),
        squared_error = error^2
      )
  }

  mets <- list(
    s7_full = extract_treat(fit.s7_full),
    s7_miss = extract_treat(fit.s7_miss),
    s8_full = extract_treat(fit.s8_full),
    s8_miss = extract_treat(fit.s8_miss)
  )

  # Coverage metrics — one row per model × variant × time × type
  coverage <- bind_rows(lapply(names(mets), function(k) {
    parts <- strsplit(k, "_")[[1]]
    mets[[k]] |>
      summarise(
        .by = c(time, type),
        ci_lwr = quantile(contrast, 0.05),
        ci_upr = quantile(contrast, 0.95)
      ) |>
      inner_join(true_effects) |>
      mutate(
        coverage = if_else(true >= ci_lwr & true <= ci_upr, 1, 0),
        ci_width = ci_upr - ci_lwr,
        model = parts[1],
        variant = parts[2]
      )
  }))

  # Error summaries per model × variant, plus miss-vs-full deltas per model
  summarise_err <- function(df, model, variant) {
    df |>
      summarise(
        .by = c(time, type),
        est = mean(contrast),
        bias = mean(error),
        abs_dev = mean(abs_dev),
        sq_err = mean(squared_error)
      ) |>
      mutate(model = model, variant = variant)
  }

  err_long <- bind_rows(
    summarise_err(mets$s7_full, "s7", "full"),
    summarise_err(mets$s7_miss, "s7", "miss"),
    summarise_err(mets$s8_full, "s8", "full"),
    summarise_err(mets$s8_miss, "s8", "miss")
  ) |>
    inner_join(true_effects)

  # Per-model miss-vs-full contrast: does missingness shift estimates?
  error_comp <- err_long |>
    pivot_wider(
      id_cols = c(time, type, model, true),
      names_from = variant,
      values_from = c(est, bias, abs_dev, sq_err)
    ) |>
    mutate(
      delta_est = est_miss - est_full,
      delta_bias = bias_miss - bias_full,
      delta_abs = abs_dev_miss - abs_dev_full,
      delta_sq = sq_err_miss - sq_err_full
    )

  # Estimate/error draws, tagged with model × variant
  est_error_dist <- bind_rows(
    mets$s7_full |> mutate(model = "s7", variant = "full"),
    mets$s7_miss |> mutate(model = "s7", variant = "miss"),
    mets$s8_full |> mutate(model = "s8", variant = "full"),
    mets$s8_miss |> mutate(model = "s8", variant = "miss")
  )

  out <- list(
    coverage = coverage,
    error_comp = error_comp,
    est_error_dist = est_error_dist,
    average_hmc_metrics = sums
  )

  return(out)
}
