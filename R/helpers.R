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
# SUPERSEDED: kept for backwards-reproducibility of earlier runs.
# New code should use comp_fun2, which accepts an arbitrary named list of
# models (e.g. list(s7=s7, s8=s8, s9=s9, s10=s10)).

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
    fit <- mod$sample(
      data = dat,
      chains = 3,
      iter_warmup = 500,
      iter_sampling = 750,
      parallel_chains = 1,
      refresh = 0,
      show_messages = FALSE
    )
    fit$draws()
    fit$sampler_diagnostics()
    fit
  }

  fit.s7_full <- fit_model(s7, d$full$dl_itt)
  fit.s7_miss <- fit_model(s7, d$miss$dl_itt)
  fit.s8_full <- fit_model(s8, d$full$dl_itt)
  fit.s8_miss <- fit_model(s8, d$miss$dl_itt)
  on.exit(
    {
      for (f in list(fit.s7_full, fit.s7_miss, fit.s8_full, fit.s8_miss)) {
        try(file.remove(f$output_files()), silent = TRUE)
      }
    },
    add = TRUE
  )

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


# Fit and evaluate an arbitrary set of models using simulated data.
#
# mods, a named list of compiled cmdstan_model objects.
#   Names become the model labels in the output, e.g.
#   list(s7 = s7, s8 = s8, s9 = s9, s10 = s10).
# d,    a simulated data set with $full, $miss, $true_effects
#   (produced by sim_wrap2 or sim_wrap3).
#
# Each model is fit on both the full and missing variants so that
# marginal group contrasts can be compared with vs. without missingness.
# Output structure matches comp_fun (coverage, error_comp, est_error_dist,
# average_hmc_metrics) so downstream plotting code works unchanged.
comp_fun2 <- function(d, mods) {
  # True treatment differences per time x type
  true_effects <- data.frame(
    time = as.character(c(rep(1, 8), rep(2, 8))),
    type = as.character(rep(1:8, 2)),
    true = c(d$true_effects[1, ], d$true_effects[2, ])
  )

  fit_model <- function(mod, dat) {
    dat$prior_only <- 0
    fit <- mod$sample(
      data = dat,
      chains = 3,
      iter_warmup = 500,
      iter_sampling = 750,
      parallel_chains = 1,
      refresh = 0,
      show_messages = FALSE
    )
    fit$draws()
    fit$sampler_diagnostics()
    fit
  }

  variants <- c("full", "miss")
  fit_keys <- as.vector(outer(names(mods), variants, paste, sep = "_"))

  fits <- list()
  for (m in names(mods)) {
    for (v in variants) {
      fits[[paste(m, v, sep = "_")]] <- fit_model(mods[[m]], d[[v]]$dl_itt)
    }
  }
  on.exit(
    {
      for (f in fits) try(file.remove(f$output_files()), silent = TRUE)
    },
    add = TRUE
  )

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

  sums <- bind_rows(lapply(fit_keys, function(k) {
    parts <- strsplit(k, "_")[[1]]
    summarise_hmc(fits[[k]], parts[1], parts[2])
  }))

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

  mets <- setNames(
    lapply(fit_keys, function(k) extract_treat(fits[[k]])),
    fit_keys
  )

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

  err_long <- bind_rows(lapply(names(mets), function(k) {
    parts <- strsplit(k, "_")[[1]]
    summarise_err(mets[[k]], parts[1], parts[2])
  })) |>
    inner_join(true_effects)

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

  est_error_dist <- bind_rows(lapply(names(mets), function(k) {
    parts <- strsplit(k, "_")[[1]]
    mets[[k]] |> mutate(model = parts[1], variant = parts[2])
  }))

  list(
    coverage = coverage,
    error_comp = error_comp,
    est_error_dist = est_error_dist,
    average_hmc_metrics = sums
  )
}


# cor_plot ##################################################################
# A function for producing correlations between variables
# from a full observed (or synthetic) data set using no missing
# values.
#
# d, a data set with participant, treat, sex,
# age, time, analyte, value.
cor_plot <- function(d) {
  compl_temp <- d |>
    filter(!is.na(value_std)) |>
    mutate(.by = participant, n = n()) |>
    filter(n == 24) |>
    select(-value) |>
    pivot_wider(names_from = time, values_from = value_std)

  #
  resids0 <- list()
  resids1 <- list()
  resids2 <- list()
  resids_time <- list()
  sigmas <- list()
  analytes <- unique(compl_temp$analyte)
  for (i in 1:8) {
    temp <- compl_temp |>
      filter(analyte == analytes[i])
    m0 <- lm(T1 ~ sex, data = temp)
    m1 <- lm(T2 ~ T1 + sex + treat, data = temp)
    m2 <- lm(T4 ~ T1 + sex + treat, data = temp)

    resids0[[i]] <- resid(m0)
    resids1[[i]] <- resid(m1)
    resids2[[i]] <- resid(m2)

    resids_time[[i]] <- data.frame(
      analyte = rep(analytes[i], nrow(temp)),
      t1 = resid(m0),
      t2 = resid(m1),
      t4 = resid(m2)
    )

    sigmas[[i]] <- data.frame(
      analyte = analytes[i],
      t1 = sigma(m0),
      t2 = sigma(m1),
      t4 = sigma(m2)
    )
  }

  resids0 <- bind_cols(resids0)
  resids1 <- bind_cols(resids1)
  resids2 <- bind_cols(resids2)

  colnames(resids0) <- analytes
  colnames(resids1) <- analytes
  colnames(resids2) <- analytes

  p1 <- bind_rows(resids_time) |>
    summarise(.by = analyte, c = cor(t2, t4)) |>
    ggplot(aes(c, analyte)) +
    geom_point(size = 3, shape = 21, fill = colscale[1]) +
    alpha_theme() +
    theme(panel.grid.major.y = element_line(color = "gray80")) +
    labs(
      x = "Residual Correlation T2 and T4",
      y = "Analyte",
      subtitle = "Correlation of residuals between T2 and T4"
    )

  p2 <- bind_rows(sigmas) |>
    pivot_longer(names_to = "time", values_to = "sigma", cols = t1:t4) |>
    ggplot(aes(sigma, analyte)) +
    geom_point(size = 3, shape = 21, aes(fill = time)) +
    alpha_theme() +
    theme(panel.grid.major.y = element_line(color = "gray80")) +
    labs(
      x = "Residual sigma",
      y = "Analyte",
      subtitle = "Residuals sigma at each time-point",
      fill = "Time"
    ) +
    scale_fill_manual(values = colscale)

  p3 <- bind_rows(
    resids0 |>
      correlate(quiet = TRUE) |>
      pivot_longer(
        names_to = "xterm",
        values_to = "c",
        cols = -one_of("term")
      ) |>
      mutate(time = "T1"),
    resids1 |>
      correlate(quiet = TRUE) |>
      pivot_longer(
        names_to = "xterm",
        values_to = "c",
        cols = -one_of("term")
      ) |>
      mutate(time = "T2"),
    resids2 |>
      correlate(quiet = TRUE) |>
      pivot_longer(
        names_to = "xterm",
        values_to = "c",
        cols = -one_of("term")
      ) |>
      mutate(time = "T4")
  ) |>
    mutate(cor.label = if_else(is.na(c), "", as.character(round(c, 2)))) |>

    ggplot(aes(xterm, term)) +

    geom_tile(aes(fill = c)) +
    geom_text(aes(label = cor.label), size = 3) +
    scale_fill_gradient2(low = colscale[1], mid = "white", high = colscale[2]) +

    facet_wrap(~time, ncol = 3) +
    alpha_theme() +
    theme(axis.title = element_blank()) +
    labs(
      fill = "Correlation (residual)",
      subtitle = "Correlation (residuals) between analytes"
    )

  return(list(p1 = p1, p2 = p2, p3 = p3))
}
