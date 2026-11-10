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
# ad, adapt delta
fit_wrap <- function(
  model,
  DATA = d,
  fit = "fit.m1",
  iter = 2000,
  iterwu = 1000,
  chains = 4,
  pchains = 4,
  refresh = 1000,
  ad = 0.8,
  init = NULL,
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
      refresh = refresh,
      adapt_delta = ad,
      init = init # NULL -> cmdstan random init; or list/function/number/path
    ) # Number of iterations before each print
    # Save the model for easy access
    fittedmodel$save_object(file = fit.filepath) # <-- use this
  } else {
    fittedmodel <- readRDS(fit.filepath)
  }
  # Assign the model object to the environment.
  assign(fit, fittedmodel, envir = parent.frame())
}


# Generic initial values for models s7-s20
#
# cmdstan's default random inits (uniform(-2, 2) on the unconstrained scale)
# can start the 2K-dimensional cholesky_factor_corr near-singular, freezing the
# chain (step size -> 0, all draws identical). This returns an init function
# that starts in a safe region:
#   * L_Omega* (cholesky_factor_corr): identity
#   * positive scales (sigma*, psi, lambda_diag, tau_*): uniform(lo, hi)
#   * everything else (alpha, beta, imputed values, loadings, *_raw):
#     normal(0, jitter)
# Parameter names and sizes are read from a one-iteration fixed_param run, so
# sizes defined in transformed data (n_miss_y, n_lambda_lower, ...) are exact.
#
# model, a compiled cmdstanr model
# DATA, the data list used to fit the model
# jitter, sd of the normal jitter for unconstrained parameters
# lo, hi, range for positive scale parameters
# returns a function(chain_id) for the init argument of fit_wrap/$sample
make_inits <- function(model, DATA = d, jitter = 0.3, lo = 0.5, hi = 1.5) {
  skel <- model$sample(
    data = DATA,
    fixed_param = TRUE,
    chains = 1,
    iter_sampling = 1,
    init = 0, # identity correlations, unit scales, zeros elsewhere
    refresh = 0,
    show_messages = FALSE
  )
  pars <- model$variables()$parameters
  sizes <- skel$metadata()$stan_variable_sizes[names(pars)]

  function(chain_id = 1) {
    out <- lapply(names(pars), function(nm) {
      sz <- sizes[[nm]]
      if (is.null(sz)) sz <- 0 # zero-length vectors are absent from metadata
      n <- prod(sz)
      if (grepl("^L_Omega", nm)) {
        diag(sz[1])
      } else {
        x <- if (grepl("^(sigma|psi|lambda_diag|tau_)", nm)) {
          runif(n, lo, hi)
        } else {
          rnorm(n, 0, jitter)
        }
        if (pars[[nm]]$dimensions == 0) x else array(x, dim = sz)
      }
    })
    names(out) <- names(pars)
    out
  }
}


# Check for effect of divergent transitions in parameters
#
# fit, a cmdstanr fit (e.g., from fit_wrap)
# vars, a string indicating which parameter names should be used. The string
#  is used in a grep call. Use e.g., "^(psi|lambda_lower|lambda_diag|
#  sigma|sigma_pre|mu_pre)".
# nparam, how many parameters of top affected parameters?
# returns a vector of standardized differences between divergenet/non-divergent
#  transistions
ceck_div <- function(fit, vars, nparam = 15) {
  div <- fit$sampler_diagnostics(format = "df")$divergent__ # aligned row-wise with dd

  if (all(div == 0)) {
    message("No divergent transitions detected")
  } else {
    draws <- fit$draws(format = "df")
    vars <- grep(vars, colnames(dd), value = TRUE)
    scan <- sapply(vars, function(v) {
      x <- dd[[v]]
      (mean(x[div == 1]) - mean(x[div == 0])) / sd(x, na.rm = TRUE) # std. mean diff
    })
    sort(abs(scan), decreasing = TRUE)[1:nparam]
  }
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
# extra_data, a named list added to each Stan data list, e.g. list(nu = 4)
#   for the Student-t models (s12, s13, s18).
# use_inits, if TRUE, chains start from make_inits() (built per data set)
#   instead of cmdstan's random inits.
#
# Each model is fit on both the full and missing variants so that
# marginal group contrasts can be compared with vs. without missingness.
# Output structure matches comp_fun (coverage, error_comp, est_error_dist,
# average_hmc_metrics) so downstream plotting code works unchanged.
comp_fun2 <- function(d, mods, extra_data = list(), use_inits = FALSE) {
  # Parameter names present in a fit (zero-length parameters, e.g. y_impute
  # in complete data, are absent from the output)
  fit_pars <- function(fit, mod) {
    intersect(names(mod$variables()$parameters), fit$metadata()$stan_variables)
  }

  K <- ncol(d$true_effects)
  # True treatment differences per time x type
  true_effects <- data.frame(
    time = as.character(c(rep(1, K), rep(2, K))),
    type = as.character(rep(1:K, 2)),
    true = c(d$true_effects[1, ], d$true_effects[2, ])
  )

  fit_model <- function(mod, dat) {
    dat$prior_only <- 0
    dat[names(extra_data)] <- extra_data
    fit <- mod$sample(
      data = dat,
      chains = 3,
      iter_warmup = 500,
      iter_sampling = 750,
      parallel_chains = 1,
      refresh = 0,
      show_messages = FALSE,
      init = if (use_inits) make_inits(mod, DATA = dat) else NULL
    )
    # Read parameters only: generated quantities (y_rep, log_lik, ...) are
    # several hundred MB per fit and are not used here.
    fit$draws(variables = fit_pars(fit, mod))
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
    fit$summary(variables = fit_pars(fit, mods[[model]])) |>
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
  # Check distinct to help with debug
  if (length(unique(d$analyte)) != 12) stop("The function expect 12 analytes")

  compl_temp <- d |>
    filter(!is.na(value_std)) |>
    mutate(.by = participant, n = n()) |>
    # This is hardcoded and assumes 11 measurements and 3 time points
    filter(n == 36) |>
    select(-value) |>
    pivot_wider(names_from = time, values_from = value_std)

  #
  resids0 <- list()
  resids1 <- list()
  resids2 <- list()
  resids_time <- list()
  sigmas <- list()
  analytes <- unique(compl_temp$analyte)
  for (i in seq_along(analytes)) {
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
    # Tiles are too small for printed values at text width; fill shows the
    # correlation
    scale_fill_gradient2(low = colscale[1], mid = "white", high = colscale[2]) +

    facet_wrap(~time, ncol = 3) +
    alpha_theme() +
    theme(
      # alpha_theme() sets axis.title.x/y, which override axis.title
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 6),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    labs(
      fill = "Correlation (residual)",
      subtitle = "Correlation (residuals) between analytes"
    )

  return(list(p1 = p1, p2 = p2, p3 = p3))
}


# Gurka-DeBoer MetS severity score ###########################################
# Adult (20-64 y), waist-based, sex-specific non-Hispanic White equations from
# Gurka et al. (2014) Metabolism 63(2):218-25:
#   Z = intercept + b1 * waist + b2 * HDL + b3 * SBP + b4 * ln(TG) + b5 * glucose
# with waist in cm, SBP in mmHg and HDL, TG and glucose in mg/dL.
# NOTE: the slopes match Table 1 of DeBoer & Gurka (2017) Diabetes Metab Syndr
# Obes 10:65-72, but that table lists intercepts -5.4473 (m) and -7.2536 (f).
# The intercepts are now fixed to exactly match DeBoer/Gurka 2017.

# The intercept shifts the level of the score only, not treatment contrasts.
gurka_coef <- data.frame(
  sex = c("m", "f"),
  intercept = c(-5.4473, -7.2536),
  waist = c(0.0125, 0.0254),
  HDL = c(-0.0251, -0.0120),
  sys = c(0.0047, 0.0075),
  TG = c(0.8244, 0.5800),
  glucose = c(0.0106, 0.0203)
)

# Unit conversion from data units (mmol/L) to score units (mg/dL), see
# R/data-prep.R. TG enters the score as ln(TG), where the conversion becomes
# an additive offset log(88.57).
gurka_conv <- c(
  waist = 1,
  HDL = 38.6654,
  sys = 1,
  TG = 88.57,
  glucose = 18.0156
)


# Back-transform standardised values of one analyte to score units.
#
# z, values (vector or draws x participants matrix) on the standardised scale
# sc, one row of dat$scale (data_load_s7s10) for the analyte
# returns values in score units (ln(mg/dL) for TG)
gurka_natural <- function(z, sc) {
  x <- sc$mean + sc$sd * z
  if (sc$transform == "log") x + log(gurka_conv[[sc$analyte]]) else
    x * gurka_conv[[sc$analyte]]
}


# Per-participant Gurka-DeBoer score per posterior draw
#
# Uses the (observed or imputed) values in pre_wide (T1), y_wide_1 (T2) and
# y_wide_2 (T4) from models s7-s13, so all participants get a posterior score
# at all time-points. Scores from observed values are constant over draws.
#
# fit, a cmdstanr fit (e.g., fit.s13.post)
# dat, output of data_load_s7s10()
# coef, coefficient table (default gurka_coef)
# returns a long data frame with .draw, P, sex, treat, time and score
gurka_score <- function(fit, dat, coef = gurka_coef) {
  an <- names(gurka_conv)
  sc <- dat$scale[match(an, dat$scale$analyte), ]

  part <- dat$dat |>
    distinct(P, sex, treat) |>
    arrange(P) |>
    mutate(
      sex = if_else(sex < 0, "m", "f"),
      treat = dat$treat$treat[match(treat, dat$treat$treat_idx)]
    )
  # Participant-specific coefficients (participants x terms)
  b <- coef[match(part$sex, coef$sex), ]

  vars <- c(pre_wide = "T1", y_wide_1 = "T2", y_wide_2 = "T4")
  dr <- posterior::as_draws_matrix(fit$draws(names(vars)))

  out <- list()
  for (v in names(vars)) {
    # draws x participants; the intercept is added per participant column
    score <- matrix(b$intercept, nrow(dr), nrow(part), byrow = TRUE)
    for (j in seq_along(an)) {
      z <- unclass(dr[, sprintf("%s[%d,%d]", v, part$P, sc$type[j])])
      score <- score + sweep(gurka_natural(z, sc[j, ]), 2, b[[an[j]]], "*")
    }
    out[[v]] <- data.frame(
      .draw = rep(seq_len(nrow(dr)), times = nrow(part)),
      P = rep(part$P, each = nrow(dr)),
      sex = rep(part$sex, each = nrow(dr)),
      treat = rep(part$treat, each = nrow(dr)),
      time = vars[[v]],
      score = as.vector(score)
    )
  }
  bind_rows(out)
}


# Treatment contrast (krill oil - placebo) on the Gurka-DeBoer score
#
# The score is linear in waist, HDL, SBP, ln(TG) and glucose, and each of these
# is a linear back-transform of the standardised model scale. The expected
# difference in score between treatments (same baseline values and sex) is
# therefore exactly
#   sum_k b_k * sd_k * c_k * (alpha[t, k, 2] - alpha[t, k, 1])
# where c_k is the unit conversion (1 for ln(TG)). Coefficients are
# sex-specific, so the contrast is returned per sex and as a sample-weighted
# average over sex.
#
# fit, a cmdstanr fit (models s7-s13, alpha[time, type, treat])
# dat, output of data_load_s7s10()
# coef, coefficient table (default gurka_coef)
# by_component, if TRUE the contrast is returned per score component (one term
#   per analyte, summing exactly to the total contrast)
# returns a long data frame with .draw, time, sex (and component) and contrast
gurka_contrast <- function(fit, dat, coef = gurka_coef, by_component = FALSE) {
  an <- names(gurka_conv)
  sc <- dat$scale[match(an, dat$scale$analyte), ]
  a <- posterior::as_draws_matrix(fit$draws("alpha"))

  # Proportion of women in the sample, used for the average over sex
  p_f <- dat$dat |> distinct(P, sex) |> summarise(p = mean(sex > 0)) |> pull(p)

  out <- list()
  for (t in 1:2) {
    d <- sapply(
      sc$type,
      function(k)
        a[, sprintf("alpha[%d,%d,2]", t, k)] -
          a[, sprintf("alpha[%d,%d,1]", t, k)]
    )
    # Contribution of each component (draws x components) per sex
    eff <- lapply(c(m = "m", f = "f"), function(s) {
      w <- unlist(coef[coef$sex == s, an]) *
        sc$sd *
        if_else(sc$transform == "log", 1, gurka_conv[an])
      sweep(d, 2, w, "*")
    })
    for (j in if (by_component) seq_along(an) else 0) {
      pick <- function(x) if (j == 0) rowSums(x) else x[, j]
      out[[paste(t, j)]] <- data.frame(
        .draw = seq_len(nrow(a)),
        time = c("T2", "T4")[t],
        component = if (j == 0) "total" else an[j],
        m = pick(eff$m),
        f = pick(eff$f),
        average = (1 - p_f) * pick(eff$m) + p_f * pick(eff$f)
      )
    }
  }
  res <- bind_rows(out) |>
    pivot_longer(c(m, f, average), names_to = "sex", values_to = "contrast")
  if (!by_component) res <- select(res, -component)
  res
}


# Continuous MetS score from external reference values #######################
# Sum of (signed, weighted) z-scores of the MetS criteria relative to sex- and
# age-specific reference means and SDs, as in Guseman et al. (2018) Acad
# Pediatr 18(5):589-92. The reference values are supplied by the user, e.g.
# from a population survey (HUNT4, SCAPIS, ...).
#
# ref, a data frame with columns analyte, mean and sd, and optionally sex
#   ("m"/"f") and age_lwr, age_upr (inclusive age band). Each participant must
#   match exactly one row per analyte. Values are in data units: waist (cm),
#   sys/dias/map (mmHg), HDL and glucose (mmol/L), and TG as mean and SD of
#   ln(TG, mmol/L) (for ln(mg/dL) subtract log(88.57) from the mean).
# w, named weights incl. sign. Analytes: waist, HDL, TG, glucose, sys, dias and
#   map (mean arterial pressure, dias + (sys - dias) / 3).

# Reference rows (one per participant) for one analyte
#
# part, participant lookup (dat$part from data_load_s7s10())
# r, the rows of ref for one analyte
cmets_match_ref <- function(part, r) {
  if (is.null(part))
    stop(
      "dat$part is missing: source R/data-prep.R and re-run data_load_s7s10()"
    )
  if (nrow(r) == 0) stop("No reference rows for an analyte in w")
  i <- sapply(seq_len(nrow(part)), function(p) {
    ok <- rep(TRUE, nrow(r))
    if ("sex" %in% names(r)) ok <- ok & r$sex == part$sex[p]
    if ("age_lwr" %in% names(r))
      ok <- ok & part$age[p] >= r$age_lwr & part$age[p] <= r$age_upr
    if (sum(ok) != 1)
      stop(sprintf(
        "Participant P = %d matches %d reference rows for %s",
        part$P[p],
        sum(ok),
        r$analyte[1]
      ))
    which(ok)
  })
  r[i, ]
}


# Natural-scale value (data units) of an analyte, incl. derived map
#
# get_z, function(analyte) returning standardised values (vector or matrix)
# scale, dat$scale
cmets_natural <- function(a, get_z, scale) {
  if (a == "map") {
    dias <- cmets_natural("dias", get_z, scale)
    return(dias + (cmets_natural("sys", get_z, scale) - dias) / 3)
  }
  sc <- scale[scale$analyte == a, ]
  sc$mean + sc$sd * get_z(a)
}


# Per-participant continuous MetS score per posterior draw
#
# Uses the (observed or imputed) values in pre_wide (T1), y_wide_1 (T2) and
# y_wide_2 (T4) from models s7-s13.
#
# fit, a cmdstanr fit (e.g., fit.s13.post)
# dat, output of data_load_s7s10()
# ref, reference values (see above)
# w, named weights incl. sign
# returns a long data frame with .draw, P, sex, treat, time and score
cmets_score <- function(
  fit,
  dat,
  ref,
  w = c(waist = 1, HDL = -1, TG = 1, glucose = 1, map = 1)
) {
  part <- dat$part
  type_idx <- setNames(dat$scale$type, dat$scale$analyte)

  vars <- c(pre_wide = "T1", y_wide_1 = "T2", y_wide_2 = "T4")
  dr <- posterior::as_draws_matrix(fit$draws(names(vars)))

  out <- list()
  for (v in names(vars)) {
    # draws x participants
    get_z <- function(a)
      unclass(dr[, sprintf("%s[%d,%d]", v, part$P, type_idx[[a]])])
    score <- 0
    for (a in names(w)) {
      r <- cmets_match_ref(part, ref[ref$analyte == a, ])
      x <- cmets_natural(a, get_z, dat$scale)
      z <- sweep(sweep(x, 2, r$mean, "-"), 2, r$sd, "/")
      score <- score + w[[a]] * z
    }
    out[[v]] <- data.frame(
      .draw = rep(seq_len(nrow(dr)), times = nrow(part)),
      P = rep(part$P, each = nrow(dr)),
      sex = rep(part$sex, each = nrow(dr)),
      treat = rep(part$treat, each = nrow(dr)),
      time = vars[[v]],
      score = as.vector(score)
    )
  }
  bind_rows(out)
}


# Treatment contrast (krill oil - placebo) on the continuous MetS score
#
# The score is linear in the model scale, so the expected contrast for
# participant p (same baseline values and sex) is
#   sum_a w_a * SD_a^marg * (alpha[t, a, 2] - alpha[t, a, 1]) / SD_a,p^ref
# (map: 2/3 of the dias and 1/3 of the sys contrast). Reference SDs depend on
# sex and age, so contrasts are averaged over the participants in the sample,
# per sex and overall. Reference means cancel.
#
# fit, a cmdstanr fit (models s7-s13, alpha[time, type, treat])
# dat, output of data_load_s7s10()
# ref, reference values (see above)
# w, named weights incl. sign
# by_component, if TRUE the contrast is returned per score component (one term
#   per element of w, summing exactly to the total contrast)
# returns a long data frame with .draw, time, sex (and component) and contrast
cmets_contrast <- function(
  fit,
  dat,
  ref,
  w = c(waist = 1, HDL = -1, TG = 1, glucose = 1, map = 1),
  by_component = FALSE
) {
  part <- dat$part
  type_idx <- setNames(dat$scale$type, dat$scale$analyte)
  a_dr <- posterior::as_draws_matrix(fit$draws("alpha"))

  out <- list()
  for (t in 1:2) {
    # Contrast on the standardised scale (draws x participants, constant over
    # participants); scale$mean cancels in the natural-scale difference
    get_z <- function(a) {
      d <- unclass(
        a_dr[, sprintf("alpha[%d,%d,2]", t, type_idx[[a]])] -
          a_dr[, sprintf("alpha[%d,%d,1]", t, type_idx[[a]])]
      )
      matrix(d, nrow(a_dr), nrow(part))
    }
    scale0 <- mutate(dat$scale, mean = 0)
    # Contribution of each component (draws x participants)
    eff <- lapply(names(w), function(a) {
      r <- cmets_match_ref(part, ref[ref$analyte == a, ])
      w[[a]] * sweep(cmets_natural(a, get_z, scale0), 2, r$sd, "/")
    })
    names(eff) <- names(w)
    if (!by_component) eff <- list(total = Reduce(`+`, eff))
    for (k in names(eff)) {
      out[[paste(t, k)]] <- data.frame(
        .draw = seq_len(nrow(a_dr)),
        time = c("T2", "T4")[t],
        component = k,
        m = rowMeans(eff[[k]][, part$sex == "m", drop = FALSE]),
        f = rowMeans(eff[[k]][, part$sex == "f", drop = FALSE]),
        average = rowMeans(eff[[k]])
      )
    }
  }
  res <- bind_rows(out) |>
    pivot_longer(c(m, f, average), names_to = "sex", values_to = "contrast")
  if (!by_component) res <- select(res, -component)
  res
}


# Per-participant values of single analytes per posterior draw
#
# Observed or imputed values in pre_wide (T1), y_wide_1 (T2) and y_wide_2
# (T4) from models s7-s13, back-transformed to the natural scale (exp() for
# log-transformed analytes, e.g. HOMA2IR).
#
# fit, a cmdstanr fit (e.g., fit.s13.post)
# dat, output of data_load_s7s10()
# analytes, analyte names as in dat$scale
# returns a long data frame with .draw, P, time, analyte, value and observed
#   (1 if the value was observed, 0 if imputed)
wide_values <- function(fit, dat, analytes) {
  part <- dat$part
  sc <- dat$scale[match(analytes, dat$scale$analyte), ]

  vars <- c(pre_wide = "T1", y_wide_1 = "T2", y_wide_2 = "T4")
  dr <- posterior::as_draws_matrix(fit$draws(names(vars)))

  # Observed indicators per participant, time and type
  obs <- bind_rows(
    dat$dat |> distinct(P, type, observed = pre_miss) |> mutate(time = "T1"),
    dat$dat |>
      transmute(P, type, observed = y_miss, time = c("T2", "T4")[time])
  )

  out <- list()
  for (v in names(vars)) {
    for (j in seq_len(nrow(sc))) {
      z <- unclass(dr[, sprintf("%s[%d,%d]", v, part$P, sc$type[j])])
      x <- sc$mean[j] + sc$sd[j] * z
      if (sc$transform[j] == "log") x <- exp(x)
      out[[paste(v, j)]] <- data.frame(
        .draw = rep(seq_len(nrow(dr)), times = nrow(part)),
        P = rep(part$P, each = nrow(dr)),
        time = vars[[v]],
        analyte = sc$analyte[j],
        type = sc$type[j],
        value = as.vector(x)
      )
    }
  }
  bind_rows(out) |>
    left_join(obs, by = c("P", "type", "time")) |>
    select(-type)
}
