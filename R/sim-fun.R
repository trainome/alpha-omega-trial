# Avoid testing code in the bottom when sourcing
test <- FALSE


# Simulate type-specific measurements conditionally on pre-intervention
# values and co-variates for K number of measurements types from
# multivariate normal distribution(s).
#
# n,    the number of observations (participants)
# K,    number of measurement types (k)
# form, linear combination of predictor variables for each k.
#   If "pre" is included, it will be expanded to pre1...preK internally,
#   with each type's pre-intervention value as a separate predictor.
# betas, a list of length 2, one p * K matrix per post-intervention
#   time-point. Parameters must contain the intercept as the first row.
# preds, predictor variables with nrow == n. Should not include pre
#   columns as these are simulated internally.
# C1,   pre-intervention correlations (upper triangular vector of
#   length (K*K - K) / 2). The triangular matrix fills by column.
# C2,   a list of length 2, one vector of post-intervention residual
#   correlations per time-point (upper triangular, length (K*K-K)/2)
# S1,   a vector of length K, pre-intervention SDs per type
# S2,   a list of length 2, one vector of length K of post-intervention
#   residual SDs per timepoint
# C_cross, a vector of length K, between-time-point residual correlations
#   within each type
# verbose, logical, will print input values and parameters to the console.
# -----------------------------------------------------------------
# Missing values. Some participants are missing in full, some have missing at
# random in a specific type.
# missing_form_t2,  formula specifying the predictor for missing values at t2
# missing_form_t4,  formula for t4 (if missing at t2 then also miss at t4)
# missing_betas_t2,    parameter values for missing values simulation
# missing_betas_t4,    parameter values for missing values simulation
# missing_y_prop, a proportion of rows (y) will be missing at random
# missing_pre_prop, a proportion of rows in pre will be missing at random.
sim_type <- function(
  n,
  K,
  form_pre,
  betas_pre,
  form,
  betas,
  preds,
  C1,
  C2,
  S1,
  S2,
  C_cross,
  missing_form_t2 = NULL,
  missing_form_t4 = NULL,
  missing_beta_t2,
  missing_beta_t4,
  missing_y_prop = 0,
  missing_pre_prop = 0,
  verbose = TRUE
) {
  # Input validation --------------------------------------------------------

  if (nrow(preds) != n) stop("The n is not equal to nrow(preds)")

  expected_len <- (K * K - K) / 2

  if (length(C1) != expected_len) stop("Length of C1 must be (K*K - K) / 2")
  if (length(S1) != K) stop("S1 must be of length K")

  if (!is.list(C2) || length(C2) != 2) stop("C2 must be a list of length 2")
  if (!is.list(S2) || length(S2) != 2) stop("S2 must be a list of length 2")
  if (!is.list(betas) || length(betas) != 2)
    stop("betas must be a list of length 2")

  for (t in 1:2) {
    if (length(C2[[t]]) != expected_len)
      stop(sprintf("C2[[%d]] must be of length (K*K - K) / 2", t))
    if (length(S2[[t]]) != K) stop(sprintf("S2[[%d]] must be of length K", t))
    if (!is.matrix(betas[[t]])) stop(sprintf("betas[[%d]] must be a matrix", t))
  }

  if (length(C_cross) != K) stop("C_cross must be of length K")
  if (any(abs(C_cross) > 1)) stop("C_cross values must be in [-1, 1]")

  # Pre-intervention covariance and draws -----------------------------------
  build_sigma <- function(C, S) {
    R <- matrix(1, K, K)
    R[upper.tri(R)] <- C # This fills by column
    R[lower.tri(R)] <- t(R)[lower.tri(R)]
    diag(S) %*% R %*% diag(S)
  }

  Sigma1 <- build_sigma(C1, S1)
  Z_pre <- MASS::mvrnorm(n, rep(0, K), Sigma1) # n x K

  desig_pre <- model.matrix(form_pre, data = preds)
  mu_pre <- desig_pre %*% betas_pre

  pre <- Z_pre + mu_pre

  colnames(pre) <- paste0("type", 1:K)

  # Design matrix ---------------------------------------------------------
  # If "pre" is in the formula, expand to pre1...preK per outcome
  form_terms <- attr(terms(form), "term.labels")
  has_pre <- "pre" %in% form_terms # True if pre is in form

  if (has_pre) {
    pre_df <- as.data.frame(pre)
    colnames(pre_df) <- paste0("pre", 1:K)
    preds_expanded <- cbind(preds, pre_df)

    # One design matrix per outcome with matching preK column
    design <- lapply(1:K, function(k) {
      form_k <- as.formula(gsub("\\bpre\\b", paste0("pre", k), deparse(form)))
      model.matrix(form_k, data = preds_expanded)
    })
  } else {
    # Shared design matrix across outcomes, replicated into list for consistency
    design <- replicate(K, model.matrix(form, data = preds), simplify = FALSE)
  }

  # Number of parameters
  p <- ncol(design[[1]])

  # Validate betas dimensions now that p is known
  for (t in 1:2) {
    if (!all(dim(betas[[t]]) == c(p, K)))
      stop(sprintf(
        "betas[[%d]] must be %d x %d (predictors x K), got %d x %d",
        t,
        p,
        K,
        nrow(betas[[t]]),
        ncol(betas[[t]])
      ))
  }

  # Conditional means: n x K per timepoint -----------------------------------

  # Design is now a list and betas a list with matrixes. This selects the
  # relevant list and row of the betas matrix.

  mu_t1 <- sapply(1:K, function(k) design[[k]] %*% betas[[1]][, k])
  mu_t2 <- sapply(1:K, function(k) design[[k]] %*% betas[[2]][, k])
  mu_mat <- cbind(mu_t1, mu_t2) # n x 2K

  # Post-intervention block covariance matrix --------------------------------
  Sigma2_t1 <- build_sigma(C2[[1]], S2[[1]])
  Sigma2_t2 <- build_sigma(C2[[2]], S2[[2]])

  # Cross-time-point covariance: same type at t1 and t2 only
  Sigma_cross <- diag(S2[[1]] * C_cross * S2[[2]])

  Sigma_full <- rbind(
    cbind(Sigma2_t1, Sigma_cross),
    cbind(t(Sigma_cross), Sigma2_t2)
  )

  # Check positive definiteness
  eigs <- eigen(Sigma_full, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigs <= 0))
    stop("Sigma_full is not positive definite - check correlation inputs")

  # Draw residuals and add conditional means ---------------------------------
  Z <- MASS::mvrnorm(n, rep(0, 2 * K), Sigma_full) # Residuals
  Y <- Z + mu_mat # n x 2K conditional means

  # Extract time-point-wise post-intervention values
  y_t1 <- Y[, 1:K]
  y_t2 <- Y[, (K + 1):(2 * K)]

  # Set type names on columns
  type_names <- paste0("type", 1:K)
  colnames(y_t1) <- type_names
  colnames(y_t2) <- type_names

  # If verbose - Validate input in console ----------------------------------
  if (verbose) {
    # Display parameters ---
    param_names <- colnames(design[[1]])
    if (has_pre) {
      # Replace pre1/pre2/... with generic "preK" for display clarity
      param_names <- gsub("pre\\d+", "preK", param_names)
    }
    for (t in 1:2) {
      b <- betas[[t]]
      rownames(b) <- param_names
      colnames(b) <- paste("Type", 1:K)
      cat(sprintf("Betas for timepoint %d:\n", t))
      print(b)
      cat("\n")
    }
  }

  # Missing-value indicators (applied to a copy; originals kept for the
  # full/no-missing dataset).
  miss_t2_ind <- rep(0L, n)
  miss_t4_ind <- rep(0L, n)

  if (!is.null(missing_form_t2)) {
    X_miss_t2 <- model.matrix(missing_form_t2, data = preds)
    p_miss_t2 <- plogis(X_miss_t2 %*% missing_beta_t2)
    miss_t2_ind <- rbinom(n, 1, p_miss_t2)
  }

  if (!is.null(missing_form_t4)) {
    X_miss_t4 <- model.matrix(missing_form_t4, data = preds)
    p_miss_t4 <- plogis(X_miss_t4 %*% missing_beta_t4)
    miss_t4_ind <- rbinom(n, 1, p_miss_t4)
    # Monotone: missing at t2 implies missing at t4
    miss_t4_ind[miss_t2_ind == 1] <- 1
  }

  y_t1_miss <- y_t1
  y_t2_miss <- y_t2
  y_t1_miss[miss_t2_ind == 1, ] <- NA
  y_t2_miss[miss_t4_ind == 1, ] <- NA

  # Collect all data in a data frame
  preds$participant <- 1:n

  build_long <- function(y_t1_m, y_t2_m, pre_m) {
    pre_df <- data.frame(pre_m) |>
      dplyr::mutate(participant = 1:n) |>
      tidyr::pivot_longer(
        cols = tidyselect::starts_with("type"),
        names_to = "type",
        values_to = "pre"
      )

    y_df <- dplyr::bind_rows(
      data.frame(y_t1_m) |>
        dplyr::mutate(participant = 1:n, time = 1) |>
        tidyr::pivot_longer(
          cols = tidyselect::starts_with("type"),
          names_to = "type",
          values_to = "y"
        ),
      data.frame(y_t2_m) |>
        dplyr::mutate(participant = 1:n, time = 2) |>
        tidyr::pivot_longer(
          cols = tidyselect::starts_with("type"),
          names_to = "type",
          values_to = "y"
        )
    )

    y_df |>
      inner_join(pre_df, by = c("participant", "type")) |>
      inner_join(preds, by = "participant")
  }

  full <- build_long(y_t1, y_t2, pre)
  miss <- build_long(y_t1_miss, y_t2_miss, pre)

  # Row-wise random missingness — applied to the miss copy only
  if (missing_y_prop > 0) {
    idx_y <- sample(seq_len(nrow(miss)), round(nrow(miss) * missing_y_prop, 0))
    miss$y[idx_y] <- NA
  }
  if (missing_pre_prop > 0) {
    idx_pre <- sample(
      seq_len(nrow(miss)),
      round(nrow(miss) * missing_pre_prop, 0)
    )
    miss$pre[idx_pre] <- NA
  }

  return(list(full = full, miss = miss))
}


# Testing the function ---------------------------------------------------

# Do not run when sourced
if (test) {
  N <- 10 # Number of participants
  K <- 4 # number of measurement types

  # The formula for the pre-intervention variables
  form_pre <- ~ 1 + age + sex

  # beta_pre is a p * K matrix of parameter (p in 1 to P)
  # values for each type (k in 1:K).
  # fmt: skip
  betas_pre <- matrix(
  c(0,    0,    0,    0, # intercepts
    0.2,    0.3,    0.4,    0.5, # age
    0.1,    0.1,    0.1,    0.1  # sex
  ),
  nrow = 3,
  ncol = K,
  byrow = TRUE
)

  # The same formula for both time-points
  form_main <- ~ 1 + pre + group + sex

  # betas matrix are p * K
  # fmt: skip
  betas_t1 <- matrix(
  c(0,        0,    0,        0, # intercepts
    0.9,    0.8,    0.7,    0.8, # pre_k effect 
    0.2,    0.3,    0.4,    0.5, # group
    0.1,    0.1,    0.1,    0.1  # sex
  ),
  nrow = 4,
  ncol = K,
  byrow = TRUE
)

  # fmt: skip
  betas_t2 <- matrix(
  c(  0,      0,      0,    0,   # intercepts
    0.9,    0.8,    0.7,    0.8, # preK effect (each column is that type's own pre)
    0.2,    0.3,    0.4,    0.5, # group
    0.1,    0.1,    0.1,    0.1  # sex
  ),
  nrow = 4,
  ncol = K,
  byrow = TRUE
)

  # Predictor variables, used in *all* sims
  preds_all <- data.frame(
    group = rep(c(0, 1), each = N / 2),
    sex = rep(c(0, 1), each = N / 2),
    age = rnorm(N)
  )

  # Vector of correlations (upper triangle, fill by column) for the pre int
  # values.
  c1 <- c(
    0.4, # r12
    0, # r13
    0.1, # r23
    -0.2, # r14
    0, # r24
    0 # r34
  )

  # A list of two (per time-point) residual correlations
  # (upper triangle, fill by row)
  c2 <- list(
    c(
      0.4, # r12
      0, # r13
      0.1, # r23
      -0.2, # r14
      0, # r24
      0 # r34
    ),
    c(
      0.4, # r12
      0, # r13
      0.1, # r23
      -0.2, # r14
      0, # r24
      0 # r34
    )
  )

  # Pre intervention SD's per type (length K)
  s1 <- c(1, 1, 1, 1)
  # Residual SD per type in each time point
  s2 <- list(
    c(0.1, 0.2, 0.3, 0.4),
    c(0.1, 0.2, 0.3, 0.4)
  )
  # Between time-point (residual) correlations within each type (lenght K)
  c_cross <- c(0.2, 0.2, 0.2, 0.2)

  #
  # missing_form_t2,  formula specifying the predictor for missing values at t2
  miss_t2 <- ~ 1 + age
  miss_t4 <- ~ 1 + age

  # Parameter values for the form complete missingness
  miss_beta_t2 <- c(-2.5, -0.5)
  miss_beta_t4 <- c(-2, -0.5)

  # Random missing proportion
  miss_prop_pre <- 0.02 # e.g., 2% missing of all rows
  miss_prop_y <- 0.05 # e.g., 5% miss of all rows post-int

  simout <- sim_type(
    n = N,
    K = K,

    form_pre = form_pre,
    betas_pre = betas_pre,

    form = form_main,
    preds = preds_all,
    betas = list(betas_t1, betas_t2),
    C1 = c1,
    S1 = s1,
    C2 = c2,
    S2 = s2,
    C_cross = c_cross,
    missing_form_t2 = miss_t2,
    missing_form_t4 = miss_t4,
    missing_beta_t2 = miss_beta_t2,
    missing_beta_t4 = miss_beta_t4,
    missing_pre_prop = miss_prop_pre,
    missing_y_prop = miss_prop_y,
    verbose = TRUE
  )
}


# Factor-based covariance builder -----------------------------------------
# Constructs a covariance matrix from a factor loading matrix, guaranteeing
# positive definiteness even with strong correlations.
#
# Lambda, a K x F matrix of factor loadings
# psi,    a K-vector of uniquenesses (residual variances, must be > 0)
# S,      a K-vector of marginal SDs per type
#
# Returns a K x K covariance matrix: diag(S) %*% R %*% diag(S)
# where R = cov2cor(Lambda %*% t(Lambda) + diag(psi)).
build_sigma_factor <- function(Lambda, psi, S) {
  K <- nrow(Lambda)
  R_raw <- Lambda %*% t(Lambda) + diag(psi)
  d <- sqrt(diag(R_raw))
  R <- R_raw / outer(d, d)
  diag(S) %*% R %*% diag(S)
}


# sim_type2 ----------------------------------------------------------------
# Same interface and output as sim_type, but covariance matrices are
# specified via factor loadings (Lambda/psi) instead of raw correlation
# vectors. This avoids non-positive-definite matrices when correlations
# are strong.
#
# n, K, form_pre, betas_pre, form, betas, preds — same as sim_type
#
# Lambda1,  K x F1 loading matrix for the pre-intervention covariance
# psi1,     K-vector of uniquenesses for the pre-intervention covariance
# S1,       K-vector of pre-intervention SDs
#
# Lambda2,  list of length 2, one K x F2 loading matrix per post-int
#           time-point
# psi2,     list of length 2, one K-vector of uniquenesses per time-point
# S2,       list of length 2, one K-vector of residual SDs per time-point
#
# C_cross,  K-vector of between-time-point residual correlations within
#           each type (same as sim_type)
#
# missing_* arguments — same as sim_type
# verbose — same as sim_type
sim_type2 <- function(
  n,
  K,
  form_pre,
  betas_pre,
  form,
  betas,
  preds,
  Lambda1,
  psi1,
  S1,
  Lambda2,
  psi2,
  S2,
  C_cross,
  missing_form_t2 = NULL,
  missing_form_t4 = NULL,
  missing_beta_t2,
  missing_beta_t4,
  missing_y_prop = 0,
  missing_pre_prop = 0,
  verbose = TRUE
) {
  # Input validation --------------------------------------------------------
  if (nrow(preds) != n) stop("n is not equal to nrow(preds)")
  if (nrow(Lambda1) != K) stop("Lambda1 must have K rows")
  if (length(psi1) != K) stop("psi1 must be of length K")
  if (any(psi1 <= 0)) stop("psi1 values must be > 0")
  if (length(S1) != K) stop("S1 must be of length K")

  if (!is.list(Lambda2) || length(Lambda2) != 2)
    stop("Lambda2 must be a list of length 2")
  if (!is.list(psi2) || length(psi2) != 2)
    stop("psi2 must be a list of length 2")
  if (!is.list(S2) || length(S2) != 2) stop("S2 must be a list of length 2")
  if (!is.list(betas) || length(betas) != 2)
    stop("betas must be a list of length 2")

  for (t in 1:2) {
    if (nrow(Lambda2[[t]]) != K)
      stop(sprintf("Lambda2[[%d]] must have K rows", t))
    if (any(psi2[[t]] <= 0)) stop(sprintf("psi2[[%d]] values must be > 0", t))
    if (length(S2[[t]]) != K) stop(sprintf("S2[[%d]] must be of length K", t))
    if (!is.matrix(betas[[t]])) stop(sprintf("betas[[%d]] must be a matrix", t))
  }

  if (length(C_cross) != K) stop("C_cross must be of length K")
  if (any(abs(C_cross) > 1)) stop("C_cross values must be in [-1, 1]")

  # Pre-intervention covariance and draws -----------------------------------
  Sigma1 <- build_sigma_factor(Lambda1, psi1, S1)
  Z_pre <- MASS::mvrnorm(n, rep(0, K), Sigma1)

  desig_pre <- model.matrix(form_pre, data = preds)
  mu_pre <- desig_pre %*% betas_pre
  pre <- Z_pre + mu_pre
  colnames(pre) <- paste0("type", 1:K)

  # Design matrix -----------------------------------------------------------
  form_terms <- attr(terms(form), "term.labels")
  has_pre <- "pre" %in% form_terms

  if (has_pre) {
    pre_df <- as.data.frame(pre)
    colnames(pre_df) <- paste0("pre", 1:K)
    preds_expanded <- cbind(preds, pre_df)

    design <- lapply(1:K, function(k) {
      form_k <- as.formula(gsub("\\bpre\\b", paste0("pre", k), deparse(form)))
      model.matrix(form_k, data = preds_expanded)
    })
  } else {
    design <- replicate(K, model.matrix(form, data = preds), simplify = FALSE)
  }

  p <- ncol(design[[1]])

  for (t in 1:2) {
    if (!all(dim(betas[[t]]) == c(p, K)))
      stop(sprintf(
        "betas[[%d]] must be %d x %d (predictors x K), got %d x %d",
        t,
        p,
        K,
        nrow(betas[[t]]),
        ncol(betas[[t]])
      ))
  }

  # Conditional means: n x K per timepoint -----------------------------------
  mu_t1 <- sapply(1:K, function(k) design[[k]] %*% betas[[1]][, k])
  mu_t2 <- sapply(1:K, function(k) design[[k]] %*% betas[[2]][, k])
  mu_mat <- cbind(mu_t1, mu_t2)

  # Post-intervention block covariance matrix --------------------------------
  Sigma2_t1 <- build_sigma_factor(Lambda2[[1]], psi2[[1]], S2[[1]])
  Sigma2_t2 <- build_sigma_factor(Lambda2[[2]], psi2[[2]], S2[[2]])

  # Cross-time-point covariance: same type at t1 and t2 only
  Sigma_cross <- diag(S2[[1]] * C_cross * S2[[2]])

  Sigma_full <- rbind(
    cbind(Sigma2_t1, Sigma_cross),
    cbind(t(Sigma_cross), Sigma2_t2)
  )

  eigs <- eigen(Sigma_full, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigs <= 0))
    stop(
      "Sigma_full is not positive definite - check C_cross relative to Lambda2/psi2"
    )

  # Draw residuals and add conditional means ---------------------------------
  Z <- MASS::mvrnorm(n, rep(0, 2 * K), Sigma_full)
  Y <- Z + mu_mat

  y_t1 <- Y[, 1:K]
  y_t2 <- Y[, (K + 1):(2 * K)]

  type_names <- paste0("type", 1:K)
  colnames(y_t1) <- type_names
  colnames(y_t2) <- type_names

  # Verbose output ----------------------------------------------------------
  if (verbose) {
    param_names <- colnames(design[[1]])
    if (has_pre) param_names <- gsub("pre\\d+", "preK", param_names)
    for (t in 1:2) {
      b <- betas[[t]]
      rownames(b) <- param_names
      colnames(b) <- paste("Type", 1:K)
      cat(sprintf("Betas for timepoint %d:\n", t))
      print(b)
      cat("\n")
    }
    cat("Implied pre-intervention correlation matrix:\n")
    print(round(cov2cor(Sigma1), 2))
    cat("\nImplied post-int residual correlation (t1):\n")
    print(round(cov2cor(Sigma2_t1), 2))
    cat("\nImplied post-int residual correlation (t2):\n")
    print(round(cov2cor(Sigma2_t2), 2))
  }

  # Missing values ----------------------------------------------------------
  miss_t2_ind <- rep(0L, n)
  miss_t4_ind <- rep(0L, n)

  if (!is.null(missing_form_t2)) {
    X_miss_t2 <- model.matrix(missing_form_t2, data = preds)
    p_miss_t2 <- plogis(X_miss_t2 %*% missing_beta_t2)
    miss_t2_ind <- rbinom(n, 1, p_miss_t2)
  }

  if (!is.null(missing_form_t4)) {
    X_miss_t4 <- model.matrix(missing_form_t4, data = preds)
    p_miss_t4 <- plogis(X_miss_t4 %*% missing_beta_t4)
    miss_t4_ind <- rbinom(n, 1, p_miss_t4)
    miss_t4_ind[miss_t2_ind == 1] <- 1
  }

  y_t1_miss <- y_t1
  y_t2_miss <- y_t2
  y_t1_miss[miss_t2_ind == 1, ] <- NA
  y_t2_miss[miss_t4_ind == 1, ] <- NA

  # Collect output -----------------------------------------------------------
  preds$participant <- 1:n

  build_long <- function(y_t1_m, y_t2_m, pre_m) {
    pre_df <- data.frame(pre_m) |>
      dplyr::mutate(participant = 1:n) |>
      tidyr::pivot_longer(
        cols = tidyselect::starts_with("type"),
        names_to = "type",
        values_to = "pre"
      )

    y_df <- dplyr::bind_rows(
      data.frame(y_t1_m) |>
        dplyr::mutate(participant = 1:n, time = 1) |>
        tidyr::pivot_longer(
          cols = tidyselect::starts_with("type"),
          names_to = "type",
          values_to = "y"
        ),
      data.frame(y_t2_m) |>
        dplyr::mutate(participant = 1:n, time = 2) |>
        tidyr::pivot_longer(
          cols = tidyselect::starts_with("type"),
          names_to = "type",
          values_to = "y"
        )
    )

    y_df |>
      inner_join(pre_df, by = c("participant", "type")) |>
      inner_join(preds, by = "participant")
  }

  full <- build_long(y_t1, y_t2, pre)
  miss <- build_long(y_t1_miss, y_t2_miss, pre)

  if (missing_y_prop > 0) {
    idx_y <- sample(seq_len(nrow(miss)), round(nrow(miss) * missing_y_prop, 0))
    miss$y[idx_y] <- NA
  }
  if (missing_pre_prop > 0) {
    idx_pre <- sample(
      seq_len(nrow(miss)),
      round(nrow(miss) * missing_pre_prop, 0)
    )
    miss$pre[idx_pre] <- NA
  }

  return(list(full = full, miss = miss))
}


# A convenience function for checking and structuring simulated data
#
# The data from the simulation function can be approximately standardized. To
# be safe, this function first performs standardization, similar to what is
# being done on the observed data. This will allow for use of the same priors
# (we expect the same e.g., intercept).
#
#
# d, a data set from the simulation function
# verbose, if TRUE prints the scaling procedure.
#
struct_sim <- function(d, verbose = TRUE) {
  # Standardize data
  dstd <- d |>
    mutate(
      group = dense_rank(group),
      sex = dense_rank(sex),
      type = dense_rank(type),
      participant = dense_rank(participant)
    ) |>
    mutate(time = paste0("time", time)) |>
    pivot_wider(names_from = time, values_from = y) |>
    pivot_longer(
      cols = c(pre, time1, time2),
      values_to = "y",
      names_to = "time"
    ) |>
    mutate(
      time = if_else(time == "pre", "T1", if_else(time == "time1", "T2", "T3"))
    ) |>
    mutate(.by = type, ystd = (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE))

  # Table the standardized and unstandardized data
  sum_tab <- dstd |>
    summarise(
      .by = type,
      m = mean(y, na.rm = TRUE),
      s = sd(y, na.rm = TRUE),
      mstd = mean(ystd, na.rm = TRUE),
      sstd = sd(ystd, na.rm = TRUE),
      n_miss = sum(is.na(y)),
      prop_miss = sum(is.na(y)) / n()
    )

  if (verbose) {
    cat("Synthetic data: Descriptive per type.\n")
    print(sum_tab)
  }

  # Organize data lists
  wide_complete <- dstd |>
    select(-y) |>
    # failsafe, average duplicates
    pivot_wider(names_from = time, values_from = ystd, values_fn = mean)

  # Complete case, single time-point
  tmp <- wide_complete |>
    select(participant:T1, y = T2) |>
    filter(if_all(everything(), complete.cases)) |>

    mutate(.by = participant, rows = n()) |>
    filter(rows == length(unique(wide_complete$type))) |>
    mutate(participant = dense_rank(participant))

  dl_compl_single <- list(
    N = nrow(tmp),
    N_P = length(unique(tmp$participant)),
    K = length(unique(tmp$type)),
    N_TX = length(unique(tmp$group)),

    sex = tmp$sex - 1.5,

    y = tmp$y,
    pre = tmp$T1,

    type = tmp$type,
    P = tmp$participant,
    treat = tmp$group,
    prior_only = 1
  )

  # Missing included (intention to treat (itt)), single time-point
  tmp <- wide_complete |>
    select(participant:T1, y = T2) |>
    mutate(
      y_miss = if_else(is.na(y), 0, 1),
      pre_miss = if_else(is.na(T1), 0, 1),
      # Adding placeholders
      y = if_else(is.na(y), -999, y),
      pre = if_else(is.na(T1), -999, T1)
    ) |>
    # Complete-case indicator (same definition as data-prep.R): 1 if the
    # participant has every outcome and every pre value observed.
    mutate(.by = participant, n_y = sum(y_miss), n_pre = sum(pre_miss)) |>
    mutate(cc = if_else(n_y == max(n_y) & n_pre == max(n_pre), 1L, 0L))

  dl_itt_single <- list(
    N = nrow(tmp),
    N_P = length(unique(tmp$participant)),
    K = length(unique(tmp$type)),
    N_TX = length(unique(tmp$group)),

    y = tmp$y,
    pre = tmp$pre,
    y_obs = tmp$y,
    pre_obs = tmp$pre,

    sex = tmp$sex - 1.5,
    sex_idx = as.integer(tmp$sex),
    cc = tmp$cc,

    y_miss = tmp$y_miss,
    pre_miss = tmp$pre_miss,

    type = tmp$type,
    P = tmp$participant,
    treat = tmp$group,
    prior_only = 1
  )

  # missing included (itt), both time-points
  tmp <- wide_complete |>
    pivot_longer(cols = c(T2, T3), values_to = "y", names_to = "time") |>

    mutate(time = dense_rank(time)) |>

    mutate(
      y_miss = if_else(is.na(y), 0, 1),
      pre_miss = if_else(is.na(T1), 0, 1),
      # Adding placeholders
      y = if_else(is.na(y), -999, y),
      pre = if_else(is.na(T1), -999, T1)
    ) |>
    # Complete-case indicator (same definition as data-prep.R): 1 if the
    # participant has every outcome and every pre value observed. Required
    # by the `cc` data variable in models s7-s10.
    mutate(.by = participant, n_y = sum(y_miss), n_pre = sum(pre_miss)) |>
    mutate(cc = if_else(n_y == max(n_y) & n_pre == max(n_pre), 1L, 0L))

  dl_itt <- list(
    N = nrow(tmp),
    N_P = length(unique(tmp$participant)),
    K = length(unique(tmp$type)),
    N_TX = length(unique(tmp$group)),

    y = tmp$y,
    pre = tmp$pre,
    y_obs = tmp$y,
    pre_obs = tmp$pre,

    sex = tmp$sex - 1.5,
    # Sex index (1/2) required by s7 and s9; sex is already dense-ranked.
    sex_idx = as.integer(tmp$sex),
    cc = tmp$cc,

    y_miss = tmp$y_miss,
    pre_miss = tmp$pre_miss,
    time = tmp$time,

    type = tmp$type,
    P = tmp$participant,
    treat = tmp$group,
    prior_only = 1
  )

  out_list <- list(
    raw_df = dstd,
    dl_compl_single = dl_compl_single,
    dl_itt_single = dl_itt_single,
    dl_itt = dl_itt
  )
}


# A wrapper for sim_type2
#
# This function is a wrapper for the sim_type2 function to allow for repeated
# simulation using slightly different correlation structures and fixed
# "true" effects. The strategy is to hard-code most parameters and add
# randomness noise to some aspects of the generative model.
#
# group1, group/treatment effects at follow-up time point 1.
# group2, group/treatment effects at follow-up time point 2.
# loadings_beta, a matrix with parameters for the beta distribution
#   for loadings. Will produce numbers between 0 and 1.
# c_corr, correlations between time-points within participant.
# N, number of participants
# K, number of
#
sim_wrap2 <- function(
  group1 = c(0.4, 0.4, rep(0, K - 2)),
  group2 = c(0.4, 0.4, rep(0, K - 2)),
  loadings_beta = matrix(c(2, 10, 1, 10, 1, 10), ncol = 2, byrow = TRUE),
  c_corr = rbeta(K, 5, 10),
  N = 148,
  K = 8
) {
  N <- N # Number of participants
  K <- K # number of measurement types

  # Provide error message when effect do not coorespond to K, number of
  # measurements.
  if (length(group1) != K)
    stop("True effects (time == post1) must have length K")
  if (length(group2) != K)
    stop("True effects (time == post2) must have length K")

  # The formula for the pre-intervention variables
  form_pre <- ~ 1 + age + sex

  # beta_pre is a p * K matrix of parameter (p in 1 to P)
  # values for each type (k in 1:K).
  # fmt: skip
  betas_pre <- matrix(
    c(rep(0, K), # intercepts
      runif(K, -0.2, 0.2),  # age
      runif(K, -0.2, 0.2)  # sex
    ),
    nrow = 3,
    ncol = K,
    byrow = TRUE
  )

  # The same formula for both time-points
  form_main <- ~ pre + group + sex

  # betas matrix are p * K
  # fmt: skip

  beta_pre <- runif(K, 0.7, 0.95)

  betas_t1 <- matrix(
    c(
      rep(0, K), # intercepts
      beta_pre, # pre_k effect
      group1, # group
      rep(0, K) # sex
    ),
    nrow = 4,
    ncol = K,
    byrow = TRUE
  )

  # fmt: skip
  betas_t2 <-  matrix(
    c(rep(0, K), # intercepts
      beta_pre,  # pre_k effect 
      group2,  # group
      rep(0, K)   # sex
    ),
    nrow = 4,
    ncol = K,
    byrow = TRUE
  )

  # Predictor variables, used in *all* sims
  preds_all <- data.frame(
    group = rep(c(0, 1), each = N / 2),
    sex = rep(c(0, 1), N / 2),
    age = rnorm(N)
  )

  # Loadings for pre-intervention loading matrix
  # Assuming that

  Lambda1 <- matrix(
    rbeta(2 * K, loadings_beta[1, 1], loadings_beta[1, 2]),
    ncol = 2,
    byrow = TRUE
  )

  # Error in the factor model
  psi1 <- 1 - rowSums(Lambda1)

  # Transform to loadings
  Lambda1 <- sqrt(Lambda1)
  # Add negative correlations
  Lambda1[1, ] <- -Lambda1[1, ]
  # Pre intervention SD's per type (length K)
  S1 <- rep(1, K)

  Lambda2a <- matrix(
    rbeta(4 * K, loadings_beta[2, 1], loadings_beta[2, 2]),
    ncol = 4,
    byrow = TRUE
  )

  Lambda2b <- matrix(
    rbeta(4 * K, loadings_beta[3, 1], loadings_beta[3, 2]),
    ncol = 4,
    byrow = TRUE
  )

  # Error in the factor model
  psi2a <- 1 - rowSums(Lambda2a)
  psi2b <- 1 - rowSums(Lambda2b)

  # Transform to loadings
  Lambda2a <- sqrt(Lambda2a)
  Lambda2b <- sqrt(Lambda2b)

  # Add some negative correlations
  Lambda2a[3, ] <- -Lambda2a[3, ]
  Lambda2b[3, ] <- -Lambda2b[3, ]

  # Pre intervention SD's per type (length K)
  s2a <- rep(1, K)
  s2b <- rep(1, K)

  Lambda2 <- list(Lambda2a, Lambda2b)
  psi2 <- list(psi2a, psi2b)
  S2 <- list(s2a, s2b)

  # Between time-point (residual) correlations within each type (length K)
  c_cross <- c_corr

  # missing_form_t2,  formula specifying the predictor for missing values at t2
  miss_t2 <- ~ 1 + age
  miss_t4 <- ~ 1 + age

  # Parameter values for the form complete missingness
  miss_beta_t2 <- c(-2, -0.5)
  miss_beta_t4 <- c(-2, -0.5)

  # Random missing proportion
  miss_prop_pre <- 0.005 # e.g., 2% missing of all rows
  miss_prop_y <- 0.02 # e.g., 5% miss of all rows post-int

  simout <- sim_type2(
    n = N,
    K = K,

    form_pre = form_pre,
    betas_pre = betas_pre,

    form = form_main,
    preds = preds_all,
    betas = list(betas_t1, betas_t2),

    Lambda1 = Lambda1,
    psi1 = psi1,
    S1 = S1,

    Lambda2 = Lambda2,
    psi2 = psi2,
    S2 = S2,

    C_cross = c_cross,
    missing_form_t2 = miss_t2,
    missing_form_t4 = miss_t4,
    missing_beta_t2 = miss_beta_t2,
    missing_beta_t4 = miss_beta_t4,
    missing_pre_prop = miss_prop_pre,
    missing_y_prop = miss_prop_y,
    verbose = FALSE
  )

  d2_full <- struct_sim(simout$full, verbose = FALSE)
  d2_miss <- struct_sim(simout$miss, verbose = FALSE)

  out <- list(
    full = d2_full,
    miss = d2_miss,
    true_effects = matrix(c(group1, group2), ncol = K, byrow = TRUE)
  )

  return(out)
}


# sim_type3 ----------------------------------------------------------------
# Heteroskedastic variant of sim_type2: residual SDs differ by treatment,
# but the correlation structure (from Lambda2, psi2, C_cross) is shared
# across treatments. Matches exactly what model s10 assumes: one L_Omega,
# treatment-specific sigmas.
#
# Differences from sim_type2:
#   - S2 argument shape: list of length 2 (per time), where each element is
#     a K x N_TX matrix. Column tx is the SD vector for treatment tx.
#   - Residuals are drawn from the correlation matrix (unit scales), then
#     rescaled per participant based on that participant's treatment.
#   - New treat_col argument picks the treatment column from preds.
#   - prop_noise, proportion of simulated values that should be added extra
#     noise.
#   - sd_noise, the added noise is drawn from a normal distribution with mean
#     0 and SD of sd_noise.
#
# Everything else -- pre-intervention draws, design matrix, mu_mat,
# missingness, and output structure -- is identical to sim_type2.
sim_type3 <- function(
  n,
  K,
  form_pre,
  betas_pre,
  form,
  betas,
  preds,
  Lambda1,
  psi1,
  S1,
  Lambda2,
  psi2,
  S2, # list of 2; each element a K x N_TX matrix
  C_cross,
  missing_form_t2 = NULL,
  missing_form_t4 = NULL,
  missing_beta_t2,
  missing_beta_t4,
  missing_y_prop = 0,
  missing_pre_prop = 0,
  prop_noise = 0,
  sd_noise = 0,
  treat_col = "group",
  verbose = TRUE
) {
  if (nrow(preds) != n) stop("n is not equal to nrow(preds)")
  if (nrow(Lambda1) != K) stop("Lambda1 must have K rows")
  if (length(psi1) != K) stop("psi1 must be of length K")
  if (any(psi1 <= 0)) stop("psi1 values must be > 0")
  if (length(S1) != K) stop("S1 must be of length K")

  if (!is.list(Lambda2) || length(Lambda2) != 2)
    stop("Lambda2 must be a list of length 2")
  if (!is.list(psi2) || length(psi2) != 2)
    stop("psi2 must be a list of length 2")
  if (!is.list(S2) || length(S2) != 2) stop("S2 must be a list of length 2")
  if (!is.list(betas) || length(betas) != 2)
    stop("betas must be a list of length 2")

  for (t in 1:2) {
    if (!is.matrix(S2[[t]]) || nrow(S2[[t]]) != K)
      stop(sprintf("S2[[%d]] must be a K x N_TX matrix", t))
    if (nrow(Lambda2[[t]]) != K)
      stop(sprintf("Lambda2[[%d]] must have K rows", t))
    if (any(psi2[[t]] <= 0)) stop(sprintf("psi2[[%d]] values must be > 0", t))
    if (!is.matrix(betas[[t]])) stop(sprintf("betas[[%d]] must be a matrix", t))
  }
  N_TX <- ncol(S2[[1]])
  if (ncol(S2[[2]]) != N_TX)
    stop("S2[[1]] and S2[[2]] must have the same number of columns (N_TX)")

  if (!(treat_col %in% colnames(preds)))
    stop(sprintf("treat_col '%s' not found in preds", treat_col))
  treat_raw <- preds[[treat_col]]
  treat_idx <- as.integer(as.factor(treat_raw)) # 1..N_TX
  if (max(treat_idx) > N_TX)
    stop("Number of treatment levels in preds exceeds S2 column count")

  if (length(C_cross) != K) stop("C_cross must be of length K")
  if (any(abs(C_cross) > 1)) stop("C_cross values must be in [-1, 1]")

  # Pre-intervention (identical to sim_type2) ---------------------------------
  Sigma1 <- build_sigma_factor(Lambda1, psi1, S1)
  Z_pre <- MASS::mvrnorm(n, rep(0, K), Sigma1)
  desig_pre <- model.matrix(form_pre, data = preds)
  mu_pre <- desig_pre %*% betas_pre
  pre <- Z_pre + mu_pre
  colnames(pre) <- paste0("type", 1:K)

  # Design matrix (identical to sim_type2) ------------------------------------
  form_terms <- attr(terms(form), "term.labels")
  has_pre <- "pre" %in% form_terms

  if (has_pre) {
    pre_df <- as.data.frame(pre)
    colnames(pre_df) <- paste0("pre", 1:K)
    preds_expanded <- cbind(preds, pre_df)
    design <- lapply(1:K, function(k) {
      form_k <- as.formula(gsub("\\bpre\\b", paste0("pre", k), deparse(form)))
      model.matrix(form_k, data = preds_expanded)
    })
  } else {
    design <- replicate(K, model.matrix(form, data = preds), simplify = FALSE)
  }

  p <- ncol(design[[1]])
  for (t in 1:2) {
    if (!all(dim(betas[[t]]) == c(p, K)))
      stop(sprintf(
        "betas[[%d]] must be %d x %d (predictors x K), got %d x %d",
        t,
        p,
        K,
        nrow(betas[[t]]),
        ncol(betas[[t]])
      ))
  }

  mu_t1 <- sapply(1:K, function(k) design[[k]] %*% betas[[1]][, k])
  mu_t2 <- sapply(1:K, function(k) design[[k]] %*% betas[[2]][, k])
  mu_mat <- cbind(mu_t1, mu_t2)

  # Correlation-scale 2K x 2K block matrix -----------------------------------
  unit_S <- rep(1, K)
  R2_t1 <- build_sigma_factor(Lambda2[[1]], psi2[[1]], unit_S)
  R2_t2 <- build_sigma_factor(Lambda2[[2]], psi2[[2]], unit_S)
  R_cross <- diag(C_cross)
  R_full <- rbind(cbind(R2_t1, R_cross), cbind(t(R_cross), R2_t2))

  eigs <- eigen(R_full, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigs <= 0))
    stop("R_full is not positive definite - check C_cross vs Lambda2/psi2")

  # Standardized residuals and per-participant scaling -----------------------
  Z_std <- MASS::mvrnorm(n, rep(0, 2 * K), R_full)

  # Row i of s_mat is c(S2[[1]][, treat_idx[i]], S2[[2]][, treat_idx[i]])
  s_t1 <- t(S2[[1]])[treat_idx, , drop = FALSE] # n x K
  s_t2 <- t(S2[[2]])[treat_idx, , drop = FALSE] # n x K
  s_mat <- cbind(s_t1, s_t2) # n x 2K

  Z <- Z_std * s_mat
  Y <- Z + mu_mat

  y_t1 <- Y[, 1:K]
  y_t2 <- Y[, (K + 1):(2 * K)]
  colnames(y_t1) <- paste0("type", 1:K)
  colnames(y_t2) <- paste0("type", 1:K)

  # Add noise to prop_noise proportion of observations (post-treatment)
  if (prop_noise > 0 && sd_noise > 0) {
    n_y1 <- length(y_t1)
    n_y2 <- length(y_t2)
    n_total <- n_y1 + n_y2
    n_noise <- round(prop_noise * n_total, 0)
    if (n_noise > 0) {
      idx <- sample.int(n_total, n_noise)
      eps <- rnorm(n_noise, 0, sd_noise)
      flat <- c(as.vector(y_t1), as.vector(y_t2))
      flat[idx] <- flat[idx] + eps
      y_t1 <- matrix(
        flat[1:n_y1],
        n,
        K,
        dimnames = list(NULL, paste0("type", 1:K))
      )
      y_t2 <- matrix(
        flat[(n_y1 + 1):n_total],
        n,
        K,
        dimnames = list(NULL, paste0("type", 1:K))
      )
    }
  }

  if (verbose) {
    cat("Heteroskedastic sim (sim_type3). Treatment-specific SDs:\n")
    for (t in 1:2) {
      cat(sprintf("Time %d (rows = type, cols = treatment):\n", t))
      print(round(S2[[t]], 3))
    }
    cat("Implied pre-intervention correlation matrix:\n")
    print(round(cov2cor(Sigma1), 2))
    cat("\nImplied post-int residual correlation (t1, unit scale):\n")
    print(round(R2_t1, 2))
    cat("\nImplied post-int residual correlation (t2, unit scale):\n")
    print(round(R2_t2, 2))
  }

  # Missing values (identical to sim_type2) -----------------------------------
  miss_t2_ind <- rep(0L, n)
  miss_t4_ind <- rep(0L, n)

  if (!is.null(missing_form_t2)) {
    X_miss_t2 <- model.matrix(missing_form_t2, data = preds)
    p_miss_t2 <- plogis(X_miss_t2 %*% missing_beta_t2)
    miss_t2_ind <- rbinom(n, 1, p_miss_t2)
  }
  if (!is.null(missing_form_t4)) {
    X_miss_t4 <- model.matrix(missing_form_t4, data = preds)
    p_miss_t4 <- plogis(X_miss_t4 %*% missing_beta_t4)
    miss_t4_ind <- rbinom(n, 1, p_miss_t4)
    miss_t4_ind[miss_t2_ind == 1] <- 1
  }

  y_t1_miss <- y_t1
  y_t2_miss <- y_t2
  y_t1_miss[miss_t2_ind == 1, ] <- NA
  y_t2_miss[miss_t4_ind == 1, ] <- NA

  # Collect output (identical to sim_type2) -----------------------------------
  preds$participant <- 1:n

  build_long <- function(y_t1_m, y_t2_m, pre_m) {
    pre_df <- data.frame(pre_m) |>
      dplyr::mutate(participant = 1:n) |>
      tidyr::pivot_longer(
        cols = tidyselect::starts_with("type"),
        names_to = "type",
        values_to = "pre"
      )
    y_df <- dplyr::bind_rows(
      data.frame(y_t1_m) |>
        dplyr::mutate(participant = 1:n, time = 1) |>
        tidyr::pivot_longer(
          cols = tidyselect::starts_with("type"),
          names_to = "type",
          values_to = "y"
        ),
      data.frame(y_t2_m) |>
        dplyr::mutate(participant = 1:n, time = 2) |>
        tidyr::pivot_longer(
          cols = tidyselect::starts_with("type"),
          names_to = "type",
          values_to = "y"
        )
    )
    y_df |>
      inner_join(pre_df, by = c("participant", "type")) |>
      inner_join(preds, by = "participant")
  }

  full <- build_long(y_t1, y_t2, pre)
  miss <- build_long(y_t1_miss, y_t2_miss, pre)

  if (missing_y_prop > 0) {
    idx_y <- sample(seq_len(nrow(miss)), round(nrow(miss) * missing_y_prop, 0))
    miss$y[idx_y] <- NA
  }
  if (missing_pre_prop > 0) {
    idx_pre <- sample(
      seq_len(nrow(miss)),
      round(nrow(miss) * missing_pre_prop, 0)
    )
    miss$pre[idx_pre] <- NA
  }

  return(list(full = full, miss = miss))
}


# sim_wrap3 ----------------------------------------------------------------
# Wrapper for sim_type3 exposing only the knobs that matter for a
# simulation sweep. Interface mirrors sim_wrap2 with one addition:
#
# sigma_ratio,  K-vector of SD ratios treat2 / treat1 per type. Values
#   equal to 1 give a homoskedastic cell; values != 1 give heteroskedastic
#   cells. Default mixes both so recovery can be tested across the
#   null and non-null regimes in one simulation.
# sigma_base,   K-vector of baseline (treat 1) SDs per type. Default rep(1, K).
#
# prop_noise, proportion of observed measurements with additional noise.
# sd_noise, SD of the additional noise.
#
# The cross-time correlation C_cross and pre-intervention structure are
# handled as in sim_wrap2.
sim_wrap3 <- function(
  group1 = c(0.4, 0.4, rep(0, K - 2)),
  group2 = c(0.4, 0.4, rep(0, K - 2)),
  loadings_beta = matrix(c(2, 10, 1, 10, 1, 10), ncol = 2, byrow = TRUE),
  c_corr = rbeta(K, 5, 10),
  sigma_base = rep(1, K),
  sigma_ratio = rep(c(1.5, 0.6, 1.0, 1.0, 1.3, 0.8, 1.0, 1.0), length.out = K),
  N = 148,
  K = 8,
  prop_noise = 0.025,
  sd_noise = 1.5
) {
  if (length(group1) != K) stop("group1 must have length K")
  if (length(group2) != K) stop("group2 must have length K")
  if (length(sigma_base) != K) stop("sigma_base must have length K")
  if (length(sigma_ratio) != K) stop("sigma_ratio must have length K")

  form_pre <- ~ 1 + age + sex

  betas_pre <- matrix(
    c(rep(0, K), runif(K, -0.2, 0.2), runif(K, -0.2, 0.2)),
    nrow = 3,
    ncol = K,
    byrow = TRUE
  )

  form_main <- ~ pre + group + sex
  beta_pre <- runif(K, 0.7, 0.95)

  betas_t1 <- matrix(
    c(rep(0, K), beta_pre, group1, rep(0, K)),
    nrow = 4,
    ncol = K,
    byrow = TRUE
  )
  betas_t2 <- matrix(
    c(rep(0, K), beta_pre, group2, rep(0, K)),
    nrow = 4,
    ncol = K,
    byrow = TRUE
  )

  preds_all <- data.frame(
    group = rep(c(0, 1), each = N / 2),
    sex = rep(c(0, 1), N / 2),
    age = rnorm(N)
  )

  # Pre-intervention loading structure (identical to sim_wrap2)
  Lambda1 <- matrix(
    rbeta(2 * K, loadings_beta[1, 1], loadings_beta[1, 2]),
    ncol = 2,
    byrow = TRUE
  )
  psi1 <- 1 - rowSums(Lambda1)
  Lambda1 <- sqrt(Lambda1)
  Lambda1[1, ] <- -Lambda1[1, ]
  S1 <- rep(1, K)

  Lambda2a <- matrix(
    rbeta(4 * K, loadings_beta[2, 1], loadings_beta[2, 2]),
    ncol = 4,
    byrow = TRUE
  )
  Lambda2b <- matrix(
    rbeta(4 * K, loadings_beta[3, 1], loadings_beta[3, 2]),
    ncol = 4,
    byrow = TRUE
  )
  psi2a <- 1 - rowSums(Lambda2a)
  psi2b <- 1 - rowSums(Lambda2b)
  Lambda2a <- sqrt(Lambda2a)
  Lambda2b <- sqrt(Lambda2b)
  Lambda2a[3, ] <- -Lambda2a[3, ]
  Lambda2b[3, ] <- -Lambda2b[3, ]

  Lambda2 <- list(Lambda2a, Lambda2b)
  psi2 <- list(psi2a, psi2b)

  # Per-treatment post-intervention SDs: K x 2 matrices (col 1 = treat 1)
  S2_t1 <- cbind(sigma_base, sigma_base * sigma_ratio)
  S2_t2 <- cbind(sigma_base, sigma_base * sigma_ratio)
  S2 <- list(S2_t1, S2_t2)

  c_cross <- c_corr

  miss_t2 <- ~ 1 + age
  miss_t4 <- ~ 1 + age
  miss_beta_t2 <- c(-2, -0.5)
  miss_beta_t4 <- c(-2, -0.5)
  miss_prop_pre <- 0.005
  miss_prop_y <- 0.02

  simout <- sim_type3(
    n = N,
    K = K,
    form_pre = form_pre,
    betas_pre = betas_pre,
    form = form_main,
    preds = preds_all,
    betas = list(betas_t1, betas_t2),
    Lambda1 = Lambda1,
    psi1 = psi1,
    S1 = S1,
    Lambda2 = Lambda2,
    psi2 = psi2,
    S2 = S2,
    C_cross = c_cross,
    missing_form_t2 = miss_t2,
    missing_form_t4 = miss_t4,
    missing_beta_t2 = miss_beta_t2,
    missing_beta_t4 = miss_beta_t4,
    missing_pre_prop = miss_prop_pre,
    missing_y_prop = miss_prop_y,
    prop_noise = prop_noise,
    sd_noise = sd_noise,
    treat_col = "group",
    verbose = FALSE
  )

  d3_full <- struct_sim(simout$full, verbose = FALSE)
  d3_miss <- struct_sim(simout$miss, verbose = FALSE)

  out <- list(
    full = d3_full,
    miss = d3_miss,
    true_effects = matrix(c(group1, group2), ncol = K, byrow = TRUE),
    true_sigma = S2 # K x 2 matrices per time, for recovery checks
  )

  return(out)
}
