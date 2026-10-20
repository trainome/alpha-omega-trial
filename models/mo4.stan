/*

Three-timepoint multilevel spline model for OGTT data with missing-data
imputation and participant-level intercepts per time point (cf. mo3.stan).

Population: a cell mean per group g (time x treatment) plus a spline curve
per group, as in mo3.

Participant level: in mo3 the participant intercept was common to all time
points. Since the basis B has no intercept column, all basis columns are zero
at min = 0, and participant deviations in fasting glucose were therefore
constrained to be equal at all time points. Here the participant intercept
at time t is decomposed into a shared intercept and a time-specific deviation:

  a_tp = a_p + delta_tp,   delta_tp ~ normal(0, tau_delta)

The shared intercept (a_p) and the time-specific spline deviations form the
same correlated block as in mo3 (M = 1 + K * N_T, non-centered via a Cholesky
factor). The time-specific intercept deviations (delta) are independent,
non-centered and share one sd (tau_delta). mo3 is the special case
tau_delta = 0.

Missing-data convention (as in mo3 and s9-s12): y_miss[i] == 1 means OBSERVED,
y_miss[i] == 0 means MISSING. n_miss = N - sum(y_miss).

*/

data {
  int<lower=1> N;                       // observations (observed + missing rows)
  vector[N] y_obs;                      // glucose (standardized); placeholder where missing
  array[N] int<lower=0, upper=1> y_miss; // observed indicator (1 = observed, 0 = missing)

  int<lower=1> K;                       // spline basis columns (ncol(B) in R)
  matrix[N, K] B;                       // bs(min, df = K) basis, intercept = FALSE

  // population grouping (time * treatment)
  int<lower=1> N_G;                     // number of population cells
  array[N] int<lower=1, upper=N_G> G;   // population cell index per obs

  // participant grouping
  int<lower=1> N_P;                     // number of participants
  array[N] int<lower=1, upper=N_P> P;   // participant id per obs

  // time grouping for the participant-level intercept deviation and spline
  int<lower=1> N_T;                     // number of time points (3 in the full data)
  array[N] int<lower=1, upper=N_T> Tt;  // time index per obs

  int<lower=0, upper=1> prior_only;     // 1 = ignore likelihood
}
transformed data {
  int M = 1 + K * N_T;                  // group-level effects: intercept + spline-per-time
  int n_miss = N - sum(y_miss);         // number of missing y values to impute
}
parameters {
  // population level model parameters
  vector[N_G] a_g;                      // cell means (0 + g)
  matrix[K, N_G] b_g;                   // spline coefficients per cell (bs:g)
  real<lower=0> sigma;

  // participant level model parameters (non-centered, correlated)
  vector<lower=0>[M] sd_p;              // sd of each group-level effect
  cholesky_factor_corr[M] L_p;          // correlations among group-level effects
  matrix[M, N_P] z_p;                   // standardized group-level effects

  // time-specific participant intercept deviations (non-centered)
  real<lower=0> tau_delta;              // sd of time-specific intercept deviations
  matrix[N_T, N_P] z_delta;             // standardized deviations

  // missing-data parameters
  vector[n_miss] y_impute;              // imputed (missing) glucose values
}
transformed parameters {
  vector[N] y;                          // full response (observed + imputed)
  matrix[N_P, M] r_p = (diag_pre_multiply(sd_p, L_p) * z_p)';
  matrix[N_T, N_P] delta = tau_delta * z_delta;
  vector[N] mu;

  // Reassemble the full y vector from observed values and imputed parameters.
  {
    int idx_y = 1;
    for (i in 1:N) {
      if (y_miss[i] == 1) {
        y[i] = y_obs[i];
      } else {
        y[i] = y_impute[idx_y];
        idx_y += 1;
      }
    }
  }

  for (n in 1:N) {
    // population: cell mean + cell-specific spline
    real eta = a_g[G[n]] + dot_product(B[n], b_g[, G[n]]);
    // participant: shared intercept + time-specific intercept deviation
    // + random spline for this observation's time block
    int off = 1 + (Tt[n] - 1) * K;            // start of this time's spline block in r_p
    row_vector[K] rb = r_p[P[n], (off + 1):(off + K)];
    mu[n] = eta + r_p[P[n], 1] + delta[Tt[n], P[n]] + dot_product(B[n], rb);
  }
}
model {
  // priors
  a_g ~ normal(0, 2);
  to_vector(b_g) ~ normal(0, 2);
  sigma ~ exponential(1);
  sd_p ~ exponential(1);
  L_p ~ lkj_corr_cholesky(1);
  to_vector(z_p) ~ std_normal();
  tau_delta ~ exponential(2);
  to_vector(z_delta) ~ std_normal();

  // weak prior on imputed values (always applied)
  y_impute ~ normal(0, 2);

  // likelihood (skipped under prior_only); ties imputed y to the model
  if (prior_only == 0) {
    y ~ normal(mu, sigma);
  }
}
generated quantities {
  corr_matrix[M] Cor_p = multiply_lower_tri_self_transpose(L_p);
  // implied correlation between a participant's intercepts at two time points
  real rho_int = square(sd_p[1]) / (square(sd_p[1]) + square(tau_delta));
  vector[N] y_pred = mu;
  array[N] real y_rep = normal_rng(mu, sigma);
}
