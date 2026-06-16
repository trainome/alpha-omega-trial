/* 

Three-timepoint multilevel spline model for OGTT data, WITH missing-data 
imputation (cf. mo2.stan, plus the imputation setup used in s9-s12.stan).

Participant set is fixed by the complete T1 data; rows for missing post-int
OGTT measurements are still included (their design -- min/B, G/Tt,
participant/P -- is known), and the missing y values become parameters
(y_impute). The final analysis therefore averages over the posterior
uncertainty in the imputed values.

Missing-data convention (as in s9-s12): y_miss[i] == 1 means OBSERVED,
y_miss[i] == 0 means MISSING. n_miss = N - sum(y_miss).

prior_only == 1 ignores the likelihood (samples imputed values from their
weak prior and all parameters from their priors only). 

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

  // time grouping for the participant-level spline
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

  // missing-data parameters
  vector[n_miss] y_impute;              // imputed (missing) glucose values
}
transformed parameters {
  vector[N] y;                          // full response (observed + imputed)
  matrix[N_P, M] r_p = (diag_pre_multiply(sd_p, L_p) * z_p)';
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
    // participant: random intercept + random spline for this observation's time block
    int off = 1 + (Tt[n] - 1) * K;            // start of this time's spline block in r_p
    row_vector[K] rb = r_p[P[n], (off + 1):(off + K)];
    mu[n] = eta + r_p[P[n], 1] + dot_product(B[n], rb);
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

  // weak prior on imputed values (always applied)
  y_impute ~ normal(0, 2);

  // likelihood (skipped under prior_only); ties imputed y to the model
  if (prior_only == 0) {
    y ~ normal(mu, sigma);
  }
}
generated quantities {
  corr_matrix[M] Cor_p = multiply_lower_tri_self_transpose(L_p);
  vector[N] y_pred = mu;
  array[N] real y_rep = normal_rng(mu, sigma);
}
