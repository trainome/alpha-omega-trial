/* Baseline (single timepoint) multilevel spline model for OGTT data.
   brms-style parameterization (cf. m6 in archive/ogtt-analysis.qmd):

     stglucose ~ 1 + bs(min, df = K) + (1 + bs(min, df = K) | participant)

   A population mean curve (intercept + spline) plus correlated participant-level
   random intercept and spline coefficients, non-centered via a Cholesky factor
   of the correlation matrix. The basis B is bs(min, df = K) WITHOUT an intercept
   column, so the random intercept is not confounded with the spline weights.
   
   Priors are selected to be less informative than previous models 
   (archive/ogtt-analysis.qmd).
   
   */

data {
  int<lower=1> N;                       // N observations
  vector[N] y;                          // glucose (standardized)

  int<lower=1> K;                       // spline basis columns (ncol(B) in R)
  matrix[N, K] B;                       // bs(min, df = K) basis, intercept = FALSE

  int<lower=1> N_P;                     // number of participants
  array[N] int<lower=1, upper=N_P> P;   // participant id per observation
  // Prior only flag 
  int<lower=0, upper=1> prior_only;     // Prior only, if 1 -> No likelihood.           
}
transformed data {
  int M = K + 1;                        // group-level effects: intercept + spline
  // This helps to build the Cholesky correlation structure
}
parameters {
  // population level model parameters
  real a;                               // intercept
  vector[K] b;                          // population spline coefficients
  real<lower=0> sigma;                  // residual sd

  // participant level model (non-centered, correlated)
  vector<lower=0>[M] sd_p;              // sd of each group-level effect
  cholesky_factor_corr[M] L_p;          // correlations among group-level effects
  matrix[M, N_P] z_p;                   // standardized group-level effects
}
transformed parameters {
  // participant deviations: rows = participants, cols = [intercept, spline_1..K]
  matrix[N_P, M] r_p = (diag_pre_multiply(sd_p, L_p) * z_p)';
  vector[N] mu;

  for (n in 1:N) {
    // participant spline coefficients = population + deviation
    // The ' indicate the transpose, of the population spline + 
    // the participant-level deviations
    row_vector[K] bc = b' + r_p[P[n], 2:M];
    // INtercepts plus sum of splines parameters times the synthetic predictor.
    mu[n] = a + r_p[P[n], 1] + dot_product(B[n], bc);
  }
}
model {
  // priors
  a ~ normal(0, 2);
  b ~ normal(0, 2);
  sigma ~ exponential(1);
  sd_p ~ exponential(1);
  L_p ~ lkj_corr_cholesky(1);
  to_vector(z_p) ~ std_normal();

  // likelihood
  if (prior_only == 0) {
  y ~ normal(mu, sigma);
  }
}
generated quantities {
  corr_matrix[M] Cor_p = multiply_lower_tri_self_transpose(L_p);
  vector[N] y_pred = mu;
  array[N] real y_rep = normal_rng(mu, sigma);
}
