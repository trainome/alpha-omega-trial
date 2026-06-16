/* Three-timepoint multilevel spline model for OGTT data.
   brms-style parameterization (see. m5 in archive/ogtt-analysis.qmd):

     stglucose ~ 0 + g + bs(min, df = K):g + (1 + bs(min, df = K):time | participant)

   Population: a cell mean per group g (0 + g) plus a spline curve per group
   (bs(min):g). Group g indexes the population cells, i.e. time x treatment

   Participant level: a random intercept plus random spline coefficients that may
   differ by time (1 + bs(min):time | participant), correlated across the whole
   block and non-centered via a Cholesky factor. The basis B is bs(min, df = K)
   WITHOUT an intercept column (set intercept = FALSE). */

data {
  int<lower=1> N;                       // observations
  vector[N] y;                          // glucose (standardized)

  int<lower=1> K;                       // spline basis columns (e.g. df = 4)
  matrix[N, K] B;                       // bs(min, df = K) basis, intercept = FALSE

  // population grouping (time x treatment)
  int<lower=1> N_G;                     // number of population cells
  array[N] int<lower=1, upper=N_G> G;   // population cell index per obs

  // participant grouping
  int<lower=1> N_P;                     // number of participants
  array[N] int<lower=1, upper=N_P> P;   // participant id per obs

  // time grouping for the participant-level spline
  int<lower=1> N_T;                     // number of time points (e.g. 3)
  array[N] int<lower=1, upper=N_T> Tt;  // time index per obs
  // Prior only flag 
  int<lower=0, upper=1> prior_only;     // Prior only, if 1 -> No likelihood.   
}
transformed data {
  int M = 1 + K * N_T;                  // group-level effects: intercept + spline-per-time
}
parameters {
  // population level
  vector[N_G] a_g;                      // cell means (0 + g)
  matrix[K, N_G] b_g;                   // spline coefficients per cell (bs:g)
  real<lower=0> sigma;

  // participant level (non-centered, correlated)
  vector<lower=0>[M] sd_p;              // sd of each group-level effect
  cholesky_factor_corr[M] L_p;          // correlations among group-level effects
  matrix[M, N_P] z_p;                   // standardized group-level effects
}
transformed parameters {
  matrix[N_P, M] r_p = (diag_pre_multiply(sd_p, L_p) * z_p)';
  vector[N] mu;

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
