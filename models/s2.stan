/* A joint model for outcome scores. 
A multivariate normal distribution is used for the participant-level outcomes
estimating a correlation structur between outcome types. 

The model uses the Cholesky factor parameterzation for stability.
See Stan User guide 
(https://mc-stan.org/docs/stan-users-guide/regression.html#multivariate-outcomes)

*/

data {
  int<lower=0> N;
  int<lower=0> N_P;
  int<lower=0> K;
  
  vector[N] y;
  vector[N] pre;
  
  array[N] int<lower=1, upper=K> type;
  array[N] int<lower=1, upper=N_P> P;
  array[N] int<lower=1, upper=2> treat;
  
  int<lower=0, upper=1> prior_only;
}
parameters {
  matrix[K, 2] alpha;
  vector[K] beta;
  vector<lower=0>[K] sigma;
  cholesky_factor_corr[K] L_Omega;
}
transformed parameters {
  matrix[N_P, K] mu_wide;
  
  for (i in 1:N) {
    mu_wide[P[i], type[i]] = alpha[type[i], treat[i]] + beta[type[i]] * pre[i];
  }
}
model {
  to_vector(alpha) ~ normal(0, 1);
  beta ~ normal(1, 0.4);
  sigma ~ exponential(1);
  L_Omega ~ lkj_corr_cholesky(2);
  
  if (prior_only == 0) {
    matrix[K, K] L_Sigma = diag_pre_multiply(sigma, L_Omega);
    
    // Build y_wide the same way
    matrix[N_P, K] y_wide;
    for (i in 1:N) {
      y_wide[P[i], type[i]] = y[i];
    }
    
    for (p in 1:N_P) {
      y_wide[p]' ~ multi_normal_cholesky(mu_wide[p]', L_Sigma);
    }
  }
}
generated quantities {
  matrix[K, K] Omega = multiply_lower_tri_self_transpose(L_Omega);
  vector[N] y_rep;
  vector[N_P] log_lik_part;

  {
    matrix[K, K] L_Sigma = diag_pre_multiply(sigma, L_Omega);
    matrix[N_P, K] y_rep_wide;

    matrix[N_P, K] y_wide;
    for (i in 1:N) {
      y_wide[P[i], type[i]] = y[i];
    }

    for (p in 1:N_P) {
      y_rep_wide[p] = multi_normal_cholesky_rng(mu_wide[p]', L_Sigma)';
      log_lik_part[p] = multi_normal_cholesky_lpdf(y_wide[p]' | mu_wide[p]', L_Sigma);
    }

    for (i in 1:N) {
      y_rep[i] = y_rep_wide[P[i], type[i]];
    }
  }
}
