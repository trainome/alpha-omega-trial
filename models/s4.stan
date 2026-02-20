/* */

data {
  int<lower=0> N;
  int<lower=0> N_part;
  int<lower=0> K;
  
  vector[N] y_obs;
  vector[N] pre_obs;
  
  array[N] int<lower=0, upper=1> y_miss;
  array[N] int<lower=0, upper=1> pre_miss;
  
  array[N] int<lower=1, upper=K> type;
  array[N] int<lower=1, upper=N_part> P;
  array[N] int<lower=1, upper=2> treat;
  
  vector[N] age;
  array[N] int<lower=1, upper=2> sex;
  
  int<lower=0, upper=1> prior_only;
}
transformed data {
  int n_miss_y = N - sum(y_miss);
  int n_miss_pre = N - sum(pre_miss);
}
parameters {
  matrix[K, 2] alpha;
  vector[K] beta;
  vector[K] beta_age;
  vector[K] beta_sex;
  
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0>[K] sigma;
  
  // Missing data parameters
  vector[n_miss_y] y_impute;
  vector[n_miss_pre] pre_impute;
  
  // Pre-value imputation model
  vector[K] mu_pre;
  vector[K] beta_age_pre;
  vector[K] beta_sex_pre;
  vector<lower=0>[K] sigma_pre;
  cholesky_factor_corr[K] L_Omega_pre;
}
transformed parameters {
  vector[N] y;
  vector[N] pre;
  
  {
    int idx_y = 1;
    int idx_pre = 1;
    
    for (i in 1:N) {
      if (y_miss[i] == 1) {
        y[i] = y_obs[i];
      } else {
        y[i] = y_impute[idx_y];
        idx_y += 1;
      }
      
      if (pre_miss[i] == 1) {
        pre[i] = pre_obs[i];
      } else {
        pre[i] = pre_impute[idx_pre];
        idx_pre += 1;
      }
    }
  }
  
  // Build wide-format matrices (order independent)
  matrix[N_part, K] mu_wide;
  matrix[N_part, K] mu_pre_wide;
  for (i in 1:N) {
    mu_wide[P[i], type[i]] = alpha[type[i], treat[i]] 
                              + beta[type[i]] * pre[i]
                              + beta_age[type[i]] * age[i]
                              + beta_sex[type[i]] * (sex[i] - 1);
    
    mu_pre_wide[P[i], type[i]] = mu_pre[type[i]] 
                                  + beta_age_pre[type[i]] * age[i]
                                  + beta_sex_pre[type[i]] * (sex[i] - 1);
  }
}
model {
  // Outcome model priors
  to_vector(alpha) ~ normal(0, 0.5);
  beta ~ normal(1, 0.2);
  beta_age ~ normal(0, 0.2);
  beta_sex ~ normal(0, 0.2);
  sigma ~ exponential(1);
  L_Omega ~ lkj_corr_cholesky(4);
  
  // Pre-value model priors
  mu_pre ~ normal(0, 1);
  beta_age_pre ~ normal(0, 0.2);
  beta_sex_pre ~ normal(0, 0.2);
  sigma_pre ~ exponential(1);
  L_Omega_pre ~ lkj_corr_cholesky(2);
  
  if (prior_only == 0) {
    // Pre-value imputation model
    matrix[K, K] L_Sigma_pre = diag_pre_multiply(sigma_pre, L_Omega_pre);
    matrix[N_part, K] pre_wide;
    for (i in 1:N) {
      pre_wide[P[i], type[i]] = pre[i];
    }
    for (p in 1:N_part) {
      pre_wide[p]' ~ multi_normal_cholesky(mu_pre_wide[p]', L_Sigma_pre);
    }
    
    // Outcome model
    matrix[K, K] L_Sigma = diag_pre_multiply(sigma, L_Omega);
    matrix[N_part, K] y_wide;
    for (i in 1:N) {
      y_wide[P[i], type[i]] = y[i];
    }
    for (p in 1:N_part) {
      y_wide[p]' ~ multi_normal_cholesky(mu_wide[p]', L_Sigma);
    }
  }
}
generated quantities {
  matrix[K, K] Omega = multiply_lower_tri_self_transpose(L_Omega);
  matrix[K, K] Omega_pre = multiply_lower_tri_self_transpose(L_Omega_pre);
  
  vector[N] y_rep;
  vector[N_part] log_lik_part;
  
  {
    matrix[K, K] L_Sigma = diag_pre_multiply(sigma, L_Omega);
    matrix[N_part, K] y_wide;
    for (i in 1:N) {
      y_wide[P[i], type[i]] = y[i];
    }
    
    for (p in 1:N_part) {
      vector[K] y_part = multi_normal_cholesky_rng(mu_wide[p]', L_Sigma);
      log_lik_part[p] = multi_normal_cholesky_lpdf(y_wide[p]' | mu_wide[p]', L_Sigma);
      
      for (i in 1:N) {
        if (P[i] == p) {
          y_rep[i] = y_part[type[i]];
        }
      }
    }
  }
  
  vector[N] y_expect;
  for (i in 1:N) {
    y_expect[i] = mu_wide[P[i], type[i]];
  }
}

