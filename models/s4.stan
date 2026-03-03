/* 
Including missing data data imputation in joint model. 
*/


data {
  int<lower=0> N;                             // Number of obs. 
  int<lower=0> N_P;                           // Number of participants
  int<lower=0> K;                             // Number of measures
  
  vector[N] y_obs;                            // Data (incl missing w. placehol)
  vector[N] pre_obs;                          // Pre-int dat (incl miss)
  
  array[N] int<lower=0, upper=1> y_miss;      // Missing data indicator (miss=0)
  array[N] int<lower=0, upper=1> pre_miss;    // Miss dat ind for pre (miss=0)
  
  array[N] int<lower=1, upper=K> type;        // Measurement type
  array[N] int<lower=1, upper=N_P> P;         // Participant 
  array[N] int<lower=1, upper=2> treat;       // Treatment index
  
  array[N] real<lower=-0.5, upper=0.5> sex;   // Sex
  
  int<lower=0, upper=1> prior_only;           // Prior only flag
}
transformed data {
  int n_miss_y = N - sum(y_miss);             // Count the number of missing
  int n_miss_pre = N - sum(pre_miss);
}
parameters {
  
  // Parameters for the main model 
  
  // Type (k)
  matrix[K, 2] alpha;
  vector[K] beta;
  vector[K] beta_sex;
  
  // Correlation matrix and sigma k
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0>[K] sigma;
  
  // Missing data parameters
  vector[n_miss_y] y_impute;
  vector[n_miss_pre] pre_impute;
  
  vector[K] mu_pre;
  vector<lower=0>[K] sigma_pre;
  cholesky_factor_corr[K] L_Omega_pre;
}
transformed parameters {
  vector[N] y;
  vector[N] pre;
  matrix[N_P, K] mu_wide;
  
  // Fill y and pre
  {
    int idx_y = 1;
    int idx_pre = 1;
    
    // If missing values, use parameter
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
  
  // Compute mu_wide using pre
  for (i in 1:N) {
    mu_wide[P[i], type[i]] = alpha[type[i], treat[i]] + 
    beta[type[i]] * pre[i] + 
    beta_sex[type[i]] * sex[i];
  }
}
model {
  // Priors
  to_vector(alpha) ~ normal(0, 1);
  beta ~ normal(1, 0.5);
  beta_sex ~ normal(0, 0.5);
  sigma ~ exponential(1);
  L_Omega ~ lkj_corr_cholesky(4);
  
  mu_pre ~ normal(0, 1);
  sigma_pre ~ exponential(1);
  L_Omega_pre ~ lkj_corr_cholesky(4);
  
  // Weak prior on imputed values
  y_impute ~ normal(0, 2);
  pre_impute ~ normal(0, 2);
  
  
  
  if (prior_only == 0) {
    // Pre-value imputation model
    matrix[K, K] L_Sigma_pre = diag_pre_multiply(sigma_pre, L_Omega_pre);
    matrix[N_P, K] pre_wide;
    for (i in 1:N) {
      pre_wide[P[i], type[i]] = pre[i];
    }
    for (p in 1:N_P) {
      pre_wide[p]' ~ multi_normal_cholesky(mu_pre, L_Sigma_pre);
    }
    
    // Outcome model
    matrix[K, K] L_Sigma = diag_pre_multiply(sigma, L_Omega);
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
  matrix[K, K] Omega_pre = multiply_lower_tri_self_transpose(L_Omega_pre);
  
  vector[N] y_rep;
  vector[N_P] log_lik_part;
  vector[N] y_expect;
  
  {
    matrix[K, K] L_Sigma = diag_pre_multiply(sigma, L_Omega);
    matrix[N_P, K] y_rep_wide;
    
    for (p in 1:N_P) {
      y_rep_wide[p] = multi_normal_cholesky_rng(mu_wide[p]', L_Sigma)';
      log_lik_part[p] = 0;
    }
    
    // Log-lik only for observed y
    for (i in 1:N) {
      y_rep[i] = y_rep_wide[P[i], type[i]];
      y_expect[i] = mu_wide[P[i], type[i]];
      
      if (y_miss[i] == 1) {
        log_lik_part[P[i]] += normal_lpdf(y[i] | mu_wide[P[i], type[i]], sigma[type[i]]);
      }
    }
  }
}

