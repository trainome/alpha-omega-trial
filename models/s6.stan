/* 
Joint model with missing data imputation, sex predictor,
time-specific intercepts, and a single 2K-dimensional MVN 
with LKJ prior on the full correlation structure.
*/

data {
  int<lower=0> N;                             // Number of obs. 
  int<lower=0> N_P;                           // Number of participants
  int<lower=0> K;                             // Number of measures
  
  vector[N] y_obs;                            // Data (incl missing w. placehol)
  vector[N] pre_obs;                          // Pre-int dat (incl miss)
  
  array[N] int<lower=0, upper=1> y_miss;      // Missing data indicator (miss=0)
  array[N] int<lower=0, upper=1> pre_miss;    // Miss dat ind for pre (miss=0)
  
  array[N] int<lower=1, upper=2> time;        // Time indicator
  array[N] int<lower=1, upper=K> type;        // Measurement type
  array[N] int<lower=1, upper=N_P> P;         // Participant 
  array[N] int<lower=1, upper=2> treat;       // Treatment index
  
  array[N] real<lower=-0.5, upper=0.5> sex;   // Sex
  
  int<lower=0, upper=1> prior_only;           // Prior only flag
}
transformed data {
  int n_miss_y = N - sum(y_miss);
  int n_miss_pre = N - sum(pre_miss);
  int M = 2 * K;
  
  // Participant observation index for generated quantities
  array[N_P, M] int obs_idx;
  array[N_P] int obs_count;
  for (p in 1:N_P) obs_count[p] = 0;
  for (i in 1:N) {
    obs_count[P[i]] += 1;
    obs_idx[P[i], obs_count[P[i]]] = i;
  }
}
parameters {
  // Time-specific intercepts: alpha[time][type, treat]
  array[2] matrix[K, 2] alpha;
  vector[K] beta;
  // Time-specific sex parameters
  matrix[K, 2] beta_sex;
  
  // Single 2K-dimensional correlation and scale
  cholesky_factor_corr[M] L_Omega;
  vector<lower=0>[M] sigma;               // sigma[1:K] = time 1, sigma[K+1:M] = time 2
  
  // Missing data parameters
  vector[n_miss_y] y_impute;
  vector[n_miss_pre] pre_impute;
  
  // Pre-value imputation model
  vector[K] mu_pre;
  vector<lower=0>[K] sigma_pre;
  cholesky_factor_corr[K] L_Omega_pre;
}
transformed parameters {
  vector[N] y;
  vector[N] pre;
  matrix[N_P, K] mu_wide_1 = rep_matrix(0.0, N_P, K); 
  matrix[N_P, K] mu_wide_2 = rep_matrix(0.0, N_P, K); 
  
  // Fill y and pre
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
  
  // Compute mu_wide for each time-point
  for (i in 1:N) {
    real mu_i = alpha[time[i]][type[i], treat[i]] +
                beta[type[i]] * pre[i] +
                beta_sex[type[i], time[i]] * sex[i];
    if (time[i] == 1)
      mu_wide_1[P[i], type[i]] = mu_i;
    else
      mu_wide_2[P[i], type[i]] = mu_i;
  }
}
model {
  // Priors
  for (t in 1:2) to_vector(alpha[t]) ~ normal(0, 1);
  beta ~ normal(1, 0.4);
  to_vector(beta_sex) ~ normal(0, 0.5);
  
  sigma ~ exponential(1);
  L_Omega ~ lkj_corr_cholesky(2); 
  mu_pre ~ normal(0, 1);
  sigma_pre ~ exponential(1);
  L_Omega_pre ~ lkj_corr_cholesky(2);
  
  // Weak prior on imputed values
  //if (prior_only == 1) {
    y_impute ~ normal(0, 2);
    pre_impute ~ normal(0, 2);
  //}
  
  if (prior_only == 0) {
    // Pre-value imputation model
    matrix[K, K] L_Sigma_pre = diag_pre_multiply(sigma_pre, L_Omega_pre);
    matrix[N_P, K] pre_wide;
    for (i in 1:N) {
      pre_wide[P[i], type[i]] = pre[i];
    }
    
    array[N_P] vector[K] pre_vec;
    for (p in 1:N_P) pre_vec[p] = pre_wide[p]';
    
    pre_vec ~ multi_normal_cholesky(mu_pre, L_Sigma_pre);
    
    // Outcome model: single 2K-dimensional MVN
    matrix[M, M] L_Sigma = diag_pre_multiply(sigma, L_Omega);
    
    matrix[N_P, K] y_wide_1;
    matrix[N_P, K] y_wide_2;
    for (i in 1:N) {
      if (time[i] == 1) {
        y_wide_1[P[i], type[i]] = y[i];
      } else {
        y_wide_2[P[i], type[i]] = y[i];
      }
    }
    
    array[N_P] vector[M] y_arr;
    array[N_P] vector[M] mu_arr;
    for (p in 1:N_P) {
      y_arr[p][1:K]      = y_wide_1[p]';
      y_arr[p][(K+1):M]  = y_wide_2[p]';
      mu_arr[p][1:K]     = mu_wide_1[p]';
      mu_arr[p][(K+1):M] = mu_wide_2[p]';
    }
    
    y_arr ~ multi_normal_cholesky(mu_arr, L_Sigma);
  }
}
generated quantities {
  // Full 2K x 2K correlation matrix
  matrix[M, M] Omega = multiply_lower_tri_self_transpose(L_Omega);
  matrix[K, K] Omega_pre = multiply_lower_tri_self_transpose(L_Omega_pre);
  
  // Extract interpretable sub-blocks
  matrix[K, K] Omega_1   = Omega[1:K, 1:K];           // Within time 1
  matrix[K, K] Omega_2   = Omega[(K+1):M, (K+1):M];   // Within time 2
  matrix[K, K] Omega_cross = Omega[1:K, (K+1):M];     // Cross-timepoint
  
  // Per-outcome cross-timepoint correlations (diagonal of cross block)
  vector[K] rho_time;
  for (k in 1:K)
    rho_time[k] = Omega[k, K + k];
  
  vector[N] y_rep;
  vector[N_P] log_lik_part;
  vector[N] y_expect;
  
  {
    matrix[M, M] L_Sigma = diag_pre_multiply(sigma, L_Omega);
    
    matrix[N_P, K] y_wide_1 = rep_matrix(0, N_P, K);
    matrix[N_P, K] y_wide_2 = rep_matrix(0, N_P, K);
    for (i in 1:N) {
      if (time[i] == 1) {
        y_wide_1[P[i], type[i]] = y[i];
      } else {
        y_wide_2[P[i], type[i]] = y[i];
      }
    }
    
    for (p in 1:N_P) {
      vector[M] mu_stacked;
      mu_stacked[1:K]     = mu_wide_1[p]';
      mu_stacked[(K+1):M] = mu_wide_2[p]';
      
      // Posterior predictive draw
      vector[M] y_part = multi_normal_cholesky_rng(mu_stacked, L_Sigma);
      
      // Log-likelihood
      vector[M] y_stacked;
      y_stacked[1:K]     = y_wide_1[p]';
      y_stacked[(K+1):M] = y_wide_2[p]';
      log_lik_part[p] = multi_normal_cholesky_lpdf(y_stacked | mu_stacked, L_Sigma);
      
      // Map back to long format
      for (j in 1:obs_count[p]) {
        int i = obs_idx[p, j];
        if (time[i] == 1) {
          y_rep[i] = y_part[type[i]];
        } else {
          y_rep[i] = y_part[K + type[i]];
        }
      }
    }
    
    for (i in 1:N) {
      if (time[i] == 1) {
        y_expect[i] = mu_wide_1[P[i], type[i]];
      } else {
        y_expect[i] = mu_wide_2[P[i], type[i]];
      }
    }
  }
}
