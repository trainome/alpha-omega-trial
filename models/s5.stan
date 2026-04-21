/* 
Including missing data data imputation in joint model, 
adding sex as a predictor. This model alos includes both time-points. Sex
is time-specific.

The var-cov matrix here is separated by time, i.e., no between time-points
correlations are modelled (or used for imputation).

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
  int n_miss_y = N - sum(y_miss);             // Count the number of missing
  int n_miss_pre = N - sum(pre_miss);
}
parameters {
  
  // Parameters for the main model 
  
  // An array per time-point, Type (k) treatments (2)
  array[2] matrix[K, 2] alpha;
  vector[K] beta;
  // Time-specific sex parameters
  matrix[K,2] beta_sex;
  
  // Correlation matrix and sigma k for each time point
  cholesky_factor_corr[K] L_Omega_1;
  cholesky_factor_corr[K] L_Omega_2;
  vector<lower=0>[K] sigma_1;
  vector<lower=0>[K] sigma_2;
  
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
  // mu matricies for both time points
  matrix[N_P, K] mu_wide_1 = rep_matrix(0.0, N_P, K); 
  matrix[N_P, K] mu_wide_2 = rep_matrix(0.0, N_P, K); 
  
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
  
  // Compute mu_wide using pre for each time-point
  for (i in 1:N) {
  real mu_i = alpha[time[i]][type[i], treat[i]] +
            beta[type[i]] * pre[i] +
            beta_sex[type[i],time[i]] * sex[i];
  if (time[i] == 1)
    mu_wide_1[P[i], type[i]] = mu_i;
  else
    mu_wide_2[P[i], type[i]] = mu_i;
}
}
model {
  // Priors
  for (t in 1:2)  to_vector(alpha[t]) ~ normal(0, 1);
  beta ~ normal(1, 0.4);
  to_vector(beta_sex) ~ normal(0, 0.5);
  sigma_1 ~ exponential(1);
  L_Omega_1 ~ lkj_corr_cholesky(2);
  sigma_2 ~ exponential(1);
  L_Omega_2 ~ lkj_corr_cholesky(2);
  
  mu_pre ~ normal(0, 1);
  sigma_pre ~ exponential(1);
  L_Omega_pre ~ lkj_corr_cholesky(2);
  
  // Weak prior on imputed values
  // if (prior_only == 1) {
  y_impute ~ normal(0, 2);
  pre_impute ~ normal(0, 2);
  // }
  
  
  
 if (prior_only == 0) {
    matrix[K, K] L_Sigma_pre = diag_pre_multiply(sigma_pre, L_Omega_pre);
    matrix[N_P, K] pre_wide;
    for (i in 1:N) {
      pre_wide[P[i], type[i]] = pre[i];
    }

    // Convert to array of vectors for vectorized MVN
    array[N_P] vector[K] pre_vec;
    for (p in 1:N_P) pre_vec[p] = pre_wide[p]';

    pre_vec ~ multi_normal_cholesky(mu_pre, L_Sigma_pre);

    // Outcome model
    matrix[K, K] L_Sigma_1 = diag_pre_multiply(sigma_1, L_Omega_1);
    matrix[K, K] L_Sigma_2 = diag_pre_multiply(sigma_2, L_Omega_2);
    matrix[N_P, K] y_wide_1;
    matrix[N_P, K] y_wide_2;
    for (i in 1:N) {
      if (time[i] == 1) {
        y_wide_1[P[i], type[i]] = y[i];
      } else {
        y_wide_2[P[i], type[i]] = y[i];
      }
    }

    array[N_P] vector[K] y_vec_1;
    array[N_P] vector[K] y_vec_2;
    array[N_P] vector[K] mu_vec_1;
    array[N_P] vector[K] mu_vec_2;
    for (p in 1:N_P) {
      y_vec_1[p] = y_wide_1[p]';
      y_vec_2[p] = y_wide_2[p]';
      mu_vec_1[p] = mu_wide_1[p]';
      mu_vec_2[p] = mu_wide_2[p]';
    }

    y_vec_1 ~ multi_normal_cholesky(mu_vec_1, L_Sigma_1);
    y_vec_2 ~ multi_normal_cholesky(mu_vec_2, L_Sigma_2);
  }
}
generated quantities {
  
  // Back transform to correlation matrix
  matrix[K, K] Omega_1 = multiply_lower_tri_self_transpose(L_Omega_1);
  matrix[K, K] Omega_2 = multiply_lower_tri_self_transpose(L_Omega_2);
  matrix[K, K] Omega_pre = multiply_lower_tri_self_transpose(L_Omega_pre);
  
  // Vectors for y_rep, log_lik and y_exp
  vector[N] y_rep;
  vector[N_P] log_lik_part;
  vector[N] y_expect;
  
  
  // Scope this part locally, L_Sigma_1 will not be saved
  // in the output.
  {
    matrix[K, K] L_Sigma_1 = diag_pre_multiply(sigma_1, L_Omega_1);
    matrix[K, K] L_Sigma_2 = diag_pre_multiply(sigma_2, L_Omega_2);
    
    matrix[N_P, K] y_rep_wide_1;
    matrix[N_P, K] y_rep_wide_2;
    
    for (p in 1:N_P) {
      y_rep_wide_1[p] = multi_normal_cholesky_rng(mu_wide_1[p]', L_Sigma_1)';
      y_rep_wide_2[p] = multi_normal_cholesky_rng(mu_wide_2[p]', L_Sigma_2)';
      log_lik_part[p] = 0;
    }
    
    for (i in 1:N) {
      if (time[i] == 1) {
        y_rep[i]    = y_rep_wide_1[P[i], type[i]];
        y_expect[i] = mu_wide_1[P[i], type[i]];
      } else {
        y_rep[i]    = y_rep_wide_2[P[i], type[i]];
        y_expect[i] = mu_wide_2[P[i], type[i]];
      }
      // Conditional log_lik for non missing data
      if (y_miss[i] == 1) {
        int    k   = type[i];
        int    p   = P[i];
        real   mu  = (time[i] == 1) ? mu_wide_1[p, k] : mu_wide_2[p, k];
        real   sig = (time[i] == 1) ? sigma_1[k]       : sigma_2[k];
        log_lik_part[p] += normal_lpdf(y[i] | mu, sig);
      }
    }
  }
}

