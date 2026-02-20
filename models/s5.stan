data {
  int<lower=0> N;
  int<lower=0> N_part;
  int<lower=0> K;
  
  vector[N] y_obs;
  vector[N] pre_obs;
  
  array[N] int<lower=0, upper=1> y_miss;
  array[N] int<lower=0, upper=1> pre_miss;
  
  array[N] int<lower=1, upper=K> score;
  array[N] int<lower=1, upper=N_part> P;
  array[N] int<lower=1, upper=2> treat;
  array[N] int<lower=1, upper=2> time;
  
  int<lower=0, upper=1> prior_only;
}
transformed data {
  int n_miss_y = N - sum(y_miss);
  int n_miss_pre = N - sum(pre_miss);
  
  // Precompute missingness indices
  array[n_miss_y] int miss_y_idx;
  array[n_miss_y] int miss_y_score;
  array[n_miss_y] int miss_y_time;
  array[n_miss_pre] int miss_pre_idx;
  array[n_miss_pre] int miss_pre_score;
  
  {
    int jy = 1;
    int jp = 1;
    for (i in 1:N) {
      if (y_miss[i] == 0) {
        miss_y_idx[jy] = i;
        miss_y_score[jy] = score[i];
        miss_y_time[jy] = time[i];
        jy += 1;
      }
      if (pre_miss[i] == 0) {
        miss_pre_idx[jp] = i;
        miss_pre_score[jp] = score[i];
        jp += 1;
      }
    }
  }
  
  // Precompute participant-row lookup
  array[N_part, 2 * K] int obs_idx;
  array[N_part] int obs_count;
  for (p in 1:N_part) obs_count[p] = 0;
  for (i in 1:N) {
    obs_count[P[i]] += 1;
    obs_idx[P[i], obs_count[P[i]]] = i;
  }
}
parameters {
  matrix[K, 2] alpha_mid;
  matrix[K, 2] alpha_post;
  vector[K] beta;
  
  cholesky_factor_corr[K] L_Omega_mid;
  cholesky_factor_corr[K] L_Omega_post;
  vector<lower=0>[K] sigma_mid;
  vector<lower=0>[K] sigma_post;
  
  // Non-centered imputation
  vector[n_miss_y] y_impute_raw;
  vector[n_miss_pre] pre_impute_raw;
  
  // Pre-value imputation model
  vector[K] mu_pre;
  vector<lower=0>[K] sigma_pre;
  cholesky_factor_corr[K] L_Omega_pre;
}
transformed parameters {
  vector[N] y;
  vector[N] pre;
  
  // Reconstruct pre (non-centered)
  {
    int idx_pre = 1;
    for (i in 1:N) {
      if (pre_miss[i] == 1) {
        pre[i] = pre_obs[i];
      } else {
        pre[i] = mu_pre[miss_pre_score[idx_pre]] 
                 + sigma_pre[miss_pre_score[idx_pre]] * pre_impute_raw[idx_pre];
        idx_pre += 1;
      }
    }
  }
  
  // Reconstruct y (non-centered)
  {
    int idx_y = 1;
    for (i in 1:N) {
      if (y_miss[i] == 1) {
        y[i] = y_obs[i];
      } else {
        if (miss_y_time[idx_y] == 1) {
          y[i] = alpha_mid[miss_y_score[idx_y], treat[i]] 
                 + beta[miss_y_score[idx_y]] * pre[i]
                 + sigma_mid[miss_y_score[idx_y]] * y_impute_raw[idx_y];
        } else {
          y[i] = alpha_post[miss_y_score[idx_y], treat[i]] 
                 + beta[miss_y_score[idx_y]] * pre[i]
                 + sigma_post[miss_y_score[idx_y]] * y_impute_raw[idx_y];
        }
        idx_y += 1;
      }
    }
  }
  
  // Build wide-format mean matrices
  matrix[N_part, K] mu_wide_mid = rep_matrix(0, N_part, K);
  matrix[N_part, K] mu_wide_post = rep_matrix(0, N_part, K);
  
  for (i in 1:N) {
    if (time[i] == 1) {
      mu_wide_mid[P[i], score[i]] = alpha_mid[score[i], treat[i]] + beta[score[i]] * pre[i];
    } else {
      mu_wide_post[P[i], score[i]] = alpha_post[score[i], treat[i]] + beta[score[i]] * pre[i];
    }
  }
}
model {
  // Priors
  to_vector(alpha_mid) ~ normal(0, 0.5);
  to_vector(alpha_post) ~ normal(0, 0.5);
  beta ~ normal(1, 0.2);
  sigma_mid ~ exponential(1);
  sigma_post ~ exponential(1);
  L_Omega_mid ~ lkj_corr_cholesky(4);
  L_Omega_post ~ lkj_corr_cholesky(4);
  
  mu_pre ~ normal(0, 1);
  sigma_pre ~ exponential(1);
  L_Omega_pre ~ lkj_corr_cholesky(2);
  
  // Non-centered priors
  y_impute_raw ~ std_normal();
  pre_impute_raw ~ std_normal();
  
  if (prior_only == 0) {
    // Pre-value model
    matrix[K, K] L_Sigma_pre = diag_pre_multiply(sigma_pre, L_Omega_pre);
    matrix[N_part, K] pre_wide = rep_matrix(0, N_part, K);
    for (i in 1:N) {
      pre_wide[P[i], score[i]] = pre[i];
    }
    for (p in 1:N_part) {
      pre_wide[p]' ~ multi_normal_cholesky(mu_pre, L_Sigma_pre);
    }
    
    // Mid-study model
    matrix[K, K] L_Sigma_mid = diag_pre_multiply(sigma_mid, L_Omega_mid);
    matrix[N_part, K] y_wide_mid = rep_matrix(0, N_part, K);
    for (i in 1:N) {
      if (time[i] == 1) y_wide_mid[P[i], score[i]] = y[i];
    }
    for (p in 1:N_part) {
      y_wide_mid[p]' ~ multi_normal_cholesky(mu_wide_mid[p]', L_Sigma_mid);
    }
    
    // Post-treatment model
    matrix[K, K] L_Sigma_post = diag_pre_multiply(sigma_post, L_Omega_post);
    matrix[N_part, K] y_wide_post = rep_matrix(0, N_part, K);
    for (i in 1:N) {
      if (time[i] == 2) y_wide_post[P[i], score[i]] = y[i];
    }
    for (p in 1:N_part) {
      y_wide_post[p]' ~ multi_normal_cholesky(mu_wide_post[p]', L_Sigma_post);
    }
  }
}
generated quantities {
  matrix[K, K] Omega_mid = multiply_lower_tri_self_transpose(L_Omega_mid);
  matrix[K, K] Omega_post = multiply_lower_tri_self_transpose(L_Omega_post);
  matrix[K, K] Omega_pre = multiply_lower_tri_self_transpose(L_Omega_pre);
  
  vector[N] y_rep;
  vector[N_part] log_lik_part;
  
  {
    matrix[K, K] L_Sigma_mid = diag_pre_multiply(sigma_mid, L_Omega_mid);
    matrix[K, K] L_Sigma_post = diag_pre_multiply(sigma_post, L_Omega_post);
    
    matrix[N_part, K] y_wide_mid = rep_matrix(0, N_part, K);
    matrix[N_part, K] y_wide_post = rep_matrix(0, N_part, K);
    for (i in 1:N) {
      if (time[i] == 1) {
        y_wide_mid[P[i], score[i]] = y[i];
      } else {
        y_wide_post[P[i], score[i]] = y[i];
      }
    }
    
    for (p in 1:N_part) {
      vector[K] y_part_mid = multi_normal_cholesky_rng(mu_wide_mid[p]', L_Sigma_mid);
      vector[K] y_part_post = multi_normal_cholesky_rng(mu_wide_post[p]', L_Sigma_post);
      
      log_lik_part[p] = multi_normal_cholesky_lpdf(y_wide_mid[p]' | mu_wide_mid[p]', L_Sigma_mid)
                      + multi_normal_cholesky_lpdf(y_wide_post[p]' | mu_wide_post[p]', L_Sigma_post);
      
      for (j in 1:obs_count[p]) {
        int i = obs_idx[p, j];
        if (time[i] == 1) {
          y_rep[i] = y_part_mid[score[i]];
        } else {
          y_rep[i] = y_part_post[score[i]];
        }
      }
    }
  }
  
  vector[N] y_expect;
  for (i in 1:N) {
    if (time[i] == 1) {
      y_expect[i] = mu_wide_mid[P[i], score[i]];
    } else {
      y_expect[i] = mu_wide_post[P[i], score[i]];
    }
  }
}

