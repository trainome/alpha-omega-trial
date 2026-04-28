/*
Model S10:
- Extension of S8: joint model with missing data imputation, sex predictor,
  time-specific intercepts, and a single 2K-dimensional MVN with LKJ prior
  on the full correlation structure.
- Residual scales are treatment-specific: sigma is indexed by time x type x treat.
  The correlation structure L_Omega is shared across treatments; only the
  marginal scales differ.
- Metabolic syndrome diagnosis is derived in generated quantities from
  the continuous biomarker model (hard threshold: >= 3 of 5 criteria).
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

  // Diagnostic criteria data
  array[5] int diag_idx = {1, 3, 5, 7, 8};  // set the right index for subscores
  vector[5] threshold_male = [-1.1, 0.684, 0.0164, 1.88, -0.432]'; // Add thresholds
  vector[5] threshold_female = [-0.233, 0.684, 0.0164, 1.88, -1.48]'; // Add thresholds
  array[5] int above = {0, 1, 1, 1, 1};      // 1 = above is abnormal

  // Participant-level sex and treat variables (constant per participant)
  array[N_P] real sex_p;
  array[N_P] int treat_p;
  for (i in 1:N) {
    sex_p[P[i]] = sex[i];
    treat_p[P[i]] = treat[i];
  }

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

  // Shared correlation; treatment-specific scales
  cholesky_factor_corr[M] L_Omega;
  // sigma[t][k, tx] = residual SD at time t, type k, treatment tx
  array[2] matrix<lower=0>[K, 2] sigma;

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
  matrix[N_P, K] pre_wide = rep_matrix(0.0, N_P, K);
  matrix[N_P, K] y_wide_1 = rep_matrix(0.0, N_P, K);
  matrix[N_P, K] y_wide_2 = rep_matrix(0.0, N_P, K);
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

  // Build wide matrices
  for (i in 1:N) {
    pre_wide[P[i], type[i]] = pre[i];

    if (time[i] == 1)
      y_wide_1[P[i], type[i]] = y[i];
    else
      y_wide_2[P[i], type[i]] = y[i];

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
  // Priors for the MVN model
  for (t in 1:2) to_vector(alpha[t]) ~ normal(0, 1);
  beta ~ normal(1, 0.4);
  to_vector(beta_sex) ~ normal(0, 0.5);

  for (t in 1:2) to_vector(sigma[t]) ~ exponential(1);
  L_Omega ~ lkj_corr_cholesky(2);
  mu_pre ~ normal(0, 1);
  sigma_pre ~ exponential(1);
  L_Omega_pre ~ lkj_corr_cholesky(2);

  // Weak prior on imputed values
  y_impute ~ normal(0, 2);
  pre_impute ~ normal(0, 2);

  // Likelihood for the MVN model
  if (prior_only == 0) {
    // Pre-value imputation model
    matrix[K, K] L_Sigma_pre = diag_pre_multiply(sigma_pre, L_Omega_pre);

    array[N_P] vector[K] pre_vec;
    for (p in 1:N_P) pre_vec[p] = pre_wide[p]';
    pre_vec ~ multi_normal_cholesky(mu_pre, L_Sigma_pre);

    // Full post-int model: per-participant L_Sigma with treat-specific scales
    for (p in 1:N_P) {
      vector[M] sig_p;
      sig_p[1:K]     = sigma[1][, treat_p[p]];
      sig_p[(K+1):M] = sigma[2][, treat_p[p]];
      matrix[M, M] L_Sigma_p = diag_pre_multiply(sig_p, L_Omega);

      vector[M] y_p;
      vector[M] mu_p;
      y_p[1:K]      = y_wide_1[p]';
      y_p[(K+1):M]  = y_wide_2[p]';
      mu_p[1:K]     = mu_wide_1[p]';
      mu_p[(K+1):M] = mu_wide_2[p]';

      y_p ~ multi_normal_cholesky(mu_p, L_Sigma_p);
    }
  }
}
generated quantities {
  matrix[M, M] Omega = multiply_lower_tri_self_transpose(L_Omega);
  matrix[K, K] Omega_pre = multiply_lower_tri_self_transpose(L_Omega_pre);

  matrix[K, K] Omega_1   = Omega[1:K, 1:K];
  matrix[K, K] Omega_2   = Omega[(K+1):M, (K+1):M];
  matrix[K, K] Omega_cross = Omega[1:K, (K+1):M];

  vector[K] rho_time;
  for (k in 1:K)
    rho_time[k] = Omega[k, K + k];

  vector[N] y_rep;
  vector[N_P] log_lik_part;
  vector[N] y_expect;

  // Per-participant, per-time diagnosis: 1 if >= 3 of 5 criteria met
  matrix[N_P, 3] diagnosed;
  // Per-participant, per-time count of criteria met (0-5)
  matrix[N_P, 3] n_criteria;

  {
    for (p in 1:N_P) {
      vector[M] sig_p;
      sig_p[1:K]     = sigma[1][, treat_p[p]];
      sig_p[(K+1):M] = sigma[2][, treat_p[p]];
      matrix[M, M] L_Sigma_p = diag_pre_multiply(sig_p, L_Omega);

      vector[M] mu_stacked;
      mu_stacked[1:K]     = mu_wide_1[p]';
      mu_stacked[(K+1):M] = mu_wide_2[p]';

      vector[M] y_part = multi_normal_cholesky_rng(mu_stacked, L_Sigma_p);

      vector[M] y_stacked;
      y_stacked[1:K]     = y_wide_1[p]';
      y_stacked[(K+1):M] = y_wide_2[p]';
      log_lik_part[p] = multi_normal_cholesky_lpdf(y_stacked | mu_stacked, L_Sigma_p);

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

  // Compute diagnosis from continuous biomarker values
  for (p in 1:N_P) {
    vector[5] thresh;
    if (sex_p[p] < 0)
      thresh = threshold_male;
    else
      thresh = threshold_female;

    for (t in 1:3) {
      int count = 0;
      for (d in 1:5) {
        int k = diag_idx[d];
        real val;
        if (t == 1)
          val = pre_wide[p, k];
        else if (t == 2)
          val = y_wide_1[p, k];
        else
          val = y_wide_2[p, k];

        if (above[d] == 1)
          count += (val >= thresh[d]);
        else
          count += (val <= thresh[d]);
      }
      n_criteria[p, t] = count;
      diagnosed[p, t] = (count >= 3) ? 1.0 : 0.0;
    }
  }
}
