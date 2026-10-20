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
  array[N] int<lower=0, upper=1> cc;          // Complete case indicator
  // Used in the log-likelihood to create comparable data sets for loo.

  int<lower=0, upper=1> prior_only;           // Prior only flag
}
transformed data {
  int n_miss_y = N - sum(y_miss);
  int M = 2 * K;

  // Diagnostic criteria data
  array[6] int diag_idx = {2, 6, 8, 10, 11, 12};
// set the right index for subscores
//
// analyte  type std_threshold_male std_threshold_female
// HDL         2            -1.10                -0.233
// TG          6             0.947                0.947
// dias        8            -0.0200              -0.0200
// glucose     10             1.88                 1.88
// sys         11             0.0638               0.0638
// waist       12            -0.432               -1.48

  vector[6] threshold_male   = [-1.1,   0.947, -0.02, 1.88, 0.0638, -0.432]';
  vector[6] threshold_female = [-0.233, 0.947, -0.02, 1.88, 0.0638, -1.48]';
  array[6] int above = {0, 1, 1, 1, 1, 1};   // 1 = above is abnormal

  // Participant-level sex and treat variables (constant per participant)
  array[N_P] real sex_p;
  array[N_P] int treat_p;
  for (i in 1:N) {
    sex_p[P[i]] = sex[i];
    treat_p[P[i]] = treat[i];
  }

  // Per-(P, k) pre-imputation bookkeeping. pre is participant-level (T1)
  // but each (P, k) appears in two rows (time=1, time=2) carrying the same
  // pre_obs/pre_miss. Index pre_impute by (P, k) so each missing pre value
  // is one parameter, not two.
  array[N_P, K] int pre_pk_miss;
  array[N_P, K] int pre_pk_idx;
  matrix[N_P, K] pre_obs_pk = rep_matrix(0.0, N_P, K);
  for (p in 1:N_P) for (k in 1:K) {
    pre_pk_miss[p, k] = 0;
    pre_pk_idx[p, k] = 0;
  }
  for (i in 1:N) {
    if (pre_miss[i] == 1) pre_obs_pk[P[i], type[i]] = pre_obs[i];
    else                  pre_pk_miss[P[i], type[i]] = 1;
  }
  int n_miss_pre_pk = 0;
  for (p in 1:N_P) for (k in 1:K) {
    if (pre_pk_miss[p, k] == 1) {
      n_miss_pre_pk += 1;
      pre_pk_idx[p, k] = n_miss_pre_pk;
    }
  }

  // Participant observation index for generated quantities
  array[N_P, M] int obs_idx;
  array[N_P] int obs_count;
  for (p in 1:N_P) obs_count[p] = 0;
  for (i in 1:N) {
    obs_count[P[i]] += 1;
    obs_idx[P[i], obs_count[P[i]]] = i;
  }

  // Per-participant mask of valid entries in the 2K stacked vector.
  // Entry j corresponds to time t=((j-1) div K)+1 and type k=((j-1) mod K)+1.
  // An entry is "valid" iff y is observed at (t, k) AND pre is observed for k.
  // Only valid entries contribute to log_lik_part in generated quantities.
  array[N_P, M] int valid_mask;
  array[N_P] int valid_count;
  array[N_P, M] int valid_idx;
  {
    array[N_P, K] int pre_observed;
    array[N_P, 2, K] int y_observed;
    for (p in 1:N_P) {
      for (k in 1:K) pre_observed[p, k] = 0;
      for (t in 1:2) for (k in 1:K) y_observed[p, t, k] = 0;
    }
    for (i in 1:N) {
      if (pre_miss[i] == 1) pre_observed[P[i], type[i]] = 1;
      if (y_miss[i] == 1)   y_observed[P[i], time[i], type[i]] = 1;
    }
    for (p in 1:N_P) {
      for (j in 1:M) { valid_mask[p, j] = 0; valid_idx[p, j] = 0; }
      for (t in 1:2) {
        for (k in 1:K) {
          int j = (t - 1) * K + k;
          if (y_observed[p, t, k] == 1 && pre_observed[p, k] == 1)
            valid_mask[p, j] = 1;
        }
      }
      int c = 0;
      for (j in 1:M) {
        if (valid_mask[p, j] == 1) {
          c += 1;
          valid_idx[p, c] = j;
        }
      }
      valid_count[p] = c;
    }
  }

  // Complete-case participant index for complete-case log-lik (loo)
  array[N_P] int cc_p = rep_array(0, N_P);
  for (i in 1:N) cc_p[P[i]] = cc[i];

  // Number of complete-case participants
  int N_CC = 0;
  for (p in 1:N_P) N_CC += cc_p[p];

  // cc_pos[p] = column of participant p in log_lik (0 if not complete case)
  // cc_id[j]  = participant id for column j
  array[N_P] int cc_pos = rep_array(0, N_P);
  array[N_CC] int cc_id;
  {
    int j = 1;
    for (p in 1:N_P)
      if (cc_p[p] == 1) { cc_pos[p] = j; cc_id[j] = p; j += 1; }
  }
}
parameters {
  // Time-specific intercepts: alpha[time][type, treat]
  array[2] matrix[K, 2] alpha;
  // Time-specific pre-value coefficients
  matrix[K, 2] beta;
  // Time-specific sex parameters
  matrix[K, 2] beta_sex;

  // Shared correlation; treatment-specific scales
  cholesky_factor_corr[M] L_Omega;
  // sigma[t][k, tx] = residual SD at time t, type k, treatment tx
  array[2] matrix<lower=0>[K, 2] sigma;

  // Missing data parameters
  vector[n_miss_y] y_impute;
  vector[n_miss_pre_pk] pre_impute_pk;

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

    for (i in 1:N) {
      if (y_miss[i] == 1) {
        y[i] = y_obs[i];
      } else {
        y[i] = y_impute[idx_y];
        idx_y += 1;
      }

      // pre is participant-level: index by (P, type) so duplicated rows
      // share the same imputed value.
      if (pre_pk_miss[P[i], type[i]] == 0) {
        pre[i] = pre_obs_pk[P[i], type[i]];
      } else {
        pre[i] = pre_impute_pk[pre_pk_idx[P[i], type[i]]];
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
                beta[type[i], time[i]] * pre[i] +
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
  to_vector(beta) ~ normal(1, 0.4);
  to_vector(beta_sex) ~ normal(0, 0.5);

  for (t in 1:2) to_vector(sigma[t]) ~ exponential(1);
  L_Omega ~ lkj_corr_cholesky(6);   // stronger LKJ regularises near-singular modes
  mu_pre ~ normal(0, 1);
  sigma_pre ~ exponential(1);
  L_Omega_pre ~ lkj_corr_cholesky(2);

  // Weak prior on imputed values
  y_impute ~ normal(0, 2);
  pre_impute_pk ~ normal(0, 2);

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

      // Marginal log-density over only the valid (y- and pre-observed) entries.
      // Participants with valid_count[p] == 0 contribute 0 -- filter in R.
      int v = valid_count[p];
      if (v == 0) {
        log_lik_part[p] = 0;
      } else {
        matrix[M, M] Sigma_p = multiply_lower_tri_self_transpose(L_Sigma_p);
        vector[v] y_sub;
        vector[v] mu_sub;
        matrix[v, v] Sigma_sub;
        for (a in 1:v) {
          int ia = valid_idx[p, a];
          y_sub[a]  = y_stacked[ia];
          mu_sub[a] = mu_stacked[ia];
          for (b in 1:v) {
            int ib = valid_idx[p, b];
            Sigma_sub[a, b] = Sigma_p[ia, ib];
          }
        }
        log_lik_part[p] = multi_normal_lpdf(y_sub | mu_sub, Sigma_sub);
      }

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

  // Complete-case log-lik for cross-model loo comparison. Named "log_lik" so
  // loo() picks it up by default. For a complete case every entry is valid, so
  // this is the full multivariate density of the participant vector.
  vector[N_CC] log_lik;
  for (p in 1:N_P)
    if (cc_p[p] == 1)
      log_lik[cc_pos[p]] = log_lik_part[p];

  // Compute diagnosis from continuous biomarker values
  // diag_idx order: {HDL=2, TG=6, dias=8, glucose=10, sys=11, waist=12}
  // BP (dias d=3 OR sys d=5) collapses into a single criterion,
  // so the maximum count is 5 (HDL, TG, BP, glucose, waist).
  for (p in 1:N_P) {
    vector[6] thresh;
    if (sex_p[p] < 0)
      thresh = threshold_male;
    else
      thresh = threshold_female;

    for (t in 1:3) {
      int count = 0;
      int bp_pass = 0;
      for (d in 1:6) {
        int k = diag_idx[d];
        real val;
        int pass;
        if (t == 1)
          val = pre_wide[p, k];
        else if (t == 2)
          val = y_wide_1[p, k];
        else
          val = y_wide_2[p, k];

        if (above[d] == 1)
          pass = (val >= thresh[d]);
        else
          pass = (val <= thresh[d]);

        if (d == 3 || d == 5) {
          if (pass == 1) bp_pass = 1;
        } else {
          count += pass;
        }
      }
      count += bp_pass;
      n_criteria[p, t] = count;
      diagnosed[p, t] = (count >= 3) ? 1.0 : 0.0;
    }
  }
}
