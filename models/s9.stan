/*
Model S9:
Extension of S7 (univariate normal per observation, missing-data imputation,
pre-value covariate, time-specific treatment effects) where the residual
SD is treatment-, time- and type-specific.

sigma is indexed by time x type x treat, so the model can represent
different within-group dispersion between treatments in the same type/time
cell.

Structure otherwise identical to S7 so the two can be compared via
log_lik_part / loo and by the same y_rep-based posterior predictive checks.

*/
data {
  int<lower=0> N;                           // Number of obs
  int<lower=0> N_P;                         // Number of participants
  int<lower=0> K;                           // Number of score types
  int<lower=0> N_TX;                        // Number of treatments

  vector[N] y_obs;                            // Data (incl missing w. placehol)
  vector[N] pre_obs;                          // Pre-int dat (incl miss)

  array[N] int<lower=0, upper=1> y_miss;      // Missing data indicator (miss=0)
  array[N] int<lower=0, upper=1> pre_miss;    // Miss dat ind for pre (miss=0)

  array[N] int<lower=1, upper=2> time;        // Time indicator
  array[N] int<lower=1, upper=K> type;        // Measurement type
  array[N] int<lower=1, upper=N_P> P;         // Participant
  array[N] int<lower=1, upper=2> treat;       // Treatment index
  array[N] int<lower=1, upper=2> sex_idx;     // Sex index (1=male, 2=female)
  array[N] int<lower=0, upper=1> cc;          // Complete case indicator
  // Used in the log-likelihood to create comparable data sets for loo.

  int<lower=0, upper=1> prior_only;         // Prior only flag
}

transformed data {
  int n_miss_y = N - sum(y_miss);             // Count the number of missing
  int n_miss_pre = N - sum(pre_miss);

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

  // Participant-level sex index (constant per participant)
  array[N_P] int sex_p_idx;
  for (i in 1:N) sex_p_idx[P[i]] = sex_idx[i];

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
  // An array per time-point, Type (k) treatments (2)
  array[2] matrix[K, N_TX] alpha;
  // Pre-score covariate parameter, separate for each type, time
  matrix[K, 2] beta;
  // Time-, type- and treatment-specific sigmas:
  // sigma[time][type, treat]
  array[2] matrix<lower=0>[K, N_TX] sigma;

  // Missing data parameters
  vector[n_miss_y] y_impute;
  vector[n_miss_pre] pre_impute;

}
transformed parameters {

  vector[N] y;
  vector[N] pre;


  // Local scope for the declared index vars here (idx_).
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

   // Linear model. Keeping it here to make mu and sigma available
   // in the gc block.
  vector[N] mu;
  vector[N] sig;

  for(i in 1:N) {
    mu[i] = alpha[time[i]][type[i], treat[i]] + beta[type[i],time[i]] * pre[i];
    sig[i]= sigma[time[i]][type[i], treat[i]];
  }

}

model {

  // Priors for alpha (across all time/type/treat)
  // Priors
  for (t in 1:2)  to_vector(alpha[t]) ~ normal(0, 1);

  for(ti in 1:2)
   for (ty in 1:K)
  beta[ty, ti] ~ normal(1, 0.4);

  for (ti in 1:2) to_vector(sigma[ti]) ~ exponential(1);

  // Weak prior on imputed values
  // if (prior_only == 1) {
  y_impute ~ normal(0, 2);
  pre_impute ~ normal(0, 2);
  // }


  if (prior_only == 0) {
    y ~ normal(mu, sig);
  }



}

generated quantities {
  vector[N] y_rep;
  vector[N_P] log_lik_part;

  // Per-participant, per-time diagnosis: 1 if >= 3 of 5 criteria met
  matrix[N_P, 3] diagnosed;
  // Per-participant, per-time count of criteria met (0-5)
  matrix[N_P, 3] n_criteria;

  // Y rep using the row indicies
  for(i in 1:N)
    y_rep[i] = normal_rng(mu[i], sig[i]);

  // Log lik using participant to correspond to MVNorm models.
  // Only paired observations (y AND pre observed) contribute.
  for (p in 1:N_P)
    log_lik_part[p] = 0;
  for (i in 1:N)
    if (y_miss[i] == 1 && pre_miss[i] == 1)
      log_lik_part[P[i]] += normal_lpdf(y[i] | mu[i], sig[i]);

  // Complete-case log-lik for cross-model loo comparison. Named "log_lik" so
  // loo() picks it up by default. For a complete case every observation is
  // valid, so this is the full density of the participant.
  vector[N_CC] log_lik;
  for (p in 1:N_P)
    if (cc_p[p] == 1)
      log_lik[cc_pos[p]] = log_lik_part[p];

  // Compute diagnosis from continuous biomarker values (incl imputed)
  {
    matrix[N_P, K] pre_wide  = rep_matrix(0.0, N_P, K);
    matrix[N_P, K] y_wide_1  = rep_matrix(0.0, N_P, K);
    matrix[N_P, K] y_wide_2  = rep_matrix(0.0, N_P, K);
    for (i in 1:N) {
      pre_wide[P[i], type[i]] = pre[i];
      if (time[i] == 1) y_wide_1[P[i], type[i]] = y[i];
      else              y_wide_2[P[i], type[i]] = y[i];
    }

    for (p in 1:N_P) {
      vector[6] thresh;
      if (sex_p_idx[p] == 1) thresh = threshold_male;
      else                   thresh = threshold_female;

      // diag_idx order: {HDL=2, TG=6, dias=8, glucose=10, sys=11, waist=12}
      // BP (dias d=3 OR sys d=5) collapses into a single criterion,
      // so the maximum count is 5 (HDL, TG, BP, glucose, waist).
      for (t in 1:3) {
        int count = 0;
        int bp_pass = 0;
        for (d in 1:6) {
          int k = diag_idx[d];
          real val;
          int pass;
          if (t == 1)      val = pre_wide[p, k];
          else if (t == 2) val = y_wide_1[p, k];
          else             val = y_wide_2[p, k];

          if (above[d] == 1) pass = (val >= thresh[d]);
          else               pass = (val <= thresh[d]);

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
}
