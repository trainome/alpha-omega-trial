/*
Model S18:
Univariate Student-t counterpart of S13. Extension of S7 (univariate
observation likelihood per row, missing-data imputation, pre-value covariate,
time-specific treatment effects and sex, as in S13) where the normal likelihood is replaced by a
univariate Student-t with nu fixed as data (nu = 4, as in S12/S13).

sigma is indexed by time x type only (as in S7 and S13, no treatment index).
Under the t likelihood, sigma is the *scale*; the marginal SD is
sqrt(nu/(nu-2)) * sigma (for nu > 2).

Compared with S11 (univariate t, time x type x treat scales) this halves the
scale parameters. Compared with S13 (MVT, full correlation, time x type
scales) it drops the among-type correlation, so S13 vs S18 isolates the value
of the correlation structure under the t likelihood.

Complete-case log_lik is exported for loo comparison against S7-S17.
Metabolic syndrome diagnosis is derived in generated quantities from the
continuous biomarker model (harmonised definition: >= 3 of the 5 criteria waist,
HDL, TG, BP, glucose; central obesity is not required).
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
  array[N] real<lower=-0.5, upper=0.5> sex;   // Sex (-0.5 = male, 0.5 = female)
  array[N] int<lower=0, upper=1> cc;          // Complete case indicator
  // Used in the log-likelihood to create comparable data sets for loo.

  int<lower=1, upper=30> nu;                  // nu parameter for the mvt-dist

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
  // Time-specific sex parameters (as in S13)
  matrix[K, 2] beta_sex;
  // Time- and type-specific scale (interpreted as the Student-t scale
  // parameter, not the marginal SD), shared across treatments as in S7/S13:
  // sigma[type, time]
  matrix<lower=0>[K, 2] sigma;

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
    mu[i] = alpha[time[i]][type[i], treat[i]] + beta[type[i],time[i]] * pre[i] +
            beta_sex[type[i], time[i]] * sex[i];
    sig[i] = sigma[type[i], time[i]];

  }

}

model {

  // Priors for alpha (across all time/type/treat)
  // Priors
  for (t in 1:2)  to_vector(alpha[t]) ~ normal(0, 1);

  for(ti in 1:2)
   for (ty in 1:K)
  beta[ty, ti] ~ normal(1, 0.4);
  to_vector(beta_sex) ~ normal(0, 0.5);

  to_vector(sigma) ~ exponential(1);


  // Weak prior on imputed values
  // if (prior_only == 1) {
  y_impute ~ normal(0, 2);
  pre_impute ~ normal(0, 2);
  // }


  if (prior_only == 0) {
    y ~ student_t(nu, mu, sig);
  }



}

generated quantities {
  vector[N] y_rep;
  vector[N_P] log_lik_part;

  // Per-participant, per-time diagnosis: 1 if >= 3 of the 5 criteria
  // (waist, HDL, TG, BP, glucose) are met (harmonised definition)
  matrix[N_P, 3] diagnosed;
  // Per-participant, per-time count of criteria met (0-5)
  matrix[N_P, 3] n_criteria;

  // Y rep using the row indicies
  for(i in 1:N)
    y_rep[i] = student_t_rng(nu, mu[i], sig[i]);

  // Log lik using participant to correspond to MVNorm models.
  // Only paired observations (y AND pre observed) contribute.
  for (p in 1:N_P)
    log_lik_part[p] = 0;
  for (i in 1:N)
    if (y_miss[i] == 1 && pre_miss[i] == 1)
      log_lik_part[P[i]] += student_t_lpdf(y[i] | nu, mu[i], sig[i]);

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

      // diag_idx order (d): {HDL=2, TG=6, dias=8, glucose=10, sys=11, waist=12}
      // BP (dias d=3 OR sys d=5) collapses into a single criterion,
      // so the maximum count is 5 (HDL, TG, BP, glucose, waist).
      for (t in 1:3) {
        int count = 0;
        int bp_pass = 0;
        int waist_pass = 0;
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
          } else if (d == 6) {
            waist_pass = pass;
          } else {
            count += pass;
          }
        }
        count += bp_pass;
        // n_criteria keeps the total count (0-5) incl. waist for reference
        n_criteria[p, t] = count + waist_pass;
        // Harmonised definition: any 3 of the 5 criteria (waist not required)
        diagnosed[p, t] = (n_criteria[p, t] >= 3) ? 1.0 : 0.0;
      }
    }
  }
}
