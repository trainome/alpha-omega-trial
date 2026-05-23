/*
Model S13:
Extension of S9: replaces the univariate normal observation likelihood with a
univariate Student-t (student_t) where nu (degrees of freedom) is type-specific. 

sigma is indexed by sex x time x type x treat (as in S9). Under the t
likelihood, this parameterizes the *scale*; the marginal SD is
sqrt(nu/(nu-2)) * sigma (for nu > 2).

Structure otherwise identical to S9 so the two can be compared via
log_lik_part / loo and the same y_rep-based posterior predictive checks.

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

  int<lower=0, upper=1> prior_only;         // Prior only flag
}

transformed data {
  int n_miss_y = N - sum(y_miss);             // Count the number of missing
  int n_miss_pre = N - sum(pre_miss);

  // Diagnostic criteria data
  array[6] int diag_idx = {2, 4, 6, 8, 9, 10};
// set the right index for subscores
//
// analyte  type std_threshold_male std_threshold_female
// HDL         2            -1.10                -0.233
// TG          4             0.947                0.947
// dias        6            -0.0200              -0.0200
// glucose     8             1.88                 1.88
// sys         9             0.0638               0.0638
// waist      10            -0.432               -1.48

  vector[6] threshold_male   = [-1.1,   0.947, -0.02, 1.88, 0.0638, -0.432]';
  vector[6] threshold_female = [-0.233, 0.947, -0.02, 1.88, 0.0638, -1.48]';
  array[6] int above = {0, 1, 1, 1, 1, 1};   // 1 = above is abnormal

  // Participant-level sex index (constant per participant)
  array[N_P] int sex_p_idx;
  for (i in 1:N) sex_p_idx[P[i]] = sex_idx[i];
}

parameters {
  // An array per time-point, Type (k) treatments (2)
  array[2] matrix[K, N_TX] alpha;
  // Pre-score covariate parameter, separate for each type, time
  matrix[K, 2] beta;
  // Sex-, time-, type- and treatment-specific scale (interpreted as the
  // Student-t scale parameter, not the marginal SD):
  // sigma[sex_idx][time][type, treat]
  array[2, 2] matrix<lower=0>[K, N_TX] sigma;

  // Degrees of freedom for the Student-t observation likelihood, treatment-
  // specific
  vector<lower=1>[K] nu;

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
  vector[N] nu_k; // type-specific 
  
  for(i in 1:N) {
    mu[i] = alpha[time[i]][type[i], treat[i]] + beta[type[i],time[i]] * pre[i];
    sig[i] = sigma[sex_idx[i], time[i]][type[i], treat[i]];
    nu_k[i] = nu[type[i]]; 
    
  }

}

model {

  // Priors for alpha (across all time/type/treat)
  // Priors
  for (t in 1:2)  to_vector(alpha[t]) ~ normal(0, 1);

  for(ti in 1:2)
   for (ty in 1:K)
  beta[ty, ti] ~ normal(1, 0.4);

  for (s in 1:2) for (ti in 1:2) to_vector(sigma[s, ti]) ~ exponential(1);
  nu ~ gamma(2, 0.1);

  // Weak prior on imputed values
  // if (prior_only == 1) {
  y_impute ~ normal(0, 2);
  pre_impute ~ normal(0, 2);
  // }


  if (prior_only == 0) {
    y ~ student_t(nu_k, mu, sig);
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
    y_rep[i] = student_t_rng(nu_k[i], mu[i], sig[i]);

  // Log lik using participant to correspond to MVNorm models.
  // Only paired observations (y AND pre observed) contribute.
  for (p in 1:N_P)
    log_lik_part[p] = 0;
  for (i in 1:N)
    if (y_miss[i] == 1 && pre_miss[i] == 1)
      log_lik_part[P[i]] += student_t_lpdf(y[i] | nu_k[i], mu[i], sig[i]);

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

      // diag_idx order: {HDL=2, TG=4, dias=6, glucose=8, sys=9, waist=10}
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
