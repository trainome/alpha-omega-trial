/* 
Model S9:
A joint model for Metabolic Syndrome Criteria subscores.
-  We model each subscore separately with data imputation on both 
post-intervention time-points. 
- Pre-values are co-variate and score-specific treatment 
effect as the focal parameters
- A metabolic syndrome diagnostic criteria indicator is used in a logistic
model.

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
  
  int<lower=0, upper=1> prior_only;         // Prior only flag
}

transformed data {
  int n_miss_y = N - sum(y_miss);             // Count the number of missing
  int n_miss_pre = N - sum(pre_miss);
}

parameters {
  // An array per time-point, Type (k) treatments (2)
  array[2] matrix[K, N_TX] alpha;  
  // Pre-score covariate parameter, separate for each type, time
  matrix[K, 2] beta;                     
  // Time- and type-specific sigmas
  matrix<lower=0>[K,2] sigma;
  
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
    sig[i]= sigma[type[i],time[i]];
  }

}

model {

  // Priors for alpha (across all time/type/treat)
  // Priors
  for (t in 1:2)  to_vector(alpha[t]) ~ normal(0, 1); 
  
  for(ti in 1:2)
   for (ty in 1:K)
  beta[ty, ti] ~ normal(1, 0.4);
  
  for(ti in 1:2)
   for (ty in 1:K)
  sigma[ty,ti] ~ exponential(1);
  
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
  
  // Y rep using the row indicies
  for(i in 1:N)
    y_rep[i] = normal_rng(mu[i], sig[i]);
  
  // Log lik using participant to correspond to MVNorm models
  for (p in 1:N_P)
    log_lik_part[p] = 0;
  for (i in 1:N)
    log_lik_part[P[i]] += normal_lpdf(y[i] | mu[i], sig[i]);
}
