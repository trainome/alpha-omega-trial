/* 
Analysis of fatty acid profiles in all available data. Imputation of post-treat
values are done in the model. Heterogenous sigma between treatments. 
*/
data {
  int<lower=0> N;
  int<lower=0> N_FA;
  int<lower=0> N_P;
  int<lower=0> N_TX;
  int<lower=0> N_TI;
  
  // Index vars
  array[N] int<lower=1, upper=N_P> P; 
  array[N] int<lower=1, upper=N_FA> FA;
  array[N] int<lower=1, upper=N_TX> TX;
  array[N] int<lower=1, upper=N_TI> TI;
  
  // Predictors
  vector[N] pre;                            // Pre-treatment 
  vector[N] y;                              // Outcome (obs values; miss:-999)
  array[N] int<lower=0, upper=1> obs;       // 1 = observed, 0 = missing
  
  // Prior flag
  int<lower=0, upper=1> prior_only;
}
transformed data {
  // Split indices into observed and missing
  int N_obs = sum(obs);
  int N_mis = N - N_obs;
  
  array[N_obs] int ii_obs;
  array[N_mis] int ii_mis;
  
  {
    int i_obs = 1;
    int i_mis = 1;
    for (n in 1:N) {
      if (obs[n] == 1) {
        ii_obs[i_obs] = n;
        i_obs += 1;
      } else {
        ii_mis[i_mis] = n;
        i_mis += 1;
      }
    }
  }
}
parameters {
  array[N_FA, N_TI, N_TX] real alpha;
  vector[N_FA] beta;
  // Heterogenous sigma
  array[N_FA, N_TX] real<lower=0> sigma;
  
  // Imputed missing outcomes
  vector[N_mis] y_mis;
}
transformed parameters {
  // Combine observed and imputed into a single vector
  vector[N] y_all;
  for (n in 1:N_obs)
    y_all[ii_obs[n]] = y[ii_obs[n]];
  for (n in 1:N_mis)
    y_all[ii_mis[n]] = y_mis[n];
}
model {
  // Linear model
  vector[N] mu;
  vector[N] sig;
  for (n in 1:N) {
    mu[n] = alpha[FA[n], TI[n], TX[n]] + beta[FA[n]] * pre[n];
    sig[n] = sigma[FA[n], TX[n]];
  }
  
  // Priors
  for (fa in 1:N_FA)
    for (ti in 1:N_TI)
      for (tx in 1:N_TX)
        alpha[fa, ti, tx] ~ normal(0, 1);
  beta  ~ normal(1, 1);
  
  for (fa in 1:N_FA)
      for (tx in 1:N_TX)
        sigma[fa, tx] ~ exponential(1);
  
  // Prior on missing values
  y_mis ~ normal(0, 2);
  
  // Likelihood — all data (observed + imputed) 
  if (prior_only == 0) {
    y_all ~ normal(mu, sig);
  }
}
generated quantities {
  array[N] real y_rep;
  vector[N_mis] y_imp; // Store posterior draws of imputed values
  vector[N_obs] log_lik; // log-likelihood on observed data
  
  if (prior_only == 0) {
    for (n in 1:N) {
      real mu_n = alpha[FA[n], TI[n], TX[n]] + beta[FA[n]] * pre[n];
      y_rep[n] = normal_rng(mu_n, sigma[FA[n], TX[n]]);
    }
    for (k in 1:N_obs) {                    
        int i = ii_obs[k];
        log_lik[k] = normal_lpdf(y[i] |
                       alpha[FA[i], TI[i], TX[i]] + beta[FA[i]] * pre[i],
                       sigma[FA[i], TX[i]]);
      }
  }
  
  // Copy imputed values for easy extraction
  y_imp = y_mis;
}
