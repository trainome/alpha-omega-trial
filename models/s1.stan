/* 
A preliminary joint model for Metabolic Syndrome Criteria subscores.
The model models each subscore separately on complete cases with pre-values 
as co-variate and score-specific treatment effect as the focal parameters

*/
data {
  int<lower=0> N;                           // Number of obs
  int<lower=0> N_P;                         // Number of participants 
  int<lower=0> K;                           // Number of score types
  int<lower=0> N_TX;                        // Number of treatments
  
  vector[N] y;                              // Outcome (post-treat scores)
  vector[N] pre;                            // Pre-treatment scores
  
  array[N] int<lower=1, upper=K> type;      // Score type 
  array[N] int<lower=1, upper=N_P> P;       // Participant ID
  array[N] int<lower=1, upper=N_TX> treat;  // 1 = control, 2 = treatment
  
  int<lower=0, upper=1> prior_only;         // Prior only flag
}

parameters {
  array[K, N_TX] real alpha;          // Intercept (mean for treat*type)
  vector[K] beta;                     // Pre-score covariate parameter
  vector<lower=0>[K] sigma;
}
transformed parameters {
   // Linear model. Keeping it here to make mu and sigma available 
   // in the gc block.
  vector[N] mu;
  vector[N] sig;
  
  for(i in 1:N) {
    mu[i] = alpha[type[i], treat[i]] + beta[type[i]] * pre[i]; 
    sig[i]= sigma[type[i]];
  }

}

model {

  // Priors
  for (ty in 1:K)
   for (tx in 1:N_TX)
      alpha[ty, tx] ~ normal(0, 1);  
  beta ~ normal(1, 0.5);
  sigma ~ exponential(1);
  
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
  
  // Log lik using participant to correspond to 
  // later models.
  for (p in 1:N_P)
    log_lik_part[p] = 0;
  for (i in 1:N)
    log_lik_part[P[i]] += normal_lpdf(y[i] | mu[i], sig[i]);
}
