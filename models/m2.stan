/* 
Complete-case analysis of fatty acid composition.

The model uses only complete cases (T1 to T4) to establish 
an ANCOVA-style model for each fatty acid.

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
  vector[N] y;                              // Outcome

  // Prior flag
  int<lower=0, upper=1> prior_only;
}

parameters {
  array[N_FA, N_TI, N_TX] real alpha;
  vector[N_FA] beta;
  vector<lower=0>[N_FA] sigma;
}

transformed parameters {

}

model {


 // Linear model
 vector[N] mu;
 vector[N] sig;
 for (n in 1:N) {
   mu[n] = alpha[FA[n], TI[n], TX[n]] + beta[FA[n]] * pre[n];
   sig[n]= sigma[FA[n]];
 }

  // Priors
  for (fa in 1:N_FA)
   for (ti in 1:N_TI)
    for (tx in 1:N_TX)
      alpha[fa, ti, tx] ~ normal(0, 1);
  beta  ~ normal(1, 1);
  sigma ~ exponential(1);



  // Likelihood
  if (prior_only == 0) {
     y ~ normal(mu, sig);
  }

 
}

generated quantities {
  
  array[N] real y_rep;
  if (prior_only == 0) {
    for (n in 1:N) {
      real mu_n = alpha[FA[n], TI[n], TX[n]] + beta[FA[n]] * pre[n];
      y_rep[n] = normal_rng(mu_n, sigma[FA[n]]);
    }
  }
}
