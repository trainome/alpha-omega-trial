/* 
Model 1. 
Drop-out analysis using baseline characteritics as predictor variables.

The model does outcome-wise prediction of droup-out per treatment group
while controlling for sex.

The parameter(s) of interest is beta_x and beta_xtx which 
gives os the log-odds of drop-out per treatment group when 
increasing the predictor by 1 standard deviation.  

*/


// 
data {
  int<lower=0> N;
  int<lower=0> N_A;
  int<lower=0> N_P;
  int<lower=0> N_TX;
  int<lower=0> N_S;
  
  // Index vars
  array[N] int<lower=1, upper=N_P> P; 
  array[N] int<lower=1, upper=N_A> A;

  
  // Predictors
  array[N] real S;                               // Sex (-0.5/0.5 coding)
  array[N] real x;                               // Outcome pred
  array[N] real TX;                              // Treatment .5 coding
  
  array[N] int<lower=0, upper=1> do4;            // Drop-out to T4
  array[N] int<lower=0, upper=1> do2;            // T2   
  
  // Prior flag
  int<lower=0, upper=1> prior_only;
  
}

parameters {
  vector[N_A] alpha4;
  vector[N_A] beta4_x;
  vector[N_A] beta4_s;
  vector[N_A] beta4_tx;
  vector[N_A] beta4_xtx;
  
  vector[N_A] alpha2;
  vector[N_A] beta2_x;
  vector[N_A] beta2_s;
  vector[N_A] beta2_tx;
  vector[N_A] beta2_xtx;
  
  
}

model {
  
  // Time-point 4
  alpha4 ~ normal(0, 0.5);
  beta4_x ~ normal(0, 0.5);
  beta4_s ~ normal(0, 0.5);
  beta4_tx ~ normal(0, 0.5);
  beta4_xtx ~ normal(0, 0.5);
  
  // Time-point 2
  alpha2 ~ normal(0, 0.5);
  beta2_x ~ normal(0, 0.5);
  beta2_s ~ normal(0, 0.5);
  beta2_tx ~ normal(0, 0.5);
  beta2_xtx ~ normal(0, 0.5);
  
  

  
  
  if (prior_only == 0) {
      for (n in 1:N) {
        do4[n] ~ bernoulli_logit(
          alpha4[A[n]] + 
          beta4_x[A[n]] * x[n] + 
          beta4_s[A[n]] * S[n] + 
          beta4_tx[A[n]] * TX[n] + 
          beta4_xtx[A[n]] * x[n] * TX[n]);
          
        do2[n] ~ bernoulli_logit(
          alpha2[A[n]] + 
          beta2_x[A[n]] * x[n] + 
          beta2_s[A[n]] * S[n] + 
          beta2_tx[A[n]] * TX[n] + 
          beta2_xtx[A[n]] * x[n] * TX[n]);
          
          
    }
  }
  

  
}

