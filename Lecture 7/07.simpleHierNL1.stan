data {
  int<lower=0> nObs;        // Total number of observations
  int<lower=0> N;
  int<lower=0> obs;   // obs as scalar
  real<lower=2> kappa;      // concentration 
  real<lower=0> alpha_omega;    // priors on Omega
  real<lower=0> beta_omega;
}  

parameters {
  real<lower=0, upper=1> omega;     // overall prior mode
  real<lower=0, upper=1> theta;     // prob. of water
}

model {
  omega ~ beta(alpha_omega,beta_omega);
  { //  a & b are local parameters and not saved.
    real alpha;
    real beta;
    alpha = omega * (kappa - 2) +1;
    beta = (1 - omega) * (kappa - 2) + 1;
    theta ~ beta(alpha, beta);           
  }
  
  // obs ~ binomial(N, theta);    comment out likelihood
}

