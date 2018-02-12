data {
  int<lower=0> nObs;  // Total number of observations
  int<lower=0> N[nObs];
  int<lower=0> obs[nObs];   // obs as scalar
  real<lower=0, upper=1> omega;  // mode as input data
  real<lower=2> kappa;  // concentration 
}  

parameters {
  real<lower=0, upper=1> p;     // prob. of water
}

transformed parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  
  alpha = omega * (kappa - 2) + 1;
  beta = (1 - omega) * (kappa - 2) + 1;
}

model {
  p ~ beta(alpha, beta);                // prior on theta
  obs ~ binomial(N, p);    
}

