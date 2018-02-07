data {
  int<lower=0> nObs;                // Total number of observations
  int<lower=0> obs;                 // scalar count of Ws
  real<lower=0> alpha;              // alpha and beta input as data
  real<lower=0> beta;               //  rather than hard-coded
}  

parameters {
  real<lower=0, upper=1> p; // prob. of water
}

model {
  p ~ beta(alpha, beta);            // prior on theta 
  obs ~ binomial(nObs, p);  // binomial likelihood function    
}

