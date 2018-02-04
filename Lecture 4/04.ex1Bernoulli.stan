data {
  int<lower=0> nObs;                // Total number of observations
  int<lower=0, upper=1> obs[nObs];  // 1D array of observations
  real<lower=0> alpha;
  real<lower=0> beta;
} 

parameters {
  real<lower=0, upper=1> p;     // prob. of water
}

model {
  p ~ beta(alpha, beta);                //prior on p
  for(n in 1:nObs) {        
    obs[n] ~ bernoulli(p);      // bernoulli likelihood function    
  }
}

