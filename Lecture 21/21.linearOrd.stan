data{
  int<lower=0> nObs;
  int<lower=0> K;
  int<lower=1> nVar;
  int<lower=1> obs[nObs];
  matrix[nObs, nVar] X;
  real<lower=0> alphaSD;
  real<lower=0> betaSD;
}

parameters{
  ordered[K-1] alpha;
  vector[nVar] beta;
}

model{
  vector[nObs] phi;
  phi = X*beta;
  
  alpha ~ normal(0 , alphaSD);
  beta ~ normal(0, betaSD);
  for ( i in 1:nObs ) {
    obs[i] ~ ordered_logistic( phi[i] , alpha);
  }
}
