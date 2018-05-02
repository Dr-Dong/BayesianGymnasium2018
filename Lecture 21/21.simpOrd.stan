data{
  int<lower=0> nObs;
  int<lower=0> K;
  int<lower=1> obs[nObs];
  real<lower=0> alphaSD;
}

parameters{
  ordered[K-1] alpha;
}

model{
  vector[nObs] phi;
  
  alpha ~ normal(0 , alphaSD);
  for ( i in 1:nObs ) {
    phi[i] = 0;
    obs[i] ~ ordered_logistic( phi[i] , alpha );
  }
}
