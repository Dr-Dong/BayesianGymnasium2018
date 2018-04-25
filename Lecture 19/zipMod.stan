data {
  int<lower=0> nObs;
  int<lower=1> nVar;
  int<lower=0> obs[nObs];
  matrix[nObs, nVar] XT;
  matrix[nObs, nVar] XL;
  real<lower=0> thetaSD;
  real<lower=0> lambdaSD;
}

parameters {
  vector[nVar] betaT;
  vector[nVar] betaL;
}

transformed parameters {
  vector[nObs] logitTheta;
  vector[nObs] logLambda;
  
  logitTheta = XT * betaT;
  logLambda = XL * betaL;
}

model {
  betaT ~ normal(0, thetaSD);
  betaL ~ normal(0, lambdaSD);
  
  for(n in 1:nObs) {
    if(obs[n] == 0)
      target += log_sum_exp(
        bernoulli_logit_lpmf(1 | logitTheta[n]),
        bernoulli_logit_lpmf(0 | logitTheta[n])
          + poisson_log_lpmf(obs[n]|logLambda[n]));
      else
        target += bernoulli_logit_lpmf(0 | logitTheta[n])
          + poisson_log_lpmf(obs[n] | logLambda[n]);
  }
}
