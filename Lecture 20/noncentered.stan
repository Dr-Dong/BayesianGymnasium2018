data {
  int<lower=1> nObs;
  int<lower=1> nPar;
  int<lower=1> nPlants;
  vector[nObs] obs;
  int<lower=1> plant[nObs];
  matrix[nPlants, 1] hVar;
  matrix[nObs, nPar] afternoon;
}

parameters {
  matrix[nPar, nPlants] z;
  matrix[1, nPar] gamma;
  vector<lower=0>[2] sigma;
  real<lower=0> tau;
  cholesky_factor_corr[2] Omega;
}

transformed parameters {
  matrix[nPlants, nPar] beta;
  
  beta = hVar*gamma + (diag_pre_multiply(sigma,Omega) * z)';
}

model {
  vector[nObs] mu;

  Omega ~ lkj_corr_cholesky(2);
  tau ~ cauchy(0,1);
  sigma ~ cauchy(0,1);
  to_vector(gamma) ~ normal(0,1);
  to_vector(z) ~ normal(0,1);

  for(n in 1:nObs) 
  mu[n] = afternoon[n] * beta[plant[n]]';

  obs ~ normal(mu, tau);
}



