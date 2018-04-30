data {
  int<lower=1> nObs;
  int<lower=1> nPlants;
  vector[nObs] obs;
  int<lower=1> plant[nObs];
  int afternoon[nObs];
}

parameters {
  vector[nPlants] bPlant;
  vector[nPlants] aPlant;
  vector[2] gamma;
  vector<lower=0>[2] sigma;
  real<lower=0> tau;
  corr_matrix[2] Omega;
}

transformed parameters {
  vector[2] beta[nPlants];
  cov_matrix[2] Sigma;
  for(p in 1:nPlants) {
    beta[p,1] = aPlant[p];
    beta[p,2] = bPlant[p];  
  }

  Sigma = quad_form_diag(Omega,sigma);
}

model {
  vector[nObs] mu;
  Omega ~ lkj_corr(2);
  tau ~ cauchy(0 , 1);
  sigma ~ cauchy(0 , 1);
  gamma ~ normal(0 , 1);
  beta ~ multi_normal( gamma , Sigma);
  
  for ( n in 1:nObs ) 
    mu[n] = aPlant[plant[n]] + bPlant[plant[n]] * afternoon[n];
  
  obs ~ normal( mu , tau );
}
