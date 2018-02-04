library(shinystan)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("YOUR DIRECTORY")

obs <- rep(c(1,0), times=c(7,3)) # our Bernoulli observations
nObs <- length(obs)              # number of observations
alpha <- 1                       # Prior for alpha
beta <- 1                        # Prior for beta
dat <- list(obs = obs, nObs=nObs, alpha=alpha, beta=beta)

mod1 <- stan(file="04.ex1Bernoulli.stan", #path to .stan file
             data=dat,
             iter=2000, # number of MCMC iterations
             chains=4,  # number of independent MCMC chains
             seed=3)    # set the seed so run is repeatable

traceplot(mod1, par="p")

traceplot(mod1, par="lp__")

print(mod1)

print(mod1, par="p")

stan_dens(mod1, par="p")

stan_dens(mod1, par="p")