library(rstan)
library(shinystan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


obs <- rep(c(1,0), times=c(7,3)) # our Bernoulli observations
nObs <- length(obs)              # number of observations
alpha <- 1                       # Prior for alpha
beta <- 1                        # Prior for beta
dat <- list(obs = obs, nObs=nObs, alpha=alpha, beta=beta)

mod1 <- stan(file="04.ex1Bernoulli.stan", #path to .stan file
             data=dat,
             iter=2000, # number of MCMC iterations
             chains=4,  # number of independent MCMC chains
             seed=3)    # se


# Traceplots 
traceplot(mod1, par="lp__")

# summary
print(mod1)

# Density plots 
stan_dens(mod1, par="p")

# Histogram plot
stan_hist(mod1, par="p")

# ShinyStan
launch_shinystan(mod1)

# Binomial data
nObs <- length(obs)
obs <- sum(obs)
alpha <- 1
beta <- 1
dat2 <-  list(nObs=nObs, obs = obs, alpha=alpha, beta=beta)

# binomial model
modBin1 <- stan(file="05.ex1Binomial.stan", #path to .stan file
                data=dat2,
                iter=2000, # number of MCMC iterations
                chains=4,  # number of independent MCMC chains
                seed=3)    # set the seed so run is repeatable

# bernoulli model results
print(mod1)  

# Binomial model results
print(modBin1)


####################### Uninformed flat prior
omega <- 0.5
kappa <- 2
dat11 <-  list(nObs=nObs, obs = obs, omega=omega, kappa=kappa)

modBin11 <- stan(file="05.binomialMode.stan", #path to .stan file
                 data=dat11,
                 iter=2000, # number of MCMC iterations
                 chains=4,  # number of independent MCMC chains
                 seed=3)    # set the seed so run is repeatable

## No likelihood Beta(1,1)
modNoLik11 <- stan(file="05.noLikBinomial.stan", #path to .stan file
                   data=dat11,
                   iter=2000, # number of MCMC iterations
                   chains=4,  # number of independent MCMC chains
                   seed=3)    # set the seed so run is repeatable

# plots
p11 <- as.matrix(modBin11, par="p")
pNoLik11 <- as.matrix(modNoLik11, par="p")
par(mfrow=c(1,2))
par(mar=c(4,3,0.1,0.5))
plot(density(pNoLik11), xlab="p", las=1, main="", ylab="")
abline(v=omega, lty=2, lwd=2.5, col="red")
plot(density(p11), xlab="p", las=1, main="", ylab="Density")
abline(v=omega, lty=2, lwd=2.5, col="red")



##################### Weakly informed
## Essentially Beta(2,2)
omega <- 0.5
kappa <- 4
dat22 <-  list(nObs=nObs, obs = obs, omega=omega, kappa=kappa)

## Beta(2,2) parameterized by mode
modBin22 <- stan(file="05.binomialMode.stan", #path to .stan file
                 data=dat22,
                 iter=2000, # number of MCMC iterations
                 chains=4,  # number of independent MCMC chains
                 seed=3)    # set the seed so run is repeatable

## No likelihood Beta(1,1)
modNoLik22 <- stan(file="05.noLikBinomial.stan", #path to .stan file
                   data=dat22,
                   iter=2000, # number of MCMC iterations
                   chains=4,  # number of independent MCMC chains
                   seed=3)    # set the seed so run is repeatable

p22 <- as.matrix(modBin22, par="p")
pNoLik22 <- as.matrix(modNoLik22, par="p")
par(mfrow=c(1,2))
par(mar=c(4,3,0.1,0.5))
plot(density(pNoLik22), xlab="p", las=1, main="", ylab="")
abline(v=omega, lty=2, lwd=2.5, col="red")
plot(density(p22), xlab="p", las=1, main="", ylab="Density")
abline(v=omega, lty=2, lwd=2.5, col="red")


################  Strongly informed

## Essentially Beta(10,10)
omega <- 0.5
kappa <- 20
dat1010 <-  list(nObs=nObs, obs = obs, omega=omega, kappa=kappa)

## Beta(2,2) parameterized by mode
modBin1010 <- stan(file="05.binomialMode.stan", #path to .stan file
                   data=dat1010,
                   iter=2000, # number of MCMC iterations
                   chains=4,  # number of independent MCMC chains
                   seed=3)    # set the seed so run is repeatable

## No likelihood Beta(1,1)
modNoLik1010 <- stan(file="05.noLikBinomial.stan", #path to .stan file
                     data=dat1010,
                     iter=2000, # number of MCMC iterations
                     chains=4,  # number of independent MCMC chains
                     seed=3)    # set the seed so run is repeatable

# plot
p1010 <- as.matrix(modBin1010, par="p")
pNoLik1010 <- as.matrix(modNoLik1010, par="p")
par(mfrow=c(1,2))
par(mar=c(4,3,0.1,0.5))
plot(density(pNoLik1010), xlab="p", las=1, main="", ylab="")
abline(v=omega, lty=2, lwd=2.5, col="red")

plot(density(p1010), xlab="p", las=1, main="", ylab="Density")
abline(v=omega, lty=2, lwd=2.5, col="red")



