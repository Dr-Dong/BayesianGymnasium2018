library(rstan)
library(shinystan)
library(car)
library(mvtnorm)
library(rethinking)
library(MASS)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("../utilityFunctions.R")

#################################
# Plot of three different priors varying in shrinkage
par(mar=c(3,3.2,0.1,0.5))
curve(dnorm(x,0,0.5), from=-3, to=3, lty=3, lwd=1.5, las=1, ann=FALSE)
curve(dnorm(x,0,1), lty=1, lwd=1.5, add=TRUE, col="blue")
curve(dnorm(x,0,10), lty=2, lwd=1.5, add=TRUE)

mtext(expression(beta), side=1, line=2)
mtext("density", side=2, line=2.23)
text(1.5,0.7, "N(0,0.5)")
text(1.8,0.3, "N(0,1)")
text(0,0.1, "N(0,10)")

#################################

#################################
# Code from last class experimenting with in vs. out of sample
# deviance. Conduct simulations of regressions with more or less 
# than the "true" number of parameters. We then compare the 
# fit of the training dataset to the test dataset.
N <- 20
kseq <- 1:5
nSims <- 500

# If you have a mac set this
cores <- 4

dev <- sapply(kseq, function (k) {
  #  print(k);
  # If you have a mac use this
  r <- mcreplicate( nSims, sim.train.test(N=N, k=k), mc.cores=cores)
  # If you don't have a mac use this
  #  r <- replicate(nSims, sim.train.test( N=N, k=k ));
  c( mean(r[1,]), mean(r[2,]), sd(r[1,]), sd(r[2,]) )
} )

#################################

#################################
# Plot of results

par(mar=c(3,3.2,0.1,0.5))
plot(1:5, dev[1,], ylim=c(45,60), 
     xlim=c(1,5.3), las=1, ann=FALSE, pch=16, col="blue")

mtext(text = "# parameters (k)", side=1, line = 2, cex=1)
mtext(text = "deviance", side=2, line = 2.2, cex=1)

points((1:5), dev[2,])
lines(kseq, dev[1,] + 2*kseq, lty=2)
for (i in kseq) {
  lines(c(i,i)+0.1, c(dev[1,i], dev[2,i]))
  text(i+0.3, dev[2,i] - (dev[2,i]-dev[1,i])+2, 
       labels = eval(round(dev[2,i]-dev[1,i],1)))
}

#################################

#################################
# Read in milk dataset

milk <- read.csv("milk.csv")
cort <- milk$neocortex.perc/100
obs <- milk$kcal.per.g
mass <- log(milk$mass)

nObs <- nrow(milk)
# priors for later
bMu <- 0; bSD <- 1; sigmaSD <- 2.5

# Set up design matrix of response variables
xMat <- as.matrix(model.matrix(~cort + mass))
xMat[1:3,]

#################################

#################################
# Intercept only model
X1 <- as.matrix(xMat[,1]) # intercept-only design matrix
nVar <- ncol(X1)
# set up the data as a list:
d1 <- list(nObs=nObs, nVar=nVar, obs=obs, X=X1, bMu=bMu, 
           bSD=bSD, sigmaSD=2.5)

mod1 <- stan(file="15.multMod.stan", data=d1, iter=2000,
             chains=4, seed=867.5309, pars=c("log_lik"))

#################################

#################################
# Extract log-likelihood as a matrix

log_lik <- t(as.matrix(mod1, "log_lik"))

#################################

#################################
# Compute log pointwise predictive density

logSumExp <- function(x, n.iter) {
  out <- log(sum(exp(x))) - log(n.iter)
  return(out)
}

lppd <- apply(log_lik, 1, logSumExp, n.iter=ncol(log_lik))

#################################

#################################
# Calculate effective number of parameters

pWAIC <- apply(log_lik, 1, var)

#################################

#################################
# Calculate expected pointwise predictive density

(elpd <- sum(lppd) - sum(pWAIC))

# and calculate WAIC
-2 * elpd

# and calculate SE of WAIC
waic_vec <- -2 * (lppd - pWAIC)
sqrt(nObs * var(waic_vec))

#################################

#################################
# ######### intercept + neocortex
X2 <- as.matrix(xMat[,1:2])
nVar <- ncol(X2)
d2 <- list(nObs=nObs, nVar=nVar, obs=obs, X=X2, bMu=bMu, bSD=bSD, sigmaSD=2.5)

mod2 <- stan(file="15.multMod.stan", data=d2, iter=2000,
             chains=4, seed=867.5309, pars=c("log_lik"))

#################################

#################################
########## intercept + mass
X3 <- as.matrix(xMat[,c(1,3)])
nVar <- ncol(X3)
d3 <- list(nObs=nObs, nVar=nVar, obs=obs, X=X3, bMu=bMu, bSD=bSD, sigmaSD=2.5)

mod3 <- stan(file="15.multMod.stan", data=d3, iter=2000,
             chains=4, seed=867.5309, pars=c("log_lik"))

#################################

#################################
#X4 <- as.matrix(xMat)
nVar <- ncol(X4)
d4 <- list(nObs=nObs, nVar=nVar, obs=obs, X=X4, bMu=bMu, bSD=bSD, sigmaSD=2.5)

mod4 <- stan(file="15.multMod.stan", data=d4, iter=2000,
             chains=4, seed=867.5309, pars=c("log_lik"))

#################################

#################################
# Calculate WAIC using loo package

library(loo)
likM1 <- extract_log_lik(mod1) # extract log likelihood
(aicM1 <- waic(likM1))

# do the same for second model
likM2 <- extract_log_lik(mod2) # extract log likelihood
(aicM2 <- waic(likM2))

# Compare the difference in WAIC's by hand
dM1M2 <- aicM1$pointwise[,1] - aicM2$pointwise[,1]
(elpdDiff <- sum(dM1M2))
(seDiff <- sd(dM1M2) * sqrt(nObs))

#################################

#################################
# Compare WAICs using compare function in loo

compare(aicM1, aicM2)

# or use compare function from rethinking package

#detach("package:loo", unload=TRUE)
library(rethinking)
(waicTab <-  compare(mod1, mod2, mod3, mod4))

#################################