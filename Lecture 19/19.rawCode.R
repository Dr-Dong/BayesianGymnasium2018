library(rstan)
library(shinystan)
library(car)
library(mvtnorm)
library(rethinking)
library(MASS)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("../utilityFunctions.R")

#####################################
# Read in data

dat <- read.csv("plantCount.csv")
names(dat)

#####################################

#####################################
# Set up design matrices

XT <- model.matrix(~moist, data=dat)
XL <- model.matrix(~nit, data=dat)
nObs <- nrow(dat)
obs <- dat$count

# Create data list
d <- list(nObs=nObs, nVar=2, obs=obs, XT=XT, XL=XL, thetaSD=1, lambdaSD=1)

# Run model
m1 <- stan(file="zipMod.stan", data=d, iter=2000, chains=4, seed=867.5309)
#####################################

#####################################
# Look at parameters for probability of true absence

print(m1, pars="betaT")

# Convert to natural scale

logistic(1.35)
logistic(1.35-0.23)

#####################################

#####################################
# Look at Poisson part of model 

print(m1, pars="betaL")

# Interpret on natural scale

exp(0.75)
exp(0.04)
exp(0.75+0.04)

#####################################

#####################################
# Plot results

par(mar=c(3,3.2,0.1,0.5))
par(mfrow=c(1,2))

# Plot theta
theta <- logistic(as.matrix(m1,"logitTheta"))
avTheta <- 1-colMeans(theta)
hdiTheta <- 1 - apply(theta, 2, HDI, credMass=0.95)

x <-  dat$moist
y <- ifelse(dat$count == 0, 0, 1)
plot(x, y,  type="p", ann=FALSE, las=1, pch=16)
mtext("moisture", side=1, line=2)
mtext("Pr(presence)", side=2, line=2.2)

polygon(c(x, rev(x)), c(hdiTheta[1, ], rev(hdiTheta[2,])),
        col="#50505050")
lines(x,avTheta, lwd=2)

# plot lambda
lambda <- exp(as.matrix(m1, "logLambda"))
avLam <- colMeans(lambda)
hdiLam <- apply(lambda, 2, HDI, credMass=0.95)

x <-  dat$nit
y <- dat$count
plot(x, y,  type="p", ann=FALSE, las=1, pch=16)
mtext("Nitrogen", side=1, line=2)
mtext("Abundance", side=2, line=2.2)

polygon(c(x, rev(x)), c(hdiLam[1, ], rev(hdiLam[2,])),
        col="#50505050")
lines(x,avLam, lwd=2)
#####################################

#####################################
#
#####################################

#####################################
#
#####################################

#####################################
#
#####################################

#####################################
#
#####################################

#####################################
#
#####################################

#####################################
#
#####################################

#####################################
#
#####################################

#####################################
#
#####################################

#####################################
#
#####################################

#####################################
#
#####################################

#####################################
#
#####################################

#####################################
#
#####################################

#####################################
#
#####################################

#####################################
#
#####################################

#####################################
#
#####################################

