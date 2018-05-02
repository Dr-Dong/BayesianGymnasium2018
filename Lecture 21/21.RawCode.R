library(rstan)
library(shinystan)
library(car)
library(mvtnorm)
library(rethinking)
library(MASS)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("../utilityFunctions.R")
setwd("~/Dropbox/BayesClass/2018 Class/Lecture 21")

##################################
#Read in data and set up

dat <- read.csv("trolley.csv")
names(dat)
dim(dat)

obs <- dat$response     # response variable
action <- dat$action    # action or no
intent <- dat$intention # intention or no
contact <- dat$contact  # contact or no

##################################

##################################
# Plot results

par(mar=c(3,3.2,0.1,0.5))
plot(table(obs), type="h", xlim=c(1,7), ann=FALSE, las=1)
mtext("response", side=1, line=2)

##################################

##################################
# Plot cumulative data

echo=c(-(1:4))}
par(mar=c(3,3.2,0.1,0.5))
par(mfrow=c(1,3))
plot(table(obs), type="h", xlim=c(1,7), ann=FALSE, las=1)
mtext("response", side=1, line=2)

# discrete proportion of each response variable k
pK <- table(obs)/nrow(dat)

# convert to cumulative proportions
cumPK <- cumsum(pK)
plot(1:7, cumPK, type="b", las=1)
mtext("response", side=1, line=2)
mtext("Cum. proportion", side=2, line=2.2, cex=0.8)

##################################

##################################
# Logit function

logit <- function(x) {
  log(x/(1-x))
}

(lco <- round(logit(cumPK), 4))

##################################

##################################
# Add in log cumulative odds

par(mar=c(3,3.2,0.1,0.5))
par(mfrow=c(1,3))
plot(table(obs), type="h", xlim=c(1,7), ann=FALSE, las=1)
mtext("response", side=1, line=2)

plot(1:7, cumPK, type="b", las=1)
mtext("response", side=1, line=2)
mtext("Cum. proportion", side=2, line=2.2, cex=0.8)

plot(1:7, lco, type="b", las=1)
mtext("response", side=1, line=2)
mtext("log-cumulative odds", side=2, line=2.2, cex=0.8)

##################################

##################################
# Cumulative proportions

par(mar=c(3,3.2,0.1,0.5))

plot(cumPK, type="h", xlim=c(1,7), ylim=c(0,1), ann=FALSE, las=1)

par(new=TRUE)
plot(cumPK, type="b", axes=FALSE, ylim=c(0,1))
segments(x0=1+0.05, x1=1+0.05, y0=0, y1=cumPK[1], col="blue", lwd=2)

for(i in 2:7)
  segments(x0=i+0.05, x1=i+0.05, y0=cumPK[i-1], y1=cumPK[i],
           col="blue", lwd=2)
mtext("response", side=1, line=2)
mtext("cumulative proportion", side=2, line=2.2)

##################################

##################################
# Set up data for stan model and run 
# note stan model takes a LONG time

nObs <- nrow(dat)
K <- max(obs)

# function to set initial values so the ordering is correct
inits <- function() {
  list(alpha = sort(runif(K-1, -2,2)))
}

d <- list(nObs=nObs, K=K, obs=obs, alphaSD=5)
# m1 <- stan(file="21.simpOrd.stan", data=d, iter=2000, chains=4,
#            seed=867.5309, init=inits)

print(m1, pars="alpha")

##################################

##################################
# Convert our alphas to probabilities

alpha <- as.data.frame(m1, pars="alpha")
(probs <- logistic(colMeans(alpha)))

##################################

##################################
# Subtract 0.5 from each of our alphas and 
# compare results

# original estimates: 
(pk <- round(dordlogit(1:K, phi = 0, colMeans(alpha)), 3))
# average outcome
sum(pk*1:7)

# subtracting 0.5 from each:
(pk <- round(dordlogit(1:K, phi = 0.5, colMeans(alpha)), 3))
sum(pk*1:7)

##################################

##################################
# Ordered logistic regression with predictors
# model takes a LONG time to run
X <- model.matrix(~action*intent + intent*contact)[,-1]
nVar <- ncol(X)

inits2 <- function() {
  list(alpha = sort(runif(K-1, -2,2)), 
       beta = runif(nVar, -2,2))
}

d2 <- list(nObs=nObs, K=K, nVar=nVar, obs=obs, X=X, alphaSD=5,
           betaSD=1)
# m2 <- stan(file="21.linearOrd.stan", data=d2, iter=1000, chains=4,
#            seed=867.5309, init=inits2)

# print results

print(m2, pars=c("alpha", "beta"))

##################################

##################################
# Extract results

post <- as.data.frame(m2, pars=c("alpha", "beta"))
names(post) <- c("a1", "a2", "a3", "a4", "a5", "a6",
                 "bA", "bI", "bC", "bAI", "bIC")

##################################

##################################
# plot results of ordered logistic as counterfactual
# plots
par(mar=c(3,3.2,1.2,0.5))
par(mfrow=c(1,3))
xseq <-seq(0.9, 0.1, length=7)

# empty plot for when A=C=0
plot(0:1,0:1, type="n", xlim=c(0,1), ylim=c(0,1), xaxp=c(0,1,1), las=1)
mtext("intention", side=1, line=2)

A <- 0  # value for action
C <- 0  # value for contact
I <- 0:1 # value for interaction

for(i in 1900:2000) {
  pz <- post[i,]
  ak <- as.numeric(pz[1:6])
  phi <- pz$bA*A + pz$bI*I + pz$bC*C + pz$bAI*A*I + pz$bIC*I*C
  pk <- pordlogit(1:6, a=ak, phi=phi)
  
  for(k in 1:6) lines(I, pk[,k], col="#6495ed30")
}
abline(h=0:1, lty=2)
mtext(concat("action=", A, ", contact=",C))
text(x=xseq, y=c(0.03, 0.12, 0.2, 0.35, 0.53, 0.7, 0.9), 
     labels = 1:7, col="blue", cex=1.1)


# empty plot for A=1, C=0
plot(0:1,0:1, type="n", xlim=c(0,1), ylim=c(0,1), xaxp=c(0,1,1), las=1)
mtext("intention", side=1, line=2)

A <- 1  # value for action
C <- 0  # value for contact
I <- 0:1

for(i in 1900:2000) {
  pz <- post[i,]
  ak <- as.numeric(pz[1:6])
  phi <- pz$bA*A + pz$bI*I + pz$bC*C + pz$bAI*A*I + pz$bIC*I*C
  pk <- pordlogit(1:6, a=ak, phi=phi)
  
  for(k in 1:6) lines(I, pk[,k], col="#6495ed30")
}
abline(h=0:1, lty=2)
mtext(concat("action=", A, ", contact=",C))
text(x=xseq, y=c(0.1, 0.26, 0.36, 0.53, 0.69, 0.8, 0.93), 
     labels = 1:7, col="blue", cex=1.1)


# empty plot for A=0, C=1
plot(0:1,0:1, type="n", xlim=c(0,1), ylim=c(0,1), xaxp=c(0,1,1), las=1)
mtext("intention", side=1, line=2)

A <- 0  # value for action
C <- 1  # value for contact
I <- 0:1

for(i in 1900:2000) {
  pz <- post[i,]
  ak <- as.numeric(pz[1:6])
  phi <- pz$bA*A + pz$bI*I + pz$bC*C + pz$bAI*A*I + pz$bIC*I*C
  pk <- pordlogit(1:6, a=ak, phi=phi)
  
  for(k in 1:6) lines(I, pk[,k], col="#6495ed30")
}
abline(h=0:1, lty=2)
mtext(concat("action=", A, ", contact=",C))
text(x=xseq, y=c(0.13, 0.35, 0.45, 0.56, 0.69, 0.8, 0.93), 
     labels = 1:7, col="blue", cex=1.1)

##################################

##################################
#
##################################

##################################
#
##################################

##################################
#
##################################

##################################
#
##################################

##################################
#
##################################

##################################
#
##################################

##################################
#
##################################

##################################
#
##################################

##################################
#
##################################

##################################
#
##################################

##################################
#
##################################

##################################
#
##################################

##################################
#
##################################

