#setwd("~/Dropbox/BayesClass/2018 Class/Lecture 18")

library(rstan)
library(shinystan)
library(car)
library(mvtnorm)
library(rethinking)
library(MASS)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("../utilityFunctions.R")

#############################
## Simulate binomial sampling
set.seed(3)
trillium <- rbinom(2000, 1000, 1/1000)
c(mean(trillium), var(trillium))

#############################


#############################
# read in data

kline <- read.csv("kline.csv")
head(kline)

#############################


#############################
## sort dataframe by pop size.
dat <- kline[order(kline$pop),]
obs <- dat$tools
contact <- ifelse(dat$contact=="low", 0, 1)
pop <- log(dat$pop)
nObs <- nrow(dat)
# make design matrix
predMat <- model.matrix(~pop*contact)

# make design matrix for new simulated data
newPop <- rep(seq(6, 13, length=30), 2)
newCont <- rep(0:1, each=30)
newMat <- model.matrix(~newPop*newCont)
nNew <- nrow(newMat)

#############################


#############################
##Specify design matrices
X <- predMat
newX <- newMat

# declare data list
d1 <- list(nObs=nObs, nVar=ncol(X), nNew=nNew, nVarN=ncol(newX),
           obs=obs, X=X, newX=newX, b0SD=10, bMu=0, bSD=1)

# run model
m1 <- stan(file="poissonMod.stan", data=d1, iter=2000, chains=4,
           seed=867.5309, pars="lambda", include=FALSE)

# extract betas and print result
beta <- as.matrix(m1, "beta")
print(m1,"beta")

#############################


#############################
#plot results
par(mar=c(3,3,0.1,0.5))
plotBeta <- plot(m1, pars="beta")
plotBeta + theme(text=element_text(family="ArialMT"))

#############################


#############################
#Calculate lambdas
lamLo <- exp(beta[,1] + beta[,2]*8)
lamHi <- exp(beta[,1] + beta[,3] + (beta[,2] + beta[,4])*8)

# compare difference
diff <- lamHi - lamLo
round(sum(diff>0)/length(diff),3)

#############################


#############################
#paired plot
pairs(m1, pars="beta")

#############################


#############################
#1. A model without an interaction:

X <- predMat[,1:3]
newX <- newMat[,1:3]

d2 <- list(nObs=nObs, nVar=ncol(X), nNew=nNew, nVarN=ncol(newX),
           obs=obs, X=X, newX=newX, b0SD=10, bMu=0, bSD=1)

m2 <- stan(file="poissonMod.stan", data=d2, iter=2000, chains=4,
           seed=867.5309, pars="lambda", include=FALSE)

beta1 <- as.matrix(m1, "beta")

#############################


#############################
#2. A model with just population size

X <- predMat[,1:2]

newX <- newMat[,1:2]

d3 <- list(nObs=nObs, nVar=ncol(X), nNew=nNew, nVarN=ncol(newX),
           obs=obs, X=X, newX=newX, b0SD=10, bMu=0, bSD=1)

m3 <- stan(file="poissonMod.stan", data=d3, iter=2000, chains=4,
           seed=867.5309, pars="lambda", include=FALSE)

#############################


#############################
#3. A model with just contact rate

X <- predMat[,c(1,3)]
newX <- newMat[,c(1,3)]

d4 <- list(nObs=nObs, nVar=ncol(X), nNew=nNew, nVarN=ncol(newX),
           obs=obs, X=X, newX=newX, b0SD=10, bMu=0, bSD=1)

m4 <- stan(file="poissonMod.stan", data=d4, iter=2000, chains=4,
           seed=867.5309, pars="lambda", include=FALSE)

#############################


#############################
#Compare models

compare(m1, m2, m3, m4)

#############################


#############################
# Plot results side by side

par(mar=c(3,3.2,0.1,0.5))
par(mfrow=c(1,2))

# interaction model
lam <- as.matrix(m1, pars="newLam")
avLam <- exp(colMeans(lam))
hdiLam <- exp(apply(lam, 2, HDI, credMass=0.95))

avLo <- avLam[newCont==0]
avHi <- avLam[newCont==1]

hdiLo <- hdiLam[, newCont==0]
hdiHi <- hdiLam[, newCont==1]

x <- newPop[1:30]
plot(x, avLo, type="n", ann=FALSE, las=1)
mtext("log Pop", side=1, line=2)
mtext("Total tools", side=2, line=2.2)

polygon(c(x, rev(x)), c(hdiHi[1, ], rev(hdiHi[2,])), 
        col="#88CCEE50")
lines(x,avHi, lwd=2, col="#88CCEE")

polygon(c(x, rev(x)), c(hdiLo[1, ], rev(hdiLo[2,])), 
        col="#50505080")
lines(x,avLo, lwd=2, lty=2)


lam <- as.matrix(m2, pars="newLam")
avLam <- exp(colMeans(lam))
hdiLam <- exp(apply(lam, 2, HDI, credMass=0.95))

avLo <- avLam[newCont==0]
avHi <- avLam[newCont==1]

hdiLo <- hdiLam[, newCont==0]
hdiHi <- hdiLam[, newCont==1]

plot(x, avLo, type="n", ann=FALSE)
mtext("log Pop", side=1, line=2)

polygon(c(x, rev(x)), c(hdiHi[1, ], rev(hdiHi[2,])), 
        col="#88CCEE50")
lines(x,avHi, lwd=2, col="#88CCEE")

polygon(c(x, rev(x)), c(hdiLo[1, ], rev(hdiLo[2,])), 
        col="#50505080")
lines(x,avLo, lwd=2, lty=2)

#############################


#############################
# Centered model

cX <- model.matrix(~scale(pop, scale = FALSE) * contact)
cnewX <- model.matrix(~scale(newPop, scale=FALSE) * newCont)

dc <- list(nObs=nObs, nVar=ncol(cX), nNew=nNew, nVarN=ncol(cnewX),
           obs=obs, X=cX, newX=cnewX, b0SD=10, bMu=0, bSD=1)

mc <- stan(file="poissonMod.stan", data=dc, iter=2000, chains=4,
           seed=867.5309, pars="lambda", include=FALSE)

cBeta <- as.matrix(mc, "beta")
print(mc, "beta")

# paired plot
pairs(mc, pars="beta")

#############################


#############################
# Simulate dataset sampling for 30 days
yD <- rpois(30, 5)

#############################


#############################
# Simulate dataset sampling for a week

yW <- rpois(4, 3*7)

#############################


#############################
# Set up data

obs <- c(yD, yW)

expose <- c(rep(1,30), rep(7, 4))

state <- c(rep(0,30), rep(1, 4))

#############################


#############################
# Run model

dex <- list(nObs=length(obs), expose=expose, state=state,
            obs=obs)

mEx <- stan(file="exposeMod.stan", data=dex, iter=2000, chains=4,
            seed=867.5309, pars="lambda", include=FALSE)

pars <- as.data.frame(mEx, c("alpha", "beta"))

GA <- quantile(exp(pars$alpha), c(0.025, 0.5, 0.975))
TN <- quantile(exp(rowSums(pars)), c(0.025, 0.5, 0.975))

# look at results
rbind(GA,TN)

#############################

