library(rstan)
library(shinystan)
library(car)
library(mvtnorm)
library(rethinking)
library(MASS)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("../utilityFunctions.R")
setwd("~/Dropbox/BayesClass/2018 Class/Lecture 17")
#########################
#read in data

clutch <- read.csv("ClutchData.csv")
clutch[1:5,]
o <- order(clutch$Clutch)
clutch <- clutch[o,]
#########################

#########################
#Set up data and run 4 models

obs <- clutch$Mort
N <- clutch$Clutch
nObs <- length(obs)
treat <- clutch$Treat
predMat <- as.matrix(model.matrix(~Clutch*Treat, data=clutch))
bMu <- 0; bSD<-5

#########################

#########################
#1. Linear regression with clutch size as a predictor.

X <- predMat[,1:2]
d1 <- list(nObs=nObs, nVar=ncol(X), obs=obs, N=N, X=X,
           bMu=bMu, bSD=bSD)
m1 <- stan(file="logitMod.stan", data=d1, iter=2500, chains=4,
           seed=867.5309)

#########################

#########################
# 2. Linear regression with predator access as a predictor.

X <- predMat[,c(1,3)]
d2 <- list(nObs=nObs, nVar=ncol(X), obs=obs, N=N, X=X,
           bMu=bMu, bSD=bSD)
m2 <- stan(file="logitMod.stan", data=d2, iter=2500, chains=4,
           seed=867.5309)

#########################

#########################
# 3. Linear regression with clutch size and predator access as  predictors.
```{r, message=FALSE, warning=FALSE, cache=TRUE, verbose=FALSE}
X <- predMat[,c(1:3)]
d3 <- list(nObs=nObs, nVar=ncol(X), obs=obs, N=N, X=X,
           bMu=bMu, bSD=bSD)
m3 <- stan(file="logitMod.stan", data=d3, iter=2500, chains=4,
           seed=867.5309)

#########################

#########################
# 4. Linear regression with an interaction between clutch size and predator access.

X <- predMat
d4 <- list(nObs=nObs, nVar=ncol(X), obs=obs, N=N, X=X,
           bMu=bMu, bSD=bSD)
m4 <- stan(file="logitMod.stan", data=d4, iter=2500, chains=4,
           seed=867.5309)

#########################

#########################
# Compare results using WAIC
(comp <- compare(m1,m2,m3,m4))

#########################

#########################
# Model averaging

wts <- round(comp@output$weight[1:2],2)
betaM3 <- as.matrix(m3, pars="beta")
betaM3 <- cbind(betaM3, rep(0, nrow(betaM3)))
betaM4 <- as.matrix(m4, pars="beta")

avgBeta <- betaM4*wts[1] + betaM3*wts[2]
Mean <- colMeans(avgBeta)
StDev <- apply(avgBeta, 2, sd)
hdi <- apply(avgBeta, 2, HDI, credMass = 0.95)
x <- as.matrix(round(cbind(Mean, StDev, t(hdi)),3))
colnames(x) <- c("Mean", "SD", "0.025","0.975")
x

#########################

#########################
# Look at results 

noPred <- avgBeta[,2] + avgBeta[,4]*0
pred <- avgBeta[,2] + avgBeta[,4]*1
quantile(logistic(noPred), probs=c(0.025,0.5,0.975))
quantile(logistic(pred), probs=c(0.025,0.5,0.975))
#########################

#########################
# No predators clutch of 1
(np1 <- round(mean(logistic((avgBeta[,1]))),3))

# No predators clutch of 2
(np2 <- (round(mean(logistic((avgBeta[,1]+noPred))),3)))

# with predators clutch of 1
(p1 <- round(mean(logistic(avgBeta[,1]+avgBeta[,3])),3))

# with predators clutch of 2
(p2 <- round(mean(logistic(avgBeta[,1]+avgBeta[,3] + pred)),4))

np2-np1 # no predators
p2-p1   # predators

#########################

#########################
# extract the p's, the expected log-odds of the probabilities 
#of each observation, and model average them.
propObs <- obs/N
logitP3 <- as.matrix(m3, "p")
logitP4 <- as.matrix(m4, "p")
avgLogitP <- logitP4*wts[1]+logitP3*wts[2]

pz <- logistic(avgLogitP)
mnP <- colMeans(pz)

hdiP <- apply(pz, 2, HDI, credMass=0.95)

plot((obs/N) ~ N, ylim=c(0,1), las=1, type="n", ann=FALSE)
mtext("P(mortality)", side=2, line=2.2)
mtext("Clutch size", side=1, line=2)
# no predators
x <- N[treat==0]
y <- propObs[treat==0]
mnNP <- mnP[treat==0]
hdiNP <- hdiP[,treat==0]

points(x, y)
polygon(c(x, rev(x)), c(hdiNP[1,], rev(hdiNP[2,])), col="#50505050")
lines(x, mnNP, col="black")

# predators
x <- N[treat==1]
y <- propObs[treat==1]
mnPr <- mnP[treat==1]
hdiPr <- hdiP[,treat==1]

points(x, y, col="blue")
polygon(c(x, rev(x)), c(hdiPr[1,], rev(hdiPr[2,])), col="#88CCEE50", border="blue")
lines(x, mnPr, col="blue")
#########################

