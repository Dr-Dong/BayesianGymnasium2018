library(rstan)
library(shinystan)
library(car)
library(mvtnorm)
library(rethinking)
library(MASS)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/Dropbox/BayesClass/2018 Class/Lecture 16")

#############################
#Logistic link plot

par(mar=c(3,3.2,0.1,0.5))
par(mfrow=c(1,2))

curve(2.3*x, from=-1, to=1, ylim=c(-4,4), col="cornflowerblue", lwd=2, las=1)
mtext(text = "x", side=1, line = 2, cex=1)
mtext(text = "log-odds", side=2, line = 2, cex=1)
abline(h=seq(-4, 4, by=1), col="#50505050")

curve(logistic(2.3*x), from=-1, to=1, ylim=c(0,1), ann=FALSE, las=1, lwd=2)
abline(h=logistic(seq(-4,4,by=1)), col="#50505050")
mtext(text = "x", side=1, line = 2, cex=1)
mtext(text = "probability", side=2, line=2.4, cex=1)
#############################

#############################
# log-link plot

par(mar=c(3,3.2,0.1,0.5))
par(mfrow=c(1,2))
curve(2.3*x, from=-1, to=1, ylim=c(-3,3), col="cornflowerblue", 
      lwd=2, las=1)
mtext(text = "x", side=1, line = 2, cex=1)
mtext(text = "log(abundance)", side=2, line = 2, cex=1)
abline(h=seq(-3, 3, by=1), col="#50505050")

curve(exp(2.3*x), from=-1, to=1, ylim=c(0,10), ann=FALSE, 
      las=1, lwd=2)
abline(h=exp(seq(-3,3,by=1)), col="#50505050")
mtext(text = "x", side=1, line = 2, cex=1)
mtext(text = "Abundance count", side=2, line=2.2, cex=1)

#############################

