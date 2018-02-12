results='hide'}
library(rstan)
library(shinystan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("../utilityFunctions.R")


nObs <- 5
N <- c(2, 3, 3, 2, 4)
obs <- c(1, 3, 3, 2, 4)
omega <- 0.5
kappa <- 2

dat <- list(nObs=nObs, N=N, obs=obs, omega=omega, kappa=kappa)

mod1 <- stan(file="06.binomialMode.stan", #path to .stan file
             data=dat,
             iter=2000, # number of MCMC iterations
             chains=4,  # number of independent MCMC chains
             seed=3, 
             verbose = FALSE) # get rid of stupid messages and warnings 

p <- as.matrix(mod1, "p") # Extract posterior of parameter p


#######
pDens <- as.data.frame(density(p, from=0, to=1.2, n=1024, adjust=2)[1:2])

# Use which.max to identify the index for the highest density value
MAP <- pDens$x[which.max(pDens$y)]
plot(pDens, xlim = c(0.4,1), xlab="p", main="", las=1, type="l")
abline(v=c(mean(p), median(p), MAP), lty=1:3, 
       col=c("black", "blue", "red"), lwd=3)

#########
#What is the probability the proportrion of water is less than 0.75 or between 0.8-0.9?

sum(p < 0.75)/length(p) # less than 0.75
int <- pDens[pDens$x < 0.75,]
sum(p > 0.8 & p < 0.9)/sum(p)
plot(pDens, xlim = c(0.25,1), xaxs = "i", yaxs = "i", main="",
     las = 1, type = "l", bty="l")
polygon(x=c(int$x, rev(int$x)), y=c(rep(0,dim(int)[1]), rev(int$y)), col="blue")

#######
# quantiles
quantile(p, probs=c(0.025, 0.975))

###############

round(quantile(p, c(0.25, 0.75)), 2)
plotInterval(p, probs=c(0.25,0.75), HDI=FALSE, col="blue")

###############


round(quantile(p,0.75),2)

####################

HDI(p, credMass=0.5)

#####################

plotInterval(p, HDI=TRUE, interval=0.5, col="blue")

#######################

# Stan model with much fewer posterior iterations
mod200 <- stan(file="06.binomialMode.stan", #path to .stan file
               data=dat,
               iter=100, # number of MCMC iterations
               chains=4,  # number of independent MCMC chains
               seed=3, 
               verbose = FALSE) # get rid of stupid messages and warnings 

p200 <- as.matrix(mod200, "p")
HDI(p200,0.5) # 200 iterations
HDI(p,0.5) # 4000 iterations


#########################
# Plot it out
par(mfrow=c(1,2))
par(mar=c(3,3,0.1,0.5))

# using a grey50 with 50% transparency
plot(p200, pch=16, col="#50505050",las=1, ylim=c(0.5,1))
plot(p, pch=16, col="#50505050",las=1, ylim=c(0.5,1))


##########################
# Generated quantities model
set.seed(2)
nObs <- 20
N <- round(runif(nObs, 4,10),digits=0)
obs <- rbinom(nObs, size = N, prob=0.73)
omega <- 0.5
kappa <- 4
dat1 <-  list(nObs=nObs, N=N, obs = obs, omega=omega, kappa=kappa)

modGQ <- stan(file="06.binomialGQ.stan", #path to .stan file
              data=dat1,
              iter=2000, # number of MCMC iterations
              chains=4,  # number of independent MCMC chains
              seed=3, 
              verbose = FALSE) # get rid of stupid messages and warnings 

yNew <- as.matrix(modGQ, "yNew")

##############################

obsMat <- t(replicate(obs,n = nrow(yNew), simplify=TRUE))
mean(apply(obsMat > yNew, 2, mean))

##############################

yAvg <- apply(yNew,2, median)
plot(yAvg, obs, las=1, pch=16)
abline(a=0,b=1)



