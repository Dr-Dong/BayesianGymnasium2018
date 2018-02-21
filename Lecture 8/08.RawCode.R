library(rstan)
library(shinystan)
library(car)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("../utilityFunctions.R")


##############################
# Simulation code

omega <- 0.3 # population mode (avg. infection rate across sites)
kappa <- 20  # Concentration (how variable sites are)
nObs <- 60   # total number of sites

a <- omega*(kappa-2) +1       # Beta shape parameter
b <- (1-omega)*(kappa-2) +1   # Beta shape parameter

set.seed(5) # for repeatability 
# Total n. of newts sampled at a site. Sites have 5, 10, 25, or 40 newts 
N <- as.integer(rep(c(5, 10, 25, 40), each=15)) 
trueInfect <- rbeta(nObs, a, b)       # true infection prevalance at a site
infect <- rbinom(nObs, N, trueInfect) # simulated number of infected newts
# Site indicator (not really necessary here but included for generality)
site <- 1:nObs 

newtDat <- data.frame(site=site, N=N, trueInfect=trueInfect, 
                      infect=infect, propInf=infect/N)

########################

########################
# Read in data

newtDat <- read.csv("newtDat.csv")
head(newtDat)

#########################

#########################
# Set up the data

nObs    <- nrow(newtDat)
nSites  <- max(newtDat$site)
site    <- newtDat$site
N       <- newtDat$N
obs     <- newtDat$infect
alpha   <- 2                  # weak beta
beta    <- 2                  # priors
sigma   <- 2.5                # regularizing scale prior

dat <- list(nObs=nObs, nSites=nSites, site=site, N=N, obs=obs,
            alpha=alpha, beta=beta, sigma=sigma)

#############################

#############################
# First model: Complete pooling

modCP <- stan(file="completePooling.stan", data=dat, iter=2000,  
              chains=4, seed=3, verbose = FALSE)

thetaCP <- as.matrix(modCP, pars="theta")

# Plot posterior distribution of theta

par(mar=c(3,3,0.1,0.5))
plot(density(thetaCP, adj=2), main="", xlab="", ylab="",las=1)
mtext(text = expression(paste(theta)), side=1, line = 2)
mtext(text = "P(theta|y)", side=2, line = 2.1)

################################

################################
# No pooling model

modNP <- stan(file="noPooling.stan", data=dat, iter=2000,  
              chains=4, seed=3, verbose = FALSE)

thetaNP <- as.matrix(modNP, pars="theta")
medianNP <- apply(thetaNP, 2, median)

# Plot out posterior distributions of no-pooling thetas

plot(modNP, pars="theta")

###################################

###################################
# Partial-pooling model

modPP <- stan(file="partialPooling.stan", data=dat, iter=2000,  
              chains=4, seed=3, verbose = FALSE)

omega   <- as.matrix(modPP, pars="omega")
kappa   <- as.matrix(modPP, pars="kappa")

thetaPP <- as.matrix(modPP, pars="theta")
medianPP <- apply(thetaPP, 2, median)


#######
# plot results

par(mar=c(3,3,0.1,0.5))

# Display no-pooling estimates of infection for each pond
plot(newtDat$propInf, col="blue", pch=16, xaxt="n", las=1, xlab="", ylab="") 
axis(1, at=c(1, 15, 30, 45, 60))
mtext(text = "site", side=1, line = 2, cex=1)
mtext("p(infection)", side=2, line=2.2, cex=1)

## Draw vertical dividers between site sample sizes
abline(v=c(15.5, 30.5, 45.5))
text(8, 0, "N=5")
text(22.5, 0, "N=10")
text(37.5, 0, "N=25")
text(53, 0, "N=40")

# Add in partial-pooling estimates
points(medianPP)

# Add in median of population-level mode
abline(h=median(omega), lty=3, lwd=2)

# Connect no-pooling estimates to partial-pooling estimates
for(i in 1:length(medianNP)) {
  arrows(x0=i, x1=i, y0=newtDat$propInf[i], y1=medianPP[i], length = 0.05)
}

##################################
# Plot infection proportions as a whole for east TN

# calculate a & b for each posterior iteration and turn them into a matrix
par(mar=c(3,3,0.1,0.5))
a <- omega * (kappa - 2) + 1           
b <- (1-omega) * (kappa - 2) + 1
propTN <- cbind(a,b)

# Plot empty plot
plot(0:1,0:1, type="n", ylim=c(0,3.2)) # plot empty plot
mtext(text = "Probability of infection", side=1, line = 2, cex=1)

# Overlay infection probability densities for each posterior draw
for (n in 1:nrow(propTN)) {
  curve(dbeta(x, propTN[n,1], propTN[n,2]), add=TRUE, col="#50505020")
}

# Overlay median posterior probability density
medianTN <- apply(propTN,2,median)  
curve(dbeta(x,medianTN[1], medianTN[2]), add=TRUE, col="cornflowerblue", lwd=3)

##################################

##################################

#Calculate error for no-pooling and partial pooling models

par(mar=c(3,4,0.1,0.5))
errorNP <- abs(medianNP - newtDat$trueInfect)
errorPP <- abs(medianPP - newtDat$trueInfect)

# calculate the average error for each sample size
avgErrNP <- aggregate(errorNP~N, FUN=mean)$errorNP
avgErrPP <- aggregate(errorPP~N, FUN=mean)$errorPP

# plot the no pooling error in blue
plot(errorNP, col="blue", pch=16, xaxt="n", las=1, xlab="", ylab="")
axis(1, at=c(1, 15, 30, 45, 60))
mtext(text = "site", side=1, line = 2, cex=1)
mtext("absolute error", side=2, line=2.5, cex=1)


# denote the different sample sizes 
abline(v=c(15.5, 30.5, 45.5))
text(8, 0.275, "N=5")
text(22.5, 0.275, "N=10")
text(37.5, 0.275, "N=25")
text(53, 0.275, "N=40")

# plot the partial pooling error estimates
points(errorPP)

# Plot the average error as horizontal bars

segments(x0=c(1, 16, 31, 46), x1=c(14, 29, 44, 59), y0=avgErrPP, 
         y1=avgErrPP, lty=2, lwd=2)
segments(x0=c(1,16,31,46), x1=c(14, 29, 44, 59), y0=avgErrNP, 
         y1=avgErrNP, lty=1, lwd=2, col="blue")

#################################