library(rstan)
library(shinystan)
library(car)
library(xtable)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("../utilityFunctions.R")

###############################
# Read in Data
milk <- read.csv("milkFull.csv")
unique(milk$clade)

# Set up data for stan
obs <- milk$kcal
x   <- as.integer(milk$clade)
nObs <- nrow(milk)
nVar <- max(x)
aMu <- 0.6         
aSD <- sigmaSD <- 10

dat <- list(nObs=nObs, nVar=nVar, obs=obs, x=x, aMu=aMu, aSD=aSD, sigmaSD=sigmaSD)

# run model
intMod <- stan(file="13.interceptMod.stan", data=dat, iter=2000, chains=4, seed=867.5309)

# Look at summary
round(summary(intMod, pars=c("alpha", "sigma"),  probs = c(0.025, 0.5, 0.975))$summary,2)

###################################

###################################
# Interaction dataset

rugged <- read.csv("RuggedGDP.csv")
rugged <- rugged[order(rugged$rugged),]
afr <- rugged[rugged$africa == 1, ] # African dataset
nafr <- rugged[rugged$africa == 0, ] # non-African dataset

# order both subsets
afr <- afr[order(afr$rugged),]
nafr <- nafr[order(nafr$rugged),]


# African linear regression
afrDat <- list(nObs=nrow(afr), obs=log(afr$GDP), xvar=afr$rugged, aSD=20, bSD=1, sigmaSD=10)

afrMod <- stan(file="13.uniMod.stan", data=afrDat, iter=2000, chains=4, seed=867.5309)
afrMu <- as.matrix(afrMod, "mu")

# NonAfrican linear regression
nafrDat <- list(nObs=nrow(nafr), obs=log(nafr$GDP), xvar=nafr$rugged, aSD=20, bSD=1, sigmaSD=10)

nafrMod <- stan(file="13.uniMod.stan", data=nafrDat, iter=2000, chains=4, seed=867.5309)
nMu <- as.matrix(nafrMod, "mu")
###################################

###################################
# Plot univariate models
par(mar=c(3,3.2,0.1,0.5))
par(mfrow=c(1,2))

### AFRICA
# Mean & HDI 
afrHDI <- apply(afrMu,2, HDI, credMass=0.95)
afrMean <- colMeans(afrMu)

# Make an empty plot
x <- afrDat$xvar
y <- afrDat$obs
plot(x, y, type="n", las=1, bty="l")
mtext(text = "Ruggedness", side=1, line = 2, cex=1)
mtext(text = "log(GDP)", side=2, line = 2.2, cex=1)

# plot uncertainty interval in mu as a polygon
polygon(x=c(x, rev(x)), y=c(afrHDI[1,], rev(afrHDI[2,])), col="#50505080", border="black")

# plot the data points and mean regression line
points(x, y, pch=1, col="blue")
lines(afrMean~x, col="black", lwd=2)
text(3, 9.7, "Africa", font=2)

### NONAFRICA
nHDI <- apply(nMu,2, HDI, credMass=0.95)
nMean <- colMeans(nMu)

# Make an empty plot
x <- nafrDat$xvar
y <- nafrDat$obs
plot(x, y, type="n", las=1, bty="l")
mtext(text = "Ruggedness", side=1, line = 2, cex=1)

# plot uncertainty interval in mu as a polygon
polygon(x=c(x, rev(x)), y=c(nHDI[1,], rev(nHDI[2,])), col="#50505080", border="black")

# plot the data points and mean regression line
points(x, y, pch=1, col="black")
lines(nMean~x, col="black", lwd=2)
text(3, 10.9, "not Africa", font=2)

###################################

###################################
## Run multivariate model

X <- cbind(rugged$rugged, rugged$africa)

fullDat <- list(nObs=nrow(rugged), nVar=ncol(X), obs=log(rugged$GDP), X=X, aMu=0, aSD=20, bMu=0, bSD=1, sigmaSD=10)

fullMod <- stan(file="13.multMod.stan", data=fullDat, iter=2000,chains=4, seed=867.5309)
mu <- as.matrix(fullMod, "mu")
muHDI <- apply(mu, 2, HDI, credMass=0.95)
muMn <- colMeans(mu)
  
##################################
# Plot results

par(mar=c(3,3.2,0.1,0.5))
par(mfrow=c(1,1))
### AFRICA
# Mean & HDI 
afrHDI <- muHDI[,rugged$africa==1]
afrMean <- muMn[rugged$africa==1]

### not AFRICA
# Mean & HDI
nHDI <- muHDI[,rugged$africa==0]
nMean <- muMn[rugged$africa==0]

# Make an empty plot
x <- nafrDat$xvar
y <- nafrDat$obs
plot(x, y, type="n", las=1, bty="l", ylim=c(6,11))
mtext(text = "Ruggedness", side=1, line = 2, cex=1)
mtext(text = "log(GDP)", side=2, line = 2.2, cex=1)

# plot uncertainty interval in mu as a polygon
polygon(x=c(x, rev(x)), y=c(nHDI[1,], rev(nHDI[2,])), col="#50505080", border="black")

# plot the data points and mean regression line
points(x, y, pch=1, col="black")
lines(nMean~x, col="black", lwd=2)

### AFRICA
x <- afrDat$xvar
y <- afrDat$obs

# plot uncertainty interval in mu as a polygon
polygon(x=c(x, rev(x)), y=c(afrHDI[1,], rev(afrHDI[2,])), col="#88CCEE80", border="#88CCEE")

# plot the data points and mean regression line
points(x, y, pch=1, col="blue")
lines(afrMean~x, col="blue", lwd=2)

###################################

###################################
# transformed parameters {
#   vector[nObs] mu;
#   vector[nObs] gamma;
#   
#   gamma = betaR + betaAR*A;
#   mu = alpha + to_vector(gamma * R') + betaA*A;
# }

###################################

###################################
# African linear regression
X <- model.matrix(~rugged*africa, data=rugged)[,2:4]

intDat <- list(nObs=nrow(rugged), obs=log(rugged$GDP),R=X[,"rugged"], A=X[,"africa"], aMu=0, aSD=20, bMu=0, bSD=1, sigmaSD=10)

intxnMod <- stan(file="13.intrxnMod.stan", data=intDat, iter=2000, chains=4, seed=867.5309)

mu <- as.matrix(intxnMod, "mu")
muHDI <- apply(mu, 2, HDI, credMass=0.95)
muMn <- colMeans(mu)

###################################

###################################
## Plot out results
par(mar=c(3,3.2,0.1,0.5))
par(mfrow=c(1,2))
### AFRICA
# Mean & HDI
afrHDI <- muHDI[,rugged$africa==1]
afrMean <- muMn[rugged$africa==1]

### not AFRICA
# Mean & HDI
nHDI <- muHDI[,rugged$africa==0]
nMean <- muMn[rugged$africa==0]

### AFRICA
x <- afr$rugged
y <- log(afr$GDP)
plot(x, y, type="n", las=1, bty="l")
mtext(text = "Ruggedness", side=1, line = 2, cex=1)
mtext(text = "log(GDP)", side=2, line = 2.2, cex=1)

# plot uncertainty interval in mu as a polygon
polygon(x=c(x, rev(x)), y=c(afrHDI[1,], rev(afrHDI[2,])), col="#50505080", border="black")

# plot the data points and mean regression line
points(x, y, pch=1, col="blue")
lines(afrMean~x, col="black", lwd=2)

### NOT AFRICA
# Make an empty plot
x <- nafr$rugged
y <- log(nafr$GDP)
plot(x, y, type="n", las=1, bty="l")
mtext(text = "Ruggedness", side=1, line = 2, cex=1)

# plot uncertainty interval in mu as a polygon
polygon(x=c(x, rev(x)), y=c(nHDI[1,], rev(nHDI[2,])), col="#50505080", border="black")

# plot the data points and mean regression line
points(x, y, pch=1, col="black")
lines(nMean~x, col="black", lwd=2)
###################################

###################################
round(summary(intxnMod, pars=c("alpha", "betaR", "betaA", "betaAR"), probs = c(0.025, 0.5, 0.975))$summary,2)

###################################

###################################
slopes <- as.matrix(intxnMod, pars=c("betaR", "betaAR"))
gammaA <- slopes[,"betaR"] + slopes[,"betaAR"]*1
gammaNA <- slopes[,"betaR"] + slopes[,"betaAR"]*0
mean(gammaA)
mean(gammaNA)

###################################

###################################
par(mar=c(3,3.2,0.1,0.5))
par(mfrow=c(1,2))
plot(density(gammaNA), main="", xlim=c(-0.5,0.6), col="black",
     lwd=2, las=1)
lines(density(gammaA), col="blue",lwd=2)
mtext(text = expression(gamma), side=1, line = 2, cex=1.2)
mtext(text = expression(p(gamma)), side=2, line = 1.9, cex=1.2)
text(0.3,4, "Africa")
text(0.05, 5, "Not Africa")
plot(density(gammaA-gammaNA), main="", xlim=c(-0.2,0.8),
     col="red",lwd=2, lty=1, las=1)
mtext(text = expression(paste(Delta," ", gamma)), side=1, line = 2, cex=1.2)

###################################

