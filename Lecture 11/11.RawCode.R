library(rstan)
library(shinystan)
library(car)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("../utilityFunctions.R")

#######################
# Read in reef data

reefs <- read.csv("urchinDat.csv")
head(reefs)

########################

########################
# run chiton model
nObs <- nrow(reefs)
urchins <- reefs$urchins
chitons <- reefs$chitons
seaStars <- reefs$seaStars

### Chiton model
o <- order(chitons)
chitonDat <- list(nObs=nObs, obs=urchins[o], xvar=as.vector(scale(chitons[o])), 
                  aSD=10, bSD=1, sigmaSD=10)

chitMod <- stan(file="uniMod.stan", data=chitonDat, iter=2000, 
                chains=4, seed=3)

# extract posterior estimates of alpha, beta, and mu
chitPar <- as.matrix(chitMod, pars=c("alpha", "beta", "mu"))
###########################

###########################
# Run sea star model

# Sea star model
o <- order(seaStars)
starDat <- list(nObs=nObs, obs=urchins[o], xvar=as.vector(scale(seaStars[o])), 
                aSD=10, bSD=1, sigmaSD=10)

starMod <- stan(file="uniMod.stan", data=starDat, iter=2000, 
                chains=4, seed=3)

# extract posterior estimates of alpha, beta, and mu
starPar <- as.matrix(starMod, pars=c("alpha", "beta", "mu"))

##########################

##########################
# Look graphically
# at results
par(mar=c(3,3.2,0.1,0.5))
par(mfrow=c(1,2))
# Mean & HDI for chitons
chitHDI <- apply(chitPar,2, HDI, credMass=0.95)
chitMean <- colMeans(chitPar)

# Make an empty plot
x <- chitonDat$xvar
y <- chitonDat$obs
plot(x, y, type="n", las=1, bty="l")

mtext(text = "Urchin density", side=2, line = 2.2, cex=1)
mtext(text = "Chiton density", side=1, line = 2, cex=1)

# plot uncertainty interval in mu as a polygon
polygon(x=c(x, rev(x)), y=c(chitHDI[1, -c(1:2)], 
                            rev(chitHDI[2, -c(1:2)])), col="#50505080", border="grey80")

# plot the data points and mean regression line
points(x, y, pch=16, col="red")
abline(a=chitMean[1], b=chitMean[2], col="red", lwd=2)

### Plot seastar resutls
starHDI <- apply(starPar,2, HDI, credMass=0.95)
starMean <- colMeans(starPar)
# Make an empty plot
x <- starDat$xvar
y <- starDat$obs
plot(x, y, type="n", las=1, bty="l")

mtext(text = "Urchin density", side=2, line = 2.2, cex=1)
mtext(text = "Sea star density", side=1, line = 2, cex=1)

# plot uncertainty interval in mu as a polygon
polygon(x=c(x, rev(x)), y=c(starHDI[1, -c(1:2)], 
                            rev(starHDI[2, -c(1:2)])), col="#50505080", border="grey80")

# plot the data points and mean regression line
points(x, y, pch=16, col="blue")
abline(a=starMean[1], b=starMean[2], col="blue", lwd=2)

# print summaries
round(summary(chitMod, pars=c("alpha", "beta"))$summary,2)

round(summary(starMod, pars=c("alpha", "beta"))$summary,2)
#########################

#########################
# Run multiple regression model
dat <- list(nObs=nObs, nVar=2, obs=urchins, x1=as.vector(scale(chitons)), 
            x2 = as.vector(scale(seaStars)), aSD=10, bSD=1, sigmaSD=10)

multMod <- stan(file="multiMod.stan", data=dat, iter=2000,
                chains=4, seed=867.5309, pars="mu", 
                include = FALSE)

# print results
round(summary(multMod, pars=c("alpha", "beta"))$summary,2)

####################################

###################################
# plot results
plotMult <- plot(multMod, pars=c("alpha", "beta", "sigma"), ci_level=0.5)
plotMult + theme(text=element_text(family="ArialMT"))
####################################

####################################
# Extract results for both chitons and sea stars
oCH <- order(chitons)
muCH <- as.matrix(multMod, pars="muCH")
chitHDI <- apply(muCH,2, HDI, credMass=0.95)
chitMean <- colMeans(muCH)[oCH]

oSS <- order(seaStars)
muSS <- as.matrix(multMod, pars="muSS")
starHDI <- apply(muSS,2, HDI, credMass=0.95)
starMean <- colMeans(muSS)[oSS]
#####################################

#####################################
# Make a counterfactual plot
par(mar=c(3,3.2,0.1,0.5))
par(mfrow=c(1,2))

# Make an empty plot
x <- chitonDat$xvar
y <- chitonDat$obs
plot(x, y, type="n", las=1, bty="l")

mtext(text = "Urchin density", side=2, line = 2.2, cex=1)
mtext(text = expression(paste("Chiton density | sea stars = ", bar(SS))), 
      side=1, line = 2, cex=1)

# plot uncertainty interval in mu as a polygon
polygon(x=c(x, rev(x)), y=c(chitHDI[1,oCH], 
                            rev(chitHDI[2,oCH])), col="#50505080", border="grey80")

# plot the data points and mean regression line
lines(x, chitMean, col="blue", lwd=2)


### Plot seastar resutls
x <- starDat$xvar
y <- starDat$obs
plot(x, y, type="n", las=1, bty="l")

mtext(text = "Urchin density", side=2, line = 2.2, cex=1)
mtext(text = expression(paste("Sea star density | chitons = ", bar(CH))), 
      side=1, line = 2, cex=1)

# plot uncertainty interval in mu as a polygon
polygon(x=c(x, rev(x)), y=c(starHDI[1,oSS], 
        rev(starHDI[2,oSS])), col="#50505080", border="grey80")

# plot the data points and mean regression line
lines(x, starMean, col="red", lwd=2)

#####################################

#####################################
# Read in milk data
milk <- read.csv("milk.csv")
head(milk)

mass <- log(milk$mass)
ncp <- milk$neocortex.perc
kcal <- milk$kcal.per.g
nObs <- nrow(milk)
#####################################

#####################################
# Run first model
# neoCortex
o <- order(ncp)
ncpDat <- list(nObs=nObs, obs=kcal[o], xvar=as.vector(scale(ncp[o])), 
               aSD=10, bSD=1, sigmaSD=10)

ncpMod <- stan(file="uniMod.stan", data=ncpDat, iter=2000, 
               chains=4, seed=867.5309)

# extract posterior estimates of alpha, beta, and mu
ncpPar <- as.matrix(ncpMod, pars=c("alpha", "beta", "mu"))
########################################

#######################################
# Run model for mass
# mass
o <- order(mass)
massDat <- list(nObs=nObs, obs=kcal[o], xvar=as.vector(scale(mass[o])),
                aSD=10, bSD=1, sigmaSD=10)

massMod <- stan(file="uniMod.stan", data=massDat, iter=2000,
                chains=4, seed=867.5309)

# extract posterior estimates of alpha, beta, and mu
massPar <- as.matrix(massMod, pars=c("alpha", "beta", "mu"))
#########################################

#########################################
# Plot

par(mar=c(3,3.2,0.1,0.5))
par(mfrow=c(1,2))
# Mean & HDI for NCP
ncpHDI <- apply(ncpPar,2, HDI, credMass=0.95)
ncpMean <- colMeans(ncpPar)

# Make an empty plot
x <- ncpDat$xvar
y <- ncpDat$obs
plot(x, y, type="n", las=1, bty="l", xaxt="n")
at <- seq(-2, 1.5, by=0.5)
axis(1, at=at, labels=round(at*sd(ncp) + mean(ncp)))
mtext(text = "kCal per g", side=2, line = 2.2, cex=1)
mtext(text = "% Neocortex", side=1, line = 2, cex=1)

# plot uncertainty interval in mu as a polygon

polygon(x=c(x, rev(x)), y=c(ncpHDI[1, -c(1:2)], 
      rev(ncpHDI[2, -c(1:2)])), col="#50505050", border="grey80")

# plot the data points and mean regression line
abline(a=ncpMean[1], b=ncpMean[2], col="red", lwd=2)


### Plot mass resutls
massHDI <- apply(massPar,2, HDI, credMass=0.95)
massMean <- colMeans(massPar)
# Make an empty plot
x <- massDat$xvar
y <- massDat$obs
plot(x, y, type="n", las=1, bty="l", xaxt="n")
at <- seq(-2, 1.5, by=1)
axis(1, at=at, labels=round(at*sd(mass) + mean(mass),1))
#mtext(text = "kCal per g", side=2, line = 2.2, cex=1)
mtext(text = "log(mass)", side=1, line = 2, cex=1)

# plot uncertainty interval in mu as a polygon
polygon(x=c(x, rev(x)), y=c(massHDI[1, -c(1:2)], 
      rev(massHDI[2, -c(1:2)])), col="#50505050", border="grey80")

# plot the data points and mean regression line
abline(a=massMean[1], b=massMean[2], col="blue", lwd=2)

#####################################

#####################################
# Print results
round(summary(ncpMod, pars=c("alpha", "beta"))$summary,2)

round(summary(massMod, pars=c("alpha", "beta"))$summary,2)
######################################

######################################
# Run both together
dat <- list(nObs=nObs, nVar=2, obs=kcal, x1=as.vector(scale(ncp)), 
            x2 = as.vector(scale(mass)), aSD=10, bSD=1, sigmaSD=10)

milkMod <- stan(file="multiMod.stan", data=dat, iter=2000,
                chains=4, seed=867.5309, pars="mu", include=FALSE)

round(summary(milkMod, pars=c("alpha", "beta"))$summary,2)
#######################################

####################################
# Extract parameters
oNCP <- order(ncp)
muNCP <- as.matrix(milkMod, pars="muCH")
ncpHDI <- apply(muNCP,2, HDI, credMass=0.95)
ncpMean <- colMeans(muNCP)[oNCP]

oM <- order(mass)
muM <- as.matrix(milkMod, pars="muSS")
massHDI <- apply(muM,2, HDI, credMass=0.95)
massMean <- colMeans(muM)[oM]

#########################################

#########################################
# Make a counterfactual plot

par(mar=c(3,3.2,0.1,0.5))
par(mfrow=c(1,2))

# Make an empty plot
x <- ncpDat$xvar
y <- ncpDat$obs
plot(x, y, type="n", las=1, bty="l", xaxt="n")
at <- seq(-2, 1.5, by=0.5)
axis(1, at=at, labels=round(at*sd(ncp) + mean(ncp)))
mtext(text = "kCal per g", side=2, line = 2.2, cex=1)
mtext(text = expression(paste("% neocortex | mass = ", bar(M))), 
      side=1, line = 2, cex=1)

# plot uncertainty interval in mu as a polygon
polygon(x=c(x, rev(x)), y=c(ncpHDI[1,oNCP], 
                            rev(ncpHDI[2,oNCP])), col="#50505080", border="grey80")

# plot the data points and mean regression line
lines(x, ncpMean, col="blue", lwd=2)


### Plot seastar resutls
x <- massDat$xvar
y <- massDat$obs
plot(x, y, type="n", las=1, bty="l", xaxt="n")
at <- seq(-2, 1.5, by=1)
axis(1, at=at, labels=round(at*sd(mass) + mean(mass),1))
mtext(text = expression(paste("log(mass) | % neocortex = ", bar(ncp))), 
      side=1, line = 2, cex=1)

# plot uncertainty interval in mu as a polygon
polygon(x=c(x, rev(x)), y=c(massHDI[1,oM], 
                            rev(massHDI[2,oM])), col="#50505080", border="grey80")

# plot the data points and mean regression line
lines(x, massMean, col="red", lwd=2)