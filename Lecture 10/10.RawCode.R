setwd("~/Dropbox/BayesClass/2018 Class/Lecture 10")

library(rstan)
library(shinystan)
library(car)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("../utilityFunctions.R")

##############################

# Read in data
algae <- read.csv("algae.csv")
names(algae)

##############################

##############################
## Run first linear model

# set up data
datLM <- list(nObs=nrow(algae), BM=algae$biomass, terpenes=algae$terpenes,
              aMean=150, aSD=100, bSD=10, sigmaSD=10)

# run model
modLM <- stan(file="10.modLM.stan", data=datLM, iter=2000, chains=4, seed=3)

# extract posterior estimates
parLM <- as.matrix(modLM, pars=c("alpha", "beta", "sigma"))

################################

################################
## plot marginal density plots of parameters

par(mar=c(3,3,0.1,0.5))
par(mfrow=c(1,3))

plotInterval(parLM[,"beta"], HDI=TRUE, interval=0.95, xlims=c(1, 5),
             col="cornflowerblue", yOffset=0.01)
mtext(expression(paste(bold(beta))), side=1, line=2, cex=1.2)

plotInterval(parLM[,"alpha"], HDI=TRUE, interval=0.95, xlims=c(-100, 60),
             col="cornflowerblue", yOffset=0.01)
mtext(expression(paste(bold(alpha))), side=1, line=2, cex=1.2)

plotInterval(parLM[,"sigma"], HDI=TRUE, interval=0.95, xlims=c(8, 16),
             col="cornflowerblue", yOffset=0.01)
mtext(expression(paste(bold(sigma))), side=1, line=2, cex=1.2)

## print textual results
print(modLM, pars=c("alpha", "beta", "sigma"), digits.summary=2)

#################################

#################################
## Print out correlation between alpha and beta
round(cor(parLM[,1:2]), 3)

## plot joint posterior of alpha and beta
par(mfrow=c(1,1))
par(mar=c(3,3.2,0.1,0.5))
col <- "#50505010"

# Plot alpha and beta
plot(parLM[,1:2], pch=16, col=col, las=1, ylim=c(2,5),
     xlim=c(-75,45), bty="l")

# Add ellipses
dataEllipse(as.matrix(parLM[,1:2]),level=c(0.95), add=TRUE, labels=FALSE,
            plot.points=FALSE, center.pch=FALSE, col=c(col,"#006DCC"))

mtext(text = "beta", side=2, line=2.3, cex=1.2)
mtext(text = "alpha", side=1, line=2, cex=1.2)

##################################

##################################
# Plot a subset of posterior parameter estimates
par(mar=c(3,3.2,0.1,0.5))

# Make an empty plot
plot(algae, type="n",xlim=c(0, 60), ylim=c(-65,200), las=1, bty="l", xlab="", ylab = "")

# plot regression lines by looping through the first 200 posterior
# iterations of "alpha"
for(i in 1:200) {
  abline(a=parLM[i,"alpha"], b=parLM[i,"beta"], col="#50505040")
}
points(algae, pch=16, col="red")
mtext(text = "Biomass", side=2, line = 2.3, cex=1)
mtext(text = "Terpenes", side=1, line = 2, cex=1)

####################################

####################################
## Center the data

# center (but not scale terpenes)
terpC <- as.vector(scale(algae$terpenes, scale=FALSE))
datC <- list(nObs=dim(algae)[1], BM=algae$biomass,
             terpenes=terpC, aMean=150, aSD=100, bSD=10, sigmaSD=10)

modC <- stan(file="10.modLM.stan", data=datC, iter=2000, chains=4, seed=3)

# extract posterior estimates
parC <- as.matrix(modC, pars=c("alpha", "beta"))
mu <- as.matrix(modC, pars="mu")
newBM <- as.matrix(modC, pars="newBM")

######################################

######################################
## Print results
print(modC, pars=c("alpha", "beta", "sigma"), digits_summary = 2)


## Plot joint posterior of alpha and beta again

par(mar=c(3,3.2,0.1,0.5))

# Plot alpha and beta
plot(parC, pch=16, col=col, las=1, ylim=c(2,5),
     xlim=c(145, 155), bty="l")
dataEllipse(parC,level=c(0.025, 0.5 ,0.95), add=TRUE, labels=FALSE,
            plot.points=FALSE, center.pch=FALSE, col=c(col,"#006DCC"))
mtext(text = expression(bold(alpha)), side=1, line=2, cex=1.2)
mtext(text = expression(bold(beta)), side=2, line=2.5, cex=1.2, las=1)
text(146.5, 4.8, expression(paste(rho, " = -0.004")), 
     adj=c(0.5, 0.5),  cex=1.2)

#########################################

#########################################
## Build up regression plot on part at a time

par(mfrow=c(1,1))
par(mar=c(3,3.2,0.1,0.5))

# Make an empty plot
plot(terpC, algae$biomass, type="n", ylim=c(100,200), las=1, bty="l", xaxt="n",
     ylab="")

# back transform the x-axis
axis(1, at=seq(-8, 6, by=2), 
     labels=round((seq(-8, 6, by=2)+mean(algae$terpenes))))
mtext(text = "Biomass", side=2, line = 2.3, cex=1)
mtext(text = "Terpenes", side=1, line = 2, cex=1)

# plot the data points and mean regression line
points(terpC, algae$biomass, pch=16, col="red")
abline(a=mean(parC[,"alpha"]), b=mean(parC[,"beta"]))

######################################

######################################
## plot posterior density for terpene 7
par(mar=c(3,3.1,0.1,0.5))
plot(density(mu[,7]), las=1, main="")
mtext(text = "mu | terpenes=46", side=1, line=2)

#######################################

#######################################
## replot regression with 95% uncertainty around mu

par(mar=c(3,3.2,0.1,0.5))
muHDI <- apply(mu,2, HDI, credMass=0.95)

# Make an empty plot
plot(terpC, algae$biomass, type="n", ylim=c(100,200), las=1, 
     bty="l", xaxt="n", xaxs="i", ylab="")

# back transform the x-axis
axis(1, at=seq(-8, 6, by=2), 
     labels=round((seq(-8, 6, by=2)+mean(algae$terpenes))))
mtext(text = "Biomass", side=2, line = 2.3, cex=1)
mtext(text = "Terpenes", side=1, line = 2, cex=1)

# plot uncertainty interval in mu as a polygon
polygon(x=c(terpC, rev(terpC)), y=c(muHDI[1,], rev(muHDI[2,])), 
        col="#50505080", border="grey80")

# plot the data points and mean regression line
points(terpC, algae$biomass, pch=16, col="red")
abline(a=mean(parC[,"alpha"]), b=mean(parC[,"beta"]))

############################################

############################################
## calculate 95% uncertainty in "data" and plot
bmHDI <- apply(newBM,2, HDI, credMass=0.95)

par(mar=c(3,3.2,0.1,0.5))
# Make an empty plot
plot(terpC, algae$biomass, type="n", ylim=c(100,200), las=1, 
     bty="l", xaxt="n", xaxs="i", xlab=FALSE)

# back transform the x-axis
axis(1, at=seq(-8, 6, by=2), 
     labels=round((seq(-8, 6, by=2)+mean(algae$terpenes))))
mtext(text = "Biomass", side=2, line = 2.3, cex=1)
mtext(text = "Terpenes", side=1, line = 2, cex=1)

# plot uncertainty interval in mu as a polygon
polygon(x=c(terpC, rev(terpC)), y=c(bmHDI[1,], rev(bmHDI[2,])), 
        col="lightblue", border="blue")

# plot uncertainty interval in mu as a polygon
polygon(x=c(terpC, rev(terpC)), y=c(muHDI[1,], rev(muHDI[2,])), 
        col="#50505080", border="grey80")

# plot the data points and mean regression line
points(terpC, algae$biomass, pch=16, col="red")
abline(a=mean(parC[,"alpha"]), b=mean(parC[,"beta"]))

#######################################