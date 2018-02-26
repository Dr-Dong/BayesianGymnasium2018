setwd("~/Dropbox/BayesClass/2018 Class/Lecture 9")

library(rstan)
library(shinystan)
library(car)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("../utilityFunctions.R")

#################################
## Simulate random walk in stadium

pos <- replicate(20, runif(100,-1, 1)) # simulate positions 
cumPos <- t(apply(pos,1,cumsum)) # calculate cumulative position at each step
cumPos <- cbind(rep(0,100), cumPos) # add initial step


### Plot
par(mar=c(3,3,0.1,0.5))
plot(1:100, cumPos[,21], xlim=c(0,20), type="n", las=1, axes=FALSE, xaxs="i")

axis(1, at=seq(0,20, by=5))
axis(2, at=seq(-10,10, by=5), las=1)
mtext(text = "Step", side=1, line = 2)
mtext(text = "position", side=2, line = 2.1)

for(i in 1:nrow(cumPos)) {
  lines(0:20, cumPos[i,], col="#50505070")
}


### Plot density
par(mar=c(3,3,0.1,0.5))
plot(density(cumPos[,16],adj=1.5, from=-10, to=10), main="", las=1, xlab="",ylab="")
axis(1, at=seq(-10,10, length=5), las=1)
mtext(text = "position at step 16", side=1, line = 2)
#################


##########################
# Simulate multiplicative gene interactions
prod(1 + runif(10, 0,0.1))

par(mar=c(3,3,0.1,0.5))
growth <- replicate(1000, prod(1+runif(12,0,0.1)))
plot(density(growth), main="",las=1)


## Log scale
par(mar=c(3,3,0.1,0.5))
logBig <- replicate(1000, log(prod(1+runif(12,0,0.5))))
plot(density(logBig), main="",las=1)
##############################


###############################
## Seaweed simulation function
# You give it the total number of desired observations, the intercept, 
# the slope, and the standard deviation and it makes a dataframe for you. 
# Play with it, changing parameters as you wish to see how the model
# differs. Also, in the script, terpenes is the x variable. I didn't
# include arguments to change that, but it would be easy by changing
# the mean=50 & sd=3 to whatever you want.

seaweedSim <- function(nObs, alpha, beta, sigma) {
  terpenes <- round(rnorm(nObs, mean=50, sd=3), digits=2) 
  error <- rnorm(nObs, mean=0, sd=sigma)
  biomass <- alpha + beta * terpenes + error
  out <- data.frame(terpenes, biomass)
  out <- out[order(terpenes),]
  return(out)
}

set.seed(20)
algae <- seaweedSim(nObs=50, alpha=0, beta=3, sigma=12)

########################

#########################
## Can also read in data as .csv
algae <- read.csv("algae.csv")
head(algae)
summary(algae)
##########################


#################################
## Cauchy vs Normal plot
par(mar=c(3,3,0.1,0.5))
# normal(0,5) curve
curve(dnorm(x, 0, 5), from=0, to=20, las=1)
# half-Cauchy(0,5) curve
curve(dcauchy(x, 0, 5), add=TRUE, col="blue")
# normal(0,10) curve
curve(dnorm(x, 0, 10), add=TRUE, col="red")
# half-Cauchy(0,10) curve
curve(dcauchy(x, 0, 10), add=TRUE, col="cornflowerblue")

mtext(text = expression(bold(sigma)), side=1, line = 2)
text(13.5, 0.075, "Normal (0,5)", font=1,cex=1, col="black", adj=c(0, 0.5))
text(13.5, 0.0675, "Cauchy (0,5)", font=1,cex=1, col="blue", adj=c(0, 0.5))
text(13.5, 0.06, "Normal (0,10)", font=1,cex=1, col="red", adj=c(0, 0.5))
text(13.5, 0.0525, "Cauchy (0,10)", font=1,cex=1, col="cornflowerblue", 
     adj=c(0, 0.5))

#############################


################################
# Run model

dat <- list(nObs=dim(algae)[1], BM=algae$biomass, muMean=150, 
            muSD=30, sigmaSD=10)

intMod <- stan(file="09.modMean.stan", data=dat, iter=2000, chains=4, seed=3)

parMod <- as.data.frame(intMod, pars=c("mu", "sigma"))


## Print results
print(intMod, pars=c("mu", "sigma"), digits.summary=2)

##################################


##################################
# Plot marginal densities
par(mar=c(3,3,0.15,0.5))
par(mfrow=c(1,2))

plotInterval(parMod$mu, HDI=TRUE, credMass=0.95, xlims=c(140, 160),
             col="blue", yOffset=0.01)
mtext(expression(paste(bold(mu))), side=1, line=2, cex=1.2)

plotInterval(parMod$sigma, HDI=TRUE, credMass=0.95, xlims=c(10, 25), 
             col="blue", yOffset=0.01)
mtext(expression(paste(bold(sigma))), side=1, line=2, cex=1.2)

####################################


####################################
# Plot joint posterior
dev.off()
par(mar=c(3,3,0.1,0.5))
col <- "#50505020"
plot(parMod, pch=16, col=col, las=1, ylim=c(10,25), 
     xlim=c(140,160), bty="l")
dataEllipse(as.matrix(parMod),level=c(0.25,0.5,0.95), add=TRUE, labels=FALSE,
            plot.points=FALSE, center.pch=FALSE, col=c(col,"#006DCC"))
mtext(text = expression(paste(sigma)), side=2, line=2.2, cex=1.2, las=1)
mtext(text = expression(paste(mu)), side=1, line=2, cex=1.2)

###################################


###################################
## Estimate biomass for the Dictyota population

par(mar=c(3,3.2,0.1,0.5))
# plot empty plot
plot(0:1,0:1, type="n", xlim=c(100, 200), ylim=c(0,0.035), las=1, bty="l") 
mtext(text = "Estimated biomass", side=1, line = 2, cex=1)

# Overlay posterior biomass densities
for (n in 1:nrow(parMod)) {
  curve(dnorm(x, parMod[n,1], parMod[n,2]), add=TRUE, col="#50505010")
}

# Overlay median posterior probability density
medBM <- apply(parMod,2,median)
curve(dnorm(x,medBM[1], medBM[2]), add=TRUE, col="cornflowerblue", lwd=3)

##############################


################################
# Look at Correlation between parameters
cor(parMod)

################################