setwd("~/Dropbox/BayesClass/2018 Class/Lecture 12")
library(rstan)
library(shinystan)
library(car)
library(xtable)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("../utilityFunctions.R")

#############################
## simulation of height and leg lengths

set.seed(2)
N <- 100                  # number of individuals
height <- rnorm(N, 10, 2) # simulated heights

legProp <- runif(N, 0.4, 0.5) # leg as proportion of height

leftLeg <- height*legProp + rnorm(N, 0, 0.02) # sim left leg
rightLeg <- height*legProp + rnorm(N, 0, 0.02)# sim right leg

###############################

###############################
# run model...will take a while!
dat <- list(nObs=N, nVar=2, obs=height, X=cbind(leftLeg, rightLeg),
            aMu=10, aSD=10, bMu=2, bSD=10, sigmaSD=10)

legsMod <- stan(file="12.multMod.stan", data=dat, iter=2000,
                chains=4, seed=867.5309, pars="mu", include
                =FALSE)
################################

################################
# Plot results and print summary
plotLegs <- plot(legsMod, pars=c("alpha", "beta", "sigma"), ci_level=0.5)
plotLegs + theme(text=element_text(family="ArialMT"))

round(summary(legsMod, pars=c("alpha", "beta", "sigma"),  
              probs = c(0.025, 0.5, 0.975))$summary,2)

##################################

##################################
# Plot joint posterior of right nad left leg

par(mar=c(3,3.2,0.1,0.5))
col <- "#50505010"

legs <- as.matrix(legsMod, pars="beta")

# Plot alpha and beta
plot(legs, pch=16, col=col, las=1, xlim=c(-3,12), 
     ylim=c(-10,5), bty="l", ann=FALSE)

# Plot confidence ellipses
dataEllipse(legs,level=c(0.95), add=TRUE, labels=FALSE,
            plot.points=FALSE, center.pch=FALSE, col=c(col,"#006DCC"))

mtext(text=expression(beta[R]), side=2, line=2.1, cex=1.2, las=1)
mtext(text=expression(beta[L]), side=1, line=2, cex=1.2)

#################################

#################################
# Compute posterior distribution of the sums of beta
# and plot
par(mar=c(3,3.2,0.1,0.5))
plotInterval(rowSums(legs), HDI=TRUE, interval=0.95, xlims=c(1.8, 2.4),
             col="cornflowerblue", yOffset=0.1)
mtext(expression(paste(beta[L] + beta[R])), side=1, line=2.1, cex=1.2)

###################################

###################################
# Fit regression with just one of the legs

dat1 <- list(nObs=N, nVar=2, obs=height, xvar=leftLeg, 
             aMu=10, aSD=10, bMu=2, bSD=10, sigmaSD=10)

legMod <- stan(file="12.uniMod.stan", data=dat1, iter=2000, chains=4, 
               seed=2, pars=c("alpha","beta", "sigma"))

round(summary(legMod, probs = c(0.025, 0.5, 0.975))$summary,2)

#####################################

#####################################
# function to simulate collinearity for 1 iter
collSim <- function (rho=0.9) {
  x <- rnorm(N, mean=rho*leftLeg, sd=(1-rho^2)*var(leftLeg))
  m <- lm(height ~ leftLeg + x)
  return(sqrt(diag(vcov(m))[2]))
}

# function to replicate collSim iter times for a given rho
repColSim <- function (rho=0.9, iter=100) {
  stdev <- replicate(iter, collSim(rho))
  return(mean(stdev))

# Run the functions over a sequence of rhos    
  rhoSeq <- seq(0, 0.99, by=0.01)  # sequence of correlations
  stddev <- sapply(rhoSeq, function(z) {repColSim(rho=z, iter=100)})  
}

#####################################

#####################################
# Plot SD of betas as a function of rho
par(mar=c(3,3.6,0.1,0.5))
plot(rhoSeq, stddev, type="l", bty="l", ann=FALSE, las=1, col="red", lwd=2)
mtext(text=expression(bar(SD[beta])), side=2, line=2.1, cex=1.2)
mtext(text = expression(rho), side=1, line=2, cex=1.2)

######################################

######################################
# Read in plant data

plants <- read.csv("plants.csv")
head(plants)

obs <- plants$h1
# create a design matrix with initial height centered but not scaled
X <- cbind(as.vector(scale(plants$h0, scale=FALSE)), 
           plants$treats, plants$fungus)

#######################################

#######################################
# Run Bayesian model for plant dataset A
dat <- list(nObs=nrow(plants), nVar=dim(X)[2], obs=obs, X=X, 
            aMu=10, aSD=10, bMu=0, bSD=10, sigmaSD=10)

plantMod <- stan(file="12.multMod.stan", data=dat, iter=2000,
                 chains=4, seed=2, pars=c("alpha","beta", 
                 "sigma"))

round(summary(plantMod, pars=c("alpha", "beta", "sigma"), 
              probs = c(0.025, 0.5, 0.975))$summary,2)

########################################

########################################
# Run model excluding fungus

X1 <- cbind(as.vector(scale(plants$h0, scale=FALSE)), 
            plants$treats)

datNF <- list(nObs=nrow(plants), nVar=dim(X1)[2], obs=obs, X=X1, 
              aMu=10, aSD=10, bMu=0, bSD=10, sigmaSD=10)

modNF <- stan(file="12.multMod.stan", data=datNF, iter=2000, 
              chains=4, seed=2, pars=c("alpha","beta", "sigma"))

round(summary(modNF, pars=c("alpha", "beta", "sigma"), 
              probs = c(0.025, 0.5, 0.975))$summary,2)

####################################

####################################
# Read in milk data
milk <- read.csv("milkFull.csv")
unique(milk$clade)
###################################

###################################
# Make dummy variable for Streppsirrhines
(strep <- ifelse(milk$clade == "Strepsirrhine", 1, 0))

# Make dummy variable for new world monkeys
nwm <- ifelse(milk$clade == "New World Monkey", 1, 0)

# Dummy variable for old world monkeys
owm <- ifelse(milk$clade == "Old World Monkey", 1, 0)

clade <- cbind(nwm, owm, strep)
###################################

###################################
#Can also use model.matrix

X <- model.matrix(~clade, data=milk)
X <- as.matrix(X) # needed to get rid of the attribute list

####################################

###################################
# Fit the model
nObs <- nrow(milk)
nVar <- ncol(clade)
obs <- milk$kcal
X   <- clade
aMu <- bMu <- 0.6         
aSD <- sigmaSD <- 10
bSD <- 1

dat <- list(nObs=nObs, nVar=nVar, obs=obs, X=X, aMu=aMu, aSD=aSD, 
            bMu=bMu, bSD=bSD, sigmaSD=sigmaSD)

milkMod <- stan(file="12.multMod.stan", data=dat, iter=2000,
                chains=4, seed=867.5309, pars="mu", include=FALSE)

round(summary(milkMod, pars=c("alpha", "beta", "sigma"),  
              probs = c(0.025, 0.5, 0.975))$summary,2)

###################################

###################################
# Extract Posterior estimates
post <- as.matrix(milkMod, pars=c("alpha","beta"))

# bind together the intercept, and then for beta 2-4 use 
# Sweep function to add them to the intercept
mu <- cbind(post[,1], sweep(post[,2:4], MARGIN=1, STATS = post[,1], FUN='+'))

# Calculate column means for mu
means <- colMeans(mu)

# Calculate  column SD
SD <- apply(mu, 2, sd)

# Calculate 95% HDI
muHDI <- apply(mu, 2, HDI, credMass=0.95)

# put everything together
sumTab <- data.frame(Mean=means, SD=SD, lower.95=muHDI[1,],
                     upper.95=muHDI[2,], row.names=c("Ape", "NWM", "OWM", "Strep"))

#########################################

#########################################
# Plot results as marginal plot and  contrast plot
par(mar=c(3,3.2,0.1,0.5))
par(mfrow=c(1,2))
plot(density(mu[,2], adj=2), las=1, col="blue", lwd=2, main="", 
     xlim=c(0.5, 1))
lines(density(mu[,3], adj=2), col="black", lwd=2)
mtext(text = "Marginal milk energy", side=1, line = 2, cex=1)
text(0.6,8, "New World\nMonkeys", col="blue")
text(0.9,7, "Old World\nMonkeys", col="black")

# Contrast plot
dif <- mu[,3]-mu[,2]
plotInterval(dif, HDI = TRUE, interval = 0.95, xlims=c(-0.2,0.3), col="#50505080")
mtext(text = "Joint milk energy difference", side=1, line = 2, cex=1)

# Calculate proportion of times mu 3 < mu 2
round(sum(mu[,3]<mu[,2])/nrow(mu),2)
###################################

###################################
# Run multiple intercepts model
obs <- milk$kcal
x   <- as.integer(milk$clade)
nObs <- nrow(milk)
nVar <- max(x)
aMu <- 0.6         
aSD <- sigmaSD <- 10

dat <- list(nObs=nObs, nVar=nVar, obs=obs, x=x, aMu=aMu, aSD=aSD, 
            sigmaSD=sigmaSD)

intMod <- stan(file="12.interceptMod.stan", data=dat, iter=2000,
               chains=4, seed=867.5309)

round(summary(intMod, pars=c("alpha", "sigma"),  
              probs = c(0.025, 0.5, 0.975))$summary,2)
#################################