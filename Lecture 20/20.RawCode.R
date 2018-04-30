library(rstan)
library(shinystan)
library(car)
library(mvtnorm)
library(rethinking)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("../utilityFunctions.R")


#################################
# Set up bee parameters
a <- 3.5        # average morning nectaring time
b <- -1       # average diff in afternoon time
sigma.a <- 1    # std dev. in intercept
sigma.b <- 0.5  # std dev. in slope
rho <- -0.7     # Correlation between intercept and slope

#################################

#################################
# build vector of two means
Mu <- c(a, b)

# build VCV matrix directly

cov.ab <- sigma.a * sigma.b * rho
Sigma <- matrix(c(sigma.a^2, cov.ab, cov.ab, sigma.b^2), 
                ncol=2)

# Or decompose it 
sigmas <- c(sigma.a, sigma.b)
Omega <- matrix(c(1, rho, rho, 1), nrow=2)
Sigma <- diag(sigmas) %*% Omega %*% diag(sigmas)

#################################

#################################
# Simulate data
n.plants <- 20 # no. of plants
library(MASS)
set.seed(5)
plantEffects <- mvrnorm(n.plants, Mu, Sigma)
aPlant <- plantEffects[,1]
bPlant <- plantEffects[,2]

#################################

#################################
# Visualize data in contour plot
par(mar=c(3,3.3,0.1,0.5))

plot(aPlant, bPlant, col="cornflowerblue", las=1, ann=FALSE, 
     pch=16)
mtext(expression(paste("Intercepts ",alpha)), side=1, line=2, 
      cex=1.1)
mtext(expression(paste("Slopes ",beta)), side=2, line=2.2, cex=1.1)

library(ellipse)
c.int <- c(0.05, 0.25, 0.5, 0.75, 0.95)
for(i in c.int) {
  lines(ellipse(Sigma, centre=as.vector(Mu), level=i), 
        col="#50505080")
}

#################################

#################################
# Simulate bees nectaring 
n.visits <- 10
afternoon <- rep(0:1, n.visits*n.plants/2)
plantID <- rep(1:n.plants, each=n.visits)

mu <- aPlant[plantID] + bPlant[plantID]*afternoon 
tau <- 0.5  # std dev. w/in plant

set.seed(2)
nectar <- rnorm(n.visits * n.plants, mu, tau)
dat <- data.frame(plant=plantID, afternoon=afternoon, 
                  nectar=nectar)
head(dat)

#################################

#################################
# LKJ prior plot
par(mar=c(3,3.1,0.1,0.5))

plot(density(rlkjcorr(1e4, K=2, eta=4)[,1,2]), lty=3, lwd=2, adj=1,
     ann=FALSE, las=1, col="blue")
mtext("Correlation", side=1, line=2)

lines(density(rlkjcorr(1e4, K=2, eta=2)[,1,2]), lty=1, lwd=2)
lines(density(rlkjcorr(1e4, K=2, eta=1)[,1,2]), lty=2, lwd=2, 
      col="red")

text(-0.13, 0.39, expression(paste(eta, " = 1")), adj=c(0,0), 
     cex=1.1, col="red")
text(-0.13, 0.6, expression(paste(eta, " = 2")), adj=c(0,0), cex=1.1)

text(-0.13, 0.9, expression(paste(eta, " = 4")), adj=c(0,0), 
     cex=1.1, col="blue")
#################################

#################################
# Set up data and run model
obs <- nectar
nObs <- nrow(dat)
d <- list(nObs=nObs, obs=obs, nPlants=max(dat$plant), plant=dat$plant, afternoon=dat$afternoon)

m1 <- stan(file = "varSlopesMod.stan", data = d, iter = 2000, 
           chains = 4, seed = 4)

post <- extract(m1, c("aPlant", "bPlant", "gamma", "Sigma"))

#################################

#################################
# Look at posterior distribution of correlations between
# slope and intercept
par(mar=c(3,3.1,0.1,0.5))
Omega <- as.matrix(m1, pars="Omega")
Rho <- Omega[,2]

plot(density(Rho), lty=1, lwd=2, ann=FALSE, las=1, col="blue", xlim=c(-1,1))
mtext("Correlation", side=1, line=2)

dens(rlkjcorr(1e4, K=2, eta=2)[,1,2], add=TRUE, lwd=2,lty=2, adj=2)
text(0,1, "prior")
text(-0.25,1.85, "posterior", col="blue")

#################################

#################################
# Plot results
par(mar=c(3,3.1,0.1,0.5))
# compute unpooled estimates from data
a1 <- sapply(1:n.plants, function(i) mean(nectar[plantID==i & 
                                                   afternoon == 0]))
b1 <- sapply(1:n.plants, function(i) mean(nectar[plantID==i & 
                                                   afternoon == 1]))-a1

# Extract posterior means of partially-pooled estimates
aPP <- colMeans(post$aPlant)
bPP <- colMeans(post$bPlant)

# plot both and connect with lines
plot(a1, b1, ann=FALSE, pch=16, col="cornflowerblue", 
     xlim=c(min(a1)-0.1, max(a1)+0.1), ylim=c(min(b1)-0.1, 
                                              max(b1)+0.1))
points(aPP,bPP, pch=1)

mtext(expression(paste("Intercepts ",alpha)), side=1, line=2, 
      cex=1.1)
mtext(expression(paste("Slopes ",beta)), side=2, line=2, cex=1.1)

for(n in 1:n.plants) lines(c(a1[n], aPP[n]), c(b1[n], bPP[n]))

### now superimpose contours of population
Mu.est <- colMeans(post$gamma)
Sigma.est <- colMeans(post$Sigma)

for(l in c.int) {
  lines(ellipse(Sigma.est, centre=Mu.est, level=l), col="#50505070")
}

#################################

#################################
# Run non-centered model

afternoon <- model.matrix(~dat$afternoon)
nPar <- ncol(afternoon)
hVar <- as.matrix(rep(1, max(dat$plant)))
d <- list(nObs=nObs, nPar=nPar, nPlants=max(dat$plant), obs=obs, 
          plant=dat$plant, hVar=hVar, afternoon=afternoon)

m2 <- stan(file = "noncentered.stan", data = d, iter = 2000, 
           chains = 4, seed = 4, pars="z", include=FALSE)

#################################

#################################
#
#################################

#################################
#
#################################

#################################
#
#################################

#################################
#
#################################

#################################
#
#################################

#################################
#
#################################

#################################
#
#################################

#################################
#
#################################

