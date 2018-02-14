library(rstan)
library(shinystan)
library(car)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("../utilityFunctions.R")

################### set up data for first hierarchical model

# Low certainty in omega, high kappa
alpha_omega <- 2
beta_omega <- 2
kappa <- 100

N <- 6
obs <- 5
nObs <- length(N)
datEx1 <- list(alpha_omega=alpha_omega, beta_omega=beta_omega, 
               kappa=kappa, nObs=nObs, N=N, obs=obs)

# No likelihood
modExNL1 <- stan(file="07.simpleHierNL1.stan", data=datEx1, 
                 iter=2000, chains=4, seed=3, verbose = FALSE)
parNL <- as.matrix(modExNL1, pars=c("theta","omega"))

par(mar=c(3,3,0.1,0.5))
plot(density(parNL[,"omega"], adj=2), main="", xlab="", ylab="",las=1)
mtext(text = expression(paste(omega)), side=1, line = 2)

par(mar=c(3,3,0.1,0.5))
col <- "#50505008" # specify nice grey color with transparancy
plot(parNL, pch=16, col=col, las=1)
dataEllipse(parNL,level=c(0.25,0.5,0.95), add=TRUE, labels=FALSE, 
            plot.points=FALSE, center.pch=FALSE, col=c(col,"#006DCC"))
mtext(text = expression(paste(omega)), side=2, line = 2.2)
mtext(text = expression(paste(theta)), side=1, line = 2)


# With likelihoood
modEx1 <- stan(file="07.simpleHier1.stan", data=datEx1, iter=2000, 
               chains=4, seed=3, verbose = FALSE)
parL <- as.matrix(modEx1, pars=c("theta","omega"))


par(mar=c(3,3,0.1,0.5))
plot(density(parL[,"omega"], adj=2), main="", xlab="", ylab="",las=1)
lines(density(parNL[,"omega"], adj=2),col="blue",lty=2)
mtext(text = expression(paste(omega)), side=1, line = 2)


par(mar=c(3,3,0.1,0.5))
par(mfrow=c(1,2))
plot(parL, pch=16, col=col, las=1)
dataEllipse(parL,level=c(0.25,0.5,0.95), add=TRUE, labels=FALSE, 
            plot.points=FALSE, center.pch=FALSE, col=c(col,"#006DCC"))
mtext(text = expression(paste(omega)), side=2, line = 2.2)
mtext(text = expression(paste(theta)), side=1, line = 2)

plot(density(parL[,"theta"], adj=2), main="", xlab="", ylab="",las=1)
mtext(text = expression(paste(theta)), side=1, line = 2)


##################### High certainty in omega, low kappa

a_omega <- 20
b_omega <- 20
kappa <- 5

N <- 6
obs <- 5
datEx2 <- list(a_omega=a_omega, b_omega=b_omega, kappa=kappa,
               N=N, obs=obs)

#### no likelihood
modExNL2 <- stan(file="07.simpleHierNL1.stan", data=datEx2,
                 iter=2000, chains=4, seed=3, verbose = FALSE)
parNL2 <- as.matrix(modExNL2, pars=c("theta","omega"))


par(mfrow=c(1,2))
par(mar=c(3,3,0.1,0.5))
plot(density(parNL2[,"omega"], adj=2), main="", xlab="", ylab="",
     las=1,ylim=c(0,4.5))
mtext(text = expression(paste(omega)), side=1, line = 2)

par(mar=c(3,3,0.1,0.5))
plot(density(parNL2[,"theta"], adj=2), main="", xlab="", ylab="",
     las=1, ylim=c(0,4.5))
mtext(text = expression(paste(theta)), side=1, line = 2)


par(mar=c(3,3,0.1,0.5))
col <- "#50505010"
plot(parNL2, pch=16, col=col, las=1, ylim=c(0,1))
dataEllipse(parNL2,level=c(0.25,0.5,0.95), add=TRUE, labels=FALSE, 
            plot.points=FALSE, center.pch=FALSE, col=c(col,"#006DCC"))
mtext(text = expression(paste(omega)), side=2, line = 2.2)
mtext(text = expression(paste(theta)), side=1, line = 2)


######### with likelihood
modEx2 <- stan(file="07.simpleHier1.stan", data=datEx2,
               iter=2000, chains=4, seed=3, verbose = FALSE)
parL2 <- as.matrix(modEx2, pars=c("theta","omega"))


par(mar=c(3,3,0.1,0.5))
plot(density(parL2[,"omega"], adj=2), main="", xlab="", ylab="",las=1)
lines(density(parNL2[,"omega"], adj=2),col="blue",lty=2)
mtext(text = expression(paste(omega)), side=1, line = 2)


par(mar=c(3,3,0.1,0.5))
par(mfrow=c(1,2))

plot(parL2, pch=16, col=col, las=1, xlim=c(0,1), ylim=c(0,1))
dataEllipse(parL2,level=c(0.25,0.5,0.95), add=TRUE, labels=FALSE, 
            plot.points=FALSE, center.pch=FALSE, col=c(col,"#006DCC"))
mtext(text = expression(paste(omega)), side=2, line = 2.2)
mtext(text = expression(paste(theta)), side=1, line = 2)

plot(density(parL2[,"theta"], adj=2), main="", xlab="", ylab="",las=1)
mtext(text = expression(paste(theta)), side=1, line = 2)

dev.off()
par(mar=c(3,3,0.1,0.5))
plot(0:1,0:1, xaxt="n", yaxt="n", type="n", bty="l", xlab="", ylab="")
polygon(x=c(0,1,1,0), y=c(c(0,0.4,0,0)))
polygon(x=c(0,1,1,0), y=c(c(1,1,1,0.6)))
mtext(text = "parameters", side=2, line = 2.2)
mtext(text = "data", side=1, line = 2)
text(0.3,0.9, "wide data", font=2,cex=1.5)
text(0.6,0.1, "long data", font=2, cex=1.5)

