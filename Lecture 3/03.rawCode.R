########### Code for Lecture 3. Note that some of it may be
## slightly different in class. 

N <- 10 # number of globe tosses
obs <- c("W", "W", "L", "W", "L", "L", "W", "W", "W", "W") # observed data
y <- ifelse(obs == "W", 1, 0) # convert to 1's & 0's

#################

# calculate bernoulli likelihood

peez <- seq(0,1, length=10000) # sequence of candidate p's

# vector of likelihood for the fixed data for each candidate theta
liks <- numeric(length = length(peez))

for(i in 1:length(peez)) {
  temp <- dbinom(y, 1, peez[i])
  liks[i] <- prod(temp)
}

# plot it
par(mar=c(3.5,5,0.1,0.5)) # set margins
plot(peez, liks, type="l", lwd=3, col="blue", las=1, ylim=c(0,0.0025), 
  frame.plot=FALSE, xlab="", ylab="")
mtext(text = "Likelihood", side=2, line = 4)
mtext(text = "p's", side=1, line = 2, font=3)

(maxLik <- peez[liks == max(liks)])


##########
#Log likelihoods
logLiks <- numeric(length(peez)) # vector of log-likelihoods

for(i in 1:length(peez)) {
  temp <- dbinom(y, 1, peez[i], log=TRUE)
  logLiks[i] <- sum(temp)
}

par(mar=c(3.5,3.4,0.1,0.5)) # set margins
plot(peez, logLiks, type="l", lwd=3, col="blue", las=1, ylim=c(-30, -5), 
  frame.plot=FALSE, xlab="", ylab="")
mtext(text = "log-Likelihood", side=2, line = 2.5)
mtext(text = "p's", side=1, line = 2, font=3)

##########################

#Thinking about the beta distribution

par(mfrow = c(1,2))
par(mar =  c(4,4,0.1,0.5))
curve(dbeta(x, shape1 = 1, shape2 = 1), las=1, ylab = "p(p|a,b)", xlab="p",
  col="blue", lwd=3)

curve(dbeta(x, shape1 = 7, shape2 = 3), las=1, ylab="", xlab="p", 
  col="blue", lwd=3)

################################

# small kappa
par(mfrow=c(1,2))
par(mar=c(4,4,0.1,0.5))

mu <- 0.7 # expected mean proportion of water 
omega <- 0.7 # proportion of water as the mode
kappa <- 5 # low confidence in central measures

# plot in terms of mean:
a <- mu * kappa
b <- (1 - mu) * kappa

curve(dbeta(x, shape1=a, shape2=b),las=1, ylab="p(P|a, b)", xlab="p", 
  col="blue", lwd=2) 
abline(v = mu, col="red", lwd=2, lty=3)
text(0.76, 1.5, expression(mu), cex=2)

# Plotting in terms of mode:
a <- omega * (kappa - 2) + 1
b <- (1 - omega) * (kappa - 2) + 1

curve(dbeta(x, shape1=a, shape2=b),las=1, ylab="", xlab="p", 
  col="blue", lwd=2) 
abline(v=omega, col="red", lwd=2, lty=3)
text(0.7, 1.5, expression(omega), cex=2)

######################################
# large kappa

ar(mfrow=c(1,2))
par(mar=c(4,4,0.1,0.5))

mu <- 0.7 # expected mean proportion of water 
omega <- 0.7 # proportion of water as the mode
kappa <- 50 # low confidence in central measures

# plot in terms of mean:
a <- mu * kappa
b <- (1 - mu) * kappa

curve(dbeta(x, shape1=a, shape2=b),las=1, ylab="P(p|a, b)", xlab="p", 
  col="blue", lwd=2) 
abline(v = mu, col="red", lwd=2, lty=3)
text(0.6, 5.5, expression(mu), cex=2)
# Plotting in terms of mode:
a <- omega * (kappa - 2) + 1
b <- (1 - omega) * (kappa - 2) + 1

curve(dbeta(x, shape1=a, shape2=b),las=1, ylab="", xlab="p", 
  col="blue", lwd=2) 
abline(v=omega, col="red", lwd=2, lty=3)


text(0.6, 5.5, expression(omega), cex=2)