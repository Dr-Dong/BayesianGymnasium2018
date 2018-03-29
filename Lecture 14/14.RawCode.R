
# Run this to install rethinking package
install.packages(c('devtools','coda','mvtnorm','loo'))
library(devtools)
install_github("rmcelreath/rethinking")


######################
# Make up some data

library(rethinking)

x <- c(37, 35.5, 34.5, 41.5, 55.5, 61, 53.5)
y <- c(438, 452, 612, 521, 752, 871, 1350)
par(mar=c(3,3.2,0.1,0.5))
plot(x,y, pch=16, col="blue", las=1)
######################


######################
#Create models of increasing complexity

m1 <- lm(y ~ x) # straight line
m2 <- lm(y ~ x + I(x^2)) # quadratic
m3 <- lm(y ~ x + I(x^2) + I(x^3)) # cubic
m4 <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4)) # 4th degree
m5 <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5)) # 5th degree
m6 <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6)) # 6th degree

######################


######################
# Predict new data

# Sequence from 20 to 65, by 0.1
x1 <- seq(20,65,0.1) 
y2 <- predict(m2, new=data.frame(x=x1))
y3 <- predict(m3, new=data.frame(x=x1))
y4 <- predict(m4, new=data.frame(x=x1))
y5 <- predict(m5, new=data.frame(x=x1))
y6 <- predict(m6, new=data.frame(x=x1))
######################


######################
# Plot results in 6-panel figure

par(mfrow=c(3, 2))
par(mar=c(3,3.2,0.1,0.5))

# straight line
plot(x,y, pch=16, col="blue", las=1)
abline(m1)
text(40,1300,expression(paste(R^2, "= 0.48")), cex=1.5)

# Quadratic
plot(x,y, pch=16, col="blue", las=1)
lines(x1,y2)
text(40,1300,expression(paste(R^2, "= 0.53")), cex=1.5)

# Cubic
plot(x,y, pch=16, col="blue", las=1)
lines(x1, y3)
text(40,1300,expression(paste(R^2, "= 0.67")), cex=1.5)

# x^4
plot(x,y, pch=16, col="blue", ylim=c(300,1500), las=1)
lines(x1,y4)
text(40,1400,expression(paste(R^2, "= 0.81")), cex=1.5)

# x^5
plot(x,y, pch=16, col="blue", ylim=c(0,1800), las=1)
lines(x1,y5)
text(40,1700,expression(paste(R^2, "= 0.99")), cex=1.5)

# x^6
plot(x,y, pch=16, col="blue", ylim=c(0,1500), las=1)
lines(x1,y6)
text(40,1400,expression(paste(R^2, "= 1")), cex=1.5)

######################

######################
# Cross-validate and plot results for linear and 
# polynomials

par(mar=c(3,3.2,0.1,0.5))
par(mfrow=c(1,2))

# Polynomial regression
plot(x,y, pch=16, col="blue", ylim=c(-500,2000), las=1)
for(i in 1:length(x)){
  x1 <- x[-i]; y1 <- y[-i]
  m <- lm(y1 ~ x1 + I(x1^2) + I(x1^3) + I(x1^4) + I(x1^5))
  xz <- seq(20,65,0.1)
  yz <- predict(m, new=data.frame(x1=xz))
  lines(xz, yz)
}

# linear regression
plot(x,y, pch=16, col="blue", ylim=c(-500,2000), las=1)
for(i in 1:length(x)){
  x1 <- x[-i]; y1 <- y[-i]
  m <- lm(y1 ~ x1)
  xz <- seq(20,65,0.1)
  yz <- predict(m, new=data.frame(x1=xz))
  lines(xz, yz)
}
######################

######################
# Plot K-L divergence

P <- c(0.3,0.7) # true target probs
ws <- seq(0.01,0.99, by=0.001) # probs for wrinkly spreader
Q <- cbind(ws, rev(ws)) # matrix of candidate values

KL <- function (P, Q) { # function to calculate divergence
  kl <- vector(length=length(P)) 
  for(i in 1:length(P)) {
    kl[i] <- P[i]*(log(P[i])-log(Q[i]))
  }
  return(sum(kl))
}

kl <- apply(Q, 1, KL, P=P) # calculate divergence for each row of Q
par(mar=c(3,3.2,0.1,0.5))
plot(Q[,1], kl, las=1, col="cornflowerblue", type="l", lwd=3)
mtext(text = expression(q[ws]), side=1, line=2)
mtext(text = expression(D[KL]), side=2, line=2.1)

abline(v=0.3, lty=3, lwd=2) # where Q=P
text(0.45, 1.5, "Q = P")
######################


######################
#
N <- 100
kseq <- 1:5
nSims <- 500

# If you have a mac set this
cores <- 4

dev <- sapply(kseq, function (k) {
  #  print(k);
  # If you have a mac use this
  r <- mcreplicate( nSims, sim.train.test(N=N, k=k), mc.cores=cores)
  # If you don't have a mac use this
  #  r <- replicate(nSims, sim.train.test( N=N, k=k ));
  c( mean(r[1,]), mean(r[2,]), sd(r[1,]), sd(r[2,]) )
} )
```

```{r, fig.height=2.75, fig.width=5}
par(mar=c(3,3.2,0.1,0.5))
plot(1:5, dev[1,], ylim=c(min(dev[1:2,])-5, max(dev[1:2,]+10)), 
     xlim=c(1,5.1), las=1, ann=FALSE, pch=16, col="blue")

mtext(text = "# parameters", side=1, line = 2, cex=1)
mtext(text = "deviance", side=2, line = 2.2, cex=1)

points((1:5)+0.1, dev[2,])
for (i in kseq) {
  IN <- dev[1,i] + c(-1,1)*dev[3,i]
  OUT <- dev[2,i] + c(-1,1)*dev[4,i]
  lines(c(i,i), IN, col="blue")
  lines(c(i,i)+0.1, OUT)  
}

text(1.8, 55, "in", col="blue")
text(2.3, 59, "out")
######################