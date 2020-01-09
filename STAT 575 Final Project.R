library(energy)
library(mlbench)
library(mvtnorm)

d <- 2
n <- 30
m <- 1000
alpha <- .1
(cvk <- qnorm(1-(alpha/2), d*(d+2), sqrt(8*d*(d+2)/n)))
(cvs <- qchisq(1-alpha, (d*(d+1)*(d+2))/6))

sigma <- matrix(c(1,0,0,1), ncol=2)
x <- rmvnorm(n, mean=c(0,0), sigma=sigma, method="chol")

as.integer(abs(skew(x,2,20)) >= cvs)
as.integer(abs(kur(x,2,20)) >= cvk)
abs(skew(x,2,5))
abs(kur(x,2,5))

#Algorithm for sample multivariate normal skewness
zbar <- numeric(d)
bz <- matrix(NA, nrow=n, ncol=n)

skew <- function(z,d,n){
  c <- cov(z)
  ic <- solve(c)
  for(k in 1:d){
    zbar[k] <- mean(z[,k])
  }
  for(i in 1:n){
    tZ <- t(z[i,]-zbar)
    cZ <- tZ %*% ic
    for(j in 1:n){
      Z <- z[j,]-zbar
      bz[i,j] <- (cZ %*% Z)^3
    }
  }
  b1 <- sum(bz)/n^2
  return(n*b1/6)
}
(s <- skew(x,d,n))

#Algorithm for sample multivariate normal kurtosis
xbar <- numeric(d)
bx <- numeric(n)

kur <- function(x,d,n){
  c <- cov(x)
  ic <- solve(c)
  for(i in 1:d){
    xbar[i] <- mean(x[,i])
  }
  for(j in 1:n){
    X <- x[j,]-xbar
    tX <- t(x[j,]-xbar)
    bx[j] <- (tX %*% ic %*% X)^2
  }
  return(mean(bx))
}
(k <- kur(x,d,n))

#Multivariate Normality Tests
epsilon <- seq(0,1,.05)
N <- length(epsilon)
y <- numeric(2*n)
Y <- matrix(y, ncol=2)
skewness <- kurtosis <- shapwilk <- energy <- numeric(m)
m.skewness <- m.kurtosis <- m.shapwilk <- m.energy <- numeric(N)

for(i in 1:N){ #for each epsilon
  e <- epsilon[i]
  for(j in 1:m){ #for each replicate
    isigma <- sample(c(1,10), replace=TRUE,
                    size=n, prob=c(1-e,e))
    for(k in 1:n){ #creating the multivariate distribution
      sigma <- matrix(c(isigma[k],0,0,isigma[k]),ncol=2)
      Y[k,1:2] <- rmvnorm(1,mean=c(0,0),
                   sigma=sigma)
    }
    skewness[j] <- as.integer(abs(skew(Y,d,n)) >= cvs)
    kurtosis[j] <- as.integer(abs(kur(Y,d,n)) >= cvk)
    shapwilk[j] <- as.integer(
      shapiro.test(Y)$p.value <= alpha)
    energy[j] <- as.integer(
      mvnorm.etest(Y, R=200)$p.value <= alpha)
  }
  m.skewness[i] <- mean(skewness)
  m.kurtosis[i] <- mean(kurtosis)
  m.shapwilk[i] <- mean(shapwilk)
  m.energy[i] <- mean(energy)
  print(c(epsilon[i], m.skewness[i], m.kurtosis[i], 
          m.shapwilk[i], m.energy[i]))
}

#Results
data.frame(epsilon = epsilon, skewness = m.skewness, kurtosis = m.kurtosis,
           shapiro.wilk = m.shapwilk, energy = m.energy)

#Plot the empirical estimates of power
plot(epsilon, m.skewness, ylim=c(0,1), type="l",
     xlab = "epsilon", ylab = "power")
lines(epsilon, m.kurtosis, lty=2, col="red")
lines(epsilon, m.shapwilk, lty=3, col="blue")
lines(epsilon, m.energy, lty=4, col = "green")
legend("topright", 1, col = c("black", "red", "blue", "green"),
       c("skewness", "kurtosis", "S-W", "energy"),
       lty = c(1,2,3,4), inset = .02)


##Scratch Work##

for(j in 1:m){ 
  isigma <- sample(c(1,10), replace=TRUE,
                   size=n, prob=c(1-e,e))
  for(k in 1:n){ 
    sigma <- matrix(c(isigma[k],0,0,isigma[k]),ncol=2)
    Y[k,1:2] <- rmvnorm(1,mean=c(0,0),
                    sigma=sigma)
  }
}

#Checking means between two multivariate normal dist random generations
z1 <- rmvnorm(20, mean=c(0,0), sigma=matrix(c(1,0,0,1), ncol=2))
data.frame(mean(z1[,1]), mean(z1[,2]))

z2 <- numeric(20)
Z2 <- matrix(z2, ncol=2)
for(k in 1:10){ 
  Z2[k,1:2] <- rmvnorm(1,mean=c(0,0),
                      sigma=matrix(c(1,0,0,1), ncol=2))
}
data.frame(mean(Z2[,1]), mean(Z2[,2]))

#Testing sample function
mean <- sample(c(1,10), replace=TRUE, 
                size=10, prob=c(.5,.5))
sample(1:ncol(mean), 1)
mean2 <- seq(0,10,1)
rnorm(10, mean2, 1)

(isigma <- sample(c(1,10), replace=TRUE, 
               size=1, prob=c(.5,.5)))

sigma <- matrix(c(isigma,0,0,isigma), ncol=2)
(samp <- replicate(5, rmvnorm(1, c(0,0), sigma)))

xbar[1] <- mean(x[,1])
xbar[2] <- mean(x[,2])

c <- cov(x)
ic <- solve(c)
X <- x[2,]-xbar
tX <- t(x[2,]-xbar)
(bx[2] <- (tX %*% ic %*% X))

for(j in 1:n){
  X <- x[j,]-xbar
  tX <- t(x[j,]-xbar)
  bx[j] <- (tX %*% ic %*% X)^2
}

(xbar(x,2))

mean(x[,1])
xbar[1] <- mean(x[,1])


sigma <- matrix(c(1,0,0,1), ncol=2)
z <- rmvnorm(n, mean=c(0,0), sigma=sigma, method="chol")
bZ <- matrix(NA, nrow=n, ncol=n)
c <- cov(z)
ic <- solve(c)
for(k in 1:d){
  zbar[k] <- mean(z[,k])
}
for(i in 1:n){
  tZ <- t(z[i,]-zbar)
  cZ <- tZ %*% ic
  for(j in 1:n){
    Z <- z[j,]-zbar
    bz[i,j] <- cZ %*% Z
    bZ[i,j] <- (bz[i,j])^3
  }
}
(b1 <- sum(bZ)/n^2)
n*b1/6
cvs
