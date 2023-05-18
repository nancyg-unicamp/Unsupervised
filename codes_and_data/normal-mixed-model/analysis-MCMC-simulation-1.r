rm(list = ls())
rm(list=ls(all=TRUE))
rm(list=ls(all.names=TRUE))

library(mixtools)
library(arm)
library(lattice)
library(splines2)
library(splines)
library(statmod)
library(mixtools)
require(mvtnorm)
require(truncnorm)
require(invgamma)
require(monomvn)
library(coda)

MCMC <- function(y,X,n,L,n.bases){
  
  ## prior hyperparameter
  tau0.2 <- 100
  tau1.2 <- 100
  a0 <- 0.01
  b0 <- 0.01
  shape.sigma <- (a0+n/2)
  Sigma.beta <- diag(2.5^2,dim(X)[2])
  mu.beta <- rep(0,dim(X)[2])
  
  set.seed(292102)
  ## initial values gamma
  kmeans <- kmeans(y,2)$cluster
  
  gamma.post1 <- matrix(0,nrow=L,ncol=n)
  gamma.post1[1,] <- rbinom(n,1,0.5)
    
  ## initial values phi
  phi.post1 <- array(0,dim=c(L,dim(X)[2])) 
  
  ## initial values for etas, sig2, sig2phi
  mu0.post1 <- mu1.post1 <- sig2.post1 <- sig2phi.post1 <-  matrix(0,nrow=L)
  
  mu0.post1[1,1] <- 0
  mu1.post1[1,1] <- 3
  
  ## initial values for sigma
  sig2.post1[1,1] <- var(y)
  
  ## initial value for mu_phi
  muphi.post1 <- matrix(0,nrow=L,dim(X)[2])
  
  set.seed(12345)
  start.time <- Sys.time()
  for (r in  2:L){
    
    ## update gamma 
    dgauss0 <- dnorm((y-mu0.post1[r-1])/sqrt(sig2.post1[r-1,]))
    dgauss1 <- dnorm((y-mu1.post1[r-1])/sqrt(sig2.post1[r-1,]))
    regre <- X%*%phi.post1[r-1,]
    p <- exp(regre)/(1+exp(regre))
    ratio <- (dgauss1*p)/(dgauss0*(1-p)+dgauss1*p)
    gamma.post1[r,] <- rbinom(n,size=1,prob=ratio)
    
    ## update eta0
    sigeta <- (sum(gamma.post1[r,]==0)/sig2.post1[r-1,]+1/tau0.2)^(-1)
    mueta <- (sigeta/sig2.post1[r-1,])*sum(y[gamma.post1[r,]==0])
    mu0.post1[r,1] <- rnorm(1,mueta,sqrt(sigeta))
    
    ### update eta1
    sigeta <- (sum(gamma.post1[r,]==1)/sig2.post1[r-1,]+1/tau1.2)^(-1)
    mueta <- (sigeta/sig2.post1[r-1,])*sum(y[gamma.post1[r,]==1])
    mu1.post1[r,1] <- rtruncnorm(1, a=mu0.post1[r], b=Inf, mean = mueta, sd = sqrt(sigeta))
    
    ### update sigma
    mu.y <- sum((y-mu0.post1[r,]*(1-gamma.post1[r,])-mu1.post1[r,]*gamma.post1[r,])^2)
    rate <- b0+mu.y/2
    sig2.post1[r,1] <- rinvgamma(1, shape=shape.sigma, rate=rate)
    
    ### update thetas e phis
    ## prior.df = 1 is Cauchy
    ## prior.df = Inf is Normal (default)
    ## Otherwise prior is t-Student
    ## link: logit/probit
    
    fit <- bayesglm (gamma.post1[r,] ~ X-1, 
                     family=binomial(link="logit"),
                     prior.scale=2.5, prior.df=Inf)
    phi.post1[r,] <- fit$coefficients
  }
  return(list(theta=phi.post1,mu0=mu0.post1,mu1=mu1.post1,sig2=sig2.post1,gamma=gamma.post1))
} ## end funcion MCMC

n <- 100 ## 300, 500
a0 <- b0 <- 0.01
tau0.2 <- 100
L <- 15000
M <- 100

## T,S matrices
nk1 <- 7
knots1 <- seq(16,60,length=nk1)
t <- seq(16,60,1)
bs.x1 <- bs(t, knots = knots1)
s.spline1 <- matrix(as.numeric(bs.x1), nr = nrow(bs.x1)) 
n.bases1 <- dim(s.spline1)[2]

nk2 <- 5
knots2 <- seq(0,1,length=nk2)
s <- seq(0,1,length=30) 
bs.x2 <- bs(s, knots = knots2)
s.spline2 <- matrix(as.numeric(bs.x2), nr = nrow(bs.x2)) 
n.bases2 <- dim(s.spline2)[2]

for (l in 1:M){
  print(l)
  y <- read.table(paste0("y",l,".dat"))
  y <- as.matrix(y)
  
  ## R matrix
  eeg1 <- read.table(paste0("eeg1-",l,".dat"))
  eeg1 <- as.matrix(eeg1)
  
  eeg2 <- read.table(paste0("eeg2-",l,".dat"))
  eeg2 <- as.matrix(eeg2)
  
  R1 <- array(0,dim=c(n,1,n.bases1))
  for (i in 1:n){
    for (k in 1:n.bases1){
      R1[i,1,k] <- s.spline1[,k]%*%eeg1[i,]
    }
  }
  
  R2 <- array(0,dim=c(n,1,n.bases2))
  for (i in 1:n){
    for (k in 1:n.bases2){
      R2[i,1,k] <- s.spline2[,k]%*%eeg2[i,]
    }  
  }
  
  R.i <- matrix(0,n,n.bases1+n.bases2)
  for (i in 1:n){
    R.i[i,1:n.bases1] <- R1[i,1,]    
    R.i[i,(n.bases1+1):(n.bases1+n.bases2)] <- R2[i,1,]
  }
  X <- R.i
  
  start_time <- Sys.time()
   post.means <- MCMC(y,X,n,L,n.bases) 
  end_time <- Sys.time()
 }

