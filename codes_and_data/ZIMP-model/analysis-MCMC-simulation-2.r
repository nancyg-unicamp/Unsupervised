rm(list = ls())
rm(list=ls(all=TRUE))
rm(list=ls(all.names=TRUE))

library(arm)
library(lattice)
library(splines2)
library(splines)
library(statmod)
library(mixtools)
require(mvtnorm)
# Library for sampling from Truncated Normal distribution
require(truncnorm)
require(invgamma)
require(monomvn)
library(mvnfast)
library(coda)
library(tidyverse)
library(fda)
library(latex2exp)

MCMC <- function(y,X,n,L,n.bases){
  ## prior hyperparameter 
   a0.1 <- 0.01 ## gamma shape class 1 
   b0.1 <- 0.01 ## gamma rate

   a0.2 <- 0.03 ## gamma parameters class 2
   b0.2 <- 0.01

  ## initial values gamma

   gamma.post1 <- array(0,dim=c(L,3,n))
   gamma.post1[1,,] <- rmultinom(n,1,c(1/3,1/3,1/3))

  ## initial values phi
   phi1.post1 <- array(0,dim=c(L,n.bases)) 
   phi2.post1 <- array(0,dim=c(L,n.bases))

 ## initial values for etas, sig2, sig2phi
   sig2phi1.post1 <-  matrix(0,nrow=L)
   sig2phi2.post1 <-  matrix(0,nrow=L)

 ## mus = lambdas
  lambda1.post1 <- lambda2.post1 <-  rep(0,L)
  ratio <- array(0,dim=c(L,3,n))

  lambda1.post1[1] <- 0.01 
  lambda2.post1[1] <- 0.05 

 ## initial values for sigma_phi
  sig2phi1.post1[1,] <- 0.05 #var(phi1.post1[1,])
  sig2phi2.post1[1,] <- 0.05 #var(phi2.post1[1,])

 ## initial value for mu_phi

  muphi1.post1 <- matrix(0,nrow=L,n.bases)
  muphi1.post1[1,] <- 0 

  muphi2.post1 <- matrix(0,nrow=L,n.bases)
  muphi2.post1[1,] <- 0 

  set.seed(292102)
  start.time <- Sys.time()
  for (r in  2:L){
  
    regre1 <- X%*%as.vector(t(phi1.post1[r-1,]))
    regre2 <- X%*%as.vector(t(phi2.post1[r-1,]))
  
    ratio.1 <- exp(-lambda1.post1[r-1])*(lambda1.post1[r-1]^y)*exp(regre1)
    ratio.2 <- exp(-lambda2.post1[r-1])*(lambda2.post1[r-1]^y)*exp(regre2)
  
    ratio[r,1,] <- ifelse(y==0,1,0)/(ifelse(y==0,1,0) + ratio.1 + ratio.2)
    ratio[r,2,] <- ratio.1/(ifelse(y==0,1,0) + ratio.1 + ratio.2)
    ratio[r,3,] <- ratio.2/(ifelse(y==0,1,0) + ratio.1 + ratio.2)
  
    gamma.post1[r,,] <- sapply(1:n,function(x){rmultinom(1,size=1,prob=c(ratio[r,1,x],ratio[r,2,x],ratio[r,3,x]))})
  
    ## update lambda0 = mu0
    a1.1 <- a0.1 + sum(y*gamma.post1[r,2,])
    b1.1 <- b0.1 + sum(gamma.post1[r,2,])
    lambda1.post1[r] <- rgamma(1,a1.1,rate=b1.1)
  
    ## update lambda1 = mu1
    a1.2 <- a0.2 + sum(y*gamma.post1[r,3,])
    b1.2 <- b0.2 + sum(gamma.post1[r,3,])
    lambda2.post1[r] <- rgamma(1,a1.2,rate=b1.2)
  
    ### update thetas e phis
  
    thetaphi1 <- bayesglm (gamma.post1[r,2,] ~ X-1, 
                   family=binomial(link="logit"),
                  prior.scale=2.5, prior.df=Inf)

    coefs1 <-  thetaphi1$coef
    phi1.post1[r,] <- coefs1  
  
   thetaphi <- bayesglm (gamma.post1[r,3,] ~ X-1, 
                         family=binomial(link="logit"),
                         prior.scale=2.5, prior.df=Inf)
  
   coefs <-  thetaphi$coef
   phi2.post1[r,] <- coefs 
 
   } 
   return(list(phi1=phi1.post1,phi2=phi2.post1,mu0=lambda1.post1,mu1=lambda2.post1,gamma=gamma.post1))
}

n <- 150
L <- 15000
M <- 100

## montando matrix T,S
set.seed(12345)

nk1 <- 13
nk2 <- 13
n <- 150

t <- seq(0,10,length.out=256)
Z1 <- matrix(rnorm(n*nk1,-0.1,1),nrow=n,ncol=nk1)
U1 <- matrix(runif(nk1*nk1,0,1),nrow=nk1,ncol=nk1)
C1 <- Z1%*%U1

s <- seq(0,1,length.out=256)
Z2 <- matrix(rnorm(n*nk2,0,1),nrow=n,ncol=nk2)
U2 <- matrix(runif(nk2*nk2,0,1),nrow=nk2,ncol=nk2)
C2 <- Z2%*%U2

mknots1 <- seq(0,10,length=nk1-1)
mknots2 <- seq(0,1,length=nk2-1)

#interior knots
mknots1 <- mknots1[-c(1,length(mknots1))]
mknots2 <- mknots2[-c(1,length(mknots2))]
mknots2 <- mknots2[-c(1,length(mknots2))]


#regressao b-splines knots equi-spaced''

# matrix b-splines 
B1 <- bs(t, knots = mknots1)
B2 <- bs(s, knots = mknots2)

# functional covariates
X1 <- C1%*%t(B1)/10
X2 <- C2%*%t(B1)

ind <- rbinom(n,1,0.3)

for (i in 1:n){
  X1[i,] <- X1[i,] + ind[i]*dnorm(t,6,2)
  X2[i,] <- X2[i,] + ind[i]*dnorm(s,0.5,1)
}

eeg1 <- X1
eeg2 <- X2

#interior knots
mknots1 <- seq(0,10,length = 7)
mknots2 <- seq(0,1,length = 5)

mknots1 <- mknots1[-c(1,length(mknots1))]
mknots2 <- mknots2[-c(1,length(mknots2))]

#regressao b-splines knots equi-espacados''

# matrix b-splines 
B1 <- bs(t, knots = mknots1)
B2 <- bs(s, knots = mknots2)
n.bases1 <- dim(B1)[2]
n.bases2 <- dim(B2)[2]

dim.beta <- (n.bases1+n.bases2)

n.bases <-n.bases1+n.bases2

## montando matrix R

R1 <- array(0,dim=c(n,1,n.bases1))
for (i in 1:n){
  for (k in 1:n.bases1){
    R1[i,1,k] <- B1[,k]%*%eeg1[i,]
  }
}

R2 <- array(0,dim=c(n,1,n.bases2))
for (i in 1:n){
  for (k in 1:n.bases2){
    R2[i,1,k] <- B2[,k]%*%eeg2[i,]
  }  
}

R.i <- matrix(0,n,n.bases)
for (i in 1:n){
  R.i[i,1:n.bases1] <- R1[i,1,]    
  R.i[i,(n.bases1+1):n.bases] <- R2[i,1,]
}
X <- R.i

for (l in 1:M){
  print(l)
  y <- read.table(paste0("y",l,".dat"))
  y <- as.matrix(y)
  
  start_time <- Sys.time()
    post.means <- MCMC(y,X,n,L,n.bases) 
  end_time <- Sys.time()
}
