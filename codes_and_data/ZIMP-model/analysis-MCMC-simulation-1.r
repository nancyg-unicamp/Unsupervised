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

 #eta.lm <- lm(y~gamma.post1[1,])$coef
  lambda1.post1[1] <- 0.01 
  lambda2.post1[1] <- 0.05 

 ## initial values for sigma_phi
  sig2phi1.post1[1,] <- 0.05
  sig2phi2.post1[1,] <- 0.05 

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
  
   coefs2 <-  thetaphi$coef
   phi2.post1[r,] <- coefs2 
  } 
   return(list(phi1=phi1.post1,phi2=phi2.post1,mu0=lambda1.post1,mu1=lambda2.post1,gamma=gamma.post1))
}

n <- 100 ## 300,500
L <- 15000
M <- 100

## montando matrix T,S
nk1 <- 7
knots1 <- seq(16,60,length=nk1)
t1 <- seq(16,60,1)
bs.x1 <- bs(t1, knots = knots1)
s.spline1 <- matrix(as.numeric(bs.x1), nr = nrow(bs.x1))   
n.bases1 <- dim(s.spline1)[2]

nk2 <- 5
knots2 <- seq(0,1,length=nk2)
t2 <- seq(0,1,length=30) 
bs.x2 <- bs(t2, knots = knots2)
s.spline2 <- matrix(as.numeric(bs.x2), nr = nrow(bs.x2))   
n.bases2 <- dim(s.spline2)[2]

dim.beta <- (n.bases1+n.bases2)

set.seed(12345)
for (l in 1:M){
  print(l)
  y <- read.table(paste0("y",l,".dat"))
  y <- as.matrix(y)
  
  ## eegs
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
  post.means <- MCMC(y,X,n,L,dim.beta) 
  end_time <- Sys.time()
}

