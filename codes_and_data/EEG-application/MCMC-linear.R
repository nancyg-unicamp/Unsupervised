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
require(truncnorm)
library(truncnorm)
require(invgamma)
require(monomvn)
library(mvnfast)
library(coda)

MCMC <- function(y,Z,R.i,n,n.eeg,n.bases,M){
  ## priors for hyperparameters
  tau0.2 <- 100
  tau1.2 <- 100
  a0 <- 0.01
  b0 <- 0.01
  shape.sigma <- (a0+n/2)
  shape.sigma2phi <- a0 + (n.eeg*n.bases)/2 
  sig.muphi <- 10e-6 
  T.star <- diag(10e-6,3)
  T.til <- solve(solve(T.star)+t(Z)%*%Z)
  theta.star <- rep(0,dim(Z)[2])
  T.til <- solve(solve(T.star)+t(Z)%*%Z)
  
  ## initial values gamma
  gamma.post1 <- matrix(0,nrow=M,ncol=n)
  normalEM <- normalmixEM(y,2)
  class <- apply(normalEM$posterior,1,which.max)
  gamma.post1[1,] <- kmeans(y, centers = sort(kmeans(y, centers = 2)$centers))$cluster-1 #class-1
  
  ## initial values phi
  phi.post1 <- array(0,dim=c(M,n.eeg*n.bases)) 
  
  ## initial values theta
  theta.post1 <- matrix(0,nrow=M,ncol=dim(Z)[2])
  theta.post1[1,] <- glm(gamma.post1[1,]~Z-1,family=binomial(link=logit))$coef
  
  ## initial values for etas, sig2, sig2phi
  mu0.post1 <- mu1.post1 <- sig2.post1 <- matrix(0,nrow=M)
  
  #eta.lm <- lm(y~gamma.post1[1,])$coef
  mu0.post1[1] <- -1 
  mu1.post1[1] <- 9 
  
  ## initial values for sigma
  sig2.post1[1,] <- 20 
  
  ## initial value for mu_phi
  
  muphi.post1 <- matrix(0,nrow=M,n.bases*n.eeg)
  muphi.post1[1,] <- 0 
  
  set.seed(12345)
  for (r in  2:M){
    
    ## update gamma
    dgauss0 <- dnorm((y-mu0.post1[r-1])/sqrt(sig2.post1[r-1,]))
    dgauss1 <- dnorm((y-mu1.post1[r-1])/sqrt(sig2.post1[r-1,]))
    regre <- Z%*%(theta.post1[r-1,])+R.i%*%as.vector(t(phi.post1[r-1,]))
    p <- exp(regre)/(1+exp(regre))
    ratio <- (dgauss1*p)/(dgauss0*(1-p)+dgauss1*p)
    gamma.post1[r,] <- rbinom(n,size=1,prob=ratio)
    
    ## update mu0
    sigeta <- (sum(gamma.post1[r,]==0)/sig2.post1[r-1,]+1/tau0.2)^(-1)
    mueta <- (sigeta/sig2.post1[r-1,])*sum(y[gamma.post1[r,]==0])
    mu0.post1[r] <- rnorm(1,mueta,sqrt(sigeta))
    
    ### update mu1
    sigeta <- (sum(gamma.post1[r,]==1)/sig2.post1[r-1,]+1/tau1.2)^(-1)
    mueta <- (sigeta/sig2.post1[r-1,])*sum(y[gamma.post1[r,]==1])
    mu1.post1[r] <- rtruncnorm(1, a=mu0.post1[r], b=Inf, mean = mueta, sd = sqrt(sigeta))
    
    ### update sigma
    mu.y <- sum((y[gamma.post1[r,]==0]-mu0.post1[r,])^2)+sum((y[gamma.post1[r,]==1]-mu1.post1[r,])^2)
    rate <- b0+mu.y/2
    sig2.post1[r] <- rinvgamma(1, shape=shape.sigma, rate=rate)
    
    ### update thetas e phis
    ## prior.df = 1 is Cauchy
    ## prior.df = Inf is Normal (default)
    ## Otherwise prior is t-Student
    ## link: logit/probit
    
    thetaphi <- bayesglm(gamma.post1[r,] ~ Z[,1]+Z[,2]+Z[,3]+R.i-1, 
                         family=binomial(link="logit"),
                         prior.scale=2.5, prior.df=Inf)
    
    phi.post1[r,] <- thetaphi$coef[(dim(Z)[2]+1):(dim(theta.post1)[2]+dim(phi.post1)[2])]
    theta.post1[r,] <- thetaphi$coef[1:dim(Z)[2]]
  }
  return(list(sig2=sig2.post1,gamma=gamma.post1,mu0=mu0.post1,mu1=mu1.post1,theta=theta.post1,phi=phi.post1))
}

n <- 96
n.eeg <- 14
n.freqs <- 45

M <- 15000
ML <- (M/2)+1

#load eeg data
load(file = "data/eeg-data.rda")
eeg_tnsr <- eeg

# load scalar data - response+scalar covariates
scalars <- read.table("data/scalars.dat")
y <- scalars$HAMD
Z <- cbind(rep(1,n),scalars$chronicity,scalars$sex)

## normalizing eegs 
eeg.tnsr2 <- array(0,dim=c(n.eeg,n.freqs,n))
for (i in 1:n)
  eeg.tnsr2[,,i] <- (eeg_tnsr[,,i]-mean(eeg_tnsr[,,i]))/sd(eeg_tnsr[,,i])
eeg_tnsr <- eeg.tnsr2

## S matrix
knots1 <- c(16,16,16,16,18,20,22,25,30,35,40,45,50,55,60,60,60,60) ##seq(17,59)
x1 <- seq(16,60)
bs.matrix <- spline.des(knots1, x1, ord = 4, derivs=0, outer.ok = TRUE,
                        sparse = FALSE)
s.spline1 <- bs.matrix$design
n.bases <- dim(s.spline1)[2]

## R matrix
R <- array(0,dim=c(n,n.eeg,n.bases))
for (i in 1:n){
  for (j in 1:n.eeg){
    R[i,j,] <- t(s.spline1)%*%eeg_tnsr[j,,i]
  }
}

R.sum <- array(0,dim=c(n.eeg,n.bases,n.bases))
for (j in 1:n.eeg){
  for (i in 1:n){
    R.sum[j,,] <- R.sum[j,,]+R[i,j,]%*%t(R[i,j,])
  }}

R.i <- matrix(0,96,n.eeg*n.bases)
for (i in 1:n){
  for (j in 1:n.eeg){
    R.i[i,((j-1)*n.bases+1):((j-1)*n.bases+n.bases)] <- R[i,j,]    
  }
}

R.isum <- matrix(0,n.eeg*n.bases,n.eeg*n.bases)
for (i in 1:n){
  R.isum <- R.isum+R.i[i,]%*%t(R.i[i,])
}

post.means <- MCMC(y,Z,R.i,n,n.eeg,n.bases,M)
