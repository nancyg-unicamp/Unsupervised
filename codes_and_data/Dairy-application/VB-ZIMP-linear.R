rm(list = ls())
rm(list=ls(all=TRUE))
rm(list=ls())

library(boot)
library(arm)
library(lattice)
library(splines2)
library(splines)
library(statmod)
library(mixtools)
require(mvtnorm)
require(invgamma)
require(monomvn)
library(mvnfast)
library(coda)
library(pracma)
library(tidyverse)
library(fda)
library(latex2exp)
library(readxl)


g <- function(zeta,regre1,x.i,mu2,sig2,qmu.beta2){

  t2 <- dnorm(zeta,mu2,sig2)+10e-10
  tt<-sqrt(t(x.i)%*%x.i)
  regre2 <- sqrt(t(x.i)%*%x.i)*zeta + t(x.i)%*%qmu.beta2
  g <- log(1+exp(regre1)+exp(regre2))*t2
  g<- ifelse(is.infinite(g),10e10,g)
 return(g)
}

g1 <- function(eta,x.i,mu1,sig1,mu2,sig2,qmu.beta1,qmu.beta2){
  
  regre1 <- sqrt(t(x.i)%*%x.i)*eta + t(x.i)%*%qmu.beta1
  
  t1 <- dnorm(eta,mu1,sig1)+10e-1
  f1 <- function(zeta,regre1,x.i,mu2,sig2,qmu.beta2) {g(zeta,regre1,x.i,mu2,sig2,qmu.beta2)}
  g1.1<-integrate(f1,mu2-3*sig2,mu2+3*sig2,regre1,x.i,mu2,sig2,qmu.beta2,stop.on.error = F)$value
 
  g1<- g1.1*t1
  

  return(g1)
}

### integrating in p


aprox.GS <- function(x.i,mu1,sig1,mu2,sig2,qmu.beta1,qmu.beta2){
  
  f2<- function(eta,x.i,mu1,sig1,mu2,sig2,qmu.beta1,qmu.beta2) {g1(eta,x.i,mu1,sig1,mu2,sig2,qmu.beta1,qmu.beta2)}

  g2 <-  integrate(f2,mu1-3*sig1,mu1+3*sig1,x.i,mu1,sig1,mu2,sig2,qmu.beta1,qmu.beta2,stop.on.error = F)$value 
  
  return(g2)
}


elbo.f <- function(n,y,X,a1,b1,a2,b2,alpha.i,zeta1,zeta2,phi1,phi2,dim.beta,Vq.beta1,Vq.beta2,Sigma.beta1,Sigma.beta2,qmu.beta1,qmu.beta2,mu.beta1,mu.beta2){
  ###ELBO
  N <- n
  p0<- ifelse(y==0,1,0)
  p1 <- (-phi1/zeta1+y*(-log(zeta1)+digamma(phi1)))
  p2 <- (-phi2/zeta2+y*(-log(zeta2)+digamma(phi2)))
  
  E.0 <- sum(alpha.i[,1]*p0+alpha.i[,2]*p1+alpha.i[,3]*p2)
  
  #### calculating E1
  E.1 <- 0
  for (i in 1:N){
    x.i <- as.vector(X[i,])
    D <- diag(1,dim.beta)
    D[,1] <- x.i/sqrt(sum(x.i^2))
    gs <- gramSchmidt(D)
    S <- gs$Q
    V.alpha1 <- t(S)%*%Vq.beta1%*%(S)
    V.alpha2 <- t(S)%*%Vq.beta2%*%(S)
    
    ## calculates E(log (1+ exp(Xbeta1)+ exp(Xbeta2)) w.r.t. beta
    mu1<-(t(S)%*%qmu.beta1)[1]
    mu2<-(t(S)%*%qmu.beta2)[1]
    sig1<-sqrt(V.alpha1[1,1])
    sig2<-sqrt(V.alpha2[1,1])
    logq.beta <- aprox.GS(x.i,mu1,sig1,mu2,sig2,qmu.beta1,qmu.beta2)
    
    
    E.1 <- (E.1 + alpha.i[i,1]
            + alpha.i[i,2]*(t(x.i)%*%qmu.beta1)
            + alpha.i[i,3]*(t(x.i)%*%qmu.beta2) - logq.beta)
  }
  
  ## calculating E.2 thru E.5   

  E.2 <- (a1-1)*(-log(zeta1)+digamma(phi1))-b1*phi1/zeta1
  E.3 <- (a2-1)*(-log(zeta2)+digamma(phi2))-b2*phi2/zeta2
  
  E.4 <- (-(dim.beta/2)*log(2*pi) - 0.5*log(det(Sigma.beta1))-0.5*sum(diag(solve(Sigma.beta1)%*%Vq.beta1))
  -0.5*(t(qmu.beta1-mu.beta1)%*%solve(Sigma.beta1)%*%(qmu.beta1-mu.beta1))) 
  
  E.5 <- (-(dim.beta/2)*log(2*pi) - 0.5*log(det(Sigma.beta2))-0.5*sum(diag(solve(Sigma.beta2)%*%Vq.beta2))
  -0.5*(t(qmu.beta2-mu.beta2)%*%solve(Sigma.beta2)%*%(qmu.beta2-mu.beta2))) 
  
  #### Calculating F1 thru F5
  
  z.0<-ifelse(alpha.i[,1]==0,0,alpha.i[,1]*log(alpha.i[,1]))
  z.1<-ifelse(alpha.i[,2]==0,0,alpha.i[,2]*log(alpha.i[,2]))
  z.2<-ifelse(alpha.i[,3]==1,0,(1-alpha.i[,3])*log(1-alpha.i[,3]))
  
  F.1 <-sum(z.0*p0+z.1+z.2)
 
  F.2 <- -lgamma(phi1)+phi1*log(zeta1)+(phi1-1)*(-log(zeta1)+digamma(phi1))-phi1/zeta1
  
  F.3 <- -lgamma(phi2)+phi2*log(zeta2)+(phi2-1)*(-log(zeta2)+digamma(phi2))-phi2/zeta2
  
  F.4 <- -(dim.beta/2)*log(2*pi)-0.5*log(det(Vq.beta1)) - (dim.beta/2)
  F.5 <- -(dim.beta/2)*log(2*pi)-0.5*log(det(Vq.beta2)) - (dim.beta/2)
  
  elbo <- E.0+E.1+E.2+E.3+E.4+E.5-F.1-F.2-F.3-F.4-F.5
  return(elbo)
}

VB <- function(y,X,n,dim.beta){
  
  N <- n

  ## prior hyperparameter
  a1<-  0.01
  a2 <- 0.02
  b1 <- b2 <- 0.01
  mu.beta1 <- rep(0,dim(X)[2])
  Sigma.beta1 <- diag(1^2,dim(X)[2])
  mu.beta2 <- rep(0,dim(X)[2])
  Sigma.beta2 <- diag(1^2,dim(X)[2])

  ### initial values
  zeta1 <- 0.1  
  zeta2 <- 0.1 
  phi1 <-  0.1 
  phi2 <-  0.1
  
  ### beta - initial values
  qmu.beta1 <- rep(0,dim.beta)
  qmu.beta2 <- rep(0,dim.beta)
  qmu.beta1[1] <- 1  ## initial value for intercept
  qmu.beta2[1] <- 1  ## initial value for intercept
  Vq.beta1 <- diag(0.1,dim.beta)
  Vq.beta2 <- diag(0.1,dim.beta)
  
  #### calculating elbo.0
  alpha.i <-matrix(c(rep(0.2,n),rep(0.4,n),rep(0.4,n)),nrow=n,ncol=3)
  
  elbo.ant <- elbo.f(n,y,X,a1,b1,a2,b2,alpha.i,zeta1,zeta2,phi1,phi2,dim.beta,Vq.beta1,Vq.beta2,Sigma.beta1,Sigma.beta2,qmu.beta1,qmu.beta2,mu.beta1,mu.beta2)

  ind <- 0
  alphai.t <- matrix(0,ncol=3,nrow=n) ### variational of gamma_i
  
  ## updating alpha
  for (i in 1:N){
    x.i <- as.vector(X[i,])

        
        alphai.t[i,1] <- ifelse(y[i]==0,1,0)
        alphai.t[i,2] <- y[i]*(-log(zeta1)+digamma(phi1))-phi1/zeta1 + x.i%*%qmu.beta1 ## alpha_1
        alphai.t[i,3] <- y[i]*(-log(zeta2)+digamma(phi2))-phi2/zeta2 + x.i%*%qmu.beta2 ## alpha_2 
  }   
    ind0 <- ifelse(y==0,1,0)

    alpha.i[,1] <- ind0/(ind0+exp(alphai.t[,2])+exp(alphai.t[,3]))
    alpha.i[,2] <- (exp(alphai.t[,2])/(ind0+exp(alphai.t[,2])+exp(alphai.t[,3])))

    alpha.i[,3] <- (exp(alphai.t[,3])/(ind0+exp(alphai.t[,2])+exp(alphai.t[,3])))

  
  ## variational parameters  of lambda0 and lambda1
  phi1 <- a1+sum((1-alpha.i[,2])*y)
  zeta1 <- b1+sum((1-alpha.i[,2]))
  
  phi2 <- a1+sum(alpha.i[,3]*y)
  zeta2 <- b1+sum(alpha.i[,3])
  
  ## prior.df = 1 if cauchy; prior.df = 7 if t-student; prior.df = Inf if Normal
  fit1 <- bayesglm(alpha.i[,2] ~ X-1,family=binomial(link="logit"),prior.scale=2.5, prior.df=Inf)
  fit2 <- bayesglm(alpha.i[,3] ~ X-1,family=binomial(link="logit"),prior.scale=2.5, prior.df=Inf)
  
  regre1 <- X%*%fit1$coefficients
  regre2 <- X%*%fit2$coefficients
  
  pp1 <- invlogit(regre1)
  zz1 <- regre1 + (alpha.i[,2] - pp1)/(pp1*(1-pp1))
  Sigma1.z <- diag((fit1$weights^(-1)),n) ##diag(c((((1+exp(regre))^2)/exp(regre))),n)
  
  Vq.beta1 <- solve(t(X)%*%solve(Sigma1.z)%*%X + solve(Sigma.beta1))
  qmu.beta1 <- Vq.beta1%*%(t(X)%*%solve(Sigma1.z)%*%zz1+solve(Sigma.beta1)%*%mu.beta1)
  
  pp2 <- invlogit(regre2)
  zz2 <- regre2 + (alpha.i[,3] - pp2)/(pp2*(1-pp2))
  Sigma2.z <- diag((fit2$weights^(-1)),n) ##diag(c((((1+exp(regre))^2)/exp(regre))),n)
  
  Vq.beta2 <- solve(t(X)%*%solve(Sigma2.z)%*%X + solve(Sigma.beta2))
  qmu.beta2 <- Vq.beta2%*%(t(X)%*%solve(Sigma2.z)%*%zz2+solve(Sigma.beta2)%*%mu.beta2)
  
  elbo.atual <- elbo.f(n,y,X,a1,b1,a2,b2,alpha.i,zeta1,zeta2,phi1,phi2,dim.beta,Vq.beta1,Vq.beta2,Sigma.beta1,Sigma.beta2,qmu.beta1,qmu.beta2,mu.beta1,mu.beta2)
  
  while (elbo.ant < elbo.atual-0.00001){ 
    elbo.ant <- elbo.atual
    
    ind <- ind+1
    sum.B <- 0
    
    dif <- sum_param.ant <- 10
    while (abs(dif) > 5) { 
     
      ## updating alpha
       for (i in 1:N){
         x.i <- as.vector(X[i,])
         
         alphai.t[i,1] <- ifelse(y[i]==0,1,0)
         alphai.t[i,2] <- y[i]*(-log(zeta1)+digamma(phi1))-phi1/zeta1 + x.i%*%qmu.beta1 ## alpha_1
         alphai.t[i,3] <- y[i]*(-log(zeta2)+digamma(phi2))-phi2/zeta2 + x.i%*%qmu.beta2 ## alpha_2 
       }   
      ind0 <- ifelse(y==0,1,0)
      
      alpha.i[,1] <- ind0/(ind0+exp(alphai.t[,2])+exp(alphai.t[,3]))
      alpha.i[,2] <- (exp(alphai.t[,2])/(ind0+exp(alphai.t[,2])+exp(alphai.t[,3])))
      alpha.i[,3] <- (exp(alphai.t[,3])/(ind0+exp(alphai.t[,2])+exp(alphai.t[,3])))
     
      
      
      ##  variational parameters for lambda0 and lambda1
      phi1 <- a1+sum(alpha.i[,2]*y)
      zeta1 <- b1+sum(alpha.i[,2])
      
      phi2 <- a2+sum(alpha.i[,3]*y)
      zeta2 <- b2+sum(alpha.i[,3])
      
      ## prior.df = 1 if cauchy; prior.df = 7 if t-student; prior.df = Inf if Normal
      fit1 <- bayesglm(alpha.i[,2] ~ X-1,family=binomial(link="logit"),prior.scale=2.5, prior.df=Inf)
      
      regre1 <- X%*%fit1$coefficients
      pp1 <- invlogit(regre1)
      zz1 <- regre1 + (alpha.i[,2] - pp1)/(pp1*(1-pp1))
      Sigma1.z <- diag((fit1$weights^(-1)),n) 
      Vq.beta1 <- solve(t(X)%*%solve(Sigma1.z)%*%X + solve(Sigma.beta1))
      qmu.beta1 <- Vq.beta1%*%(t(X)%*%solve(Sigma1.z)%*%zz1+solve(Sigma.beta1)%*%mu.beta1)
 
      fit2 <- bayesglm(alpha.i[,3] ~ X-1,family=binomial(link="logit"),prior.scale=2.5, prior.df=Inf)
      
      regre2 <- X%*%fit2$coefficients
      pp2 <- invlogit(regre2)
      zz2 <- regre2 + (alpha.i[,3] - pp2)/(pp2*(1-pp2))
      Sigma2.z <- diag((fit2$weights^(-1)),n) 
      Vq.beta2 <- solve(t(X)%*%solve(Sigma2.z)%*%X + solve(Sigma.beta2))
      qmu.beta2 <- Vq.beta2%*%(t(X)%*%solve(Sigma2.z)%*%zz2+solve(Sigma.beta2)%*%mu.beta2)
      
      
      vec.par <- matrix(c(c(alpha.i),c(qmu.beta1),c(qmu.beta2),c(Vq.beta1),c(Vq.beta2),zeta1,zeta2,phi1,phi2))
      sum_param.atual <- sqrt(t(vec.par)%*%vec.par)
      dif <- sum_param.ant-sum_param.atual
      print(dif)
      sum_param.ant <- sum_param.atual
    }
    
    elbo.atual <- elbo.f(n,y,X,a1,b1,a2,b2,alpha.i,zeta1,zeta2,phi1,phi2,dim.beta,Vq.beta1,Vq.beta2,Sigma.beta1,Sigma.beta2,qmu.beta1,qmu.beta2,mu.beta1,mu.beta2)
    
    print(paste0(elbo.ant,elbo.atual))
  }
  return(list(alpha.i=alpha.i,zeta1=zeta1,zeta2=zeta2,phi1=phi1,phi2=phi2,qmu.beta1=qmu.beta1,qmu.beta2=qmu.beta2,Vq.beta1=Vq.beta1,Vq.beta2=Vq.beta2))
} ## end funcion VB


#### Loading Functional covariates 
load(file = 'data/funct-data.rda')

#### Loading Sacalar covariates 
load(file = 'data/covar.total.rda')

#### Loading Response variable
load(file= 'data/y.total.rda')

###### Select Learning and Test Set 


n.learn <- 208
n.test <- length(y.total) - n.learn

set.seed(12345)
aleat <- sample(1:258,n.test,replace=FALSE)

y <- y.total[-aleat]
y.test <- y.total[aleat]


covar <- as.matrix(covar.total[-aleat,],rownames=FALSE)
covar.test<-as.matrix(covar.total[aleat,],rownames=FALSE)

eeg.vacas  <-  eeg.vacas.total[-aleat,,]
eeg.vacas.test <- eeg.vacas.total[aleat,,]




n <- N <- length(y)
n.eeg <- 4
n.freqs <- 31





## matrix S
knots1 <- c(rep(1,4),seq(5,25,length.out=4),rep(31,4))
x1 <- seq(1,n.freqs)
bs.matrix <- spline.des(knots1, x1, ord = 4, derivs=0, outer.ok = TRUE,
                        sparse = FALSE)
s.spline1 <- bs.matrix$design 
n.bases <- dim(s.spline1)[2]

eeg_tnsr <- array(0,dim=c(n,n.eeg,n.freqs))

eeg_tnsr<-eeg.vacas

##  matrix R
R <- array(0,dim=c(n,n.eeg,n.bases))
for (i in 1:n){
  for (j in 1:n.eeg){
    R[i,j,] <- t(s.spline1)%*%eeg_tnsr[i,j,]
  }
}

## R_{i,j,}%*%t(R_{i,j,})
R.sum <- array(0,dim=c(n.eeg,n.bases,n.bases))
for (j in 1:n.eeg){
  for (i in 1:n){
    R.sum[j,,] <- R.sum[j,,]+R[i,j,]%*%t(R[i,j,])
  }}

R.i <- matrix(0,n,n.eeg*n.bases)
for (i in 1:n){
  for (j in 1:n.eeg){
    R.i[i,((j-1)*n.bases+1):((j-1)*n.bases+n.bases)] <- R[i,j,]    
  }
}

X1 <- R.i

## Design matrix Z
Z <- cbind(rep(1,n),covar)

X <- cbind(Z,X1)
dim.beta<- dim(X)[2]

### Start VB
  
  set.seed(12345)
  start_time <- Sys.time()
   post.means <- VB(y,X,n,dim.beta) 
  end_time <- Sys.time()
  tempo <- end_time - start_time

  gamma <- as.matrix(post.means$alpha.i)
  beta1 <- as.matrix(post.means$qmu.beta1)
  beta2 <- as.matrix(post.means$qmu.beta2)
  Vq.beta1 <- as.matrix(post.means$Vq.beta1)
  Vq.beta2 <- as.matrix(post.means$Vq.beta2)
  zeta1 <- post.means$zeta1
  zeta2 <- post.means$zeta2
  phi1 <- post.means$phi1
  phi2 <- post.means$phi2

  

save(beta1,file="beta1-VB-logit-linear-cauchy-elbo-dif.rda")
save(Vq.beta1,file="Vq.beta1-VB-logit-linear-cauchy-elbo-dif.rda")
save(beta2,file="beta2-VB-logit-linear-cauchy-elbo-dif.rda")
save(Vq.beta1,file="Vq.beta2-VB-logit-linear-cauchy-elbo-dif.rda")
save(zeta1,file="zeta1-VB-logit-linear-cauchy-elbo-dif.rda")
save(zeta2,file="zeta2-VB-logit-linear-cauchy-elbo-dif.rda")
save(phi1,file="phi1-VB-logit-linear-cauchy-elbo-dif.rda")
save(phi2,file="phi2-VB-logit-linear-cauchy-elbo-dif.rda")
save(gamma,file="gamma-VB-logit-linear-cauchy-elbo-dif.rda")

save(tempo,file="tempo-VB-logit-linear-cauchy-elbo-dif.rda")

