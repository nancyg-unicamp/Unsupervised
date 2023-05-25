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
require(truncnorm)
require(invgamma)
require(monomvn)
library(mvnfast)
library(coda)
library(pracma)

### integrating in p
aprox.GS1 <- function(x.i,qmu.beta,V.alpha,S){     #,rel.tol = 1e-10
  z <- integrate(log.GS1,-Inf,Inf,x.i,qmu.beta,V.alpha,S,stop.on.error = F)
  return(z$value)
}

log.GS1 <- function(alpha,x.i,qmu.beta,V.alpha,S){
  regre <- sqrt(t(x.i)%*%x.i)*alpha + t(x.i)%*%qmu.beta
  aprox <- log(invlogit(regre)+0.00000000000000000000001)*(dnorm(alpha,(t(S)%*%qmu.beta)[1],sqrt(V.alpha[1,1]))+10e-10) 
  return(aprox)
}

aprox.GS2 <- function(x.i,qmu.beta,V.alpha,S){     
  z <- integrate(log.GS2,-Inf,Inf,x.i,qmu.beta,V.alpha,S,stop.on.error = F)
  return(z$value)
}

log.GS2 <- function(alpha,x.i,qmu.beta,V.alpha,S){
  regre <- sqrt(t(x.i)%*%x.i)*alpha + t(x.i)%*%qmu.beta
  aprox <- log(1-invlogit(regre)+0.00000000000000000000001)*(dnorm(alpha,(t(S)%*%qmu.beta)[1],sqrt(V.alpha[1,1]))+10e-10) #,0)
  return(aprox)
}

elbo.f <- function(a0,b0,tau0.2,tau1.2,alpha.i,y,m,dim.beta,X,Vq.beta,s2,n,Sigma.beta,qmu.beta,mu.beta){
  ###ELBO
  Aq.sig2 <- a0+n/2
  Bq.sig2 <- 0.5*sum(alpha.i*((y-m[2])^2+s2[2]))+0.5*sum((1-alpha.i)*((y-m[1])^2+s2[1]))
  E.0 <- (-(n/2)*log(2*pi)) - (Aq.sig2/(b0+Bq.sig2))*Bq.sig2 + (n/2)*(log(b0+Bq.sig2) - digamma(Aq.sig2))  
  
  #### calculating E1
  E.1 <- 0
  for (i in 1:n){
    x.i <- as.vector(X[i,])
    D <- diag(1,dim.beta)
    D[,1] <- x.i/sqrt(sum(x.i^2))
    gs <- gramSchmidt(D)
    S <- gs$Q
    V.alpha <- t(S)%*%Vq.beta%*%(S)
    
    ## calculates E(log invlogit(Xbeta)) w.r.t. beta
    logq.beta1 <- aprox.GS1(x.i,qmu.beta,V.alpha,S)
    
    ## calculates E(log 1-invlogit(Xbeta)) w.r.t. beta
    logq.beta2 <- aprox.GS2(x.i,qmu.beta,V.alpha,S)
    
    E.1 <- E.1 + alpha.i[i]*logq.beta1 + (1-alpha.i[i])*logq.beta2  
  }
  
  ## calculating E.2 a E.5 and Fs  
  E.2 <- -0.5*log(2*pi*(tau0.2)) - 0.5*(1/tau0.2)*(m[1]^2+s2[1])
  E.3 <- -0.5*log(2*pi*(tau1.2)) - 0.5*(1/tau1.2)*(m[2]^2+s2[2])
  
  E.4 <- a0*log(b0) - lgamma(a0) + (a0+1)*(log(b0+Bq.sig2) - digamma(a0+n/2)) - b0*(a0+n/2)*(1/(b0+Bq.sig2))

  E.5 <- -(dim.beta/2)*log(2*pi) - 0.5*log(det(Sigma.beta))-0.5*sum(diag(solve(Sigma.beta)%*%Vq.beta))
  -0.5*(t(qmu.beta-mu.beta)%*%solve(Sigma.beta)%*%(qmu.beta-mu.beta)) 
  
  F.1 <-sum(alpha.i*log(alpha.i)+(1-alpha.i)*log(1-alpha.i))
  
  F.2 <- -0.5*log(2*pi)-0.5-0.5*log(s2[1])
  
  F.3 <- -0.5*log(2*pi)-0.5-0.5*log(s2[2])
  
  F.4 <- (a0+n/2)*log(b0+Bq.sig2)-lgamma(a0+n/2)+(a0+n/2+1)*(log(b0+Bq.sig2)-digamma(a0+n/2))-(a0+n/2)

    F.5 <- -(dim.beta/2)*log(2*pi)-0.5*log(det(Vq.beta)) - (dim.beta/2)
  
  elbo <- E.0+E.1+E.2+E.3+E.4+E.5-F.1-F.2-F.3-F.4-F.5
  return(elbo)
}

VB <- function(y,X,n,dim.beta){
  ## prior hyperparameter
  tau0.2 <- 100 # variance for m_0
  tau1.2 <- 100 # variance for m_1
  a0 <- A <- 0.01
  b0 <- B <- 0.01
  shape.sigma <- (a0+n/2)
  mu.beta <- rep(0,dim(X)[2])
  Sigma.beta <- diag(2.5^2,dim(X)[2])
  mu.mu <- c(0,0)
  
  ### initial values
  
  ###for s2
  s2.0 <- 20 #var(y)
  
  ### for mean (mu.mu,sig2.mu,mu.beta,Sigma.beta)
  means = normalmixEM(y,2)
  sort.means <- sort(means$mu,index.return=T)
  
  v <- means$posterior[,2] #round(means$posterior[,2])
  m <- means$mu
  s2 <- (means$sigma)^2
  
  if (means$mu[1]>means$mu[2]){
    m <- means$mu[sort.means$ix]
    s2 <- (means$sigma[sort.means$ix])^2
    v <- 1-means$posterior[,2] #round(1-means$posterior[,2])
  }   
  
  ### beta - initial values
  qmu.beta <- rep(0,dim.beta)
  qmu.beta[1:3] <- c(-0.1,-0.3,-0.4)
  Vq.beta <- diag(0.1,dim.beta)
  m <- c(-1,9)
  
  ## initial values gamma
  kmeans <- kmeans(y,2)$cluster
  
  normalEM <- normalmixEM(y,2)
  class <- apply(normalEM$posterior,1,which.max)
  alpha.i <- normalEM$posterior[,2]
  
  #### elbo.0
  
  elbo.ant <- elbo.f(a0,b0,tau0.2,tau1.2,alpha.i,y,m,dim.beta,X,Vq.beta,s2,n,Sigma.beta,qmu.beta,mu.beta)
  qsig2.inv <- 1/var(y) ### 
  
  ind <- sig2 <- 0
  alphai.t <- matrix(0,ncol=2,nrow=n) ### posterior probability of c_i
  
  sum.B <- 0
  ## updating alpha
  for (i in 1:n){
    x.i <- as.vector(X[i,])
    D <- diag(1,dim.beta)
    D[,1] <- x.i/sqrt(sum(x.i^2))
    gs <- gramSchmidt(D)
    S <- gs$Q
    V.alpha <- t(S)%*%Vq.beta%*%(S)
    
    ## calculates E(log \Phi{-1}(Xbeta)) w.r.t. beta
    logq.beta1 <- aprox.GS1(x.i,qmu.beta,V.alpha,S)
    
    ## calculates E(log 1-\Phi{-1}(Xbeta)) w.r.t. beta
    logq.beta2 <- aprox.GS2(x.i,qmu.beta,V.alpha,S)
    
    alphai.t[i,1] <- exp(logq.beta2-0.5*qsig2.inv*((y[i]-m[1])^2+s2[1])/2) ## alpha_0
    alphai.t[i,2] <- exp(logq.beta1-0.5*qsig2.inv*((y[i]-m[2])^2+s2[2])/2) ## alpha_1
    
    alpha.i[i] <- alphai.t[i,2]/(alphai.t[i,1]+alphai.t[i,2])
    sum.B <- sum.B + (((y[i]-m[2])^2+s2[2])*alpha.i[i]) + (((y[i]-m[1])^2+s2[1])*(1-alpha.i[i]))
  }
  
  ## A and B come from prior for sigma^2
  Bqsig2 <- B + 0.5*sum.B
  qsig2.inv <-  (A+n/2)/Bqsig2 ### E[1/sigma^2]
  
  s2[1] <- 1/(qsig2.inv*sum((1-alpha.i))+1/tau0.2)
  m[1] <- (qsig2.inv*(t((1-alpha.i))%*%y)+mu.mu[1]/tau0.2)*s2[1]
  
  s2[2] <- 1/(qsig2.inv*sum((alpha.i))+1/tau0.2)
  m[2] <- (qsig2.inv*(t((alpha.i))%*%y)+mu.mu[2]/tau1.2)*s2[2]
  
  
  ## prior.df = 1 if cauchy; 
  ## prior.df = 7 if t-student; 
  ## prior.df = Inf if Normal
  ## link = logit, probit
  fit <- bayesglm(alpha.i ~ X-1,
                  family=binomial(link="logit"),prior.scale=2.5, prior.df=Inf)
  
  regre <- X%*%fit$coefficients
  pp <- invlogit(regre)
  zz<- regre + (alpha.i - pp)/(pp*(1-pp))
  Sigma.z <- diag((fit$weights^(-1)),n) 
  Vq.beta <- solve(t(X)%*%solve(Sigma.z)%*%X + solve(Sigma.beta))
  qmu.beta <- Vq.beta%*%(t(X)%*%solve(Sigma.z)%*%zz+solve(Sigma.beta)%*%mu.beta)
  
  ### E(\sigma^2) w.r.t. q(sigma^2)
  sig2 <- Bqsig2/((A+n/2)-1)
  
  elbo.atual <- elbo.f(a0,b0,tau0.2,tau1.2,alpha.i,y,m,dim.beta,X,Vq.beta,s2,n,Sigma.beta,qmu.beta,mu.beta)
  
  while (elbo.ant < elbo.atual-0.00000001){ 
    elbo.ant<-elbo.atual
    
    ind <- ind+1
    sum.B <- 0
    
    dif <- sum_param.ant <- 10
    while (abs(dif) > 0.1) { 
      sum.B <- 0
    ## updating alpha
    for (i in 1:n){
      x.i <- as.vector(X[i,])
      D <- diag(1,dim.beta)
      D[,1] <- x.i/sqrt(sum(x.i^2))
      gs <- gramSchmidt(D)
      S <- gs$Q
      V.alpha <- t(S)%*%Vq.beta%*%(S)
      
      ## calculates E(log invlogit{-1}(Xbeta)) w.r.t. beta
      logq.beta1 <- aprox.GS1(x.i,qmu.beta,V.alpha,S)
      
      ## calculates E(log 1-invlogit{-1}(Xbeta)) w.r.t. beta
      logq.beta2 <- aprox.GS2(x.i,qmu.beta,V.alpha,S)
      
      alphai.t[i,1] <- exp(logq.beta2-0.5*qsig2.inv*((y[i]-m[1])^2+s2[1])/2) ## alpha_0
      alphai.t[i,2] <- exp(logq.beta1-0.5*qsig2.inv*((y[i]-m[2])^2+s2[2])/2) ## alpha_1
      
      alpha.i[i] <- alphai.t[i,2]/(alphai.t[i,1]+alphai.t[i,2])
      sum.B <- sum.B + (((y[i]-m[2])^2+s2[2])*alpha.i[i]) + (((y[i]-m[1])^2+s2[1])*(1-alpha.i[i]))
    }
    
    ## A and B come from prior for sigma^2
    Bqsig2 <- B + 0.5*sum.B
    qsig2.inv <-  (A+n/2)/Bqsig2 ### E[1/sigma^2]
    
    s2[1] <- 1/(qsig2.inv*sum((1-alpha.i))+1/tau0.2)
    m[1] <- (qsig2.inv*(t((1-alpha.i))%*%y)+mu.mu[1]/tau0.2)*s2[1]
    
    s2[2] <- 1/(qsig2.inv*sum((alpha.i))+1/tau0.2)
    m[2] <- (qsig2.inv*(t((alpha.i))%*%y)+mu.mu[2]/tau1.2)*s2[2]
    
    ## prior.df = 1 if cauchy; 
    ## prior.df = 7 if t-student; 
    ## prior.df = Inf if Normal
    ## link = logit, probit
    fit <- bayesglm(alpha.i ~ X-1,
                    family=binomial(link="logit"),
                    prior.scale=2.5, prior.df=Inf)
    
    regre <- X%*%fit$coefficients
    pp <- invlogit(regre)
    zz<- regre + (alpha.i - pp)/(pp*(1-pp))
    Sigma.z <- diag((fit$weights^(-1)),n) 
    Vq.beta <- solve(t(X)%*%solve(Sigma.z)%*%X + solve(Sigma.beta))
    qmu.beta <- Vq.beta%*%(t(X)%*%solve(Sigma.z)%*%zz+solve(Sigma.beta)%*%mu.beta)
    
    ### E(\sigma^2) w.r.t. q(sigma^2)
    sig2 <- Bqsig2/((A+n/2)-1)
    
    vec.par <- matrix(c(c(alpha.i),Bqsig2,c(m),c(s2),c(qmu.beta),c(Vq.beta)))
    sum_param.atual <- sqrt(t(vec.par)%*%vec.par)
    dif <- sum_param.ant-sum_param.atual
    print(dif)
    sum_param.ant <- sum_param.atual
    }
    
    elbo.atual <- elbo.f(a0,b0,tau0.2,tau1.2,alpha.i,y,m,dim.beta,X,Vq.beta,s2,n,Sigma.beta,qmu.beta,mu.beta)
    
    print(paste0(elbo.ant,elbo.atual))
  }
  return(list(sig2=sig2,alpha.i=alpha.i,s2=s2,mu=m,qmu.beta=qmu.beta,Vq.beta=Vq.beta))
} ## end funcion VB


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

## S matrix
knots1 <- c(16,16,16,16,18,20,22,25,30,35,40,45,50,55,60,60,60,60) ##seq(17,59)
x1 <- seq(16,60)
bs.matrix <- spline.des(knots1, x1, ord = 4, derivs=0, outer.ok = TRUE,
                        sparse = FALSE)
s.spline1 <- bs.matrix$design
n.bases <- dim(s.spline1)[2]
dim.beta <- n.bases*n.eeg+dim(Z)[2]

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

X <- cbind(Z,R.i)
set.seed(123456)
post.means <- VB(y,X,n,dim.beta) 
