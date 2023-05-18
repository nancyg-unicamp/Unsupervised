rm(list = ls())
rm(list=ls(all=TRUE))
rm(list=ls(all.names=TRUE))

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
  regre <- sqrt(t(x.i)%*%x.i)%*%alpha
  aprox <- log(invlogit(regre)+0.00000000000000000000001)*(dnorm(alpha,(t(S)%*%qmu.beta)[1],sqrt(V.alpha[1,1]))+10e-10) #,0)
  return(aprox)
}

aprox.GS2 <- function(x.i,qmu.beta,V.alpha,S){     
  z <- integrate(log.GS2,-Inf,Inf,x.i,qmu.beta,V.alpha,S,stop.on.error = F)
  return(z$value)
}

log.GS2 <- function(alpha,x.i,qmu.beta,V.alpha,S){
  regre <- sqrt(t(x.i)%*%x.i)%*%alpha
  aprox <- log(1-invlogit(regre)+0.00000000000000000000001)*(dnorm(alpha,(t(S)%*%qmu.beta)[1],sqrt(V.alpha[1,1]))+10e-10) #,0)
  return(aprox)
}

elbo.f <- function(a0,b0,tau0.2,tau1.2,alpha.i,y,m,dim.beta,X,Vq.beta,s2,n,Sigma.beta,qmu.beta,mu.beta){
  ###ELBO
  N<-n
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
  
  ## calculating E.2 to E.5 and Fs  
   
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
  if (is.nan(elbo)){elbo<--(10^7)}
  
  return(elbo)
}

VB <- function(y,X,n,dim.beta){
  N<-n
  #### hyperparameters
  ## prior hyperparameter
  tau0.2 <- 100 # variance for m_0
  tau1.2 <- 100 # variance for m_1
  a0 <- A <- 0.01
  b0 <- B <- 0.01
  shape.sigma <- (a0+n/2)
  mu.beta <- rep(0,dim(X)[2])
  Sigma.beta <- diag(2.5^2,dim(X)[2])
  
  ## initial value of mixture means
  mu.mu <- c(0,3)
  
  ### initial values for mixture variances (are equal)
  s2.0 <- var(y)
  
  ### for mean (mu.mu,sig2.mu,mu.beta,Sigma.beta)
  set.seed(123456)
  means = normalmixEM(y,2)
  sort.means <- sort(means$mu,index.return=T)
  
  v <- means$posterior[,2]
  m <- means$mu
  s2 <- (means$sigma)^2
  
  if (means$mu[1]>means$mu[2]){
    m <- means$mu[sort.means$ix]
    s2 <- (means$sigma[sort.means$ix])^2
    v <- 1-means$posterior[,2]
  }   
  
  ### beta - initial values
  qmu.beta <- rep(0,dim.beta)
  Vq.beta <- diag(0.1,dim.beta)
  
  ## initial values alpha
  kmeans <- kmeans(y,2)$cluster
  
  normalEM <- normalmixEM(y,2)
  class <- apply(normalEM$posterior,1,which.max)
  alpha.i <- pmax(normalEM$posterior[,2],0.01)
  
  #### calcula elbo.0
  
  elbo.ant <- elbo.f(a0,b0,tau0.2,tau1.2,alpha.i,y,m,dim.beta,X,Vq.beta,s2,n,Sigma.beta,qmu.beta,mu.beta)
  
  if(is.nan(elbo.ant)){elbo.ant<--(10^7)}
  
  qsig2.inv <- 1/var(y) ### 
  
  ind <- sig2 <- 0
  alphai.t <- matrix(0,ncol=2,nrow=N) ### posterior probability of c_i
  
  sum.B <- 0
  ## updating alpha
  for (i in 1:N){
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
  qsig2.inv <-  (A+N/2)/Bqsig2 ### E[1/sigma^2]
  
  s2[1] <- 1/(qsig2.inv*sum((1-alpha.i))+1/tau0.2)
  m[1] <- (qsig2.inv*(t((1-alpha.i))%*%y)+mu.mu[1]/tau0.2)*s2[1]
  
  s2[2] <- 1/(qsig2.inv*sum((alpha.i))+1/tau0.2)
  m[2] <- (qsig2.inv*(t((alpha.i))%*%y)+mu.mu[2]/tau1.2)*s2[2]
  
  
  ## prior.df = 1 if cauchy; prior.df = 7 if t-student; prior.df = Inf if Normal
  fit <- bayesglm(alpha.i ~ X-1,family=binomial(link="logit"),prior.scale=2.5, prior.df=Inf)
  
  regre <- X%*%fit$coefficients
  pp <- invlogit(regre)
  zz<- regre + (alpha.i - pp)/(pp*(1-pp))
  Sigma.z <- diag((fit$weights^(-1)),n)
  Vq.beta <- solve(t(X)%*%solve(Sigma.z)%*%X + solve(Sigma.beta))
  qmu.beta <- Vq.beta%*%(t(X)%*%solve(Sigma.z)%*%zz+solve(Sigma.beta)%*%mu.beta)
  
  ### E(\sigma^2) w.r.t. q(sigma^2)
  sig2 <- Bqsig2/((A+N/2)-1)
  
  elbo.atual <- elbo.f(a0,b0,tau0.2,tau1.2,alpha.i,y,m,dim.beta,X,Vq.beta,s2,n,Sigma.beta,qmu.beta,mu.beta)
  
  while ( elbo.ant < elbo.atual-0.00001){ 
    elbo.ant<-elbo.atual
    
    ind <- ind+1
    sum.B <- 0
    
    dif <- sum_param.ant <- 10
    while (abs(dif) > 1) { 
      sum.B <- 0
      ## updating alpha
      for (i in 1:N){
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
      qsig2.inv <-  (A+N/2)/Bqsig2 ### E[1/sigma^2]
      
      s2[1] <- 1/(qsig2.inv*sum((1-alpha.i))+1/tau0.2)
      m[1] <- (qsig2.inv*(t((1-alpha.i))%*%y)+mu.mu[1]/tau0.2)*s2[1]
      
      s2[2] <- 1/(qsig2.inv*sum((alpha.i))+1/tau0.2)
      m[2] <- (qsig2.inv*(t((alpha.i))%*%y)+mu.mu[2]/tau1.2)*s2[2]
      
      
      ## prior.df = 1 - cauchy; 
      ## prior.df = Inf - Normal;
      ## otherwise t-student; 
      fit <- bayesglm(alpha.i ~ X-1,family=binomial(link="logit"),prior.scale=2.5, prior.df=Inf)
      
      regre <- X%*%fit$coefficients
      pp <- invlogit(regre)
      zz<- regre + (alpha.i - pp)/(pp*(1-pp))
      Sigma.z <- diag((fit$weights^(-1)),n) 
      Vq.beta <- solve(t(X)%*%solve(Sigma.z)%*%X + solve(Sigma.beta))
      qmu.beta <- Vq.beta%*%(t(X)%*%solve(Sigma.z)%*%zz+solve(Sigma.beta)%*%mu.beta)
      
      ### E(\sigma^2) w.r.t. q(sigma^2)
      sig2 <- Bqsig2/((A+N/2)-1)
      
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

n <- 150
a0 <- b0 <- 0.01
tau0.2 <- 100
M <- 100

## T,S matrices
nk1 <- 7
knots1 <- seq(0,10,length=nk1)
t <- seq(0,10,length.out=256)
bs.x1 <- bs(t, knots = knots1)
s.spline1 <- matrix(as.numeric(bs.x1), nr = nrow(bs.x1)) 
n.bases1 <- dim(s.spline1)[2]

nk2 <- 5
knots2 <- seq(0,1,length=nk2)
s <- seq(0,1,length.out=256) 
bs.x2 <- bs(s, knots = knots2)
s.spline2 <- matrix(as.numeric(bs.x2), nr = nrow(bs.x2)) 
n.bases2 <- dim(s.spline2)[2]

dim.beta <- n.bases1+n.bases2

set.seed(12345)
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
  post.means <- VB(y,X,n,dim.beta) 
  end_time <- Sys.time()
}
