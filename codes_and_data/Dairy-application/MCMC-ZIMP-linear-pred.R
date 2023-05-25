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

## Design matrix Z
Z <- cbind(rep(1,n),covar)

M <- 15000
ML <- 10000+1



## matrix S
knots1 <- c(rep(1,4),seq(5,25,length.out=4),rep(31,4))
x1 <- seq(1,n.freqs)
bs.matrix <- spline.des(knots1, x1, ord = 4, derivs=0, outer.ok = TRUE,
                        sparse = FALSE)
s.spline1 <- bs.matrix$design 
n.bases <- dim(s.spline1)[2]

eeg_tnsr <- array(0,dim=c(n,n.eeg,n.freqs))

eeg_tnsr<-eeg.vacas

## matrix R
R <- array(0,dim=c(n,n.eeg,n.bases))
for (i in 1:n){
  for (j in 1:n.eeg){
    R[i,j,] <- t(s.spline1)%*%eeg_tnsr[i,j,]
  }
}

##  R_{i,j,}%*%t(R_{i,j,})
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



## prior hyperparameters

a0.1 <- 0.01 ## gamma shape class 1 
b0.1 <- 0.01 ## gamma rate

a0.2 <- 0.05 ## gamma parameters class 2
b0.2 <- 0.01

## initial values gamma

gamma.post1 <- array(0,dim=c(M,3,n))
gamma.post1[1,,] <- rmultinom(n,1,c(1/3,1/3,1/3))

## initial values phi
phi1.post1 <- array(0,dim=c(M,n.eeg*n.bases)) 
phi2.post1 <- array(0,dim=c(M,n.eeg*n.bases))

## initial values theta
theta1.post1 <- matrix(0,nrow=M,ncol=dim(Z)[2])


theta2.post1 <- matrix(0,nrow=M,ncol=dim(Z)[2])


## initial values for lambdas
lambda1.post1 <- lambda2.post1 <-  rep(0,M)
ratio <- array(0,dim=c(M,3,n))

lambda1.post1[1] <- 0.01 
lambda2.post1[1] <- 0.05 


set.seed(12345)
start.time <- Sys.time()
for (r in  2:M){
  
  ## Update gamma

  regre1 <- Z%*%theta1.post1[r-1,]+R.i%*%as.vector(t(phi1.post1[r-1,]))
  regre2 <- Z%*%theta2.post1[r-1,]+R.i%*%as.vector(t(phi2.post1[r-1,]))
  
  
  ratio.1 <- exp(-lambda1.post1[r-1])*(lambda1.post1[r-1]^y)*exp(regre1)
  ratio.2 <- exp(-lambda2.post1[r-1])*(lambda2.post1[r-1]^y)*exp(regre2)
 
  ratio[r,1,] <- ifelse(y==0,1,0)/(ifelse(y==0,1,0) + ratio.1 + ratio.2)
  ratio[r,2,] <- ratio.1/(ifelse(y==0,1,0) + ratio.1 + ratio.2)
  ratio[r,3,] <- ratio.2/(ifelse(y==0,1,0) + ratio.1 + ratio.2)
  
  gamma.post1[r,,] <- sapply(1:n,function(x){rmultinom(1,size=1,prob=c(ratio[r,1,x],ratio[r,2,x],ratio[r,3,x]))})
  
  ## update lambda0 
  a1.1 <- a0.1 + sum(y*gamma.post1[r,2,])
  b1.1 <- b0.1 + sum(gamma.post1[r,2,])
  
  lambda1.post1[r] <- rgamma(1,a1.1,rate=b1.1)
  
  ## update lambda1 
  a1.2 <- a0.2 + sum(y*gamma.post1[r,3,])
  b1.2 <- b0.2 + sum(gamma.post1[r,3,])
  
  lambda2.post1[r] <- rgamma(1,a1.2,rate=b1.2)
  

  ### update thetas e phis
  ## prior.df = 1 is Cauchy
  ## prior.df = Inf is Normal (default)
  ## Otherwise prior is t-Student
  ## link: logit/probit
  
  thetaphi1 <- bayesglm (gamma.post1[r,2,] ~ Z+R.i-1, 
                   family=binomial(link="logit"),
                  prior.scale=2.5, prior.df=1)

  coefs1 <-  thetaphi1$coef
  theta1.post1[r,] <- coefs1[1:dim(Z)[2]]
  phi1.post1[r,] <- coefs1[(dim(Z)[2]+1):length(coefs1)]
  
  
  thetaphi2 <- bayesglm (gamma.post1[r,3,] ~ Z+R.i-1, 
                         family=binomial(link="logit"),
                         prior.scale=2.5, prior.df=1)
  
  coefs2 <-  thetaphi2$coef
  theta2.post1[r,] <- coefs2[1:dim(Z)[2]]
  phi2.post1[r,] <- coefs2[(dim(Z)[2]+1):length(coefs2)]
   
  

}
end.time <- Sys.time()
end.time-start.time

ss <- seq(ML,M,by=100)

 for (r in ss){
   write.table(t(theta1.post1[r,]),"theta1.zip.logit-linear-cauchy.dat",row.names=F,col.names=F,append=T)
   write.table(t(theta2.post1[r,]),"theta2.zip.logit-linear-cauchy.dat",row.names=F,col.names=F,append=T)
   write.table(t(lambda1.post1[r]),"lambda1.zip.logit-linear-cauchy.dat",row.names=F,col.names=F,append=T)
   write.table(t(lambda2.post1[r]),"lambda2.zip.logit-linear-cauchy.dat",row.names=F,col.names=F,append=T)
   write.table(t(gamma.post1[r,,]),"gamma.zip.logit-linear-cauchy.dat",row.names=F,col.names=F,append=T)
   write.table(t(phi1.post1[r,]),"phi1.zip.logit-linear-cauchy.dat",row.names=F,col.names=F,append=T)
   write.table(t(phi2.post1[r,]),"phi2.zip.logit-linear-cauchy.dat",row.names=F,col.names=F,append=T)
   write.table(t(ratio[r,,]),"ratio.zip.logit-linear-cauchy.dat",row.names=F,col.names=F,append=T)
 }

write.table(t(y),"y.dat",row.names=F,col.names=F,append=T)
write.table(t(y.test),"y-test.dat",row.names=F,col.names=F,append=T)

gamma.est<-apply(gamma.post1[ss,,],c(2,3),mean)

class.vaca <- apply(gamma.est,2,max)
class.vaca <- ifelse(class.vaca==gamma.est[1,],1,ifelse(class.vaca==gamma.est[2,],2,3))

theta1.est <- apply(theta1.post1[ss,],2,mean)
theta2.est <- apply(theta2.post1[ss,],2,mean)

gamma.est<-apply(gamma.post1[ss,,],c(2,3),mean)


pdf("ZIP_lambda.pdf")
par(mfrow=c(1,2))
hist(lambda1.post1[ss],xlab=expression(lambda[1]),prob=T,main="")
hist(lambda2.post1[ss],xlab=expression(lambda[2]),prob=T,main="")
dev.off()

#### estimating w_j(t)

library(depth)

dd=n.bases
ss1 <- seq(1,dd)
ss2 <- seq(dd+1,2*dd)
ss3 <- seq(2*dd+1,3*dd)
ss4 <- seq(3*dd+1,4*dd)

phi1.1 <- phi1.post1[ss,ss1]
phi2.1 <- phi1.post1[ss,ss2]
phi3.1 <- phi1.post1[ss,ss3]
phi4.1 <- phi1.post1[ss,ss4]
w1.1 <- phi1.1%*%t(s.spline1)
w2.1 <- phi2.1%*%t(s.spline1)
w3.1 <- phi3.1%*%t(s.spline1)
w4.1 <- phi4.1%*%t(s.spline1)

phi1.2 <- phi2.post1[ss,ss1]
phi2.2 <- phi2.post1[ss,ss2]
phi3.2 <- phi2.post1[ss,ss3]
phi4.2 <- phi2.post1[ss,ss4]
w1.2 <- phi1.2%*%t(s.spline1)
w2.2 <- phi2.2%*%t(s.spline1)
w3.2 <- phi3.2%*%t(s.spline1)
w4.2 <- phi4.2%*%t(s.spline1)

pdf("ZIP_w1.pdf")
par(mfrow=c(2,2))
box.func1<-fbplot(t(w1.1),plot=T,xlab="Days",ylab=TeX("$w_{11}$"))
abline(h=0)
box.func2<-fbplot(t(w2.1),plot=T,xlab="Days",ylab=TeX("$w_{21}$"))
abline(h=0)
box.func2<-fbplot(t(w3.1),plot=T,xlab="Days",ylab=TeX("$w_{31}$"))
abline(h=0)
box.func2<-fbplot(t(w4.1),plot=T,xlab="Days",ylab=TeX("$w_{41}$"))
abline(h=0)
dev.off()



pdf("ZIP_w2.pdf")
par(mfrow=c(2,2))
box.func1<-fbplot(t(w1.2),plot=T,xlab="Days",ylab=TeX("$w_{12}$"))
abline(h=0)
box.func2<-fbplot(t(w2.2),plot=T,xlab="Days",ylab=TeX("$w_{22}$"))
abline(h=0)
box.func2<-fbplot(t(w3.2),plot=T,xlab="Days",ylab=TeX("$w_{32}$"))
abline(h=0)
box.func2<-fbplot(t(w4.2),plot=T,xlab="Days",ylab=TeX("$w_{42}$"))
abline(h=0)
dev.off()

par(mfrow=c(1,1))

##########  Predictive values for new observation y[new] with eeg_tnsr[new]


n.t <- length(y.test)

## montando matrix R
R.test <- array(0,dim=c(n.t,n.eeg,n.bases))
for (i in 1:n.t){
  for (j in 1:n.eeg){
    R.test[i,j,] <- t(s.spline1)%*%eeg.vacas.test[i,j,]
  }
}

R.i.test <- matrix(0,n.t,n.eeg*n.bases)
for (i in 1:n.t){
  for (j in 1:n.eeg){
    R.i.test[i,((j-1)*n.bases+1):((j-1)*n.bases+n.bases)] <- R.test[i,j,]    
  }
}


## Design matrix Z
Z.test <- cbind(rep(1,n.t),covar.test)

y0.test<-ifelse(y.test==0,1,0)


gamma.prev  <- matrix(0,nrow=n.t,ncol=3)

for (s in ss){
  regre1 <- Z.test%*%theta1.post1[s,]+R.i.test%*%as.vector(t(phi1.post1[s,]))
  regre2 <- Z.test%*%theta2.post1[s,]+R.i.test%*%as.vector(t(phi2.post1[s,]))
      p.0 <- 1/(1+exp(regre1)+exp(regre2))
      p.1 <- exp(regre1)/(1+exp(regre1)+exp(regre2))
      p.2 <- exp(regre2)/(1+exp(regre1)+exp(regre2))
      
  gamma.prev[,1] <- gamma.prev[,1]+p.0
  gamma.prev[,2] <- gamma.prev[,2]+p.1
  gamma.prev[,3] <- gamma.prev[,3]+p.2
    }

gamma.prev <- gamma.prev/length(ss)

dg<- rgb(0,100/255,0,alpha=0.5)
mb <- rgb(123/255,104/255,238/255,alpha=0.5)
r3<- rgb(220/255,20/255,60/255,alpha=0.5)


#################  Plots 

Class <- c(rep("Pure Zero",208),rep("Low",208),rep("High",208))
Health <- c(y,y,y)
Gamma.est<- c(gamma.est[1,],gamma.est[2,],gamma.est[3,])


data <- data.frame(Health,Class,Gamma.est)
colnames(data) <-  c("Health Events", "Class", "gamma")
ggplot(data, aes(fill=Class, y=gamma, x=Health)) + 
  geom_bar(position="fill", stat="identity",alpha=0.5) +
  scale_fill_manual(name="Class",values=c("mediumslateblue","darkgreen","Red3"),
                    labels=c("High", "Low","Pure Zero")) +
  ylab("Estimated probabilities") + xlab("Number of Health Events")
  theme_bw() + theme(panel.grid.minor = element_blank())

colnames(gamma_prev.frame) <- c("Health events", "gamma[1]","gamma[2]","gamma[3]")

tidy.g<-gather(gamma_prev.frame,key="Key",value="Values")

ggplot(tidy.g, aes(fill=))



lambda.frame<-as.data.frame(cbind(c(lambda1.post1[ss],lambda2.post1[ss]),
                                  c(rep(1,length(ss)),rep(2,length(ss)))))
colnames(lambda.frame) <- c("lambda", "Intensity")
lambda.frame$Intensity <- as.factor(lambda.frame$Intensity)

ggplot(lambda.frame,aes(x=lambda)) +
  geom_histogram(data=subset(lambda.frame,Intensity == 1),
                 alpha = 0.5,bins=20,aes(y=..density.., fill=Intensity)) +
  geom_histogram(data=subset(lambda.frame,Intensity == 2),
                 alpha = 0.5,bins=20,aes(y=..density..,fill = Intensity)) +
  xlab(expression(lambda)) +
  scale_fill_manual(name="Intensity",values=c("darkgreen", "mediumslateblue"),
                    labels=c("Low","High")) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())





theta.frame <- as.data.frame(cbind(rbind(theta1.post1[ss,],theta2.post1[ss,]),
                                    c(rep(1,length(ss)),rep(2,length(ss)))))
colnames(theta.frame) <- c("Intercept","theta[1]", "theta[2]", "theta[3]", "theta[4]", "Intensity")
theta.frame$Intensity <- as.factor(theta.frame$Intensity)


ggplot(theta.frame,aes(x=Intercept)) +
  geom_histogram(data=subset(theta.frame,Intensity == 1),
                 alpha = 0.5,bins=7,aes(y=..density.., fill=Intensity)) +
  geom_histogram(data=subset(theta.frame,Intensity == 2),
                 alpha = 0.5,bins=9,aes(y=..density..,fill = Intensity))+
  xlab("Intercept")  +
  scale_fill_manual(values=c("darkgreen", "mediumslateblue")) +
  guides(fill = "none") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())


temp <- theta.frame
colnames(temp)[2:5] <- c("theta[1]", "theta[2]", "theta[3]", "theta[4]")
tidy <- gather(temp, key="Type", value = "Value", 2:5)


ggplot(tidy,aes(x=Value)) +
  geom_histogram(data=subset(tidy,Intensity == 1),
                 alpha = 0.5,bins=7,aes(y=..density.., fill=Intensity)) +
  geom_histogram(data=subset(tidy,Intensity == 2),
                 alpha = 0.5,bins=9,aes(y=..density..,fill = Intensity))  +
  xlab("") +
  scale_fill_manual(name="Intensity",values=c("darkgreen", "mediumslateblue"),
                    labels=c("Low","High")) +
  theme_bw() +  theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) +
  facet_wrap(~Type, scales = "free", nrow=2,
             labeller = label_parsed)
