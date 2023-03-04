###############################################
# Section 3: Bayesian Inference 
#
# John Henderson         
#
# PS 236B Spring 2013
###############################################

rm(list=ls(all=TRUE))
              
# 1. Linear model:
#   - known sigma and normal prior on beta

mu0 <- 1.2; tau0 <- 1.3; sigma <- 1.5; N <- 10
set.seed(1004)
y <- rnorm(mean=2.3,sd=sqrt(sigma),n=N)

tau_star <- function(n,tau0,sigma){
     return(1/((1/tau0^2)+(n/sigma^2)))}

mu_star <- function(y,n,mu0,tau0,sigma,tau_star){
     return(tau_star*((mu0/tau0^2)+
       (n*mean(y)/sigma^2)))}

ts <- tau_star(n=N,tau0,sigma)
ms <- mu_star(y,n=N,mu0,tau0,sigma,tau_star=ts)
posterior <- rnorm(mean=ms,sd=sqrt(ts),n=1000)
bhat <- mean(posterior)
   
par(mfrow=c(1,2))                          
plot(posterior,ylim=c(0,4))
abline(h=bhat,lwd=2,col='red')
abline(h=2.3,lty=1,col='pink',lwd=2)
abline(h=c(quantile(posterior,prob=c(.025,.975))),lty=2,lwd=2)   
    
# 2. Similar model with diffuse prior and unknown sigma

#sigma <- 1.5; N <- 10
#set.seed(1004)
#y <- rnorm(mean=2.3,sd=sqrt(sigma),n=N)
s2 <- 1/(N-1)*sum((y-mean(y))^2) 

invchi=((N-1)*s2)/(rchisq(1000,df=N-1))
posterior1 <- rnorm(mean=mean(y),sd=sqrt(invchi/N),n=1000)
bhat1 <- mean(posterior1)
                        
plot(posterior1,ylim=c(0,4))
abline(h=bhat1,lwd=2,col='blue')
abline(h=2.3,lty=1,col='pink',lwd=2)
abline(h=c(quantile(posterior1,prob=c(.025,.975))),lty=2,lwd=2)
     
# 3. Illustrate for normal prior, the posterior as a 'compromise' between prior and data

outs=array(0,150) 
mu_data <- mean(y)
sigma_data <- sd(y)
mu0 <- -1.2; tau0 <- 1.3; sigma <- 1.5;        

tau_star <- function(n,tau0,sigma){
     return(1/((1/tau0^2)+(n/sigma^2)))}

mu_star <- function(y,n,mu0,tau0,sigma,tau_star){
     return(tau_star*((mu0/tau0^2)+
       (n*mean(y)/sigma^2)))}
                 
for(k in 1:length(outs)){
	ts <- tau_star(n=k-1,tau0,sigma=sigma_data)
	ms <- mu_star(y=mu_data,n=k-1,mu0,tau0,sigma=sigma_data,tau_star=ts)
    outs[k] <- ms                 
}                         

plot(outs,ylim=c(mu0-.5,2.3+.5),col='red',lwd=.5)
abline(h=c(mu0,2.3),lty=2)
         
#  same plot with very different variance ration


outs=array(0,150) 
mu_data <- mean(y)
sigma_data <- 10.5 #sd(y)
mu0 <- -1.2; tau0 <- 1.3; sigma <- 1.5;        

tau_star <- function(n,tau0,sigma){
     return(1/((1/tau0^2)+(n/sigma^2)))}

mu_star <- function(y,n,mu0,tau0,sigma,tau_star){
     return(tau_star*((mu0/tau0^2)+
       (n*mean(y)/sigma^2)))}
                 
for(k in 1:length(outs)){
	ts <- tau_star(n=k-1,tau0,sigma=sigma_data)
	ms <- mu_star(y=mu_data,n=k-1,mu0,tau0,sigma=sigma_data,tau_star=ts)
    outs[k] <- ms                 
}
       
plot(outs,ylim=c(mu0-.5,2.3+.5),col='red',lwd=.5)
abline(h=c(mu0,2.3),lty=2)
               
# if there is more variability in y, then it takes many more observations to 
#  converge on the population average 
                                      
# END