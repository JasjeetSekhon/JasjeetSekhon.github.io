###############################################
# Section 4: Hierarchical Models and MCMC Estimation
#
# John Henderson         
#
# PS 236B Spring 2013
###############################################

rm(list=ls(all=TRUE))
 

# 1. Gibbs sampler; Bivariate Normal, with p=.65 correlation 
            
# condition each draw on the prior draw at t-1
set.seed(1005)
theta.gibbs.1 <- function(y,theta,p){
   rnorm(y[1]+p*(theta[2]-y[2]),1-p^2,n=1)}

theta.gibbs.2 <- function(y,theta,p){
   rnorm(y[2]+p*(theta[1]-y[1]),1-p^2,n=1)}

# use four starting chains
y1=c(2.5,2.5); y2=c(-2.5,-2.5)
y3=c(2.5,-2.5); y4=c(-2.5,2.5)

chains <- matrix(NA, 1000, 8)
chains[1,] <- c(y1, y2, y3, y4)
p <- .65;    

# center around origin; could do this elswhere
y <- c(0,0)
          
# gibbs sampler here
for(i in 2:1000){
 chains[i,c(1,3,5,7)] <- 
 c(theta.gibbs.1(y,theta=chains[i-1,1:2],p),
   theta.gibbs.1(y,theta=chains[i-1,3:4],p),                   
   theta.gibbs.1(y,theta=chains[i-1,5:6],p),
   theta.gibbs.1(y,theta=chains[i-1,7:8],p))

 chains[i,c(2,4,6,8)] <-
 c(theta.gibbs.2(y,theta=chains[i,1:2],p),
   theta.gibbs.2(y,theta=chains[i,3:4],p),
   theta.gibbs.2(y,theta=chains[i,5:6],p),   
   theta.gibbs.2(y,theta=chains[i,7:8],p))
}      

# need to discard first Bth samples, which may 
#  not have yet converged; discard 500 of 1000 samples
burnin <- 500
   
# here is each markov chain for the first 100 samples
plot(chains[,1:2],col='white',ylim=c(-3,3),xlim=c(-3,3))
for(i in 2:100){lines(lty=2,col='red',
   x=c(chains[i-1,1],chains[i,1]),
   y=c(chains[i-1,2],chains[i,2]))}

for(i in 2:100){lines(lty=2,col='blue',
   x=c(chains[i-1,3],chains[i,3]),
   y=c(chains[i-1,4],chains[i,4]))}

for(i in 2:100){lines(lty=2,col='darkgreen',
   x=c(chains[i-1,5],chains[i,5]),
   y=c(chains[i-1,6],chains[i,6]))}

for(i in 2:100){lines(lty=2,col='darkorange',
   x=c(chains[i-1,7],chains[i,7]),
   y=c(chains[i-1,8],chains[i,8]))}


# here is each markov chain after deleting burnin 
plot(chains[,1:2],col='white',ylim=c(-3,3),xlim=c(-3,3))
for(i in burnin:nrow(chains)){lines(lty=2,col='red',
   x=c(chains[i-1,1],chains[i,1]),
   y=c(chains[i-1,2],chains[i,2]))}

for(i in burnin:nrow(chains)){lines(lty=2,col='blue',
   x=c(chains[i-1,3],chains[i,3]),
   y=c(chains[i-1,4],chains[i,4]))}

for(i in burnin:nrow(chains)){lines(lty=2,col='darkgreen',
   x=c(chains[i-1,5],chains[i,5]),
   y=c(chains[i-1,6],chains[i,6]))}

for(i in burnin:nrow(chains)){lines(lty=2,col='darkorange',
   x=c(chains[i-1,7],chains[i,7]),
   y=c(chains[i-1,8],chains[i,8]))}

# can jitter sample repeatedly from these chains OR continue the markov chain if needed


# 2. Metropolis Algorith to estimate a bivariate linear model
                                                                        
rm(list=ls())
set.seed(1005); x <- rnorm(n=40); 
y <- -2.7 + 3.5*x + rnorm(n=40,sd=4); 
# s2 <- var(y)
           
# note everything is on log scale; so sum of log of independent distributions 
loglike <- function(param){
  a = param[1]; b = param[2]; sds= param[3];
  sll <- dnorm(y, a + b*x, sd = sds, log = T);
  sum_sll <- sum(sll); return(sum_sll)}

prior <- function(param){
  a <- param[1]; b = param[2]; sds = param[3];
  apr <- dnorm(a, sd=6, log = T)
  bpr <- dnorm(b, sd=6, log = T)  
  sdpr <- dunif(sds, min=0, max=20, log = T)
  return(apr+bpr+sdpr)}

poster <- function(param){
  return(loglike (param) + prior(param))}
                                    
# note some important tuning parameters;
#  - the jumping distribution is normal; 
#  - with radius 1 in the three 3-dimensional parameter space

propose <- function(param){
  return(rnorm(3,mean=param,sd=c(1,1,1)))}
       
# here is the MCMC runs 
runMCMC <- function(starts, iters){
  chain = matrix(NA,iters+1,3)
  probs = matrix(NA,iters+1,2)
  chain[1,] = starts
  for (i in 1:iters){
   prop = propose(chain[i,])
   probab = exp(poster(prop) - poster(chain[i,]))
   r=min(probab,1)
   ifs=sample(c(1,0),replace=F,
       prob=c(r,1-r),size=1)
   if(ifs==1){
     chain[i+1,] = prop
   } else if(ifs==0){
     chain[i+1,] = chain[i,]}}
  return(chain)}


starts= c(4,0,2)
iters=20000
chain = runMCMC(starts, iters)

burnin = 10000                                       

# rate of accepting proposals should be neither too high nor too low
#  - roughly between 20 to 30%
accept = 1-mean(duplicated(chain[-(1:burnin),]))
       

# check the performance of the model 
# red are true parameters
# blue are OLS paramers 

trueA=-2.7
trueB=3.5
trueSd=4

par(mfrow = c(2,3))
hist(chain[-(1:burnin),1],nclass=30, , main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnin),1]))
abline(v = trueA, col="red" )
abline(v=summary(lm(y~x))$coef[1],col='blue',lty=2,lwd=2)
hist(chain[-(1:burnin),2],nclass=30, main="Posterior of b", xlab="True value = red line")
abline(v = mean(chain[-(1:burnin),2]))
abline(v = trueB, col="red" )
abline(v=summary(lm(y~x))$coef[2],col='blue',lty=2,lwd=2)
hist(chain[-(1:burnin),3],nclass=30, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnin),3]) )
abline(v = trueSd, col="red" )
abline(v=sd(y),col='blue',lty=2,lwd=2)
plot(chain[-(1:burnin),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
abline(h = trueA, col="red" )
abline(h=summary(lm(y~x))$coef[1],col='blue',lty=2,lwd=2)
plot(chain[-(1:burnin),2], type = "l", xlab="True value = red line" , main = "Chain values of b", )
abline(h = trueB, col="red" )
abline(h=summary(lm(y~x))$coef[2],col='blue',lty=2,lwd=2)
plot(chain[-(1:burnin),3], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
abline(h = trueSd, col="red" )
abline(h=sd(y),col='blue',lty=2,lwd=2)

# for comparison:
summary(lm(y~x))
                                     

# 3. Model fit and evaluation issues 

# model performance is very important to check;
#  what to look out for?
#   - convergence in the trace plots, i.e., model convergence
#   - rejection/accept rates     
#   - sufficient burnin
#   - model fit 
#      - in and out of sample prediction; 
#      - test simulated data of a similar type for explicit fit
#
#  what to amend
#   - change prior distributions
#   - add further parameter structure (avoid unneeded complexities or dependencies)
#   - vary steps and step size in the proposal search; change symmetric jump distribution (or Metropolis-Hastings)
#   - sample many more iterations and burnin; many more chains 
#
#   BUT if the model is very hard to estimate using MCMC do not many dependencies; may need to pursue deterministic solutions

     
# low accept rate, model is not converged   
#  may need to enlarge how far the posterior is explored 
rm(list=ls())

set.seed(1005); x <- rnorm(n=500); 
y <- -2.7 + 3.5*x + rnorm(n=500,sd=4); s2 <- var(y)

loglike <- function(param){
  a = param[1]; b = param[2]; sds= param[3];
  sll <- dnorm(y, a + b*x, sd = sds, log = T);
  sum_sll <- sum(sll); return(sum_sll)}

prior <- function(param){

  a <- param[1]; b = param[2]; sds = param[3];
  apr <- dnorm(a, sd=6, log = T)
  bpr <- dnorm(b, sd=6, log = T)  
  sdpr <- dunif(sds, min=0, max=20, log = T)
  return(apr+bpr+sdpr)}

poster <- function(param){
  return(loglike (param) + prior(param))}

propose <- function(param){
  return(rnorm(3,mean=param,sd=c(.25,.25,.5)))}

runMCMC <- function(starts, iters){
  chain = matrix(NA,iters+1,3)
  probs = matrix(NA,iters+1,2)
  chain[1,] = starts
  for (i in 1:iters){
   prop = propose(chain[i,])
   probab = exp(poster(prop) - poster(chain[i,]))
   r=min(probab,1)
   ifs=sample(c(1,0),replace=F,
       prob=c(r,1-r),size=1)
   if(ifs==1){
     chain[i+1,] = prop
   } else if(ifs==0){
     chain[i+1,] = chain[i,]}}
  return(chain)}


starts= c(4,0,2)
iters=200
chain = runMCMC(starts, iters)

burnin = 100
accept = 1-mean(duplicated(chain[-(1:burnin),]))


trueA=-2.7
trueB=3.5
trueSd=4

par(mfrow = c(2,3))
hist(chain[-(1:burnin),1],nclass=30, , main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnin),1]))
abline(v = trueA, col="red" )
abline(v=summary(lm(y~x))$coef[1],col='blue',lty=2,lwd=2)
hist(chain[-(1:burnin),2],nclass=30, main="Posterior of b", xlab="True value = red line")
abline(v = mean(chain[-(1:burnin),2]))
abline(v = trueB, col="red" )
abline(v=summary(lm(y~x))$coef[2],col='blue',lty=2,lwd=2)
hist(chain[-(1:burnin),3],nclass=30, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnin),3]) )
abline(v = trueSd, col="red" )
abline(v=sd(y),col='blue',lty=2,lwd=2)
plot(chain[-(1:burnin),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
abline(h = trueA, col="red" )
abline(h=summary(lm(y~x))$coef[1],col='blue',lty=2,lwd=2)
plot(chain[-(1:burnin),2], type = "l", xlab="True value = red line" , main = "Chain values of b", )
abline(h = trueB, col="red" )
abline(h=summary(lm(y~x))$coef[2],col='blue',lty=2,lwd=2)
plot(chain[-(1:burnin),3], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
abline(h = trueSd, col="red" )
abline(h=sd(y),col='blue',lty=2,lwd=2)




# too high accept rate, model is not converged 
#  may need to shrink how far the posterior is explored           
rm(list=ls())

set.seed(1005); x <- rnorm(n=500); 
y <- -2.7 + 3.5*x + rnorm(n=500,sd=4); s2 <- var(y)

loglike <- function(param){
  a = param[1]; b = param[2]; sds= param[3];
  sll <- dnorm(y, a + b*x, sd = sds, log = T);
  sum_sll <- sum(sll); return(sum_sll)}

prior <- function(param){

  a <- param[1]; b = param[2]; sds = param[3];
  apr <- dnorm(a, sd=6, log = T)
  bpr <- dnorm(b, sd=6, log = T)  
  sdpr <- dunif(sds, min=0, max=20, log = T)
  return(apr+bpr+sdpr)}

poster <- function(param){
  return(loglike (param) + prior(param))}

propose <- function(param){
  return(rnorm(3,mean=param,sd=c(.1,.1,.1)))}

runMCMC <- function(starts, iters){
  chain = matrix(NA,iters+1,3)
  probs = matrix(NA,iters+1,2)
  chain[1,] = starts
  for (i in 1:iters){
   prop = propose(chain[i,])
   probab = exp(poster(prop) - poster(chain[i,]))
   r=min(probab,1)
   ifs=sample(c(1,0),replace=F,
       prob=c(r,1-r),size=1)
   if(ifs==1){
     chain[i+1,] = prop
   } else if(ifs==0){
     chain[i+1,] = chain[i,]}}
  return(chain)}

starts= c(4,0,2)
iters=500
chain = runMCMC(starts, iters)

burnin = 250
accept = 1-mean(duplicated(chain[-(1:burnin),]))

trueA=-2.7
trueB=3.5
trueSd=4

par(mfrow = c(2,3))
hist(chain[-(1:burnin),1],nclass=30, , main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnin),1]))
abline(v = trueA, col="red" )
abline(v=summary(lm(y~x))$coef[1],col='blue',lty=2,lwd=2)
hist(chain[-(1:burnin),2],nclass=30, main="Posterior of b", xlab="True value = red line")
abline(v = mean(chain[-(1:burnin),2]))
abline(v = trueB, col="red" )
abline(v=summary(lm(y~x))$coef[2],col='blue',lty=2,lwd=2)
hist(chain[-(1:burnin),3],nclass=30, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnin),3]) )
abline(v = trueSd, col="red" )
abline(v=sd(y),col='blue',lty=2,lwd=2)
plot(chain[-(1:burnin),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
abline(h = trueA, col="red" )
abline(h=summary(lm(y~x))$coef[1],col='blue',lty=2,lwd=2)
plot(chain[-(1:burnin),2], type = "l", xlab="True value = red line" , main = "Chain values of b", )
abline(h = trueB, col="red" )
abline(h=summary(lm(y~x))$coef[2],col='blue',lty=2,lwd=2)
plot(chain[-(1:burnin),3], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
abline(h = trueSd, col="red" )
abline(h=sd(y),col='blue',lty=2,lwd=2)


# very high accept rate; but beat the problem with many more iterations 
# metropolis
rm(list=ls())

set.seed(1005); x <- rnorm(n=500); 
y <- -2.7 + 3.5*x + rnorm(n=500,sd=4); s2 <- var(y)

loglike <- function(param){
  a = param[1]; b = param[2]; sds= param[3];
  sll <- dnorm(y, a + b*x, sd = sds, log = T);
  sum_sll <- sum(sll); return(sum_sll)}

prior <- function(param){

  a <- param[1]; b = param[2]; sds = param[3];
  apr <- dnorm(a, sd=6, log = T)
  bpr <- dnorm(b, sd=6, log = T)  
  sdpr <- dunif(sds, min=0, max=20, log = T)
  return(apr+bpr+sdpr)}

poster <- function(param){
  return(loglike (param) + prior(param))}

propose <- function(param){
  return(rnorm(3,mean=param,sd=c(.1,.1,.1)))}

runMCMC <- function(starts, iters){
  chain = matrix(NA,iters+1,3)
  probs = matrix(NA,iters+1,2)
  chain[1,] = starts
  for (i in 1:iters){
   prop = propose(chain[i,])
   probab = exp(poster(prop) - poster(chain[i,]))
   r=min(probab,1)
   ifs=sample(c(1,0),replace=F,
       prob=c(r,1-r),size=1)
   if(ifs==1){
     chain[i+1,] = prop
   } else if(ifs==0){
     chain[i+1,] = chain[i,]}}
  return(chain)}


starts= c(4,0,2)
iters=50000
chain = runMCMC(starts, iters)

burnin = 25000
accept = 1-mean(duplicated(chain[-(1:burnin),]))



trueA=-2.7
trueB=3.5
trueSd=4

par(mfrow = c(2,3))
hist(chain[-(1:burnin),1],nclass=30, , main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnin),1]))
abline(v = trueA, col="red" )
abline(v=summary(lm(y~x))$coef[1],col='blue',lty=2,lwd=2)
hist(chain[-(1:burnin),2],nclass=30, main="Posterior of b", xlab="True value = red line")
abline(v = mean(chain[-(1:burnin),2]))
abline(v = trueB, col="red" )
abline(v=summary(lm(y~x))$coef[2],col='blue',lty=2,lwd=2)
hist(chain[-(1:burnin),3],nclass=30, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnin),3]) )
abline(v = trueSd, col="red" )
abline(v=sd(y),col='blue',lty=2,lwd=2)
plot(chain[-(1:burnin),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
abline(h = trueA, col="red" )
abline(h=summary(lm(y~x))$coef[1],col='blue',lty=2,lwd=2)
plot(chain[-(1:burnin),2], type = "l", xlab="True value = red line" , main = "Chain values of b", )
abline(h = trueB, col="red" )
abline(h=summary(lm(y~x))$coef[2],col='blue',lty=2,lwd=2)
plot(chain[-(1:burnin),3], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
abline(h = trueSd, col="red" )
abline(h=sd(y),col='blue',lty=2,lwd=2)



# one benefit of alternative/pre-packaged approaches (jags/bugs, mcmcpack, etc) is that these 
#  search over tuning parameters without the user needing to be aware of these parameters

# the cost is that users are not made aware, and may suffer from estimating poor models without
#  doing the diagnostics 

                                 
# END