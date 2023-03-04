library(Matching)
library(rbounds) 
load(file=url("http://sekhon.berkeley.edu/causalinf/data/section10.RData"))

# using Luke Keele rbounds package    
?psens
psens(results.edu, GammaInc=.1, Gamma=2)

delta <- seq(1.39,3,.01)
lambda <- ((1.38 * seq(1.39,3,.01)) - 1)/(seq(1.39,3,.01) -1.38)

plot(lambda,delta, ylim=c(1.2,max(delta)), xlim=c(1,10), type="l")
abline(v=1.38, lty=2)
abline(h=1.38, lty=2)

cbind(delta,lambda)                          