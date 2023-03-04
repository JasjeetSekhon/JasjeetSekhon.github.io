### PS 236 Section, October 2008   ###
### Randomization inference        ###


## We want to create a function that tests the hypothesis of no treatment effect
##   under treatment assignment according to fixed margins

## 'resp' is a vector of responses ('r' in the section slides)
## 'TR' is a vector of the actual treatment assignment ('Z' in the slides)
## 'nT' is the number of treated units ('m' in the section slides---only 1
##   strata)
## 'nC' is the number of control units
## 'runs' is the number of simulations

randInf <- function(resp, TR, runs=10000){
  nT <- sum(TR)
  nC <- length(TR) - nT
  Z <- sapply(1:runs, function(x){
    trPrime <- sample(c(rep(0,nC),rep(1,nT)),(nC+nT))
    sum(resp == trPrime)
    })
  correct <- sum(resp == TR)
  cat("p-value of null hypothesis of no effect:", sum(Z >= correct)/runs)
  }