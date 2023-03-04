### PS 236 Problem Set 1   ###
### Fall 2008              ###
### Solutions              ###


## Question 2              ##

## We want to create a function that tests the hypothesis of no treatment effect
##   under treatment assignment according to binomial randomization

## 'resp' is a vector of responses ('r' in the section slides)
## 'TR' is a vector of the actual treatment assignment ('Z' in the slides)
## 'runs' is the number of simulations

randInf <- function(resp, TR, runs=10000){
  Z <- sapply(1:runs, function(x){
    trPrime <- rbinom(length(TR), 1, 0.5)
    sum(resp == trPrime)
    })
  correct <- sum(resp == TR)
  cat("p-value of null hypothesis of no effect:", sum(Z >= correct)/runs,"\n")
  }

## Using Mike's guesses and the actual assignment
randInf(c(1,0,1,1,0,0,0,1,1,0), c(1,0,1,0,0,0,1,0,1,0))


## Question 3              ##

## 'p' is the probability of treatment assignment when treatment group is too
##   small
## 'runs' is the number of simulations
## 'wind' is the distance from the 0.5 proportion that the treated group is
##   permitted
## 'n' is the number of units in each simulation
## 'nT' is a vector of the results for each simulation

tAssign <- function(p=0.5, runs=1000, wind=0.001){
  z <- sapply(1:runs, function(x){
    n <- round(runif(1, 100, 1000))
    nT <- rbinom(1,1,0.5)
    sapply(2:n, function(y){
      fracT <- sum(nT)/length(nT)
      if (fracT < 0.5) {
        ti <- rbinom(1,1,p)
        nT <<- c(nT, ti)
        }
      if (fracT > 0.5) {
        ti <- rbinom(1,1,(1-p))
        nT <<- c(nT, ti)
        }
      if (fracT == 0.5) {
        ti <- rbinom(1,1,0.5)
        nT <<- c(nT, ti)
        }
      })
    sum(nT)/length(nT)
    })
  inRange <- sum(abs(z - 0.5) <= wind)/runs
  paste("Fraction of sims within",wind,"window:",inRange)
  }

## A more parsimonious version of the above code
tAssign <- function(p=0.5, runs=1000, wind=0.001){
  z <- sapply(1:runs, function(x){
    n <- round(runif(1, 100, 1000))
    nT <- rbinom(1,1,0.5)
    sapply(2:n, function(y){
      fracT <- sum(nT)/length(nT)
      ti <- rbinom(1,1,((fracT < 0.5)*p+(fracT > 0.5)*(1-p) + (fracT == 0.5)*0.5))
      nT <<- c(nT, ti)
      })
    sum(nT)/length(nT)
    })
  inRange <- sum(abs(z - 0.5) <= wind)/runs
  paste("Fraction of sims within",wind,"window:",inRange)
  }

## Perform the simulations for different values of p
Z <- sapply(seq(0.5,1,0.05), function(x) tAssign(p=x))
