## PS 236 Section, September 2008 ###
### OLS Monte Carlos               ###

## Starting with an aside, a function looks for a variable first in its
##   environment, then in subsequently higher environments (e.g., nested
##   functions), then globally
a <- 15

b <- function(x) {
  a*x
  }
  
b2 <- function(x) {
  a <- 3
  a*x
  }
  
b(3)
b2(3)

## Notice that R won't find variables defined in a function.
b3 <- function(x) {
  a2 <- 3
  a2*x
  }

b3(3)
a2

## We can use a double arrow to define a variable in a function globally
b3 <- function(x) {
  a2 <<- 3
  a2*x
  }

b3(3)
a2

## This can be useful for debugging your code, just like inserting
##   'print("Hello")' at different locations to identify problem areas

## Now to the main event!
## We will generate data, then run 'Monte Carlo' simulations to see
##   how accurate our parameter estimates are
## Monte Carlos are just simulations of random events

## First, load the MASS library to access the 'mvrnorm' function
library(MASS)

## Now, set the 'seed' of the random generator in R
set.seed(1027)

runif(5)

## Whenever you set the seed to the same number, you will get the
##  exact same results from random number draws
## To see this, reset the seed and take 'new' draws
set.seed(1027)

runif(5)

## Let's do a simulation of OLS
## In the end, we want to know how many runs in the simulation include
##    the correct treatment effect in the confidence interval
## We'll have 'N" observations with two covariates
## There will be 'runs' simulations
## 'means' will be the means of the covariates
## 'sigma' will be the covariance matrix between the covariates
##   that we specify
## NB: 'sigma' must be symmetric
## 'b' will be the covariate effects
## 'status' will be treatment and it will be assigned randomly
## 'Tr' will be the treatment effect
## 'eVar' will be the variance of epsilon (it has mean 0)
## 'alpha' will be the level of the confidence interval

## Now, let's write these into a function using 'sapply'
sims <- function(means=c(1,5), sigma=matrix(c(1,0.5,0.5,1), nrow=2),
  b=c(2,1), Tr=10, eVar=1, N=150, runs=1000, alpha=0.05){
  level <- 1 - alpha/2
  cVal <- qnorm(level)
  # Notice the '<<-' in the next line. It defines the 'z' variable globally
  z <<- sapply(1:runs, function(x){
    # First, generate the true values
    X <- mvrnorm(n=N, mu=means, Sigma=sigma)
    e <- rnorm(N, 0,eVar)
    status <- round(runif(N),0)
    status[1] <- ifelse(sum(status) == 0, 1, status[1])
    Y <- X %*% b + status*Tr + e
    # Now, estimate the parameters
    # First without covariates, but with treatment status
    nCV <- lm(Y ~ status)
    ci1 <- abs(nCV$coefficients[2] - Tr)/summary(nCV)$coefficients[2,2]
    inCI1 <- ifelse(ci1 < cVal, 1, 0)
    # Now with covariates
    CV <- lm(Y ~ status + X)
    ci2 <- abs(CV$coefficients[2] - Tr)/summary(CV)$coefficients[2,2]
    inCI2 <- ifelse(ci2 < cVal, 1, 0)
    Z <- c(nCV$coefficients[2], inCI1, CV$coefficients[2], inCI2)
    names(Z) <- c("No covar estimate", "No covar correct coverage",
    "With covars estimate", "With covars correct coverage")
    return(Z)
    } )
  rowMeans(z)
}

set.seed(1027)
sims()
Z

## Let's make 'status' correlated with the error
## Specifically, if the error is positive, then 'status' will be 1
endSims <- function(means=c(1,5), sigma=matrix(c(1,0.5,0.5,1), nrow=2),
  b=c(2,1), Tr=10, eVar=1, N=150, runs=1000, alpha=0.05){
  level <- 1 - alpha/2
  cVal <- qnorm(level)
  z <- sapply(1:runs, function(x){
    # First, generate the true values
    X <- mvrnorm(n=N, mu=means, Sigma=sigma)
    e <- rnorm(N, 0,eVar)
    # THE FOLLOWING LINE IS THE ONLY DIFFERENCE FROM 'sims'
    status <- ifelse(e > 0, 1, 0)
    status[1] <- ifelse(sum(status) == 0, 1, status[1])
    Y <- X %*% b + status*Tr + e
    # Now, estimate the parameters
    # First without covariates, but with treatment status
    nCV <- lm(Y ~ status)
    ci1 <- abs(nCV$coefficients[2] - Tr)/summary(nCV)$coefficients[2,2]
    inCI1 <- ifelse(ci1 < cVal, 1, 0)
    # Now with covariates
    CV <- lm(Y ~ status + X)
    ci2 <- abs(CV$coefficients[2] - Tr)/summary(CV)$coefficients[2,2]
    inCI2 <- ifelse(ci2 < cVal, 1, 0)
    Z <- c(cor(e, status), nCV$coefficients[2], inCI1,
      CV$coefficients[2], inCI2)
    names(Z) <- c("Treatment-error correlation",
      "No covar estimate", "No covar correct coverage",
      "With covars estimate", "With covars correct coverage")
    return(Z)
    } )
  rowMeans(z)
}

set.seed(1027)
endSims()
## Comprehension check: Does matching fix this problem?


## Let's make treatment correlated with the covariates, but not the
##   error
## Comprehension check: What results do you expect? ##
corSims <- function(means=c(1,5),
  sigma=matrix(c(1,0.5,0.25,0.5,1,0.5, 0.25, 0.5,1), nrow=3),
  b=c(2,1), Tr=10, eVar=1, N=150, runs=1000, alpha=0.05){
  level <- 1 - alpha/2
  cVal <- qnorm(level)
  z <- sapply(1:runs, function(x){
    # First, generate the true values
    # THE FOLLOWING LINE IS DIFFERENT FROM 'sims'
    X <- mvrnorm(n=N, mu=c(0, means), Sigma=sigma)
    e <- rnorm(N, 0,eVar)
    # THE FOLLOWING THREE LINES ARE DIFFERENT FROM 'sims'
    status <- ifelse(X[,1] > 0, 1, 0)
    status[1] <- ifelse(sum(status) == 0, 1, status[1])
    X <- X[,c(2:3)]
    Y <- X %*% b + status*Tr + e
    # Now, estimate the parameters
    # First without covariates, but with treatment status
    nCV <- lm(Y ~ status)
    ci1 <- abs(nCV$coefficients[2] - Tr)/summary(nCV)$coefficients[2,2]
    inCI1 <- ifelse(ci1 < cVal, 1, 0)
    # Now with covariates
    CV <- lm(Y ~ status + X)
    ci2 <- abs(CV$coefficients[2] - Tr)/summary(CV)$coefficients[2,2]
    inCI2 <- ifelse(ci2 < cVal, 1, 0)
    Z <- c(nCV$coefficients[2], inCI1, CV$coefficients[2], inCI2)
    names(Z) <- c("No covar estimate", "No covar correct coverage",
    "With covars estimate", "With covars correct coverage")
    return(Z)
    } )
  rowMeans(z)
}

set.seed(1027)
corSims()