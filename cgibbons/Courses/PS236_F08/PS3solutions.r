#### PS 236 Problem Set 3 ####
#### Fall 2008            ####
#### Solution code        ####

library(Matching)
setwd("C:/Users/Charles/Documents/R")


### Problem 2

## For all parts, define the inputs as:
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
## 'approach' takes four values:
##   1 is standard OLS
##   2 is matching
##   3 is matching followed by OLS
##   4 is matching with an OLS bias adjustment

sims <- function(means, sigma, approach,
  b=c(2,-2,4,1,-5), Tr=5, eVar=1, N=250, runs=1000, alpha=0.05){
  level <- 1 - alpha/2
  cVal <- qnorm(level)
  z <- sapply(1:runs, function(x){
    # First, generate the true values
    X <- mvrnorm(n=N, mu=means, Sigma=sigma)
    e <- rnorm(N, 0,eVar)
    status <- round(runif(N),0)
    status[1] <- ifelse(sum(status) == 0, 1, status[1])
    Y <- X %*% b + status*Tr + e
    # Now, estimate the parameters
    # With OLS
    if (approach == 1) {
      CV <- lm(Y ~ status + X)
      ci <- abs(CV$coefficients[2] - Tr)/summary(CV)$coefficients[2,2]
      inCI <- ifelse(ci < cVal, 1, 0)
    }
    if (approach == 2) {
      CV <- Match(Y=Y, Tr=status, X=X)
      ci <- abs(CV$est - Tr)/CV$se
      inCI <- ifelse(ci < cVal, 1, 0)
    }
    if (approach == 3) {
      mCV <- Match(Y=Y, Tr=status, X=X)
      Y <- mCV$mdata$Y
      status <- mCV$mdata$Tr
      X <- mCV$mdata$X
      CV <- lm(Y ~ status + X)
      ci <- abs(CV$coefficients[2] - Tr)/summary(CV)$coefficients[2,2]
      inCI <- ifelse(ci < cVal, 1, 0)
    }
    if (approach == 4) {
      CV <- Match(Y=Y, Tr=status, X=X, Z=X, BiasAdjust=TRUE)
      ci <- abs(CV$est - Tr)/CV$se
      inCI <- ifelse(ci < cVal, 1, 0)
    }
    Z <- c(inCI)
    names(Z) <- c("Correct coverage %")
    return(Z)
    } )
  mean(z)
}

## Part a
set.seed(1027)
means <- c(round(runif(5,-20,20)))
a <- matrix(runif(25,-1,1), ncol=5)
sigma <- a %*% t(a)
sims(means=means, sigma=sigma, approach=1)

## Part b
set.seed(1027)
means <- c(round(runif(5,-20,20)))
a <- matrix(runif(25,-1,1), ncol=5)
sigma <- a %*% t(a)
sims(means=means, sigma=sigma, approach=2)

## Part c
set.seed(1027)
means <- c(round(runif(5,-20,20)))
a <- matrix(runif(25,-1,1), ncol=5)
sigma <- a %*% t(a)
sims(means=means, sigma=sigma, approach=3)

## Part d
set.seed(1027)
means <- c(round(runif(5,-20,20)))
a <- matrix(runif(25,-1,1), ncol=5)
sigma <- a %*% t(a)
sims(means=means, sigma=sigma, approach=4)


### Problem 3

## Load the data
X <- read.csv("ExactMatching.csv")

## Part a
# Perform OLS
summary(lm(Outcome ~ Treatment + X1 + X2, data=X))

## Part b
# Divide the sample into treatment and control
treated <- X[X$Treatment == 1,]
control <- X[X$Treatment == 0,]

# Create estimates of the treatment effect for a treated individual and similar
#   controls
pairEst <-sapply(1:dim(treated)[1], function(x){
  # Record the covariate values for the treated observation
  x1 <- treated[x,3]
  x2 <- treated[x,4]
  # Find the mean of the outcome (column 1) for controls with the same covariate
  #   values
  contMean <- mean(control[control[,3] == x1 & control[,4] == x2,1])
  # Calculate the difference between the treated and mean of the controls
  treated[x,1] - contMean
  })

# Find the number of dropped treated observations
sum(is.na(pairEst))

# Calculate the estimated treatment effect
mean(pairEst, na.rm=TRUE)

# Check your answer
exMatched <- Match(Y=X$Outcome, Tr=X$Treatment, X=X[,c(3,4)], exact=TRUE)
summary(exMatched)

## Part c
matched <- Match(Y=X$Outcome, Tr=X$Treatment, X=X[,c(3,4)])
summary(matched)