### PS 236 Section, November 2008     ###
### The Bootstrap                     ###

setwd("C:/Users/Charles/Documents/Grad school/Third Year/Political Science 236")
load("GM.RData")
output <- GM[,1]
materials <- GM[,4]
capital <- GM[,3]
labor <- GM[,2]

## You are hired by General Motors to study their production line. Specifically,
##   GM believes that production exhibits "increasing returns to scale,"
##   under which total average cost per vehicle would be lower at higher levels
##   of production.
## You assume that production involves three factors---materials, labor, and
##   capital---and follows a Cobb-Douglas production function:
##     Y = M^(alpha) * L^(beta) * K^(gamma).
## Increasing returns to scale exist if and only if
##     alpha + beta + gamma > 1.    (1)
## You run the regression
##     y = alpha*m + beta*l + gamma*k,
##   where the lowercase letters represent the log of their capital counterparts.
## You test hypothesis (1) by bootstrapping.

bootOLS <- function(Y, X, nSims=1000, Alpha=0.05){
  boots <- sapply(1:nSims, function(x){
    # We need to choose observatons from our dataset randomly for each simulation.
    #   By default, the number of observations chosen will be the number in the
    #   dataset (i.e., dim(X)[1]). We sample with replacement
    obs <- sample(1:dim(X)[1], replace=TRUE)
    # Now, we make a new X and Y based just on these rows
    Xb <- X[obs,]
    Yb <- Y[obs]
    # Do a regression using these new subsamples
    bootReg <- lm(Yb ~ Xb)
    # And keep the estimated coefficients
    bootB <- bootReg$coefficients
  })
  # Let's do a confidence interval for our GM scenario.
  # alpha, beta, and gamma correspond to the second through fourth coefficient,
  #   so let's sum them. Note that the coefficients for a given simulation
  #   appear as a column and each row is a different coefficient.
  sums <- colSums(boots[2:dim(boots)[1],])
  # The estimate is the mean of this result
  est <- mean(sums)
  # The standard error of the results
  se <- sd(sums)
  # The Alpha critical value is used to find the confidence interval
  ci <- quantile(sums, c(Alpha/2, (1-Alpha/2)))
  cat("    Estimate:", est, "\n    Standard Error:", se, "\n   ", Alpha*100, "% Confidence Interval: (", ci[1], ",", ci[2],")\n")
  list(est, se, ci)
}

bootGM <- bootOLS(Y=output, X=cbind(materials, labor, capital))
olsGM <- lm(output ~ cbind(materials, labor, capital))
summary(olsGM)

  
  
## Hypothesis tests
## Case 1
## This is a relatively straightforward application of the bootstrap

case1 <- function(x, z, nSims=1000, null=0){
  diffs <- sapply(1:nSims, function(y){
    # Sample from the samples, keeping each sample's size the same
    xb <- x[sample(1:length(x), replace=TRUE)]
    zb <- z[sample(1:length(z), replace=TRUE)]
    # Take the difference of the subsamples
    mean(xb) - mean(zb)
  })
  # Count how often results appear in the smaller tail (0 < p < 0.5) and double
  #   it. Think about why the following line makes sense (Hint: Draw a picture
  #   of what the line is counting).
  bT <- 2*ifelse(mean(x) - mean(z) > null, sum(diffs < null)/length(diffs),
    sum(diffs >= null)/length(diffs))
  # Compute the analytical p-value
  aT <- 2*(1 - pnorm(abs(mean(x) - mean(z) - null)/sqrt(var(x)/length(x) + var(z)/length(z))))
  cat("    Bootstrapped p-value:", bT, "\n    Analytical p-value:", aT, "\n")
  c(bT,aT)
}
x <- rnorm(50, mean=0, sd=1)
z <- rnorm(75, mean=5, sd=2)
case1(x=x,z=z,null=-5)
t.test(x,z,mu=-5)$p.val



## Case 2
## Here, we are essentially asking, "how often is the difference between the
##   bootstrapped statistic and the full sample statistic bigger than the difference
##   between the full sample statistic and the null hypothesis?"

case2 <- function(x, z, nSims=1000, null=0){
  # Find the full sample statistic
  fs <- mean(x) - mean(z)
  diffs <- sapply(1:nSims, function(y){
    xb <- x[sample(1:length(x), replace=TRUE)]
    zb <- z[sample(1:length(z), replace=TRUE)]
    # Take the difference of the subsamples
    mean(xb) - mean(zb)
  })
  # See how often the bootstrapped difference minus the full sample difference
  #   is bigger than the full sample difference minus the null hypothesis in
  #   the smaller tail. Again, think about what the next line means.
  bT <- 2*ifelse(mean(x) - mean(z) > null, sum(diffs - fs >= fs - null)/length(diffs),
    sum(diffs - fs < fs - null)/length(diffs))
  # Compute the analytical p-value
  aT <- 2*(1 - pnorm(abs(mean(x) - mean(z) - null)/sqrt(var(x)/length(x) + var(z)/length(z))))
  cat("    Bootstrapped p-value:", bT, "\n    Analytical p-value:", aT, "\n")
  c(bT,aT)
}
x <- rnorm(50, mean=0, sd=1)
z <- rnorm(75, mean=5, sd=2)
case2(x=x,z=z,null=-5)
t.test(x,z,mu=-5)$p.val



## Case 3: Studentized version
## This case is very much like case 2, except now we standardize the distances

case3 <- function(x, z, nSims=1000, null=0){
  # Studentize the null
  Snull <- null/sqrt(var(x)/length(x) + var(z)/length(z))
  # Find the full sample statistic
  fs <- (mean(x) - mean(z) - null)/sqrt(var(x)/length(x) + var(z)/length(z))
  diffs <- sapply(1:nSims, function(y){
    xb <- x[sample(1:length(x), replace=TRUE)]
    zb <- z[sample(1:length(z), replace=TRUE)]
    # Take the studentized difference of the subsamples
    (mean(xb) - mean(zb) - null)/sqrt(var(xb)/length(xb) + var(zb)/length(zb))
  })
  # Note on the following line that the null is incorporated into the test
  #   statistics already and thus does not appear separately
  bT <- 2*ifelse(fs > 0, sum(diffs - fs >= fs)/length(diffs),
    sum(diffs - fs < fs)/length(diffs))
  aT <- 2*(1 - pnorm(abs(mean(x) - mean(z) - null)/sqrt(var(x)/length(x) + var(z)/length(z))))
  cat("    Bootstrapped p-value:", bT, "\n    Analytical p-value:", aT, "\n")
  c(bT,aT)
}
x <- rnorm(50, mean=0, sd=1)
z <- rnorm(75, mean=2, sd=2)
case3(x=x,z=z,null=-2)
t.test(x,z,mu=-2)$p.val

## Case 4: Same distributions
## This is just like Case 2, except we make our draws differently
## Null of no difference in means assumed
case4 <- function(x, z, nSims=1000){
  fs <- mean(x) - mean(z)
  # Assume (as under the null) that all the observations came from the same
  #   process and join them together
  allObs <- c(x,z)
  diffs <- sapply(1:nSims, function(y){
    # Take draws from the united set of observations, keeping the x and z
    #   vectors the same length.
    xb <- allObs[sample(1:length(allObs), size=length(x), replace=TRUE)]
    zb <- allObs[sample(1:length(allObs), size=length(z), replace=TRUE)]
    mean(xb) - mean(zb)
  })
  bT <- 2*ifelse(mean(x) - mean(z) > 0, sum(diffs >= fs)/length(diffs),
    sum(diffs < fs)/length(diffs))
  aT <- 2*(1 - pnorm(abs(mean(x) - mean(z))/sqrt(var(x)/length(x) + var(z)/length(z))))
  cat("    Bootstrapped p-value:", bT, "\n    Analytical p-value:", aT, "\n")
  c(bT,aT)
}
x <- rnorm(50, mean=0, sd=1)
z <- rnorm(75, mean=0, sd=2)
case4(x=x,z=z)
t.test(x,z,mu=0)$p.val



