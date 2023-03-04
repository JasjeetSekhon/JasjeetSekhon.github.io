#################################################################
## compute "panel-corrected" standard errors a la Beck and Katz
## simon jackman, dept of political science, stanford university
## june 2003
##
## edits by Jasjeet Sekhon.  http://jseekhon.fas.harvard.edu
## Feb, 2005
#################################################################

#################################################################
## function pcse
##
## inputs:
## lmobj: an object return by lm, e.g., foo <- lm(y ~ x1 + x2)
## group: unit id
##
## returns:  a k-by-k variance covariance matrix, where k is the
## number of predictors in lmobj (including the intercept)
## the "panel-corrected" standard errors are the square root of
## the diagonal of this matrix
#################################################################

pcse <- function(lmobj,group){
  e <- resid(lmobj)                                 ## ols residuals
  id <- unique(group)                               
  n <- length(id)
  J <- length(e)/n
  u <- matrix(NA,J,n)
  for(i in 1:n){
    u[,i] <- e[group==id[i]]
  }
  Sigma <- crossprod(u)/J                           ## MLE
#  cat("Cross-Unit Error Var-Covar Matrix (MLE)\n")
#  print(Sigma)

  X <- model.matrix(lmobj)                          ## X from lmobj
  k <- dim(X)[2]
  V1 <- matrix(0,k,k)
  V2 <- matrix(0,k,k)
  for(i in 1:n){                                    ## loop over units
    oki <- group==id[i]
    V1 <- V1 + crossprod(X[oki,])                   ## accumulate x-products
    for(j in 1:n){                                  ## loop again for x-correlations
      okj <- group==id[j]
      V2 <- V2 + (Sigma[i,j] * t(X[oki,])%*%X[okj,]) ## the "middle matrix"
    }
  }
  iV1 <- solve(V1)
  V <- iV1%*%V2%*%iV1
  return(V)
}


