#### PS 236 Problem Set 1 ####
#### Fall 2008            ####
#### Solution code        ####


## Load the following STATA file
## First, load the 'foreign' package
library(foreign)

## And the 'xtable' package that makes tables for LaTeX
library(xtable)

## Now, set the working directory
setwd("C:/Users/Charles/R")

## Load the data and store it as 'X' ##
X <- read.dta("ownership.dta")


### Part a                ###
colSums(apply(X, 2, is.na))
summary(X)

attach(X)

pov <- ifelse(X$poverty < 100, 1, NA)
pov <- ifelse(X$ poverty >= 100, 0, pov)

ownInd <- ifelse(ownershd %in% c("Owned free and clear",
  "Owned with mortgage or loan"), 1, NA)
ownInd <- ifelse(ownershd %in% c("No cash rent", "With cash rent"), 0,
  ownInd)
ownInd <- as.factor(ownInd)

detach(X)


### Part b                ###
sum(is.na(ownInd))


### Part c                ###
povOwn <- tapply(pov, ownInd, mean)
xtable(rbind(povOwn[[1]], povOwn[[2]]))


### Part d                ###
##  Note that 'na.rm' means remove NAs when performing calculations. It is
##    needed here because of the NAs in ownInd
mean(pov[ownInd == 1 & X$statefip == "California"], na.rm=TRUE)
mean(pov[ownInd == 1 & X$statefip == "Oregon"], na.rm=TRUE)


### Part e                ###
incOwn <- tapply(X$hhincome, ownInd, function(x){
  z <- c(mean(x), sd(x), quantile(x))
  names(z) <- c("Mean", "St. Dev.", "0%", "25%", "50%", "75%", "100%")
  return(z)
  } )
xtable(rbind(incOwn[[1]], incOwn[[2]]))
  
  
### Part f                ###
incOwnPov <-tapply(X$hhincome[pov == 1], ownInd[pov == 1], function(x){
  z <- c(mean(x), sd(x), quantile(x))
  names(z) <- c("Mean", "St. Dev.", "0%", "25%", "50%", "75%", "100%")
  return(z)
  } )
xtable(rbind(incOwnPov[[1]], incOwnPov[[2]]))
  
  
### Part g                ###
smed <- tapply(X$hhincome, X$statefip, median)
xtable(as.matrix(smed[!is.na(smed)],ncol=1))


