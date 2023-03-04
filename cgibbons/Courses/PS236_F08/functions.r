### PS 236 Section, September 2008 ###
### Writing Programs               ###


## Basic function writing
## name <- function(x) { ...stuff... }
## 'name' is the name of the function and is called for data y 'name(y)'
## 'x' is the argument that you will pass to the function

## Let's create a new squared function
sqrd <- function(x) {
  x*x
  }

sqrd(4)
sixteen <- sqrd(4)
sixteen

## How about raising one variable to the power of another?
## Let's make the default power be 2 (i.e., squaring)
pow <- function(a,b=2) {
  a^b
  }

x <- 10
pow(x)
y <- 4
pow(x,y)

## If you associate the values you want to call with their variable name,
##   you can place them in any order
pow(b=y, a=x)

## Let's rewrite the function to give us the base, the power, and the result
pow2 <- function(a,b=2) {
  list(a,b,a^b)
  }

## Why a list rather than just 'c()'?
x <- c(1:10)
pow2(x)

## Load the following STATA file
##   A 3% random sample of the 2006 American Community Survey by the Census
##   Bureau from the IPUMS database
## First, load the 'foreign' package
## If you don't have this package, R should download it for you
library(foreign)

## Now, set the working directory --- note the direction of the slashes
## BE SURE TO SET TO YOUR COMPUTER'S WORKING DIRECTORY (where you stored the
##   data file).
setwd("C:/Users/Charles/R")

## Load the data and store it as 'X' ##
X <- read.dta("ownership.dta")

## Create a poverty indicator variable
attach(X)
pov <- ifelse(X$poverty < 100, 1, NA)
pov <- ifelse(X$ poverty >= 100, 0, pov)
detach(X)

## Let's write a function to find the poverty rate for a given state
povRate <- function(x) {
  mean(pov[X[,2] == x])
  }

povRate("California")

## Running that program for every state would be really time consuming. . .
## Let's automate it

## First, using a for loop

## A short illustration of a for loop
## for (i in 'set that i should belong to') { ... stuff ... }
w <- NA
for (i in 1:10) {
  w <- c(w, i)
  }
  
rates <- function(x) {
  w <- data.frame(cbind(NA,NA))
  names(w) <- c("State", "Poverty Rate")
  for (i in levels(x)) {
    z <- round(mean(pov[X[,2] == i])*100,2)
    z <- ifelse(z == "NaN", NA, z)
    w <- rbind(w, c(i,z))
    }
  w[,2] <- as.numeric(w[,2])
  w <- w[is.na(w[,2]) == 0,]
  w
  }
  
ratesA <- rates(X[,2])

## Now using 'sapply'
## 'sapply' vectorizes a function over a set of variables

## A short illustration of 'sapply'
## sapply('set that x should belong to', function(x) { ... stuff ... })
sapply(1:10, function(x){
  x
  } )

ratesSA <- function(x) {
  w <- sapply(levels(x), function(y) {
    z <- round(mean(pov[X[,2] == y])*100,2)
    z <- ifelse(z == "NaN", NA, z)
    z
    } )
  w <- w[is.na(w) == 0]
  }
  
ratesB <- ratesSA(X[,2])

## Now even simpler!
## Since 'state' is a factor, we can use 'tapply'
## tapply(variable of interest, factor(s) to subdivide, function to use)
ratesT <- round(tapply(pov, X[,2], mean)*100,2)
ratesT <- ratesT[!is.na(ratesT)]

## We can use 'tapply' on two factors
twoFact <- round(tapply(pov, list(X[,2],X[,3]), mean)*100,2)

## Get rid of rows that are all 'NA'
for(i in dim(twoFact)[1]:1) {
  if (sum(is.na(twoFact[i,])) == dim(twoFact)[2]) {
    twoFact <- twoFact[-c(i),]
  }
}
  
## Get rid of columns that are all 'NA'
for(i in dim(twoFact)[2]:1) {
  if (sum(is.na(twoFact[,i])) == dim(twoFact)[1]) {
    twoFact <- twoFact[,-c(i)]
  }
}
