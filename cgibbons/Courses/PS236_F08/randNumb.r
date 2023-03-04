### PS 236 Section, September 2008 ###
### Random Variable R Code ###

## You can use R to generate random numbers
## 'unif' gives the uniform distribution
## 'norm' gives the normal distribution
## 'r' generates random numbers
## 'p' gives quantile of a given value of the distribution
## 'q' gives the value for a given quantile of the distribution

## Let's take 1000 draws from a standard normal distribution
r <- rnorm(1000)

## Let's use these draws to create 1000 draws from a N(10,5) distribution
## If P is N(0,1), then Q ~ N(m,s) is distributed as sP + m
r10 <- 5*r + 10

## Note: Of course, you could have used R to create Q directly

## Let's get the summary statistics to see if it worked
mean(r10)
sd(r10)
sd(r10)^2
var(r10)
max(r10)
min(r10)
median(r10)

## What if we want the qth quantile?
q <- 0.5
sort(r10)[round(q*length(r10))]

## Or, more accurately
quantile(r10, probs=.1*c(0:10))

## Let's pick a value from r10 at random
v <- sample(r10,1)

## What quartile is v in?
sum(r10 <= v)/length(r10)