#
# Example code for the difference pathology
#

#set the random number seed so we get the same answer every time
set.seed(39291)

#Let's generate some data to replicate the pathology I talked about.

#the number of observations:
obs  <- 10000

#draw some random normal errors with mean 0 and standard deviation equal to 1
u <- rnorm(obs,0,1);

#y equals:
y <- 5 + u;

#let's create delta y
dely1 <- y[1:(obs-1)] - y[2:obs]

#lets create the variables we are going to work with, only keep 950 of
#the obs, just to make sure we don't runoff the end of the index
dy1 <- dely1[1:(obs-50)];
dy2 <- dely1[2:(obs-49)];
dy3 <- dely1[3:(obs-48)];
dy4 <- dely1[4:(obs-47)];
dy5 <- dely1[5:(obs-46)];
dy6 <- dely1[6:(obs-45)];

print(summary(lm(dy1 ~ dy2 + dy3 +dy4 + dy5 + dy6)));

lm2 <- lm(dy1 ~ dy2 + dy3)
summary(lm2)

Y  <- dy1
X1 <- dy2
X2 <- dy3

base  <- var(X1)*var(X2)-cov(X1,X2)^2
b1  <- (cov(X1,Y)*var(X2) - cov(X2,Y)*cov(X1,X2))/base
b2  <- (cov(X2,Y)*var(X1) - cov(X1,Y)*cov(X1,X2))/base

a1  <- mean(Y) - b1*mean(X1) - b2*mean(X2);
cat("\n*****************************************\n")
cat("\nResults from doing it by hand\n")
cat("intercept :",a1,"\n")
cat("dy2       :",b1,"\n")
cat("dy3       :",b2,"\n")
cat("\n")
cat("cov(X1,Y) :",cov(X1,Y),"\n")
cat("cov(X2,Y) :",cov(X2,Y),"\n")
cat("cov(X1,X2):",cov(X1,X2),"\n")
cat("\n*****************************************\n")



