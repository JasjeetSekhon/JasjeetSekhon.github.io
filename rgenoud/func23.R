# 23 functions as defined in
# Xin Yao, Yong Liu, and Guangming Lin
# Evolutionary Programming Made Faster
# IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL. 3, NO. 2, JULY 1999

testfuncs <- list()
testfuncs[[1]] <- function(x) { sum(x^2) }

testfuncs[[2]] <- function(x) { sum(abs(x)) + prod(abs(x)) }

testfuncs[[3]] <- function(x) {
  n <- length(x)
  sum(sapply(1:n, function(i){ sum(x[1:i])^2 }))
}

testfuncs[[4]] <- function(x) { max(abs(x)) }

testfuncs[[5]] <- function(x) {
  n <- length(x)-1
  sum(sapply(1:n, function(i){ 100*(x[i+1]-x[i]^2)^2 + (x[i]-1)^2 }))
}

testfuncs[[6]] <- function(x) { sum(floor(x+0.5)^2) }

testfuncs[[7]] <- function(x) {
  n <- length(x)-1
  sum(sapply(1:n, function(i){ i*x[i]^4 })) + runif(1,0,1)
}

testfuncs[[8]] <- function(x) {
  n <- length(x)
  -sum(sapply(1:n, function(i){ x[i]*sin(sqrt(abs(x[i]))) }))
}

testfuncs[[9]] <- function(x) {
  n <- length(x)
  sum(sapply(1:n, function(i){ x[i]^2 - 10*cos(2*pi*x[i]) + 10 }))
}

testfuncs[[10]] <- function(x) {
  n <- length(x)
  -20*exp(-0.2*sqrt(sum(sapply(1:n, function(i){ x[i]^2 }))/30)) -
    exp(sum(sapply(1:n, function(i){ cos(2*pi*x[i]) }))/30) + 20 + exp(1)
}

testfuncs[[11]] <- function(x) {
  n <- length(x)
  sum(sapply(1:n, function(i){ x[i]^2 }))/4000 -
    prod(sapply(1:n, function(i){ cos(x[i]/sqrt(i)) })) + 1
}

testfuncs[[12]] <- function(x) {
  n <- length(x)
  y <- 1 + (x+1)/4;
  u <- function(z,a,k,m) {
    ifelse(z>a, k*(z-a)^m, ifelse(z>=-a, 0, k*(-z-a)^m))
  }
  pi/n *(10*sin(pi*y[1])^2 +
      sum(sapply(1:(n-1), function(i){ (y[i]-1)^2*(1+10*sin(pi*y[i+1])^2) })) +
      (y[n]-1)^2) +
    sum(u(x,10,100,4))
}

testfuncs[[13]] <- function(x) {
  n <- length(x)
  u <- function(z,a,k,m) {
    ifelse(z>a, k*(z-a)^m, ifelse(z>=-a, 0, k*(-z-a)^m))
  }
  0.1 *(sin(pi*3*x[1])^2 +
      sum(sapply(1:(n-1), function(i){ (x[i]-1)^2*(1+sin(3*pi*x[i+1])) })) +
      (x[n]-1)^2)*(1+sin(2*pi*x[n])^2) +
    sum(u(x,5,100,4))
}

testfuncs[[14]] <- function(x) {
  a <- rbind(
    c(-32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32,
      -32, -16, 0, 16, 32, -32, -16, 0, 16, 32),
    c(-32, -32, -32, -32, -32, -16, -16, -16, -16, -16, 0, 0, 0, 0, 0,
       16, 16, 16, 16, 16, 32, 32, 32, 32, 32))
  1/(1/500 + sum(1/(1:25 + apply((x-a)^6,2,sum))))
}

testfuncs[[15]] <- function(x) {
  a <- c(0.1957, 0.1947, 0.1735, 0.1600, 0.0844, 0.0627, 0.0456, 0.0342,
         0.0323, 0.0235, 0.0246)
  b <- 1/c(0.25, 0.5, 1, 2, 4, 6, 8, 10, 12, 14, 16)
  sum(sapply(1:11, function(i){
    (a[i] - x[1]*(b[i]^2+b[i]*x[2])/(b[i]^2+b[i]*x[3]+x[4]))^2 }))
}

testfuncs[[16]] <- function(x) {
  4*(x[1]^2 - x[2]^2 + x[2]^4) - 2.1*x[1]^4 + x[1]^6/3 + x[1]*x[2]
}

testfuncs[[17]] <- function(x) {
  (x[2] - 5.1/(4*pi^2)*x[1]^2 + 5/pi*x[1] - 6)^2 +
    10*(1-1/(8*pi))*cos(x[1]) + 10
}

testfuncs[[18]] <- function(x) {
  (1+(x[1]+x[2]+1)^2*(19-14*x[1]+3*x[1]^2-14*x[2]+6*x[1]*x[2]+3*x[2]^2))*
  (30+(2*x[1]-3*x[2])^2*(18-32*x[1]+12*x[1]^2+48*x[2]-36*x[1]*x[2]+27*x[2]^2))
}

testfuncs[[19]] <- function(x) {
  n <- 3;
  d <- c(1, 1.2, 3, 3.2)
  a <- rbind(
    c(3, 10, 30),
    c(0.1, 10, 35),
    c(3, 10, 30),
    c(0.1, 10, 35)
  )
  p <- rbind(
    c(0.3689, 0.1170, 0.2673),
    c(0.4699, 0.4387, 0.7470),
    c(0.1091, 0.8732, 0.5547),
    c(0.038150, 0.5743, 0.8828)
  )
  -sum(sapply(1:4, function(i){ d[i]*
    exp(-sum(sapply(1:n, function(j){ a[i,j]*(x[j]-p[i,j])^2 }))) }))
}

testfuncs[[20]] <- function(x) {
  n <- 6;
  d <- c(1, 1.2, 3, 3.2)
  a <- rbind(
    c(10, 3, 17, 3.5, 1.7, 8),
    c(0.05, 10, 17, 0.1, 8, 14),
    c(3, 3.5, 1.7, 10, 17, 8),
    c(17, 8, 0.05, 10, 0.1, 14)
  )
  p <- rbind(
    c(0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886),
    c(0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991),
    c(0.2348, 0.1415, 0.3522, 0.2883, 0.3047, 0.6650),
    c(0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381)
  )
  -sum(sapply(1:4, function(i){ d[i]*
    exp(-sum(sapply(1:n, function(j){ a[i,j]*(x[j]-p[i,j])^2 }))) }))
}

shekel <- function(x,m) {
  d <- c(0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5)
  a <- rbind(
    c(4, 4, 4, 4),
    c(1, 1, 1, 1),
    c(8, 8, 8, 8),
    c(6, 6, 6, 6),
    c(3, 7, 3, 7),
    c(2, 9, 2, 9),
    c(5, 5, 3, 3),
    c(8, 1, 8, 1),
    c(6, 2, 6, 2),
    c(7, 3.6, 7, 3.6)
  )
  -sum(sapply(1:m, function(i){ 1/(crossprod(x-a[i,]) + d[i]) }))
}

testfuncs[[21]] <- function(x) { shekel(x,5) }

testfuncs[[22]] <- function(x) { shekel(x,7) }

testfuncs[[23]] <- function(x) { shekel(x,10) }

# testNparms:  number of parameters in the argument vector for each function
testNparms <- c(rep(30,13), 2,4,2,2,2,4,6,4,4,4)

# bounds specified for each function
testbounds <- list(
  matrix(c(-100,100), testNparms[1], 2, byrow=TRUE),
  matrix(c(-10,10), testNparms[2], 2, byrow=TRUE),
  matrix(c(-100,100), testNparms[3], 2, byrow=TRUE),
  matrix(c(-100,100), testNparms[4], 2, byrow=TRUE),
  matrix(c(-30,30), testNparms[5], 2, byrow=TRUE),
  matrix(c(-100,100), testNparms[6], 2, byrow=TRUE),
  matrix(c(-1.28,1.28), testNparms[7], 2, byrow=TRUE),
  matrix(c(-500,500), testNparms[8], 2, byrow=TRUE),
  matrix(c(-5.12,5.12), testNparms[9], 2, byrow=TRUE),
  matrix(c(-32,32), testNparms[10], 2, byrow=TRUE),
  matrix(c(-600,600), testNparms[11], 2, byrow=TRUE),
  matrix(c(-50,50), testNparms[12], 2, byrow=TRUE),
  matrix(c(-50,50), testNparms[13], 2, byrow=TRUE),
  matrix(c(-65.536,65.536), testNparms[14], 2, byrow=TRUE),
  matrix(c(-5,5), testNparms[15], 2, byrow=TRUE),
  matrix(c(-5,5), testNparms[16], 2, byrow=TRUE),
  matrix(c(-5,10,0,15), testNparms[17], 2, byrow=TRUE),
  matrix(c(-2,2), testNparms[18], 2, byrow=TRUE),
  matrix(c(0,1), testNparms[19], 2, byrow=TRUE),
  matrix(c(0,1), testNparms[20], 2, byrow=TRUE),
  matrix(c(0,10), testNparms[21], 2, byrow=TRUE),
  matrix(c(0,10), testNparms[22], 2, byrow=TRUE),
  matrix(c(0,10), testNparms[23], 2, byrow=TRUE)
)
