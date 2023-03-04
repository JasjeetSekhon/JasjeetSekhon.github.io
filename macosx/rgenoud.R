library(rgenoud)

#maximize a bivariate normal mixture which looks like a claw
biclaw <- function(xx) {
  mNd2 <- function(x1, x2, mu1, mu2, sigma1, sigma2, rho)
    {
      z1 <- (x1-mu1)/sigma1;
      z2 <- (x2-mu2)/sigma2;
      w <- (1.0/(2.0*pi*sigma1*sigma2*sqrt(1-rho*rho))) ;
      w <- w*exp(-0.5*(z1*z1 - 2*rho*z1*z2 + z2*z2)/(1-rho*rho)) ;
      as.double(w);
    }
  x1 <- xx[1]+1;
  x2 <- xx[2]+1;
  
  y <- (0.5*mNd2(x1,x2,0.0,0.0,1.0,1.0,0.0) +
        0.1*(mNd2(x1,x2,-1.0,-1.0,0.1,0.1,0.0) +
             mNd2(x1,x2,-0.5,-0.5,0.1,0.1,0.0) +
             mNd2(x1,x2,0.0,0.0,0.1,0.1,0.0) +
             mNd2(x1,x2,0.5,0.5,0.1,0.1,0.0) +
             mNd2(x1,x2,1.0,1.0,0.1,0.1,0.0)));
  
  as.double(y);
}

foo <- function()
  {
    biclaw1 <- genoud(biclaw, nvars=2, P9=50, pop.size=2000, max=TRUE, wait.generations=100,
                      max.generations=100, hard.generation.limit=TRUE,
                      MemoryMatrix=FALSE);
  }

print(system.time(foo()))

