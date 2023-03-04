ks.Dstat <- function(x, y)
  {
    n.x <- length(x)
    n.y <- length(y)
    n <- n.x + n.y
                  
    w <- c(x, y)
    z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
    z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]        
    
    return( max(abs(z)) )
  } #ks.Dstat


x <- rnorm(100, 0, 2)
y <- rnorm(100,0,1)

ks.Dstat(x, y)
ks.test(x, y)

plot(ecdf(x), lty=1)
lines(ecdf(y), lty=2)

