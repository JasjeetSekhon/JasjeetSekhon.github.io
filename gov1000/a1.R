par(mfrow=c(2,1))

x  <- seq(-4,4,by=.1)
plot(x,dnorm(x,0,1),type="l",ylab="Density",xlab="some variable", ylim=c(0,.4))

x  <- seq(-4,4,by=.1)
plot(x,dnorm(x,0,1.5),type="l",ylab="Density",xlab="some variable", ylim=c(0,.4))
