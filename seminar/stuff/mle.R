library(ggplot2)
load("section3.RData")
plot.data <- na.omit(plot.data)

pdf("./loess_regression.pdf", width=6)
ggplot(plot.data, aes(x=idh, y = lula1st.pct)) + geom_point(size=.7) + geom_smooth(size = 1.2) + scale_x_continuous("Human Development Index") + scale_y_continuous("Lula Vote Share") + stat_smooth(method="lm", se=TRUE, size  = 1.2, color="red")
dev.off()


summary(lm(lula1st.pct~1, data = plot.data))

#Fix sigma to be 1
log.lik <- c()
mus <-  seq(0,1,.001)
for(i in 1:length(mus)){
 log.lik <- c(log.lik,-1/2*sum(plot.data$lula1st.pct-mus[i])^2)
}
mus[log.lik==max(log.lik)]

ggplot(data.frame(log.lik=log.lik, mus = mus), aes(x=mus, y = log.lik)) + geom_line()


model <- (lm(lula1st.pct~bolsa.pct, data = plot.data))
summary(model)

#Specify the likelihood function
ols.lf1 <- function(beta, y=plot.data$lula1st.pct, X=model.matrix(model)) {
  sigma2 <- .01459264 #this is the sigma^2 from the regression
  if (sigma2 <= 0) return(NA)
  n <- nrow(X)
  e <- y - X%*%beta                                  # t() = matrix transpose
  logl <- ((-n/2)*log(2*pi)) - ((n/2)*log(sigma2)) - ((t(e)%*%e)/(2*sigma2))
  return(logl) # since optim() does minimisation by default.
}

library(lattice)
mle.plt <- expand.grid(b0=seq(.1,.5,.005), b1 =seq(.4,1,.005))
mle.plt$lik <- apply(mle.plt, 1, ols.lf1)
wireframe(lik~b0 * b1, mle.plt, aspect = c(1,.5), drape = T)
ggplot(mle.plt, aes(x=b0, y = b1, z = lik)) + geom_tile(aes(fill = lik)) + stat_contour() 

ols.lf1(c(.297072, .746038), plot.data$lula1st.pct, model.matrix(model))

