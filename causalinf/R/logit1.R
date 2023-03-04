
logit <- function(XB)
  {
    exp.XB <- exp(XB)
    pred <- exp.XB/(1+exp.XB)
    return(pred)
  }

linear.predictor <- seq(-5,5,by=.1)
probs <- logit(linear.predictor)
png(file="logit_func.png")
plot(linear.predictor, probs, type="l", xlab="X*Beta", ylab="Predicted Probability", ylim=c(-.1,1.1))

lines(linear.predictor,lm(probs~linear.predictor)$fitted, type="l", col="red")

dev.off()
