library(ggplot2)
library(boot)
load("~/Dropbox/projects/ps236b/ipw/blattman.RData")

data <- data[order(data$abd),]

logit <- function(lp){
  1/(1+exp(-lp))
}

log.like <- function(beta, y, x){
  x.0 <- x[y==0,]
  x.1 <- x[y==1,]
 -1 * (sum(log(logit(x.1%*%beta))) + sum(log(1-logit(x.0%*%beta)))) #optim minimizes, hence the -1
}

#Propensity score, start with linear probability
propensity.lp <- lm(abd~C_ach+C_akw+C_ata+C_kma+C_oro+C_pad+C_paj+age, data = data)
design.matrix <- model.matrix(propensity.lp)
starting.values <- coef(propensity.lp)
logit.results <- optim(par = starting.values, fn = log.like, method = "BFGS", hessian=T, y=data$abd, x=design.matrix)

#coefficients
logit.results$par
coef(glm(abd~C_ach+C_akw+C_ata+C_kma+C_oro+C_pad+C_paj+age, data = data, family = binomial(link = "logit")))
#variance covariance matrix
solve(logit.results$hessian)
#standard error
sqrt(diag(solve(logit.results$hessian)))

data$prop.score <- logit(design.matrix%*%logit.results$par)
ggplot(data, aes(x=prop.score, color=as.factor(abd))) + geom_density(size =2)

ggplot(data, aes(x=prop.score, color=as.factor(abd))) + geom_histogram(size =1, aes(fill=as.factor(abd)))


#estimate average treatment effect
#IPW estimator
mean(((data$educ*data$abd)/data$prop.score) - ((data$educ*(1-data$abd))/(1-data$prop.score)))

#use the bootstrap to get variances
ipw.boot <- function(data, dv, indices){
  data <- data[indices,]
  dv <- data[,dv]
  mean(((dv*data$abd)/data$prop.score) - ((dv*(1-data$abd))/(1-data$prop.score)))
}
boot.ci(boot(data = data, ipw.boot, R = 500, dv = "educ"), type ="perc")



###Normalize weights to 1
sum((data$abd/data$prop.score))^(-1) *sum((data$educ*data$abd)/data$prop.score) - sum(((1-data$abd)/(1-data$prop.score)))^(-1) * sum((data$educ*(1-data$abd))/(1-data$prop.score))

##IPW ATE
#Estimate outcome equation
data$resp.cont <- predict(lm(educ~ age + fthr_ed + mthr_ed + orphan96+ hh_fthr_frm + hh_size96, data = data, subset=abd==0), newdata=data)
data$resp.tr <-predict(lm(educ~age + fthr_ed + mthr_ed + orphan96+ hh_fthr_frm + hh_size96, data = data, subset = abd==1), newdata = data)

mean(((data$age*data$abd)/data$prop.score) - ((data$age*(1-data$abd))/(1-data$prop.score)) - ((data$abd - data$prop.score)/(data$prop.score*(1-data$prop.score))) * ((1-data$prop.score) * data$resp.tr + data$prop.score*data$resp.cont))
