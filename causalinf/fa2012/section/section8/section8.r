rm(list=ls(all=TRUE))
load("~/Dropbox/projects/ps236/hw/hw3/hw3.RData")
load("~/Dropbox/projects/ps236/section/section8/plots.RData")


###
#Bootstrap the Mean and Median
###

##
#Mean
##

#Choose the number of bootstrap samples
n.b <- 500

#create the matrices to collect the bootstrap estimates
means.treat <- matrix(,nrow=n.b)
means.contr <- matrix(,nrow=n.b)
meds.treat <- matrix(,nrow=n.b)
meds.contr <- matrix(,nrow=n.b)

y.treat <- hw3.data$vote.pop[hw3.data$treat=="client"]
y.contr <- hw3.data$vote.pop[hw3.data$treat=="pub.pol"]

#loop n.b. times
for (i in 1:n.b){
  #Be sure to sample WITH REPLACEMENT
  y.bs.treat <-y.treat[sample(1:length(y.treat),length(y.treat),replace=TRUE)]
  y.bs.contr <-y.contr[sample(1:length(y.contr),length(y.contr),replace=TRUE)]
  #Collect the boot strap estimates
  means.treat[i] <- mean(y.bs.treat)
  means.contr[i] <- mean(y.bs.contr)
  meds.treat[i] <- median(y.bs.treat)
  meds.contr[i] <- median(y.bs.contr)
}

#Bootstrapped Standard Errors of the Mean
sd(means.treat)
sd(means.contr)

#what about using the standard formulas?
sd(y.treat)/sqrt(length(y.treat))
sd(y.contr)/sqrt(length(y.contr))

#Standard errors of the median
sd(meds.treat)
sd(meds.contr)

#Let's plot the distribution
hist(meds.treat, breaks="FD",col="blue",lwd=2)
hist(meds.contr, breaks="FD",col="red")


#Hypothesis testing
#Be sure to check which side you're on
mean(y.treat)-mean(y.contr)
mean.diff.bs <- means.treat-means.contr
length(sum(mean.diff.bs[mean.diff.bs<=0]))/n.b * 2 #multiply by 2 to adjust for two sided hypothesis testing
hist(means.treat-means.contr, breaks="FD", col="red")

###
#Bootstrapped Percentile Confidence Intervals
###
diff.means <- matrix(NA,nrow=n.b)
diff.medians <- matrix(NA,nrow=n.b)
diff.means.better <- matrix(NA,nrow=n.b)
for (i in 1:n.b){
  y.bs.treat <-y.treat[sample(1:length(y.treat),length(y.treat),replace=TRUE)]
  y.bs.contr <-y.contr[sample(1:length(y.contr),length(y.contr),replace=TRUE)]
  diff.means[i] <- mean(y.bs.treat) - mean(y.bs.contr)
  diff.means.better[i] <- (mean(y.bs.treat) - mean(y.bs.contr))/sqrt(var(y.bs.treat)/length(y.bs.treat) + var(y.bs.contr)/length(y.bs.contr)) - (mean(y.treat)-mean(y.contr))/sqrt(var(y.treat)/length(y.treat) + var(y.contr)/length(y.contr))
  diff.medians[i] <- median(y.bs.treat) - median(y.bs.contr)
}
hist(diff.means, breaks="FD",col="blue",lwd=2)
hist(diff.medians, breaks="FD",col="blue",lwd=2)
hist(diff.means.better, breaks="FD",col="blue",lwd=2)


#use the quantile function to return the confidence interval
quantile(diff.means,probs=c(.025,.975))
t.test(y.treat,y.contr)
quantile(diff.medians,probs=c(.025,.975))
#For better boostrap intervals do the following
quantile(diff.means.better, probs=c(.025,.975))
(mean(y.treat)-mean(y.contr)) - quantile(diff.means.better,.975)*sd(diff.means)
(mean(y.treat)-mean(y.contr)) - quantile(diff.means.better,.025)*sd(diff.means)


#Use the boot package
library(boot)
#Create a function
diff.means.stat <- function(data, indx){
  bs.data <- data[indx,]
  diff.mean <- mean(bs.data$vote.pop[bs.data$treat=="client"]) - mean(bs.data$vote.pop[bs.data$treat=="pub.pol"])
  var <- var(bs.data$vote.pop[bs.data$treat=="client"])/length(bs.data$vote.pop[bs.data$treat=="client"]) + var(bs.data$vote.pop[bs.data$treat=="pub.pol"])/length(bs.data$vote.pop[bs.data$treat=="pub.pol"])
  c(diff.mean, var)
}
#Use the boot command
boot.results <- boot(data = hw3.data, statistic = diff.means.stat, R=1000)
plot(boot.results)
#For bootstrapped confidence inervals
boot.ci(boot.results)

###
#Bootstrap Regression Model
###


#estimate the model using the "true" data
lm.model <- lm(vote.pop~treat+as.factor(block),data=hw3.data)
#we will sample from the residuals
residuals <- residuals(lm.model)

n.b <- 500
betas.bs <- matrix(nrow=9,ncol=n.b)
for (i in 1:n.b){
  #sample the residuals
  e.star <- residuals[sample(1:length(residuals),replace=TRUE)]
  #Use model.matrix to pull out the design matrix
  #Use coef to pull out the vector of model coefficients
  Y.star <- model.matrix(lm.model) %*% coef(lm.model) + e.star
  #re-estimate the model
  beta.star <- coef(lm(Y.star ~ hw3.data$treat + as.factor(hw3.data$block)))
  betas.bs[,i] <- beta.star
}

#Bias
hist(betas.bs[2,], breaks="FD",col="blue",lwd=2)
abline(v=coef(lm.model)[2],lwd=3,col="red")

coef(lm.model)[2]
mean(betas.bs[2,])

#standard error
sd(betas.bs[2,])

#confidence interval
quantile(betas.bs[2,],probs=c(.025,.975))

#variance-covariance matrix
var(betas.bs[2,])
cov(betas.bs[2,],betas.bs[1,])


#Difference between the first and second coefficient
hist(betas.bs[2,]-betas.bs[1,], breaks="FD",col="blue",lwd=2)


##
###Bootstrap Lowess
###


attach(fd.occs.inc)

bw=.6
#Lowess function requires a "band width"
lo.rain.inc <- lowess(x=stdized.yearly,y=std.basket,f=bw)
plot(stdized.std.monthly,std.basket,cex=.05,col="gray80",ylim=c(-.1,.1),xlab="Standardized Rainfall",ylab="Agricultural Income")
#To save on time, i'll keep the number of the bootstrap samples small
n.b <- 50                    
lo.b.x <- matrix(nrow=nrow(fd.occs.inc),ncol=n.b)
lo.b.y <- matrix(nrow=nrow(fd.occs.inc),ncol=n.b)
for (i in 1:n.b){
 bs.index <-  sample(1:nrow(fd.occs.inc),replace=TRUE)
 x <- stdized.yearly[bs.index]
 y <- std.basket[bs.index]
 lo.b <- lowess(x=x,y=y,f=bw)
 lines(lo.b, col="gray50")
 lo.b.x[,i] <- lo.b$x
 lo.b.y[,i] <- lo.b$y
}
lines(lo.rain.inc,col=2,lwd=2)
#need to apply a quantile to every vector of bootstrapped Y's corresponding to each value of X 

conf.int.b <- apply(lo.b.y,1,quantile,probs=c(.025,.975))

plot(stdized.std.monthly,std.basket,cex=.05,col="gray80",ylim=c(-.1,.1),xlab="Standardized Rainfall",ylab="Agricultural Income")
lines(lo.rain.inc,col=2,lwd=2)
lines(lo.rain.inc$x,conf.int.b[1,],lty=2,col=1)
lines(lo.rain.inc$x,conf.int.b[2,],lty=2,col=1)

detach(fd.occs.inc)

