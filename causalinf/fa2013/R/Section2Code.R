# Section 2 Code-Regression

# First run the following experiment to show that beta.hat=tau.hat

t=rnorm(25,1,1)
c=rnorm(25,0,4)
test=t.test(t,c)
outcomes=c(t,c)
treat=c(rep(1,25),rep(0,25))
data=cbind(treat,outcomes)
data=as.data.frame(data)
colnames(data)=c("Treat","Outcome")
model=lm(Outcome~Treat,data=data)
summary(model)[2,4]

tau.hat=test$estimate[1]-test$estimate[2]

beta.hat=summary(model)$coef[2,1]

tau.hat

beta.hat


# OLS with error in X's

X=rnorm(50,3,1)
Y=X+rnorm(50, 10, 4)
data=cbind(X,Y)
data=as.data.frame(data)
model=lm(Y~X,data)

pdf("XsWithoutError.pdf")
plot(X,Y,xlim=c(0,6),ylim=c(0,25),pch=16)
abline(v=0)
abline(h=0)
abline(model,col="blue")
dev.off()


pdf("XsWithError.pdf")
plot(X,Y,xlim=c(0,6),ylim=c(0,25),pch=16)
abline(v=0)
abline(h=0)
abline(model,col="blue")
X.with.error=X+rnorm(50,0,1.5)
points(X.with.error,Y,col="gray",pch=16)
data2=cbind(X.with.error,Y)
data2=as.data.frame(data2)
model.with.X.error=lm(Y~X.with.error,data2)
abline(model.with.X.error,col="green")
legend(0.2,25,c("No Error","Error in X"),col=c("blue","green"),lty=c(1,1))
dev.off()

estimates=rep(0,10000)
for(i in 1:10000){

X.with.error=X+rnorm(50,0,2)
data2=cbind(X.with.error,Y)
data2=as.data.frame(data2)
model.with.X.error=lm(Y~X.with.error,data2)
estimates[i]=summary(model.with.X.error)$coef[2,1]
}

pdf("XErrorHistogram.pdf")
hist(estimates,xlim=c(-summary(model)$coef[2,1]-1,summary(model)$coef[2,1]+1),freq=FALSE, col="darkslategray1", main="Histogram of Estimates with Error in X")
abline(v=summary(model)$coef[2,1],col="blue")
text(summary(model)$coef[2,1]-0.2, 0.3,"True Value",srt=90,col="blue")
dev.off()



pdf("YsWithError.pdf")
plot(X,Y,xlim=c(0,6),ylim=c(0,25),pch=16)
abline(v=0)
abline(h=0)
abline(model,col="blue")
Y.with.error=Y+rnorm(50,0,2)
points(X,Y.with.error,col="gray",pch=16)
data3=cbind(X.with.error,Y)
data3=as.data.frame(data3)
model.with.Y.error=lm(Y.with.error~X,data3)
abline(model.with.Y.error,col="red")
legend(0.2,25,c("No Error","Error in Y"),col=c("blue","red"),lty=c(1,1))
dev.off()


estimates=rep(0,10000)
for(i in 1:10000){

Y.with.error=Y+rnorm(50,0,2)
data3=cbind(X.with.error,Y)
data3=as.data.frame(data3)
model.with.Y.error=lm(Y.with.error~X,data3)
estimates[i]=summary(model.with.Y.error)$coef[2,1]
}


pdf("YErrorHistogram.pdf")
hist(estimates,xlim=c(-summary(model)$coef[2,1]-1,summary(model)$coef[2,1]+1),freq=FALSE, col="darkslategray1", main="Histogram of Estimates with Error in Y")
abline(v=summary(model)$coef[2,1],col="blue")
text(summary(model)$coef[2,1]-0.2, 0.3,"True Value",srt=90,col="blue")
dev.off()

# No bias here





X=rnorm(100,0,1)
Z=rnorm(100,0,1)
Y=X+Z+rnorm(100, 0, 1)
data=cbind(Y,X,Z)
data=as.data.frame(data)
model=lm(Y~X+Z, data)


estimates=rep(0,10000)
for(i in 1:10000){

	Z.with.error=Z+rnorm(100, 0, 3)
	data=cbind(Y,X,Z.with.error)
	data=as.data.frame(data)
	model.with.Z.error=lm(Y~X+Z.with.error,data)
	estimates[i]=summary(model.with.Z.error)$coef[2,1]
	}

pdf("ZErrorHistogram.pdf")
	
hist(estimates,xlim=c(-1.5,1.5),freq=FALSE, main="Histogram of Estimates with Error in Z",col="darkslategray1")
abline(v=summary(model)$coef[2,1],col="blue")
text(summary(model)$coef[2,1]-0.05, 4,"True Value",srt=90,col="blue")

dev.off()