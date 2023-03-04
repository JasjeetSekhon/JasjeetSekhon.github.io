################################
# OLS examples
################################

#install.packages("ggplot2")
library(ggplot2)
setwd("~/Documents/Slides_OLS") #working directory

set.seed(3)

n=200
x1 = rnorm(n,mean=10,1)
#x2 = rnorm(n,mean=2,10)
epsilon = rnorm(n,0,2)
y = 10+5*x1+epsilon

### mesurment error:
noise = rnorm(n,0,2)
x1_noise = x1+noise

### OLS estimation:
lm(y~x1_noise)


### Plot:
pdf(file="ols1.pdf",width=5,height=4)
ggplot(data.frame(x1,x1_noise,y))+theme_bw()+

  geom_point(aes(x=x1,y=y),col="blue",alpha=0.8)+  
geom_smooth(aes(x=x1,y=y),col="blue",method="lm",fill="blue",alpha=0.3)+

  geom_point(aes(x=x1_noise,y=y),col="red",alpha=0.5)+
geom_smooth(aes(x=x1_noise,y=y),col="red",method="lm",fill="red",alpha=0.3)+
  labs(title="measurment error with mean 0")
dev.off()

### mesurment error 2:
noise = rnorm(n,5,2)
x1_noise = x1+noise

### Plot:
pdf(file="ols2.pdf",width=5,height=4)
ggplot(data.frame(x1,x1_noise,y))+theme_bw()+
  
  geom_point(aes(x=x1,y=y),col="blue",alpha=0.8)+  
  geom_smooth(aes(x=x1,y=y),col="blue",method="lm",fill="blue",alpha=0.3)+
  
  geom_point(aes(x=x1_noise,y=y),col="red",alpha=0.5)+
  geom_smooth(aes(x=x1_noise,y=y),col="red",method="lm",fill="red",alpha=0.3)+
  labs(title="measurment error with mean 5")
dev.off()

### The relathionship between the measurment error and the bias of the estimator:

R=200
beta.vec1 = rep(999,R)
noiseSE.vec = seq(0,15,length=R)

for (i in c(1:R)){
  noise = rnorm(n,0,noiseSE.vec[i])
  x1_noise0 = x1+noise
  beta.vec1[i] = coef(lm(y~x1,data=data.frame(x1=x1_noise0,y=y)))[2]
}

pdf(file="ols3.pdf",width=5,height=4)
ggplot(data.frame(beta=beta.vec1,SD_noise = noiseSE.vec),aes(x=SD_noise,y=beta))+
  geom_point()+ylim(0,5)+theme_bw()+labs(x="Variance of noise")
dev.off()

### Changes in the expecation of the measurment error:
R=200
beta.vec2 = rep(999,R)
noiseE.vec = seq(0,100,length=R)

for (i in c(1:R)){
  noise = rnorm(n,noiseE.vec[i],1)
  x1_noise0 = x1+noise
  beta.vec2[i] = coef(lm(y~x1,data=data.frame(x1=x1_noise0,y=y)))[2]
}

pdf(file="ols4.pdf",width=5,height=4)
ggplot(data.frame(beta=beta.vec2,SD_noise = noiseE.vec),aes(x=SD_noise,y=beta))+
  geom_point()+ylim(0,5)+theme_bw()+labs(x="Expectation of noise")
dev.off()

cat("Cehck",mean(beta.vec1),"\n")
cat("Cehck",mean(beta.vec2),"\n")


### The relathionship between the measurment error and the distribution of the noise:

R=200
beta.vec3 = rep(999,R)
noise3.vec = seq(0.01,10,length=R)

for (i in c(1:R)){
  noise = rexp(n,noise3.vec[i])
  x1_noise0 = x1+noise
  beta.vec3[i] = coef(lm(y~x1,data=data.frame(x1=x1_noise0,y=y)))[2]
}

pdf(file="ols5.pdf",width=5,height=4)
ggplot(data.frame(beta=beta.vec3,SD_noise = noise3.vec),aes(x=SD_noise,y=beta))+
  geom_point()+ylim(0,5)+theme_bw()+labs(x="Rate parameter (variance = 1/Rate)")
dev.off()

###################################
# Measurment error in "y"
###################################

R=2000
beta.vec4 = rep(999,R)
noise4.vec = seq(0.01,10,length=R)

x1 = rnorm(n,mean=10,1)
epsilon = rnorm(n,0,2)
y = 10+5*x1+epsilon

for (i in c(1:R)){
  noise = rnorm(n,3,2)
  y_noise = y*noise
  beta.vec4[i] = coef(lm(y~x1,data=data.frame(x1=x1,y=y_noise)))[2]
}

hist(beta.vec4,breaks=40)
abline(v=15,col="red4",lwd=2)
abline(v=5,col="green4",lwd=2)

