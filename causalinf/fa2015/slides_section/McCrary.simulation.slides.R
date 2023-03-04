#######################################################################
# McCrary (2006) section
#######################################################################

library(ggplot2)
library(rdd)

set.seed(12345)

setwd("~/Dropbox/causalinf.private/fa2015/slides/slides8_mccrary")

### Example 1
b=0.8
c=0
x=rnorm(1000,mean=0,sd=1)
g = floor((x-c)/b)*b+b/2+c

pdf(file="fig_histogram1.pdf",width=4.5,height=4.5)
par(cex=0.7)
barplot(table(g),col="lightblue",las=1)
grid()
dev.off()

### Constructing X and Y:
b=0.05
c=0
x=rnorm(1000,mean=0,sd=1)
g = floor((x-c)/b)*b+b/2+c
X = seq(min(g),max(g)+b,by=b)
n = length(x)

Y = rep(NA,length(X))
for (j in c(1:length(X))){
  Y[j] = (1/(n*b))*sum(abs(g-X[j])<0.1^6)
}

# Compare the Y from my code to the bins from the function "DCdensity"
Y.DCdensity=DCdensity(x,cutpoint=0,ext.out=TRUE,bin=b,plot=FALSE)$data[,2]

pdf(file="fig_histogram2.pdf",width=4.5,height=4.5)
ggplot(data.frame(Y.me=Y,Y.DCdensity=Y.DCdensity),
       aes(x=Y.me,y=Y.DCdensity))+
  geom_point()+geom_line()+
  labs(x="My code Y (bin height)",
       y="DCdensity Y (bin height)" )
dev.off()  

#################################
# Smoothing
#################################
x<-runif(1000,-1,1)
x<-x+2*(runif(1000,-1,1)>0 & x<0)

a=DCdensity(x,cutpoint=0,ext.out=TRUE)

### Bandwidth:
#h=0.2 #arbitrary #works
h=a$bw
b=a$binsize

c=0
g = floor((x-c)/b)*b+b/2+c
X = seq(min(g),max(g)+b,by=b)

Y = rep(NA,length(X))
for (j in c(1:length(X))){
  Y[j] = (1/(n*b))*sum(abs(g-X[j])<0.1^6)
}


# Example the cut-point: r=0
r=0

K = apply(matrix((X-r)/h,ncol=1),1,function(x){return(max(0,1-abs(x)))})  
X1 = X-r
lm.kernel = lm(Y~X1,weights=K)
f.density = coef(lm.kernel)[1]

### Density function given bandwidth:

f.kernel.lm = function(r,side,bandwidth,X,Y,cutpoint){
  
  K = apply(matrix((X-r)/bandwidth,ncol=1),1,function(x){return(max(0,1-abs(x)))})  
  X1 = (X-r)
  
  # side (below or above "c"):
  if (side=="left"){
    K1 = K[X<=cutpoint & r<=cutpoint]
    X1 = X1[X<=cutpoint & r<=cutpoint]
    Y1 = Y[X<=cutpoint & r<=cutpoint]
  }
  if (side=="right"){
    K1 = K[X>=cutpoint & r>=cutpoint]
    X1 = X1[X>=cutpoint & r>=cutpoint]
    Y1 = Y[X>=cutpoint & r>=cutpoint]
  }
  
  lm.kernel = lm(Y1~X1,weights=K1)
  return(coef(lm.kernel)[1])
}

### calculating the smoother for the values of "r":

f.density1 = rep(NA,length(x))

max.x = length(x)
for (i in c(1:max.x)){
  cat("Iteration:  ",i," out of  ",max.x,"\n")
  if (x[i]<c){
    side0="left"
  }
  if (x[i]>=c){
    side0="right"
  }
  
  f.density1[i] = f.kernel.lm(x[i],side=side0,bandwidth=h,X=X,Y=Y,cutpoint=c)
}

### Figure;

pdf(file="fig_mccrary1.pdf")
ggplot()+
  #geom_point(aes(x=x,y=f.density1))+
  geom_line(aes(x=x[x<0],y=f.density1[x<0]),lwd=1)+
  geom_line(aes(x=x[x>0],y=f.density1[x>0]),lwd=1)+
  geom_vline(x=0,col="red",lty=2,lwd=1)+
  geom_point(aes(x=X,y=Y),col="blue")+
  labs(x="Running (forcing) variable",
       y="Estimated density")
  
dev.off()

pdf(file="fig_mccrary2.pdf")
DCdensity(x,cutpoint=0,ext.out=FALSE)
dev.off()

###########################################################
# Point estimate function

f.estimate <- function(cutpoint,X,Y){
  r=cutpoint
  X.left = X[X<=cutpoint]
  Y.left = Y[X<=cutpoint]
  K.left = apply(matrix((X.left-r)/h,ncol=1),1,function(x){return(max(0,1-abs(x)))})  
  X1.left = X.left-r
  lm.kernel.left = lm(Y.left~X1.left,weights=K.left)
  
  X.right = X[X>=cutpoint]
  Y.right = Y[X>=cutpoint]
  K.right = apply(matrix((X.right-r)/h,ncol=1),1,function(x){return(max(0,1-abs(x)))})  
  X1.right = X.right-r
  lm.kernel.right = lm(Y.right~X1.right,weights=K.right)
  theta <- log(coef(lm.kernel.right)[1])-log(coef(lm.kernel.left)[1])
  return(theta)
}

theta.obs <- f.estimate(cutpoint=0,X=X,Y=Y)

# Check:
a=DCdensity(x,cutpoint=0,ext.out=TRUE,plot=FALSE)
a$theta

###########################################################
# Bootstrap vs. analytical SE of theta estimate


B=1000
boot.matrix <- matrix(sample(c(1:length(Y)),length(Y)*B,replace=TRUE),
                      ncol=B,
                      nrow=length(Y))
theta.boot <- rep(NA,B)
for (j in c(1:B)){
  index = boot.matrix[,j]
  Xb <- X[index]
  Yb <- Y[index]
  theta.boot[j] <- f.estimate(cutpoint=0,X=Xb,Y=Yb)
}

# Comparison:
a=DCdensity(x,cutpoint=c,ext.out=TRUE,plot=FALSE)
a$se
sd(theta.boot,na.rm=TRUE)

# check normality:
pdf(file="fig_mccrary3.pdf")
ggplot(data.frame(theta=theta.boot-f.estimate(cutpoint=0,X=X,Y=Y)),aes(x=theta))+
  geom_density()+
  geom_density(aes(x=rnorm(B,mean=0,sd=a$se)),col="red")+
  xlim(-1.2,1.2)
dev.off()







