
rm(list=ls(all=TRUE))
set.seed(12345)
n=1000

setwd("~/Dropbox/causalinf.private/fa2014/Slides/Slides_qqplot")

### Figure 1
treat = rnorm(n,mean=0,sd=1)
control = rnorm(n,mean=0,sd=1)

pdf("fig1.pdf",width=4.5,height=4.5)
par(cex=0.7)
qqplot(treat,control,col="blue")
abline(a=0,b=1,lwd=2,col="red3") 
dev.off()

### Figure 2
treat = rnorm(n,mean=1,sd=1)
control = rnorm(n,mean=0,sd=1)

pdf("fig2.pdf",width=4.5,height=4.5)
par(cex=0.7)
qqplot(treat,control,col="blue")
abline(a=0,b=1,lwd=2,col="red3") 
dev.off()

### Figure 3
treat = rnorm(n,mean=0,sd=2)
control = rnorm(n,mean=0,sd=1)

pdf("fig3.pdf",width=4.5,height=4.5)
par(cex=0.7)
qqplot(treat,control,col="blue")
abline(a=0,b=1,lwd=2,col="red3") 
dev.off()

### Figure 4
treat = runif(n,0,1)
control = runif(n,0,1)

pdf("fig4.pdf",width=4.5,height=4.5)
par(cex=0.7)
qqplot(treat,control,col="blue")
abline(a=0,b=1,lwd=2,col="red3") 
dev.off()

### Figure 5
treat = runif(n,0,1)
control = runif(n,-0.5,1.5)

pdf("fig5a.pdf",width=4.5,height=4.5)
par(cex=0.7)
qqplot(treat,control,col="blue")
abline(a=0,b=1,lwd=2,col="red3") 
dev.off()

pdf("fig5b.pdf",width=4.5,height=4.5)
par(cex=0.7)
boxplot(treat,control,names=c("Treat","Control"))
dev.off()

### Figure 6
treat = rexp(n,1)
control = rnorm(n,mean=1,sd=1)

pdf("fig6a.pdf",width=4.5,height=4.5)
par(cex=0.7)
qqplot(treat,control,col="blue")
abline(a=0,b=1,lwd=2,col="red3") 
dev.off()

pdf("fig6b.pdf",width=4.5,height=4.5)
par(cex=0.7)
boxplot(treat,control,names=c("Treat","Control"))
dev.off()











