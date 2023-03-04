
#install.packages("locfit")

rm(list=ls(all=TRUE))
library(locfit)
library(sampleSelection)
library(foreign)
library(ggplot2)
library(Hmisc)
library(rdd)
library(texreg)

setwd("")

###############################
# Kernel functions
###############################

K.tcub = function(t){
  return((1-abs(t)^3)^3*(abs(t)<1))
}

pdf("kernel2.pdf",width=4.5,height=4)
par(cex=0.7)
plot(seq(-1.5,1.5,length=1000),K.tcub(seq(-1.5,1.5,length=1000)),las=1,
     xlab="t",ylab="Kernel",main="Tri-Cube Kernel")
dev.off()

K.epan = function(t){
  return((1-t^2)*(abs(t)<1))
}

pdf("kernel3.pdf",width=4.5,height=4)
par(cex=0.7)
plot(seq(-1.5,1.5,length=1000),K.epan(seq(-1.5,1.5,length=1000)),las=1,
     xlab="t",ylab="Kernel",main="Epanechnikov Kernel")
dev.off()

K.gauss = function(t){
  return((1/sqrt(2*pi)*exp(-t^2/2)))
}

pdf("kernel4.pdf",width=4.5,height=4)
par(cex=0.7)
plot(seq(-3.5,3.5,length=1000),K.gauss(seq(-3.5,3.5,length=1000)),las=1,
     xlab="t",ylab="Kernel",main="Gaussian Kernel")
abline(v=c(-1,1),lty=2,col="red3",lwd=2)
dev.off()

### One more figure
X=seq(5,20,length=1000)

h=3
t = (7-X)/h
plot(K.gauss(t))


###############################
# Smoother for a distribution
###############################

n=10000
x=rnorm(n,mean=100,sd=10)
a <- hist(x,breaks=100,main=NA,col="lightblue",las=1)

plot(x,las=1)
hist(x,breaks=100,main=NA,col="lightblue",las=1,xlab="X")

with(a,plot(mids,counts,las=1,ylab="Counts",xlab="Mid-bins",col="blue3"))
plot(fit,add=TRUE)


data0 = data.frame(counts=a$counts,bin=a$mids)
fit <- locfit(counts~bin,data=data0,alpha=c(0.1^5,5),deg=1,kern="tcub")
predict(fit,newdata=data0)


plot(density(x,n=length(x)))
lines(density(x)$x,density(x)$y)


############################
# Lee data
############################

### Load data from data directory
rd <- read.dta("./RDReplication.dta")

use <- rd$Use == 1
close.25 <- abs(rd$DifDPct) < .25 & !is.na(rd$DifDPct)
close.5 <- abs(rd$DifDPct) < .5 & !is.na(rd$DifDPct)
close1 <- abs(rd$DifDPct) < 1 & !is.na(rd$DifDPct)
close2 <- abs(rd$DifDPct) < 2 & !is.na(rd$DifDPct)
close5 <- abs(rd$DifDPct) < 5 & !is.na(rd$DifDPct)
close10 <- abs(rd$DifDPct) < 10 & !is.na(rd$DifDPct)
close25 <- abs(rd$DifDPct) < 25 & !is.na(rd$DifDPct)

cov = c("PrvTrmsD","DWinPrv","DPctPrv")

rd$dinc <- rd$DWinPrv == 1 & !is.na(rd$DWinPrv)
rd$rinc <- rd$DWinPrv == 0 & !is.na(rd$DWinPrv)
rd$inc.margin <- ifelse(rd$dinc, rd$DifDPct, -rd$DifDPct)

### Notes:
# DifDPPrv has less missing values than DPctPrv


########################
# Polynomial regression 
########################
d = rd[close10==1  & use,]

d$DifDPct2 = d$DifDPct^2
d$DifDPct3 = d$DifDPct^3
d$DifDPct4 = d$DifDPct^4
d$DifDPct5 = d$DifDPct^5

d$right = (d$DifDPct>0)*1

d$rDifDPct=d$DifDPct*d$right
d$rDifDPct2 = d$right*d$DifDPct2
d$rDifDPct3 = d$right*d$DifDPct3
d$drDifDPct4 = d$right*d$DifDPct4
d$drDifDPct5 = d$right*d$DifDPct5

lm2 = lm(DifDPPrv~right+DifDPct+DifDPct2+
           DifDPct:right+DifDPct2:right,data=d)

lm3 = lm(DifDPPrv~right+DifDPct+DifDPct2+DifDPct3+
           DifDPct:right+DifDPct2:right+DifDPct3:right,data=d)

lm4 = lm(DifDPPrv~right+DifDPct+DifDPct2+DifDPct3+DifDPct4+
           DifDPct:right+DifDPct2:right+DifDPct3:right+DifDPct4:right,data=d)

lm5 = lm(DifDPPrv~right+DifDPct+DifDPct2+DifDPct3+DifDPct4+DifDPct5+
           DifDPct:right+DifDPct2:right+DifDPct3:right+DifDPct4:right+DifDPct5:right,data=d)

texreg(list(lm2,lm3,lm4,lm5),custom.model.names = paste("Polynomials",c(2:5)))


### Bins: ###
# Calculating means by bins
n.bins=40
#breaks0 = unique(c(seq(min(d$DifDPct),0,length=n.bins/2),seq(0,max(d$DifDPct),length=n.bins/2)))
breaks0 = seq(-10,10,by=0.5)
bins0 = cut(d$DifDPct,breaks=breaks0)
bins = tapply(d$DifDPct,bins0,mean)
bin.DifDPPrv = tapply(d$DifDPPrv,bins0,function(x){return(mean(x,na.rm=TRUE))})

######################
### Plots:
######################

### Fig poly 2
xx = lm2$model$DifDPct[order(lm2$model$DifDPct)]
yy = lm2$fitt[order(lm2$model$DifDPct)]

pdf("fig_poly2.pdf",width=4.5,height=4.5)
par(cex=0.7)
plot(xx[xx<0],yy[xx<0],las=1,col="red4",xlab="Runing variable (Dem Margin of victory)",
     ylab="Dem vote share previous elections",type="l",xlim=c(min(xx),max(xx)),ylim=c(min(yy),max(yy)),
     lwd=2,main="Polynomial of degree 2")
lines(xx[xx>0],yy[xx>0],col="blue4",lwd=2)
abline(v=0,lwd=2,col=1,lty=2)
points(bins,bin.DifDPPrv,type="p",las=1,pch=16,col=ifelse(bins<=0,"green4","purple3"))
dev.off()


### Fig poly 3
xx = lm3$model$DifDPct[order(lm3$model$DifDPct)]
yy = lm3$fitt[order(lm3$model$DifDPct)]

pdf("fig_poly3.pdf",width=4.5,height=4.5)
par(cex=0.7)
plot(xx[xx<0],yy[xx<0],las=1,col="red4",xlab="Runing variable (Dem Margin of victory)",
     ylab="Dem vote share previous elections",type="l",xlim=c(min(xx),max(xx)),ylim=c(min(yy),max(yy)),
     lwd=2,main="Polynomial of degree 3")
lines(xx[xx>0],yy[xx>0],col="blue4",lwd=2)
abline(v=0,lwd=2,col=1,lty=2)
points(bins,bin.DifDPPrv,type="p",las=1,pch=16,col=ifelse(bins<=0,"green4","purple3"))
dev.off()


### Fig poly 4
xx = lm4$model$DifDPct[order(lm4$model$DifDPct)]
yy = lm4$fitt[order(lm4$model$DifDPct)]

pdf("fig_poly4.pdf",width=4.5,height=4.5)
par(cex=0.7)
plot(xx[xx<0],yy[xx<0],las=1,col="red4",xlab="Runing variable (Dem Margin of victory)",
     ylab="Dem vote share previous elections",type="l",xlim=c(min(xx),max(xx)),ylim=c(min(yy),max(yy)),
     lwd=2,main="Polynomial of degree 4")
lines(xx[xx>0],yy[xx>0],col="blue4",lwd=2)
abline(v=0,lwd=2,col=1,lty=2)
points(bins,bin.DifDPPrv,type="p",las=1,pch=16,col=ifelse(bins<=0,"green4","purple3"))
dev.off()

### Fig poly 5
xx = lm5$model$DifDPct[order(lm5$model$DifDPct)]
yy = lm5$fitt[order(lm5$model$DifDPct)]

pdf("fig_poly5.pdf",width=4.5,height=4.5)
par(cex=0.7)
plot(xx[xx<0],yy[xx<0],las=1,col="red4",xlab="Runing variable (Dem Margin of victory)",
     ylab="Dem vote share previous elections",type="l",xlim=c(min(xx),max(xx)),ylim=c(min(yy),max(yy)),
     lwd=2,main="Polynomial of degree 5")
lines(xx[xx>0],yy[xx>0],col="blue4",lwd=2)
abline(v=0,lwd=2,col=1,lty=2)
points(bins,bin.DifDPPrv,type="p",las=1,pch=16,col=ifelse(bins<=0,"green4","purple3"))
dev.off()




### Replicating Lee (2008) figures:
d0 <- read.dta("./table3group100.dta")

for (window in c(1,0.5,0.25,0.1)){
  d=d0[d0$use==1,]
  d = d[!is.na(d$difdemshare),]
  d = d[abs(d$difdemshare)<window,]
  
  d$right = (d$difdemshare>0)*1
  
  poly.runing = poly(d$difdemshare,degree=4)
  x.poly = data.frame(poly.runing,poly.runing*d$right)
  
  lm4.lee = lm(demshareprev~right+difdemshare+difdemshare2+difdemshare3+difdemshare4+
                 difdemshare:right+difdemshare2:right+difdemshare3:right+difdemshare4:right,data=d)
  
  # figure
  xx = lm4.lee$model$difdemshare[order(lm4.lee$model$difdemshare)]
  yy = lm4.lee$fitt[order(lm4.lee$model$difdemshare)]
  
  pdf(paste("fig_lee",100*window,".pdf",sep=""),width=4.5,height=4.5)
  par(cex=0.7)
  plot(xx[xx<0],yy[xx<0],las=1,col="red4",xlab="Runing variable (Dem Margin of victory)",
       ylab="Dem vote share previous elections",type="l",xlim=c(min(xx),max(xx)),ylim=c(min(yy),max(yy)),
       lwd=2,main=paste("Window around the cut off : ","[",-window,",",window,"]",sep=""))
  lines(xx[xx>0],yy[xx>0],col="blue4",lwd=2)
  abline(v=0,lwd=2,col=1,lty=2)
  dev.off()
}






########################
# Looking on bins plots 
########################

d = rd[close10==1  & use,]

fit0 <- locfit(DifDPPrv~DifDPct,data=d[d$DifDPct<0,],alpha=c(0.1^5,0.5),deg=1,kern="tcub")
fit1 <- locfit(DifDPPrv~DifDPct,data=d[d$DifDPct>0,],alpha=c(0.1^5,0.5),deg=1,kern="tcub")

plot(fit0,get.data=TRUE)

loess0 = predict(fit0,newdata=d[d$DifDPct<0,])
loess1 = predict(fit1,newdata=d[d$DifDPct>0,])

xx = d$DifDPct[d$DifDPct<0]
plot(xx[order(xx)],loess0[order(xx)],type="l",xlim=c(-1,1),ylim=c(min(loess0),max(loess1)))

xx = d$DifDPct[d$DifDPct>0]
lines(xx[order(xx)],loess1[order(xx)],type="l",xlim=c(-1,1))

########################
# Looking on bins plots 
########################

d = rd[close10==1 & use,]

# Calculating means by bins
n.bins=40
#breaks0 = unique(c(seq(min(d$DifDPct),0,length=n.bins/2),seq(0,max(d$DifDPct),length=n.bins/2)))
breaks0 = seq(-10,10,by=0.5)
bins0 = cut(d$DifDPct,breaks=breaks0)
bin.DPctPrv = tapply(d$DPctPrv,bins0,function(x){return(mean(x,na.rm=TRUE))})
bin.DifDPPrv = tapply(d$DifDPPrv,bins0,function(x){return(mean(x,na.rm=TRUE))})
bins = tapply(d$DifDPct,bins0,mean)

bins00 = cut(d$inc.margin,breaks=breaks0,right=FALSE)
bin.inc.density = tapply(d$inc.margin,bins00,length)


plot(bins,bin.DPctPrv,type="p",las=1,pch=16,col=ifelse(abs(bins)<=2,"red3","green4"))
abline(v=0,col="blue3",lwd=2)

plot(bins,bin.DifDPPrv,type="p",las=1,pch=16,col=ifelse(abs(bins)<=2,"red3","green4"))
abline(v=0,col="blue3",lwd=2)

plot(bins,bin.density,type="p",las=1,pch=16,col=ifelse(abs(bins)<=2,"red3","green4"))
abline(v=0,col="blue3",lwd=2)


################################
# Previous vote share margin
################################
d = rd[close10==1 & use,]

bandwidth = IKbandwidth(d$DifDPct,d$DifDPPrv,cutpoint=0)
#CCT.bandwidth=rdbwselect(y=d$DifDPPrv,x=d$DifDPct,c=0)

n.bins=40
breaks0 = seq(-10,10,by=0.5)
bins0 = cut(d$DifDPct,breaks=breaks0)
bin.DifDPPrv = tapply(d$DifDPPrv,bins0,function(x){return(mean(x,na.rm=TRUE))})
bins = tapply(d$DifDPct,bins0,mean)

fit0 <- locfit(DifDPPrv~DifDPct,data=d[d$DifDPct<0,],alpha=c(0.1^5,bandwidth),deg=1,kern="tcub")
fit1 <- locfit(DifDPPrv~DifDPct,data=d[d$DifDPct>0,],alpha=c(0.1^5,bandwidth),deg=1,kern="tcub")

loess0.fit = predict(fit0,newdata=d[d$DifDPct<0,],se.fit=TRUE)
loess1.fit = predict(fit1,newdata=d[d$DifDPct>0,],se.fit=TRUE)

loess0=loess0.fit$fit
loess1=loess1.fit$fit

loess0.se=loess0.fit$se.fit
loess1.se=loess1.fit$se.fit

pdf("fig_loess1.pdf",width=4.5,height=4.5)
par(cex=0.7)
plot(bins,bin.DifDPPrv,type="p",las=1,pch=16,col=ifelse(bins<=0,"red3","green4"),
     ylim = c(min(d$DifDPct,min(loess0-loess0.se)),max(d$DifDPct,max(loess1+loess1.se))),
     xlab="Democratic margin",ylab="Democratic margin in previous election",cex.lab=1.2)
xx = d$DifDPct[d$DifDPct>0]
lines(xx[order(xx)],loess1[order(xx)],type="l",xlim=c(-1,1),col="green4",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(loess1[order(xx)]+loess1.se[order(xx)],rev(loess1[order(xx)]-loess1.se[order(xx)])),
        col=rgb(0, 0.5, 0,0.2),border=NA)

xx = d$DifDPct[d$DifDPct<0]
lines(xx[order(xx)],loess0[order(xx)],type="l",xlim=c(-1,1),col="red3",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(loess0[order(xx)]+loess0.se[order(xx)],rev(loess0[order(xx)]-loess0.se[order(xx)])),
        col=rgb(0.5, 0, 0,0.2),border=NA)
abline(v=0,col="blue3",lwd=2, lty=2)
dev.off()


### Not accurate due to "span" parameter...

dp = data.frame(DifDPct=d$DifDPct,Treat=(d$DifDPct>0)*1,
                loess = c(loess0,loess1),DifDPPrv=d$DifDPPrv)

### plot ith span=0.5
span0=0.5
pdf("fig_loess2.pdf",width=4.5,height=4.5)
ggplot(data=dp, aes(x = DifDPct, y = DifDPPrv)) + 
  geom_smooth(data=subset(dp,Treat==0),aes(x = DifDPct, y = DifDPPrv),
              method = "loess", span=span0,se = TRUE, size = 1,col="blue",fill="blue",alpha=0.3)+
geom_smooth(data=subset(dp,Treat==1),aes(x = DifDPct, y = DifDPPrv),
            method = "loess",span=span0, se = TRUE, size = 1,col="red",fill="red",alpha=0.3)+
theme_bw()+ # White background
geom_vline(x=0,lty=2)+
  labs(y="Democratic margin in previous election",x="Democratic margin")
dev.off()

### plot ith span=0.25
span0=0.25
pdf("fig_loess3.pdf",width=4.5,height=4.5)
ggplot(data=dp, aes(x = DifDPct, y = DifDPPrv)) + 
  geom_smooth(data=subset(dp,Treat==0),aes(x = DifDPct, y = DifDPPrv),
              method = "loess", span=span0,se = TRUE, size = 1,col="blue",fill="blue",alpha=0.3)+
  geom_smooth(data=subset(dp,Treat==1),aes(x = DifDPct, y = DifDPPrv),
              method = "loess",span=span0, se = TRUE, size = 1,col="red",fill="red",alpha=0.3)+
  theme_bw()+ # White background
  geom_vline(x=0,lty=2)+
  labs(y="Democratic margin in previous election",x="Democratic margin")
dev.off()

### plot ith span=0.1
span0=0.1
pdf("fig_loess4.pdf",width=4.5,height=4.5)
ggplot(data=dp, aes(x = DifDPct, y = DifDPPrv)) + 
  geom_smooth(data=subset(dp,Treat==0),aes(x = DifDPct, y = DifDPPrv),
              method = "loess", span=span0,se = TRUE, size = 1,col="blue",fill="blue",alpha=0.3)+
  geom_smooth(data=subset(dp,Treat==1),aes(x = DifDPct, y = DifDPPrv),
              method = "loess",span=span0, se = TRUE, size = 1,col="red",fill="red",alpha=0.3)+
  theme_bw()+ # White background
  geom_vline(x=0,lty=2)+
  labs(y="Democratic margin in previous election",x="Democratic margin")
dev.off()



######################################
# Bootstrap confidence interval
######################################

d = rd[close10==1 & use,]

bandwidth = IKbandwidth(d$DifDPct,d$DifDPPrv,cutpoint=0)

n.bins=40
breaks0 = seq(-10,10,by=0.5)
bins0 = cut(d$DifDPct,breaks=breaks0)
bin.DifDPPrv = tapply(d$DifDPPrv,bins0,function(x){return(mean(x,na.rm=TRUE))})
bins = tapply(d$DifDPct,bins0,mean)

fit0 <- locfit(DifDPPrv~DifDPct,data=d[d$DifDPct<0,],alpha=c(0.1^5,bandwidth),deg=1,kern="tcub")
fit1 <- locfit(DifDPPrv~DifDPct,data=d[d$DifDPct>0,],alpha=c(0.1^5,bandwidth),deg=1,kern="tcub")

loess0 = predict(fit0,newdata=d[d$DifDPct<0,],se.fit=FALSE)
loess1 = predict(fit1,newdata=d[d$DifDPct>0,],se.fit=FALSE)

#### SE calculation using bootstrap: - Ask Jas!!!
L0 = 200  #boootstrap replications

boot0 = matrix(NA,nrow=length(loess0.fit),ncol=L0)
boot1 = matrix(NA,nrow=length(loess1.fit),ncol=L0)

for (j in c(1:L0)){
  if(j%%10==0){cat("Iteration",j,"\n")}
  id0 = sample(c(1:nrow(d)),nrow(d),replace=TRUE)
  d0 = d[id0,]
  
  bandwidth0 = IKbandwidth(d0$DifDPct,d0$DifDPPrv,cutpoint=0)
  
  fit0 <- locfit(DifDPPrv~DifDPct,data=d0[d0$DifDPct<0,],alpha=c(0.1^5,bandwidth0),deg=1,kern="tcub")
  fit1 <- locfit(DifDPPrv~DifDPct,data=d0[d0$DifDPct>0,],alpha=c(0.1^5,bandwidth0),deg=1,kern="tcub")
  
  boot0[,j] = predict(fit0,newdata=d[d$DifDPct<0,])
  boot1[,j] = predict(fit1,newdata=d[d$DifDPct>0,])
}

min.ci0 = apply(boot0,1,function(x){return(quantile(x,prob=0.025))})
max.ci0 = apply(boot0,1,function(x){return(quantile(x,prob=0.975))})

min.ci1 = apply(boot1,1,function(x){return(quantile(x,prob=0.025))})
max.ci1 = apply(boot1,1,function(x){return(quantile(x,prob=0.975))})

### Generating the figure:

pdf("fig_loess5.pdf",width=4.5,height=4)
par(cex=0.7)
plot(bins,bin.DifDPPrv,type="p",las=1,pch=16,col=ifelse(bins<=0,"red3","green4"),
     ylim = c(min(d$DifDPct,min(min.ci0)),max(d$DifDPct,max(max.ci1))),
     xlab="Democratic margin",ylab="Democratic margin in previous election",cex.lab=1.2)
xx = d$DifDPct[d$DifDPct>0]
lines(xx[order(xx)],loess1[order(xx)],type="l",xlim=c(-1,1),col="green4",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(max.ci1[order(xx)],rev(min.ci1[order(xx)])),
        col=rgb(0, 0.5, 0,0.2),border=NA)

xx = d$DifDPct[d$DifDPct<0]
lines(xx[order(xx)],loess0[order(xx)],type="l",xlim=c(-1,1),col="red3",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(max.ci0[order(xx)],rev(min.ci0[order(xx)])),
        col=rgb(0.5, 0, 0,0.2),border=NA)

abline(v=0,col="blue3",lwd=2, lty=2)
dev.off()


##################################################
# A comparison figure: Bootstrap vs Asymptotic SE
##################################################

loess0.fit = predict(fit0,newdata=d[d$DifDPct<0,],se.fit=TRUE)
loess1.fit = predict(fit1,newdata=d[d$DifDPct>0,],se.fit=TRUE)

loess0=loess0.fit$fit
loess1=loess1.fit$fit

loess0.se=loess0.fit$se.fit
loess1.se=loess1.fit$se.fit

pdf("fig_loess6.pdf",width=5.5,height=4.5)
### Figure analytical
par(mfrow=c(1,2),cex=0.7)
plot(bins,bin.DifDPPrv,type="p",las=1,pch=16,col=ifelse(bins<=0,"red3","green4"),
     ylim = c(min(d$DifDPct,min(min.ci0)),max(d$DifDPct,max(max.ci1))),
     xlab="Democratic margin",ylab="Democratic margin in previous election",cex.lab=1.2)
xx = d$DifDPct[d$DifDPct>0]
lines(xx[order(xx)],loess1[order(xx)],type="l",xlim=c(-1,1),col="green4",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(loess1[order(xx)]+loess1.se[order(xx)],rev(loess1[order(xx)]-loess1.se[order(xx)])),
        col=rgb(0, 0.5, 0,0.2),border=NA)

xx = d$DifDPct[d$DifDPct<0]
lines(xx[order(xx)],loess0[order(xx)],type="l",xlim=c(-1,1),col="red3",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(loess0[order(xx)]+loess0.se[order(xx)],rev(loess0[order(xx)]-loess0.se[order(xx)])),
        col=rgb(0.5, 0, 0,0.2),border=NA)
abline(v=0,col="blue3",lwd=2, lty=2)


### Figure bootstrap

plot(bins,bin.DifDPPrv,type="p",las=1,pch=16,col=ifelse(bins<=0,"red3","green4"),
     ylim = c(min(d$DifDPct,min(min.ci0)),max(d$DifDPct,max(max.ci1))),
     xlab="Democratic margin",ylab="",cex.lab=1.2)
xx = d$DifDPct[d$DifDPct>0]
lines(xx[order(xx)],loess1[order(xx)],type="l",xlim=c(-1,1),col="green4",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(max.ci1[order(xx)],rev(min.ci1[order(xx)])),
        col=rgb(0, 0.5, 0,0.2),border=NA)

xx = d$DifDPct[d$DifDPct<0]
lines(xx[order(xx)],loess0[order(xx)],type="l",xlim=c(-1,1),col="red3",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(max.ci0[order(xx)],rev(min.ci0[order(xx)])),
        col=rgb(0.5, 0, 0,0.2),border=NA)

abline(v=0,col="blue3",lwd=2, lty=2)
dev.off()

##################################################
# Binary outcome and local logistic regression
##################################################

d = rd[close10==1 & use,]

bandwidth = IKbandwidth(d$DifDPct,d$DWinPrv,cutpoint=0)
#CCT.bandwidth=rdbwselect(y=d$DifDPPrv,x=d$DifDPct,c=0)

n.bins=40
breaks0 = seq(-10,10,by=0.5)
bins0 = cut(d$DifDPct,breaks=breaks0)
bin.DWinPrv = tapply(d$DWinPrv,bins0,function(x){return(mean(x,na.rm=TRUE))})
bins = tapply(d$DifDPct,bins0,mean)


### Using a linear probability model

fit0 <- locfit(DWinPrv~DifDPct,data=d[d$DifDPct<0,],alpha=c(0.1^5,bandwidth),
               deg=1,kern="tcub")
fit1 <- locfit(DWinPrv~DifDPct,data=d[d$DifDPct>0,],alpha=c(0.1^5,bandwidth),
               deg=1,kern="tcub")

loess0.fit = predict(fit0,newdata=d[d$DifDPct<0,],se.fit=TRUE)
loess1.fit = predict(fit1,newdata=d[d$DifDPct>0,],se.fit=TRUE)

loess0=loess0.fit$fit
loess1=loess1.fit$fit

loess0.se=loess0.fit$se.fit
loess1.se=loess1.fit$se.fit

pdf("fig_loess_logit1.pdf",width=4.5,height=4.5)
par(cex=0.7)
plot(bins,bin.DWinPrv,type="p",las=1,pch=16,col=ifelse(bins<=0,"red3","green4"),
     ylim = c(min(min(loess0-loess0.se)),max(max(loess1+loess1.se))),
     xlab="Democratic margin",ylab="Democratic won in previous election",cex.lab=1.2)
xx = d$DifDPct[d$DifDPct>0]
lines(xx[order(xx)],loess1[order(xx)],type="l",xlim=c(-1,1),col="green4",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(loess1[order(xx)]+loess1.se[order(xx)],rev(loess1[order(xx)]-loess1.se[order(xx)])),
        col=rgb(0, 0.5, 0,0.2),border=NA)

xx = d$DifDPct[d$DifDPct<0]
lines(xx[order(xx)],loess0[order(xx)],type="l",xlim=c(-1,1),col="red3",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(loess0[order(xx)]+loess0.se[order(xx)],rev(loess0[order(xx)]-loess0.se[order(xx)])),
        col=rgb(0.5, 0, 0,0.2),border=NA)
abline(v=0,col="blue3",lwd=2, lty=2)
dev.off()


### Using a Logistic regression
fit0 <- locfit(DWinPrv~DifDPct,data=d[d$DifDPct<0,],alpha=c(0.1^5,bandwidth),
               deg=1,kern="tcub",family="binomial",link="logit")
fit1 <- locfit(DWinPrv~DifDPct,data=d[d$DifDPct>0,],alpha=c(0.1^5,bandwidth),
               deg=1,kern="tcub",family="binomial",link="logit")

loess0.fit = predict(fit0,newdata=d[d$DifDPct<0,],se.fit=TRUE)
loess1.fit = predict(fit1,newdata=d[d$DifDPct>0,],se.fit=TRUE)

loess0=loess0.fit$fit
loess1=loess1.fit$fit

loess0.se=loess0.fit$se.fit
loess1.se=loess1.fit$se.fit

pdf("fig_loess_logit2.pdf",width=4.5,height=4.5)
par(cex=0.7)
plot(bins,bin.DWinPrv,type="p",las=1,pch=16,col=ifelse(bins<=0,"red3","green4"),
     ylim = c(min(min(loess0-loess0.se)),max(max(loess1+loess1.se))),
     xlab="Democratic margin",ylab="Democratic won in previous election",cex.lab=1.2)
xx = d$DifDPct[d$DifDPct>0]
lines(xx[order(xx)],loess1[order(xx)],type="l",xlim=c(-1,1),col="green4",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(loess1[order(xx)]+loess1.se[order(xx)],rev(loess1[order(xx)]-loess1.se[order(xx)])),
        col=rgb(0, 0.5, 0,0.2),border=NA)

xx = d$DifDPct[d$DifDPct<0]
lines(xx[order(xx)],loess0[order(xx)],type="l",xlim=c(-1,1),col="red3",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(loess0[order(xx)]+loess0.se[order(xx)],rev(loess0[order(xx)]-loess0.se[order(xx)])),
        col=rgb(0.5, 0, 0,0.2),border=NA)
abline(v=0,col="blue3",lwd=2, lty=2)
dev.off()

###############################################################
# Bootstrap confidence interval - Logistic local regression
###############################################################

d = rd[close10==1 & use,]

bandwidth = IKbandwidth(d$DifDPct,d$DWinPrv,cutpoint=0)

n.bins=40
breaks0 = seq(-10,10,by=0.5)
bins0 = cut(d$DifDPct,breaks=breaks0)
bin.DWinPrv = tapply(d$DWinPrv,bins0,function(x){return(mean(x,na.rm=TRUE))})
bins = tapply(d$DifDPct,bins0,mean)

fit0 <- locfit(DWinPrv~DifDPct,data=d[d$DifDPct<0,],alpha=c(0.1^5,bandwidth),
               deg=1,kern="tcub",family="binomial",link="logit")
fit1 <- locfit(DWinPrv~DifDPct,data=d[d$DifDPct>0,],alpha=c(0.1^5,bandwidth),
               deg=1,kern="tcub",family="binomial",link="logit")


loess0 = predict(fit0,newdata=d[d$DifDPct<0,],se.fit=FALSE)
loess1 = predict(fit1,newdata=d[d$DifDPct>0,],se.fit=FALSE)

#### SE calculation using bootstrap: - Ask Jas!!!
L0 = 600  #boootstrap replications

boot0 = matrix(NA,nrow=length(loess0),ncol=L0)
boot1 = matrix(NA,nrow=length(loess1),ncol=L0)

for (j in c(1:L0)){
  if(j%%10==0){cat("Iteration",j,"\n")}
  id0 = sample(c(1:nrow(d)),nrow(d),replace=TRUE)
  d0 = d[id0,]
  
  bandwidth0 = IKbandwidth(d0$DifDPct,d0$DWinPrv,cutpoint=0)
  
  fit0 <- locfit(DWinPrv~DifDPct,data=d0[d0$DifDPct<0,],alpha=c(0.1^5,bandwidth0),
                 deg=1,kern="tcub",family="binomial",link="logit")
  fit1 <- locfit(DWinPrv~DifDPct,data=d0[d0$DifDPct>0,],alpha=c(0.1^5,bandwidth0),
                 deg=1,kern="tcub",family="binomial",link="logit")
  
  boot0[,j]  = predict(fit0,newdata=d[d$DifDPct<0,])
  boot1[,j]  = predict(fit1,newdata=d[d$DifDPct>0,])
}

min.ci0 = apply(boot0,1,function(x){return(quantile(x,prob=0.025))})
max.ci0 = apply(boot0,1,function(x){return(quantile(x,prob=0.975))})

min.ci1 = apply(boot1,1,function(x){return(quantile(x,prob=0.025))})
max.ci1 = apply(boot1,1,function(x){return(quantile(x,prob=0.975))})

### Generating the figure:

pdf("fig_loess_logit3.pdf",width=4.5,height=4)
par(cex=0.7)
plot(bins,bin.DWinPrv,type="p",las=1,pch=16,col=ifelse(bins<=0,"red3","green4"),
     ylim = c(min(min(min.ci0)),max(max(max.ci1))),
     xlab="Democratic margin",ylab="Democratic won in previous election",cex.lab=1.2)
xx = d$DifDPct[d$DifDPct>0]
lines(xx[order(xx)],loess1[order(xx)],type="l",xlim=c(-1,1),col="green4",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(max.ci1[order(xx)],rev(min.ci1[order(xx)])),
        col=rgb(0, 0.5, 0,0.2),border=NA)

xx = d$DifDPct[d$DifDPct<0]
lines(xx[order(xx)],loess0[order(xx)],type="l",xlim=c(-1,1),col="red3",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(max.ci0[order(xx)],rev(min.ci0[order(xx)])),
        col=rgb(0.5, 0, 0,0.2),border=NA)

abline(v=0,col="blue3",lwd=2, lty=2)
dev.off()

##################################################################
# MP's for sale example
##################################################################

rm(list = ls())

setwd("")

library(foreign)
library(ggplot2)
library(Hmisc)
library(plyr)
#library(equivalence)
library(dummies)

balance.plot0 =1
if (balance.plot0==1){
  
  # Function to create figure with both summary statistics and plots of p-values for different variables across different groups
  # Author: Rocio Titiunik
  # Date: September 5th, 2008
  # Version: 1.0
  
  # Version to distribute publicly
  
  # NOTE: This function is *very far* from being a real function. As you'll see, there are many parameters set inside the function
  #       that should be arguments to the function instead. I'm currently working on a version that is fully flexible and does not
  #       set any parameters inside the function. But for now, this is it.
  
  # About the arguments of the function:
  
  # 'results':     a matrix whose rows are different variables; whose first two columns contain the means for treated and control;
  #                and whose remaining columns have the pvalues to be plotted for every variable
  
  # 'title':       title of the overall graph
  
  # at1, at2,at3:  scalars which indicates where to locate the three differents groups (mean treatment, mean controls, graph area) in the figure area
  
  # xlim1 :         the left limit of the x-axis; right limit is always set to 1
  
  # textsize:      scalar indicating the size of text in the figure
  
  # legend:        logical indicating whether the legend should be included
  
  # legendx:       scalar indicating the x-coordinate of the legend's location
  
  # legendy:       scalar indicating the y-coordinate of the legend's location
  
  # parcex:        scalar setting cex parameter
  
  
  
  plot.pval <- function(results, title=NULL, legend,legendx=0.7,legendy=3, textsize=0.9, parcex=0.8, at1=-0.35, at2=-0.15, at3=-0.9,xlim1=-0.85) {
    
    
    # set values of different parameters
    xlim = c(xlim1,1); pchset = c(21,24,22,23); pchcolset = c("blue","yellow","red","darkgreen")
    
    # set margins and letter size
    par(cex=parcex, mai = c(0.5, 0.35, 1.1, 0.35))
    
    # set number of rows 
    ny = nrow(results)
    
    # create the empty figure
    if(!is.null(title))  plot(x=NULL,axes=F, xlim=xlim, ylim=c(1,ny),xlab="",ylab="", main=title)
    if(is.null(title))   plot(x=NULL,axes=F, xlim=xlim, ylim=c(1,ny),xlab="",ylab="")
    
    # add the 0, 0.05 and 0.1 vertical lines
    abline(v=c(0,0.05,0.1),lty=c(1,4,4), lwd=c(1,2,2))
    axis(side=1,at=c(0,0.05,0.1,1),tick=TRUE, las=2, cex.axis=0.7)
    
    # add labels on top of the three areas of the graph
    axis(side=3,at=at1,labels="Mean\nTreated",tick=FALSE, padj=0.5,cex.axis=textsize)
    axis(side=3,at=at2,labels="Mean\nControl",tick=FALSE, padj=0.5,cex.axis=textsize)
    axis(side=3,at=0.5,labels="P-values",tick=FALSE, padj=0.5,cex.axis=textsize)
    
    # Fill the figure with the information which is inside the 'results' matrix
    # First, add the p-values as points
    for(i in 4:ncol(results)) points(results[,i],ny:1, pch = pchset[i-4+1], col = pchcolset[i-4+1], bg = pchcolset[i-4+1])
    
    # Second, add each variable name and the means for treated and control
    for(i in 1:ny) {
      text(at3,ny-i+1,results[i,1],adj = 0,cex=textsize) # variable name
      text(at1,ny-i+1,results[i,2], cex=textsize)        # treatment mean
      text(at2,ny-i+1,results[i,3], cex=textsize)        # control mean
    }
    
    # Add dotted horizontal lines every two variables to make it prettier
    for(i in seq(2,by=2,length.out=floor((ny-1)/2))) abline(h = i+0.5, lty = 3)
    
    # Add legend
    if(legend) legend(x=legendx, y=legendy, c("T-test","Wilcoxon","KS"), pch=pchset, pt.bg = pchcolset, cex=0.8)
  }
  
}


# loading the data
dd <- read.dta("./mps.dta")
dd$treated = factor((dd$wobl_margin>=0)*1)

setwd("~/Dropbox/causalinf.private/fa2014/Slides/Slides_loess")
#################
# Balance around a window according to: equivalence and KS tests 

### Covariates:
#generating age variable
dd$age = dd$xxyod-dd$xxyob
dd$treat = (dd$wobl_margin>0)*1

# Note: variable "labour" - party indicator - not included in the covariates...

cov.occ = c("xxoc_local_politics","xxoc_civil_serv","xxoc_local_politics",
            "xxoc_business","xxoc_white_collar","xxoc_journalist")

cov.all = c("xxfemale","xxyob","xxyod","age","wobl_race_effective_number_of_ca","wobl_race_electorate",
            "prev_attempts","wobl_race_turnout","pre_wobl","labour")


f.stat = function(x1,x0){
  #sd.bias = summary(x1-x0)[4]/sd(c(x1,x0),na.rm=T) # need to correct
  ks = ks.test(x1,x0)$p.value
  ttest = t.test(x1,x0)$p.value
  wilcox = wilcox.test(x1,x0)$p.value
  
  return(c(mean(x1,na.rm=T),mean(x0,na.rm=T),ttest,wilcox,ks))
}

dw = dd[abs(dd$wobl_margin)<=0.1,]

### Figure 1 - no key covariates, looks balance
tab=round(t(mapply(f.stat,as.list(dw[dw$treat==1,c(cov.all,cov.occ)]),
                   as.list(dw[dw$treat==0,c(cov.all,cov.occ)]))),dig=3)
tab = cbind(rownames(tab),tab)

pdf("fig_loess_mps1.pdf",width=8,height=4.5)
plot.pval(tab,legend=TRUE)
dev.off()

###########################################################################################
# Covariate continuous at cut-point plot
###########################################################################################
library(rdd)
library(locfit)
d0 = dd
d = d0[d0$labour==0,]

### Previous vote share margin
bandwidth = IKbandwidth(d$wobl_margin,d$pre_wobl,cutpoint=0)
#CCT.bandwidth=rdbwselect(y=d$DifDPPrv,x=d$DifDPct,c=0)

n.bins=40
breaks0 = seq(-0.1,.1,by=0.01)
bins0 = cut(d$wobl_margin,breaks=breaks0)
bin.prv = tapply(d$pre_wobl,bins0,function(x){return(mean(x,na.rm=TRUE))})
bins = tapply(d$wobl_margin,bins0,function(x){return(mean(x,na.rm=TRUE))})

fit0 <- locfit(pre_wobl~wobl_margin,data=d[d$wobl_margin<0,],alpha=c(0.1^5,bandwidth),deg=1,kern="tcub")
fit1 <- locfit(pre_wobl~wobl_margin,data=d[d$wobl_margin>0,],alpha=c(0.1^5,bandwidth),deg=1,kern="tcub")

loess0.fit = predict(fit0,newdata=d[d$wobl_margin<0,],se.fit=TRUE)
loess1.fit = predict(fit1,newdata=d[d$wobl_margin>0,],se.fit=TRUE)

loess0=loess0.fit$fit
loess1=loess1.fit$fit

loess0.se=loess0.fit$se.fit
loess1.se=loess1.fit$se.fit

pdf("fig_loess_mps2_nonlabour.pdf",width=4.5,height=4.5)
par(cex=0.7)
plot(bins,bin.prv,type="p",las=1,pch=16,col=ifelse(bins<=0,"red3","green4"),
     ylim = c(min(d$DifDPct,min(loess0-loess0.se)),max(d$DifDPct,max(loess1+loess1.se))),
     xlab="Vote margin",ylab="Vote margin in previous election",cex.lab=1.2)
xx = d$wobl_margin[d$wobl_margin>0]
lines(xx[order(xx)],loess1[order(xx)],type="l",xlim=c(-1,1),col="green4",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(loess1[order(xx)]+loess1.se[order(xx)],rev(loess1[order(xx)]-loess1.se[order(xx)])),
        col=rgb(0, 0.5, 0,0.2),border=NA)

xx = d$wobl_margin[d$wobl_margin<0]
lines(xx[order(xx)],loess0[order(xx)],type="l",xlim=c(-1,1),col="red3",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(loess0[order(xx)]+loess0.se[order(xx)],rev(loess0[order(xx)]-loess0.se[order(xx)])),
        col=rgb(0.5, 0, 0,0.2),border=NA)
abline(v=0,col="blue3",lwd=2, lty=2)
dev.off()


########################################
d = dd
#d$prev_attempts = (d$prev_attempts>0)*1

### Previous attempts
bandwidth = IKbandwidth(d$wobl_margin,d$prev_attempts,cutpoint=0)
#CCT.bandwidth=rdbwselect(y=d$DifDPPrv,x=d$DifDPct,c=0)

n.bins=40
breaks0 = seq(-0.1,0.1,by=0.01)
bins0 = cut(d$wobl_margin,breaks=breaks0)
bin.prv = tapply(d$prev_attempts,bins0,function(x){return(mean(x,na.rm=TRUE))})
bins = tapply(d$wobl_margin,bins0,function(x){return(mean(x,na.rm=TRUE))})

fit0 <- locfit(prev_attempts~wobl_margin,data=d[d$wobl_margin<0,],alpha=c(0.1^5,bandwidth),
               deg=1,kern="tcub",family="binomial")
fit1 <- locfit(prev_attempts~wobl_margin,data=d[d$wobl_margin>0,],alpha=c(0.1^5,bandwidth),
               deg=1,kern="tcub",family="binomial")

loess0.fit = predict(fit0,newdata=d[d$wobl_margin<0,],se.fit=TRUE,family="binomial")
loess1.fit = predict(fit1,newdata=d[d$wobl_margin>0,],se.fit=TRUE,family="binomial")

loess0=loess0.fit$fit
loess1=loess1.fit$fit

loess0.se=loess0.fit$se.fit
loess1.se=loess1.fit$se.fit

pdf("fig_loess_mps3.pdf",width=4.5,height=4.5)
par(cex=0.7)
plot(bins,bin.prv,type="p",las=1,pch=16,col=ifelse(bins<=0,"red3","green4"),
     ylim = c(0,1.5),
     xlab="Vote margin",ylab="Previous attempts",cex.lab=1.2)
xx = d$wobl_margin[d$wobl_margin>0]
lines(xx[order(xx)],loess1[order(xx)],type="l",xlim=c(-1,1),col="green4",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(loess1[order(xx)]+loess1.se[order(xx)],rev(loess1[order(xx)]-loess1.se[order(xx)])),
        col=rgb(0, 0.5, 0,0.2),border=NA)

xx = d$wobl_margin[d$wobl_margin<0]
lines(xx[order(xx)],loess0[order(xx)],type="l",xlim=c(-1,1),col="red3",lwd=2)

polygon(c(xx[order(xx)],rev(xx[order(xx)])),
        c(loess0[order(xx)]+loess0.se[order(xx)],rev(loess0[order(xx)]-loess0.se[order(xx)])),
        col=rgb(0.5, 0, 0,0.2),border=NA)
abline(v=0,col="blue3",lwd=2, lty=2)
dev.off()









