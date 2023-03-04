
######################################################
# Wilcoxon sign rank test

rm(list=ls(all=TRUE))
set.seed(12345)
setwd("")

#install.packages("PairedData")
library(ggplot2)
library(PairedData)
data(ChickWeight)
attach(ChickWeight)
head(ChickWeight)

sum_Z_is_c_si <- (OpenRange>Confinement)*1 #calculating c_si:
d_si <- rank(abs(OpenRange-Confinement)) #calculating d_si:
statistic <- sum(d_si*sum_Z_is_c_si)

wilcox.test(OpenRange,Confinement,paired=TRUE)

### Permutation diftribution under the null
L = 10000
Y = ChickWeight[,c(2,3)]

stat0 <- rep(999,L)
for (i in c(1:L)){
  OpenRange0 <- rep(999,10)
  Confinement0 <- rep(999,10)
  
  for(j in c(1:10)){
    id0 <- sample(c(2,3),1)
    OpenRange0[j] <- Y[j,c(2,3) %in% id0]
    Confinement0[j] <- Y[j,!c(2,3) %in% id0]
  }
  
  sum_Z_is_c_si0 <- (OpenRange0>Confinement0)*1 #calculating c_si:
  d_si0 <- rank(abs(OpenRange0-Confinement0)) #calculating d_si:
  stat0[i] <- sum(d_si0*sum_Z_is_c_si0)
}

pdf(file="fig_sign1.pdf",width=4.5,height=4.5)
ggplot(data=data.frame(stat0),aes(x=stat0))+geom_histogram(fill="blue")+
  theme_bw()+geom_vline(xintercept = statistic,colour="red",size=1)
dev.off()

### P-value
min(sum(statistic<=stat0)/L,sum(statistic>=stat0)/L)*2



###############################
# Confidence sets
###############################


rm(list=ls(all=TRUE))

t=c (12,12,12.9,13.6,16.6,17.2,17.5,18.2,19.1,19.3,19.8,20.3,20.5,20.6,21.3,21.6,22.1)
c=c(5,5.4,6.1,10.9,11.8,12,12.3,14.8,15,16.8,17.2,17.2,17.4,17.5,18.5,18.7,18.7,19.2)


### calculate a one-sided confidence interval:
L = 500
tau.gride = seq(-10,30,length=L)
pv.gride = rep(999,length(tau.gride))

for (j in c(1:length(tau.gride))){
  #if (j %% 5==0){cat("Iteration: ",j,"\n")}
  pv.gride[j] = wilcox.test(t-tau.gride[j],c,exact=FALSE,alternative="greater")$p.value
}  

pdf("fig_ci1.pdf",width=4.5,height=4.5)
par(cex=0.7)
plot(tau.gride,pv.gride,las=1,pch=16,col="blue2",
     xlab="Treatment effect",ylab="P-value")
abline(h=0.05,lty=1,col="red2",lwd=2)

abline(v=tau.gride[tau.gride==min(tau.gride[pv.gride>0.05])],lty=3,col="black",lwd=2)

arrows(x0=tau.gride[tau.gride==min(tau.gride[pv.gride>0.05])]+5,y0=0.4,
       x1=tau.gride[tau.gride==min(tau.gride[pv.gride>0.05])]+10,y1=0.4,code=2,
       length=.2,lwd=2)

tau.min = round(tau.gride[tau.gride==min(tau.gride[pv.gride>0.05])],dig=2)
dev.off()

tau.min

t.test(t,c,alternative="greater")


### Two sided confidence set:

### calculate a one-sided confidence interval:
L = 2000
tau.gride = seq(-5,15,length=L)
pv.gride = rep(999,length(tau.gride))

for (j in c(1:length(tau.gride))){
  #if (j %% 5==0){cat("Iteration: ",j,"\n")}
  pv.gride[j] = wilcox.test(t-tau.gride[j],c,exact=FALSE,alternative="two.sided")$p.value
}  

pdf("fig_ci2.pdf",width=4.5,height=4.5)
par(cex=0.7)
plot(tau.gride,pv.gride,las=1,pch=16,col="blue2",
     xlab="Treatment effect",ylab="P-value")
abline(h=0.05,lty=1,col="red2",lwd=2)

abline(v=tau.gride[tau.gride==min(tau.gride[pv.gride>0.05])],lty=3,col="black",lwd=2)

abline(v=tau.gride[tau.gride==max(tau.gride[pv.gride>0.05])],lty=3,col="black",lwd=2)
dev.off()


