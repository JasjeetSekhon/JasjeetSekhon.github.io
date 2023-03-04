rm(list=ls())           
load(file=url("http://sekhon.berkeley.edu/causalinf/data/section5.RData"))

#Whats the naive resul`t?
#positive with a T-stat of 3
summary(lm(PMDB.win.04~treat,data=data))

#Plot the density of the forcing variable
plot(density(data$vote.margin),col="purple", lwd=2, main="Density of the Vote Margin")
abline(v=0,lwd=2)

hist(data$vote.margin,col="purple", lwd=2,breaks=15)
abline(v=0,lwd=2)

####
###Binning
####

#number of bins
bins=100

#create bin midpoints using the cut function
bins.losers <- levels(cut(data$vote.margin[data$treat==0],bins/2))
bin.midpoints.losers <- (as.numeric(as.numeric(sub(pattern="\\((.*),.*",x=bins.losers,replacement="\\1"))) + as.numeric(sub(pattern=".*,(.*).",x=bins.losers,replacement="\\1")))/2
bins.winners <- levels(cut(data$vote.margin[data$treat==1],bins/2))
bin.midpoints.winners <- (as.numeric(as.numeric(sub(pattern="\\((.*),.*",x=bins.winners,replacement="\\1"))) + as.numeric(sub(pattern=".*,(.*).",x=bins.winners,replacement="\\1")))/2


#GENDER
#mean within bins for losers, use the "cut " function
bin.means.losers <- tapply(data$male.loser[data$treat==0],cut(data$vote.margin[data$treat==0],bins/2),mean)
#mean within bins for winners
bin.means.winners <- tapply(data$male.winner[data$treat==1],cut(data$vote.margin[data$treat==1],bins/2),mean)
#Plot
plot(1,type="n",ylim=c(.8,1),xlim=c(-.3,.3),xlab="2000 Vote Margin", ylab="P(Male Candidate)")
points(bin.midpoints.losers,bin.means.losers,pch=19,col="blue")
points(bin.midpoints.winners,bin.means.winners,pch=19,col="red")
abline(v=0,lwd=2)
#Add a loess line
lines(loess.smooth(bin.midpoints.losers,bin.means.losers),col="blue",lwd=2)
lines(loess.smooth(bin.midpoints.winners,bin.means.winners),col="red",lwd=2)

#Population in 2000
#number of bins
bins=100
#mean within bins for losers, use the "cut " function
bin.means.losers <- tapply(log(data$pop.2000[data$treat==0]),cut(data$vote.margin[data$treat==0],bins/2),mean)
#mean within bins for winners
bin.means.winners <- tapply(log(data$pop.2000[data$treat==1]),cut(data$vote.margin[data$treat==1],bins/2),mean)
#Plot
plot(1,type="n",ylim=c(8.9,9.6),xlim=c(-.3,.3),xlab="2000 Vote Margin", ylab="Population")
points(bin.midpoints.losers,bin.means.losers,pch=19,col="blue")
points(bin.midpoints.winners,bin.means.winners,pch=19,col="red")
abline(v=0,lwd=2)
#Add a loess line
lines(loess.smooth(bin.midpoints.losers,bin.means.losers),col="blue",lwd=2)
lines(loess.smooth(bin.midpoints.winners,bin.means.winners),col="red",lwd=2)


#Number of Cattle in 1985
#number of bins
bins=100
#mean within bins for losers, use the "cut " function
bin.means.losers <- tapply(log(data$cattle.1985[data$treat==0]+1),cut(data$vote.margin[data$treat==0],bins/2),mean)
#mean within bins for winners
bin.means.winners <- tapply(log(data$cattle.1985[data$treat==1]+1),cut(data$vote.margin[data$treat==1],bins/2),mean)
#Plot
plot(1,type="n",ylim=c(8.9,10.5),xlim=c(-.3,.3),xlab="2000 Vote Margin", ylab="Log Number of Cattle in 1985")
points(bin.midpoints.losers,bin.means.losers,pch=19,col="blue")
points(bin.midpoints.winners,bin.means.winners,pch=19,col="red")
abline(v=0,lwd=2)
#Add a loess line
lines(loess.smooth(bin.midpoints.losers,bin.means.losers),col="blue",lwd=2)
lines(loess.smooth(bin.midpoints.winners,bin.means.winners),col="red",lwd=2)




#PROBABILITY OF WINNING IN 2004
#number of bins
bins=100
#mean within bins for losers, use the "cut " function
bin.means.losers <- tapply(data$PMDB.win.04[data$treat==0],cut(data$vote.margin[data$treat==0],bins/2),mean)
#mean within bins for winners
bin.means.winners <- tapply(data$PMDB.win.04[data$treat==1],cut(data$vote.margin[data$treat==1],bins/2),mean)
#Plot
plot(1,type="n",ylim=c(0,.6),xlim=c(-.3,.3),xlab="2000 Vote Margin", ylab="Probability of Winning in 2004")
points(bin.midpoints.losers,bin.means.losers,pch=19,col="blue")
points(bin.midpoints.winners,bin.means.winners,pch=19,col="red")
abline(v=0,lwd=2)
#Add a loess line
lines(loess.smooth(bin.midpoints.losers,bin.means.losers),col="blue",lwd=2)
lines(loess.smooth(bin.midpoints.winners,bin.means.winners),col="red",lwd=2)

#Pct Vote share in 2004
#number of bins
#mean within bins for losers, use the "cut " function
bin.means.losers <- tapply(data$PMDB.vote.share.04[data$treat==0],cut(data$vote.margin[data$treat==0],bins/2),mean)
#mean within bins for winners
bin.means.winners <- tapply(data$PMDB.vote.share.04[data$treat==1],cut(data$vote.margin[data$treat==1],bins/2),mean)
#Plot
plot(1,type="n",ylim=c(0,.5),xlim=c(-.5,.5),xlab="2000 Vote Margin", ylab="Vote Share in 2004")
points(bin.midpoints.losers,bin.means.losers,pch=19,col="blue")
points(bin.midpoints.winners,bin.means.winners,pch=19,col="red")
abline(v=0,lwd=2)

#Add a loess line
lines(loess.smooth(bin.midpoints.losers,bin.means.losers),col="blue",lwd=2)
lines(loess.smooth(bin.midpoints.winners,bin.means.winners),col="red",lwd=2)



####Estimation
##Simple rectangular kernel
#whats our bandwidth? 
h <- .03
#trim the data
data.trim <- data[data$vote.margin<h & data$vote.margin>-h,]
#estimate LATE
late.trim <- mean(data.trim$PMDB.win.04[data.trim$treat==1]) - mean(data.trim$PMDB.win.04[data.trim$treat==0])
#estimate SE using Neyman formula
se.trim <- sqrt(var(data.trim$PMDB.win.04[data.trim$treat==1])/length(data.trim$PMDB.win.04[data.trim$treat==1]) +  var(data.trim$PMDB.win.04[data.trim$treat==0])/length(data.trim$PMDB.win.04[data.trim$treat==0]))
late.trim
se.trim

#what about a covariate?
late.trim <- mean(log(data.trim$pop.2000[data.trim$treat==1])) - mean(log(data.trim$pop.2000[data.trim$treat==0]))
se.trim <- sqrt(var(log(data.trim$pop.2000[data.trim$treat==1]))/length(data.trim$gini.2000[data.trim$treat==1]) +  log(var(data.trim$pop.2000[data.trim$treat==0]))/length(data.trim$pop.2000[data.trim$treat==0]))
late.trim
se.trim


##Local Linear Regression
#whats our bandwidth?
h <- .03
#trim the data
data.trim <- data[data$vote.margin<h & data$vote.margin>-h,]
#estimate two regressions
left.reg <- lm(PMDB.win.04~vote.margin,data=data.trim[data.trim$vote.margin<0,])
right.reg <- lm(PMDB.win.04~vote.margin,data=data.trim[data.trim$vote.margin>=0,])
coef(left.reg)[1] - coef(right.reg)[1] 
#For inference, use the combined model
summary(lm(PMDB.win.04~treat + vote.margin + vote.margin*treat, data = data.trim))

h <- seq(from = .01, to = .2, by = .01)



##CROSS VALIDATION
data.window <- data[data$vote.margin>-.5 & data$vote.margin<.5, ]

cv <- rep(NA,times=length(h))
for(j in 1:length(h)){
  y.hat <- rep(NA,times=nrow(data.window))
  for (i in 1:nrow(data.window)){
    c <- data.window$vote.margin[i]
    if(c<0){
      x.bw <- data$vote.margin[data$vote.margin< c & data$vote.margin>(c- h[j])]
      y.bw <- data$PMDB.vote.share.04[data$vote.margin< c & data$vote.margin>(c- h[j])]
    }
    if(c>=0){
      x.bw <- data$vote.margin[data$vote.margin> c & data$vote.margin<(c+ h[j])]
      y.bw <- data$PMDB.vote.share.04[data$vote.margin> c & data$vote.margin<(c+ h[j])]
    }
    x.bw <- x.bw - c
    y.hat[i] <- coef(lm(y.bw~x.bw))[1] + c
  }
    cv[j] <- sum((data.window$PMDB.vote.share.04 - y.hat)^2)/nrow(data.window)
    print(j)
}


plot(h,cv,type="l",lwd=2,col="red")
points(h[which(cv==min(cv))], min(cv),pch = 19, cex=2,col="blue")

h.opt <- h[which(cv==min(cv))]
#trim the data
data.trim <- data[data$vote.margin<h.opt & data$vote.margin>-h.opt,]
#For inference, use the combined model
summary(lm(PMDB.win.04~treat + vote.margin + vote.margin*treat, data = data.trim))
