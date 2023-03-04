##################################
# P-score example
##################################
rm(list=ls(all=TRUE))
library(xtable)

setwd("")
load("welder.RData")
d = data

x = d[,c("age","smoker","black")]
treat = d$treat

#######################################
# Balance prior to matching
#######################################

f.stat = function(x1,x0){
  #sd.bias = summary(x1-x0)[4]/sd(c(x1,x0),na.rm=T) # need to correct
  ks = ks.test(x1,x0)$p.value
  ttest = t.test(x1,x0)$p.value
  wilcox = wilcox.test(x1,x0)$p.value
  
  return(c(mean(x1,na.rm=T),mean(x0,na.rm=T),ttest,wilcox,ks))
}

tab=round(t(mapply(f.stat,as.list(x[treat==1,]),as.list(x[treat==0,]))),dig=3)
colnames(tab) = c("Ave. Treat","Ave. control","T-test","Wilcoxon","KS")
xtable(tab,captioin="Pre-matching: balance table",dig=3)

#######################################
# Estimating P-score
#######################################

ps.model1 <- glm(treat~(.),data=data.frame(treat=treat,x),family=binomial(link=logit))
ps.model2 <- glm(treat~(.)^2,data=data.frame(treat=treat,x),family=binomial(link=logit))

ps1 <- ps.model1$fit
ps2 <- ps.model2$fit

pdf("fig1_ps.pdf",width=5,height=4)
par(mfrow=c(1,2),cex=0.7)
plot(density(ps1[treat==0]),col="blue3",las=1,main="No interactions",ylim=c(0,4.8),
     xlab="P-score",lwd=2)
lines(density(ps1[treat==1]),col="red4",lwd=2)
legend("topleft",lty=rep(1,2),col=c("red4","blue3"),legend=c("Treatment","Control"),lwd=2,bty="n")

plot(density(ps2[treat==0]),col="blue3",las=1,main="With interactions",ylim=c(0,4.5),
     xlab="P-score",lwd=2)
lines(density(ps2[treat==1]),col="red4",lwd=2)
legend("topleft",lty=rep(1,2),col=c("red4","blue3"),legend=c("Treatment","Control"),lwd=2,bty="n")
dev.off()

pdf("fig2_ps.pdf",width=5,height=4)
par(mfrow=c(1,2),cex=0.7)
boxplot(ps1[treat==1],ps1[treat==0],col=c("red4","blue3"),names=c("Treatment","Control"),las=1,
        ylab="P-score",main="No interactions")
boxplot(ps2[treat==1],ps2[treat==0],col=c("red4","blue3"),names=c("Treatment","Control"),las=1,
        ylab="P-score",main="With interactions")
dev.off()

par(mfrow=c(1,1),cex=0.7)
plot(density(ps1[treat==0]),col=rgb(red= 0, green = 0, blue = 150, maxColorValue = 255),
     las=1,main="No interactions",ylim=c(0,4.8),xlim=c(-0.1,1),
     xlab="P-score",lwd=2,lty="solid",bty="l")
polygon(density(ps1[treat==0]), col = rgb(red = 0, green = 0, blue = 150,
                                          alpha = 100, maxColorValue = 255), lty = "solid")

lines(density(ps1[treat==1]),lwd=2,
      col = rgb(red = 150, green = 0, blue = 0,alpha = 200, maxColorValue = 255))
polygon(density(ps1[treat==1]), col = rgb(red = 150, green = 0, blue = 0,
                                          alpha = 100, maxColorValue = 255))

#######################################
# Matching with replacment
#######################################

### Model 1

ps1.t = ps1[treat==1]
ps1.c = ps1[treat==0]
nt <- length(ps1.t)
d2 = function(x1,x2){return((x1-x2)^2)}
ps1.c.match <-index.control <- as.list(rep(999,nt))
for (i in c(1:nt)){
  ps1.c.match[i] = ps1.c[which.min(d2(ps1.t[i],ps1.c))]
  index.control[i] = which.min(d2(ps1.t[i],ps1.c))
}
cat("Check for ties",length(index.control)==length(unlist(index.control)),"\n")
cat("If TRUE there are no ties","\n")
index.control1 = unlist(index.control)
ps1.c.match = unlist(ps1.c.match)

ks.test(ps1.t,ps1.c.match)

### Model 2

ps2.t = ps2[treat==1]
ps2.c = ps2[treat==0]

ps2.c.match <-index.control <- as.list(rep(999,nt))
for (i in c(1:nt)){
  ps2.c.match[i] = ps2.c[which.min(d2(ps2.t[i],ps1.c))]
  index.control[i] = which.min(d2(ps2.t[i],ps1.c))
}

cat("Check for ties",length(index.control)==length(unlist(index.control)),"\n")
cat("If TRUE there are no ties","\n")

index.control2 = unlist(index.control)
ps2.c.match = unlist(ps2.c.match)

ks.test(ps2.t,ps2.c.match)
#####################################
# P-score balance after matching
#####################################

pdf("fig11_ps.pdf",width=5,height=4)
par(mfrow=c(1,2),cex=0.7)
plot(density(ps1[treat==0]),col="blue3",las=1,main="No interactions",ylim=c(0,4.8),
     xlab="P-score",lwd=2)
lines(density(ps1[treat==1]),col="red4",lwd=2)
lines(density(ps1.c.match),col="blue3",lty=2,lwd=2)
legend("topleft",lty=c(1,1,2),col=c("red4","blue3","blue3"),
       legend=c("Treatment","Control","Matched Control"),lwd=2,bty="n")

plot(density(ps2[treat==0]),col="blue3",las=1,main="With interactions",ylim=c(0,4.8),
     xlab="P-score",lwd=2)
lines(density(ps2[treat==1]),col="red4",lwd=2)
lines(density(ps2.c.match),col="blue3",lty=2,lwd=2)
legend("topleft",lty=c(1,1,2),col=c("red4","blue3","blue3"),
       legend=c("Treatment","Control","Matched Control"),lwd=2,bty="n")
dev.off()


pdf("fig22_ps.pdf",width=5.5,height=4)
par(mfrow=c(1,2),cex=0.7)
boxplot(ps1[treat==1],ps1.c.match,ps1[treat==0],col=c("red4","lightblue","blue3"),
        names=c("Treatment","\n Match \n control","Control"),las=1, ylab="P-score",main="No interactions")

boxplot(ps2[treat==1],ps2.c.match,ps2[treat==0],col=c("red4","lightblue","blue3"),
        names=c("Treatment","\n Match \n control","Control"),las=1, ylab="P-score",
        main="With interactions")
dev.off()


#############################################
# Which observations in the control is used?
#############################################

nc.match1 = length(unique(ps1.c.match))
index1 <- c(rep(0,nt-nc.match1),index.control1)
nc.match2 = length(unique(ps2.c.match))
index2 <- c(rep(0,nt-nc.match2),index.control2)

pdf("fig3_ps.pdf",width=5.5,height=4)
par(mfrow=c(1,2),cex=0.7)
hist(index1,breaks=length(index1),las=1,xlab="Observation index \n Control",
     ylim=c(0,12),col="lightblue",main="No interactions")
hist(index2,breaks=length(index2),las=1,xlab="Observation index \n Control",
     ylim=c(0,12),col="lightblue",main="With interactions")
dev.off()


####################################
# Matching P-score example - univariate matching
# PSID example 
# The effect of college on wife logwage
####################################

rm(list=ls(all=TRUE))
setwd("")

library(AER)
data(Mroz)
d = Mroz

d$lfp = (d$lfp=="yes")*1
d$hc = (d$hc=="yes")*1
d$wc = (d$wc=="yes")*1

x = d[,c("age","k5","k618","lfp","hc","inc")]
treat = d$wc
y = d$lwg

f.stat = function(x1,x0){
  #sd.bias = summary(x1-x0)[4]/sd(c(x1,x0),na.rm=T) # need to correct
  ks = ks.test(x1,x0)$p.value
  ttest = t.test(x1,x0)$p.value
  wilcox = wilcox.test(x1,x0)$p.value
  
  return(c(mean(x1,na.rm=T),mean(x0,na.rm=T),ttest,wilcox,ks))
}

tab=round(t(mapply(f.stat,as.list(x[treat==1,]),as.list(x[treat==0,]))),dig=3)
colnames(tab) = c("Ave. Treat","Ave. control","T-test","Wilcoxon","KS")
xtable(tab,captioin="Pre-matching: balance table",dig=3)

### Regression estimation:
lm.wc = lm(lwg~(.),data=d)
summary(lm.wc)

lm.wc2 = lm(lwg~(.),data=d[d$lfp==1,])
summary(lm.wc2)

texreg(list(lm.wc,lm.wc2),custom.model.names = c("Full sample","Only working women"),dig=4)

### P-score balance prior to matching:

ps.model1 <- glm(treat~(.),data=data.frame(treat=treat,x),family=binomial(link=logit))
ps.model2 <- glm(treat~(.)^2,data=data.frame(treat=treat,x),family=binomial(link=logit))

ps1 <- ps.model1$fit
ps2 <- ps.model2$fit

pdf("fig1_psid.pdf",width=5,height=4)
par(mfrow=c(1,2),cex=0.7)
boxplot(ps1[treat==1],ps1[treat==0],col=c("red4","blue3"),names=c("Treatment","Control"),las=1,
        ylab="P-score",main="No interactions")
boxplot(ps2[treat==1],ps2[treat==0],col=c("red4","blue3"),names=c("Treatment","Control"),las=1,
        ylab="P-score",main="With interactions")
dev.off()

### Matching - no interactions

ps1.t = ps1[treat==1]
ps1.c = ps1[treat==0]
nt <- length(ps1.t)
d2 = function(x1,x2){return((x1-x2)^2)}
ps1.c.match <-index.control <- as.list(rep(999,nt))
for (i in c(1:nt)){
  ps1.c.match[i] = ps1.c[which.min(d2(ps1.t[i],ps1.c))]
  index.control[i] = which.min(d2(ps1.t[i],ps1.c))
}
cat("Check for ties",length(index.control)==length(unlist(index.control)),"\n")
cat("If TRUE there are no ties","\n")
index.control1 = unlist(index.control)
ps1.c.match = unlist(ps1.c.match)

ks.test(ps1.t,ps1.c.match)

### P-score after matching
pdf("fig2_psid.pdf",width=5.5,height=4)
par(mfrow=c(1,1),cex=0.7)
boxplot(ps1[treat==1],ps1.c.match,ps1[treat==0],col=c("red4","lightblue","blue3"),
        names=c("Treatment","\n Match \n control","Control"),las=1, ylab="P-score",main="No interactions")
dev.off()


### Estimation: Differnce in means:
t.test(y[treat==1],y[index.control1])

dc = rbind(d[index.control1,],d[treat==1,])
lm.wc = lm(lwg~(.),data=dc)
summary(lm.wc)

lm.wc2 = lm(lwg~(.),data=dc[dc$lfp==1,])
summary(lm.wc2)

texreg(list(lm.wc,lm.wc2),custom.model.names = c("Full sample \\ (After matching)",
                                                 "Only working women \\ (After matching)"),dig=4)


### Check the matching worked (multiple replacment)
summary(as.numeric(table(index.control1)))


############################################################
# The relation between husband education and wife education
############################################################

with(d,table(wc,hc))
with(dc,table(wc,hc))

























