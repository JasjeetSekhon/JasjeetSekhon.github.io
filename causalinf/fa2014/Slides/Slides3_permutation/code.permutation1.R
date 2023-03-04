
####################################
# Code section 3
####################################

###############################
# Example difference in means
###############################
setwd("")

rm(list=ls(all=TRUE))

set.seed(13)
x1 = rexp(1000,rate=0.6)
x2 = rexp(1000,rate=0.5)

# Observed statistic
statistic.obs = mean(x2)-mean(x1)
mean(x2)-mean(x1)

x = c(x1,x2)
t = c(rep(1,length(x2)),rep(0,length(x1)))

f.permute = function(){
  id = sample(c(1:length(x)),length(x2))
  t0= rep(0,length(x))
  t0[id]=1
  statistic0 = mean(x[t0==1])-mean(x[t0==0])
  return(statistic0)
}

stat.permutation = replicate(10000,f.permute())

pdf(file="figure1.pdf",width=5.5,height=4)
par(cex=0.8)
hist(stat.permutation,breaks=seq(-0.5,0.5,length=30),col="lightblue",las=1,
     main="Permutation distribution",xlab="Value of test statistic")
abline(v=statistic.obs,col="red4",lwd=2)
text(0.25,1500,"Observed value",col="red4",cex=1.2,lwd=2)
dev.off()


pdf(file="figure2.pdf",width=5.5,height=4)
par(cex=0.8)
hist(stat.permutation,breaks=seq(-0.5,0.5,length=30),col="lightblue",las=1,
     main="Permutation distribution",xlab="Value of test statistic")
abline(v=statistic.obs,col="red4",lwd=2)
text(0.25,1500,"Observed value",col="red4",cex=1.2,lwd=2)

abline(v=-statistic.obs,col="green4",lwd=2)
text(-0.25,1500,"Hypothetical \n observed value",col="green4",cex=1.2,lwd=2)
dev.off()


pdf(file="figure3.pdf",width=5,height=3)
par(cex=0.8)
hist(abs(stat.permutation),breaks=seq(-0.01,0.5,length=50),col="lightblue",las=1,
     main="Permutation distribution",xlab="Value of test statistic")
abline(v=statistic.obs,col="red4",lwd=2)
text(0.3,1000,"Observed value",col="red4",cex=1.2,lwd=2)
dev.off()

######################################################
# Wilcoxon rank sum test
q = rank(x)

f.permute = function(){
  id = sample(c(1:length(x)),length(x2))
  t0= rep(0,length(x))
  t0[id]=1
  statistic0 = sum(q[t0==1])
  return(statistic0)
}
stat.permutation = replicate(10000,f.permute())
statistic.obs = sum(q[t==1])

pdf(file="figure4.pdf",width=5.5,height=4)
par(cex=0.8)
hist(stat.permutation,col="lightblue",las=1,
     breaks=seq(min(min(stat.permutation),statistic.obs)-10,max(max(stat.permutation),statistic.obs)+10,length=200),
     main="Permutation distribution",xlab="Value of test statistic")
abline(v=statistic.obs,col="red4",lwd=2)
text(955000,180,"Observed value",col="red4",cex=1.2,lwd=2)
dev.off()

### two sided P-value
min((sum(stat.permutation<=statistic.obs)+1)/(length(stat.permutation)+1),
    (sum(stat.permutation>=statistic.obs)+1)/(length(stat.permutation)+1))


############################
# KS examples
############################
setwd("")

set.seed(16)
x=rnorm(50,mean=2,sd=1)
y=rnorm(100,mean=2,sd=2)

pdf("figure0_ks.pdf",width=4.5,height=4)
par(cex=0.7)
boxplot(x,y,col=c("red2","blue2"),names=c("X","Y"),las=1)
dev.off()

t.test(x,y)
wilcox.test(x,y)
ks.test(x,y,exact=TRUE)
ks.test(x,y)$stat

################
### plot:
################

### The emperical distribution:
Fx = function(x0,x){
  return(sum(x<=x0)/length(x))
}

Fy = function(y0,y){
  return(sum(y<=y0)/length(y))
}

A = seq(min(c(x,y)),max(c(x,y)),length=1000)
Fx1 = apply(matrix(A,ncol=1),1,Fx,x=x)
Fy1 = apply(matrix(A,ncol=1),1,Fy,y=y)

# adjusting the emperical CDF:

pdf("figure1_ks.pdf",width=5,height=4.5)
par(cex=0.7,cex.lab=1.3)
plot(A,Fx1,col="red2",xlab="w",ylab="F(w)",las=1)
points(A,Fy1,col="blue2")

## adding the KS statistic:
ind = which(abs(Fx1-Fy1)==max(abs(Fx1-Fy1)))
statistic.obs = round(abs(Fx1-Fy1)[ind][1],dig=3)
# add arrow:
arrows(x0=A[ind][1],y0=Fx1[ind][1],x1=A[ind][1],Fy1[ind][1],lwd=3,code=3)
text(-1,0.9,paste("KS statistic: \n D =  ",statistic.obs),cex=1.3,lwd=2)
legend("bottomright",col=c("blue2","red2"),lty=rep(1,2),legend=c("Y CDF","X CDF"),lwd=2,
       cex=1.3)
dev.off()

### Permutation distribution:

treat = c(rep(1,length(x)),rep(0,length(y)))
outcome = c(x,y)

f.ks = function(){
  treat0 = sample(treat,length(treat),replace=FALSE)
  
  x0 = outcome[treat0==1]
  y0 = outcome[treat0==0]
  
  A0 = seq(min(c(x0,y0)),max(c(x0,y0)),length=1000)
  Fx11 = apply(matrix(A0,ncol=1),1,Fx,x=x0)
  Fy11 = apply(matrix(A0,ncol=1),1,Fy,y=y0)
  
  ind0 = which(abs(Fx11-Fy11)==max(abs(Fx11-Fy11)))
  statistic = round(abs(Fx11-Fy11)[ind0][1],dig=4)
  return(statistic)
}

ks.permutation1 = replicate(1000,f.ks())


pdf("figure2_ks.pdf",width=5,height=4.5)
hist(ks.permutation1,breaks=seq(0,0.4,length=30),col="lightblue",
     main="KS test - Binary variables")
abline(v=statistic.obs,col="red4",lwd=2)
pv1 = sum(statistic.obs<ks.permutation1)/length(ks.permutation1)
text(0.28,170,paste("One sided P-value: ",pv1),lwd=2,cex=1)
dev.off()

pdf("figure2b_ks.pdf",width=5,height=4.5)
hist(ks.permutation1,breaks=seq(0,0.4,length=30),col="lightblue",
     main="KS test - Binary variables")
abline(v=statistic.obs,col="red4",lwd=2)
pv1 = sum(statistic.obs<ks.permutation1)/length(ks.permutation1)
text(0.28,170,paste("One sided P-value: ",pv1),lwd=2,cex=1)
abline(v=0.02,col="green4",lwd=2)
dev.off()


######################################################
# Example with binary variables
######################################################

set.seed(14)
x=rbinom(50,size=1,prob=0.5)
y=rbinom(100,size=1,prob=0.7)

t.test(x,y)
wilcox.test(x,y)
ks.test(x,y)

### Emperical CDF:

A = seq(-1,1.2,length=1000)
Fx1 = apply(matrix(A,ncol=1),1,Fx,x=x)
Fy1 = apply(matrix(A,ncol=1),1,Fy,y=y)

# adjusting the emperical CDF:
pdf("figure3_ks.pdf",width=5,height=4.5)
par(cex=0.7,cex.lab=1.3)
plot(A,Fx1,col="red2",xlab="t",ylab="F(t)",type="l",lwd=2,las=1)
lines(A,Fy1,col="blue2",lwd=2)

## adding the KS statistic:
ind = which(abs(Fx1-Fy1)==max(abs(Fx1-Fy1)))
statistic.obs = round(abs(Fx1-Fy1)[ind][1],dig=3)
# add arrow:
arrows(x0=A[ind][50],y0=Fx1[ind][50],x1=A[ind][50],Fy1[ind][50],lwd=2,code=3)
text(-0.5,0.9,paste("KS statistic: ",statistic.obs),cex=1.5,lwd=2)
legend("bottomright",col=c("blue2","red2"),lty=rep(1,2),legend=c("Y CDF","X CDF"),lwd=2,
       cex=1.3)
dev.off()

ks.permutation2 = replicate(1000,f.ks())

pdf("figure4_ks.pdf",width=5,height=4.5)
hist(ks.permutation2,breaks=seq(0,0.4,length=40),col="lightblue",
     main="KS test - Binary variables")
abline(v=statistic.obs,col="red4",lwd=2)
pv2 = sum(statistic.obs<ks.permutation2)/length(ks.permutation2)
text(0.3,100,paste("One sided P-value: ",pv2),lwd=2,cex=1)
dev.off()











