##################################
# P-score and Matching example
##################################
rm(list=ls(all=TRUE))

library(xtable)
library(texreg)
library(foreign)
library(weights)

setwd("~/Dropbox/causalinf.private/fa2014/Slides/Slides_matching2")
data = read.dta("nswre74.dta")
d = data

cov=c("age","ed","black","hisp","married","nodeg","re74","re75")
x = d[,c("age","ed","black","hisp","married","nodeg","re74","re75")]
treat = d$treat

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

# balance plot:
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

tab = cbind(rownames(tab),tab)
par(mfrow=c(1,1),cex=0.7)
plot.pval(tab,legend=TRUE)

#############################################################################
# Example of matching without replacment using Mahalanobis distance matrix
#############################################################################
rm(list=ls(all=TRUE))
source("smahal.R")

### data without the experimental controls with CPS 3
d1 = read.dta("cps1re74.dta")

### Balance tables: CPS-1 vs NSW treatment group ###
x1 = d1[,c("age","ed","black","hisp","married","nodeg","re74","re75")]
treat1 = d1$treat

# Calculating the Mahalanobis distance matrix using the function "smahal.R":
distance.matrix = smahal(treat1,x1)

library("optmatch")
pair.match <-pairmatch(distance.matrix, controls=1,data=d1)

# The matched data set:
dm = d1[is.na(pair.match)=="FALSE",]

x1.m = dm[,c("age","ed","black","hisp","married","nodeg","re74","re75")]
treat1.m = dm$treat

f.stat = function(x1,x0){
  #sd.bias = summary(x1-x0)[4]/sd(c(x1,x0),na.rm=T) # need to correct
  ks = ks.test(x1,x0)$p.value
  ttest = t.test(x1,x0)$p.value
  wilcox = wilcox.test(x1,x0)$p.value
  
  return(c(mean(x1,na.rm=T),mean(x0,na.rm=T),ttest,wilcox,ks))
}

tab1=round(t(mapply(f.stat,as.list(x1.m[treat1.m==1,]),as.list(x1.m[treat1.m==0,]))),dig=3)
colnames(tab1) = c("Ave. Treat","Ave. control","T-test","Wilcoxon","KS")
xtable(tab1,captioin="Pre-matching: balance table",dig=3)

### Using Jas's "Match" function:
match = Match(Tr=treat1,X=x1,replace=FALSE,Weight=2)
dm2 = d1[c(match$index.treated,match$index.control),]
x2.m = dm2[,c("age","ed","black","hisp","married","nodeg","re74","re75")]
treat2.m = dm2$treat

tab2=round(t(mapply(f.stat,as.list(x2.m[treat2.m==1,]),as.list(x2.m[treat2.m==0,]))),dig=3)
colnames(tab2) = c("Ave. Treat","Ave. control","T-test","Wilcoxon","KS")
xtable(tab2,captioin="Pre-matching: balance table",dig=3)

tab.compare = cbind(tab2[,2:3],tab1[,2:3])
xtable(tab.compare,captioin="Balance table: Comparing both matching procedures",dig=3)

################################################################
rm(list=ls(all=TRUE))

### data without the experimental controls with CPS 3
d1 = read.dta("cps1re74.dta")

f.stat = function(x1,x0){
  #sd.bias = summary(x1-x0)[4]/sd(c(x1,x0),na.rm=T) # need to correct
  ks = ks.test(x1,x0)$p.value
  ttest = t.test(x1,x0)$p.value
  wilcox = wilcox.test(x1,x0)$p.value
  
  return(c(mean(x1,na.rm=T),mean(x0,na.rm=T),ttest,wilcox,ks))
}

### Balance tables: CPS-1 vs NSW treatment group ###
x1 = d1[,c("age","ed","black","hisp","married","nodeg","re74","re75")]
treat1 = d1$treat
tab=round(t(mapply(f.stat,as.list(x1[treat1==1,]),as.list(x1[treat1==0,]))),dig=3)
colnames(tab) = c("Ave. Treat","Ave. control","T-test","Wilcoxon","KS")
xtable(tab,captioin="CPS-1: balance table",dig=3)

### Balance tables: CPS-1 vs NSW control group ###
x0.nsw = d[d$treat==0,c("age","ed","black","hisp","married","nodeg","re74","re75")]
x0.cps = d1[d1$treat==0,c("age","ed","black","hisp","married","nodeg","re74","re75")] 

tab=round(t(mapply(f.stat,as.list(x0.nsw),as.list(x0.cps))),dig=3)
colnames(tab) = c("Ave. NSW","Ave. CPS","T-test","Wilcoxon","KS")
xtable(tab,caption="",dig=3)

### Regression adjustment:
lm1 = lm(re78~(.),data=d1)
summary(lm1)

###########################
### Estimating the P-score
###########################
x1 = d1[,c("age","ed","black","hisp","married","nodeg","re74","re75")]
ps.model1 <- glm(treat~(.),data=data.frame(treat=d1$treat,x1),family=binomial(link=logit)); summary(ps.model1)
#ps.model2 <- lm(treat~(.)^2,data=data.frame(treat=d3$treat,x3))
#anova(ps.model2,ps.model1)
ps <- ps.model1$fit

pdf("fig1_match.pdf",width=5,height=4)
par(mfrow=c(1,1),cex=0.7)
plot(density(ps[d1$treat==0]),col=rgb(red= 0, green = 0, blue = 150, maxColorValue = 255),
     las=1,main="Distribution of the P-score \n Treatment vs Control",ylim=c(0,6),
     xlab="P-score",lwd=2,lty="solid",bty="l")
polygon(density(ps[d1$treat==0]), col = rgb(red = 0, green = 0, blue = 150,
                                            alpha = 100, maxColorValue = 255), lty = "solid")
lines(density(ps[d1$treat==1]),lwd=2,
      col = rgb(red = 150, green = 0, blue = 0,alpha = 200, maxColorValue = 255))
polygon(density(ps[d1$treat==1]), col = rgb(red = 150, green = 0, blue = 0,
                                            alpha = 100, maxColorValue = 255))

legend("topright",lty=rep(1,2),col=c(rgb(red= 150, green = 0, blue = 0, maxColorValue = 255),
                                     rgb(red= 0, green = 0, blue = 150, maxColorValue = 255)),
       legend=c("Treatment","Control"),lwd=2,bty="n")
dev.off()


pdf("fig2_match.pdf",width=5,height=4)
par(mfrow=c(1,1),cex=0.7)
boxplot(ps[d1$treat==1],ps[d1$treat==0],col=c("red4","blue3"),names=c("Treatment","Control"),las=1,
        ylab="P-score",
        main="Distribution of the P-score \n NSW treatment vs CPS-1 control")
dev.off()

###########################
### IPW
###########################

# defining the weights:
w = rep(999,dim(d1)[1])
w[d1$treat==1] = 1/ps[d1$treat==1]
w[d1$treat==0] = 1/(1-ps[d1$treat==0])
summary(w)

### Regression adjustment CPS-1, using IPSW: 
lm2 = lm(re78~(.),data=d1,weights=w);summary(lm2)

### Experimental benchmark:
lm0 = lm(re78~(.),data=d); summary(lm0)

### Regression adjustmen with CPS-1t:
lm1 = lm(re78~(.),data=d1);summary(lm1)

screenreg(list(lm0,lm1,lm2))

texreg(list(lm0,lm1,lm2),
       custom.model.names = c("Regression NSW","Regression CPS-1", "IPSW"))


### Checking the balance after weighting by the P-score:

f.stat.weights = function(x1,x0,t0,weights){
  #sd.bias = summary(x1-x0)[4]/sd(c(x1,x0),na.rm=T) # need to correct
  ttest0 = wtd.t.test(x1,x0,weight=weights[t0==1],weighty=weights[t0==0])
  ttest.pv = ttest0$coefficients[3]
  mean.x1 = ttest0$additional[2]
  mean.x0 = ttest0$additional[3]
  return(c(mean.x1,mean.x0,ttest.pv))
}

#### Overlap condition:
ind = ps>=min(ps[d1$treat==1]) & ps<=max(ps[d1$treat==1])
#ind = !is.na(ps)

# Weighted
tab1=round(t(mapply(f.stat.weights,as.list(x1[treat1==1 & ind,]),as.list(x1[treat1==0 & ind,]),
                    MoreArgs=list(t0=treat1[ind],weights=w[ind]))),dig=3)
colnames(tab1) = c("Ave. Treat","Ave. control","T-test")
xtable(tab1,captioin="After IPSW: balance table",dig=3)

# Not weighted
tab=round(t(mapply(f.stat,as.list(x1[treat1==1  & ind,]),as.list(x1[treat1==0 & ind,]))),dig=3)
colnames(tab) = c("Ave. Treat","Ave. control","T-test","Wilcoxon","KS")

### Regression check

anova(lm(as.formula(paste("treat~",paste(cov,collapse="+"))),data=d1,weights=w),
      lm(as.formula(paste("treat~","1")),data=d1,weights=w))

anova(lm(as.formula(paste("treat~",paste(cov,collapse="+"))),data=d1),
      lm(as.formula(paste("treat~","1")),data=d1))


##########################################
### Smith and Todd 2005 diff-in-diff
##########################################
rm(list=ls(all=TRUE))
d = read.dta("nswre74.dta")
cov=c("age","ed","black","hisp","married","nodeg","re74","re75")
x = d[,c("age","ed","black","hisp","married","nodeg","re74","re75")]
treat = d$treat


# diff-in-diff
with(d[d$treat==1,],mean(re78-re74)) - with(d[d$treat==0,],mean(re78-re74))
with(d[d$treat==1,],mean(re78-re75)) - with(d[d$treat==0,],mean(re78-re75))

# Regression adjustments
d$diff74 = with(d,re78-re74)
d$diff75 = with(d,re78-re75)

lm01 = lm(diff74~(.),data=d[,c("diff74","treat")]) #DinD
lm02 = lm(diff75~(.),data=d[,c("diff75","treat")]) #DinD

lm1 = lm(re78~(.),data=d[,c("re78","treat",cov)])

lm2 = lm(diff74~(.),data=d[,c("diff74","treat",cov)]) 
lm3 = lm(diff75~(.),data=d[,c("diff75","treat",cov)])


screenreg(list(lm01,lm02,lm1,lm2,lm3))

texreg(list(lm01,lm02,lm1,lm2,lm3),
       custom.model.names = c("DinD re74","DinD re75","Regression NSW",
                              "DinD re74", "DinD re75"))


###########################
### P-score 
###########################
d3 = read.dta("cps3re74.dta")

### Balance tables: CPS-1 vs NSW treatment group ###
x3 = d3[,c("age","ed","black","hisp","married","nodeg","re74","re75")]
treat3 = d3$treat

ps.model3 <- glm(treat~(.),data=data.frame(treat=d3$treat,x3),family=binomial(link=logit)); summary(ps.model3)
#ps.model2 <- lm(treat~(.)^2,data=data.frame(treat=d3$treat,x3))
#anova(ps.model2,ps.model1)
ps <- ps.model3$fit

# defining the weights:
w = rep(999,dim(d3)[1])
w[d3$treat==1] = 1/ps[d3$treat==1]
w[d3$treat==0] = 1/(1-ps[d3$treat==0])
summary(w)

### limiting to observations with overlap in P-score of treated units
ind = ps>=min(ps[d3$treat==1]) & ps<=max(ps[d3$treat==1])

w.caliper = rep(999,dim(d3)[1])
w.caliper[ind & d3$treat==1] = ps[ind & d3$treat==1]
w.caliper[ind & d3$treat==0] = 1-ps[ind & d3$treat==0]
w.caliper <- w.caliper[w.caliper!=999]
summary(1/w.caliper)

lm4 = lm(re78~(.),data=d3[ind,],weights=w.caliper);summary(lm4)
lm4.check = lm(re78~(.),data=d3[ind,]);summary(lm4.check)

### Simple weighted difference in means, and not weighted difference in means
sum(d3$re78*d3$treat)/sum(d3$treat/ps) - sum(d3$re78*(1-d3$treat))/sum((1-d3$treat)/(1-ps))
sum(d3$re78*d3$treat)/sum(d3$treat) - sum(d3$re78*(1-d3$treat))/sum(1-d3$treat)

wtd.t.test(d3$re78[d3$treat==1],d3$re78[d3$treat==0],weight = ps[d3$treat==1],weighty = 1-ps[d3$treat==0])
wtd.t.test(d3$re78[d3$treat==1],d3$re78[d3$treat==0])

#######################################
# Matching with replacment
#######################################

ps.t = ps[d3$treat==1]
ps.c = ps[d3$treat==0]
nt <- length(ps.t)
d2 = function(x1,x2){return((x1-x2)^2)}
ps.c.match <-index.control <- as.list(rep(999,nt))
for (i in c(1:nt)){
  ps.c.match[i] = ps.c[which.min(d2(ps.t[i],ps.c))]
  index.control[i] = which.min(d2(ps.t[i],ps.c))
}
cat("Check for ties",length(index.control)==length(unlist(index.control)),"\n")
cat("If TRUE there are no ties","\n")
index.control = unlist(index.control)
ps.c.match = unlist(ps.c.match)

ks.test(ps.t,ps.c.match)
boxplot(ps.t,ps.c)

# data after matching:
dm3 = rbind(d3[d3$treat==0,][index.control,],d3[d3$treat==1,])






#############################################
# Which observations in the control is used?
#############################################

nc.match = length(unique(ps.c.match))
index <- c(rep(0,length(d3$treat[d3$treat==0])-nc.match),index.control)

#pdf("fig3_ps.pdf",width=5.5,height=4)
par(mfrow=c(1,2),cex=0.7)
hist(index,breaks=length(index),las=1,xlab="Observation index \n Control",
     col="lightblue",main="All observations")
hist(index[index>0],breaks=length(index[index>0]),las=1,xlab="Observation index \n Control",
     col="lightblue",main="Observations in matched data")
#dev.off()

summary(as.numeric(table(index.control)))

### Regression adjustment:
lm1.m = lm(re78~(.),data=dm3); summary(lm1.m)

sum(dm3$re78*dm3$treat)/sum(dm3$treat) - sum(dm3$re78*(1-dm3$treat))/sum(1-dm3$treat)




###############################################################
# Checking equivalence of weighting and matching
###############################################################

### Code not good!!

#weights0 = rep(999,dim(d3[d3$treat==0,])[1])
#d3$weights0[d3$treat==1]=1
#weights0[c(1:length(weights0)) %in% index[index>0]] = index


###############################################################
# CPS - 1
###############################################################

### data without the experimental controls with CPS 3
d1 = read.dta("cps1re74.dta")





###########################
### Estimating the P-score
###########################
x1 = d1[,c("age","ed","black","hisp","married","nodeg","re74","re75")]
ps.model1 <- glm(treat~(.),data=data.frame(treat=d1$treat,x1),family=binomial(link=logit)); summary(ps.model1)
#ps.model2 <- lm(treat~(.)^2,data=data.frame(treat=d3$treat,x3))
#anova(ps.model2,ps.model1)
ps <- ps.model1$fit


### P-score balance

par(mfrow=c(1,1),cex=0.7)
boxplot(ps[d1$treat==1],ps[d1$treat==0],col=c("red4","blue3"),names=c("Treatment","Control"),las=1,
        ylab="P-score",
        main="Distribution of the P-score \n Treatment vs Control")

### looking only on the overlap:
ind = min(ps[d1$treat==1])<ps

boxplot(ps[d1$treat==1 & ind],ps[d1$treat==0 & ind],col=c("red4","blue3"),names=c("Treatment","Control"),las=1,
        ylab="P-score",
        main="Distribution of the P-score \n Treatment vs Control")

summary(lm(re78~(.),data=d1[ind,]))

###########################
### P-score and IPW
###########################
summary(1/ps)
lm2 = lm(re78~(.),data=d1,weights=1/ps);summary(lm2)

### limiting to observations with overlap in P-score of treated units
ind = ps>min(ps[d1$treat==1]) & ps<max(ps[d1$treat==1])
lm4 = lm(re78~(.),data=d1[ind,],weights=1/ps[ind]);summary(lm4)

summary(1/ps[ind])

### limiting to observations with overlap in P-score
ind = ps>min(ps[d1$treat==1]) & ps<max(ps[d1$treat==0])
lm3 = lm(re78~(.),data=d1[ind,],weights=1/ps[ind]);summary(lm3)

summary(1/ps[ind])









### Balance test:
#source("/media/yotam/Windows8_OS/Users/Yotam/Dropbox/Class.test/R.code/class.test.functions_new.R")

### Observed test statistic (RMSE)
#stat.obs0 = lm(treat~(.)^2,data=data.frame(x,treat))
#stat.obs1.nomatch = rmse(stat.obs0$fitt,treat)

### Permutation distribution
#times = 500
#sim1.nomatch = apply(matrix(c(1:times),ncol=1),1,f.class,x=x,nt=length(treat[treat==1]),method="OLS",
#                     criterion="RMSE")
#par(mfrow=c(1,1),cex=0.7)
#hist(sim1.nomatch,breaks=seq(min(sim1.nomatch),max(sim1.nomatch),length.out=20),col="lightblue",
#     xlim=c(min(sim1.nomatch,stat.obs1.nomatch)-0.01,max(sim1.nomatch)+0.01),
#     xlab="RMSE",main="Permutation distribution of the RMSE \n 1,000 permutations")
#abline(v=stat.obs1.nomatch,col="red4",lwd=2)









