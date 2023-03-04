##################################
# P-score and Matching example
##################################
rm(list=ls(all=TRUE))

library(xtable)
library(texreg)
library(foreign)
library(weights)
library(Matching)

setwd("") # link to the data in the section notes
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
plot.pval(tab,legend=TRUE)

### data without the experimental controls with CPS 3
d1 = read.dta("cps1re74.dta")

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

### Matching with replacement on P-score
match=Match(Tr=treat1,X=ps,replace=TRUE,ties=FALSE)
index = c(match$index.treated,match$index.control)
dm = d1[index,]
lm3 = lm(re78~(.),data=dm);summary(lm3) 

screenreg(list(lm0,lm1,lm2,lm3))

texreg(list(lm0,lm1,lm2,lm3),
       custom.model.names = c("Regression NSW","Regression CPS-1", "IPSW","Matching on P-score"))

### Checking the balance after weighting by the P-score:

f.stat.weights = function(x1,x0,t0,weights){
  #sd.bias = summary(x1-x0)[4]/sd(c(x1,x0),na.rm=T) # need to correct
  ttest0 = wtd.t.test(x1,x0,weight=weights[t0==1],weighty=weights[t0==0])
  ttest.pv = ttest0$coefficients[3]
  mean.x1 = ttest0$additional[2]
  mean.x0 = ttest0$additional[3]
  return(c(mean.x1,mean.x0,ttest.pv))
}

# Weighted
ind = !is.na(ps)
tab1=round(t(mapply(f.stat.weights,as.list(x1[treat1==1 & ind,]),as.list(x1[treat1==0 & ind,]),
                    MoreArgs=list(t0=treat1[ind],weights=w[ind]))),dig=3)
colnames(tab1) = c("Ave. Treat","Ave. control","T-test")
xtable(tab1,captioin="After IPSW: balance table",dig=3)

### Distribution of the weights
summary(w[treat1==0])
summary(w[treat1==1])

######################################################################
# Limiting the sample - imposing overlap: No p-score \approx 1 or 0
######################################################################

#### Overlap condition:
#ind = ps>=min(ps[d1$treat==1]) & ps<=max(ps[d1$treat==1])
ind = ps>=0.1 & ps<=0.9 # the second limitatin has no meaning as the max P-score is 0.48

# weighted
tab1=round(t(mapply(f.stat.weights,as.list(x1[treat1==1 & ind,]),as.list(x1[treat1==0 & ind,]),
                    MoreArgs=list(t0=treat1[ind],weights=w[ind]))),dig=3)
colnames(tab1) = c("Ave. Treat","Ave. control","T-test")
xtable(tab1,captioin="After IPSW: balance table",dig=3)

### Regression adjustment CPS-1, using IPSW: 
lm5 = lm(re78~(.),data=d1[ind,],weights=w[ind]);summary(lm5)

### Distribution of the weights
summary(w[treat1==0 & ind])
summary(w[treat1==1 & ind])

texreg(list(lm0,lm1,lm5,lm3),
       custom.model.names = c("Regression NSW","Regression CPS-1", "IPSW","Matching on P-score"))

##########################################
# Limiting the sample with a different overlap condition - imposing overlap!
##########################################

#### Overlap condition:
ind = ps>=min(ps[d1$treat==1]) & ps<=max(ps[d1$treat==0])

# weighted
tab1=round(t(mapply(f.stat.weights,as.list(x1[treat1==1 & ind,]),as.list(x1[treat1==0 & ind,]),
                    MoreArgs=list(t0=treat1[ind],weights=w[ind]))),dig=3)
colnames(tab1) = c("Ave. Treat","Ave. control","T-test")
xtable(tab1,captioin="After IPSW: balance table",dig=3)

### Regression adjustment CPS-1, using IPSW: 
lm5 = lm(re78~(.),data=d1[ind,],weights=w[ind]);summary(lm5)

### Distribution of the weights
summary(w[treat1==0 & ind])
summary(w[treat1==1 & ind])

pdf("fig3_match.pdf",width=5,height=4)
par(mfrow=c(1,1),cex=0.7)
boxplot(ps[d1$treat==1 & ind],ps[d1$treat==0  & ind],col=c("red4","blue3"),names=c("Treatment","Control"),las=1,
        ylab="P-score",
        main="Distribution of the P-score \n NSW treatment vs CPS-1 control")
dev.off()


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








