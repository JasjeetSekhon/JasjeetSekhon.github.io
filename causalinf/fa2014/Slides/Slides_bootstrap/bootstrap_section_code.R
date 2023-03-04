##################################
# P-score and Matching example
##################################
rm(list=ls(all=TRUE))

library(xtable)
library(foreign)
library(Matching)

# Link to the data see section notes
setwd("")
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
xtable(tab,captioin="Experimental balance table",dig=3)

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
    xlim = c(xlim1,1); pchset = c(21,24,22,23); pchcolset = c("blue2","darkgreen","orange3","red2")
    
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

pdf("fig_boot1.pdf",width=4.5,height=4.5)
plot.pval(tab,legend=TRUE)
dev.off()


### Regression adjustment
lm1 = lm(re78~(.),data=d); summary(lm1)

### Bootstrap:

B=1000
beta.boot = rep(NA,B)
for (b in c(1:B)){
  if(b %% 50==0){cat("Iteration: ",b,"\n")}
  index = sample(rownames(d),length(rownames(d)),replace=TRUE)
  d0 = d[index,]
  beta.boot[b] = coef(lm(re78~(.),data=d0))[2]
}

pdf("fig_boot2.pdf",width=4.5,height=4.5)
par(cex=0.7)
plot(density(beta.boot),las=1, main="Bootstrap distribution of treatment effect",lwd=2,
     xlab=expression(beta),ylab="Density")
abline(v = quantile(beta.boot,prob=c(0.025,0.975)),col="red4",lty=2,lwd=2)
dev.off()
#abline(v = confint(lm1)[2,],col="blue4",lty=2,lwd=2)


### Parametric bootstrap:

n = dim(d)[1]
beta.boot.parm = rep(NA,B)
for (b in c(1:B)){
  if(b %% 50==0){cat("Iteration: ",b,"\n")}
  epsilon.b = sample(lm1$res,n,replace=TRUE)
  y.b = as.matrix(cbind(rep(1,n),treat,x,x[,"age"]^2)) %*% matrix(coef(lm1),ncol=1) +epsilon.b
  beta.boot.parm[b] = coef(lm(y.b~(.),data=d))[2]
}

par(cex=0.7)
plot(density(beta.boot),las=1, main="Bootstrap distribution of beta",lwd=2)
abline(v = quantile(beta.boot,prob=c(0.025,0.975)),col="red4",lty=2,lwd=2)
#abline(v = quantile(beta.boot.parm,prob=c(0.025,0.975)),col="orange4",lty=2,lwd=2)
#abline(v = confint(lm1)[2,],col="blue4",lty=2,lwd=2)

tab = rbind(quantile(beta.boot,prob=c(0.025,0.975)),
            quantile(beta.boot.parm,prob=c(0.025,0.975)),
            confint(lm1)[2,])
rownames(tab) = c("Non-parametric","Parametric","Analytical")
xtable(tab)


##########################################
# The "boot" package
##########################################
library(boot)

f.lm = function(data,index){
  return(coef(lm(re78~(.),data=data[index,]))[2])
}

boot0 = boot(data=d,statistic=f.lm,R=1000)

tab = rbind(quantile(beta.boot,prob=c(0.025,0.975)),
            quantile(beta.boot.parm,prob=c(0.025,0.975)),
            confint(lm1)[2,],
            quantile(boot0$t,prob=c(0.025,0.975)),
)
rownames(tab) = c("Non-parametric","Parametric","Analytical","Non-parametric bootstrap using boot")
xtable(tab)

pdf("fig_boot3.pdf",width=4.5,height=4.5)
par(cex=0.7)
plot(boot0)
dev.off()


ci0 = boot.ci(boot0)


















