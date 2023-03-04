
rm(list=ls(all=TRUE))
library(ggplot2)

setwd("")

d1 = read.csv("growth_1963_1997_cln.csv",header=TRUE,fill=TRUE)
d2 = read.csv("growth_1997_2013_cln.csv",header=TRUE,fill=TRUE)
d2 = d2[-1,]

growth = c(d1[d1$Area=="California",c(-1,-2)],
d2[d2$Area=="California",c(-1,-2)])
growth = unlist(growth)

year=c(1964:2013)


pdf(file="cv_example.pdf",width=4.5,height=4.5)
ggplot(data.frame(year,growth),aes(x=year,y=growth))+geom_point()+
  stat_smooth(method="lm")+theme_bw()
dev.off()



####################################################################
# Example:
####################################################################

f.lag = function(lag,y){ # generating lags:
  X = matrix(999,ncol=lag,nrow=length(y))
  for (j in c(1:lag)){
    X[,j] = c(rep(NA,j),y[1:(length(y)-j)])
  }
  colnames(X) = rep("",lag)
  return(X)
}

data0 = data.frame(y=growth,f.lag(3,y=growth))
colnames(data0) = c("y",paste("lag",c(1:3),sep=""))

# model 1
lm1 = lm(y~lag1,data=data0)
summary(lm1)

# model 2
lm2 = lm(y~lag1+lag2,data=data0)
summary(lm2)

# model 3
lm3 = lm(y~lag1+lag2+lag3,data=data0)
summary(lm3)

### Summary of results:
library(texreg)
screenreg(list(lm1,lm2,lm3))
texreg(list(lm1,lm2,lm3),digits=3,
       custom.coef.names=c("Intercept",paste("Lag",c(1:3))))

## just for fun :)
plotreg(lm3)


#### Wald test (F-test)
# note: In order to have the same data across models we need to drop the first three observaions:
data0.anova = data0[-c(1:3),]
lm1.anova = lm(y~lag1,data=data0.anova)
lm2.anova = lm(y~lag1+lag2,data=data0.anova)
lm3.anova = lm(y~lag1+lag2+lag3,data=data0.anova)

anova(lm1.anova,lm2.anova)
anova(lm1.anova,lm3.anova)
anova(lm2.anova,lm3.anova)

####################################################################
# Iterative CV
####################################################################

rmse = function(x1,x2){
  return(sqrt(mean((x1-x2)^2)))
}

for (m in c(5,10,15)){
  # "m" the minimal number of periods:
  rmse1 <- rmse2 <- rmse3 <- rmse4 <- rmse5 <- rep(999,length(growth)-m-1)
  
  
  for (t in c(1:(length(growth)-m-2))){
    cat("Period: ",t,"\n")
    ar1 = ar(growth[1:(m+t-1)],order.max=1,aic=FALSE,demean=FALSE)
    rmse1[t] = rmse(predict(ar1,newdata=growth[(m+t):length(growth)],n.ahead=length(growth[(m+t):length(growth)]))$pred,
                    growth[(m+t):length(growth)])
    
    ar2 = ar(growth[1:(m+t-1)],order.max=2,aic=FALSE,demean=FALSE)
    rmse2[t] = rmse(predict(ar2,newdata=growth[(m+t):length(growth)],n.ahead=length(growth[(m+t):length(growth)]))$pred,
                    growth[(m+t):length(growth)])
    
    ar3 = ar(growth[1:(m+t-1)],order.max=3,aic=FALSE,demean=FALSE)
    rmse3[t] = rmse(predict(ar3,newdata=growth[(m+t):length(growth)],n.ahead=length(growth[(m+t):length(growth)]))$pred,
                    growth[(m+t):length(growth)])
  }
  
  assign(paste("tab",m,sep=""),matrix(unlist(lapply(list(rmse1,rmse2,rmse3),mean)),ncol=3)) 
}


### Results:
library(xtable)

tab = round(rbind(tab5,tab10,tab15),dig=3)
colnames(tab) = paste("Model",c(1:3))
rownames(tab) = paste("M =",c(5,10,15))

xtable(tab,align=rep("c",4),dig=3)



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

#### Choosing which model is correct using cross-validation:

rmse = function(x1,x2){
  return(mean((x1-x2)^2))
}

L0=100 # number of repetitions
rmse.model.1 <- rmse.model.2 <- converge1<-converge2 <-rep(NA,L0)
a = data.frame(treat=treat,x)
for (j in c(1:L0)){
  id = sample(c(1:dim(d)[1]),round(dim(d)[1]*0.5))
  
  ps.model1 <- glm(treat~(.),data=a[id,],family=binomial(link=logit))
  ps.model2 <- glm(treat~(.)^2,data=a[id,],family=binomial(link=logit))
  
  rmse.model.1[j] = rmse(predict(ps.model1,newdata=a[-id,],type="response"),a$treat[-id])
  
  rmse.model.2[j] = rmse(predict(ps.model2,newdata=a[-id,],type="response"),a$treat[-id])
}

tab = rbind(c(mean(rmse.model.1),mean(rmse.model.2)),c(median(rmse.model.1),median(rmse.model.2)))
colnames(tab)=c("Model 1","Model 2")
rownames(tab)=c("Mean","Median")
xtable(tab)

pdf(file="cv_example2.pdf",width=4.5,height=4.5)
par(cex=0.7)
boxplot(rmse.model.1,rmse.model.2,las=1,names=c("Model 1","Model 2"))
dev.off()








