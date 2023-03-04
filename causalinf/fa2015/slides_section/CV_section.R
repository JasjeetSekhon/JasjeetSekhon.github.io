
rm(list=ls(all=TRUE))
library(ggplot2)
library(ggthemes)

setwd("/home/yotam/Dropbox/causalinf.private/fa2015/slides/slides6_CV")

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
# P-score example: PSID
##################################
library(AER)
library(reshape)

data("PSID1982", package = "AER")
d = PSID1982

rmse = function(x1,x2){return(mean((x1-x2)^2))}

covnames = colnames(d)[! colnames(d) %in% c("wage","union")]
d$tr = (d$union=="yes")*1

L0=100 # number of repetitions
rmse.ols <- rmse.probit1 <- rmse.probit2<-rmse.logit1 <- rmse.logit2 <- converge1<-converge2 <-rep(NA,L0)
data = d[,c(covnames,"tr")]
for (j in c(1:L0)){
  id = sample(c(1:dim(d)[1]),round(dim(d)[1]*0.5))
  
  ols <- lm(tr~(.),data=data[id,])
  logit1 <- glm(tr~(.),data=data[id,],family=binomial(link="logit"))
  logit2 <- glm(tr~(.)^2,data=data[id,],family=binomial(link="logit"))
  probit1 <- glm(tr~(.),data=data[id,],family=binomial(link="probit"))
  probit2 <- glm(tr~(.)^2,data=data[id,],family=binomial(link="probit"))
  
  rmse.ols[j] = rmse(predict(ols,newdata=data[-id,],type="response"),d$tr[-id])
  rmse.logit1[j] = rmse(predict(logit1,newdata=data[-id,],type="response"),d$tr[-id])  
  rmse.logit2[j] = rmse(predict(logit2,newdata=data[-id,],type="response"),d$tr[-id])
  rmse.probit1[j] = rmse(predict(probit1,newdata=data[-id,],type="response"),d$tr[-id])  
  rmse.probit2[j] = rmse(predict(probit2,newdata=data[-id,],type="response"),d$tr[-id])
}

dp <- data.frame(ols=rmse.ols,logit1=rmse.logit1,logit2=rmse.logit2,probit1=rmse.probit1,probit2=rmse.probit2) 
dp = melt(dp,variable_name="model")
dp = rename(dp,c(value="rmse"))

pdf("union_pscore_cv1.pdf")
ggplot(dp,aes(y=rmse,x=model,fill=model))+
  geom_boxplot()+
  labs(x="",y="RMSE",title="RMSE of predicting union membership")+
  theme_bw()+ 
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  scale_fill_economist()+
  guides(fill=guide_legend(title="Model type: "))+ # adding legend title
  theme(legend.position="bottom")
dev.off()

pdf("union_pscore_cv2.pdf")
ggplot(dp[dp$model %in% c("ols","logit1","probit1"),],aes(y=rmse,x=model,fill=model))+
  geom_boxplot()+
  labs(x="",y="RMSE",title="RMSE of predicting union membership")+
  theme_bw()+ 
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  scale_fill_economist()+
  guides(fill=guide_legend(title="Model type: "))+ # adding legend title
  theme(legend.position="bottom")
dev.off()

tapply(dp$rmse,dp$model,mean)

#######################################
# P-score example: PSID - K-folds CV
#######################################

folds.num <- 10
d$fold <- cut(c(1:dim(d)[1]),breaks=folds.num)
levels(d$fold) = c(1:folds.num)

rmse.ols <- rmse.probit1 <- rmse.probit2<-rmse.logit1 <- rmse.logit2 <- converge1<-converge2 <-rep(NA,folds.num)
data = d[,c(covnames,"tr")]
for (j in c(1:folds.num)){
  id = which(d$fold!=j)
  
  ols <- lm(tr~(.),data=data[id,])
  logit1 <- glm(tr~(.),data=data[id,],family=binomial(link="logit"))
  logit2 <- glm(tr~(.)^2,data=data[id,],family=binomial(link="logit"))
  probit1 <- glm(tr~(.),data=data[id,],family=binomial(link="probit"))
  probit2 <- glm(tr~(.)^2,data=data[id,],family=binomial(link="probit"))
  
  rmse.ols[j] = rmse(predict(ols,newdata=data[-id,],type="response"),d$tr[-id])
  rmse.logit1[j] = rmse(predict(logit1,newdata=data[-id,],type="response"),d$tr[-id])  
  rmse.logit2[j] = rmse(predict(logit2,newdata=data[-id,],type="response"),d$tr[-id])
  rmse.probit1[j] = rmse(predict(probit1,newdata=data[-id,],type="response"),d$tr[-id])  
  rmse.probit2[j] = rmse(predict(probit2,newdata=data[-id,],type="response"),d$tr[-id])
}

dp <- data.frame(ols=rmse.ols,logit1=rmse.logit1,logit2=rmse.logit2,probit1=rmse.probit1,probit2=rmse.probit2) 
dp = melt(dp,variable_name="model")
dp = rename(dp,c(value="rmse"))

tapply(dp$rmse,dp$model,mean)


#############################################
# LASSO example: PSID - Need to complete
#############################################
library(AER)
library(reshape)
library(glmnet)

data("PSID1982", package = "AER")
d = PSID1982

y = as.matrix(d$wage)
x = d[,covnames]
x$occupation<- (x$occupation=="white")*1
x$industry <- (x$industry=="yes")*1
x$south <- (x$south=="yes")*1
x$smsa <- (x$smsa=="yes")*1
x$ethnicity <- (x$ethnicity=="afam")*1
x$married <- (x$married=="yes")*1
x$gender <- (x$gender=="male")*1
x = as.matrix(x)

### Split-sample CV
#folds.num <- 10
#fold <- cut(c(1:dim(d)[1]),breaks=folds.num)
#levels(fold) = c(1:folds.num)

#rmse.vec <-rep(NA,folds.num)
#for (j in c(1:folds.num)){
#  id = which(fold!=j)
#  lasso <- glmnet(x,y,family="gaussian",alpha=1,lambda=lambda.vec[])
#}


### K-folds CV








predict(lasso,s=10)

lasso$beta




##################################
# P-score example: welder
##################################
rm(list=ls(all=TRUE))
library(xtable)

load("~/Dropbox/causalinf.private/fa2015/slides/slides5_pscore/welder.RData")
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








