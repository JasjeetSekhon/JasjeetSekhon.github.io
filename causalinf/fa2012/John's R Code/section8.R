rm(list = ls())

library(stringr)

# Code for matching with a cluster and using plot.pval

# 1. matching using a cluster of computers or chips 
#  - see http://www-personal.umich.edu/~titiunik/multiple_cpus.html for more details on 
#   matching with multiple chips 

water=read.csv(file=url("http://sekhon.berkeley.edu/causalinf/data/cross_section_wfl.csv"))

treat=water$treat
covars=water[,which(str_sub(names(water),1,5)=='y1990')]

#exclude covars with less than 3 unique values
indm=array(TRUE,ncol(covars))
for(j in 1:ncol(covars)){
sun=sort(unique(covars[,j]))
	if(length(sun)<4 & sum(sun)!=1){
		indm[j]=FALSE
	}
}

covars=covars[,indm]

mortality=water$y1999.tasatot
infect=water$y1999.tminfec

water.data=as.data.frame(cbind(treat,covars,mortality,infect))


# estimate a pscore on 10 covars
set.seed(1003)
indx=sort(sample(c(1:ncol(covars)),size=10,replace=F))

model=as.formula(paste("treat~",paste(names(covars[indx]),collapse="+")))   

# use the linear predictor
pscore=glm(model,data=water.data,family=binomial(link=logit))$linear
#pscore=glm(model,data=water.data,family=binomial(link=logit))$fitted

# plot pscore
qqplot(pscore[treat==1],pscore[treat==0],ylim=c(-5,1),xlim=c(-5,1))
abline(a=0,b=1,lty=2)

# orthogonalize covariates 
orth = covars
for(i in 1:ncol(orth)){
  orth[,i]=lm(covars[,i]~pscore)$residuals
}
   
Xmat=cbind(pscore,orth)


# using the cluster to match with genmatch

#Using the cluster
source(url("http://sekhon.berkeley.edu/causalinf/R/setupCode.R"))
library(snow)
library(Matching)
clusterCreate()  

full=as.formula(paste("treat~",paste(names(covars),collapse="+")))   

genout = GenMatch(Tr=treat, X=Xmat, BalanceMatrix=covars, estimand="ATT", M=1, pop.size=1000, wait.generations=8, max.generations=10, hard.generation.limit=T,
cluster = cl)

mout = Match(Tr=treat, X=Xmat, estimand="ATT", Weight.matrix=genout)

mb=MatchBalance(full, data = water.data, match.out = mout, nboots=50)   

clusterShutdown()


# 2. using 'plot.pval' function to make pretty balance plots
#  - see http://www-personal.umich.edu/~titiunik/R/graph.pval.public.R for more details on 
#   the original function 

source(url("http://sekhon.berkeley.edu/causalinf/R/balStatPlot.R"))

# Note that: length(mb[[1]])==dim(covars)[2]
#  - or covariates in MatchBalance must equal that to be plotted here

plot.pval(mb, covariates = names(covars))

# END