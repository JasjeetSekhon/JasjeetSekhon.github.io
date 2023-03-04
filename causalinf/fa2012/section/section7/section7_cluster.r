setwd("/home/fdhidalgo/projects/ps236/")
library(Matching)
load("hw6.RData")

covar <- data[,names(data)[4:17]]
pscore.fmla <- as.formula(paste("foxnews2000~",paste(names(covar),collapse="+")))
unmatched.bal <- MatchBalance(pscore.fmla,data=data, nboots=50)

pscore <- glm(pscore.fmla,data = data,family = binomial(link = logit))

lp <- pscore$linear.predictor
orth.covar <- covar
# Orthogonalize covariates
for(i in 1:ncol(orth.covar)){
  orth.covar[,i]<-lm(orth.covar[,i]~lp)$residuals
}


match.data <- cbind(lp, orth.covar)


#Genetic Matching


#Using the cluster
source("/home/fdhidalgo/projects/airwaves/rcode/setupCode.R")
library(snow)
clusterCreate()

genout <-  GenMatch(Tr=data$foxnews2000, X=match.data, BalanceMatrix=covar, estimand="ATT", M=1, pop.size=5000, wait.generations=8, cluster = cl)

match.gen <-  Match( Tr=data$foxnews2000, X=match.data, estimand="ATT", Weight.matrix=genout)
MatchBalance(pscore.fmla, data = data, match.out = match.gen, nboots=50)

save(genout, file = "genmatch_results.RData")

clusterShutdown()
quit(save="no")


#Command to use on the server
#nohup Rnotify.sh section6_cluster.r 4 fdhidalgo@berkeley.edu > section6.out
