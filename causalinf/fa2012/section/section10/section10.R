rm(list=ls(all=TRUE))
library(foreign)
library(Matching)

#load in Mark Huberty's plot code
source("balStatPlot.R")

load("./section10.RData")

covar <- data[,c(2:15)]
pscore.fmla <- as.formula(paste("abd~",paste(names(covar),collapse="+")))
abd <- data$abd
pscore <- glm(pscore.fmla, data = data, family = binomial(link = logit))


lp <- pscore$linear.predictor
bal.data <- covar
match.data <- cbind(lp,covar)
unmatched.bal <- MatchBalance(pscore.fmla,data=bal.data)

##PSCORE
match.pscore <- Match(Tr=abd,X=match.data$lp,M=1,estimand="ATT")
pscore.bal <- MatchBalance(pscore.fmla, data=bal.data, match.out = match.pscore)

###MANHALOBIS DISTANCE
match.mahn <- Match(Tr=abd, X=match.data,M=1,estimand="ATT",Weight=2)
mahn.bal <- MatchBalance(pscore.fmla, data = bal.data, match.out = match.mahn)
plot.pval(mahn.bal, covariates = names(bal.data))


###MANHALOBIS DISTANCE with Exact Matching on Age
exact <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE,FALSE, FALSE)
match.mahn.exact <- Match(Tr=abd, X=match.data,M=1,estimand="ATT",Weight=2,exact = exact)
mahn.exact.bal <- MatchBalance(pscore.fmla, data = bal.data, match.out = match.mahn.exact)
plot.pval(mahn.exact.bal, covariates = names(bal.data))


###MANHALOBIS DISTANCE with Exact Matching on Region and Caliper Set to .5
exact <- c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,FALSE, FALSE)
match.cal.exact <- Match(Tr=abd, X=match.data,M=1,estimand="ATT",Weight=2,exact = exact, caliper = .5)
cal.exact.bal <- MatchBalance(pscore.fmla, data = bal.data, match.out = match.cal.exact)
plot.pval(cal.exact.bal, covariates = names(bal.data))
#extract dropped observations
dropped.obs <- match.data[match.cal.exact$index.dropped,]



###MANHALOBIS DISTANCE with Caliper on Age
caliper <- c(100, 100, 100, 100, 100, 100, 100, 100, 100, .1, 100, 100, 100, 100, 100)
match.mahn.cal <- Match(Y=data$educ, Tr=abd, X=match.data,M=1,estimand="ATC",Weight=2, caliper = caliper)
mahn.cal.bal <- MatchBalance(pscore.fmla, data = bal.data, match.out = match.mahn.cal)
plot.pval(mahn.cal.bal, covariates = names(bal.data))
Dropped.obs <- match.data[match.mahn.cal$index.dropped,]



#GenMatch
genout <- GenMatch(Tr=abd,X=match.data,BalanceMatrix=match.data,estimand="ATT",pop.size=500)
match.gen <- Match(Y=data$educ,Tr=abd, X=match.data,M=1,estimand="ATT",Weight.matrix=genout)
gen.bal <- MatchBalance(pscore.fmla,match.out=match.gen,data=bal.data)
plot.pval(gen.bal, covariates = names(bal.data))


#Prioritize age by dividing its pvalue by 2
# the loss function
myLoss = function(x) {
	x[10] <- x[10]/2		
	return(sort(x))
}

gen.custom.out <- GenMatch(Tr=abd,X=match.data,BalanceMatrix=match.data,estimand="ATT",pop.size=500, loss = myLoss)
match.gen.custom <- Match(Y=data$educ,Tr=abd, X=match.data,M=1,estimand="ATT",Weight.matrix=gen.custom.out)
gen.custom.bal <- MatchBalance(pscore.fmla,match.out=match.gen.custom,data=bal.data)
plot.pval(gen.custom.bal, covariates = names(bal.data))

#Got good balance!
summary(match.gen.custom)



