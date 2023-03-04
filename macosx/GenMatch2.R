library(Matching)
gc()
gc()
set.seed(38913)

data(lalonde)
attach(lalonde)

#The covariates we want to match on
X = cbind(age, educ, black, hisp, married, nodegr, u74, u75, re75, re74);

Xbig <- rbind(X, X)
dim(Xbig)

Ybig <- c(treat, treat)
length(Ybig)

#Let's call GenMatch() to find the optimal weight to give each
#covariate in 'X' so as we have achieved balance on the covariates in
#'BalanceMat'. This is only an example so we want GenMatch to be quick
#to the population size has been set to be only 15 via the 'pop.size'
#option.

foo <- print(system.time(GenMatch(Tr=Ybig, X=Xbig, BalanceMatrix=Xbig, estimand="ATE", M=1,
                                  pop.size=16, max.generations=10, wait.generations=1)))
cat("GENMATCH: 1",foo,"\n")

foo <- print(system.time(GenMatch(Tr=Ybig, X=Xbig, BalanceMatrix=Xbig, estimand="ATE", M=1,
                                  pop.size=16, max.generations=10, wait.generations=1)))
cat("GENMATCH: 2",foo,"\n")

foo <- print(system.time(GenMatch(Tr=Ybig, X=Xbig, BalanceMatrix=Xbig, estimand="ATE", M=1,
                                  pop.size=16, max.generations=10, wait.generations=1)))
cat("GENMATCH: 3",foo,"\n")
