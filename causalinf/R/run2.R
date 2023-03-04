#library(Matching, lib="~/libs")
source("GenMatchMP1.R")

#example(GenMatch)

data(lalonde)
attach(lalonde)

#The covariates we want to match on
X = cbind(age, educ, black, hisp, married, nodegr, u74, u75, re75, re74);
X <- rbind(X,X,X)
treat <- c(treat,treat,treat)
Y <- c(re78,re78,re78)

genout <- GenMatch(Tr=treat, X=X, estimand="ATE",
                   pop.size=100, max.generations=3, wait.generations=3,
                   hard.generation.limit=TRUE)


