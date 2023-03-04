# Example of using GenMatch with lalonde2.RData.  This is modified
# from the example included in the GenMatch help page.
#
# For output see http://sekhon.berkeley.edu/causalinf/R/GenMatch1.Rout
#

library(Matching)

#download the data file from the web.  It is better to copy it to your
#computer so you don't have to keep doing this!!
load(url("http://sekhon.berkeley.edu/causalinf/R/lalonde2.RData"))
#if the file is located locally
#load("lalonde2.RData")

#let's attach the dataset
attach(lalonde2)

# Let's print out the variable names:

#> names(lalonde2)
# [1] "treat"     "age"       "education" "black"     "hispan"    "married"  
# [7] "nodegree"  "re74"      "re75"      "re78"      "u74"       "u75"


#The covariates we want to match on
X = cbind(age, education, black, hispan, married, nodegree, u74, u75, re75, re74);

#The covariates we want to obtain balance on
#ARE THESE ENOUGH??  WHAT ELSE MAY WE WANT TO LOOK AT?
BalanceMat <- cbind(age, education, black, hispan, married, nodegree, u74, u75, re75, re74,
                    I(re74*re75));

#Let's call GenMatch() to find the optimal weight to give each
#covariate in 'X' so as we have achieved balance on the covariates in
#'BalanceMat'. For details on 'pop.size', 'max.generations' and
#'wait.generations' see the GenMatch help and the genoud help.  These
#three are the most important 'GenMatch' specific options aside from
#'BalanceMatrix'.

###############
# WARNING 1: THESE SETTINGS ARE FOR A QUICK RUN.  FOR EXAMPLE, A LARGER POP.SIZE IS NEEDED.
# WARNING 2: NO PSCORE IS BEING USED IN THIS EXAMPLE
# WARNING 3: NO ORTHOGONALIZATION IS BEING DONE (THERE IS AFTER ALL NO PSCORE)
###############

genout <- GenMatch(Tr=treat, X=X, BalanceMatrix=BalanceMat, estimand="ATT", M=1,
                   pop.size=100, wait.generations=1)

#The outcome variable
Y=re78/1000;

# Now that GenMatch() has found some weights, let's estimate
# our causal effect of interest using those weights
mout <- Match(Y=Y, Tr=treat, X=X, estimand="ATT", Weight.matrix=genout)
summary(mout)

#                        
#Let's determine if balance has actually been obtained on the variables of interest
#WHAT ELSE MAY BE OF INTEREST?                        
mb <- MatchBalance(treat~age +education+black+ hispan+ married+ nodegree+ u74+ u75+
                   re75+ re74+ I(re74*re75),
                   match.out=mout, nboots=1000)



