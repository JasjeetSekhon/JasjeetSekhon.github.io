library(Matching)

data(lalonde)
attach(lalonde)

#The covariates we want to match on
X = cbind(age, educ, black, hisp, married, nodegr, u74, u75, re75, re74)

#The covariates we want to obtain balance on
BalanceMat <- cbind(age, educ, black, hisp, married, nodegr, u74, u75, re75, re74,
                    I(re74*re75))

#
#Let's call GenMatch() to find the optimal weight to give each
#covariate in 'X' so as we have achieved balance on the covariates in
#'BalanceMat'. This is only an example so we want GenMatch to be quick
#so the population size has been set to be only 16 via the 'pop.size'
#option. This is *WAY* too small for actual problems.
#For details see http://sekhon.berkeley.edu/papers/MatchingJSS.pdf.
#

my.fitfunc <- function(matches, BM)
  {
    # my.fitfunc is a function to calculate balance by the mean
    # standardized difference in the eQQ plot for each variable.
    # Minimize the maximum of these differences across variables.
    # Lexical optimization is conducted.  This is the same as
    # fit.func="qqmean.max"

    # note: "matches" has three columns:
    # column 1: index of treated obs
    # column 2: index of control obs
    # column 3: weights for matched-pairs

    # note: BM is the BalanceMatrix the user passed into GenMatch
    
    index.treated <- matches[,1]
    index.control <- matches[,2]
    
    nvars <- ncol(BM)
    qqsummary   <- c(rep(NA,nvars))
    
    for (i in 1:nvars)
      {    
        
        qqfoo <- qqstats(BM[,i][index.treated], BM[,i][index.control], standardize=TRUE)
        qqsummary[i] <- qqfoo$meandiff
      } #end of for loop
    
    return(sort(qqsummary, decreasing=TRUE))
  }

genout <- GenMatch(Tr=treat, X=X, BalanceMatrix=BalanceMat, estimand="ATT", M=1,
                   pop.size=16, max.generations=10, wait.generations=1,
                   fit.func=my.fitfunc)

#The outcome variable
Y=re78/1000

#
# Now that GenMatch() has found the optimal weights, let's estimate
# our causal effect of interest using those weights
#
mout <- Match(Y=Y, Tr=treat, X=X, estimand="ATT", Weight.matrix=genout)
summary(mout)

#                        
#Let's determine if balance has actually been obtained on the variables of interest
#                        
mb <- MatchBalance(treat~age +educ+black+ hisp+ married+ nodegr+ u74+ u75+
                   re75+ re74+ I(re74*re75),
                   match.out=mout, nboots=500, ks=TRUE, mv=FALSE)
