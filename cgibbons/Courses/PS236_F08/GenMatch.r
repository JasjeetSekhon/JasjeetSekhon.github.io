### PS 236 Section, November 2008     ###
### GenMatch                          ###

## GenMatch was developed by Jas to achieve better performance in multivariate
##   matching. We are going to follow Jas' Journal of Statisitcal Software
##   article listed on the syllabus to provide a comparison of Mahalanobis
##   distance matching and p-score matching to GenMatch

## Load the matching package and the LaLonde data that comes with the package
library(Matching)
data(lalonde)

## Examine the LaLonde data
?lalonde
summary(lalonde)

## Save some LaLonde data objects
Y <- lalonde$re78        # Outcome: post-treatment income
Tr <- lalonde$treat      # Treatment status
re74 <- lalonde$re74     # Income in 1974 (pre-treatment)


## First attempt at creating a propensity score model
##   Since no link is specified the default for the binomial family is logit
glm1 <- glm(treat~age + educ + black + hisp + married + nodegr + re74 + re75,
  family=binomial, data=lalonde)

## One-to-one matching with replacement using the p-score estimated above
##   Remember, no Y is needed
##   ATT is default estimand
rr1 <- Match(Y=NULL, Tr=Tr, X=glm1$fitted);

## Let's check for balance
## MatchBalance doesn't actually do a function (the __ ~ __ + __ . . . notation)
##   but it's easy for R to use.
##   The 'I()' function tells R to create the variable inside. For example,
##   lm(Y ~ X + Z) would estimate Y = X * a + Z * b, whereas
##   lm(Y ~ I(X + Z)) would estimate Y = (X + Z) * c
##   Note also that we check for balance on many things that we didn't include
##   in the p-score or matching model. This serves to provide a good check---
##   if we achieve balance on a covariate that we didn't ask for balance on,
##   that's a good sign.
MatchBalance(
  # Specify what to look for balance on
  treat~age + I(age^2) + educ + I(educ^2) + black +
  hisp + married + nodegr + re74 + I(re74^2) + re75 + I(re75^2) +
  u74 + u75 + I(re74*re75) + I(age*nodegr) + I(educ*re74) + I(educ*re75),
  # Specify that matching output that gives the matched dataset
  match.out=rr1,
  # Specify the number of bootstrap simulations in the KS test
  nboots=1000,
  # Specify the data source for the covariates in the "model" above
  data=lalonde)
  
## In interpreting these results, we would like the after-matching means and
##   the ratio of variances of treatment and control to be the same. The output
##   of MatchBalance provides these raw statistics ('mean treatment', 'mean
##   control', 'var ratio (Tr/Co)'). Also are the mean, median, and maximum
##   differences between between the actual ('raw') values in the QQ-plot and
##   the standardized ('eCDF') version. These figures give the statistics for
##   the distance between the points in the following QQ-plot and the 45-degree
##   line (re74 given as a covariate example).
## Let's use a (somewhat) conrete
##   illustration of what a QQ-plot does: if there are ten treated and 20
##   controls in your matched sample, the QQ-plot orders re74 for each group
##   and plots the first treated against the second control, the second treated
##   against the fourth control, etc. This plots the first decile value of
##   treated against the first decile of control, etc. The number of breaks
##   (n-tiles) is based upon the group with the fewest observations.
##   MatchBalance takes the difference between these two values (NB: NOT the
##   matched pairs) and calculates the mean, median, and maximum among them.
qqplot(re74[rr1$index.control], re74[rr1$index.treated])
abline(coef=c(0,1), col=2)

## In summary, we want the standardized mean difference ('std mean diff') to
##   be close to 0, the QQ-plot statistics to be close to 0, the variance
##   ratio to be near 1, the t-test and KS p-values to be near 1 (implying that
##   you cannot reject the null hypothesis of mean or distribution equality
##   with high certainty). The bootrapped KS value is preferred in general
##   because it deals well with distributions that have point masses.

## The matching result above is quite poor. Let's try to improve it by using
##   a better p-score model and GenMatch.

## This is the p-score model that Dehejia and Wahba use
## Note that GenMatch doesn't require a p-score, but it can be helpful.
dw.pscore <- glm(treat~age + I(age^2) + educ + I(educ^2) + black +
  hisp + married + nodegr + re74 + I(re74^2) + re75 +
  I(re75^2) + u74 + u75, family=binomial, data=lalonde)

## Create a matrix of variables to match on (X) and for GenMatch to test
##   and maximize balance on (BalanceMatrix)
attach(lalonde)
X <- cbind(age, educ, black, hisp, married, nodegr, re74, re75, u74, u75)
BalanceMatrix <- cbind(age, I(age^2), educ, I(educ^2), black,
  hisp, married, nodegr, re74 , I(re74^2), re75,
  I(re75^2), u74, u75, I(re74*re75), I(age*nodegr),
  I(educ*re74), I(educ*re75))
detach(lalonde)

## Let's do GenMatch
gen1 <- GenMatch(
  # Treatment and the matching variables are specified as before, but note
  #   that Y cannot be inputted (rather than just being optional as above).
  Tr=Tr,
  # Unlike the J. Stat. Software article, we'll match on the Dehejia-Wahba
  #   p-score as well
  X=cbind(dw.pscore$linear.predictors, X),
  # Give the matrix of covariates for which we desire balance
  BalanceMatrix=BalanceMatrix,
  # GenMatch runs in iterations, creating 'populations' of potential weighting
  #   matrices (the 'W' matrix in the equation of page 6 of the J. Stat.
  #   Software article). The higher the 'pop.size', the more potential matrices
  #   that will be examined in an iteration, but the longer it will take. Jas'
  #   article uses 10000; we'll use 1000, while 100 is the default. Note that
  #   this parameter is perhaps the most important to achieving good balance
  #   ('GenMatch' converges to the true optimum asymtotically as a function of
  #   population size, not generations (see below)).
  pop.size=1000,
  # 'max.generations' specifies that maximum number of iterations that will run.
  #   Again, more is better, but we'll use 1000 (default is 100). Note that this
  #   is a "soft" limit and is only binding if 'hard.generation.limit' is set
  #   to TRUE (default is FALSE).
  max.generations=1000,
  # 'wait.generations' is the number of generations run where the value of the
  #   loss function that we are minimizing (e.g., smallest p-value) doesn't
  #   change ('tolerance' sets how big "no change" can be) before
  #   GenMatch determines that balance is optimal. We keep the default of 4.
  wait.generations=4
  # By default, GenMatch runs non-exact, no caliper matching with replacement
  #   and ties, but these defaults can be changed.
  # See '?GenMatch' for more details, but these are the parameters that you
  #   should focus on for the most part.
  # Lastly, you should note that you'll get generatino-by-generation results
  #   reporting of some key statistics. You can control the amount of
  #   reporting using 'print.level'.
)

## GenMatch only finds the optimal weighting matrix---it doesn't give any
##   statistics. Once GenMatch has this matrix, we can supply it to
##   'Match' and 'MatchBalance' to get our results. We'll look at the same
##   balance statistics as above.
mgen1 <- Match(Y=Y, Tr=Tr, X=X, Weight.matrix=gen1)
MatchBalance(treat~age + I(age^2) + educ + I(educ^2) + black +
  hisp + married + nodegr + re74 + I(re74^2) + re75 + I(re75^2) +
  u74 + u75 + I(re74*re75) + I(age*nodegr) + I(educ*re74) + I(educ*re75),
  data=lalonde, match.out=mgen1, nboots=1000)

## Now, the results and objects returned by 'Match' can be used as we have
##   already seen.

## Get our estimated ATT. Recall that the experimental benchmark for the
##   subsample that includes data on re74 is $1794.
summary(mgen1)

