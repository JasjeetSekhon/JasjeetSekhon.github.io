rm(list = ls())
# the matching library
library(Matching)

# GenMatch help file
?GenMatch
# we have to pass in Tr and X, and have many more options

# the lalonde data set
data(lalonde)

names(lalonde)
attach(lalonde)

# we have to define what we want to match on (what we want to find weights for in our weight matrix)
X = cbind(age, educ, black, hisp, married, nodegr, u74, u75, re75, re74)

# and what we want to balance on (what we include in our loss function)
BalanceMat = cbind(age, educ, black, hisp, married, nodegr, u74, u75, re75, re74, I(re74*re75))

# pop.size defines the number of randomly generated units are in our population.  This is a very important parameter, and for publishable quality work you should set this to 1000 or higher
# max.generations is the maximum number of generations you want GenMatch to run (although, for this to be a hard limit, you have to set hard.generation.limit = TRUE)
# wait.generations is the number of generations you "wait" for improvement

genout = GenMatch(Tr = treat, X = X, BalanceMatrix = BalanceMat, estimand = "ATT", pop.size = 1000, max.generations = 5, wait.generations = 4, hard.generation.limit = TRUE)

names(genout)
genout$value
genout$par
genout$Weight.matrix
diag(genout$Weight.matrix)
genout$matches

# once we get our genmatch results, we run match passing in the weight matrix from genmatch
mout = Match(Tr = treat, X = X, estimand = "ATT", Weight.matrix = genout)

# and finally we can check the balance
mb = MatchBalance(treat ~ age + educ + black + hisp + married + nodegr + u74 + u75 + re75+ re74+ I(re74*re75), match.out=mout, nboots=500)

# if we ever get stuck in the middle of our code, we can restart genmatch from where it ends by using
# the weight matrix.  We can pass in a set of weights that it will include in its population.
# We do this using "starting.values"
genout1.5 = GenMatch(Tr = treat, X = X, BalanceMatrix = BalanceMat, estimand = "ATT", starting.values = diag(genout$Weight.matrix))

# note that if we get stuck and we can't save the object, we can still use the output.  If we have:
#GENERATION: 13
#Lexical Fit..... 2.566608e-01  2.566608e-01  2.658496e-01  2.917272e-01  3.173158e-01  3.173158e-01  5.177252e-01  7.817784e-01  7.817784e-01  8.367785e-01  8.599837e-01  9.792672e-01  9.999906e-01  9.999906e-01  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  
#unique......... 69, #Total UniqueCount: 943
#var 1:
#best............ 8.978087e+02
#mean............ 8.699059e+02
#variance........ 1.580914e+04
#var 2:
#best............ 4.602277e+02
#mean............ 4.651478e+02
#variance........ 2.543910e+03
#var 3:
#best............ 1.170968e+02
#mean............ 1.234293e+02
#variance........ 2.385922e+03
#var 4:
#best............ 3.698448e+02
#mean............ 3.785874e+02
#variance........ 2.891896e+03
#var 5:
#best............ 5.325505e+01
#mean............ 7.211983e+01
#variance........ 1.032834e+04
#var 6:
#best............ 4.427314e+02
#mean............ 4.488160e+02
#variance........ 2.900735e+03
#var 7:
#best............ 2.205042e+01
#mean............ 5.536643e+01
#variance........ 1.743466e+04
#var 8:
#best............ 5.995381e+02
#mean............ 5.978292e+02
#variance........ 1.238149e+04
#var 9:
#best............ 4.084588e+02
#mean............ 4.111752e+02
#variance........ 2.970262e+03
#var 10:
#best............ 6.915549e+02
#mean............ 6.856579e+02
#variance........ 3.822394e+03

# then the weights are given by the "best line", so if the code crashes during generation 14, but we
# have output from 13, we can reconstruct the weights as:
lastWeights = c(8.978087e+02, 4.602277e+02, 1.170968e+02, 3.698448e+02, 5.325505e+01, 4.427314e+02, 2.205042e+01, 5.995381e+02, 4.084588e+02, 6.915549e+02)

# and we could then set our starting.values as starting.values = lastWeights, and GenMatch will
# pick up where it left off because this trial should have good fit given it was
# the best fitting member of the last generation.

#Note that we also have the other usual settings such as "caliper" and "exact".  It is important to note that if you use a caliper in GenMatch, you need to send the same caliper in to Match.  Same for any other settings that you pass in to GenMatch that affect that matches, such as "estimand".

genout2 = GenMatch(Tr = treat, X = X, estimand = "ATE", caliper = 0.2)
mout2 = Match(Tr = treat, X = X, estimand = "ATE", caliper = 0.2)
mb = MatchBalance(treat ~ age + educ + black + hisp + married + nodegr + u74 + u75 + re75+ re74+ I(re74*re75), match.out=mout2, nboots=500)



###################################################################################################
# The Loss Function

# a very important aspect of GenMatch is the loss function
# the default is to minimize the maximum discrepency.  We could imagine that this
# might not be the ideal loss function for all situations.

# first, let's try to write our own function that mimics the default

# to write your own loss function, you must have a function that takes in a vector of p-values
# and returns either a scalar (only one number), or a vector of p-values that will be sorted
# lexically


# this function takes the p-values passed to it by GenMatch and sorts them and returns the minimum one
myMin = function(x) {
	cat("My Loss Function receives:", x, "\n")
	cat("My Loss Function returns:", sort(x), "\n")
	return(sort(x))
}

# note here that I am setting the seeds for GenMatch.  GenMatch uses many different seeds,
# so to get the exact same answer each time (since it is a random search algorithm), you need
# to set a seed, as well as pass in a unif.seed and an int.seed value.
set.seed(555) 
# first, let's see what GenMatch's default returns
m1 = GenMatch(Tr = treat, X = X, unif.seed=85, int.seed=5, pop.size = 16)

# Weight at Solution
# X[ 1] :	5.347816e-02
# X[ 2] :	4.364047e+02
# X[ 3] :	3.697330e+02
# X[ 4] :	1.809479e+02
# X[ 5] :	5.990888e+02
# X[ 6] :	4.307825e+02
# X[ 7] :	4.570632e+02
# X[ 8] :	4.457457e+02
# X[ 9] :	1.867992e+02
# X[10] :	7.324183e+02


# now, let's try our own "minimum" loss function
# note, though, that it is only returning the sorted vector, and genmatch will automatically
# lexically sort this.
set.seed(555) 
m2 = GenMatch(Tr = treat, X = X, unif.seed=85, int.seed=5, pop.size = 16, loss = myMin)

# Weight at Solution
# X[ 1] :	5.347816e-02
# X[ 2] :	4.364047e+02
# X[ 3] :	3.697330e+02
# X[ 4] :	1.809479e+02
# X[ 5] :	5.990888e+02
# X[ 6] :	4.307825e+02
# X[ 7] :	4.570632e+02
# X[ 8] :	4.457457e+02
# X[ 9] :	1.867992e+02
# X[10] :	7.324183e+02


# so, we see that we get the same answer.

# but, if we only returned the minimum pvalue, then the answer would change a little.

myMin = function(x) {
	cat("My Loss Function receives:", x, "\n")
	cat("My Loss Function returns:", min(x), "\n")
	return(min(x))
}

set.seed(555) 
m2.5 = GenMatch(Tr = treat, X = X, unif.seed=85, int.seed=5, pop.size = 16, loss = myMin)

# X[ 1] :	1.556710e+00
# X[ 2] :	2.994715e+02
# X[ 3] :	5.330026e+02
# X[ 4] :	2.213101e+02
# X[ 5] :	6.044632e+02
# X[ 6] :	5.787384e+02
# X[ 7] :	2.460252e+02
# X[ 8] :	3.684978e+02
# X[ 9] :	2.029297e+02
# X[10] :	5.695304e+02

# We could also write a function that does the median.  Note here that we only return the scalar.
myMedian = function(x) {
	cat("My Loss Function receives:", x, "\n")
	cat("My Loss Function returns:", median(x), "\n")
	return(median(x))
}

set.seed(555) 
m3 = GenMatch(Tr = treat, X = X, unif.seed=85, int.seed=5, pop.size = 16, loss = myMedian)


# Weight at Solution
# X[ 1] :	5.150578e+02
# X[ 2] :	2.477110e+02
# X[ 3] :	3.794624e+02
# X[ 4] :	4.518964e+02
# X[ 5] :	8.125151e+01
# X[ 6] :	9.354612e+02
# X[ 7] :	5.518071e+02
# X[ 8] :	1.947243e+02
# X[ 9] :	2.277983e+02
# X[10] :	5.163755e+02

# and, since this is a standard built in function, we can check if it is the same
set.seed(555) 
m4 = GenMatch(Tr = treat, X = X, unif.seed=85, int.seed=5, pop.size = 16, loss = median)

# Weight at Solution
# X[ 1] :	5.150578e+02
# X[ 2] :	2.477110e+02
# X[ 3] :	3.794624e+02
# X[ 4] :	4.518964e+02
# X[ 5] :	8.125151e+01
# X[ 6] :	9.354612e+02
# X[ 7] :	5.518071e+02
# X[ 8] :	1.947243e+02
# X[ 9] :	2.277983e+02
# X[10] :	5.163755e+02

# we could also do something more complex, such as make it so that we don't allow weights that make any of the variables worse in terms of balance.

# first, we have to figure out what the balance is initially.  If we choose our fit.func to be the
# default of "pvals", then our balance will be a vector of t.tests and ks.tests.  If we choose a different fit.func, such as qqmean.max, then the vector returned by genmatch will be different.  For now, we will stick with "pvals".

# t.test will conduct a t.test for the difference in means
initialize = function(X) {
	initial = NULL
	for(i in 1:ncol(X)) {
		initial = c(initial, t.test(X[,i][treat == 1], X[,i][treat == 0])$p.value)
	}
	# ks.test will conduct a ks.test for difference in distributions
	for(i in 1:ncol(X)) {
		initial = c(initial, ks.test(X[,i][treat == 1], X[,i][treat == 0])$p.value)
	}
	return(initial)
}

# note that it is important to pass in the Balance matrix, since these are the variables that
# GenMatch is using for it's measures of fitness.
initial = initialize(BalanceMat)

# now we can write our own loss function

monoBigger = function(x) {
	p.vals = x
	cat("X my Loss Gets:", x, "\n")
	if(sum(x < initial) > 0) {
		p.vals = rep(0, length(initial))	
	}
	cat("X my Loss Returns:", p.vals, "\n")	
	return(sort(p.vals))
}

set.seed(555) 
m5 = GenMatch(Tr = treat, X = X, BalanceMatrix = BalanceMat, unif.seed=85, int.seed=5, loss = monoBigger)
# looks like we can't find any thing that actually monotomically increases 11 (the ones in the
# balance matrix) variables, we can tell
# by the fact that all the weights are 1 (at least in this seed with this population)

# X[ 1] :	1.000000e+00
# X[ 2] :	1.000000e+00
# X[ 3] :	1.000000e+00
# X[ 4] :	1.000000e+00
# X[ 5] :	1.000000e+00
# X[ 6] :	1.000000e+00
# X[ 7] :	1.000000e+00
# X[ 8] :	1.000000e+00
# X[ 9] :	1.000000e+00
# X[10] :	1.000000e+00

# but maybe we can on a fewer variables:
myX = cbind(age, educ, black)
initial = initialize(myX)

initial
# 0.26594435 0.15016935 0.64735738 0.74809856 0.06287261 1.00000000

# note here that I am not passing in the Balance matrix, so GenMatch will
# automatically set the balance matrix equal to the myX matrix (which only contains 3 variables)
set.seed(555) 
m6 = GenMatch(Tr = treat, X = myX, unif.seed=85, int.seed=5, loss = monoBigger)

# here we found some units that improved weights
# X[ 1] :	8.555239e+02
# X[ 2] :	1.297584e+02
# X[ 3] :	3.964410e+01

# so, perhaps this is a reason that we don't want to use monotonic increases as a loss function
# all the time.



detach(lalonde)