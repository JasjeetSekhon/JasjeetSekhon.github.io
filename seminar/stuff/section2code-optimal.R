load("./section3data.RData")
names(data)
data

library(MatchIt)
?matchit

# Optimal Pair Matching
match.1 = matchit(treat ~ age + black + smoker, data = data, method = "optimal", distance = "logit", ratio = 1)

names(match.1)
match.1$match.matrix
match.1$subclass

# how do we created a matched data set?
index.treated = c(1:nrow(data))[as.numeric(rownames(match.1$match.matrix))]
index.treated

index.control = c(1:nrow(data))[as.numeric(match.1$match.matrix)]
index.control

match.1.matched = rbind(data[index.treated,], data[index.control, ])

boxplot(age ~ treat, data = data)
boxplot(age[c(index.treated, index.control)] ~ treat[c(index.treated, index.control)], data = data)

# Optimal Multiple Control Matching

# by default, the ratio is 1:1 matching.  You can change it with ratio
match.temp = matchit(treat ~ age + black + smoker, data = data, method = "optimal", distance = "mahalanobis", ratio = 2)

#however we don't have enough data to do 2:1 matching.  Ideally we would do variable control matching, so it would put either 0, 1, or 2 controls for each treated, however this function doesn't seem to be working in their code.


# the code for this doesn't work!  I'm still trying to figure out
# what King et. al did to screw it up.
# match.temp.2 = matchit(treat ~ age + black + smoker, data = data, method = "optimal", distance = "mahalanobis", min.controls = 0, max.controls = 2)

# we can do it in fullmatch from optmach, but it is a terrible output!  Ask me if you want to do this, I have some ideas
# on how to deal with this format.
match.temp.3 = fullmatch(mdist(treat ~ age + black + smoker, data = data), min.controls = 0, max.controls = 2)

# Full Optimal Matching
match.2 = matchit(treat ~ age + black + smoker, data = data, method = "full", distance = "mahalanobis")

# full matching won't return a match.matrix, so it is much more difficult to create a matched data set, but still possible.
# the matches are stored in the subclasses - aka strata, so from each strata we grab the treated ones and the control ones.

num.strata = max(match.2$subclass)

for(i in 1:num.strata) {
	print(match.2$subclass[match.2$subclass == i])
}

my.matches = NULL
strata = NULL


for(i in 1:num.strata) {
	index.strata = as.numeric(names(match.2$subclass)[match.2$subclass == i])
	my.matches = rbind(my.matches, cbind(index.strata[which(treat[index.strata] == 1)], index.strata[which(treat[index.strata] == 0)]))
	strata = c(strata, rep(i, nrow(cbind(index.strata[which(treat[index.strata] == 1)], index.strata[which(treat[index.strata] == 0)]))))
}

colnames(my.matches) = c("index.treated", "index.control")
rownames(my.matches) = strata

my.matches
match.2$weights

# if you want to reweight the data, you simply reweight it by the $weights returned from matchit
mean.age.treated = mean(age[treat == 1]*match.2$weights[treat == 1])
mean.age.control = mean(age[treat == 0]*match.2$weights[treat == 0])


# Non-bipartite Matching
library(nbpMatching)
?nonbimatch

dm = matrix(c(0, 106, 119, 231, 110, 101, 106, 0, 207, 126, 192, 68, 119, 207, 0, 156, 247, 25, 231, 126, 156, 0, 34, 67, 110, 192, 247, 34, 0, 212, 101, 68, 25, 67, 212, 0), nrow = 6, ncol = 6)

# their code requires that you use their function "distancematrix" to prepare your matrix for the matching
dm = distancematrix(dm)

nbm.1 = nonbimatch(dm)
nbm.1
nbm.1$halves
total.dist = sum(nbm.1$halves[, 3])
total.dist

# odd number of observations:
dm.2 = dm[1:5, 1:5]
dm.2 = distancematrix(dm.2)
#notice that it destroys the names... we can fix this by adding the sink ourselves before calling distancematrix
dm.2 = dm[1:5, 1:5]
dm.2 = cbind(dm.2, rep(0, 5))
dm.2 = rbind(dm.2, rep(0, 6))
dm.2 = distancematrix(dm.2)


nbm.2 = nonbimatch(dm.2)
nbm.2