load(file=url("http://sekhon.berkeley.edu/causalinf/data/section4.RData")) 

#data

means = NULL
for(i in c(2:4)) {
	means <- rbind(means, tapply(data[,i], as.factor(data$treat), mean))
}
means <- as.data.frame(means)
names(means) <- c("control", "treat")
row.names(means) <- c("age", "black", "smoker")


# glm refers to the generalized linear model
# we see that we have to pass in a formula, a family with link function
# and we can pass in a data set (which allows us to not use the $ sign to access the variables in our data set)
est <- glm(treat ~ age + black + smoker, family = binomial(link = logit), data = data)

summary(est)
names(est)

pscore <- est$fitted.values
pscore

# the distributions:
boxplot(pscore ~ data$treat, main = "Pscores by Treated and Control")

plot(density(pscore[data$treat == 1]), col = "blue", xlim = c(0, 1), main = "Pscores by Treated and Control")
points(density(pscore[data$treat == 0]), type = "l", col = "red")
legend("topleft", c("treated", "control"), col = c(4, 2), lty = c(1, 1))

par(mfrow = c(1, 3))
plot(density(data$age[data$treat == 1]), col = "blue", main = "Age by Treated and Control")
points(density(data$age[data$treat == 0]), type = "l", col = "red")
legend("topleft", c("treated", "control"), col = c(4, 2), lty = c(1, 1))

plot(density(data$black[data$treat == 1]), col = "blue", main = "Race by Treated and Control")
points(density(data$black[data$treat == 0]), type = "l", col = "red")
legend("topleft", c("treated", "control"), col = c(4, 2), lty = c(1, 1))

plot(density(data$smoker[data$treat == 1]), col = "blue", main = "Smoker by Treated and Control")
points(density(data$smoker[data$treat == 0]), type = "l", col = "red")
legend("topleft", c("treated", "control"), col = c(4, 2), lty = c(1, 1))
par(mfrow = c(1, 1))

# so, now lets add the pscore to our data set
data$pscore = pscore
data

# so, what does a pscore look like for an individual?
data[5, ]

# here we see that for observation number 5, the individual is a 35 white smoker (who received treatment).  Their pscore is ~ 0.65.  This means that the person has about a 65% chance of being assigned to treatment conditional on their observable covariates.

data[33, ]

# here we see that observation number 33, a 36 year old white smoker (who did not receive treatment) has a pscore of ~ 0.64.

mean(data$pscore[data$treat == 1])
mean(data$pscore[data$treat == 0])

# matching
library(Matching)
?Match

# here we are matching with replacement without a caliper
m1 <- Match(Tr = data$treat, X = data$pscore, estimand = "ATT")
names(m1)

m1rounded = Match(Tr = data$treat, X = round(data$pscore, 2), estimand = "ATT")


# lets look at our matched data along with the weights
cbind(data[c(m1$index.treated, m1$index.control), ], m1$weights)

# what does the density look like now?
boxplot(pscore[c(m1$index.treated, m1$index.control)] ~ c(data$treat[m1$index.treated], data$treat[m1$index.control]))

par(mfrow=c(1,1))
plot(density(pscore[m1$index.treated]), col = "blue", xlim = c(0, 1), main = "Pscores by Treated and Control \n Post Matching")
points(density(pscore[m1$index.control]), type = "l", col = "red")
legend("topleft", c("treated", "control"), col = c(4, 2), lty = c(1, 1))

par(mfrow = c(1, 3))
plot(density(data$age[m1$index.treated]), col = "blue", main = "Age by Treated and Control \n After Matching")
points(density(data$age[m1$index.control]), type = "l", col = "red")
legend("topleft", c("treated", "control"), col = c(4, 2), lty = c(1, 1))

plot(density(data$black[m1$index.treated]), col = "blue", main = "Race by Treated and Control \n After Matching")
points(density(data$black[m1$index.control]), type = "l", col = "red")
legend("topleft", c("treated", "control"), col = c(4, 2), lty = c(1, 1))

plot(density(data$smoker[m1$index.treated]), col = "blue", main = "Smoker by Treated and Control \n After Matching")
points(density(data$smoker[m1$index.control]), type = "l", col = "red")
legend("topleft", c("treated", "control"), col = c(4, 2), lty = c(1, 1))
par(mfrow = c(1, 1))

# and calculate the att
att <- mean(data$dpc[m1$index.treated] * m1$weights) - mean(data$dpc[m1$index.control] * m1$weights)

# lets let Match calculate the ATT for us
m2 <- Match(Y = data$dpc, Tr = data$treat, X = data$pscore, estimand = "ATT")

# a summary
summary(m2)

m2$est
att

m3 <- Match(Y = data$dpc, Tr = data$treat, X = data$pscore, estimand = "ATT", caliper = 0.1)

summary(m3)

cbind(data[c(m3$index.treated, m3$index.control), ], m3$weights)

# which observations were dropped?
data[m3$index.dropped, ]

m3$est
att

# notice that our ATT has changed... this is because we dropped observations.  This means our estimand isnt estimating what we think it is... but what is it measuring?

# do the distributions look better?
boxplot(pscore[c(m3$index.treated, m3$index.control)] ~ c(data$treat[m3$index.treated], data$treat[m3$index.control]))

plot(density(pscore[m3$index.treated]), col = "blue", xlim = c(0, 1), main = "Pscores by Treated and Control - Post Matching")
points(density(pscore[m3$index.control]), type = "l", col = "red")
legend("topleft", c("treated", "control"), col = c(4, 2), lty = c(1, 1))


par(mfrow = c(1, 3))
par(mfrow = c(3, 1))
plot(density(data$age[m3$index.treated]), col = "blue", main = "Age by Treated and Control \n After Matching")
points(density(data$age[m3$index.control]), type = "l", col = "red")
legend("topleft", c("treated", "control"), col = c(4, 2), lty = c(1, 1))

plot(density(data$black[m3$index.treated]), col = "blue", main = "Race by Treated and Control \n After Matching")
points(density(data$black[m3$index.control]), type = "l", col = "red")
legend("topleft", c("treated", "control"), col = c(4, 2), lty = c(1, 1))

plot(density(data$smoker[m3$index.treated]), col = "blue", main = "Smoker by Treated and Control \n After Matching")
points(density(data$smoker[m3$index.control]), type = "l", col = "red")
legend("topleft", c("treated", "control"), col = c(4, 2), lty = c(1, 1))
par(mfrow = c(1, 1))

####Another example

library(foreign)
data <- read.dta("http://sekhon.berkeley.edu/causalinf/data/BalanceTest(n147).dta")

covar <- data[,c("poverty", "pop2000", "tariq", "lelev", "iso", "lnn", "garrison","reb")]
pscore.fmla <- as.formula(paste("treat~",paste(names(covar),collapse="+")))
treat <- data$treat
pscore <- glm(pscore.fmla, data = data,family = binomial(link = logit))
lp <- pscore$linear.predictor


bal.data <- covar
match.data <- cbind(lp,covar)

unmatched.bal <- MatchBalance(pscore.fmla,data=bal.data, ks=FALSE)

bal.stats <- data.frame(covar= names(bal.data))
for (i in 1:8){
  bal.stats$std.diff[i] <- abs(as.numeric((unmatched.bal$BeforeMatching)[[i]][1]))
  bal.stats$diff.means.p[i] <- (as.numeric((unmatched.bal$BeforeMatching)[[i]][7]))
  bal.stats$var.ratios[i] <- (as.numeric((unmatched.bal$BeforeMatching)[[i]][8]))
}


dotchart(bal.stats$std.diff[order(bal.stats$std.diff)],labels=bal.stats$covar[order(bal.stats$std.diff)],main="Standardized Difference",color="black",pch=19)
dotchart(bal.stats$diff.means.p[order(bal.stats$diff.means.p)],labels=bal.stats$covar[order(bal.stats$diff.means.p)],main="Difference in Means P Value",color="black",pch=19)
dotchart(bal.stats$var.ratios[order(bal.stats$var.ratios)],labels=bal.stats$covar[order(bal.stats$var.ratios)],main="Variance Ratios",color="black",pch=19)


match.pscore <- Match(Tr=treat,X=match.data$lp,M=1,estimand="ATT")
pscore.bal <- MatchBalance(pscore.fmla, data=bal.data, match.out = match.pscore, ks=FALSE)

bal.stats <- data.frame(covar= names(bal.data))
for (i in 1:8){
  bal.stats$bm.std.diff[i] <- abs(as.numeric((pscore.bal$BeforeMatching)[[i]][1]))
  bal.stats$bm.diff.means.p[i] <- (as.numeric((pscore.bal$BeforeMatching)[[i]][7]))
  bal.stats$bm.var.ratios[i] <- (as.numeric((pscore.bal$BeforeMatching)[[i]][8]))
  bal.stats$am.std.diff[i] <- abs(as.numeric((pscore.bal$AfterMatching)[[i]][1]))
  bal.stats$am.diff.means.p[i] <- (as.numeric((pscore.bal$AfterMatching)[[i]][7]))
  bal.stats$am.var.ratios[i] <- (as.numeric((pscore.bal$AfterMatching)[[i]][8]))
}

dotchart(bal.stats$am.std.diff[order(bal.stats$bm.std.diff)],labels=bal.stats$covar[order(bal.stats$bm.std.diff)],main="Standardized Difference",color="red",pch=19, ylim=c(12,28))
points(bal.stats$bm.std.diff[order(bal.stats$bm.std.diff)],1:8,col="black",pch=19)

dotchart(bal.stats$bm.diff.means.p[order(bal.stats$bm.diff.means.p)],labels=bal.stats$covar[order(bal.stats$bm.diff.means.p)],main="Difference in Means P Value",color="black",pch=19,xlim=c(0,1))
points(bal.stats$am.diff.means.p[order(bal.stats$bm.diff.means.p)],1:8,col="red",pch=19)
dotchart(bal.stats$bm.var.ratios[order(bal.stats$bm.var.ratios)],labels=bal.stats$covar,main="Variance Ratios",color="black",pch=19)
points(bal.stats$am.var.ratios[order(bal.stats$bm.var.ratios)],1:8,col="red",pch=19)
