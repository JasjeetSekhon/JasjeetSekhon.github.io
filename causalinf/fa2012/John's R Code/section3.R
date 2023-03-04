###############################################
# Section 3: Permutation Inference
#
# John Henderson         
#
# PS 236A / Stat 239A Fall 2012
###############################################

rm(list=ls(all=TRUE))
#Change path to correct directory


load(file=url("http://sekhon.berkeley.edu/causalinf/data/section3.RData"))
data <- data[order(data$block),]

y <- data$complaints.dummy
treat <- data$tr.complaints
blocks <- data$block


#How many possible treatment assignments?
choose(2,1)^18 * choose(4,2)^12

#Lets create a treatment assignment function
treat.assign <- function(treat,blocks=NA){
  if(length(unique(blocks))==1){
    treat.vector <- sample(treat)
  }
  else{
    #randomize within blocks using tapply
    treat.vector <- tapply(treat,blocks,sample)
  #tapply returns a list, to turn into a vector, use "unlist"
  treat.vector <- unlist(treat.vector)
  }
  return(treat.vector)
}
 
set.seed(1005)
omega <- replicate(5000,treat.assign(treat,blocks))
#We only want unique treatment assignments, so let us get rid of duplicates. Specifying "margin=2" keeps unique columns
omega <- unique(omega,MARGIN=2)


#First let us check the results using standard methods
summary(glm(y~treat + blocks, family=binomial(link=logit)))
summary(lm(y~treat + blocks))    

#Randomization inference
# First let us create a function that computes the difference-in-means; then 
#  let us do an extension of Fisher's exact test with blocks 
#  (taking into account the block structure of the experiment)

# 1. Mean Differences 
diff.means <- function(y,treat,blocks){
  #use tapply to calculate the within block averages
  ave.treated <- tapply(y[treat==1],blocks[treat==1],mean)
  ave.control <- tapply(y[treat==0],blocks[treat==0],mean)
  wt <- tapply(y,blocks,length)/length(y)
  test.stat <- sum(wt*(ave.treated-ave.control))
  return(test.stat)
}

#what's the true test stat?
true.test.stat <- diff.means(y,treat,blocks)

#Now let us compute the entire distribution of test-statistics under the null hypothesis
#create a matrix to hold the test statistics
test.stat.dist <- matrix(nrow=ncol(omega))
#loop over omega, caculate the test stat every iteration
for (i in 1:ncol(omega)){
  treat.fake <- omega[,i]
  fake.test.stat <- diff.means(y,treat.fake,blocks)
  test.stat.dist[i] <- fake.test.stat
}

#Let's plot the randomization distribution
plot(density(test.stat.dist),col="blue", main="Randomization Distribution")
#where does the true test statistic fall?
abline(v=true.test.stat,col="red",lwd=2)


#what's the p-value?
sum(test.stat.dist>=true.test.stat)/ncol(omega)
            

# 2. Mantel-Haenszel 
                    
# 2x2xS Mantel-Haenszel statistic => analogous to Fisher's exact test with S>1
mantel.stat <- function(y,treat,blocks){
  #use tapply to calculate the within block averages
  sum.treated <- tapply(y[treat==1],blocks[treat==1],sum)
  test.stat <- sum(sum.treated)
  return(test.stat)
}

#what's the true test stat?
true.test.stat <- mantel.stat(y,treat,blocks)   
   
#Now let us compute the entire distribution of test-statistics under the null hypothesis
#create a matrix to hold the test statistics
test.stat.dist <- matrix(nrow=ncol(omega))
#loop over omega, caculate the test stat every iteration
for (i in 1:ncol(omega)){
  treat.fake <- omega[,i]
  fake.test.stat <- mantel.stat(y,treat.fake,blocks)
  test.stat.dist[i] <- fake.test.stat
}

#Let's plot the randomization distribution
plot(density(test.stat.dist),col="blue", main="Randomization Distribution")
#where does the true test statistic fall?
abline(v=true.test.stat,col="red",lwd=2)

sum(test.stat.dist>=true.test.stat)/ncol(omega)                


# 3. Wilcoxon Rank Sum Statistic

####
#Let's look at turnout now
#####
y <-(data$turnout08)

#create a rank sum function
strat.ranksum <- function(y,treat,blocks){
  #rank within strata
  ranks <- unlist(tapply(y,blocks,rank))
  #sum ranks of treated units
  ranksum <- sum(ranks[treat==1])
  return(ranksum)
}

#what's the true or observed test statistic?
true.test.stat <- strat.ranksum(y,treat,blocks)

#Now let us compute the entire distribution of test-statistics under the null hypothesis
#create a matrix to hold the test statistics
test.stat.dist <- matrix(ncol=ncol(omega))
#loop over omega, caculate the test stat every iteration
for (i in 1:ncol(omega)){
  treat.fake <- omega[,i]
  fake.test.stat <- strat.ranksum(y,treat.fake,blocks)
  test.stat.dist[i] <- fake.test.stat
}

#Let's plot the randomization distribution
plot(density(test.stat.dist),col="blue", main="Randomization Distribution")
#where does the true test statistic fall?
abline(v=true.test.stat,col="red",lwd=2)

#what's the p-value? [Note one-sided; just care about rarity in one direction]
1-sum(test.stat.dist>true.test.stat)/ncol(omega)

 
# 4. Permutation with Covariate Adjustment

#what's the dispersion in Y before adjustment
sd(y)

#grab the residuals
y <- lm(turnout08~turnout2006+urban+num.reg.voters+pctPresParty2006+as.factor(block),data=data)$residuals
#what's the dispersion before adjustment?
sd(y)


#what's the true or observed test statistic?
true.test.stat <- strat.ranksum(y,treat,blocks)

#Now let us compute the entire distribution of test-statistics under the null hypothesis
#create a matrix to hold the test statistics
test.stat.dist <- matrix()
#loop over omega, caculate the test stat every iteration
for (i in 1:ncol(omega)){
  treat.fake <- omega[,i]
  fake.test.stat <- strat.ranksum(y,treat.fake,blocks)
  test.stat.dist[i] <- fake.test.stat
}

#Let's plot the randomization distribution
plot(density(test.stat.dist),col="blue", main="Randomization Distribution")
#where does the true test statistic fall?
abline(v=true.test.stat,col="red",lwd=2)

#what's the p-value?
1-sum(test.stat.dist>true.test.stat)/ncol(omega)   

# END section3.R