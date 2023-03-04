rm(list=ls(all=TRUE))
#Change path to correct directory
load("/Users/ElGuapo/Dropbox/projects/ps236/section/section3/section3.RData")
data <- data[order(data$block),]

y <- data$complaints.dummy
treat <- data$tr.complaints
blocks <- data$block



#How many possible treatment assignments?
choose(2,1)^18 * choose(4,2)^12

#Let's create a treatment assignment function
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

omega <- replicate(5000,treat.assign(treat,blocks))
#We only want unique treatment assignments, so let's get rid of duplicates. Specifying "margin=2" keeps unique columns
omega <- unique(omega,MARGIN=2)


#First let's check the results using standard methods
summary(glm(y~treat + blocks, family=binomial))

#Randomization inference
##First let's create a function that computes the difference-in-means
##(taking into account the block structure of the experiment)
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

#Now let's compute the entire distribution of test-statistics under the null hypothesis
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

#Now let's compute the entire distribution of test-statistics under the null hypothesis
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

#what's the p-value?
sum(test.stat.dist<=true.test.stat)/ncol(omega)


#####
###With Covariance adjustment
#####


#what's the dispersion in Y before adjustment
sd(y)

#grab the residuals
y <- lm(turnout08~turnout2006+urban+num.reg.voters+pctPresParty2006+block,data=data)$residuals
#what's the dispersion before adjustment?
sd(y)


#what's the true or observed test statistic?
true.test.stat <- strat.ranksum(y,treat,blocks)

#Now let's compute the entire distribution of test-statistics under the null hypothesis
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
sum(test.stat.dist<=true.test.stat)/ncol(omega)

