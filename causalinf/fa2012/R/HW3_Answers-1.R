rm(list=ls())
load(file=url("http://sekhon.berkeley.edu/causalinf/data/hw2data.RData"))

#Problem 4
##Part a
reg.model = lm(vote.pop~treat+as.factor(block), data = data)
reg.estimate = coef(reg.model)[2]

diff.means = function(y,treat,blocks=NA){
  if(length(unique(blocks)==1)){
    ave.treated = mean(y[treat=="client"])
    ave.control = mean(y[treat=="pub.pol"])
    test.stat = ave.treated-ave.control
  }
else{
  #use tapply to calculate the within block averages
  ave.treated = tapply(y[treat=="client"],blocks[treat=="client"],mean)
  ave.control = tapply(y[treat=="pub.pol"],blocks[treat=="pub.pol"],mean)
  wt = tapply(y,blocks,length)/length(y)
  test.stat = sum(wt*(ave.treated-ave.control))
}
  return(test.stat)
}

itt.estimate = diff.means(data$vote.pop,data$treat,data$block)

#Problem 4
##Part b
#Let's create a treatment assignment function
treat.assign = function(treat,blocks=NA){
  if(length(unique(blocks))==1){
    treat.vector = sample(treat)
  }
  else{
    #randomize within blocks using tapply
    treat.vector = tapply(treat,blocks,sample)
  #tapply returns a list, to turn into a vector, use "unlist"
  treat.vector = unlist(treat.vector)
  }
  return(treat.vector)
}

omega = replicate(5000,treat.assign(data$treat,data$block))
#We only want unique treatment assignments, so let's get rid of duplicates. Specifying "margin=2" keeps unique columns
omega = unique(omega,MARGIN=2)


#what's the true test stat?
true.test.stat = diff.means(data$vote.pop,data$treat,data$block)

#Now let's compute the entire distribution of test-statistics under the null hypothesis
#create a matrix to hold the test statistics
test.stat.dist = matrix(nrow=ncol(omega))
#loop over omega, caculate the test stat every iteration
for (i in 1:ncol(omega)){
  treat.fake = omega[,i]
  fake.test.stat = diff.means(data$vote.pop,treat.fake,data$block)
  test.stat.dist[i] = fake.test.stat
}

#Let's plot the randomization distribution
plot(density(test.stat.dist),col="blue", main="Randomization Distribution")
#where does the true test statistic fall?
abline(v=true.test.stat,col="red",lwd=2)


#what's the p-value?
 2* min(sum(test.stat.dist>=true.test.stat)/ncol(omega),1-(sum(test.stat.dist>=true.test.stat)/ncol(omega)))

#create a rank sum function
strat.ranksum = function(y,treat,blocks){
  #rank within strata
  ranks = unlist(tapply(y,blocks,rank))
  #sum ranks of treated units
  ranksum = sum(ranks[treat=="client"])
  return(ranksum)
}

#observed test statistic?
true.test.stat = strat.ranksum(data$vote.pop, data$treat, data$block)

#Now let's compute the entire distribution of test-statistics under the null hypothesis
#create a matrix to hold the test statistics
test.stat.dist = matrix(ncol=ncol(omega))
#loop over omega, caculate the test stat every iteration
for (i in 1:ncol(omega)){
  treat.fake = omega[,i]
  fake.test.stat = strat.ranksum(data$vote.pop, treat.fake, data$block)
  test.stat.dist[i] = fake.test.stat
}

#Let's plot the randomization distribution
plot(density(test.stat.dist),col="blue", main="Randomization Distribution")
#where does the true test statistic fall?
abline(v=true.test.stat,col="red",lwd=2)


#what's the p-value?
 2* min(sum(test.stat.dist>=true.test.stat)/ncol(omega),1-(sum(test.stat.dist>=true.test.stat)/ncol(omega)))

       
#Problem 4
#Part c

omega = replicate(5000,treat.assign(data$treat,data$block))
#We only want unique treatment assignments, so let's get rid of duplicates. Specifying "margin=2" keeps unique columns
omega = unique(omega,MARGIN=2)


val.noreject = c()
for(j in seq(-1,1,by=.01)){
  #create new outcome vector with hypothesis subtracted out from treated units
  vote.pop.new = ifelse(data$treat=="client", data$vote.pop-j, data$vote.pop) 
                                        #observed test statistic?
  true.test.stat = strat.ranksum(vote.pop.new, data$treat, data$block)
                                        #Now let's compute the entire distribution of test-statistics under the null hypothesis
                                        #create a matrix to hold the test statistics
  test.stat.dist = matrix(ncol=ncol(omega))
                                        #loop over omega, caculate the test stat every iteration
  for (i in 1:ncol(omega)){
    treat.fake = omega[,i]
    fake.test.stat = strat.ranksum(vote.pop.new, treat.fake, data$block)
    test.stat.dist[i] = fake.test.stat
  }

                                        #what's the p-value?
  p.val = 2* min(sum(test.stat.dist>=true.test.stat)/ncol(omega),1-(sum(test.stat.dist>=true.test.stat)/ncol(omega)))
  #if the null hypothesis is not rejected, keep value
  if(p.val>.05) val.noreject = c(val.noreject, j)
  print(j)
}

#The 95% confidence interval is
range(val.noreject)
  
                     
# Problem 4
# Part d
     
reg.model = lm(vote.pop~treat+as.factor(block), data = data) 
y = lm(vote.pop~reg.voters+as.factor(block),data=data)$residuals   

#whats the true or observed test statistic?
true.test.stat = strat.ranksum(y,data$treat,data$block) 

#Now let us compute the entire distribution of test-statistics 
 #under the null hypothesis
#create a matrix to hold the test statistics
test.stat.dist = matrix()
#loop over omega, caculate the test stat every iteration
# reuse omega from above
for (i in 1:ncol(omega)){
  treat.fake = omega[,i]
  fake.test.stat = strat.ranksum(y,treat.fake,data$block)
  test.stat.dist[i] = fake.test.stat
}      

#Lets plot the randomization distribution
plot(density(test.stat.dist),col="blue", main="Randomization Distribution")
#where does the true test statistic fall?
abline(v=true.test.stat,col="red",lwd=2)

#whats the p-value?
1-sum(test.stat.dist>true.test.stat)/ncol(omega)   



       


