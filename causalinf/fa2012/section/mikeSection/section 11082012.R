par(mfrow = c(2,1))
plot(seq(0,14,.001),dgamma(seq(0,14,.001),3,1), 
main = "Density of a Gamma(3,1) R.V.", type = "l", ylim = c(0,.3), xlim = c(0,14), xlab = "", ylab = "")
median = replicate(10000,median(rgamma(1000,3,1)))
hist(median, freq = FALSE, breaks = 100, 
 main = "Density of the sample median of 
   1000 draws from a Gamma(3,1)")
lines(density(median))

par(mfrow = c(2,1))
sample = rgamma(1000,3,1)
hist(sample, breaks = 50, freq = FALSE)
curve(dgamma(x,3,1),add = TRUE)
bootstrapmedians = replicate(10000, median(sample(sample, replace = TRUE)))
hist(bootstrapmedians, freq = FALSE,breaks = 100, 
 main = "Histogram of medians from a bootstrap sample
 with overlay of density from 1000 draws from a 
 Gamma(3,1)")
lines(density(median))

sort(median)[c(251,9750)]
sort(bootstrapmedians)[c(251,9750)]

par(mfrow = c(2,1))
hist(median, freq = FALSE, breaks = 100, 
 main = "Density of the sample median of 
   1000 draws from a Gamma(3,1)", xlab = "95% of values within (2.550,2.802)", xlim = c(2.4,2.9))
lines(density(median))
hist(bootstrapmedians, breaks = 100, freq = FALSE,
  xlab = "95% bootstrap CI (2.560,2.799)", xlim = c(2.4,2.9))

#Some bootstrap examples.

#Generate data
x = log(abs(rnorm(1000,.5,1)))^2
hist(x, freq = FALSE, breaks = 100)

#Estimate sampling distribution of sample mean.
#Taking 10000 bootstrap samples.
bootsamps = replicate(10000, mean(sample(x, replace = TRUE)))
hist(bootsamps, freq = FALSE, breaks = 100)

#We can look at the mean of our sample and the 
#mean of our bootstrap samples
mean(x)
mean(bootsamps)

#We can get a sense for the standard error of the 
#sample mean by taking 
sd(bootsamps)

#We can get a 95% bootstrap confidence interval:
sort(bootsamps)[c(251,9750)]

#What if we had 1000 more data points 
y = rexp(1000,1/5)

#Do these histograms look the same?
par(mfrow = c(2,1))
hist(x, freq = FALSE, breaks = 100, xlim = c(0,40))
hist(y, freq = FALSE, breaks = 100, xlim = c(0,40))

#They look different.  Is this just chance?
#Compute the KS statistic
kstat = ks.test(x,y)$statistic

#Bootstrap to get pvalue
bootks = replicate(1000, ks.test(sample(c(x,y),1000,replace = TRUE), sample(c(x,y),1000,replace = TRUE))$statistic)

par(mfrow = c(1,1))
hist(bootks, freq = FALSE, xlim = c(0,.6))
abline(v = kstat, col= "red")

#Pvalue?
sum(bootks >= kstat)/1000
