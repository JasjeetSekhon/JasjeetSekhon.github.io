#set.seed(112312)

#
# hypothesis test code for the null hypothesis not equal to zero for all cases EXCEPT for case 4
#

obs  <- 100
sims <- 1000
x1 <- rnorm(obs, 1, 1)
x2 <- rnorm(obs, 0, 1)
#x2 <- rt(obs, df=0.1)

x <- c(x1,x2)

null <- 1
fs <- mean(x1)-mean(x2)
fs.std <- fs/sqrt(var(x1)+var(x2))
fs.std.null <- (fs-null)/sqrt(var(x1)+var(x2))
fs.ttest <- t.test(x1,x2,mu=null)
print(fs.ttest)
anaT <- 2*(1 - pnorm(abs(mean(x1) - mean(x2) - null)/sqrt(var(x1)/length(x1) + var(x2)/length(x2))))
cat("analytical:",anaT,"\n")

bs <- vector(length=sims, mode="numeric")

count1 <- 0
count2 <- 0
count3 <- 0
count4 <- 0

for (s in 1:sims)
  {
    indx1 <- sample(1:obs, replace=TRUE);
    sx1 <- x1[indx1]

    indx2 <- sample(1:obs, replace=TRUE);
    sx2 <- x2[indx2]

    bs[s] = mean(sx1)-mean(sx2)

    #based on the CI, probably the simplist way
    if (fs < null)
      {
        if(bs[s] >= null)
          count1 = count1+1
      } else {
        if(bs[s] < null)
          count1 = count1+1
      }

    #or this way, taking off the null
    if (fs < null)
      {
        if(bs[s]-fs < fs-null)
          count2 = count2+1
      } else {
        if(bs[s]-fs >= fs-null)
          count2 = count2+1
      }

    #studentized test.  Tests means
    sxx1 = sx1 -mean(x1) +mean(x)
    sxx2 = sx2 -mean(x2) + mean(x)
    sxx.sd <- (mean(sxx1)-mean(sxx2))/sqrt(var(sxx1)+var(sxx2))
    if (fs < null)
      {    
        if (sxx.sd < fs.std.null)
          count3 <- count3+1
      } else {
        if (sxx.sd >= fs.std.null)
          count3 <- count3+1
      }

    #or yet this way. 
    indx3 <- sample(1:(obs*2), replace=TRUE)
    sx <- x[indx3]
    sxx1 <- sx[1:obs]
    sxx2 <- sx[(obs+1):(obs*2)]
    if (fs < 0)
      {    
        if ( (mean(sxx1)-mean(sxx2)) < fs )
          count4 = count4+1
      } else {
        if ( (mean(sxx1)-mean(sxx2)) >= fs )
          count4 = count4+1
      }
  }

cat("Based on the CI.  Probably the simplist way.. ASL 1:",(count1/sims)*2,"\n")

cat("Take off the null (the full sample value).... ASL 2:",(count2/sims)*2,"\n")

cat("studentized test version.  Testing the means. ASL 3:",(count3/sims)*2,"\n")

cat("This is very close to a PERMUTATION test..... ASL 4:",(count4/sims)*2,"\n")


#cat("\n\ndrawing histogram...\n")
#hist(bs, main="Histogram of Boostrap Samples", xlab="mean(x_{b})-mean(z_{b})")
#abline(v=0, col="red")
