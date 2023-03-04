#case 4 now does the KS test

obs  <- 100
sims <- 1000
x1 <- rnorm(obs, 0, 1)
x2 <- rnorm(obs, 0, 1)
#x2 <- rt(obs, df=0.1)

x <- c(x1,x2)

fs <- mean(x1)-mean(x2)
fs.std <- fs/sqrt(var(x1)+var(x2))
fs.ttest <- t.test(x1,x2)
print(fs.ttest)

fs.ks <- ks.test(x1, x2)
print(fs.ks)

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
    if (fs < 0)
      {
        if(bs[s] >= 0)
          count1 = count1+1
      } else {
        if(bs[s] < 0)
          count1 = count1+1
      }

    #or this way, taking off the null
    if (fs < 0)
      {
        if(bs[s]-fs < fs)
          count2 = count2+1
      } else {
        if(bs[s]-fs >= fs)
          count2 = count2+1
      }

    #studentized test.  Tests means
    sxx1 = sx1 -mean(x1) +mean(x)
    sxx2 = sx2 -mean(x2) + mean(x)
    sxx.sd <- (mean(sxx1)-mean(sxx2))/sqrt(var(sxx1)+var(sxx2))
    if (fs < 0)
      {    
        if (sxx.sd < fs.std)
          count3 <- count3+1
      } else {
        if (sxx.sd >= fs.std)
          count3 <- count3+1
      }

    #or yet this way. Tests if F=G (the entire distrubtion)
    indx3 <- sample(1:(obs*2), replace=TRUE)
    sx <- x[indx3]
    sxx1 <- sx[1:obs]
    sxx2 <- sx[(obs+1):(obs*2)]
    #this will always be true for the KS D stat > 0
    if (ks.test(sxx1, sxx2)$stat >= fs.ks$stat) {
          count4 = count4+1
      }
  }

cat("Based on the CI.  Probably the simplist way.. ASL 1:",(count1/sims*2),"\n")

cat("Take off the null (the full sample value).... ASL 2:",(count2/sims*2),"\n")

cat("studentized test version.  Testing the means. ASL 3:",(count3/sims*2),"\n")

cat("bootstrap ks test ........................... ASL 4:",(count4/sims),"\n")


cat("\n\ndrawing histogram...\n")
hist(bs, main="Histogram of Boostrap Samples", xlab="mean(x_{b})-mean(z_{b})")
abline(v=0, col="red")
