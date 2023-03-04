obs  <- 100
sims <- 500
mc   <- 200

rmean1 <- 0
rvar1  <- 1
rmean2 <- 0
rvar2  <- 1

bs <- vector(length=sims, mode="numeric")
mc.results <- matrix(0, nrow=mc, ncol=5)

for (m in 1:mc)
  {

    x1 <- rnorm(obs, rmean1, rvar1)
                                        #x2 <- rnorm(obs, .3, 1)
    x2 <- rnorm(obs, rmean2, rvar2)
    x <- c(x1,x2)

    fs <- mean(x1)-mean(x2)
    fs.std <- fs/sqrt(var(x1)+var(x2))
    fs.ttest <- t.test(x1,x2)

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
      } #sims
    mc.results[m,1] = (count1/sims)/.5
    mc.results[m,2] = (count2/sims)/.5
    mc.results[m,3] = (count3/sims)/.5
    mc.results[m,4] = (count4/sims)/.5
    mc.results[m,5] = t.test(x1, x2)$p.val
  }#close mc

print(mc.results)

cat("\n")
alpha = .1
cat(alpha,"NOMINAL COVERAGE RESULTS:\n")
cat("case 1:", mean(mc.results[,1] < .1), "\n")
cat("case 2:", mean(mc.results[,2] < .1), "\n")
cat("case 3:", mean(mc.results[,3] < .1), "\n")
cat("case 4:", mean(mc.results[,4] < .1), "\n")
cat("t-test:", mean(mc.results[,5] < .1), "\n")
