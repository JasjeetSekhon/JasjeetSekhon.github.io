library(Matching)
example(ks.boot)

source("ksw.code2.R")

set.seed(3198)

mat <- cbind(lalonde$re75[rr$index.treated], lalonde$re75[rr$index.control],xweights=rr$weights)
ndata <- NULL
w <- 0
count <- 0
for(i in 1:nrow(mat))
  {
    cat("i:",i,"\n")
    count <- count+1
    w <- w+mat[i,3]
    cat("w:",w,"\n")

    if(w>(1-.000001))
      {
        if (count==1)
          {
            foo <- mat[i,]
          } else {
            cat("count:",count,"i:",i,"\n")
            foo <- apply(mat[(i-count+1):i,], 2, mean)
            foo[3] <- sum(mat[(i-count+1):i,3])
          }

        cat("rbind",i,"\n")
        cat(foo)
        ndata <- rbind(ndata,foo)
        w <- 0
        count <- 0
      }
  }

ks.test(ndata[,1], ndata[,2])
ks.w(ndata[,1], ndata[,2])
ks.boot(ndata[,1], ndata[,2],nboots=500)

ks.boot(lalonde$re75[rr$index.treated],
        lalonde$re75[rr$index.control],
        nboots=500)

ks.boot.w(lalonde$re75[rr$index.treated],
          lalonde$re75[rr$index.control],
          Tr.weights=rr$weights,
          Co.weights=rr$weights,
          nboots=500
          )

ks.w(ndata[,1], ndata[,2])

ks.w(ndata[,1], ndata[,2])
ks.test(ndata[,1], ndata[,2])

ks.w(ndata[,1], ndata[,2], alternative="less")
ks.test(ndata[,1], ndata[,2], alternative="less")

ks.w(ndata[,1], ndata[,2], alternative="greater")
ks.test(ndata[,1], ndata[,2], alternative="greater")


pvals <- c()
pvals.test <- c()
dstats <- c()
dstats.actual <- c()

for(i in 1:1000){
  print(i)
  
  tr <- rnorm(1000)
  
  con <- rnorm(1000)

  ks <- ks.test(tr,con)
  dstats.actual <- c(dstats.actual, ks$statistic)
  
  ksboot<-ks.boot.w(tr,con,nboots=500)

  kstest <- ks.test(tr,con)
  
  pvals<-c(pvals,ks$p.value)
  pvals.test <- c(pvals.test,kstest$p.value)    
  dstats<-c(dstats,ksboot$ks)
  
  cat("dstat:",mean(dstats),"\n")
  cat("dstat actual:",mean(dstats.actual),"\n")
  cat("test coverage:", sum(pvals.test<.05)/length(pvals.test),"\n")  
  cat("boot coverage:", sum(pvals<.05)/length(pvals), "\n")

  cat("\n")
}

#print(sum(pvals<.05)/length(pvals))

mean(dstats)

