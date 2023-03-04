ks.w<-function(x,y,xweights=rep(1,length(x)),yweights=rep(1,length(y)),
               alternative = c("two.sided", "less", "greater")){
  alternative <- match.arg(alternative)  
  ecdfw<-function(x,weights,val){
    ecdfw<-sum(weights[x<=val])/sum(weights)
    return(ecdfw)
  }#end ecdf
  x<-x[!is.na(x)]
  y<-y[!is.na(y)]
  z<-unique(sort(c(x,y)))  
  v<-1   # Initializes v
  for (i in 1:length(z))v[i]<-ecdfw(x,xweights,z[i])-ecdfw(y,yweights,z[i])  
  if(alternative=="greater")
    {
      ks.w<-abs(max(v))
    } else if(alternative=="less")
      {
        ks.w<-abs(min(v))
      } else {
        ks.w<-max(abs(v))
      }
  return (ks.w)
}


ks.boot.w  <- function(Tr, Co,
                       nboots=1000,
                       alternative = c("two.sided", "less", "greater"),
                       Tr.weights=rep(1,length(Tr)),
                       Co.weights=rep(1,length(Co)),
                       print.level=0)
  {
    alternative <- match.arg(alternative)
    tol <- sqrt(.Machine$double.eps)
    Tr <- Tr[!is.na(Tr) & !is.na(Tr.weights)]
    Co <- Co[!is.na(Co) & !is.na(Co.weights)]

    w    <- c(Tr, Co)
    w.weights <- c(Tr.weights,Co.weights)
    obs  <- length(w)
    n.x <- length(Tr)
    n.y <- length(Co)
    cutp <- n.x
    ks.boot.pval <- NULL
    bbcount <- 0

    if (nboots < 10)
      {
        nboots  <- 10
        warning("At least 10 'nboots' must be run; seting 'nboots' to 10")
      }

    if (nboots < 500)
      warning("For publication quality p-values it is recommended that 'nboots'\n be set equal to at least 500 (preferably 1000)") 
    
    fs.ks  <- ks.w(Tr, Co, xweights=Tr.weights, yweights=Co.weights,
                   alternative=alternative)        

    if (print.level > 0)
      cat("ks.boot:", alternative, "test\n")
        for (bb in 1:nboots)
          {
            if (print.level > 1)
              cat("s:", bb, "\n")
            
#            sindx  <- sample(1:obs, obs, replace=TRUE, prob=w.weights)
            sindx  <- sample(1:obs, obs, replace=TRUE)
            
            X1tmp <- w[sindx[1:cutp]]
            X2tmp <- w[sindx[(cutp+1):obs]]

            X1.w <- w.weights[sindx[1:cutp]]
            X2.w <- w.weights[sindx[(cutp+1):obs]]
            
            s.ks <- ks.w(X1tmp, X2tmp, xweights=X1.w, yweights=X2.w,
                         alternative=alternative)                    
            
            if (s.ks >= (fs.ks - tol) )
              bbcount  <- bbcount + 1
          }
    ks.boot.pval  <- bbcount/nboots
    
    ret  <- list(ks.boot.pvalue=ks.boot.pval, ks=fs.ks, nboots=nboots)
    class(ret)  <- "ks.boot"
    
    return(ret)
  } #end of ks.boot.w
