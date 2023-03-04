library(Matching)

load(file="a2.Rdata")

set.seed(23)
#set.seed(3827287)
#set.seed(372719)

indx <- 1:185

nrow <- nrow(foo)

indx2 <- sample(186:nrow, 1000)
#indx2 <- 186:nrow

#full
#Estimate...  2401.7 
#AI SE......  910.37 
#T-stat.....  2.6382 
#p.val......  0.0083353 


dta <- foo[c(indx,indx2),]
dta$u74 <- as.real(dta$re74==0)
dta$u75 <- as.real(dta$re75==0)

     glm1  <- glm(treat~age + I(age^2) + education + I(education^2) + black +
                  hispan + married + nodegree + re74  + I(re74^2) + re75 + I(re75^2) +
                  u74 + u75, family=binomial, data=dta)

     #
     #save data objects
     #
     X  <- glm1$fitted
     Y  <- dta$re78
     Tr  <- dta$treat

     #
     # one-to-one matching with replacement (the "M=1" option).
     # Estimating the treatment effect on the treated (the "estimand" option defaults to ATT).
     #
     rr  <- Match(Y=Y, Tr=Tr, X=X, M=1);
     summary(rr)

     #
     # Let's check for balance
     # 'nboots' and 'nmc' are set to small values in the interest of speed.
     # Please increase to at least 500 each for publication quality p-values.  
     mb  <- MatchBalance(treat~age + I(age^2) + education + I(education^2) + black +
                         hispan + married + nodegree + re74  + I(re74^2) + re75 + I(re75^2) +
                         u74 + u75, data=dta, match.out=rr, nboots=500)

print(names(dta))
lalonde2 <- dta
save(lalonde2, file="lalonde2b.RData")
