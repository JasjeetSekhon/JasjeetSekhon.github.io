# Replication of GenMatch exp2 raw bias (non-cluster version).
# Bias for some regression adjustments is also presented
#
# internal note: compare with musil:~/xchg/Match/GenMatch/mc1/lalonde3reg1.R

#
# incorrect ps and bias adj model
# skewed data (based on the lalonde sample)
#

library(Matching)
data(lalonde)

sims <- 1000
estimand <- "ATE"
obs <- nrow(lalonde)

#setDATA seed
set.seed(2810192)

attach(lalonde)
propensity  <- glm(treat~ I(age^2) + I(educ^2) + black +
                   hisp + married + nodegr + I(re74^2) + I(re75^2) +
                   u74 + u75, family=binomial, data=lalonde)

X <- cbind(rep(1, length(treat)),
           propensity$linear.pred,
           I(log(age)^2),
           I(log(educ)^2),
           I(log(re74+0.01)^2),
           I(log(re75+0.01)^2))

propensity.coeffs <- propensity$coeff

#let's create some arbritrary weights
propensity.coeffs <- as.matrix(c(
                                 1+00,  #(Intercept)        
                                 .5,    #linear.pred
                                 .01, #age
                                 -.3,  #educ
                                 -0.01, #I(re74^2)          
                                 0.01   #I(re75^2)          
                                 ))

mu = X %*% propensity.coeffs
Tr.pred <- exp(mu)/(1+exp(mu))
#print(summary(Tr.pred))

TreatmentEffect <- 1000
TreatmentReal <- matrix(nrow=obs, ncol=1)

#return objects
raw <- matrix(nrow=sims)
reg1.est <- matrix(nrow=sims)
reg2.est <- matrix(nrow=sims)
reg3.est <- matrix(nrow=sims)
reg4.est <- matrix(nrow=sims)
reg1.se <- matrix(nrow=sims)
reg2.se <- matrix(nrow=sims)
reg3.se <- matrix(nrow=sims)
reg4.se <- matrix(nrow=sims)

results <- function(s)
  {
    cat("\n")
    
    cat("RAW bias:", mean(raw[1:s]), "\n")
    cat("RAW RMSE:", sqrt(mean(raw[1:s]^2)), "\n")    
    
    cat("\n")
    
    cat("reg1 bias:", mean(reg1.est[1:s]), "\n")
    cat("reg1 RMSE:", sqrt(mean(reg1.est[1:s]^2)), "\n")

    cat("\n")    

    cat("reg2 bias:", mean(reg2.est[1:s]), "\n")
    cat("reg2 RMSE:", sqrt(mean(reg2.est[1:s]^2)), "\n")

    cat("\n")

    cat("reg3 bias:", mean(reg3.est[1:s]), "\n")
    cat("reg3 RMSE:", sqrt(mean(reg3.est[1:s]^2)), "\n")

    cat("\n")

    cat("reg4 bias:", mean(reg4.est[1:s]), "\n")
    cat("reg4 RMSE:", sqrt(mean(reg4.est[1:s]^2)), "\n")

    cat("\n")        
  }

for (s in 1:sims)
  {
    cat("\nS:", s, "\n")
    
    for(i in 1:obs)
      {
        TreatmentReal[i] = sample(0:1, 1, prob=c(1-Tr.pred[i],Tr.pred[i]))
      }

    Y <- I(TreatmentEffect*TreatmentReal) + .1*exp(.7*log(re74+0.01) + .7*log(re75+0.01)) + rnorm(obs, 0, 10)
    #in hundreds of dollars
    Y <- Y

    raw[s] = TreatmentEffect - (mean(Y[TreatmentReal==1])-mean(Y[TreatmentReal==0])) 
    cat("raw[",s,"]:", raw[s], "\n")

    a1 <- lm(Y~TreatmentReal + age + educ + black + hisp + married +nodegr + re74 + re75 + u74 + u75)

    #Dehejia model, page 21 fn a
    a2 <- lm(Y~TreatmentReal + age + I(age^2) + educ + black + hisp + married +nodegr + re74 +re75)

    re74b <- re74
    re74b[re74==0] <- 0.01
    re75b <- re75
    re75b[re75==0] <- 0.01

    a3 <- lm(Y~TreatmentReal + age + educ + black + hisp + married +nodegr + I(log(re74b)) + I(log(re75b)))

    a4 <- lm(Y~TreatmentReal + age + I(age^2) + educ + black + hisp + married +nodegr + I(log(re74b)) + I(log(re75b)))    

    reg1.est[s] <- TreatmentEffect - a1$coef[2] 
    reg1.se[s] <- summary(a1)$coeff[2,2]

    reg2.est[s] <- TreatmentEffect - a2$coef[2] 
    reg2.se[s] <- summary(a2)$coeff[2,2]

    reg3.est[s]<- TreatmentEffect - a3$coef[2]
    reg3.se[s] <- summary(a3)$coeff[2,2]

    reg4.est[s]<- TreatmentEffect - a4$coef[2] 
    reg4.se[s] <- summary(a4)$coeff[2,2]            


    cat("reg1.est[",s,"]:",reg1.est[s],"\n")
    cat("reg2.est[",s,"]:",reg2.est[s],"\n")
    cat("reg3.est[",s,"]:",reg3.est[s],"\n")
    cat("reg4.est[",s,"]:",reg4.est[s],"\n")

  } #end of sims

results(s)
