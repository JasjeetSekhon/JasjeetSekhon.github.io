library(Matching)

#dta <- read.table(file="election04.raw", header=TRUE)

dta <- read.table(file="http://sekhon.berkeley.edu/causalinf/R/election04.raw",
                  header=TRUE)

glm1 <- glm(b04pc~b00pc + b00pc_sq + d96pc1 + v_change + income + hispanic + etouch + b00pc_e + b00pcsq_e, data=dta)
summary(glm1)

glm2 <- glm(b04pc~etouch, data=dta)
summary(glm2)

rr <- Match(Y=dta$b04pc, Tr=dta$etouch, X=dta$b00pc, estimand="ATT")
summary(rr)
