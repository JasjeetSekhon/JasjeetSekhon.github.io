library(Matching)
options(width=150)

#use this if we have already downloaded the data
#dta <- read.table(file="election04.raw", header=TRUE)

#get the data from the web
dta <- read.table(file="http://sekhon.berkeley.edu/causalinf/R/election04.raw",
                  header=TRUE)

#
# For codebook see http://sekhon.berkeley.edu/causalinf/R/codebook_election04.raw.txt
#

cbind(as.data.frame(dta$county), dta$etouch)

glm1 <- glm(b04pc~etouch, data=dta)
summary(glm1)

glm2 <- glm(b04pc ~ etouch + b00pc + b00pc_sq + d96pc1 + v_change +
            income + hispanic +  b00pc_e + b00pcsq_e, data=dta)
summary(glm2)

glm3 <- glm(b04pc ~ etouch + hispanic, data=dta)
summary(glm3)

Y  <- dta$b04pc
Tr <- dta$etouch
X  <- dta$hispanic

rr <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT")
summary(rr)

MatchBalance(Tr~dta$hispanic, match.out=rr)

qqplot(dta$hispanic[Tr==1], dta$hispanic[Tr==0],
       ylim=c(0,.6), xlim=c(0,.6))
abline(0, 1, col=2)

qqplot(dta$hispanic[rr$index.treated], dta$hispanic[rr$index.control],
       ylim=c(0,.6), xlim=c(0,.6))
abline(0, 1, col=2)

mdta <- rbind(dta[rr$index.treated,], dta[rr$index.control,])

cbind(as.data.frame(mdta$county), mdta$etouch, mdta$hispanic)

cname.tr <- as.data.frame(dta$county[rr$index.treated])
cname.co <- as.data.frame(dta$county[rr$index.control])
treated <- dta$hispanic[rr$index.treated]
control <- dta$hispanic[rr$index.control]
dif <- dta$hispanic[rr$index.treated] - dta$hispanic[rr$index.control]

cbind(cname.tr, treated,  cname.co, control, dif)
