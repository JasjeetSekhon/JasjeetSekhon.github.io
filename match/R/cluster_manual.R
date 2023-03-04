#
# Example of cluster GenMatch run run using ssh
#

# The situation: there are two machines: musil and quetelet and
# quetelet is a dual chip machine.

library(snow)
library(Matching)

#we could either create the connection here, outside of GenMatch, or
#let GenMatch do it.  But let's do it here so we don't have to keep
#setting it up.

cl <- makeCluster(c("musil","quetelet","quetelet"), type="SOCK")


data(lalonde)
attach(lalonde)

#The covariates we want to match on
X = cbind(age, educ, black, hisp, married, nodegr, u74, u75, re75, re74);

genout <- GenMatch(Tr=treat, X=X, cluster=cl)

#Manually stop the cluster because we manually set it up.  Remember
#that the ssh connections are still alive.
stopCluster(cl)
