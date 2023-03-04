#
# Example of cluster GenMatch run run using ssh
#

# The situation: there are two machines: musil and quetelet and
# quetelet is a dual chip machine.

library(parallel) 
library(Matching)

#we could either create the connection here, outside of GenMatch, or
#let GenMatch do it.  But let's do it here so we don't have to keep
#setting it up.

cl <- makeCluster(c("musil","quetelet","quetelet"), type="SOCK")

#This is a silly example because the "sin" function is so quick to evaluate
genout <- genoud(sin, 1, cluster=cl)

#Manually stop the cluster because we manually set it up.  Remember
#that the ssh connections are still alive.
stopCluster(cl)
