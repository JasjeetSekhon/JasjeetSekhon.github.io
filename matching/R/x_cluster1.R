#
# Example of cluster GenMatch run run using ssh tunnel.

# The situation: there are four dual chips machines: x1, x2, x3, x4.
# We logged into x1 and want to use all four of them.

# First, execute the following UNIX command on the same machine which
# will run this file---i.e., in this example x1.  These commands only
# have to be run once.  They will stay alive until killed by someone
# or the system is rebooted.

# ssh -f -N -o keepalive=yes -R 10187:localhost:10187 x2
# ssh -f -N -o keepalive=yes -R 10187:localhost:10187 x3
# ssh -f -N -o keepalive=yes -R 10187:localhost:10187 x4

#Note that we will not run the command for x1 because we assume that
#we are on x1.  If this assumption is false, we need to execute this
#command for x1 also.  Here is a description of the options:

# -f: background the ssh job
# -N: do not execute a command, just setup the ssh tunnel
# -o keepalive=yes: keep this connection alive even when inactive
# -R 10187:localhost:10187: map local port 10187 to the port 10187 on .
#                           the remote machine.  This is the default
#                           port used by the cluster.

library(snow)
library(rgenoud, lib="~/lib/R")
library(Matching, lib="~/lib/R")

#We set the master to localhost because of the ssh tunnel described
#above.  The port has to be set equal to the port which was tunneled
#in the ssh command.

setDefaultClusterOptions(master="localhost", port=10187)

#We could either create the connection here, outside of GenMatch, or
#let GenMatch do it.  But let's do it here so we don't have to keep
#setting it up. It is recommended that you setup ssh so it doesn't require a password to
#login for more instruction on this see:
#http://atoms.vuse.vanderbilt.edu/Elements/CompMeth/Comp_Methods_ssh2.html

cl <- makeSOCKcluster(c("x1","x1","x2","x2","x3","x3","x4","x4"))

data(lalonde)
attach(lalonde)

#The covariates we want to match on
X = cbind(age, educ, black, hisp, married, nodegr, u74, u75, re75, re74);

genout <- GenMatch(Tr=treat, X=X, estimand="ATE",
                   pop.size=50, max.generations=3, wait.generations=2,
                   hard.generation.limit=TRUE,
                   cluster=cl)
#Total run time : 0 hours 0 minutes and 5 seconds

genout <- GenMatch(Tr=treat, X=X, estimand="ATE",
                   pop.size=50, max.generations=3, wait.generations=2,
                   hard.generation.limit=TRUE)
#Total run time : 0 hours 0 minutes and 33 seconds


#Manually stop the cluster because we manually set it up.  Remember
#that the ssh connections are still alive.

stopCluster(cl)
