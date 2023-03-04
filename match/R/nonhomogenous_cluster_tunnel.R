#
# Example of cluster GenMatch run run using ssh tunnel with non-homogenous layout.

# The situation: there are five dual chips machines: lapo, x1, x2, x3, x4.  We
# logged into lapo and want to use all five of them.  Also, lapo has a
# different layout than the other machines.

#First, because we have a nonhomogenous layout we need to add the following
#two lines to our .cshrc file on all of the machines (if we are using another
#shell, we need to add the lines to the approriate shell).
#
#set path=($path <path to snow library>)
#setenv R_SNOW_LIB <path to library directory in which snow is>
#
#for lapo these two lines are:
#set path=($path /usr/local/lib/R/library/snow)
#setenv R_SNOW_LIB /usr/local/lib/R/library
#
#And for the x machines the lines are:
#set path=($path /usr/lib64/R/library/snow)
#setenv R_SNOW_LIB /usr/lib64/R/library
#
#For additional details see:
#http://www.stat.uiowa.edu/~luke/R/cluster/cluster.html#Section:InhomogeneousSystems

# Second, execute the following UNIX command on the same machine which
# will run this file---i.e., in this example lapo.  These commands only
# have to be run once.  They will stay alive until killed by someone
# or the system is rebooted.

# ssh -f -N -o keepalive=yes -R 10187:localhost:10187 x1
# ssh -f -N -o keepalive=yes -R 10187:localhost:10187 x2
# ssh -f -N -o keepalive=yes -R 10187:localhost:10187 x3
# ssh -f -N -o keepalive=yes -R 10187:localhost:10187 x4

#Note that we will not run the command for lapo because we assume that
#we are on lapo.  If this assumption is false, we need to execute this
#command for lapo also.  Here is a description of the options:

# -f: background the ssh job
# -N: do not execute a command, just setup the ssh tunnel
# -o keepalive=yes: keep this connection alive even when inactive
# -R 10187:localhost:10187: map local port 10187 to the port 10187 on .
#                           the remote machine.  This is the default
#                           port used by the cluster.

library(snow)
library(Matching)

#We set the master to localhost because of the ssh tunnel described
#above.  The port has to be set equal to the port which was tunneled
#in the ssh command.

setDefaultClusterOptions(master="localhost", port=10187)

#we could either create the connection here, outside of GenMatch, or
#let GenMatch do it.  But let's do it here so we don't have to keep
#setting it up.

cl <- makeSOCKcluster(c("lapo","lapo",
                        "x1.hmdc.harvard.edu","x1.hmdc.harvard.edu",
                        "x2.hmdc.harvard.edu","x2.hmdc.harvard.edu",
                        "x3.hmdc.harvard.edu","x3.hmdc.harvard.edu",
                        "x4.hmdc.harvard.edu","x4.hmdc.harvard.edu"))

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
