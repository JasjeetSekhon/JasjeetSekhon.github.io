## Code sample for setting up a .R script to invoke the right
## routines based on whether the code is run on the cluster or not
## 6 February 2010

## This code consists of two functions that should both be called.
## clusterCreate sets up the cluster
## clsuterShutdown  shuts the cluster down and stops MPI if necessary
## They operate as a pair; clusterCreate should be called at the start of
## parallelized code, and clusterShutdown at the end.

## NOTE: This definitely works on Linux and OSX. It should work on
## Windows even if Sys.info() isn't enabled. Please report if that doesn't happen.

## The basic specification:
## Load the basic clustering libraries for using SNOW and MPI clusters
## Check whether the system name matches the pscluster name
## If so, use the Rmpi routine to discover the cluster invoked via salloc / orterun
## Else set up a socket cluster for a dual-core CPU
## Note that in either case, the cluster object is named cl, and hence is portable
## <code to run on the cluster cl>
## When done, shut the cluster down
## If on the pscluster, run quit MPI for a clean shutdown



## Load some libraries
library(snow)


##########################################
## Cluster setup
##########################################

## Check whether you are on the pscluster or your local machine
## The cluster is identified by the names of the nodes
## All node names are contained here

clusterCreate <- function(){
  thisnode <- as.list(Sys.info())$nodename
  cl.nodenames <- c("pscluster.polisci.berkeley.edu",
                    "pscluster",
                    "n0000",	
                    "n0001",
                    "n0002",
                    "n0003",
                    "n0004",
                    "n0005",
                    "n0006",
                    "n0007",
                    "n0008",
                    "n0009",
                    "n0010",
                    "n0011",
                    "n0012"
                    )
  oncluster <- thisnode %in% cl.nodenames
  rm(thisnode, cl.nodenames)
  
  if(oncluster==TRUE){
    ## Then I'm on pscluster
    library(Rmpi)
    
    ## Run the Rmpi routines to discover the cluster and invoke it
    mpirank <- mpi.comm.rank(0)
    
    if(mpirank==0){ ## Is this the master node
      cat("Launching master, mpirank=",mpirank,"\n")
      makeMPIcluster()
      
    }else{ ## or a worker node                                                                                                                
      
      cat("Launching slave with, mpirank=",mpirank,"\n")
      sink(file="/dev/null")
      slaveLoop(makeMPImaster());mpi.finalize();q()
      
    }
    
    
    
    ## Grab the cluster that was just created, and assign 
    ## it to an object for later use
    ## Notice the global variable definition
    cl <<- getMPIcluster()
    
  }else{
    ## I guess I'm on my local machine
    
    ## Set up a cluster that uses both cores in a dual-core CPU. If you have
    ## a different number of cores (1, 3, 4, etc), use that number of repetitions
    ## of "localhost". You can also see Jas Sekhon's NCPUS code.
    
    setDefaultClusterOptions(master="localhost")
    cl <<- makeCluster(c("localhost","localhost"),
                       type="SOCK"
                       )
    
  }
}





clusterShutdown <- function(cluster=cl){

 
  thisnode <- as.list(Sys.info())$nodename
  cl.nodenames <- c("pscluster.polisci.berkeley.edu",
                    "pscluster",
                    "n0000",	
                    "n0001",
                    "n0002",
                    "n0003",
                    "n0004",
                    "n0005",
                    "n0006",
                    "n0007",
                    "n0008",
                    "n0009",
                    "n0010",
                    "n0011",
                    "n0012"
                    )
  oncluster <- thisnode %in% cl.nodenames
  rm(thisnode, cl.nodenames)
  
#############################################
  ## Cluster shutdown
#############################################
  ## Shut the cluster down
  stopCluster(cl)
  
  ## If you are on the pscluster, stop MPI for a clean shutdown
  if(oncluster==TRUE){
    
    mpi.quit()
    
  }
}

## Usage
## \source(clusterSetup.R)
## library(snow)
## clusterCreate()
## Your code, running against a cluster called "cl"
## clusterShutdown()
