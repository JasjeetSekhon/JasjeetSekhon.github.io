#multi processor GenMatch

#make sure that Matching is loaded
library("Matching")
library("snow")

GenMatch.orig <- GenMatch

GenMatch <- function(..., balance=TRUE, ncpus=NULL, port=NULL, options=defaultClusterOptions)
{
  cl <- makeLocalSOCKcluster(options=options, ncpus=ncpus, port=port)
#  Tr <- match.arg(Tr)
#  gm <- GenMatch.orig(Tr=Tr)
  gm <- GenMatch.orig(..., balance=balance, cluster=cl)
  stopCluster(cl)
  return(gm)
}


#
# Socket Implementation
#

#**** allow user to be different on different machines
#**** allow machinse to be selected from a hosts list
GENnewSOCKnode <- function(machine = "localhost", UseLocalShell=FALSE, ...,
                        options = defaultClusterOptions) {
    # **** allow some form of spec here
    # **** make sure options are quoted
    options <- addClusterOptions(options, list(...))
    port <- getClusterOption("port", options)
    scriptdir <- getClusterOption("scriptdir", options)
    if (getClusterOption("homogeneous")) {
        script <- file.path(scriptdir, "RSOCKnode.sh")
        rlibs <- paste(getClusterOption("rlibs", options), collapse = ":")
        rprog <- getClusterOption("rprog", options)
    }   
    else {
        script <- "RunSnowNode RSOCKnode.sh"
        rlibs <- NULL
        rprog <- NULL
    }
    rshcmd <- getClusterOption("rshcmd", options)
    user <- getClusterOption("user", options)
    env <- paste("MASTER=", getClusterOption("master", options),
                 " PORT=", port,
                 " OUT=", getClusterOption("outfile", options),
                 sep="")
    if (! is.null(rprog))
        env <- paste(env, " RPROG=", rprog, sep="")
    if (! is.null(rlibs))
        env <- paste(env, " R_LIBS=", rlibs, sep="")

    if(!UseLocalShell)
      {
        foo <- paste(rshcmd, "-l", user, machine, "env", env, script)
#        cat("rsh snow command:\n")
#        print(foo)
        system(foo)
      } else{
        foo <- paste("env", env, script)
        cat("local snow command:\n")
        print(foo)
        system(foo)        
      }
    
    con <- socketConnection(port = port, server=TRUE, blocking=TRUE,
                            open="a+b")
    structure(list(con = con, host = machine), class = "SOCKnode")
}

GENmakeSOCKcluster <- function(names, ..., options = defaultClusterOptions,
                               UseLocalShell=FALSE) {
    if (! exists("serialize") && ! require(serialize))
        stop("the `serialize' package is needed for SOCK clusters.")
    options <- addClusterOptions(options, list(...))
    cl <- vector("list",length(names))
    for (i in seq(along=cl))
        cl[[i]] <- GENnewSOCKnode(names[[i]], options = options,
                               UseLocalShell)
    class(cl) <- c("SOCKcluster")
    cl
}

makeLocalSOCKcluster <- function(ncpus=NULL, options=defaultClusterOptions, port=NULL, ...)
{
  #calculate the number of CPUs
  if (is.null(ncpus))
    {
      sysname <- Sys.info()[[1]]
      if (sysname=="Darwin" || sysname=="openbsd")
        {
          ncpus <- system("sysctl -n hw.ncpu", intern=TRUE)
          ncpus <- as.real(ncpus)
          name <- rep("localhost", ncpus)
        } else if (sysname=="Linux"){
          ncpus <- system("grep processor /proc/cpuinfo |wc |awk '{print $1}'", intern=TRUE)
          ncpus <- as.real(ncpus)
          name <- rep("localhost", ncpus)          
        } else if (sysname=="SunOS")
          {
            ncpus <- system("/usr/sbin/psrinfo | wc |awk '{print $1}'", intern=TRUE)
            ncpus <- as.real(ncpus)
            name <- rep("localhost", ncpus)                      
          } else {
          warning("Automatic CPU detection is not supported on this operating system.  Please set the number of processors 'cpus' manually.  I will continuing assuming two cpus.", immediate.=TRUE)
          name <- rep("localhost", 2)
        }
    } else {
      name <- rep("localhost", ncpus)
    }

  if (is.null(port))
    {
      port <- 60001+round(runif(1,0,10000))
    }
  cat("Using", ncpus, "local processors on port", port,"\n")          
  
  return(GENmakeSOCKcluster(names=name, options=options, UseLocalShell=TRUE, port=port,
                            ...))
}#end of makeLocalSOCKcluster
