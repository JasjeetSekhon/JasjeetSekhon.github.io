########################################
# balance statistic plots
########################################


# this code comes from Mark Huberty, he combined Rocio's plotting code and Charlie's balance table code

# source("http://durkheim.berkeley.edu/code/balstatPlotv0.3.R")

#################################################
#################################################
## Balance Statistics Output with Plots
## Mark Huberty
## 24 October 2009
## v0.

## This function provides a one-line input to create plots of the balance statistics
## (p-values for the t-test for difference in means, and the ks statistic for difference
## in distributions) for the output of the MatchBalance function in the Matching() package for R. 

## This function combines code from the Balance Statistics functions written by Charlie Gibbons
## (http://cgibbons.berkeley.edu/Courses/PS236_F08/balanceTable.r) and the balance statistics
## plotting functions written by Rocio Titiunik (http://are.berkeley.edu/~rocio/R/graph.pval.public.R).
## The changes here are to make those two functions cooperate.

## BETA.
## v0.3 added error messaging for input values; added means and title handling
## v0.2 corrected for column mismatch that made means incorrect.

#################################################
#################################################

## Note: KS test statistics appear as 'NA' for dichotomous variables. Also,
## 'BM' is before matching and 'AM' is after matching.

## Input values
# covariates:    the list of variable names corresponding to the BalanceMatrix given to GenMatch.
#                This can either be the BalanceMatrix (data type "matrix") or a vector of covariate
#                names (data type "charatcter"). The latter is useful if BalanceMatrix contains
#                categorical vectors with n categories, which R converts to n dummy variables for
#                matching.

# mb.results:    the output of MatchBalance

# title:         optional string containing the title for the plot

# means:         logical value indicating whether to include the before- and after-matching
#                covariate means on the plot. This can be useful


## Output values
# 'results':     a matrix whose rows are different variables; whose first two columns contain the means
#                for treated and control; and whose remaining columns have the pvalues to be plotted
#                for every variable. In this function, "results" takes the output from the initial
#                balanceMat function and is used as the data object for plotting.



## Parameter values. These can be modified, but it's recommended to leave them alone.
# at1, at2,at3:  scalars which indicates where to locate the three differents groups
#                (mean treatment, mean controls, graph area) in the figure area

# xlim1 :        the left limit of the x-axis; right limit is always set to 1

# textsize:      scalar indicating the size of text in the figure

# legend:        logical indicating whether the legend should be included

# legendx:       scalar indicating the x-coordinate of the legend's location

# legendy:       scalar indicating the y-coordinate of the legend's location

# parcex:        scalar setting cex parameter


###############################################
## BEGIN PLOTTING CODE
###############################################

plot.pval <- function(mb.results, covariates, title=NULL, means=TRUE, legend=TRUE,
                      legendx=0.15,legendy=2.2, textsize=0.9, parcex=0.8, at1=-0.35,
                      at2=-0.15, at3=-0.9,xlim1=-0.85){


  ## Check to make sure "covariates" is either a character vector or a matrix
  if(!(is.character(covariates)) & !(is.matrix(covariates))){
    stop("Error: The covariates must be either a character vector or a matrix with variable names")
  }
    
                           
                                    
  balanceMat <- function(mb.results, covariates){

    ## Calculate the number of covariates
    ## If "covariates" is a matrix, take the number of columns; otherwise, take the number of
    ## elements in the vector
    n <- ifelse(class(covariates)=="matrix", dim(covariates)[2], length(covariates))

    ## Check to make sure that n is the same as the # of balance statistics
    bsCount <- length(mb.results$BeforeMatching)
    if(n != bsCount){stop("Error: the number of covariates does not match the size of the MatchBalance Results.")}

    ## Determine how the covariate names are provided, and then grab them
    if(class(covariates)=="matrix") rnames <- dimnames(covariates)[[2]]
    if(class(covariates)=="character") rnames <- covariates

    ## Construct the matrix of statistics from the MatchBalance data and attach it to the covariate names
    z <- t(sapply(1:n, function(x){
      c(rnames[x],
        round(mb.results$AfterMatching[[x]]$mean.Tr,3),
        round(mb.results$AfterMatching[[x]]$mean.Co,3),
        round(mb.results$BeforeMatching[[x]]$tt$p.value,2),
        round(mb.results$AfterMatching[[x]]$tt$p.value,2),
        round(mb.results$BeforeMatching[[x]]$tt$statistic,2),
        round(mb.results$AfterMatching[[x]]$tt$statistic,2),
        ifelse(is.null(mb.results$BeforeMatching[[x]]$ks$ks.boot.pvalue) ==
               0,round(mb.results$BeforeMatching[[x]]$ks$ks.boot.pvalue,2),
               NA),
        
        ifelse(is.null(mb.results$AfterMatching[[x]]$ks$ks.boot.pvalue) ==
               0, round(mb.results$AfterMatching[[x]]$ks$ks.boot.pvalue,2),
               NA))
    }))
    z <- as.data.frame(z)
    ##return(z)
    
    z[,2:9] <- apply(z[,2:9], 2, function(x){as.numeric(x)})
    mat <- z[,2:9]

    ## Apply the correct column names
    names(mat)<- c("Mean Tr.",
                   "Mean Con.",
                   "BM t p-value",
                   "AM t p-value",
                   "BM t stat",
                   "AM t stat",
                   "BM KS p-value",
                   "AM KS p-value")
    ## Apply the correct row names
    dimnames(mat)[[1]] <- z[,1]
    mat
  }
 

  ## Take the function above and apply it to the data supplied in the command
  results <- balanceMat(mb.results, covariates)
  
  ## set values of different parameters
  # pchset is the shape of the symbols
  # pchcolset is the color of the symbols
  xlim = c(xlim1,1); pchset = c(21,24,21,24); pchcolset = c("blue","blue", "red", "red")

  ## set margins and letter size
  par(cex=parcex, mai = c(0.5, 0.35, 1.1, 0.35))

  ## set number of rows to plot
  ny = nrow(results)

  ## create the empty figure
  if(!is.null(title))  plot(x=NULL,axes=F, xlim=xlim, ylim=c(1,ny),xlab="",ylab="", main=title)
  if(is.null(title))   plot(x=NULL,axes=F, xlim=xlim, ylim=c(1,ny),xlab="",ylab="")
  
  ## add the 0, 0.05 and 0.1 vertical lines
  abline(v=c(0,0.05,0.1),lty=c(1,4,4), lwd=c(1,2,2))
  axis(side=1,at=c(0,0.05,0.1,1),tick=TRUE, las=2, cex.axis=0.7)

  ## add labels on top of the three areas of the graph
  if(means==TRUE) axis(side=3,at=at1,labels="Mean\nTreated",tick=FALSE, padj=0.5,cex.axis=textsize)
  if(means==TRUE) axis(side=3,at=at2,labels="Mean\nControl",tick=FALSE, padj=0.5,cex.axis=textsize)
  axis(side=3,at=0.5,labels="P-values",tick=FALSE, padj=0.5,cex.axis=textsize)

  ## Fill the figure with the information which is inside the 'results' matrix
  ## Add the p-values of the t-statistics as points
  for(i in 3:4) points(results[,i],ny:1, pch = pchset[i-3+1], col = pchcolset[i-3+1], bg = pchcolset[i-3+1])
  ## Add the p-values of the ks statistics as points
  for(i in 7:8) points(results[,i],ny:1, pch = pchset[i-5+1], col = pchcolset[i-5+1], bg = pchcolset[i-5+1])

      
  ## Second, add each variable name and the means for treated and control
  for(i in 1:ny) {
    text(at3,ny-i+1,dimnames(results)[[1]][i],adj = 0,cex=textsize) # variable name
    if(means==TRUE) text(at1,ny-i+1,results[i,1], cex=textsize) # treatment mean
    if(means==TRUE) text(at2,ny-i+1,results[i,2], cex=textsize) # control mean
  }

  ## Add dotted horizontal lines every two variables to make it prettier
  for(i in seq(2,by=2,length.out=floor((ny-1)/2))) abline(h = i+0.5, lty = 3)

  ## Add legend
  # to move the legend change the x and y location
  if(legend) legend(x=-1,y=0.5,
                    c(colnames(results)[3:4],
                      colnames(results)[7:8]),
                    pch=pchset,
                    pt.bg =
                    pchcolset,
                    cex=0.8,
                    ncol=2,
                    xpd=NA
                    )

   return(results)
}


#########################################
## END PLOTTING CODE
#########################################

## USE
## plot.pval(mb.results, covariates, ...optional arguments)
