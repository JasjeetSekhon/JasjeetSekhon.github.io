### Creating Tables of Balance Statistics ###
### Political Science 236, UC Berkeley    ###
### Fall 2008                             ###

## 'balanceTable' creates LaTeX tabular output that can be placed within a table
##   environment.
## 'balanceMatrix' creates a matrix of balance statistics.

## When performing 'MatchBalance()', you must supply a formula to examine
##   balance. For the functions below, 'bal.out' must be the output from a
##   "MatchBalance' with a formula that takes the form treatment ~ 'covariates',
##   where 'covariates' is the other object that enters the functions below.


## The following creates inputs for a tabular object in
##   LaTeX. Column 1 is the covariate name, then there are before and after
##   matching columns, three in each. The first statistic is the t-statistic,
##   then the t-test p-value, then the ks bootstrap value. 'covariates' is your
##   balance matrix and 'bal.out' is the 'MatchBalance()' results.
balanceTable <- function(covariates, bal.out){
  cat("\\begin{center} \\begin{tabular}{lrrrcrrr} \\hline \\hline \n")
  cat("Covariate	&\\multicolumn{3}{c}{Before matching}",
    "&&\\multicolumn{3}{c}{After matching}", "\\", "\\", "\\cline{2-4} \\cline{6-8} \n",
    sep="")
  cat("& \\emph{t} stat. & \\emph{t p}-value &KS \\emph{p}-value &&",
    "\\emph{t} stat. & \\emph{t p}-value &KS \\emph{p}-value", "\\", "\\", "\n",
    sep="")
  z <- sapply(1:dim(covariates)[2], function(x){
    cat(dimnames(covariates)[[2]][x], "&", round(bal.out$BeforeMatching[[x]]$tt$statistic,2), "&",
    round(bal.out$BeforeMatching[[x]]$tt$p.value,2), "&",
    ifelse(is.null(bal.out$BeforeMatching[[x]]$ks$ks.boot.pvalue) == 0,
      round(bal.out$BeforeMatching[[x]]$ks$ks.boot.pvalue,2), "---"), "&&",
    round(bal.out$AfterMatching[[x]]$tt$statistic,2), "&",
    round(bal.out$AfterMatching[[x]]$tt$p.value,2), "&",
    ifelse(is.null(bal.out$AfterMatching[[x]]$ks$ks.boot.pvalue) == 0,
      round(bal.out$AfterMatching[[x]]$ks$ks.boot.pvalue,2), "---"), "\\", "\\", "\n",
      sep="")
  })
  cat("\\end{tabular} \\end{center} \n")
}

## Note: KS test statistics appear as 'NA' for dichotomous variables. Also,
##   'BM' is before matching and 'AM' is after matching.
balanceMat <- function(covariates, bal.out){
  z <- t(sapply(1:dim(covariates)[2], function(x){
    c(dimnames(covariates)[[2]][x], round(bal.out$BeforeMatching[[x]]$tt$statistic,2),
    round(bal.out$BeforeMatching[[x]]$tt$p.value,2),
    ifelse(is.null(bal.out$BeforeMatching[[x]]$ks$ks.boot.pvalue) == 0,
      round(bal.out$BeforeMatching[[x]]$ks$ks.boot.pvalue,2), NA),
    round(bal.out$AfterMatching[[x]]$tt$statistic,2),
    round(bal.out$AfterMatching[[x]]$tt$p.value,2),
    ifelse(is.null(bal.out$AfterMatching[[x]]$ks$ks.boot.pvalue) == 0,
      round(bal.out$AfterMatching[[x]]$ks$ks.boot.pvalue,2), NA))
  }))
  z <- as.data.frame(z)
  z[,2:7] <- apply(z[,2:7], 2, function(x){as.numeric(x)})
  mat <- z[,2:7]
  dimnames(mat)[[2]] <- c("BM t stat", "BM t p-value", "BM KS stat",
    "AM t stat", "AM t p-value", "AM KS stat")
  dimnames(mat)[[1]] <- z[,1]
  mat
}