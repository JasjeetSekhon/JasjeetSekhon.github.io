### PS 239 Section, January 2009              ###
### Time series bootstrap (from PS 236 exam)  ###

### Load and arrange the data
dta <- read.table(file="http://sekhon.berkeley.edu/causalinf/midterm/FearonLaitin2003.asc", header=TRUE)
attach(dta)
X <- cbind(warl,lpopl1,lmtnest,ncontig,Oil,nwstate,instab,polity2l,ethfrac,relfrac)
X <- as.data.frame(X)
detach(dta)
gdpen <- dta$gdpen
gdpenl <- dta$gdpenl
war <- dta$onset
war <- ifelse(war > 1, 1, war)
warl <- dta$warl
summary(X)

## Simple way of dealing with NAs --- exclude any row with an NA in any covariate
##   Note: For the first term below, 'apply' is summing each row (since we
##   specify '1') in X. If the row contains an 'NA', the sum is 'NA'. Hence, we
##   are keeping only those rows that do not contain any NAs in X. Then we check
##   for NAs in GDP and its lagged value for each country-year.
Xrm <- X[is.na(apply(X,1,sum)) == 0 & is.na(gdpen) == 0 & is.na(gdpenl) == 0,]
gdpenrm <- gdpen[is.na(apply(X,1,sum)) == 0 & is.na(gdpen) == 0 & is.na(gdpenl) == 0]
gdpenlrm <- gdpenl[is.na(apply(X,1,sum)) == 0 & is.na(gdpen) == 0 & is.na(gdpenl) == 0]
country <- dta$cname[is.na(apply(X,1,sum)) == 0 & is.na(gdpen) == 0 & is.na(gdpenl) == 0]

## Now, estimate the betas on this data set
lm2 <- glm(gdpenrm ~ gdpenlrm + as.matrix(Xrm))

## Notice that the betas do not change when we discarded the NAs (since glm
##   discards NAs anyway)
lm1  <- glm(gdpen ~ warl + gdpenl + lpopl1 + lmtnest + ncontig + Oil + nwstate
  + instab + polity2l + ethfrac + relfrac, data=dta)
summary(lm1)
summary(lm2)



## How the bootstrap works:
##   Take a given country
##   Find the GDP 'gdp', covariates 'gdpLag' and 'covars', and residuals 'resFS'
##     for that country
##   Shuffle the residuals of that country with replacement
##   Leave the first lagged GDP measure as-is
##   Calculate the first bootstrapped GDP by multiplying the first lagged GDP
##     and the other X covariates by the full-sample coefficients and add the
##     shuffled residual
##   Calculate subsequent bootstrapped GDPs by following the same procedure, but
##     use the previous GDP calculation as the new lagged GDP
##   Save the bootstrapped GDP and join the new bootstrapped lagged GDP to the
##     other X covariates
##   Repeat for all countries
##   Join together into a bootstrapped data set
##   Calculate the bootstrapped betas
##   Repeat many times

tsOLSBootC <- function(gdp, gdpLag, covars, country, nBoots=1000, alpha=0.05){
  # Make country a factor and remove unused levels
  country <- as.factor(country)
  country <- country[,drop=TRUE]
  # 'covars' must be a matrix
  covars <- as.matrix(covars)
  # Estimate full sample coefficients and residuals
  fs <- glm(gdp ~ gdpLag + covars)
  betaFS <- fs$coefficients
  resFS <- fs$residuals
  # Estimating bootstrapped coefficients
  Z <- sapply(1:nBoots, function(run){
    # Create a starting template for the bootrapped sample with enough columns
    #   to cover GDP, lagged GDP, and the other covariates
    bootSample <- rep(NA, times=(1 + 1 + dim(covars)[2]))
    # Create the bootstrapped sample country-by-country
    for (x in levels(as.factor(country))){
      # Isolate covariates for country x
      ctyCov <- covars[country == x, ]
      # Find the first lagged GDP for that country
      firstLag <- gdpLag[country == x][1]
      # Get residuals for that country
      ctyRes <- resFS[country == x]
      # Shuffle the residuals with replacement
      randRes <- sample(ctyRes, replace=TRUE)
      # Multiply a one (for the intercept), the first lag, and first set of
      #   covariates times the full-sample betas and add a shuffled residual
      bootCtyGDP <- crossprod(c(1, firstLag, ctyCov[1,]), betaFS) + randRes[1]
      # Now, for subsequent time periods
      for (i in 2:dim(ctyCov)[1]){
        # Like above, but use previous period's (i.e., period i-1) lagged GDP
        #   with period i covariates
        nextBootGDP <- crossprod(c(1, bootCtyGDP[(i-1)], ctyCov[i,]), betaFS) + randRes[i]
        # String bootstrapped Country GDPs together
        bootCtyGDP <- c(bootCtyGDP, nextBootGDP)
      }
      # Create a country data set of bootstrapped GDP, bootstrapped lagged GDP,
      #   and other covariates
      #   Note: Using '-length()' excludes the last entry, since we don't use
      #   the last year's GDP for a lag (no obersvations after the last)
      ctyBoot <- cbind(bootCtyGDP, c(firstLag, bootCtyGDP[-length(bootCtyGDP)]),
        ctyCov)
      # String country data sets together
      bootSample <- rbind(bootSample, ctyBoot)
    }
    # Remove the started NA row
    bootSample <- bootSample[-c(1),]
    # Calculate bootstrapped betas
    bootLM <- glm(bootSample[,1] ~ bootSample[, -c(1)])
    betaBoot <- bootLM$coefficients
  })
  # Create confidence interval output
  CI <- cbind(apply(Z,1,function(x){quantile(x, (alpha/2))}), rowMeans(Z),
    betaFS, apply(Z,1,function(x){quantile(x, (1-alpha/2))}))
  CI <- round(CI, 4)
  dimnames(CI)[[2]] <- c("Lower bound", "Estimate", "Full sample", "Upper bound")
  dimnames(CI)[[1]] <- c("Intercept", "Lagged GDP", dimnames(covars)[[2]])
  cat("Bootstrapped Estimates and ", (1-alpha)*100,"% Confidence Intervals, OLS coefficients \n", sep="")
  CI
}

bootC <- tsOLSBootC(gdpenrm, gdpenlrm, Xrm, country, nBoots=1000)


## Now, let's change 'tsOLSBootC' in only one way: add a residual chosen at
##   random from all residuals from all countries, rather than choosing it only
##   from the country's own residuals (change highlighted below)

tsOLSBoot <- function(gdp, gdpLag, covars, country, nBoots=1000, alpha=0.05){
  country <- as.factor(country)
  country <- country[,drop=TRUE]
  covars <- as.matrix(covars)
  fs <- glm(gdp ~ gdpLag + covars)
  betaFS <- fs$coefficients
  resFS <- fs$residuals
  Z <- sapply(1:nBoots, function(run){
    bootSample <- rep(NA, times=(1 + 1 + dim(covars)[2]))
    for (x in levels(as.factor(country))){
      ctyCov <- covars[country == x, ]
      firstLag <- gdpLag[country == x][1]
      randRes <- sample(resFS, size=dim(ctyCov)[1], replace=TRUE)
      bootCtyGDP <- crossprod(c(1, firstLag, ctyCov[1,]), betaFS) + randRes[1]
      for (i in 2:dim(ctyCov)[1]){
        nextBootGDP <- crossprod(c(1, bootCtyGDP[(i-1)], ctyCov[i,]), betaFS) + randRes[i]
        bootCtyGDP <- c(bootCtyGDP, nextBootGDP)
      }
      ctyBoot <- cbind(bootCtyGDP, c(firstLag, bootCtyGDP[-length(bootCtyGDP)]),
        ctyCov)
      bootSample <- rbind(bootSample, ctyBoot)
    }
    bootSample <- bootSample[-c(1),]
    bootLM <- glm(bootSample[,1] ~ bootSample[, -c(1)])
    betaBoot <- bootLM$coefficients
  })
  CI <- cbind(apply(Z,1,function(x){quantile(x, (alpha/2))}), rowMeans(Z),
    betaFS, apply(Z,1,function(x){quantile(x, (1-alpha/2))}))
  CI <- round(CI, 4)
  dimnames(CI)[[2]] <- c("Lower bound", "Estimate", "Full sample", "Upper bound")
  dimnames(CI)[[1]] <- c("Intercept", "Lagged GDP", dimnames(covars)[[2]])
  cat("Bootstrapped Estimates and ", (1-alpha)*100,"% Confidence Intervals, OLS coefficients \n", sep="")
  CI
}

bootAR <- tsOLSBoot(gdpenrm, gdpenlrm, Xrm, country, nBoots=1000)