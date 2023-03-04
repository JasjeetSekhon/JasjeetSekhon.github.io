#
# Modifying the Fearon and Laitin Model 1 (Table 1) to estimate the effect on per capita income
#


dta <- read.table(file="http://sekhon.berkeley.edu/qe/FearonLaitin2003.asc", header=TRUE)
#dta <- read.table(file="FearonLaitin2003.asc", header=TRUE)

#if you prefer to use the stata data set.
#library(foreign)
#dta <- read.dta(file="FearonLaitin2003.dta")

# a few variables which may be of interest:
# "ccode" a country code (numeric)
# "country" country name, long form
# "cname"  country name, short form
# "cmark" an indicator variable for when a new country starts (1 for the first observation for a new country)
# "year" the observation year

#recode all values for the dependant variables for logistic regressions, so that obs that are > 1 are set equal to 1
indx  <- dta$onset > 1 
dta$onset[indx]  <- 1

indx  <- dta$ethonset > 1 
dta$ethonset[indx]  <- 1

indx  <- dta$emponset > 1 
dta$emponset[indx]  <- 1

indx  <- dta$cowonset > 1 
dta$cowonset[indx]  <- 1


#a model for per capital income based on model 1 in table 1 of Fearon
#and Laitin 2003, except that we are changing what the dependent
#variable is.  For replication of Table 1 see the 'FearonLaitin_replication1.R' file.

#note that in the Fearon and Laitin dataset income per capial at time
#t-1 has more significant digits then at t, to make the number of
#significanticant digits equation we do the following:
dta$gdpenl <-  round(dta$gdpenl, 3)

# to see that gdpenl is the log of gdpen
#cbind(dta$ccode, dta$year, dta$gdpen, dta$gdpenl)

#Per capital income model
lm1  <- glm(gdpen ~ warl + gdpenl +lpopl1 +lmtnest +ncontig +Oil +nwstate +instab +polity2l +ethfrac +relfrac
, 
             data=dta) 
summary(lm1)

