#
# R version of replicating Table 1 of Fearon and Latin 2003 in R
#

# a few variables which may be of interest:
# "ccode" a country code (numeric)
# "country" country name, long form
# "cname"  country name, short form
# "cmark" an indicator variable for when a new country starts (1 for the first observation for a new country)
# "year" the observation year

# all of the other variable names you need to know are defined in the
# regressions below (just compare with variable names here with the
# names in Table 1 in the article to figure out what they are.
# Compare cofficients etc).

dta <- read.table(file="http://sekhon.berkeley.edu/qe/FearonLaitin2003.asc", header=TRUE)
#dta <- read.table(file="FearonLaitin2003.asc", header=TRUE)

#if you prefer to use the stata data set.
#library(foreign)
#dta <- read.dta(file="FearonLaitin2003.dta")

#recode all values for the dependant variables, so that obs that are > 1 are set equal to 1
indx  <- dta$onset > 1 
dta$onset[indx]  <- 1

indx  <- dta$ethonset > 1 
dta$ethonset[indx]  <- 1

indx  <- dta$emponset > 1 
dta$emponset[indx]  <- 1

indx  <- dta$cowonset > 1 
dta$cowonset[indx]  <- 1

#Model 1
glm1  <- glm(onset ~ warl + gdpenl +lpopl1 +lmtnest +ncontig +Oil +nwstate +instab +polity2l +ethfrac +relfrac, 
             data=dta, family=binomial) 
summary(glm1)


#Model 2, second > .049999
subset.indx <- dta$second > .049999
glm2  <- glm(ethonset ~ warl + gdpenl + lpopl1 + lmtnest + ncontig + Oil + nwstate + instab + polity2l + ethfrac + relfrac,
             data=dta, family=binomial, subset=subset.indx) 
summary(glm2)

#Model 3
glm3  <- glm(onset ~ warl + gdpenl + lpopl1 + lmtnest + ncontig + Oil + nwstate + instab + anocl + deml + ethfrac + relfrac, 
             data=dta, family=binomial) 
summary(glm3)


#Model 4
glm4  <- glm(emponset ~ empwarl+ empgdpenl +emplpopl+ emplmtnest+ empncontig+ Oil+ nwstate+ instab+ empethfrac,
             data=dta, family=binomial) 
summary(glm4)


#Model 5
glm5  <- glm(cowonset ~ cowwarl + gdpenl + lpopl1 + lmtnest + ncontig + Oil + nwstate + instab + anocl + deml + ethfrac + relfrac,
             data=dta, family=binomial) 
summary(glm5)


