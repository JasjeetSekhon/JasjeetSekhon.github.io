data=data.frame(c(1:12),c("Black","Brown","Black","Polar","Polar","Brown","Black","Black","Polar","Brown","Brown","Polar"),c(0,1,0,0,1,1,0,1,1,0,0,1),c(22,10,27,26,13,14,33,16,10,25,26,12))

colnames(data)=c("Bear","Type","Treat","CampSiteRaids")

t.test(data[data$Treat==1,]$CampSiteRaids,data[data$Treat==0,]$CampSiteRaids)

model=lm(CampSiteRaids~Treat,data)

summary(model)

require("sandwich")
require("lmtest")
model$newse<-vcovHC(model)
coeftest(model,model$newse)


model$newse<-vcovHC(model,type="HC1")
coeftest(model,model$newse)


data$Polar=as.numeric(data$Type=="Polar")

model=lm(CampSiteRaids~Treat+Polar,data)
model$newse<-vcovHC(model,type="HC1")
coeftest(model,model$newse)

# Now add in an interactoin for Polar

model=lm(CampSiteRaids~Treat+Polar+Treat*Polar,data)
model$newse<-vcovHC(model,type="HC1")
coeftest(model,model$newse)

# Now look at just the treatment effect for polar bears

model=lm(CampSiteRaids~Treat,data[data$Type=="Polar",])
model$newse<-vcovHC(model,type="HC1")
coeftest(model,model$newse)

# So this treatment effect equals the treatment effect from the previous model plus the interaction term. Now look at the treatment effect for non-polar bears

model=lm(CampSiteRaids~Treat,data[data$Type!="Polar",])
model$newse<-vcovHC(model,type="HC1")
coeftest(model,model$newse)

# So this is just beta_1 from the previous model.


# To get the adjusted treatment effect for the entire sample, the correct model is 

data$Polar=data$Polar-mean(data$Polar) # Converts to mean deviations

model=lm(CampSiteRaids~Treat,data[data$Type!="Polar",])
model$newse<-vcovHC(model,type="HC1")
coeftest(model,model$newse)

# The estimator is beta_1, and it is consistent and should have less variance then tau_hat

# Use the same robust standard errors as Stata

data[order(data$Type),]

mean(data[data$Type=="Black",]$CampSiteRaids)
mean(data[data$Type=="Brown",]$CampSiteRaids)
mean(data[data$Type=="Polar",]$CampSiteRaids)

data.original=data

for(i in unique(data$Type)){

data[data$Type==i,]$CampSiteRaids=data[data$Type==i,]$CampSiteRaids-mean(data[data$Type==i,]$CampSiteRaids)

}

data[order(data$Type),]

for(i in unique(data$Type)){

data[data$Type==i,]$Treat=data[data$Type==i,]$Treat-mean(data[data$Type==i,]$Treat)

}

data[order(data$Type),]

model=lm(CampSiteRaids~Treat,data)
model$newse<-vcovHC(model,type="HC1")
coeftest(model,model$newse)

model=lm(CampSiteRaids~Treat+Type,data.original)
model$newse<-vcovHC(model,type="HC1")
coeftest(model,model$newse)

summary(lm(CampSiteRaids~Treat,data),df=8)

summary(lm(CampSiteRaids~Treat+Type,data.original))


omega=function(residuals, diaghat, df){return(12/(12-4)*residuals^2)}

model=lm(CampSiteRaids~Treat,data)
model$newse<-vcovHC(model,omega=omega)
coeftest(model,model$newse,df=8)

model=lm(CampSiteRaids~Treat+Type,data.original)
model$newse<-vcovHC(model,omega=omega)
coeftest(model,model$newse,df=8)




# To get the estimator that is guarenteed to be unbiased, control for the interactions.

data.original$Brown=as.numeric(data.original$Type=="Brown")

data.original$Brown=data.original$Brown-mean(data.original$Brown) # Converts to mean deviations

omega=function(residuals, diaghat, df){return(12/(12-6)*residuals^2)}

model=lm(CampSiteRaids~Treat+Type+Polar*Treat+Brown*Treat,data.original)
model$newse<-vcovHC(model,omega=omega)
coeftest(model,model$newse,df=6)





# Here is the solution code to the Social Pressure paper



social=read.csv("VotingExperiment.csv")

# First create dummy variables for each of the treatments. Since they are factors, you have to use the numbers that represent them (Civic Duty=1, Hawthorne=3, Neighbors=4, and Self=5)

social$CivicDuty=as.numeric(as.numeric(social$treatment)==1)
social$Hawthorne=as.numeric(as.numeric(social$treatment)==3)
social$Neighbors=as.numeric(as.numeric(social$treatment)==4)
social$Self=as.numeric(as.numeric(social$treatment)==5)

# Now change voted into a binary variable so the function lm() can handle it

social$voted=as.numeric(social$voted=="Yes")

# Next, install Drew's function, which will allow you to cluster the standard errors by household

robust.se <- function(model, cluster){
 require(sandwich)
 require(lmtest)
 M <- length(unique(cluster))
 N <- length(cluster)
 K <- model$rank
 dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
 uj <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
 rcse.cov <- dfc * sandwich(model, meat = crossprod(uj)/N)
 rcse.se <- coeftest(model, rcse.cov)
 return(list(rcse.cov, rcse.se))
}

# Now compute the robust standard errors, clustering by household. The first object that will be produced is the Cov matrix. The second is the results from the regression.

model=lm(voted~CivicDuty+Hawthorne+Neighbors+Self,social)
robust.se(model,social$hh_id)




# Now we will use a model using fixed effects at the neighborhood level. We don't want to estimate parameters for all 10,000 of the neighborhoods, so we will take out the mean deviations for each neighborhood. There are many ways to do this, but some take a very long time. However, there is a way to do it very quickly.

means=ave(social$voted,social$cluster)

social$voted=social$vote-means

# Now we need to change the treatment assingments for each person from (0,1) to the mean deviations 

means=ave(social$CivicDuty,social$cluster)
social$CivicDuty=social$CivicDuty-means

means=ave(social$Hawthorne,social$cluster)
social$Hawthorne=social$Hawthorne-means

means=ave(social$Neighbors,social$cluster)
social$Neighbors=social$Neighbors-means

means=ave(social$Self,social$cluster)
social$Self=social$Self-means

# Now we need to change the degrees of freedom. We normally have 344084-4-1 degrees of freedom, since we have four treatments. Now we our estimating 9,999 new parameters, so we will have 344084-4-1-9999=334080 degrees of freedom

 robust.se <- function(model, cluster){
 require(sandwich)
 require(lmtest)
 M <- length(unique(cluster))
 N <- length(cluster)
 K <- model$rank
 dfc <- (M/(M - 1)) * ((N - 1)/(N - K - 9999)) # (K=4+1)
 uj <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
 rcse.cov <- dfc * sandwich(model, meat = crossprod(uj)/N)
 rcse.se <- coeftest(model, rcse.cov, df=334080)
 return(list(rcse.cov, rcse.se))
}

# Now we can get our results

model=lm(voted~CivicDuty+Hawthorne+Neighbors+Self,social)
robust.se(model,social$hh_id)


# Next, we can just add in the controls that they list in their table. We will have to decrease our degrees of freedom by 5.

 robust.se <- function(model, cluster){
 require(sandwich)
 require(lmtest)
 M <- length(unique(cluster))
 N <- length(cluster)
 K <- model$rank
 dfc <- (M/(M - 1)) * ((N - 1)/(N - K - 9999)) # (K=4+1)
 uj <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
 rcse.cov <- dfc * sandwich(model, meat = crossprod(uj)/N)
 rcse.se <- coeftest(model, rcse.cov, df=334075)
 return(list(rcse.cov, rcse.se))
}

# Now change are 5 controls to from "yes"/"no" to 0,1

social$g2000=as.numeric(social$g2000=="yes")
social$g2002=as.numeric(social$g2002=="yes")
social$p2000=as.numeric(social$p2000=="yes")
social$p2002=as.numeric(social$p2002=="yes")
social$p2004=as.numeric(social$p2004=="Yes")

# Now calculate the mean deviations for these controls

means=ave(social$g2000,social$cluster)
social$g2000=social$g2000-means

means=ave(social$g2002,social$cluster)
social$g2002=social$g2002-means

means=ave(social$p2000,social$cluster)
social$p2000=social$p2000-means

means=ave(social$p2002,social$cluster)
social$p2002=social$p2002-means

means=ave(social$p2004,social$cluster)
social$p2004=social$p2004-means

model=lm(voted~CivicDuty+Hawthorne+Neighbors+Self+g2000+g2002+p2000+p2002+p2004,social)
robust.se(model,social$hh_id)


# Now we should do a model with interactions. This will require us to add 5*4 interaction terms, so e will have to adjust the degrees of freedom accordingly.

 robust.se <- function(model, cluster){
 require(sandwich)
 require(lmtest)
 M <- length(unique(cluster))
 N <- length(cluster)
 K <- model$rank
 dfc <- (M/(M - 1)) * ((N - 1)/(N - K - 9999-20)) # (K=4+1)
 uj <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
 rcse.cov <- dfc * sandwich(model, meat = crossprod(uj)/N)
 rcse.se <- coeftest(model, rcse.cov, df=334055)
 return(list(rcse.cov, rcse.se))
}

model=lm(voted~CivicDuty+Hawthorne+Neighbors+Self+g2000+g2002+p2000+p2002+p2004+CivicDuty*g2000+CivicDuty*g2002+CivicDuty*p2000+CivicDuty*p2002+CivicDuty*p2004+Hawthorne*g2000+Hawthorne*g2002+Hawthorne*p2000+Hawthorne*p2002+Hawthorne*p2004+Neighbors*g2000+Neighbors*g2002+Neighbors*p2000+Neighbors*p2002+Neighbors*p2004+Self*g2000+Self*g2002+Self*p2000+Self*p2002+Self*p2004,social)
robust.se(model,social$hh_id)


# Lastly, we can add in all the other control variables. The first step is to turn everything into numbers.

social$sex=as.numeric(social$sex=="female")
social$g2004=as.numeric(social$g2004=="yes")

# Next we should take the mean deviations of our variables (using the mean in each neighborhood)

means=ave(social$sex,social$cluster)
social$sex=social$sex-means

means=ave(social$yob,social$cluster)
social$yob=social$yob-means

means=ave(social$g2004,social$cluster)
social$g2004=social$g2004-means

means=ave(social$hh_size,social$cluster)
social$hh_size=social$hh_size-means

means=ave(social$numberofnames,social$cluster)
social$numberofnames=social$numberofnames-means

means=ave(social$p2004_mean,social$cluster)
social$p2004_mean=social$p2004_mean-means

means=ave(social$g2004_mean,social$cluster)
social$g2004_mean=social$g2004_mean-means

 robust.se <- function(model, cluster){
 require(sandwich)
 require(lmtest)
 M <- length(unique(cluster))
 N <- length(cluster)
 K <- model$rank
 dfc <- (M/(M - 1)) * ((N - 1)/(N - K - 9999-55)) # (K=4+1)
 uj <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
 rcse.cov <- dfc * sandwich(model, meat = crossprod(uj)/N)
 rcse.se <- coeftest(model, rcse.cov, df=334013)
 return(list(rcse.cov, rcse.se))
}

model=lm(voted~CivicDuty+Hawthorne+Neighbors+Self+g2000+g2002+p2000+p2002+p2004+sex+yob+g2004+hh_size+numberofnames+p2004_mean+g2004_mean+CivicDuty*g2000+CivicDuty*g2002+CivicDuty*p2000+CivicDuty*p2002+CivicDuty*p2004+CivicDuty*sex+CivicDuty*yob+CivicDuty*g2004+CivicDuty*hh_size+CivicDuty*numberofnames+CivicDuty*p2004_mean+CivicDuty*g2004_mean+Hawthorne*g2000+Hawthorne*g2002+Hawthorne*p2000+Hawthorne*p2002+Hawthorne*p2004+Hawthorne*sex+Hawthorne*yob+Hawthorne*g2004+Hawthorne*hh_size+Hawthorne*numberofnames+Hawthorne*p2004_mean+Hawthorne*g2004_mean+Neighbors*g2000+Neighbors*g2002+Neighbors*p2000+Neighbors*p2002+Neighbors*p2004+Neighbors*sex+Neighbors*yob+Neighbors*g2004+Neighbors*hh_size+Neighbors*numberofnames+Neighbors*p2004_mean+Neighbors*g2004_mean+Self*g2000+Self*g2002+Self*p2000+Self*p2002+Self*p2004+Self*sex+Self*yob+Self*g2004+Self*hh_size+Self*numberofnames+Self*p2004_mean+Self*g2004_mean,social)
robust.se(model,social$hh_id)


