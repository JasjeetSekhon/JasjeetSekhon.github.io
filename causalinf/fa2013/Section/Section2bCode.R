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

data$Polar=data$Polar-mean(data$Polar)




# Now look at just the treatment effect for polar bears

model=lm(CampSiteRaids~Treat,data[data$Type=="Polar",])
model$newse<-vcovHC(model,type="HC1")
coeftest(model,model$newse)

# So this treatment effect equals the treatment effect from the previous model plus the interaction term.




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

omega=function(residuals, diaghat, df){return(12/(12-6)*residuals^2)}

model=lm(CampSiteRaids~Treat+Type+Polar*Treat+Brown*Treat,data.original)
model$newse<-vcovHC(model,omega=omega)
coeftest(model,model$newse,df=6)









