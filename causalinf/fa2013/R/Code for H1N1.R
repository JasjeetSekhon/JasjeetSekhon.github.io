library(Matching)
data=read.csv("Final_H1N1.csv")

# data[which(data==".")]=NA

length(which(data$ECMODAYS>0))
length(which(data$RRTDAYS>0))
length(which(data$VASODAYS>0))
length(which(data$STEROIDS=="1"))

which(is.na(data$ECMODAYS)==TRUE)

data$BMI[data$BMI=="."]=NA

controls=with(data[-c(63,64),],cbind(VENTDAYS, RRTDAYS, VASODAYS, APACHEII, APACHEIII, AGE, GENDER, PREG, OTHERCOM, as.numeric(BMI)>35, FLUSYNDROME, STEROIDS))

gen=GenMatch(Tr=data$ECMODAYS[-c(63,64)]>0, X=controls, estimand="ATT")

mat=Match(Y=rep(0,62),Tr=data$ECMODAYS[-c(63,64)]>0, X=controls,Weight.matrix=gen, estimand="ATT")

MatchBalance(ECMODAYS==TRUE ~ PREG + AGE, data=data[,-c(63,64)], match.out=mat, na.omit=TRUE)