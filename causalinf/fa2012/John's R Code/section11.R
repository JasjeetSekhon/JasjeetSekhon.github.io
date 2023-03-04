# A simple monte carlo showing IV in the potential outcomes framework

library(sem)

bias <- matrix(nrow=5000)
for (i in 1:5000){

  #Generate potential outcomes
  #potential outcome when instrument and treatment are 0
r_0_0 <- runif(100)*10
  #Potential outcome when treatment is 0, but instrument is 1
r_0_1 <- r_0_0 #because of exclusion restriction
#Potential outcome when instrument and treatment are 1
r_1_1 <- r_0_0 + rnorm(100,mean=1)
#Potential outcome when treatment is 1 and instrument is 0
r_1_0 <- r_1_1


#plot(density(r_1_1),col="red")
#lines(density(r_0_0),col="blue")


#Create 2 treatment vectors, one if the instrument is 1, and the other if it is 0
gamma <- rnorm(100)
t_0 <- ifelse((gamma)>.5,1,0)
#strong instrument
t_1 <- ifelse((gamma + rnorm(100,1.5,.5)) > .5, 1,0)
#weak instrument
#t_1 <- ifelse((gamma + rnorm(100,.5,.5)) > .5, 1,0)


#how many compliers?
sum((t_1-t_0)==1)
#How many defiers?
sum((t_1-t_0)==-1)
#how many always takers?
sum((t_1==1&t_0==1))
#how many never takers?
sum((t_1==0&t_0==0))


#remove defiers
t_0[(t_1-t_0)==-1] <- 0

#Let's make compliers different from the over all population
r_1_1[(t_1-t_0)==1] <- r_1_1[(t_1-t_0)==1] + rnorm(length(r_1_1[(t_1-t_0)==1]),2,1)



#Generate Instrument
z <- matrix(0,nrow=100)
z[sample(1:100,50)] <- 1

#Generate realized treatment and outcome vectors
t <- ifelse(z==1,t_1,t_0)
r <- ifelse(t==1,r_1_1, r_0_0)

#What is the reduced form estimates?
mean(r[z==1])-mean(r[z==0])

#What is the effect of the instrument on treatment?
mean(t[z==1]) - mean(t[z==0])

#What is the Wald estimate?
(mean(r[z==1])-mean(r[z==0]))/(mean(t[z==1]) - mean(t[z==0]))

#what is the truth?
#average treatment effect of treatment on compliers
 (mean(r_1_1[(t_1-t_0)==1]) - mean(r_0_0[(t_1-t_0)==1]))
#average treatment effect
#True average treatment effect
mean(r_1_1 - r_0_0 )

#What is the bias?
bias[i] <- (mean(r_1_1[(t_1-t_0)==1]) - mean(r_0_0[(t_1-t_0)==1])) - ((mean(r[z==1])-mean(r[z==0]))/(mean(t[z==1]) - mean(t[z==0])))


}

plot(density(bias),col="red")
abline(v=0,lwd=2)

#Use two stage least squares to generate estimates and standard errors
data <- data.frame(r=r,t=t,z=z)
summary(tsls(r~t,~z,data=data))
