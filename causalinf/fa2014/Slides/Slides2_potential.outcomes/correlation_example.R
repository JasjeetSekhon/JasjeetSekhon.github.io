
rm(list=ls(all=TRUE))
set.seed(12345)

#############################################
# Option 1: R=10000
#############################################


### Parameters:
m=4
R=1000

n.vec = c(c(5:20),seq(21,100,by=5))
cov.real1 <- cov.approx1 <- rep(999,length(n.vec))

for (i in c(1:length(n.vec))){
  
  N = n.vec[i]
  cat("sample size:  ",N,"\n")
  
  ## analytical:
  cov.real1[i] <- (m/N)*((m-1)/(N-1)-(m/N))
  
  ### Simulation:
  
  z1<-z2<-rep(999,R)
  for (j in c(1:R)){
    id.treat = sample(c(1:N),m,replace=FALSE)
    treat0 = rep(0,N)
    treat0[id.treat]=1
    z1[j] = treat0[1]
    z2[j] = treat0[2]
  }
  
  cov.approx1[i] <- cov(z1,z2)
}

par(cex=0.7)
plot(n.vec,cov.real1,las=1,col="blue2",pch=16)
points(n.vec,cov.approx1,las=1,col="red2",pch=16)

#############################################
# Option 2: R=500000s
#############################################


### Parameters:
m=4
R=500000

n.vec = c(c(5:20),seq(21,100,by=5))
cov.real2 <- cov.approx2 <- rep(999,length(n.vec))

for (i in c(1:length(n.vec))){
  
  N = n.vec[i]
  cat("sample size:  ",N,"\n")
  
  ## analytical:
  cov.real2[i] <- (m/N)*((m-1)/(N-1)-(m/N))
  
  ### Simulation:
  
  z1<-z2<-rep(999,R)
  for (j in c(1:R)){
    id.treat = sample(c(1:N),m,replace=FALSE)
    treat0 = rep(0,N)
    treat0[id.treat]=1
    z1[j] = treat0[1]
    z2[j] = treat0[2]
  }
  
  cov.approx2[i] <- cov(z1,z2)
}

pdf("figure_cov.pdf",width=4.5,height=4)
par(cex=0.7)
plot(n.vec,cov.approx1,las=1,col="green4",pch=16,
     ylab="Covariance",xlab="Sample size (N)")
points(n.vec,cov.approx2,las=1,col="red3",pch=16)
points(n.vec,cov.real2,las=1,col="blue3",pch=16)
dev.off()

#ggplot(data.frame(N=n.vec,approx1=cov.approx1,approx2=cov.approx2,cov.real=cov.real1))+
#  geom_point(aes(x=N,y=approx1,col="blue"))



### Correlation between the mean in the control and the mean in th treatment group:
set.seed(1)

### Parameters:
N=4
m=2
R=1000000

# outcome:
y = rnorm(N,mean=0,sd=1)

## analytical:
# ???
set.seed(12345)
### Simulation:

mean_t<-mean_c<-rep(999,R)
for (i in c(1:R)){
  id.treat = sample(c(1:N),m,replace=FALSE)
  treat0 = rep(0,N)
  treat0[id.treat]=1
  mean_c[i] = mean(y[treat0==0])
  mean_t[i] = mean(y[treat0==1])
}

cov(mean_c,mean_t)
var(mean_c)

diff.mean = mean_t - mean_c
var(diff.mean)


#############################################
# Difference in means variance:
#############################################

set.seed(12345)

m=2
N=4
p = m/N
cov.z = (m/N)*((m-1)/(N-1)-(m/N))
var.z = p*(1-p)

y = rnorm(N,mean=5,sd=1)

### analytical variance:
sum_a2 = sum(y^2)

sum_aa = rep(999,N)
for (i in c(1:N)){
  sum_aa[i] = sum(y[i]*y[-i])
}

sum_aa = sum(sum_aa)

cov.diff = -(1/m)*(1/(N-m))*(sum_a2*var.z+sum_aa*cov.z)
var.z1 = var.z*sum_a2*(1/m^2)
var.z0 = var.z*sum_a2*(1/(N-m)^2)

var.diff = var.z1+var.z0-2*cov.diff


##### Simulation:
R = 100000
z = c(rep(1,m),rep(0,N-m))

mean.diff = rep(999,R)
for (i in c(1:R)){
  
  z0 = sample(z,N)
  mean.diff[i] = mean(y[z0==1]) - mean(y[z0==0])
}

var(mean.diff)



#### Example - Nick:

sum_a2 = sum(y^2)

sum_aa = rep(999,N)
for (i in c(1:N)){
  sum_aa[i] = sum(y[i]*y[-i])
}

sum_aa = sum(sum_aa)

cov.diff = -(1/m)*(1/(N-m))*(sum_a2*var.z+sum_aa*cov.z)
var.z1 = var.z*sum_a2*(1/m^2)
var.z0 = var.z*sum_a2*(1/(N-m)^2)

var.diff = var.z1+var.z0-2*cov.diff


##### Simulation:
m=3
N=12

set.seed(12345)
y = rnorm(N,mean=0,sd=1)
y

R = 20000
z = c(rep(1,m),rep(0,N-m))

mean.diff<-mean.t <-mean.c <- rep(999,R)
for (i in c(1:R)){
  z0 = sample(z,N)
  
  mean.c[i] = mean(y[z0==0])
  mean.t[i] = mean(y[z0==1])
  mean.diff[i] = mean(y[z0==1]) - mean(y[z0==0])
}
cov(mean.c,mean.t)
var(mean.diff)

### Analytical
sigma2 = mean(y^2) - mean(y)^2

sigma2*N^2/((N-1)*m*(N-m))
(N/(N-1))*(sigma2/(N-m)+sigma2/m)

var(mean.c)+var(mean.t)-2*cov(mean.c,mean.t)

(sigma2/(N-m))*((m)/(N-1))+sigma2/m*((N-m)/(N-1))-2*cov(mean.c,mean.t)




