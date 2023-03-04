###############################################
# Section 2: MLE
#
# John Henderson         
#
# PS 236B Spring 2013
###############################################

rm(list=ls(all=TRUE))
#Change path to correct directory
      
# Multivariate Regresion Using Newtons method for MLE: alpha, beta, sigma
                                          
# randomly sampling from the following process
set.seed(1005)
x=rnorm(mean=1.2,sd=1.2,n=1000)
alpha=-2.3
beta=3.44
sigma=(2.2)
y=alpha+beta*x+rnorm(mean=0,sd=(sigma),n=1000) 

plot(x,y)   
        
# all the parameters combined first in pars  
loglike=function(pars,y,x){
	alpha=pars[1]
	beta=pars[2]
	sigma=pars[3]
	n=length(y)                                                    
	X=cbind(1,x) 
	b=c(alpha,beta)
	llik=-(n/2)*log(sigma^2)-(n/2)*log(2*pi)-(1/(2*(sigma^2)))*t(y-X%*%b)%*%(y-X%*%b)	
	return(-llik)	          
}
   

# 1. Newton 
#   - optimize using Newton method of analytical derivatives     
#   - analytical derivatives...                                                                         
jacobian=function(y,x,alpha,beta,sigma){
	b=c(alpha,beta)
	X=cbind(1,x)
	n=length(y)
	drvB=c((1/(sigma^2))*((t(X)%*%y-(t(X)%*%X)%*%b)))
	drvS=-(n/(2*(sigma^2)))+(1/(2*(sigma^4)))*t(y-X%*%b)%*%(y-X%*%b)
    drv=c(drvB,drvS)
	return(drv)
}              

hessian=function(y,x,alpha,beta,sigma){
	n=length(y)
	X=cbind(1,x)
	b=c(alpha,beta)
	q1=-((t(X)%*%X)/sigma^2)  
	q2=-(t(X)%*%y)/((sigma^2))^2+((t(X)%*%X)%*%b)/((sigma^4))
	qa=rbind(q1,c(q2))  
	cb=(n/(2*((sigma^4))))-(1/((sqrt(sigma^6))))*((t(y-(X)%*%b))%*%(y-(X)%*%b))     
	qa=cbind(qa,c(q2,cb))            
	return(-qa)	
}
       
 
# numerical optimization using Newton's method with jacobian and hessian matrices 
# a. choose starting values        
alpha=-1
beta=.5   
sigma=(.01)^2

# b. while L'(Theta) !=0     
eps=.00001
steps=.1

counts=0     
pars=c(alpha,beta,sigma)                                                     

# want to find vector where jacobian is zero, i.e. maximized 
while(any(abs(c(jacobian(x=x,y=y,alpha=alpha,beta=beta,sigma=sigma)))>eps)){
	counts=counts+1	
	alpha=pars[1]
	beta=pars[2]
   	sigma=pars[3]
	difs=c(pars-(solve(-hessian(x=x,y=y,alpha=alpha,beta=beta,sigma=sigma))%*%jacobian(x=x,y=y,alpha=alpha,beta=beta,sigma=sigma)))
    pars=difs
 	pars[3]=-pars[3]

}   
   	
    
# 2. Numerical optimization using optim
#   - optimize using optim in R 
alpha=-1
beta=.5   
sigma=(2.1)^2

pars=c(alpha,beta,sigma)     
optim(par=pars,fn=loglike,x=x,y=y,method = c("BFGS"),hessian=T)
    


# 3. Genetic Search using genoud
#   - optimize using genoud in R   
#alpha=-1
#beta=.5   
#sigma=(0.01)^2

genoud(fn=loglike,nvars=3,max=F,pop.size=100,y=y,x=x)
       
############ Invariance on Sigma ############
# Note there is an invariance problem since we're are searching for sigma; and -sigma^2 = sigma^2
#  - often best to climb upwards from 0 since sigma > 0; but must make sure to tell the Newton function to explicitely climb  
#  - also note that bad starting points can seriously impact the optimum found 

alpha=-100
beta=500   
sigma=(200.1)^2

pars=c(alpha,beta,sigma)     
optim(par=pars,fn=loglike,x=x,y=y,method = c("BFGS"),hessian=T)
   
  

# newton   

alpha=-100
beta=500   
sigma=(200.1)^2

eps=.00001
steps=.1

counts=0     
pars=c(alpha,beta,sigma)                                                     

# want to find vector where jacobian is zero, i.e. maximized 
while(any(abs(c(jacobian(x=x,y=y,alpha=alpha,beta=beta,sigma=sigma)))>eps)){
	counts=counts+1	
	alpha=pars[1]
	beta=pars[2]
   	sigma=pars[3]
	difs=c(pars-(solve(-hessian(x=x,y=y,alpha=alpha,beta=beta,sigma=sigma))%*%jacobian(x=x,y=y,alpha=alpha,beta=beta,sigma=sigma)))
    pars=difs
 	pars[3]=pars[3]

}   
  

# END