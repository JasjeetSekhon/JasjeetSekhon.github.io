###############################################
# Functions and MLE
#
# John Henderson
#   
# PS 236A/ Stat 239A Fall 2012
###############################################
 

# 1. function writing & optimization for MLE for OLS and logit regression 

# all functions have the following syntax
# no fun1 can be accessed in the remaining R sesseion

fun1<-function(stuff){
	# run some computations    
	new_stuff<-stuff+stuff

	# return statement assigns an object as output
	return(new_stuff)	
}

    
fun2<-function(stuff,...){
	# ... lets you pass through objects to other internal functions
	
	# you can write internal functions, only available to fun2
	rowmean<-function(data){
		rmean_out=apply(data,MARGIN=1,FUN=mean)
		return(rmean_out)
	}	
				                       
	add_more=sum(c(rmean_out))
	new_stuff<-stuff+add_more
     
	# return statement assigns an object as output
	# return things as lists, with variable names   
	return(list("a"=new_stuff,"b"=new_stuff^2))
}          



mle_ols<-function(y,x){
	llik<-function(y,x,a,b,s=1){
        resid=(y-a-x*b)/s
        llk=log(((1/sqrt(2*pi))*exp(-((resid)^2)/2))/s)
	}
                                       
	lower_grids<-function(x,y){ 	   
		eps=1e-4      
		eps_new=10
		eps_old=5
		b=10
		while(abs(eps_new-eps_old)>eps){
			eps_old=eps_new
			a=mean(y)-b*mean(x)
            b=b-.001
			eps_new=sum(llik(a=a,b=b,x=x,y=y))
			
		}            
		a=mean(y)-b*mean(x)      
		return(a,b)
	}
	
	return(c(a,b))
}

lm(y~x)
mle_ols(y,x)
                 
# write your own mle function using optim()


 
