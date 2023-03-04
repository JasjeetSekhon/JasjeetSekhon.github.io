     
cosponsorEM=function(y,party=party,eps=.1,better=FALSE){
     
	# party has to be in the direction of 1 for Dems or Left party    

    # liklihood function in dimension i=1,2,...,n
	llik.i=function(pars,y,beta,delta,lam_a,lam_g){	
		alpha=pars[1]
		gamma=pars[2]
		mu=(1/(1+exp(- ((-(alpha-beta)^2)-gamma-delta))))
		llk=(sum(y*log(mu)+(1-y)*log(1-mu))
			-lam_a*sum(alpha^2)
			-lam_g*sum(gamma^2)
			 #-lam_b*sum(beta^2)-lam_d*sum(delta^2)
			)   	
		return(-llk)
	}       

	# liklihood functions in dimension j=1,2,...,m
	llik.j=function(pars,y,alpha,gamma,lam_b,lam_d,lam_r1,lam_r2){	
		beta=pars[1]
		delta=pars[2]
		mu=(1/(1+exp(- ((-(alpha-beta)^2)-gamma-delta) )))
		llk=(sum(y*log(mu)+(1-y)*log(1-mu))
			#-lam_a*sum(alpha^2)-lam_g*sum(gamma^2)
			 -lam_b*sum(beta^2)
			 -lam_d*sum(delta^2)
			)   	
		return(-llk)
	}
	
	# function to create starting values sweeping over naive model estimates
	sweepStart=function(y=y,rest=F,party=party){
                        
		# party has to be in the direction of 1 for Dems or Left party
		bills=array(NA,ncol(y))
		for(j in 1:length(bills)){
			bills[j]=mean(party[which(y[,j]==1)])/(1/sum(y[,j]))	
		}

		members=array(NA,nrow(y))
		for(i in 1:length(members)){
			members[i]=mean(bills[which(y[i,]==1)])	
		}             

		bills[is.na(bills)]=mean(bills[!is.na(bills)])
		members[is.na(members)]=mean(members[!is.na(members)])

		# converge the sweep
		set.seed(1005)             
		members=(members-mean(bills))/sd(bills)
		bills=jitter((bills-mean(bills))/sd(bills),factor=2)  # small jitter to climb in probability
        

        #bills_old=bills  
		#members_old=members
		#difs=sum((members_old-members)^2,na.rm=T)
		difs=100
		while(difs>1){   

			for(j in 1:length(bills)){
				bills[j]=mean(members[which(y[,j]==1)])
			} 
		  	
			bills[is.na(bills)]=mean(bills[!is.na(bills)])             
			bills=jitter((bills-mean(bills))/sd(bills),factor=2)    
		
			for(i in 1:length(members)){
				members[i]=mean(bills[which(y[i,]==1)])	
			}            
			members[is.na(members)]=mean(members[!is.na(members)])
			members_old=members      
			
			difs=sum((members_old-members)^2,na.rm=T)    
		}

		if(rest==T){		
			gams=-(rowMeans(y)/rowVars(y))
			gams[which(is.na(gams))]=mean(gams[which(!is.na(gams))])    
			gams=(gams-mean(gams))/sd(gams)
			dels=-(colMeans(y)/colVars(y))
			dels=(dels-mean(dels))/sd(dels) 
			return(list("members"=members,"bills"=bills,"gams"=gams,"dels"=dels))
		}else{
			set.seed(1005)
			return(list("members"=members,"bills"=bills,"gams"=rnorm(n=nrow(y)),"dels"=rnorm(n=ncol(y))))
		}	
	}     
	
	
	
	betterStart=function(y=y,rest=F,party=party,m,n){
      	
		starts=sweepStart(y,party,rest=F)
		starts_alpha=alpha=starts$members
		starts_beta=beta=starts$bills
		starts_gamma=gamma=starts$gams
		starts_delta=delta=starts$dels
		
		# gain an advantage by splitting y in parts 
		set.seed(1005)
		# colsplit 
		col_indx=sample(c(1:5),replace=T,size=m)
		row_indx=sample(c(1:5),replace=T,size=n)
				
		iters=1
		for(j in 1:5){  
			iters=iters+1
			ins=initial(y=y[,which(col_indx==j)],
				alpha,gamma,
				beta[which(col_indx==j)],delta[which(col_indx==j)],
				lam_a,lam_g,lam_b,lam_d,
				n=n,m=length(which(col_indx==j)),party=party) 

			starts_alpha=(ins$alpha + starts_alpha) 
			starts_beta[which(col_indx==j)]=(ins$beta)
			starts_gamma=(ins$gamma + starts_gamma)
			starts_delta[which(col_indx==j)]=(ins$delta)
			
			alpha=starts_alpha/iters
			gamma=starts_gamma/iters
			delta=starts_delta
			beta=starts_beta
		}			
		
		iters=1
		for(i in 1:5){  
			iters=iters+1
			ins=initial(y=y[which(row_indx==i),],
				alpha[which(row_indx==i)],gamma[which(row_indx==i)],
				beta,delta,
				lam_a,lam_g,lam_b,lam_d,
				n=length(which(row_indx==i)),m=m,party=party[which(row_indx==i)]) 

			starts_alpha[which(row_indx==i)]=ins$alpha
			starts_beta=(ins$beta + starts_beta)
			starts_gamma[which(row_indx==i)]=ins$gamma
			starts_delta=(ins$delta + starts_delta)
			
			alpha=starts_alpha
			gamma=starts_gamma
			delta=starts_delta/iters
			beta=starts_beta/iters
		}
			
		return(list('alpha'=alpha,'beta'=beta,'delta'=delta,'gamma'=gamma))				
	}

	initial=function(y,alpha,gamma,beta,delta,lam_a,lam_g,lam_b,lam_d,n,m,party){
	  	
		alpha_old=alpha
		gamma_old=gamma
		beta_old=beta
		delta_old=delta


		for(i in 1:n){
			pars=c(alpha[i],gamma[i])
			pars1=optim(pars,llik.i,y=y[i,],beta=beta,delta=delta,lam_a=lam_a,lam_g=lam_g)$par
			alpha[i]=pars1[1]
			gamma[i]=pars1[2]     
		}                    

		if(!is.null(party)){
			if(mean(alpha[party==-1])>mean(alpha[party==1])){
				alpha=alpha*-1
			}                  
		}          	

		for(i in 1:m){
			pars=c(beta[i],delta[i])
			pars1=optim(pars,llik.j,y=y[,i],alpha=alpha,gamma=gamma,lam_b=lam_b,lam_d=lam_d)$par
			beta[i]=pars1[1]
			delta[i]=pars1[2]
		}
		return(list('alpha'=alpha,'beta'=beta,'delta'=delta,'gamma'=gamma))
	}


    converge=function(y,difs,eps,alpha,gamma,beta,delta,lam_a,lam_g,lam_b,lam_d,party){
		counts=0
		while(difs>eps){
        	print(difs)
			alpha_old=alpha
			gamma_old=gamma
			beta_old=beta
			delta_old=delta

			counts=counts+1   

			for(i in 1:n){
				pars=c(alpha[i],gamma[i])
				pars1=optim(pars,llik.i,y=y[i,],beta=beta,delta=delta,lam_a=lam_a,lam_g=lam_g)$par
				alpha[i]=pars1[1]
				gamma[i]=pars1[2]     
			}                    

			if(!is.null(party)){
				if(mean(alpha[party==-1])>mean(alpha[party==1])){
		 			alpha=alpha*-1
				}                  
			}          	

			for(i in 1:m){
				pars=c(beta[i],delta[i])
				pars1=optim(pars,llik.j,y=y[,i],alpha=alpha,gamma=gamma,lam_b=lam_b,lam_d=lam_d)$par
				beta[i]=pars1[1]
				delta[i]=pars1[2]
			}     

			difs=sum(c((alpha_old-alpha)^2))    
      
		}
	    return(list('alpha'=alpha,'beta'=beta,'delta'=delta,'gamma'=gamma,"counts"=counts,"difs"=difs))
	}       

	# 0. initialize starting values  
	n=nrow(y); m=ncol(y)
  
    lam_b=lam_a=.5
	lam_d=lam_g=0
    
    # normal priors on spatial terms; floating/diffuse priors on valence terms
 
	# 1. initialize optimization
	
	if(better==TRUE){            
		starts=betterStart(y=y,rest=F,party=party,m,n)
		alpha_old=alpha=starts$alpha
		beta_old=beta=starts$beta
		gamma_old=gamma=starts$gamma
		delta_old=delta=starts$delta
                     
	} else if(better==FALSE){
		starts=sweepStart(y,party,rest=F)
		alpha_old=alpha=starts$members
		beta_old=beta=starts$bills
		gamma_old=gamma=starts$gams
		delta_old=delta=starts$dels
	}
    
	inits=initial(y,alpha,gamma,beta,delta,lam_a,lam_g,lam_b,lam_d,n,m,party) 

    alpha=inits$alpha
    gamma=inits$gamma
    beta=inits$beta
    delta=inits$delta
     
	difs=sum(c((alpha_old-alpha)^2))
	eps=eps#.1
    #difs=100

    

    # 2. running optimization while loop until convergence
	runs=converge(y,difs,eps,alpha,gamma,beta,delta,lam_a,lam_g,lam_b,lam_d,party)
	   	                             
	return(list('alpha'=runs$alpha,'beta'=runs$beta,'delta'=runs$delta,'gamma'=runs$gamma))
}    
          


