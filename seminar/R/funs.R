#  John Henderson 
#    <jahenderson@berkeley.edu> 
#    <http://jahenderson.com>
#  Amplify: funs.R                    
#  Version: August 3, 2011         
                   
# funs.R: list of called functions by amplify

# qSolve is a function that solves a specific quadratic equation to conduct sensitivity analysis
#   used in: ampFisher, ampMantel, ampRanksum, ampIVPermute
qSolve=function(lambda,n,m,r){
	if(lambda!=1){
		a=lambda-1
		b=-((lambda-1)*(m+r)+n)
		c=lambda*r*m

		E=c((-b+sqrt(b^2-4*a*c))/(2*a),(-b-sqrt(b^2-4*a*c))/(2*a))
		E=E[which(max(0,r+m-n) <= round(E,digits=10) & round(E,digits=10) <= min(r,m))]
	} else if(lambda==1){
		E=(r*m)/n
	}
	V=1/(1/(E)+1/(r-E)+1/(m-E)+1/(n-r-m+E))

	return(list(E=E,V=V))
} # END qSolve 
     

# countTable is a function that creates the necessary summary 2x2 table from a vector of treated and control outcomes
#   used in: ampFisher, ampMantel
countTable=function(treated,control,strata_t=NULL,strata_c=NULL){

	if(is.null(strata_t) & is.null(strata_c)){
		u_strata=1
	} else {
		strata=c(strata_t,strata_c)
		u_strata=sort(unique(strata))
	}

	if(length(u_strata)==1){
		n=length(c(treated,control))
		m=length(treated)
		r=length(which(c(treated,control)==1))
		t=length(which(treated==1))
	} else if(length(u_strata)>1){
		n=m=r=t=c()
		for(k in 1:length(u_strata)){
			n=c(n,length(c(treated,control)[which(strata==u_strata[k])]))
			m=c(m,length(treated[which(strata[1:length(treated)]==u_strata[k])]))  
			r=c(r,length(which(c(treated,control)[which(strata==u_strata[k])]==1)))
			t=c(t,length(which(treated[which(strata[1:length(treated)]==u_strata[k])]==1)))
		}
	}
	return(invisible(list(n=n,m=m,r=r,t=t)))
} # END countTable  


# arrayTable is a function that creates the a table from summary data of treated and control outcomes in a 2x2 table 
#   used in: ampMantel
arrayTable=function(n=NULL,m=NULL,r=NULL,t=NULL,out.type='tableArray'){

	if(any(is.null(n) | is.null(m) | is.null(r) | is.null(t))){
		stop("Need all input vectors (n, m, r, t) for a 2x2xS contingency table.")
	}
	if(length(n)!=length(m) | length(n)!=length(r) | length(n)!=length(t)){
		stop("Data vectors (n, m, r, t) are not all of the same length.")
	}	

	s=length(n)
	caten=c()
	tableList=tableArray=list()

	for(k in 1:s){
		tableList[[k]]=rbind(c(t[k],r[k]-t[k]),c(m[k]-t[k],t[k]-m[k]-r[k]+n[k]))	
		rownames(tableList[[k]])=c('R1','R0')
		colnames(tableList[[k]])=c('T1','T0')
		caten=c(caten,tableList[[k]][1,1],tableList[[k]][1,2],tableList[[k]][2,1],tableList[[k]][2,2])
	}

	tableArray=array(caten, 
					 dim=c(2,2,s),
					 dimnames=list(Tr=c('1','0'),R=c('1','0'),Strata=c(as.character(1:s)))
					 )

	if(out.type=='tableArray'){
		return(tableArray)
	} else if(out.type=='tableList'){
		return(tableList)
	}
} #END arrayTable   
    

# dataList is a function that creates a list of data from an array table
#   used in: ampMantel
dataList=function(tableArray,s=NULL){

	caten=c(tableArray)
	s=length(caten)/4
	listData=matrix(NA,4,s)
	for(k in 1:s){
		listData[1,k]=sum(caten[(1+4*(k-1)):(4+4*(k-1))])
		listData[2,k]=sum(caten[c(1,3)+4*(k-1)])
		listData[3,k]=sum(caten[c(1,2)+4*(k-1)])
		listData[4,k]=caten[1+4*(k-1)]
	}
	colnames(listData)=c(as.character(1:s))
	rownames(listData)=c('n','m','r','t')
	return(listData)
} #END dataList


# uSolve is a function that computes statistics to compute approximate pvalues for unmatched data with non-binary outcomes
#   used in: ampRankSum, ampIVPermute
uSolve=function(E=E,V=V,n=n,m=m,k=k,qsort=qsort,ties=F){
	# recall that qsort is the sorted ranks of the catenated treated and control vectors of responses 
	if(ties==F){
		if(k==0){
			q1=0
			w1=0
			q0=sum(qsort[(k+1):n])/(n-k)
			w0=sum((qsort[(k+1):n]-q0)^2)/(n-k-1)
			mu=sum(qsort)*m/n
			sigma=sqrt((sum((qsort-sum(qsort)/n)^2))*(m*(n-m))/(n*(n-1)))
		} else if(k==n){
			q1=sum(qsort[1:k])/k
			w1=sum((qsort[1:k]-q1)^2)/(k-1)
			q0=0
			w0=0
			mu=sum(qsort)*m/n
			sigma=sqrt((sum((qsort-sum(qsort)/n)^2))*(m*(n-m))/(n*(n-1)))
		} else if(k==1){
			q1=sum(qsort[1:k])/k
			w1=0
			q0=sum(qsort[(k+1):n])/(n-k)
			w0=sum((qsort[(k+1):n]-q0)^2)/(n-k-1)
			mu=E*q1+(m-E)*q0
			sigma=sqrt((w1-w0)*E-(E^2+V)*((w1)/(k)+(w0)/(n-k))+(m*(n-k-m+2*E)*w0)/(n-k)+V*((q1-q0)^2))
		} else if(k==(n-1)){
			q1=sum(qsort[1:k])/k
			w1=sum((qsort[1:k]-q1)^2)/(k-1)
			q0=sum(qsort[(k+1):n])/(n-k)
			w0=0
			mu=E*q1+(m-E)*q0
			sigma=sqrt((w1-w0)*E-(E^2+V)*((w1)/(k)+(w0)/(n-k))+(m*(n-k-m+2*E)*w0)/(n-k)+V*((q1-q0)^2))
		} else if(k>1 & k<(n-1)){
			q1=(sum(qsort[1:k])/k)
			q0=(sum(qsort[(k+1):n])/(n-k))
			w1=sum((qsort[1:k]-q1)^2)/(k-1)
			w0=sum((qsort[(k+1):n]-q0)^2)/(n-k-1)
			mu=E*q1+(m-E)*q0
			sigma=sqrt((w1-w0)*E-(E^2+V)*((w1)/(k)+(w0)/(n-k))+(m*(n-k-m+2*E)*w0)/(n-k)+V*((q1-q0)^2))
		}
	} else if(ties==T){
		mu=(E*n+m*(n-k+1))/2
		sigma=sqrt(((((n+1)*(2*k-n)+2*m*(n-k+1))*E-(n+2)*E^2)/12 + (m*(n-k-m)*(n-k+1))/12 + (V*(n-1)*(n+2/3))/4))
	}
	return(list(mu=mu,sigma=sigma))
} #END uSolve	
  

# hlSolve is a function that computes Hodges-Lehmann estimates for matched and unmatched non-binary data
#   used in: ampRankSum
hlSolve=function(R,Z,t=NULL,matched=F,inc=.005){
	
	m=length(which(Z==1))
	n=length(R)
	invisible(library(stringr))  
	dg=abs(str_length((inc))-3)
	
	#for now, 1 to 1 matching only; 
	#matched == F is rank sum statistic 
	#matched == T is sum sign rank statistic 
	
	if(matched==F & is.null(t)){
		t=m*(n+1)/2
	}
	if(matched==T & is.null(t)){
		t=m*(n/2+1)/4
	}
	
	p_tau=m_tau=0
	np0=np1=nm0=nm1=1
	
	while(sign(np0)==sign(np1) & sign(nm0)==sign(nm1)){
		if(matched==F){
			np0=sum(rank((R-p_tau*Z))[which(Z==1)]) - t
			p_tau=p_tau+inc
			np1=sum(rank((R-p_tau*Z))[which(Z==1)]) - t
			
			nm0=sum(rank((R-m_tau*Z))[which(Z==1)]) - t
			m_tau=m_tau-inc
			nm1=sum(rank((R-m_tau*Z))[which(Z==1)]) - t
			
		} else if(matched==T){
			np0=sum(rank(abs(((R-p_tau*Z)[Z==1])-((R-p_tau*Z)[Z==0])))[which(((R-p_tau*Z)[Z==1])>((R-p_tau*Z)[Z==0]))]) - t
			p_tau=p_tau+inc
			np1=sum(rank(abs(((R-p_tau*Z)[Z==1])-((R-p_tau*Z)[Z==0])))[which(((R-p_tau*Z)[Z==1])>((R-p_tau*Z)[Z==0]))]) - t
			
			nm0=sum(rank(abs(((R-m_tau*Z)[Z==1])-((R-m_tau*Z)[Z==0])))[which(((R-m_tau*Z)[Z==1])>((R-m_tau*Z)[Z==0]))]) - t
			m_tau=m_tau-inc
			nm1=sum(rank(abs(((R-m_tau*Z)[Z==1])-((R-m_tau*Z)[Z==0])))[which(((R-m_tau*Z)[Z==1])>((R-m_tau*Z)[Z==0]))]) - t
		}
	}
	
	if(abs(nm1) + abs(nm0) < abs(np1) + abs(np0)){
		tau=-median(c(-c(m_tau,m_tau+inc),c(p_tau-inc,p_tau)))
	} else if(abs(nm1) + abs(nm0) > abs(np1) + abs(np0)){
		tau=median(c(-c(m_tau,m_tau+inc),c(p_tau-inc,p_tau)))
	} else{
		tau=0
	}
	
	return(tau)	
} #END hlSolve
       

# ivSolve solves for an additive treatment effect in an IV setup 
#   used in: ampIVPermute
ivSolve=function(R,Z,X,t=NULL,matched=F,inc=inc,W=NULL){
	inc_ini=inc
	m=length(which(Z==1))
	n=length(R)
	
	if(matched==F & is.null(t)){
		t=m*(n+1)/2
	}
	if(matched==T & is.null(t)){
		t=m*(n/2+1)/4
	}
	
	p_tau=m_tau=0
	np0=np1=nm0=nm1=1
	count=0
	while(sign(np0)==sign(np1) & sign(nm0)==sign(nm1)){
		count=count+1
		if(matched==F){
			np0=sum((rank(R-p_tau*X))[which(Z==1)]) - t
			p_tau=p_tau+inc
			np1=sum((rank(R-p_tau*X))[which(Z==1)]) - t
			
			nm0=sum((rank(R-m_tau*X))[which(Z==1)]) - t
			m_tau=m_tau-inc
			nm1=sum((rank(R-m_tau*X))[which(Z==1)]) - t
			
		} else if(matched==T){
			if(is.null(W)){
				np0=sum(rank(abs(((R-p_tau*X)[Z==1])-((R-p_tau*X)[Z==0])))[which(((R-p_tau*X)[Z==1])>((R-p_tau*X)[Z==0]))]) - t
				p_tau=p_tau+inc
				np1=sum(rank(abs(((R-p_tau*X)[Z==1])-((R-p_tau*X)[Z==0])))[which(((R-p_tau*X)[Z==1])>((R-p_tau*X)[Z==0]))]) - t
			
				nm0=sum(rank(abs(((R-m_tau*X)[Z==1])-((R-m_tau*X)[Z==0])))[which(((R-m_tau*X)[Z==1])>((R-m_tau*X)[Z==0]))]) - t
				m_tau=m_tau-inc
				nm1=sum(rank(abs(((R-m_tau*X)[Z==1])-((R-m_tau*X)[Z==0])))[which(((R-m_tau*X)[Z==1])>((R-m_tau*X)[Z==0]))]) - t
			} else if(!is.null(W)){
				np0=sum(wtd.rank(normwt=F,abs(((R-p_tau*X)[Z==1])-((R-p_tau*X)[Z==0])),weights=W[1:(length(W)/2)])[which(((R-p_tau*X)[Z==1])>((R-p_tau*X)[Z==0]))]) - t
				p_tau=p_tau+inc
				np1=sum(wtd.rank(normwt=F,abs(((R-p_tau*X)[Z==1])-((R-p_tau*X)[Z==0])),weights=W[1:(length(W)/2)])[which(((R-p_tau*X)[Z==1])>((R-p_tau*X)[Z==0]))]) - t
				
				nm0=sum(wtd.rank(normwt=F,abs(((R-m_tau*X)[Z==1])-((R-m_tau*X)[Z==0])),weights=W[1:(length(W)/2)])[which(((R-m_tau*X)[Z==1])>((R-m_tau*X)[Z==0]))]) - t
				m_tau=m_tau-inc
				nm1=sum(wtd.rank(normwt=F,abs(((R-m_tau*X)[Z==1])-((R-m_tau*X)[Z==0])),weights=W[1:(length(W)/2)])[which(((R-m_tau*X)[Z==1])>((R-m_tau*X)[Z==0]))]) - t
				
			}
		}
		if(count==500 | count==1000 | count==1500 | count==2000){
			inc=10*inc
		}
		if(count>4000){
			warning("Estimate not converging. Stopping run.")
			break				
		}
	}
	inc=inc_ini
	if(abs(nm1) + abs(nm0) < abs(np1) + abs(np0)){
		tau=-median(c(-c(m_tau,m_tau+inc),c(p_tau-inc,p_tau)))
	} else if(abs(nm1) + abs(nm0) > abs(np1) + abs(np0)){
		tau=median(c(-c(m_tau,m_tau+inc),c(p_tau-inc,p_tau)))
	} else{
		tau=0
	}
	
	return(tau)	
} #END ivSolve 


# permInf permutes rank sum test statistics under a null of H_0: B = B_0   
#   used in: ampIVPermute
permInf=function(R,X,Z,beta_iv,t_stat,matched,nsamp){
	set.seed(1005)
	tout=array(NA,nsamp)
	for(k in 1:nsamp){
		zt=sample(Z)
		if(matched==F){
			tout[k]=sum(rank(R-beta_iv*X)[zt==1])
		} else if(matched==T){
			tout[k]=sum((rank(abs((R-beta_iv*X)[zt==1] - (R-beta_iv*X)[zt==0])))[which((R-beta_iv*X)[zt==1]>(R-beta_iv*X)[zt==0])])
		}
	}
	return(tout)
} #END permInf

#END