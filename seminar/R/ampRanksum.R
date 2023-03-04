#  John Henderson 
#    <jahenderson@berkeley.edu> 
#    <http://jahenderson.com>
#  Amplify: ampRanksum.R                    
#  Version: August 3, 2011

# ampRanksum

# The function 'ampRanksum' implements a two paramter sensitivity analysis for unmatched groups with binary treatments 
#  and non-binary (continuous or categorical) responses.  This sensitivity analysis is based on an amplification approach 
#  of the one parameter model developed in Rosenbaum and Silber (2009) and Rosenbaum (2002).  The function uses Wilcoxons 
#  rank sum statistic, and allows the researcher to report either probability values or Hodges-Lehmann point estimates of an 
#  additive treatment effect (See Rosenbaum (2002) for more details).
#
# 'ampRanksum' takes five arguments: 
#   treated:  a vector of non-binary responses for treated units 
#   control:  a vector of non-binary responses for control units 
#   Gamma:  vector of sensitivity values to check for the association between U and treatment 
#   Delta:  vector of sensitivity values to check for the association between U and response 
#   digits:  number of rounding digits for pvalues
#
# 'ampRanksum' reports a table of (approximate) probability values or Hodges-Lehmann (HL) point estimates for a treatment 
#  effect at each level of Gamma and Delta.  Note that the HL point estimates take quite a bit longer to compute.  

ampRanksum=function(treated=NULL,control=NULL,Gamma=c(1:5),Delta=c(1:5),digits=5){
	               
	if(is.null(treated) | is.null(control)){ 	 
		stop("Missing input data.")	
	}                                            
	if(any(is.na(c(treated,control)))){ 		 
		stop("Missing data in treated or control.")
	}  
	if(length(unique(c(treated,control)))==2){
	 	stop("Outomes appear to be binary.")
	} 
	if(length(treated)==length(control)){ 	 
		warning("Length of treated and control are equal, make sure these are unmatched.")
	}     
	
	if(mean(treated)>=mean(control)){
	 	tr=treated
   		ct=control
		sgn=1
	} else if(mean(treated)<mean(control)){
	 	tr=control
		ct=treated
		sgn=-1   
	}
	
	q=rank(c(tr,ct))
	t=sum(q[1:length(tr)])
	m=length(tr)
	n=length(q)
	qsort=sort(q,decreasing=TRUE)
	
	R=c(tr,ct)
	Z=c(rep(1,length(tr)),rep(0,length(ct)))
	
	pvals=hl=lambda=matrix(0,length(Gamma),length(Delta)) 
	wilx=wilcox.test(x=tr,y=ct,exact=F,alternative='greater')$p.value
	
	for(i in 1:length(Gamma)){ 
		for(j in 1:length(Delta)){ 
			
			lambda[i,j]=(Delta[j]*Gamma[i]+1)/(Delta[j]+Gamma[i])
			
			if(lambda[i,j]==1){
					pvals[i,j]=wilx
					hl[i,j]=hlSolve(R=R,Z=Z,t=NULL,matched=F,inc=.005)*sgn
			} else {
				mu=sigma=array(NA,n+1)
				for(k in 0:n){
					qS=qSolve(lambda[i,j],n=n,m=m,r=k)
					uS=uSolve(E=qS$E,V=qS$V,n=n,m=m,k=k,qsort=qsort)
					mu[k+1]=uS$mu
					sigma[k+1]=uS$sigma
				} #k-loop
				
				if(lambda[i,j]>1){
					t_max=max(mu)
						pvals[i,j]=1-pnorm(min((t-mu)/sigma))
				} else if(lambda[i,j]<1){
					t_max=min(mu)		
						pvals[i,j]=1-pnorm(max((t-mu)/sigma))
				}
				hl[i,j]=hlSolve(R=R,Z=Z,t=t_max,matched=F,inc=.005)*sgn
			} # else lambda 
		} # j-loop
	} # i-loop
	
	options(scipen=3)
	colnames(pvals)=colnames(hl)=as.character(round(Delta[1:length(Delta)],digits=3))
	rownames(pvals)=rownames(hl)=as.character(round(Gamma[1:length(Gamma)],digits=3))
	mu=mean(treated)-mean(control)
	return(list(mu=mu,Gamma=Gamma,Delta=Delta,hl=hl,pvals=round(as.data.frame(pvals),digits=digits)))  
 	
} 

#END ampRanksum