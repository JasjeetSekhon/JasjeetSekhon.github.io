rm(list = ls())
library(stringr)
library(MASS)
library(Matching)
       
source(url("http://sekhon.berkeley.edu/causalinf/R/balStatPlot.R"))

# problem 3
# (a)

# load water for life data
water=read.csv(file=url("http://sekhon.berkeley.edu/causalinf/data/cross_section_wfl.csv"))

treat=water$treat
covars=water[,which(str_sub(names(water),1,5)=='y1990')]

#exclude covars with less than 3 unique values
indm=array(TRUE,ncol(covars))
for(j in 1:ncol(covars)){
sun=sort(unique(covars[,j]))
	if(length(sun)<4 & sum(sun)!=1){
		indm[j]=FALSE
	}
}

covars=as.matrix(covars[,indm])     
k=dim(covars)[2]
covars=cbind(covars,covars[,1]*covars[,3],covars[,1]*covars[,10],covars[,3]*covars[,10],covars[,1]*covars[,3]*covars[,10],covars[,12]^3,covars[,21]^2)                                
colnames(covars)[(k+1):dim(covars)[2]]=c(paste('v',as.character(1:length((k+1):dim(covars)[2])),sep=''))
             

mortality=water$y1999.tasatot
infect=water$y1999.tminfec
          
water.data=as.data.frame(cbind(treat,covars,mortality,infect))
          

# write a MH matching function

mahalMatching=function(treat,covars,estimand='ATT'){

	# contruct a MH distance matrix 
	S=ginv(cov(covars)) #solve(t(as.matrix(covars))%*%as.matrix(covars))	
	D=matrix(0,nrow(covars),nrow(covars))
	for(i in 1:nrow(covars)){ 
		D[i,]=mahalanobis(covars, as.numeric(c(covars[i,])), cov=S, inverted=TRUE)
	}	  	
	  
	# match on D distance  
		 # matching w replacement; ties broken by deterministic weights; 1 to 1, etc
	indx.control=which(treat==0)
	indx.treated=which(treat==1)	

	Dtr=D[indx.treated,indx.control]
	Dct=D[indx.control,indx.treated]  

	wgts=indx.matched.tr=indx.matched.ct=list()		
		
	if(estimand=='ATT'){
		for(i in 1:nrow(Dtr)){
			indx.matched.ct[[i]]=indx.control[which(Dtr[i,]==min(Dtr[i,]))]
			indx.matched.tr[[i]]=array(indx.treated[i],length(which(Dtr[i,]==min(Dtr[i,]))))
			wgts[[i]]=array(1/length(which(Dtr[i,]==min(Dtr[i,]))),length(which(Dtr[i,]==min(Dtr[i,]))))      	
		}		
	} else if(estimand=='ATC'){
		for(i in 1:nrow(Dct)){
			indx.matched.tr[[i]]=indx.treated[which(Dct[i,]==min(Dct[i,]))]
			indx.matched.ct[[i]]=array(indx.control[i],length(which(Dct[i,]==min(Dct[i,]))))
			wgts[[i]]=array(1/length(which(Dct[i,]==min(Dct[i,]))),length(which(Dct[i,]==min(Dct[i,]))))      	
		}		
	} else if(estimand=='ATE'){
		for(i in 1:nrow(Dtr)){
		   	indx.matched.ct[[i]]=indx.control[which(Dtr[i,]==min(Dtr[i,]))]
		   	indx.matched.tr[[i]]=array(indx.treated[i],length(which(Dtr[i,]==min(Dtr[i,]))))
		   	wgts[[i]]=array(1/length(which(Dtr[i,]==min(Dtr[i,]))),length(which(Dtr[i,]==min(Dtr[i,]))))      	
		}
		n=i  
		for(i in 1:nrow(Dct)){
			indx.matched.tr[[n+i]]=indx.treated[which(Dct[i,]==min(Dct[i,]))]
			indx.matched.ct[[n+i]]=array(indx.control[i],length(which(Dct[i,]==min(Dct[i,]))))
			wgts[[n+i]]=array(1/length(which(Dct[i,]==min(Dct[i,]))),length(which(Dct[i,]==min(Dct[i,]))))      	
		}	   		
	}  
	 
	indx.matched.ct=unlist(indx.matched.ct)
	indx.matched.tr=unlist(indx.matched.tr)
	weights=unlist(wgts)
	
	return(list('index.treated'=indx.matched.tr,'index.control'=indx.matched.ct,'weights'=weights))  	
}                         

mout=mahalMatching(treat,covars,estimand='ATT')

# check using Match
mout_a=Match(Tr=treat,X=covars,Weight=2) 

# not the same ... ? 
# check with mout1

mb_a=MatchBalance(treat~covars,match.out=mout_a,nboots=500)
     
# plot balance
plot.pval(mb_a,colnames(covars))


# part (b)
# set seeds
set.seed=1005
unif.seed=105
int.seed=10005

# match to utilize domains and starting values options here         
    
                                        
genout_b=GenMatch(Tr=treat,X=covars,BalanceMatrix=covars,pop.size=100,
	hard.generation.limit=T,wait.generations=50,max.generations=100,
	unif.seed=unif.seed,int.seed=int.seed)
 
mout_b=Match(Tr=treat,X=covars,Weight.matrix=genout_b) 
mb_b=MatchBalance(treat~covars,match.out=mout_b,nboots=500) 
balouts_b=plot.pval(mb_b,colnames(covars))                                       
                                                       

# identify most imbalanced on standardized differences 
sds=array(NA,ncol(covars))
for(j in 1:ncol(covars)){
	sds[j]=sd(covars[,j])
}             
   
stddif=(balouts_b[,1]-balouts_b[,2])/sds
qnt=quantile(abs(stddif),prob=.75)       
  
# upweight imbalanced using domain restrictions
domains=cbind(array(0,ncol(covars)),array(1000,ncol(covars)))      
domains[which(abs(stddif)>qnt),1]=250
domains[which(abs(stddif)<qnt),2]=750
            
# restart *close* to where genout_b ended
starting.values=diag(genout_b$Weight.matrix)
           
genout_bb=GenMatch(Tr=treat,X=covars,BalanceMatrix=covars,pop.size=100,
	hard.generation.limit=T,wait.generations=50,max.generations=100,
	starting.values=starting.values,Domains=domains,
	unif.seed=unif.seed,int.seed=int.seed)
 
mout_bb=Match(Tr=treat,X=covars,Weight.matrix=genout_bb) 
mb_bb=MatchBalance(treat~covars,match.out=mout_bb,nboots=500) 
balouts_bb=plot.pval(mb_bb,colnames(covars))
   
pre_mortality=water$y1990.tasatot
pre_infect=water$y1990.tminfec
       
# pre-privatiziation total mortality                                                                           
# before
qqplot(pre_mortality[which(treat==1)],pre_mortality[which(treat==0)],xlim=c(0,25),ylim=c(0,25))
abline(a=0,b=1,lty=2)

# after                                 
qqplot(pre_mortality[mout_bb$index.treated],pre_mortality[mout_bb$index.control],xlim=c(0,25),ylim=c(0,25))
abline(a=0,b=1,lty=2)     
                         

# pre-privatiziation water-bourne disease mortality                                                                           
# before
qqplot(pre_infect[which(treat==1)],pre_infect[which(treat==0)],xlim=c(0,25),ylim=c(0,25))
abline(a=0,b=1,lty=2)

# after                                 
qqplot(pre_infect[mout_bb$index.treated],pre_infect[mout_bb$index.control],xlim=c(0,25),ylim=c(0,25))
abline(a=0,b=1,lty=2)
  
# part (c)
ate=Match(Y=mortality,Tr=treat,X=covars,Weight.matrix=genout_bb)$est
summary(Match(Y=mortality,Tr=treat,X=covars,Weight.matrix=genout_bb))
#ate=(Match(Y=infect,Tr=treat,X=covars,Weight.matrix=genout_bb))
  
# plot unit treatment effects
plot(density(from=min(mortality[mout_bb$index.treated]-mortality[mout_bb$index.control]),
	to=max(mortality[mout_bb$index.treated]-mortality[mout_bb$index.control]),
	mortality[mout_bb$index.treated]-mortality[mout_bb$index.control]),
	xlab='Treatment Effect',main='Unit Treatment Estimate Density')
        
# part (d) 
# caliper match to drop 10% of treated units
           
# while loop to find a caliper vector that targets above imbalanced covariates, all
#  at the approximately same level
caliper=array(5,ncol(covars))   
dropped=0                            
while(dropped<.10){  
	caliper[which(abs(stddif)>qnt)]=caliper[which(abs(stddif)>qnt)]-.1

	genout_d=GenMatch(Tr=treat,X=covars,BalanceMatrix=covars,pop.size=10,
		hard.generation.limit=T,wait.generations=5,max.generations=10,
		unif.seed=unif.seed,int.seed=int.seed,caliper=caliper)

	mout_d=Match(Tr=treat,X=covars,Weight.matrix=genout_d,caliper=caliper)  	  
	
	dropped=length(mout_d$index.dropped)/length(mout_d$index.treated) 
	print(dropped)
}

caliper[which(abs(stddif)>qnt)]=caliper[which(abs(stddif)>qnt)]+.1         

genout_d=GenMatch(Tr=treat,X=covars,BalanceMatrix=covars,pop.size=100,
	hard.generation.limit=T,wait.generations=50,max.generations=100,
	unif.seed=unif.seed,int.seed=int.seed,caliper=caliper)
 
mout_d=Match(Tr=treat,X=covars,Weight.matrix=genout_d,caliper=caliper) 
mb_d=MatchBalance(treat~covars,match.out=mout_d,nboots=500) 
balouts_d=plot.pval(mb_d,colnames(covars))

dropped=length(mout_d$index.dropped)/length(mout_d$index.treated) 
print(dropped)

            
# part (e)
     
# maximize proportion balanced at p>.1
custom.loss=function(x,){
	pval=x
	p.1=length(which(pval>.1))/length(pval)  
	print(p.1)
	return(p.1)
}
   
genout_e=GenMatch(Tr=treat,X=covars,BalanceMatrix=covars,pop.size=100,
	hard.generation.limit=T,wait.generations=50,max.generations=100,
	unif.seed=unif.seed,int.seed=int.seed,caliper=caliper,loss=custom.loss)
 
mout_e=Match(Tr=treat,X=covars,Weight.matrix=genout_e,caliper=caliper) 
mb_e=MatchBalance(treat~covars,match.out=mout_e,nboots=500) 
balouts_e=plot.pval(mb_e,colnames(covars))

dropped=length(mout_e$index.dropped)/length(mout_e$index.treated) 
print(dropped)
 
ind=which(colnames(covars)=='y1990.tasatot' | colnames(covars)=='y1990.tminfec')

#ind=c(ind,ind+dim(covars)[2])  



# maximize proportion balanced at p>.1; reject matches if 
#  previous outcome is imbalanced at p<.1 
custom.loss=function(matches,BM){
	
	index.treated=matches[,1]
    index.control=matches[,2]   
    nvars=ncol(BM)   
	pvals=c(rep(NA,2*nvars))
	        
    for (i in 1:nvars){	            
		foo1 = t.test(BM[,i][index.treated], BM[,i][index.control])$p.value
		foo2 = ks.test(jitter(BM[,i][index.treated]), jitter(BM[,i][index.control]))$p.value
		pvals[i]=foo1
		pvals[i+nvars]=foo2
	}
	 
	if(any(pvals[ind]<.1)){
		p.1=-1e100
		print(p.1)
		print(pvals[ind])
		return(p.1)	
	} else{
		p.1=length(which(pvals>.1))/length(pvals) 
		print(p.1)
		print(pvals[ind])
		return(p.1)	
	} 


}

genout_ee=GenMatch(Tr=treat,X=covars,BalanceMatrix=covars,pop.size=100,
	hard.generation.limit=T,wait.generations=50,max.generations=100,
	unif.seed=unif.seed,int.seed=int.seed,caliper=caliper,fit.func=custom.loss)
	
mout_ee=Match(Tr=treat,X=covars,Weight.matrix=genout_ee,caliper=caliper) 
mb_ee=MatchBalance(treat~covars,match.out=mout_ee,nboots=500) 
balouts_ee=plot.pval(mb_ee,colnames(covars))

dropped=length(mout_ee$index.dropped)/length(mout_ee$index.treated) 
print(dropped)

# part (f) 
     
# END        