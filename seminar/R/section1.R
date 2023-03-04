###############################################
# Section 1: Observed Sensitivity Analysis: Before and After Matching
#
# John Henderson         
#
# PS 236B Spring 2013
###############################################

rm(list=ls(all=TRUE))
#Change path to correct directory

load(file=url("http://sekhon.berkeley.edu/seminar/data/section1.Rdata"))     

attach(data_rep)
attach(factors)    

# Repilcating Kam and Palmer (2009) "Reconsidering the Effects of Education on Political Participation"
          
# Kam and Palmer match 1 to 3 on a propensity score where there is serious lack of common support   

model=formula(college ~ yPubAffZ + yNewspaperZ + yRadioZ + yMagazineZ + yFamTalkZ + yFrTalkZ + yAdultTalkZ + 
				        ySPIDZ + yGovtOpinionZ + yGovtCrookZ + yGovtWasteZ + yTrGovtZ + yGovtSmartZ + yGovt4AllZ + 
						yLifeWishZ + yGLuckZ + yFPlansZ + yWinArgZ + yStrOpinionZ + yMChangeZ + yTrOthersZ + 
						yOthHelpZ + yOthFairZ + yKnowledgeZ + yNextSchZ + yGPAZ + ySchOfficerZ + ySchPublishZ + 
						yHobbyZ + ySchClubZ + yOccClubZ + yNeighClubZ + yRelClubZ + yYouthOrgZ + yClubLevZ + 
						yPhoneZ + yGenZ + yRaceZ + pNewspaperZ + pRadioZ + pTVZ + pMagazineZ + pLifeWishZ + pGLuckZ + 
						pFPlansZ + pWinArgZ + pStrOpinionZ + pMChangeZ + pTrOthersZ + pOthHelpZ + pOthFairZ + pSPIDZ + 
						pVoteZ + pPersuadeZ + pRallyZ + pOthActZ + pPolClubZ + pButtonZ + pMoneyZ + pGovtOpinionZ + 
						pGovtCrookZ + pGovtWasteZ + pTrGovtZ + pGovtSmartZ + pGovt4AllZ + pEmployZ + pEducHHZ + 
						pChurchOrgZ + pFratOrgZ + pProOrgZ + pCivicOrgZ + pCLOrgZ + pNeighClubZ + pSportClubZ + 
						pInfClubZ + pFarmGrZ + pWomenClubZ + pClubLevZ + pHHIncZ + pOwnHomeZ + pKnowledgeZ)

pscore=glm(model,family=binomial(link=logit),data=factors)
etahat=pscore$fitted.values   
                           
plot(density(etahat[college==1],from=min(etahat[college==1]),to=max(etahat[college==1])))
lines(density(etahat[college==0],from=min(etahat[college==0]),to=max(etahat[college==0])),lty=2)         

# 1. check balance using IPW and bootstrapping on 

Xbal=cbind(yGPA,yGen,yBlack,yRep,yKnowledge,yNextSch,
		   pVote,pPersuade,pParticipate2,pEmploy,pEducHH,
		   pEducW,pHHInc,pOwnHome,pRep,pKnowledge)
                                                       
pvals=list()
pvals[[1]]=pvals[[2]]=array(0,ncol(Xbal))

samps=matrix(0,1000,2)  
p=sum(college)/length(college) 
set.seed(1005)

for(j in 1:ncol(Xbal)){
	for(i in 1:nrow(samps)){ 
		indx=sample(1:length(college),size=length(college),replace=T)
		samps[i,1]=mean((college[indx]*Xbal[indx,j])/etahat[indx]-((1-college[indx])*Xbal[indx,j])/(1-etahat[indx]))
		samps[i,2]=mean((college[indx]*Xbal[indx,j])/.5-((1-college[indx])*Xbal[indx,j])/(1-.5))
		#mean(((college[indx]-etahat[indx])*Xbal[indx,j])/(p*(1-etahat[indx])))
	
	}         
	pvals[[1]][j]=1-pnorm(abs(mean(samps[,1])/sd(samps[,1])))       # ATE after balance
	pvals[[2]][j]=1-pnorm(abs(mean(samps[,2])/sd(samps[,2])))       # ATE before balance
}
          

plot(pvals[[1]],ylim=c(0,1))
points(pvals[[2]],ylim=c(0,1),col='red')
abline(h=.1,lty=2)      

# IF we had coverage in the data; COULD potentially use pscore to condition on via MATCHING                

# CAN still estimate using IWP; 
   # - but these balance tests and results could be sensitive to weights near 1 or 0

samps=matrix(0,1000,2)      
set.seed(1005)

for(i in 1:nrow(samps)){ 
	indx=sample(1:length(college),size=length(college),replace=T)
	samps[i,1]=mean((college[indx]*yppnscal[indx])/etahat[indx]-((1-college[indx])*yppnscal[indx])/(1-etahat[indx]))
	samps[i,2]=mean((college[indx]*yppnscal[indx])/.5-((1-college[indx])*yppnscal[indx])/(1-.5))
	#mean(((college[indx]-etahat[indx])*Xbal[indx,j])/(p*(1-etahat[indx])))

}    
                                           
1-pnorm(abs(mean(samps[,1])/sd(samps[,1])))
1-pnorm(abs(mean(samps[,2])/sd(samps[,2])))     

                                                                                         
# Alternatively use information in the propensity score to appraise the potential for confounding 
  
# Gamma 
Gamma=(median(etahat[college==1])/(1-median(etahat[college==1])))/(median(etahat[college==0])/(1-median(etahat[college==0])))         
# 50.94695        

# seems bad: what about drawing inferences from this pscore?

Ts=sum(rank(yppnscal)[college==1])
Et=sum(rank(yppnscal)*etahat)
Vt=sum((rank(yppnscal)^2)*etahat*(1-etahat))

1-pnorm(abs((Ts-Et)/sqrt(Vt)))
            
#0.01212666

# Note is on the scale of the IWP estimate above 


# unobserved confounding 

dev=abs((Ts-Et)/sqrt(Vt))

Gam=1
ps=Gam/(1+Gam)
et=sum(rank(yppnscal)*ps)
vt=sum((rank(yppnscal)^2)*ps*(1-ps))  

delt=((Ts-et)/sqrt(vt))

while(delt>dev){

	Gam=Gam+.001
	ps=Gam/(1+Gam)
	et=sum(rank(yppnscal)*ps)
	vt=sum((rank(yppnscal)^2)*ps*(1-ps))  

	delt=((Ts-et)/sqrt(vt))	
	
}
  
Gam
#[1] 2.404
  
# Thus the existing confounding due to X (as assumed to be in etahat) is equivalent to 
#  there being an unobserved factor Gamma of no more than 2.404 in magnitude 
#   - Thus, in this kind of study, IF an unobserved test is sensitivity to any Gamma
#      2.404 < then we should conclude that this is plausible given existing measures of 
#      confounding on X
  

# Perform a sensitivity test ignoring X; and assuming only U is remaining (also, note we are not matching on anything)   
        
source(file=url("http://sekhon.berkeley.edu/seminar/R/funs.R"))   
source(file=url("http://sekhon.berkeley.edu/seminar/R/ampRanksum.R"))   
                                                                        
ampRanksum(treated=yppnscal[college==1],control=yppnscal[college==0],Delta=10e100)  
                
# Gamma = 5; pvalue = .05

# suggests that even before matching on Kam and Palmer's propensity score the positive correlation is 
#   robust to plausible amounts of unoberved confounding; AND is robust to observed confounding 
#   captured in 81 X covariates 

# END


