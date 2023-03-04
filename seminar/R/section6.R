###############################################
# Section 6: Corpus Data and Scaling Documents 
#
# John Henderson         
#
# PS 236B Spring 2013
###############################################

rm(list=ls())
                                       
library(tm) 
library(RWeka)          
  
# can look at bigrams of any length (within computational bounds)
# let's look at phrases of length 3 or less
BigramTokenizer = function(x,mins=1,maxs=3){
	NGramTokenizer(x, Weka_control(min = mins, max = maxs))
}
  
load('~/Dropbox/seminar/data/section6_a.Rdata')
   
vs=VectorSource(texts)
corpus=Corpus(vs)

tdm <- DocumentTermMatrix(corpus, control = list(tokenize = BigramTokenizer))

# work on matrix object 
datas=inspect(tdm)
   
# 1. use naive classifier first 

wordscores=function(lefts,rights,datas){
	words=matrix(NA,ncol(datas),2)
	Wl=sum(lefts); Wr=sum(rights)
	words[,1]=lefts/Wl
	words[,2]=rights/Wr   
	
	# artificially eliminate degeneracies
	words[words[,1]==0,1]=min(words[words[,1]>0,1])
	words[words[,2]==0,2]=min(words[words[,2]>0,2])
	
	Pl=Wl/(Wl+Wr)  
	Plw=words[,1]/(words[,1]+words[,2])
	Prw=words[,2]/(words[,1]+words[,2])
	Ss=Prw-Plw                    
	lnwords=log(words[,1]/words[,2])
	alphas=matrix(NA,nrow(datas),2)   
	for(i in 1:nrow(datas)){
		Vw=sum(datas[i,])
		alphas[i,1]=sum(datas[i,]*Ss)/Vw
		alphas[i,2]=log(sum(words[,2]*datas[i,])/sum(words[,1]*datas[i,])) 
	}
	return(list("wordscores"=(alphas[,1]-mean(alphas[,1]))/sd(alphas[,1]),"bayesscores"=(alphas[,2]-mean(alphas[,2]))/sd(alphas[,2])))	
}        
    

#lefts=colSums(datas[which(a2<0),])
#rights=colSums(datas[which(a2>0),])
                                  
lefts=colSums(datas[which(party<0),])
rights=colSums(datas[which(party>0),])

scores=wordscores(lefts,rights,datas)
                                 

# 2. naive regression classifier   

# rather than motivate a probability model, let's maximize bivariate party prediction 
# -- k regressions, one for each word; then predict a score from a weighted average of word*coefficients              

coefsvars=coefs=array(NA,ncol(datas))

for(j in 1:ncol(datas)){
	smry=summary(lm(party~datas[,j]))
	coefs[j]=smry$coef[2,1]
	coefsvars[j]=smry$coef[2,2]
}      

wgts=(1-coefsvars) 
alpha_r1=alpha_r2=array(NA,nrow(datas))
for(j in 1:nrow(datas)){
	alpha_r1[j]=weighted.mean((datas[j,]*coefs)[which(datas[j,]>0)],w=rep(1,length(which(datas[j,]>0))))
	alpha_r2[j]=weighted.mean((datas[j,]*coefs)[which(datas[j,]>0)],w=coefsvars[which(datas[j,]>0)])	
}

# normalize at {0,1}
alpha_r1=(alpha_r1-mean(alpha_r1))/sd(alpha_r1)    
alpha_r2=(alpha_r2-mean(alpha_r2))/sd(alpha_r2)
scores[[1]]=(scores[[1]]-mean(scores[[1]]))/sd(scores[[1]])
scores[[2]]=(scores[[2]]-mean(scores[[2]]))/sd(scores[[2]])  
        
# highly correlated scales; {a,b} terms
alpha_outs=matrix(NA,length(alpha_r1),4)
alpha_outs[,1]=alpha_r1
alpha_outs[,2]=alpha_r2
alpha_outs[,3]=scores[[1]]
alpha_outs[,4]=scores[[2]]

cor(alpha_outs)
plot(c(cor(alpha_outs)),ylim=c(0,1))

# observe the linearity here     
plot.new()
par(mfrow=c(2,2))

plot(alpha_r2,alpha_r1,xlim=c(-2,2),ylim=c(-2,2)) 
abline(a=0,b=1,col='red',lty=2,lw=2)
plot(scores[[1]],alpha_r1,xlim=c(-2,2),ylim=c(-2,2))
abline(a=0,b=1,col='red',lty=2,lw=2)
plot(scores[[2]],alpha_r1,xlim=c(-2,2),ylim=c(-2,2))     
abline(a=0,b=1,col='red',lty=2,lw=2)
plot(scores[[2]],scores[[1]],xlim=c(-2,2),ylim=c(-2,2))     
abline(a=0,b=1,col='red',lty=2,lw=2)                                     
 



# 3. repeat with 1-grams to do wordfish 

load('~/Dropbox/seminar/data/section6_b.Rdata')                                       
datas=inspect(dtm) 
     
lefts=colSums(datas[which(party<0),])
rights=colSums(datas[which(party>0),])

scores=wordscores(lefts,rights,datas)

coefsvars=coefs=array(NA,ncol(datas))

for(j in 1:ncol(datas)){
	smry=summary(lm(party~datas[,j]))
	coefs[j]=smry$coef[2,1]
	coefsvars[j]=smry$coef[2,2]
}      

wgts=(1-coefsvars) 
alpha_r1=alpha_r2=array(NA,nrow(datas))
for(j in 1:nrow(datas)){
	alpha_r1[j]=weighted.mean((datas[j,]*coefs)[which(datas[j,]>0)],w=rep(1,length(which(datas[j,]>0))))
	alpha_r2[j]=weighted.mean((datas[j,]*coefs)[which(datas[j,]>0)],w=coefsvars[which(datas[j,]>0)])	
}

# normalize at {0,1}
alpha_r1=(alpha_r1-mean(alpha_r1))/sd(alpha_r1)    
alpha_r2=(alpha_r2-mean(alpha_r2))/sd(alpha_r2)
scores[[1]]=(scores[[1]]-mean(scores[[1]]))/sd(scores[[1]])
scores[[2]]=(scores[[2]]-mean(scores[[2]]))/sd(scores[[2]])    

      
# wordfish poisson model
source('~/Dropbox/seminar/R/Wordfish.1.3.R') 

wordouts=wordfish(input=t(datas),wordsincol=F,dir=c(1,2),tol=1e-4)         
alpha=wordouts$documents[,1]          
          
# again highly correlated
alpha_outs=matrix(NA,length(alpha_r1),5)
alpha_outs[,1]=alpha
alpha_outs[,2]=alpha_r1
alpha_outs[,3]=alpha_r2
alpha_outs[,4]=scores[[1]]
alpha_outs[,5]=scores[[2]]

cor(alpha_outs)
plot(c(cor(alpha_outs)),ylim=c(0,1))


plot.new()
par(mfrow=c(2,2))

plot(alpha_r1,alpha,xlim=c(-2,2),ylim=c(-2,2)) 
abline(a=0,b=1,col='red',lty=2,lw=2)
plot(alpha_r2,alpha,xlim=c(-2,2),ylim=c(-2,2)) 
abline(a=0,b=1,col='red',lty=2,lw=2)
plot(scores[[1]],alpha,xlim=c(-2,2),ylim=c(-2,2))
abline(a=0,b=1,col='red',lty=2,lw=2)
plot(scores[[2]],alpha,xlim=c(-2,2),ylim=c(-2,2))     
abline(a=0,b=1,col='red',lty=2,lw=2)
    

plot.new()
par(mfrow=c(2,2))

plot(alpha_r2,alpha_r1,xlim=c(-2,2),ylim=c(-2,2)) 
abline(a=0,b=1,col='red',lty=2,lw=2)
plot(scores[[1]],alpha_r1,xlim=c(-2,2),ylim=c(-2,2))
abline(a=0,b=1,col='red',lty=2,lw=2)
plot(scores[[2]],alpha_r1,xlim=c(-2,2),ylim=c(-2,2))     
abline(a=0,b=1,col='red',lty=2,lw=2)
plot(scores[[2]],scores[[1]],xlim=c(-2,2),ylim=c(-2,2))     
abline(a=0,b=1,col='red',lty=2,lw=2)

                 
# 4. how do these compare with scales from IRT models on binary choices?

# fix ordering first from different data
party_old=party

source('~/Dropbox/seminar/R/cosponsorEM.R')
      
load('~/Dropbox/seminar/R/cosponsor.Rdata')
load('~/Dropbox/seminar/data/ids.Rdata')
load('~/Dropbox/seminar/data/h_ids.Rdata')

cospo=as.matrix(cospo)[order(ids),]
#irtchoice=cosponsorEM(cospo,party)
load('~/Dropbox/seminar/data/irtchoice.Rdata')
irtalpha=irtchoice$alpha[order(ids)]   
irtalpha=(irtalpha-mean(irtalpha))/sd(irtalpha)

party=party[order(ids)] 

party_old=party_old[order(h_ids)]                 
alpha_outs=alpha_outs[order(h_ids),]
                      
alpha_outs=cbind(irtalpha,alpha_outs)
cor(alpha_outs)
            

plot.new()
par(mfrow=c(2,2))
plot(alpha_outs[,1],alpha_outs[,2],xlim=c(-2,2),ylim=c(-2,2))
abline(a=0,b=1,col='red',lty=2,lw=2)                                          
plot(alpha_outs[,1],alpha_outs[,3],xlim=c(-2,2),ylim=c(-2,2))
abline(a=0,b=1,col='red',lty=2,lw=2)
plot(alpha_outs[,1],alpha_outs[,5],xlim=c(-2,2),ylim=c(-2,2))
abline(a=0,b=1,col='red',lty=2,lw=2)


plot(density(alpha_outs[,1]),xlim=c(-2.5,2.5),ylim=c(0,.6))
lines(density(alpha_outs[,2]),xlim=c(-2.5,2.5),lty=2)
              
# plot the rank ordering
plot(rank(alpha_outs[,1]),rank(alpha_outs[,2]))


# END