###############################################
# Section 7: Topic Models
#
# John Henderson         
#
# PS 236B Spring 2013
###############################################

rm(list=ls())
                                       
library(tm) 
library(RWeka)  
library(stringr)             

library(lda) 
library(topicsmodels)       
  
# can look at bigrams of any length (within computational bounds)
# let's look at phrases of length 3 or less
BigramTokenizer = function(x,mins=1,maxs=3){
	NGramTokenizer(x, Weka_control(min = mins, max = maxs))
}
       
# 110th bills cosponsored by MCs
load('~/Dropbox/seminar/data/section6_a.Rdata')


# 0. some basic functions in R for manipulating strings        
# in base:
          
#(a) grep; searches accross objects for string matches 
poors=grep(pattern='poor',texts,ignore.case=TRUE)
poors=grep(pattern='poor',texts,ignore.case=TRUE,invert=T) # provides matches for those WITHOUT a match 
       
#(b) agrep; gives fuzzy matches based on Levenshtein distances
injury1=agrep(pattern='injury',texts,ignore.case=TRUE) 
injury2=agrep(pattern='injurya',texts,ignore.case=TRUE,max.distance=2) 

injury2=agrep(pattern=c('injurya','poor'),texts,ignore.case=TRUE,max.distance=2) 
                   
cat(c('smith','al'),sep='',file='cathere')
ad=read.table('cathere')

in2=in1=array(FALSE,length(texts))
in1[injury1]=TRUE
in2[injury2]=TRUE
         
# some non-overlap
table(in1,in2)             
     
#(c) gsub; finds and replaces patterns in strings
texts1=gsub(pattern='poor',replacement='wealthy',texts,ignore.case=TRUE)
poors=grep(pattern='poor',texts1,ignore.case=TRUE)  
                                                 
# note these mirror perl, bash, python scripting in R
#  - perl regular expressions can be used, with perl=TRUE
#  - sed, awk, tr; print and replace seamlessly, BUT inefficiently    

#(d) partial string matching:                                   
#  - agrep matches using a distance measure between two strings:
#     i.e., how many one-at-a-time character changes until the strings are identical
#  - partial matching uses a subset of characters to match exactly 
      
# char.expand; searches for unique matches using first n strings from left to right; gives string 
             
char.expand('public health servic',c(texts[1]))
char.expand('health servic',c(texts[1]))
char.expand('public health servic',c(texts[3]))  
char.expand('public health servic',c(texts[c(1,3)]))  

# limited to full string matching 
char.expand('plain york charl',c(texts[1]))              
          
# pmatch; searches for unique matches using first n strings from left to right; gives index  
pmatch('public health servic',c(texts[1]))
pmatch('public health servic',c(texts[3]))  
pmatch('public health servic',c(texts[c(1,3)]))

# charmatch; left to right truncated matches 
charmatch('pub',c(texts[1]))
charmatch('pub',c(texts[3]))  
charmatch('pub',c(texts[c(1,3)]))
             
# stringr; great package for string manipulation
??stringr
          
# str_sub; subsets strings on the character index
txt=str_sub(texts[1:100],1,5)
str_length(txt)                                                
str_length(texts)        

# (slow) way to count words [count the white spaces, one space at a time]
words=1   
m=str_length(texts[1])
j=0
while(j < m){
	j=j+1
	g=str_sub(texts[1],j,j) 
	if(g==' '){
		words=words+1
	}	
}       
   
# str_pad; adds white space at start and/or end of a string
texts[1]=str_pad(texts[1],width=5,side='both',pad=' ')                    
# str_trim; trim white space at start and/or end of a string  
texts[1]=str_trim(texts[1],side='both')                       
  
# offers many other manipulations; worth cehcking out!

# 1. managing texts in R with tm
??tm
# usefule guide to text manipulation:        
# http://cran.r-project.org/web/packages/tm/vignettes/tm.pdf

#data=read.csv('')                     

# look here to learn how to import data;
#  - can have this is a series of files, or one file index appropriately 
#  - recommendation is to work with R objects (a vector of strings)
  
# need to create an object of class corpus

vs=VectorSource(texts)  # creates a source for corpus to access
corpus=Corpus(vs)   	# creates the corpus
                                       
corpus     
         
# cleaning up the corpus ... 

# see http://cran.r-project.org/web/packages/tm/vignettes/tm.pdf for:
#   -  eliminating whitespace
#   -  removing stopwords 
#   -  converting to lower case
#   -  stemming           
# i.e., ...
#outs=Corpus(VectorSource(texts)) 
#outs=tm_map(outs, tolower)   
#outs=tm_map(outs, stripWhitespace)
#outs=tm_map(outs, stemDocument) 


# can inspect what is now in the corpus
text_data=inspect(corpus)       
# careful since inspect will print the entire full corpus 

# basic operators work on term-document matrices

tdm <- DocumentTermMatrix(corpus, control = list(tokenize = BigramTokenizer))

tdm # prints off some useful information
names(tdm) 
# ncol, nrow, dimnames   

(tdm$dimnames)[[2]][1:10] # give the word names 


# often want to work directly on the matrix object [carefully about printing too mcuh]
datas=inspect(tdm)
       
# legislator word counts
rowSums(datas)          
# word frequencies 
colSums(datas)        

# can compute distances matricex [be careful on dimensionality]                     
dmat=as.matrix(dist(datas[,1:100],upper=TRUE))                                                             
 
# find frequently occuring words in the corpus
findFreqTerms(tdm, 5)  # appeared at least 5 times
findFreqTerms(tdm, 1000)  # appeared at least 1000 times       
      
# finds words that co-occur
# findAssocs(tdm[1:20,], "public", 0.9) # fails since the document is too big

# 2. topic models using lda
# lda has its own efficient data storage functions 

# lexicalize for just 1-gram corpus analysis 
docs=lexicalize(texts)    

# topic model each line/ad; then combine these .... 
#set.seed(1005)
#num.topics <- 10
   
vocab=docs[[2]]        
docs=docs[[1]]


docs=lexicalize(texts,vocab=vocab)   

result <- lda.collapsed.gibbs.sampler(docs,
50,  ## Num clusters
vocab,
100,  ## Num iterations
0.1,
0.1,
compute.log.likelihood=TRUE)
  
# if we want n-grams use: 
library(topicmodels)

ldas=dtm2ldaformat(tdm,omit_empty = FALSE)  
vocab=ldas$vocab
docs=ldas$documents


result <- lda.collapsed.gibbs.sampler(docs,
50,  ## Num clusters
vocab,
100,  ## Num iterations
0.1,
0.1,
compute.log.likelihood=TRUE)        

??lda # for more info on the topic models

# END