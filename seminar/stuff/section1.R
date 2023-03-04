library(foreign)
library(ggplot2)
setwd("~/Dropbox/Projects/PS236b/section1/")
#Read in the data

county.codes <- c("001", "003", "005", "007", "009", "011", "013", "015", "017", "019", "021", "023", "025", "027", "029", "031", "033", "035", "037", "039", "041", "043", "045", "047", "049", "051", "053", "055", "057", "059", "061", "063", "065", "067", "069", "071", "073", "075", "077", "079", "081", "083", "085", "087", "089", "091", "093", "095", "097", "099", "101", "103", "105", "107", "109", "111", "113", "115")

###
#DOWNLOAD DATA
###

#Create a list to contain all the data
county.results <- list()
for (i in 1:length(county.codes)){
                                        #form url
  url <- paste("http://swdb.berkeley.edu/pub/data/P10/c",county.codes[i],"/c", county.codes[i], "_p10_sov_data_by_p10_svprec.dbf", sep ="")
                                        #download file
  download.file(url, paste("c", county.codes[i], "_p10_sov_data_by_p10_svprec.dbf", sep=""))
  county.results[[i]] <- read.dbf(paste("c", county.codes[i], "_p10_sov_data_by_p10_svprec.dbf", sep=""))
  names(county.results)[i] <- paste("c",county.codes[i], sep="")
  if(nrow(county.results[[i]])>0){
    county.results[[i]]$county <- paste("c",county.codes[i], sep="")
    county.results[[i]] <- county.results[[i]][,c(ncol(county.results[[i]]),1:(ncol(county.results[[i]])-1))]
  }
}


####
#CLEAN DATA
####

#Check data
n.precincts <- sapply(county.results,nrow)
#remove counties with bad data
bad.precincts <- names(n.precincts)[n.precincts==0|n.precincts==1]
county.results <- county.results[names(county.results)%in%bad.precincts==FALSE]



#First focus on county
#
#Aggregate by precinct
c001data <- county.results[[1]]
#Need to remove the "A" from SVPREC
c001data$SVPREC <- gsub("A", x = c001data$SVPREC, replacement="")
c001data <- melt(c001data, id.vars = 1:5)
c001data <- cast(c001data, county + SVPREC + ADDIST + CDDIST + SDDIST ~ ..., sum)
#Also some lines are totals, get rid of any precincts that begin with text
c001data <- c001data[ grep("^([A-Z]).*",c001data$SVPREC) * -1, ]

#Now do all counties
#Note the use of regular expressions
AggCounties <- function(data){
  print(unique(data$county))
  data$SVPREC <- gsub("A", x = data$SVPREC, replacement="")
  data <- melt(data, id.vars = 1:5)
  data <- cast(data, county + SVPREC + ADDIST + CDDIST + SDDIST ~ ..., sum)
  data <- data[grep(".*([A-Z])$",data$SVPREC)*-1,]
  data <- data[grep("^SOV.*",data$SVPREC)*-1,]
}
county.results  <- llply(county.results[1:10], AggCounties, .progress="text")
load("~/Desktop/county.results.RData")

#Put together one dataset
ExtCounties <- function(data){
  output <- data[,c("county", "SVPREC", "ADDIST", "CDDIST", "SDDIST","TOTREG","TOTVOTE")]
  gov <-  data[,grep("^GOV", names(data))]
  pr <-  data[,grep("^PR", names(data))]
  output <- cbind(output,gov,pr)  
}
statewide.results <- ldply(county.results, ExtCounties, .progress="text")


#Merge in covariates
county.covar <- list()
for (i in 1:length(county.codes)){
                                        #form url
  url <- paste("http://swdb.berkeley.edu/pub/data/P10/c",county.codes[i],"/c", county.codes[i], "_p10_voters_by_p10_srprec.dbf", sep ="")
                                        #download file
  download.file(url, paste("c", county.codes[i], "_p10_voters_by_p10_srprec.dbf", sep=""))
  county.covar[[i]] <- read.dbf(paste("c", county.codes[i], "_p10_voters_by_p10_srprec.dbf", sep=""))
  names(county.covar)[i] <- paste("c",county.codes[i], sep="")
  if(nrow(county.covar[[i]])>0){
    county.covar[[i]]$county <- paste("c",county.codes[i], sep="")
    county.covar[[i]] <- county.covar[[i]][,c(ncol(county.covar[[i]]),1:(ncol(county.covar[[i]])-1))]
  }
}

ExtCovar <- function(data){
  output <- data[,c(1, 4,5,6,7,16,17,18,19)]
}
statewide.covar <- ldply(county.covar, ExtCovar, .progress="text")



statewide.merged  <- merge(statewide.results, statewide.covar, by.x=c("county","SVPREC"), by.y=c("county", "SRPREC"), all.x=TRUE, all.y=FALSE)
#Check for problems


dim(statewide.merged[is.na(statewide.merged$HISPREP),1:10])


###
#Analysis
###
statewide.merged$turnout <- statewide.merged$TOTVOTE/statewide.merged$TOTREG
statewide.merged <- statewide.merged[statewide.merged$turnout<1,]
statewide.merged$treatment <- ifelse(statewide.merged$TOTREG>250,"control", "treatment")

gg.discont <-  ggplot(statewide.merged[statewide.merged$TOTREG<1000,], aes(x=TOTREG, y = turnout))
gg.discont + geom_point()
gg.discont + geom_point(aes(color=treatment)) 
gg.discont + geom_point(aes(color=treatment), alpha = .2) + geom_vline(xintercept=250)+ stat_smooth(method=loess, aes(fill=treatment), se=TRUE, size=1.5, color = "black")



