library(XML)
#Get a simple table from Wikipedia
url <- "http://en.wikipedia.org/wiki/List_of_deadly_earthquakes_since_1900"
earthquakes <- readHTMLTable(url)
earthquakes <- earthquakes[[1]]
names(earthquakes) <- c("date", "country", "lat", "long", "depth", "magnitude", "secondary.effects", "pde.shaking.deaths", "pde.total.deaths", "utsu.total.deaths", "emdat.total.deaths", "other.source.deaths")
earthquakes$lat <- as.numeric(as.character(earthquakes$lat))
earthquakes$long <- as.numeric(as.character(earthquakes$long))
earthquakes$magnitude <-as.numeric(gsub(".*([0-9]{1}.[0-9]{1}).*",earthquakes$magnitude, replacement="\\1"))
earthquakes$utsu.total.deaths <- as.numeric(as.character(earthquakes$utsu.total.deaths))


library(maps)
world <- data.frame(map('world', plot=FALSE)[c("x","y")])
plot(world, type="l")
points(earthquakes$long, earthquakes$lat, cex=.7, pch=16, col="red")

plot(earthquakes$magnitude, log(earthquakes$utsu.total.deaths), pch=16, cex=.7)


#What if we have a badly formed table? We can use readLines
for (j in 1:200){
  url <- paste("http://ucpay.globl.org/index.php?campus=berkeley&name=&title=&base=&overtime=&extra=&gross=&year=2009&s=gross&p=",j,sep="")
#Illustrating use of ReadLines
  html <- readLines(url)
  entry.start <- grep("orow|erow", html)
  uc.pay.temp <- data.frame(matrix(nrow=50, ncol =4))
  names(uc.pay.temp) <- c("last.name", "first.name", "title", "pay")
  for (i in 1:length(entry.start)){
    entry <- html[entry.start[i]:(entry.start[i]+11)]
    uc.pay.temp$last.name[i] <-  gsub(".*>(.*) ,.*", entry[grep("&name",entry)], replacement="\\1")
    uc.pay.temp$first.name[i] <- gsub(".*, (.*)<.*", entry[grep("&name",entry)], replacement ="\\1")
    uc.pay.temp$title[i] <- gsub(".*Browse title\">(.*)</a>.*", entry[grep("title=!",entry)], replacement ="\\1")
    pay <- gsub(".*em\">(.*)</td.*", entry[grep("\"em",entry)], replacement="\\1")
#strip out the commas and dollar signs
    uc.pay.temp$pay[i] <- as.numeric(gsub("\\$|,", pay, replacement=""))
  }
  if(j==1){
    uc.pay <- uc.pay.temp
  }
  else{
    uc.pay <- rbind(uc.pay, uc.pay.temp)
  }
  print(j)
}

hist(uc.pay$pay[uc.pay$pay<500000], breaks=30, col="blue")
pay.title <- sort(tapply(uc.pay$pay, uc.pay$title, mean), decreasing=TRUE )

pay.prof <- uc.pay[grep("^PROFESSOR", uc.pay$title),]

library(ggplot2)
ggplot(pay.prof, aes(factor(title), pay)) + geom_boxplot()

#API

library(rjson)

api <- "XXXXXXXXXXX" ###### <<<API key goes here
q <- "Egypt" # Query string, use + instead of space
records <- 500 # total number of records to return
# calculate parameter for offset
os <- 0:(records/10-1)


# read first set of data in
url <- paste ("http://api.nytimes.com/svc/search/v1/article?format=json&query=", q, "&offset=", os[1], "&fields=date&api-key=", api, sep="")
raw.data <- readLines(url, warn="F") # get them
res  <- fromJSON(raw.data) # tokenize
dat <- unlist(res$results) # convert the dates to a vector

# read in the rest via loop
for (i in 2:length(os)) {
# concatenate URL for each offset
url <- paste ("http://api.nytimes.com/svc/search/v1/article?format=json&query=", q, "&offset=", os[i], "&fields=date&api-key=", api, sep="")
raw.data <- readLines(url, warn="F")
res  <- fromJSON(raw.data)
dat <- append(dat, unlist(res$results))  # append
print(i)
}

# aggregate counts for dates and coerce into a data frame
cts <- as.data.frame(table(dat))
 
# establish date range
dat.conv <- strptime(dat, format="%Y%m%d") # need to convert dat into POSIX format for this
daterange <- c(min(dat.conv), max(dat.conv))
dat.all <- seq(daterange[1], daterange[2], by="day") # all possible days
 
# compare dates from counts dataframe with the whole data range
# assign 0 where there is no count, otherwise take count
# (take out PSD at the end to make it comparable)
dat.all <- strptime(dat.all, format="%Y-%m-%d")
freqs <- ifelse(as.character(dat.all) %in% as.character(strptime(cts$dat, format="%Y%m%d")), cts$Freq, 0)
 
plot (freqs, type="l", xaxt="n", main=paste("Search term(s):",q), ylab="# of articles", xlab="date")
axis(1, 1:length(freqs), dat.all)
lines(lowess(freqs, f=.05), col = 2, lwd=2)


#XML Example
#Use the Congress API to get member biography data
url <- "XXXXXXXXXXX" ##need to supply your own api key

trent.lott <- xmlTreeParse(url)
 trent.lott  <- xmlToList(trent.lott)
