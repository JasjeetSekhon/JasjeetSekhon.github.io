#to install rgdal on mac
#setRepositories(ind=1:2)
#install.packages("rgdal")
library(sp)
library(maps)
library(maptools)
library(rgdal)
library(RColorBrewer)
library(classInt)
library(reshape)
library(ggplot2)


#Use the builtin maps
county <- map(database="county")
map(database="county", "california")
world <- map(database="world")


#geocode
library(dismo)
fdh <- geocode('180 Noe St. San Francisco CA 94114')
map(database="county", "california,san francisco")
points(x=fdh[2], y=fdh[3], pch=16)
map(database="county", "california")
points(x=fdh[2], y=fdh[3], pch=16)


###
###Fraud in Pernambuco
###
load("/Users/ElGuapo/Dropbox/projects/ps236b/spatial/pe.fraud.RData")

qplot(pe$pct.elec.change[pe$pct.elec.change<100])
ggplot(pe[pe$electorate.55<30000,], aes(x=electorate.55, y = electorate.58)) + geom_point() + geom_smooth()  + geom_abline(intercept=0, slope = 1) 

#define the coordinate system and the projection
llCRS <- CRS("+proj=longlat +ellps=WGS84")
#use the SpatialPointsDataFrame class to tell R that it's dealing with spatial data
pe.sp <- SpatialPointsDataFrame(coords = pe[,c("long","lat")], data = pe[,c(1,4:8)], proj4string = llCRS)
class(pe.sp)
#notice the slots
str(pe.sp, max.level=2)

plot(pe.sp)
plot(pe.sp, cex = (pe.sp$pct.elec.change + 66.89)/25, pch=19)


setwd("~/Dropbox/projects/Brazil Electoral Data/maps/amc/AMC1960-2000")
#Use ReadOGR to read in shape files (as well as other formats)
brazil.LL <- readOGR(dsn = ".", layer = "9760", input_field_name_encoding = "UTF-8")
proj4string(brazil.LL) <- CRS("+proj=longlat ellps=WGS84")
#notice the slots
str(brazil.LL, max.level=2)
brazil.LL@data

plot(brazil.LL)
#create state codes
brazil.LL$uf.code <- gsub("^([0-9]{2}) .*", x = brazil.LL$NEWCOD9760, replacement="\\1")

#subset to the state of Pernambuco
pe.LL <- brazil.LL[brazil.LL$uf.code==26,]
plot(pe.LL)
pe.LL <- pe.LL[-1* grep("Fernando de Noronha", as.character(pe.LL$MIN_NOMEMU)),]

plot(pe.LL)
points(pe$long, pe$lat, pch=16)

#Find which municipalities belong to which AMC (mimimal comparable areas)
#use the overlay command, returns a vector of indices, saying which polygon each point belongs to
overlay.pe <- overlay(pe.sp, pe.LL)
#assign the polygon
pe.sp$amc <- pe.LL@data$NEWCOD9760[overlay.pe]
#for polygons with multiple observations, average within each AMC
pe.amc <- melt(pe.sp@data[is.na(pe.sp$amc)==FALSE,c("electorate.58","electorate.55", "elec.change", "pct.elec.change", "amc")], id.vars = "amc")
pe.amc <- data.frame(cast(pe.amc, amc ~variable, mean, na.rm=TRUE))
#need to create a dataset with ALL AMCs to merge back into map object
all.amc <- data.frame(amc = as.character(pe.LL@data$NEWCOD9760))
pe.amc <- merge(all.amc, pe.amc, all.x=TRUE)

#SP package merges by rownames, so both datasets need identical rownames
row.names(pe.amc) <- pe.amc$amc
row.names(pe.LL) <- as.character(pe.LL@data$NEWCOD9760)


pe.LL.fraud <- spCbind(pe.LL, pe.amc)
str(pe.LL.fraud, max.level=2)
head(pe.LL.fraud@data)

plot(pe.LL.fraud)


ncol <- 8 #number of colors
#use the Rcolor library to choose good colors
plotclr <- brewer.pal(ncol, "Oranges")
#reverse order, since big drops in registration indicates MORE fraud
plotclr <- plotclr[ncol:1]
#Create equal sized classes
class <- classIntervals(round(pe.LL.fraud$pct.elec.change), ncol, style="quantile", dataPrecision =0)
colcode <- findColours(class, plotclr, digits = 2)
plot(pe.LL.fraud, density=16, col="grey", axes=T, cex.axis=.75)
plot(pe.LL.fraud, col =colcode, add = TRUE)
title(xlab="Longitude",ylab='Latitude',cex.lab=.75,line=2.25)
title(main="Electoral Fraud in Pernambuco (1958)", sub="% Change in Registered Voters",font.sub=2)
legend(-40, -9, legend=(names(attr(colcode, "table"))), fill=attr(colcode, "palette"), cex=0.75, bty="n")

library(RgoogleMaps)
#create a bounding box (borders of the map to request) using qbbox
earth.box <- qbbox(lat = pe$lat, lon = pe$long)
#get the map from google maps
pe.google <- GetMap.bbox(earth.box$lonR, earth.box$latR, maptype="satellite")
PlotOnStaticMap(pe.google, lat = pe$lat, lon = pe$lon, FUN = points, pch =16, col="red")


#Calculate distances
#coordinates on a polygon dataset will produce centroids
centroids <- coordinates(pe.LL.fraud)
plot(pe.LL.fraud)
points(centroids, pch = 16)
recife.geocode <- geocode('Recife, PE Brazil')


#convert to radians
deg2rad <- function(deg) return(deg*pi/180)

# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Haversine formula (hf)
gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

for(i in 1:nrow(centroids)){
  print(gcd.hf(deg2rad(recife.geocode[2]), deg2rad(recife.geocode[3]), deg2rad(centroids[i,1]), deg2rad(centroids[i,2])))
}

for(i in 1:nrow(pe)){
  pe$dist[i] <- gcd.hf(deg2rad(recife.geocode[2]), deg2rad(recife.geocode[3]), deg2rad(pe$long[i]), deg2rad(pe$lat[i]))
}

ggplot(pe[pe$pct.elec.change<100,], aes(x=dist, y=pct.elec.change)) + geom_point() + stat_smooth()
