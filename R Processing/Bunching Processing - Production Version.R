# Bunching Processing.R
# Author: Mark Egge
# Last Revised: 5/12/16

# broadly speaking, given an input file with times, latitudes, and longitudes
# the script below will calcualte if given points are within 250m
# of each other at a given time
# takes ~20 minutes to run (much imrpoved over plyr, which took hours)
# may not run as written, but contains all of the code which produced the data

library(ggplot2)
library(data.table)
library(geosphere)
library(sp)
library(rgeos)
library(rgdal)

# Start the clock!
ptm <- proc.time()

working.file <- "march-16-all-rts.csv"

# define a headway threshold to be used below. Buses less than 250 meters from the
# next bus on the same line will be coded as "bunched"
headway.threshold <- 250 

options(scipen=20) # display numbers up to 20 decimal points w/o scientific notation

# the import.csv funtion that lurks in the code of every R project written by
# Artur Dubrawski's Data Mining students. I've commented out the first line which
# only reads in the first 50k lines (for testing purposes)
import.csv <- function(filename){
  #return(read.csv(filename, header=FALSE, sep=",", strip.white=TRUE, nrows=50000, fileEncoding = "UTF-8-BOM"))
  return(read.csv(filename, header=FALSE, sep=",", strip.white=TRUE, fileEncoding = "UTF-8-BOM"))
}

this.dir <- dirname(parent.frame(2)$ofile); setwd(this.dir) # set working directory to script directory

if(!is.data.frame(df)) { df <- import.csv(working.file) } # import data if not already in working env

# add headers from another file...
headers <- read.csv(file="headers.csv", header=TRUE)
names(df) <- names(headers)
remove(headers) # and cleanup workspace


# do my date conversions and column derivation before going to a table. Drop my temporary datetime column 
# before going to data.table because data.table hates POSIXc
df$dt <- strptime(df$DateTime, format='%Y-%m-%d %H:%M:%S.000')
df$Route <- substr(as.character(df$RouteNumber), 0, 2) # convert from 71A to 71, etc.
df$YMD <- format(df$dt, format='%y%m%d') 
df$M <- as.numeric(format(df$dt, format='%M'))
df$HR <- as.numeric(format(df$dt, format='%H'))
df$MOD <- (df$HR * 60) + df$M # Minute of Day
df$TIMESEG <- floor(df$MOD / 5) # define 5 minute time segments
df$latitude <- df$CurrentLatitude # copy lat long columns for later use
df$longitude <- df$CurrentLongitude

df <- subset(df, select=-c(dt, M)) # get rid of columns previously used for processing

df <- df[df$HR >= 6, ] # keep only data between 6:00 am and midnight



coordinates(df) <- c("CurrentLongitude", "CurrentLatitude") # convert to spatial data frame class,
proj4string(df) <- CRS("+proj=longlat +datum=WGS84") # set projection to WGS84
# Pittsburgh: UTM Zone 17T
res <- spTransform(df, CRS("+proj=utm +zone=17T ellps=WGS84")) # reproject from WGS84 (decimal degrees) into UTM (meters)

# given a spatial dataframe containing latitudes and longitude pairs, return TRUE
# for those rows within 250 meters of another row, FALSE otherwise
is.bunched <- function(subset) {
  coordinates(subset) <- c("CurrentLongitude", "CurrentLatitude")
  (apply(gWithinDistance(subset, subset, headway.threshold, byid=TRUE), 1, sum) > 1)
}

dt <- as.data.table(res) # convert to data table
setkeyv(dt, c('YMD', 'MOD', 'Route', 'RouteDirection')) # sort by minute of date, route, and route direction
dt <- dt[, BUNCHED := is.bunched(.SD), by=.(YMD, MOD, Route, RouteDirection)]
dt <- dt[, BUNCHED10 := ifelse(BUNCHED, 1, 0)] # code TRUE/FALSE as 1/0
dt <- dt[, DISTSEG := floor(Distance / 1000)] # divide distance along path in 1000' increments
dt <- dt[, TIMESEG := floor(df$MOD / 10)] # change time segments to ten minute intervals, just to make the vis easier

# create a distance segment and calculate departure time
setkeyv(dt, c('YMD', 'RouteNumber', 'RouteDirection', 'TripUniqueID'))
dt <- dt[, DEPTTIMESEG := min(TIMESEG), by=.(YMD, RouteNumber, RouteDirection, TripUniqueID)]

# can group by departure time segment DEPTTIMESEG or observation time segment TIMESEG
# reorder to sort by route, direction, DTIMESEG
setkeyv(dt, c('RouteNumber', 'RouteDirection', 'DEPTTIMESEG', 'DISTSEG'))
depart.average <- dt[, .(bunchmean = round(mean(BUNCHED),3)), by=.(RouteNumber, RouteDirection, DEPTTIMESEG, DISTSEG)]

names(depart.average) <- c("route", "outbound", "time", "distance", "value")
depart.average$DepartTime <- 1
depart.average <- depart.average[, outbound := ifelse(outbound=="OUTBOUND", 1, 0)] # reduce file size

# write files for observed times
for(j in 0:1) { # inbound / outbound
  routes <- levels(depart.average$route)
  for(i in 1:length(routes) ) {
    filename <- paste("mean-depart", ifelse(j==0, "inbound", "outbound"), routes[i], "tsv.txt", sep="-")
    write.table(depart.average[route==routes[i] & outbound==j,.(time, distance, value)], filename, sep="\t", row.names=FALSE)  
  }
}

# group by departure time segment by observation time segment TIMESEG
# reorder to sort by route, direction, DTIMESEG
setkeyv(dt, c('RouteNumber', 'RouteDirection', 'TIMESEG', 'DISTSEG'))
obs.average <- dt[, .(bunchmean = round(mean(BUNCHED),3)), by=.(RouteNumber, RouteDirection, TIMESEG, DISTSEG)]

names(obs.average) <- c("route", "outbound", "time", "distance", "value")
obs.average <- obs.average[, outbound := ifelse(outbound=="OUTBOUND", 1, 0)] # reduce file size

# write files for observed times
for(j in 0:1) { # inbound / outbound
  routes <- levels(obs.average$route)
  for(i in 1:length(routes) ) {
    filename <- paste("mean-obs", ifelse(j==0, "inbound", "outbound"), routes[i], "tsv.txt", sep="-")
    write.table(obs.average[route==routes[i] & outbound==j,.(time, distance, value)], filename, sep="\t", row.names=FALSE)  
  }
}


# combine the two into a single file
depart.average$DepartTime <- 1
obs.average$DepartTime <- 0
comb.average <- rbindlist(list(depart.average, obs.average), use.names=TRUE, fill=FALSE)
write.table(comb.average[, .(time, distance, value, route, outbound, DepartTime)], "all-means-tsv.txt", sep="\t", row.names=FALSE)

# melt our data to change its shape from long to wide:
#recast.average <- dcast(bunched.average, RouteNumber + RouteDirection + DEPTTIMESEG ~ DISTSEG, value.var="bunchmean", fill=0)
#write.csv(recast.average, paste("output-bunchingmeans-", working.file, sep=""), row.names=FALSE)


# take the obs calculated as bunched and convert back to a data frame
df <- as.data.frame(dt[BUNCHED==TRUE, ])

# Stop the clock
cat("Elapsed time: ", proc.time() - ptm, "seconds elapsed\n")
cat("Rows:", nrow(df))


df$dt <- strptime(df$DateTime, format='%Y-%m-%d %H:%M:%S.000')
df$HHMM <- format(df$dt, format='%I:%M %p')
write.csv(df, paste("output-", working.file, sep=""))


df$DateTime <- as.POSIXct(df$DateTime)
df$CDATETIME <- ISOdatetime(2016,03,01, format(df$DateTime, '%H'), format(df$DateTime, '%M'), 0, tz="EST")
write.csv(df[, names(df) %in% c("CDATETIME", "latitude", "longitude")], paste("output-cartodb-", working.file, sep=""), row.names=FALSE)
  
# let's graph the relative frequency of bunching by route
route.means <- dt[,list(avg = mean(BUNCHED10)), by=.(RouteNumber, RouteDirection)]
route.means$Route <- substr(as.character(route.means$RouteNumber), 0, 2)
ggplot(data = route.means, aes(x=RouteNumber, y=avg, fill=Route)) + geom_bar(stat="identity") +
  labs(title="Bunched Portion of Observations", x="Route Number", y="% of Observations") +
  facet_wrap(~RouteDirection)

# to do--find the origin point of clustering for any given trip



# junkyard

#write.table(depart.average[, .(time, distance, value, route, Outbound, DepartTime)], "departed-time-means-tsv.txt", sep="\t", row.names=FALSE)
#write.table(depart.average[route=="61C" & outbound=="1",.(time, distance, value, route, outbound, DepartTime)], "bunched-means-61c-depart-time-tsv.txt", sep="\t", row.names=FALSE)
#write.table(depart.average[route=="61D" & outbound=="OUTBOUND",.(time, distance, value)], "bunched-means-61d-depart-time-tsv.txt", sep="\t", row.names=FALSE)

#write.table(obs.average[, .(time, distance, value, route, outbound)], "observed-time-means-tsv.txt", sep="\t", row.names=FALSE)
#write.table(obs.average[route=="61C" & outbound=="OUTBOUND",.(time, distance, value)], "bunched-means-61c-observed-time-tsv.txt", sep="\t", row.names=FALSE)
#write.table(obs.average[route=="61D" & outbound=="OUTBOUND",.(time, distance, value)], "bunched-means-61d-observed-time-tsv.txt", sep="\t", row.names=FALSE)


# fill in missing values for missing minutes with East Busway Swissvale
#times <- seq(ISOdate(2016,03,01,06,00,00,tz="EST"), ISOdate(2016,03,01,23,59,59,tz="EST"), by="min")
#str.times <- format(times, format='%I:%M %p')
#missing <- str.times[!(str.times %in% df$HHMM)]
#missing.rows <- data.frame("HHMM"=missing, "RouteNumber" = "00", "RouteDirection" = "INBOUND", "latitude" = 40.42936, "longitude"= -79.91794)
#carto.df <- rbind(df[, c("HHMM", "RouteNumber", "RouteDirection", "latitude", "longitude")], missing.rows)
#carto.df <- df
#t <- as.POSIXct(carto.df$HHMM, format='%I:%M %p', origin = "2016-03-01", tz = "EST")
#carto.df$DateTime <- format(t, '%Y-%m-%d %H:%M:00')
#carto.df <- carto.df[order(carto.df$DateTime, carto.df$RouteNumber, carto.df$RouteDirection), ]
#carto.df$M <- as.numeric(format(carto.df$DateTime, '%M'))
#carto.df$HR <- as.numeric(format(carto.df$DateTime, format='%H'))
#carto.df$MOD <- (carto.df$HR * 60) + carto.df$M
#carto.df <- subset(carto.df, select=-c(M, HR))


#dist.matrix <- gDistance(subset, subset, byid=TRUE)
#min.distances <- apply(dist.matrix, 1, function(x) sort(x, partial=2)[2])
#sort(dist.matrix[1, ], partial=2)[2] #return second least value
