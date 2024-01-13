


library(sf)
library(rgdal)
library(maps)
library(dplyr)

regions <- readOGR("National_GACC_Current_20210112.shp")


coords <- data.frame(Longitude=loc.data[,1], Lattitude=loc.data[,2])
coordinates(coords) <- c("Longitude","Lattitude")
as(coords,"SpatialPoints")
proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")


proj4string(coords) <- proj4string(regions)


pointsinpoly <- over(coords,regions)

dim(loc.data)

table(pointsinpoly$OBJECTID)
plot(loc.data)
for (j in 1:10){
  points(loc.data[which(pointsinpoly$OBJECTID==j),], col=j,pch=19)
}

unique(pointsinpoly$GACCName)

library(MASS)


# START  HERE


#attach department to observations

loc.data_2 <- cbind(loc.data,pointsinpoly$OBJECTID)
loc.data_2[1:2,3] <- 2

which(is.na(loc.data_2[,3]))


for (j in 1:20){
  loc.data_2[which(is.na(loc.data_2[,3]) & loc.data_2[,1] < -100 ),3] <- loc.data_2[which(is.na(loc.data_2[,3]) & loc.data_2[,1] < -100 )+1,3]
  which(is.na(loc.data_2[,3]))
}


for (j in 1:20){
  loc.data_2[which(is.na(loc.data_2[,3]) & loc.data_2[,1] > -100 ),3] <- loc.data_2[which(is.na(loc.data_2[,3]) & loc.data_2[,1] > -100 )-1,3]
  which(is.na(loc.data_2[,3]))
}

plot(loc.data)
for (j in 1:10){
  points(loc.data_2[which(loc.data_2[,3]==j),1:2], col=j,pch=19)
}


plot(loc.data)
for (j in c(4,9)){
  points(loc.data_2[which(loc.data_2[,3]==j),1:2], col=j,pch=19)
}




# read data:

#feature engineering



#
loc.data_2 <- cbind(loc.data_2, 1:dim(loc.data_2)[1])
colnames(loc.data_2) <- c('lon', 'lat', 'dep', 'pix')

loc.data_2 <- data.frame(loc.data_2)

data_DF <- left_join(data_DF, loc.data_2, by=c('lon'='lon', 'lat'='lat') )


mean(data_train_DF[which(data_train_DF$pix==1 & data_train_DF$year==1994),]$CNT, na.rm=TRUE)

names(data_DF)
data_train_DF <- data_DF[which(data_DF$dep %in% c(4,9)), -ncol(data_DF)]

head(data_train_DF)
Cali_DF <- data_train_DF

save(data_train_DF, file='Cali.RData')
#plot(data_train_DF$lon, data_train_DF$lat)
