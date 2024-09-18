###### OBTAINING ELEVATION VALUES USING COORDINATES ######
##########################################################

#remotes::install_github("r-spatial/mapview")
#install.packages("elevatr")
#install.packages("rgdal")

library(dplyr)
library(tidyverse)
library(sf)
library(elevatr)
library(mapview)
library(tidyr)

#set working directory
setwd()

#read dataset containing columns "Latitude" and "Longitude" 
d <- read.csv("samples.csv") 

#inspect data
head(d)

summary(as.numeric(d$Latitude))
summary(as.numeric(d$Longitude))

ggplot(d, aes(Longitude, Latitude)) +
  geom_count() 

#remove missing values and set as numeric
d <- d %>% drop_na(Latitude, Longitude)

d$Latitude <- as.numeric(d$Latitude)
d$Longitude <- as.numeric(d$Longitude)

#plot coordinates on a map, here using variable "type" to colour our points. 
m <- mapview(d, xcol = "Longitude", ycol = "Latitude",
             zcol = "type", crs = 4326, grid = FALSE, alpha = 0)
m

#other customisable variables
#map.types = c("Esri.WorldShadedRelief", "CartoDB.Positron")


### Get elevations ###
coord <- sf::st_as_sf(data.frame(d$Latitude, d$Longitude),
                      coords = c("d.Longitude", "d.Latitude"), crs = 4326)

#get_elev_point gets elevation here using the EPSG:4326 projection and zoom is set to 10.
#The default zoom is z=5 but this may not be high enough resolution
#can take some time...
coord <- get_elev_point(coord, prj="EPSG:4326", src = "aws", z = 10, overwrite = TRUE)


#join elevation data with our dataframe
d1 <- as.data.frame(cbind(d, coord$elevation)) %>% 
  rename("Elevation" = "coord$elevation")

#inspect
head(d1) 
hist(d1$Elevation, main = "Elevation across samples", xlab = "Elevation (m)")

#add elevations back to original file
d2 <- d1[c("Id","Elevation")]
d3 <- merge(d, d2, by = "Id", all.x = TRUE)
head(d3)

#write.csv(d3, "samples_with_elev.csv")
