## Establishes environment
library(tidyverse)
library(ggmap)
library(ocedata)
library(RgoogleMaps)
library(maps)
library(mapproj)
library(mapdata)
library(marmap)
library(units)
library(sf)
library(ggmap)
library(patchwork)
library(raster) #may need to reinstall
library(marmap)
library(tidyverse)

## Reads in data
getwd()
lumpz <- read.csv("7/Zero lumpfish1989.csv")
lumpnz <- read.csv("7/Nonzero lumpfish1989.csv")


#TRIES TO COMBINE DF TO COMPARE
lumpz["catch"] <- c("zero")
lumpnz["catch"] <- c("nonzero")
lumpz %>% left_join(lumpnz)


#Useful options in map_data include usa, state, county, 
#world (0? in middle), and world2 (180? in middle)
?map_data()
iceland <- map_data('world', "iceland")   
head(iceland)

#the group command explains which polygons belong together, one 
#polygon for each island on the map
ggplot(data=iceland, aes(x=long, y=lat)) + 
  geom_polygon()

#Calculates ratio for long and lat (better projection)
# Djupivogur lat = 64.6568° N, long = 14.2848° W
LatitudeDeg <- 64.657  #Seattle
1/cos(LatitudeDeg/180*pi) # = 2.336

p1 <- ggplot(data=iceland, aes(x=long, y=lat)) + 
  geom_polygon(alpha=0.5)+ 
  coord_fixed(ratio=2.34) +
  geom_point(data = lumpz, aes(x=lon, y=lat), 
             color = "red")+
  geom_point(data = lumpnz, aes(x=lon, y=lat))+
  theme_minimal()
p1

bathy <- getNOAA.bathy(lon1=-28, lon2=-10, lat1=62.5, lat2=68,
                       resolution=2, keep=FALSE) 
map <- autoplot.bathy(bathy, geom=c('raster'),
                      show.legend=FALSE) +     #turn off legend
  scale_fill_etopo() +                   #special topographic colors
  theme(axis.title = element_blank()) +  #remove the axis titles
  scale_x_continuous(breaks=seq(-130,-118, 2),  #where to place the values
                     labels=paste0(seq(130,118, -2),'W'), 
                     expand = c(0, 0)) + 
  scale_y_continuous(breaks=seq(40,50,2),  #where to place the values
                     labels=paste0(seq(40,50,2),'N'),
                     expand = c(0, 0))+
  geom_point(data = lumpz, aes(x=lon, y=lat), 
             color = "red", size = '.')+
  geom_point(data = lumpnz, aes(x=lon, y=lat))+
  theme_minimal()
map

