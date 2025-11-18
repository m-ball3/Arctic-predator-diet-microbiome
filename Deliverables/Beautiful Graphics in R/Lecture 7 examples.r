#=============CREATING MAPS===========================================
#Fish 554, Beautiful Graphics in R.
#Different projections, creating a map of your study site
#By Trevor A. Branch tbranch@uw.edu, with reliance on lecture
#notes from Allan Hicks, Elizabeth Ng, Kristina Randrup. 
#Examples from: Eric Anderson, Dawn Barlow, Mike Johns, and more.
#VERSION note: most of the packages here are updated very
#  frequently, and use the latest versions of R, so it is 
#  a good idea to update R and set Rstudio to use the latest version
#  type "version" in the console to see your version (4.4.2 as of 18 Nov 2024)
#  and Tools->Global Options->General->R Version to update.
#=====================================================================
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

#===========ggplot with geom_polygon==================================
#These examples selected from Eric Anderson
#http://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html
#=====================================================================
library(tidyverse) #including ggplot2
library(ggmap)     #contains useful mapping data

#Useful options in map_data include usa, state, county, 
#world (0? in middle), and world2 (180? in middle)
usa <- map_data('usa')   
head(usa)

#the group command explains which polygons belong together, one 
#polygon for each island on the map
ggplot(data=usa, aes(x=long, y=lat, group=group)) + 
  geom_polygon()

#need to fix the relation between 1deg latitude and 1deg longitude
#at the equator this is 1:1 but nearer the poles it gets very large 
ggplot(data=usa, aes(x=long, y=lat, group=group)) + 
  geom_polygon() + 
  coord_fixed(ratio=1.25) #ratio of longitude to latitude

#1 deg latitude divided by cos(latitude in radians) gives this ratio
#The ratio varies by latitude and maps usually go across many latitudes
#therefore it is necessary to use **projections** to ensure that 
#maps are appropriately scaled. 
LatitudeDeg <- 47.606  #Seattle
1/cos(LatitudeDeg/180*pi)

LatitudeDeg <- 29.7604  #Houston differs quite a lot
1/cos(LatitudeDeg/180*pi)

#The states data
states <- map_data('state')
head(states)

#color by state
ggplot(data=states, aes(x=long, y=lat, group=group, fill=region)) + 
  geom_polygon() + theme_minimal() +
  coord_fixed(ratio=1.25) +  #ratio of longitude to latitude
  theme(legend.position='none') #remove the legend

#nice-ish map of the US showing states
#tidy up lat and long labels, change the transparency
#add some points
states <- map_data('state')
ggplot(data=states, aes(x=long, y=lat, group=group, fill=region)) + 
  coord_fixed(ratio=1.25) +  #ratio of longitude to latitude
  theme_minimal() + 
  geom_polygon(alpha=0.4, col='white') + #state boundaries
  theme(legend.position='none',
        axis.title = element_blank()) +  #remove the legend
  scale_x_continuous(breaks=seq(-70,-120,-10),  #where to place the values
                     labels=paste0(seq(70,120,10),'W')) + 
  scale_y_continuous(breaks=seq(25,50,5),  #where to place the values
                   labels=paste0(seq(25,50,5),'N')) + 
  geom_point(aes(x=-122.3321, y= 47.6062), size=2) + #add point
  geom_text(aes(x=-122.3321, y= 47.6062, label='Seattle'), 
            nudge_x=4.5, color='gray20', size=4)  #add text
  

#=============Chloropleth map with geom_polygon=======================
#Unemployment rates in Washington state during 2009 recession
#by county
#Full data set: http://datasets.flowingdata.com/unemployment09.csv
#=====================================================================
library(tidyverse) #including ggplot2
library(ggmap)     #contains useful mapping data
counties <- map_data("county")
WA_county <- subset(counties, region == "washington")
head(WA_county)

unemp <- read_csv('Data\\unemployment09WA.csv')
#change county names to lower case to match map_data names
unemp$subregion <- str_to_lower(unemp$subregion)
head(unemp)

#join the two datasets together
WA_county_unemp <- left_join(x=WA_county, y=unemp, by='subregion')
head(WA_county_unemp)

#use unemployment as the color
ggplot(data=WA_county_unemp, aes(x=long, y=lat, group=group, 
                                 fill=Unemployment)) + 
  geom_polygon(color='white') + #county outline in white
  theme_minimal() + 
  coord_fixed(ratio=1.48)  #ratio of longitude to latitude in Washington

#fix the axes, titles, change the color scales, etc. 
ggplot(data=WA_county_unemp, aes(x=long, y=lat, group=group, 
                                 fill=Unemployment)) + 
  geom_polygon(color='white') + #county outline in white
  theme_minimal() + 
  theme(legend.title = element_blank(),
        axis.title = element_blank()) + 
  coord_fixed(ratio=1.48) +  #ratio of longitude to latitude in Washington
  labs(title='Unemployment in Washington in 2009 (%)') + 
  scale_x_continuous(breaks=seq(-125,-117,2),  #where to place the values
                     labels=paste0(seq(125,117,-2),'W')) + 
  scale_y_continuous(breaks=seq(46,49,1),  #where to place the values
                     labels=paste0(seq(46,49,1),'N')) + 
  scale_fill_continuous(high = '#132B43', low = '#56B1F7')

#save the file, to ensure result is device-independent 
#I always save figures to a subdirectory "Figs"
ggsave('Figs\\WA unemployment class.png', dpi=600, width=9, height=7)

#=====Earthquakes using marmap package with bathymetry/topography=====
#Data and maps showing bathymetry and topography
#Depth/height from the 'ETOPO1' dataset. The marmap package is here:
#Pante E & Simon-Bouhet B (2013) PLoS ONE 8(9): e73051. 
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0073051
#Code idea: Dawn Barlow.
#Earthquake data downloaded 2/29/2020 from:
#https://earthquake.usgs.gov/earthquakes/search/
#=====================================================================

#installing marmap can take 30-60 seconds
#make sure R is up to date (v 4.0.x) with "version" in console
library(raster) #may need to reinstall
library(marmap)
library(tidyverse)

#get the bathymetry from NOAA. Resolution in minutes, default = 4, 
#minimum = 1. Smaller = higher resolution and longer times
#keep=TRUE writes downloaded data into a file
bathy <- getNOAA.bathy(lon1=-130, lon2=-118, lat1=40, lat2=50,
            resolution=2, keep=FALSE) 

#create a ggplot object appropriate to the bathy data object
map <- autoplot.bathy(bathy, geom=c('raster', 'contour'),
                show.legend=FALSE) +     #turn off legend
  scale_fill_etopo() +                   #special topographic colors
  theme(axis.title = element_blank()) +  #remove the axis titles
  scale_x_continuous(breaks=seq(-130,-118, 2),  #where to place the values
                     labels=paste0(seq(130,118, -2),'W'), 
                     expand = c(0, 0)) + 
  scale_y_continuous(breaks=seq(40,50,2),  #where to place the values
                     labels=paste0(seq(40,50,2),'N'),
                     expand = c(0, 0))
map

#add earthquakes to the map, all >2.5 magnitude in 2019 
#40-50N and 118-130W, downloaded from 
#https://earthquake.usgs.gov/earthquakes/search/
quakes <- read_csv(file='Data\\PNW earthquakes 2019.csv')
quakes
map + geom_point(data = quakes, 
          aes(x=longitude, y=latitude, size=mag),
          shape=21, fill='#00000044', color='black') + 
  labs(title='Earthquakes in 2019')
ggsave('Figs\\PNW earthquakes.png', width=8, height=10, dpi=600)

#=================PROJECTIONS ggplot==================================
#Plot of sea surface temperature around Antarctica 
#Data and code slightly modified from Mike Johns mjohns@pointblue.org
#Point Blue Conservation Science
#Website: http://www.michaelejohns.com/dataviz
#=====================================================================

library(tidyverse)
#read_csv does not work, so I used read.csv
data.f <- read.csv(file='7/SH20190130SST.csv') 
tail(data.f) #temp = NA when on land

#colors to use
col.p <- c("#B50012","#E60307","#F6421E","#F67F38","#FEC66E",
           "#FEFDA7","#EEF9CC","#B2DEEB","#6FA8CF","#3C65AB",
           "#2E3E96")

globe <- ggplot() +
  geom_point(data=data.f,
             aes(x=Var1, y=Var2, color=temp),
             shape=16) +
  
  #NA values for sea surface temperature are on land, hence gray
  scale_color_gradientn(colors=rev(col.p),na.value="grey30") +
  scale_x_continuous(expand = c(1,1)) +
  
  #orientation = where South Pole is, in latitude & longitude, 
  #and then the rotation (mid longitude in degrees of 
  #the map). Orthographic projection = earth as seen from 
  #infinitely far away, i.e. as a globe. 
  coord_map(projection='orthographic', 
            orientation=c(-90, -175, 0)) + 
  theme_minimal() + 
  theme(axis.title=element_blank(),
        axis.text=element_blank() )
globe #this takes some time: around 1-2 minutes on my computer

#gets sent to directory called Figs.
ggsave(filename='Figs\\SH20190130SST test.png', plot=globe, 
       width=6, height=6, dpi=300)

#=====================sf package======================================
#powerful package that allows for projections and many GIS functions
#including distance from points and lines, intersections
#and other such things only GIS software contains. 
#=====================================================================
library(units)
library(sf)
library(maps)

usa <- st_as_sf(x=maps::map('usa', plot = FALSE, fill = TRUE))
ggplot() + geom_sf(data = usa)

laea <- st_crs(x='+proj=laea +lat_0=30 +lon_0=-95') #Lambert equal area
usa <- st_transform(x=usa, crs=laea)
ggplot() + geom_sf(data = usa) + theme_minimal()

#makes the nc dataset available for North Carolina
demo(nc, ask = FALSE, echo = FALSE)

#what is inside the nc dataset?
nc

#plot chloropleth map for one column (BIR74) in the nc dataset
ggplot() + 
  geom_sf(data = nc, aes(fill = BIR74)) + #BIR74 is data column in nc
  scale_y_continuous(breaks = 34:36) +    #latitude tickmarks
  theme_minimal()


#********************ADDITIONAL CODE FOR INTEREST*********************
#not covered in lecture 

#===================ggmap using google maps terrain data==============
#for easy maps with terrain, taken from google maps
#From Eric Anderson
#http://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html
#=====================================================================
library(ggmap)
library(ggplot2)
# WARNING ----
# Paraphrased from SMEA 586 notes by Sunny Jardine
# To run this code, you need an API key from stadia maps to get 
# access to basemaps through ggmap. Do not share your API key with anyone 
# else. API keys can be obtained from Google but require a credit 
# card, this is not needed from Stadia Maps, but they limit the 
# number of requests you can make 
# https://docs.stadiamaps.com/authentication/
# Use your key here:
ggmap::register_stadiamaps("INSERT KEY HERE", write = TRUE)

randomlats <- runif(n=30, min=45, max=49)
randomlons <- runif(n=30, min=-122, max=-118)
randomdata <- data.frame(lats=randomlats, lons=randomlons)
wa_bbox <- make_bbox(data = randomdata, lat = lats, lon = lons)
wa_bbox

# grab terrain map from google
wa_terrain <- get_map(location = wa_bbox, source = 'google', 
                      maptype = 'terrain')
ggmap(wa_terrain) + 
  geom_point(data=randomdata, 
             mapping=aes(x=lons, y=lats),
             color='blue', size=4)
