#installing marmap can take 30-60 seconds
#make sure R is up to date (v 4.0.x) with "version" in console
library(raster) #may need to reinstall
library(marmap)
library(tidyverse)
library(dplyr)
library(maps)
library(ggplot2)
library(sf)

# Get Alaska polygon data correctly
alaska <- subset(map_data("world"), region == "USA")
alaska <- subset(map_data("world"), region == "USA" & long < -130 & lat > 50)


ggplot(data = alaska, aes(x = long, y = lat, group = group)) + 
  geom_polygon(fill = "lightblue", color = "black") +
  coord_fixed(1.3) +
  theme_minimal()

mean_lat <- mean(alaska$lat)
ratio <- 1 / cos(mean_lat * pi / 180)

ggplot(data = alaska, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "lightblue", color = "black") +
  coord_fixed(ratio = ratio) +
  theme_minimal()


#get the bathymetry from NOAA. Resolution in minutes, default = 4, 
#minimum = 1. Smaller = higher resolution and longer times
#keep=TRUE writes downloaded data into a file
bathy <- getNOAA.bathy(lon1=-135, lon2=-180, lat1=50, lat2=74,
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




## EXAMPLE CODE from STACKOVERFLOW


# Example sample data: your samples by region (make sure regions are lowercase to match map_data)
# Loads in data
combined_df <- read.csv("Deliverables/allrows-abundance.csv")

# Reads in shapefile
alaska_coastal_regions <- st_read("path/to/alaska_coastal_regions.shp")

# Plot 
ggplot(alaska_coastal_regions) +
  geom_sf() +
  theme_minimal()


# Join sample data to map polygons by region
map_with_samples <- left_join(states_map, sample_data, by = "region")

# Calculate centroids of each region for plotting points on the map
centroids <- map_with_samples %>%
  group_by(region) %>%
  summarise(
    centroid_long = mean(range(long)),
    centroid_lat = mean(range(lat)),
    sample_count = first(sample_count)  # sample data joined; NA if no data
  )

# Plot map polygons with fill based on sample counts (NA regions will be blank)
ggplot(map_with_samples, aes(x = long, y = lat, group = group, fill = sample_count)) + 
  geom_polygon(color = "black") +
  scale_fill_viridis_c(na.value = "grey90", option = "plasma") +
  coord_fixed(1.3) +
  theme_minimal() +
  geom_point(data = centroids, aes(x = centroid_long, y = centroid_lat, size = sample_count), color = "red") +
  labs(fill = "Sample count", size = "Sample count")
