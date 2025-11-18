#Useful options in map_data include usa, state, county, 
#world (0? in middle), and world2 (180? in middle)
alaska <- subset(map_data("world"), region == "USA")
alaska <- subset(map_data("world"), region == "USA" & long < -130 & lat > 50)
head(alaska)

#the group command explains which polygons belong together, one 
#polygon for each island on the map
ggplot(data=alaska, aes(x=long, y=lat)) + 
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