srdbData <- read.csv("srdb-data.csv")
library(dplyr)
srdbData %>% 
  select(Latitude, Longitude, Study_midyear, Rs_annual, Rs_growingseason) %>% 
  filter(!is.na(Rs_annual) | !is.na(Rs_growingseason)) ->
  srdb_filtered

# Sanity check: plot the data
library(ggplot2)
world_coords <- map_data("world")
p <- ggplot(srdb_filtered, aes(Longitude, Latitude, color = Rs_annual)) + geom_map(data = world_coords, map = world_coords, aes(long, lat, map_id = region), color = "white", fill = "white", linewidth = 0.2) + geom_point()
print(p)