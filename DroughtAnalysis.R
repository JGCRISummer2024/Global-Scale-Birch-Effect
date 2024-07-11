library(ggplot2)
library(dplyr)
library(gridExtra)


srdbData <- read.csv("srdb-data.csv") # get data
srdbData %>% 
  select(Latitude, Longitude, Study_midyear, Rs_annual, Rs_growingseason, Study_midyear) %>% 
  filter(!is.na(Rs_annual) | !is.na(Rs_growingseason)) ->
  srdb_filtered # srdb data that has data in either rs_annual, rs_growingseason, or both
rsRange <- qplot(srdb_filtered$Rs_annual, geom = "boxplot", xlab = "Annual C flux values g / cm^2")

srdb_filtered %>% 
  filter(Rs_annual >= 0 & Rs_annual <= (1100 + 1.5 * (1100 - 486.5))) ->
  srdb_filtered_no_outliers

# world map with data, colored by annual C flux
world_coords <- map_data("world")
p <- ggplot(srdb_filtered_no_outliers, aes(Longitude, Latitude, color = Rs_annual)) + geom_map(data = world_coords, map = world_coords, aes(long, lat, map_id = region), color = "white", fill = "white", linewidth = 0.2) + geom_point()

# boxplot of the years the data was collected
yearPlot <- qplot(srdb_filtered$Study_midyear, geom = "boxplot", xlab = "Year of data collection")

# annual C flux boxplot without the extreme values from before
rsRangeNoOutliers <- qplot(srdb_filtered_no_outliers$Rs_annual, geom = "boxplot", xlab = "Annual C flux values g / cm^2 without extreme values")

# plot all the graphs
grid.arrange(p, yearPlot, rsRange, rsRangeNoOutliers, nrow = 2)



