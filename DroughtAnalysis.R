srdbData <- read.csv("srdb-data.csv")
library(dplyr)
srdbData %>% 
  select(Latitude, Longitude, Study_midyear, Rs_annual, Rs_growingseason) %>% 
  filter(!is.na(Rs_annual) | !is.na(Rs_growingseason)) ->
  srdb_filtered

# Sanity check: plot the data
library(ggplot2)
p <- ggplot(srdb_filtered, aes(Longitude, Latitude, color = Rs_annual)) + geom_point()
print(p)