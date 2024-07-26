# All data in one data frame
# Import in the following order:
# srdb, spei, worldclim, modis
# data frame should have columns: long, lat, year, rs, MAT, MAP, NPP, SPEI for the current year, SPEI for the last year, and finally a boolean representing the presence of a birch effect
library(ggplot2)
library(dplyr)
library(tibble)
library(terra)
library(geodata)

# get SRDB
srdbData <- read.csv("srdb-data.csv")
srdbData %>% 
  select(Latitude, Longitude, Study_midyear, Rs_annual, Rs_growingseason, Ecosystem_type) %>% 
  filter(!is.na(Rs_annual) | !is.na(Rs_growingseason)) %>%
  # filter for only years between 1970 and 2000 since that's when the climate data is from
  filter(Study_midyear <= 2000 & Study_midyear >= 1970) ->
  srdb_filtered

# get SPEI
spei <- rast("spei12.nc")

# get worldclim
mat <- worldclim_global("tavg", "10", "worldclim_data/")
map <- worldclim_global("prec", "10", "worldclim_data/")

# Maybe average the monthly average temperatures so you get the yearly average?
# so somewhere like Thompson might have a yearly average of 5C whereas Lima might have a yearly average of like 20C
locs <- data.frame(Long = srdb_filtered$Longitude, Lat = srdb_filtered$Latitude)

matCoords <- terra::extract(mat, locs[1:2])
mapCoords <- terra::extract(map, locs[1:2])

# Remove useless ID column
matCoords <- matCoords[,-which(names(matCoords) == "ID")]
mapCoords <- mapCoords[,-which(names(mapCoords) == "ID")]

# takes the mean annual temperature and the mean annual precipitation

matMCoords <- rowMeans(matCoords, na.rm = TRUE)
mapMCoords <- rowSums(mapCoords, na.rm = TRUE)

birchEffectData <- data.frame(Longitude = srdb_filtered$Longitude, 
                              Latitude = srdb_filtered$Latitude,
                              Year = srdb_filtered$Study_midyear,
                              Annual_CO2_Flux = srdb_filtered$Rs_annual,
                              MAT = matMCoords,
                              MAP = mapMCoords
                              )

# As of 7/25/24, contains Long, Lat, Year, Flux, MAT, and MAP
head(birchEffectData)
