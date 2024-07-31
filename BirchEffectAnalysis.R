# All data in one data frame
# Import in the following order:
# srdb, spei, worldclim, modis
# data frame should have columns: long, lat, year, rs, MAT, MAP, NPP, SPEI for the current year, SPEI for the last year, and finally a boolean representing the presence of a birch effect
library(ggplot2)
library(dplyr)
library(tibble)
library(terra)
library(geodata)
library(tidyr)



# get SRDB
srdbData <- read.csv("srdb-data.csv")
srdbData %>% 
  select(Latitude, Longitude, Study_midyear, Rs_annual, Rs_growingseason, Ecosystem_type) %>% 
  # remove all NAs except for the ones in rs_growingseason
  filter(!is.na(Rs_annual)) %>% 
  filter(!is.na(Latitude)) %>% 
  filter(!is.na(Longitude)) %>%
  filter(!is.na(Study_midyear)) %>% 
  filter(Rs_annual >= 0) %>% 
  filter(Rs_annual <= 4500) %>% 
  filter(Latitude >= 0) ->
  srdb_filtered
srdb_filtered <- srdb_filtered[1:10,]
#make years easier to work with by elimiating the decimal
srdb_filtered$Study_midyear <- floor(srdb_filtered$Study_midyear)

# find growing season
# growing season is the average time from the last frost to the first one
# months where the average temperature is >0C?
# months where the average low is >0C?
# WorldClim only has minimum temp, not average min



# CACHING SPEI DATA
coords <- tibble(Longitude = srdb_filtered$Longitude, Latitude = srdb_filtered$Latitude)
coords %>% 
  distinct() ->
  coords
spei_dat_file <- paste("spei", nrow(coords), digest::digest(coords), sep = "_")
if(file.exists(spei_dat_file)) {
  message("Loading saved data ", spei_dat_file)
  spei_dat <- readRDS(spei_dat_file)
} else {
  message("Extracting SPEI data; this is slow...")
  # Documentation: https://spei.csic.es/database.html
  # Data downloaded 2024-07-17
  spei <- rast("spei12.nc")
  
  # Extract our points of interest. terra::extract() will handle making sure
  # the coordinates get mapped to the correct grid cell(s) in the data
  spei_dat <- terra::extract(spei, coords)
  saveRDS(spei_dat, file = spei_dat_file)
}
# Reshape data into a more manageable form
spei_monthly <- pivot_longer(spei_dat, -ID)
spei_monthly <- separate(spei_monthly, name, into = c("spei", "entry"), convert = TRUE)
# The SPEI data don't seem to provide 'time' explicitly in the netcdf, so
# compute it from the entries
spei_monthly$year <- ceiling(spei_monthly$entry / 12) + 1900
spei_monthly$month <- (spei_monthly$entry - 1) %% 12 + 1
spei_monthly$time <- with(spei_monthly, year + (month-1) / 12)

# Compute gsd - growing season drought
spei_monthly %>%
  arrange(ID, year, month) %>%
  group_by(ID, year) %>%
  summarise(gsd = mean(value[6:8]), .groups = "drop") ->
  gsd

# Merge back with coords
coords %>%
  mutate(ID = seq_len(nrow(coords))) %>%
  left_join(gsd, by = "ID") ->
  gsd

get_spei <- function(lon, lat, yr) {
  gsd %>% 
    dplyr::filter(Longitude == lon, Latitude == lat) %>% 
    arrange(year) -> x
  
  if(nrow(x) == 0) return(c(NA, NA, NA)) # not found
  
  y <- which(x$year == yr)
  x$gsd[c(y, y-1)]
}

message("Looking up growing season SPEI (yr 0, -1) for SRDB data...")
srdb_filtered$gsd1 <- srdb_filtered$gsd0 <- NA_real_
for(i in seq_len(nrow(srdb_filtered))) {
  spei_i <- get_spei(srdb_filtered$Longitude[i], srdb_filtered$Latitude[i], srdb_filtered$Study_midyear[i])
  srdb_filtered$gsd0[i] <- spei_i[1]
  srdb_filtered$gsd1[i] <- spei_i[2]
}

# SPEIcurr <- getSPEI()
# SPEIprev <- getSPEI(data.frame(lon = srdb_filtered$Longitude, lat = srdb_filtered$Latitude, year = floor(srdb_filtered$Study_midyear - 1)))
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

birchEffectData <- tibble(Longitude = srdb_filtered$Longitude, 
                              Latitude = srdb_filtered$Latitude,
                              Year = srdb_filtered$Study_midyear,
                              Annual_CO2_Flux = srdb_filtered$Rs_annual,
                              MAT = matMCoords,
                              MAP = mapMCoords,
                              spei_0 = srdb_filtered$gsd0,
                              spei_1 = srdb_filtered$gsd1
                              )

# As of 7/25/24, contains Long, Lat, Year, Flux, MAT, and MAP
head(birchEffectData)
