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

# lowTemp <- tribble(
#   ~city, ~minTemp, ~avgTemp,
#   "seattle", 31.3, 47.1,
#   "dc", 34.9, 58.2,
#   "new york", 32.8, 53.7,
#   "bangor", 30.8, 54.5,
#   "fargo", 29.1, 56.6,
#   "portland", 30.4, 48.3,
#   "atlanta", 28.1, 55.6,
#   "nashville", 32.7, 60.8,
#   "minneapolis", 35.7, 59.5,
#   "boise", 33, 59.9,
#   "burlington", 32.2, 58.4
# )
#plot(lowTemp$minTemp, lowTemp$avgTemp)
#tempmodel <- lm(avgTemp ~ minTemp, data = lowTemp)
#print(summary(tempmodel))
# regression shows the predicted mean temp for a month with a mean min of 28F is 52F (R^2 = 0.19)

#### get SRDB ####
srdbData <- read.csv("srdb-data.csv")
srdbData %>% 
  select(Latitude, Longitude, Study_midyear, Rs_annual, Rs_growingseason, Ecosystem_type) %>% 
  # remove all NAs except for the ones in rs_growingseason
  filter(!is.na(Rs_annual)) %>% 
  filter(!is.na(Latitude)) %>% 
  filter(!is.na(Longitude)) %>%
  filter(!is.na(Study_midyear)) %>% 
  filter(Rs_annual >= 0) %>% 
  filter(Rs_annual <= 4500) ->
  # TODO: maybe filter by latitude, or make it so it gets a different growing season for different places
  srdb_filtered
#srdb_filtered <- srdb_filtered[1:10,] # TODO: delete when done (for testing)
#make years easier to work with by elimiating the decimal
srdb_filtered$Study_midyear <- floor(srdb_filtered$Study_midyear)





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
  message("Extracting SPEI data; this will take 2+ hours...")
  # Documentation: https://spei.csic.es/database.html
  # Data downloaded 2024-07-17
  spei <- rast("spei01.nc")
  
  # Extract our points of interest. terra::extract() will handle making sure
  # the coordinates get mapped to the correct grid cell(s) in the data
  spei_dat <- terra::extract(spei, coords)
  saveRDS(spei_dat, file = spei_dat_file)
}
#### Get SPEI ####
# Reshape data into a more manageable form
spei_monthly <- pivot_longer(spei_dat, -ID)
spei_monthly <- separate(spei_monthly, name, into = c("spei", "entry"), convert = TRUE)
# The SPEI data don't seem to provide 'time' explicitly in the netcdf, so
# compute it from the entries
spei_monthly$year <- ceiling(spei_monthly$entry / 12) + 1900
spei_monthly$month <- (spei_monthly$entry - 1) %% 12 + 1
spei_monthly$time <- with(spei_monthly, year + (month-1) / 12)


# find growing season
# growing season is the average time from the last frost to the first one
# since I don't have average min temp data from WorldClim, I will use
# months that average >10C, which should roughly coincide with the last month to
# get below 28F

#### Get WorldClim ####
mat <- worldclim_global("tavg", "10", "worldclim_data/")
map <- worldclim_global("prec", "10", "worldclim_data/")

locs <- tibble(Longitude = srdb_filtered$Longitude, Latitude = srdb_filtered$Latitude)

matCoords <- terra::extract(mat, locs)
mapCoords <- terra::extract(map, locs)

# takes the mean annual temperature and the mean annual precipitation
matMCoords <- rowMeans(matCoords, na.rm = TRUE)
mapMCoords <- rowSums(mapCoords, na.rm = TRUE)

# make temperature data into a nicer format to prepare for merging into spei_monthly
matCoords %>% 
  pivot_longer(-ID) %>% 
  separate(name, into = c("wc2", "month"), sep = "tavg_", convert = TRUE) ->
  matCoords
matCoords <- matCoords[-2] # removes useless wc2 column

# attaches temp data to spei_monthly
spei_monthly %>% 
  left_join(matCoords, by = c("ID", "month")) ->
  spei_monthly


# Compute gsd - growing season drought
# this dynamically uses average temperature data to average the SPEI in months
# that average above 10C

spei_monthly %>%
  arrange(ID, year, month) %>%
  group_by(ID, year) %>%
  summarise(gsd = mean(value.x[which(value.y > 10)]), .groups = "drop") ->
  gsd
# mean(value[as.numeric(speis of months where the average is above 10C)], na.rm = TRUE)
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

#### Get NPP ####

x <- terra::rast("~/GitHub/CO2-variability-with-droughts/RenderData.tiff")
y <- terra::extract(x, 1:(360*720))
# 255 is the NA value, so we'll remove those
x_clamp <- clamp(x, lower = 1, upper = 254, values = FALSE)
# rescale the data to fit the scale on the website
x_clamp <-  x_clamp*(1950/254) + 50
net_primary_production <- terra::extract(x_clamp, locs)

#### Analysis ####
birchEffectData <- tibble(Longitude = srdb_filtered$Longitude, 
                              Latitude = srdb_filtered$Latitude,
                              Year = srdb_filtered$Study_midyear,
                              Annual_CO2_Flux = srdb_filtered$Rs_annual,
                              MAT = matMCoords,
                              MAP = mapMCoords,
                              NPP = net_primary_production$RenderData,
                              spei_0 = srdb_filtered$gsd0,
                              spei_1 = srdb_filtered$gsd1,
                              # TODO: change flag, maybe
                              spei_flag = (spei_0 >= 0 & spei_1 <= -0.5)
                              )
modelNoFlag <- lm(Annual_CO2_Flux ~ MAT + MAP + NPP + spei_0 + spei_1, data = birchEffectData)
print(car::Anova(modelNoFlag, type = "III"))
modelFlag <- lm(Annual_CO2_Flux ~ MAT + MAP + spei_0 + spei_1 + spei_flag, data = birchEffectData)
print(car::Anova(modelFlag, type = "III"))
print(birchEffectData)
