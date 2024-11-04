library(ggplot2)
library(dplyr)
library(tibble)
library(terra)
library(geodata)
library(tidyr)

#### get SRDB ####
srdbData <- read.csv("srdb-data.csv")
srdbData %>% 
  select(Latitude, Longitude, Study_midyear, Rs_annual) %>% 
  # remove all NAs except for the ones in rs_growingseason
  filter(!is.na(Rs_annual)) %>% 
  filter(!is.na(Latitude)) %>% 
  filter(!is.na(Longitude)) %>%
  filter(!is.na(Study_midyear)) %>% 
  filter(Rs_annual >= 0) %>% 
  filter(Rs_annual <= 4500) ->
  srdb_filtered
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
# get below 28F (limit of the growing season according to the USDA)

#### Get WorldClim ####

mat <- worldclim_global("tavg", "10", "worldclim_data/")
map <- worldclim_global("prec", "10", "worldclim_data/")

locs <- tibble(Longitude = srdb_filtered$Longitude, Latitude = srdb_filtered$Latitude)

matCoords <- terra::extract(mat, locs)
mapCoords <- terra::extract(map, locs)

# takes the mean annual temperature and the mean annual precipitation
matMCoords <- rowMeans(matCoords[2:13])
mapMCoords <- rowSums(mapCoords[2:13])

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

x <- terra::rast("RenderData.tiff")
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
                              spei_measured_year = srdb_filtered$gsd0,
                              spei_previous_year = srdb_filtered$gsd1,
                              birch_effect_potential = (spei_measured_year >= 0 & spei_previous_year <= -1)
                              )

modelFlag <- lm(Annual_CO2_Flux ~ MAT + MAP + NPP + birch_effect_potential, data = birchEffectData)
message(paste("Printing birch effect analysis, N=", nrow(birchEffectData)))
print(car::Anova(modelFlag, type = "III"))

# wet locations
birchEffectData %>% 
  filter(MAP > 2000) ->
  birchEffectDataWet
message(paste("Printing only wet locations (MAP > 2000mm), N=", nrow(birchEffectDataWet)))
modelFlagWet <- lm(Annual_CO2_Flux ~ MAT + MAP + NPP + birch_effect_potential, data = birchEffectDataWet)
print(car::Anova(modelFlagWet, type = "III"))

# arid/semi-arid locations excluding almost completely rainless areas
birchEffectData %>% 
  filter(MAP < 750) %>% 
  filter(MAP > 100) ->
  birchEffectDataDry
message(paste("Printing only dry locations (MAP < 750mm), N=", nrow(birchEffectDataDry)))
modelFlagDry <- lm(Annual_CO2_Flux ~ MAT + MAP + NPP + birch_effect_potential, data = birchEffectDataDry)
print(car::Anova(modelFlagDry, type = "III"))

# World map
world_coords <- map_data("world")
p <- ggplot() + geom_map(data = world_coords, map = world_coords, aes(long, lat, map_id = region), color = "gray", fill = "white") + geom_point(data = birchEffectData, aes(Longitude, Latitude, color = Annual_CO2_Flux), alpha = .5)
print(p)


# Climatology Graph
presentationLoc <- tribble(
  ~place, ~lon, ~lat,
  "Richland", -119.273, 46.282
)

presentationLoc$ID <- seq_len(nrow(presentationLoc))

tavg_coords <- terra::extract(mat, presentationLoc[2:3])
precip_coords <- terra::extract(map, presentationLoc[2:3])

monthly <- pivot_longer(tavg_coords, -ID)
monthly$month <- as.numeric(gsub("wc2.1_10m_tavg_", "", monthly$name))
monthly <- dplyr::left_join(monthly, presentationLoc)

precip_coords %>% 
  pivot_longer(-ID) ->
  precip_coords

monthly$precip <- precip_coords$value
custom_colors <- c("Temperature" = "red", "Precipitation" = "lightblue")
p <- ggplot(monthly) + 
  geom_bar(aes(x = month, y = precip, fill = "Precipitation"), stat = "identity", color = "lightblue") + 
  geom_line(aes(x = month, y = (value * 9/5) + 32, color = "Temperature"), stat = "identity", linewidth = 2) + 
  ylab("Mean Air temperature (Fahrenheit)") + 
  scale_y_continuous(sec.axis = sec_axis(transform = ~./25.4, name = "Rainfall (Inches)")) + 
  ggtitle("Climatology in Richland, WA") + 
  scale_x_discrete(name = "Month", limits = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) + 
  scale_color_manual(values = custom_colors, name = "") +
  scale_fill_manual(values = custom_colors, name = "")
print(p)
