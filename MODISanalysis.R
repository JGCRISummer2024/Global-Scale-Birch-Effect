# MODIS satellite data
## NOT SURE WHY THIS ISNT WORKING, FIGURE OUT TOMORROW
# Compute daily NPP
# We're interested in 2022 data for the TEMPEST site;
# currently using lat/lon of the NEON tower and pulling
# data for a roughly 1 km2 area around that.
# See https://rspatial.org/modis/1-introduction.html
# BBL November 2023

message("Welcome to modis.R")

# install.packages('luna', repos='https://rspatial.r-universe.dev')
library(terra)
library(luna)

# Dates of interest
start <- "2022-01-01"
end <- "2022-12-31"

# Area of interest
# I selected a single point, the coordinates of the NEON tower
# 38.890131, -76.560014, and then expand aoi slightly around it
aoi <- ext(-76.57, -76.56, 38.89, 38.90) # roughly 1 km2 around tower
aoi_vect <- terra::vect(aoi)
# EPSG:4326 - WGS 84 is lat/lon coordinate system used by GPS and others
crs(aoi_vect) <- "epsg:4326"
message("Extract area is ",
        sprintf("%4.2f", terra::expanse(aoi_vect) / 1e6), " km2")

# GPP

# SDS Name: Gpp_500m
# Description: Gross Primary Productivity
# Units: kg C/m²
# 16-bit signed integer
# Fill Value: 32761 to 32767
# No Data Value: N/A
# Valid Range: 0 to 30000
# Scale Factor: 0.0001
product <- "MOD17A2H"

tiles <- luna::getNASA(product,
                        start, end, aoi = aoi,
                        download = FALSE)

message("This will require ", length(tiles), " tiles")

download <- TRUE
if(download) {
  # Save a credentials file by
  earthdata <- list(username = "username", password = "password")
  saveRDS(earthdata, file = "modis_data/earthdata_credentials")
  
  # Read in the credentials file so we can log in and download
  credfile <- "modis_data/earthdata_credentials.RDS"
  if(file.exists(credfile)) {
    earthdata <- readRDS(credfile)
  } else {
    stop("No credentials for downloading! See comments")
  }
  
  mf <- luna::getModis(product,
                       start, end, aoi = aoi,
                       download = TRUE,
                       path = "./modis_data/",
                       username = earthdata$username,
                       password = earthdata$password)
} else {
  message("Not downloading; hope everything is here")
}

# Process downloaded files
files <- list.files("modis_data/",
                    pattern = "^MOD17A2H.*hdf$",
                    full.names = TRUE)

# Read in data files and compute spatial average
res <- list()
for(f in sort(files)) {
  message("Processing ", f)
  r <- rast(f)
  # Isolate GPP data
  gpp <- r[["Gpp_500m"]]
  # Rescale and 'clamp' the values, based on documentation:
  # https://lpdaac.usgs.gov/products/myd17a2hv061/
  gpp_clamp <- clamp(gpp / 0.0001, lower = 0, upper = 30000, values = FALSE) * 0.0001
  
  # Make sure our area of interest (aoi) is in same coordinate
  # system and then extract data for that region
  aoi_vect_trans <- project(aoi_vect, crs(gpp_clamp))
  gpp_ex <- terra::extract(gpp_clamp, aoi_vect_trans)
  
  # Compute mean in gC/m2 (per 8 day period)
  gpp_mean <- mean(gpp_ex$Gpp_500m, na.rm = TRUE) * 1000 # gC/m2
  res[[f]] <- data.frame(filename = basename(f),
                         GPP = gpp_mean)
}

results <- do.call(rbind, res)
rownames(results) <- NULL

# Separate filename into more useful data. From the help page:

# In this example of a tiled product, the filename
# MCD19A2.A2023081.h22v02.061.2023082164002.hdf indicates:
#
# MCD19A2 - Product Short Name
# A2023081 - Julian Date of Acquisition (AYYYYDDD)
# h22v02 - Tile Identifier (horizontalXXverticalYY)
# 061 - Collection Version
# 2023082164002 - Julian Date of Production (YYYYDDDHHMMSS)
# .hdf - Data Format (HDF-EOS)

library(tidyr)
results %>%
  separate(filename,
           into = c("product", "jdoa", "ti", "cv", "jdop", "df"),
           sep = "\\.") ->
  gpp

gpp$year <- as.integer(substr(gpp$jdoa, 2, 5))
gpp$day <- as.integer(substr(gpp$jdoa, 6, 8))
gpp <- gpp[gpp$year == 2022,]

# MODIS GPP is an 8-day product. Convert to daily GPP
gpp$daily_GPP <- gpp$GPP / 8

# Interpolate to every day for the entire year
gpp_interp <- data.frame(day = 1:365,
                         GPP = approx(gpp$day,
                                      gpp$daily_GPP,
                                      xout = 1:365, rule = 2)$y)


# NPP

# SDS Name: Npp_500m
# Description: Net Primary Productivity
# Units: kg C/m²/yr
# 16-bit signed integer
# Fill Value: 32761 to 32767
# No Data Value: N/A
# Valid Range: -30000 to 32700
# Scale Factor: 0.0001

# NPP tile downloaded 2023-11-02 from NASA Earthdata
# (it doesn't seem to be available through the luna package)
# MOD17A3HGF.A2022001.h12v05.061.2023035051756.hdf
files <- list.files("modis_data/",
                    pattern = "^MOD17A3HGF.*hdf$",
                    full.names = TRUE)
stopifnot(length(files) == 1)

message("Processing ", files[1])
r <- rast(files[1])

npp <- r[["Npp_500m"]]
npp_clamp <- clamp(npp / 0.0001, lower = -30000, upper = 32700, values = FALSE) * 0.0001

# Make sure our area of interest (aoi) is in same coordinate
# system and then extract data
aoi_vect_trans <- project(aoi_vect, crs(npp_clamp))
npp_ex <- terra::extract(npp_clamp, aoi_vect_trans)

npp_mean <- mean(npp_ex$Npp_500m, na.rm = TRUE) * 1000 # gC/m2
message("NPP is ", npp_mean)

# Scale the annual NPP by the GPP curve, assuming a fixed ratio
gpp_interp$NPP <- npp_mean * gpp_interp$GPP / sum(gpp_interp$GPP)

write.csv(gpp_interp, "modis_data/modis_gpp_npp_2022.csv")
message("All done")

