# SPEI
library(ggplot2)
library(dplyr)
library(tibble)
# function to get an SPEI given long, lat, and the number of months since jan 1900
getValue <- function(spatRaster, coords, time) { #spatRaster is spei
  loc <- tribble(
    ~lon, ~lat,
    coords[1], coords[2]
  )
  timeList <- terra::extract(spatRaster, loc)
  print(timeList[time])
}
coords <- tribble(
  ~place,  ~lon,         ~lat,
  "Boise City", -102.51464, 36.73183
)
coords$ID <- seq_len(nrow(coords))


# Documentation: https://spei.csic.es/database.html
# Data downloaded 2024-07-17

# These data are provided as netCDF files (http://en.wikipedia.org/wiki/NetCDF)
# Happily, the terra package can handle these!
# (There's also the `ncdf4` package for many other applications)
library(terra)
spei <- rast("spei12.nc")
getValue(spei, c(-102.75, 36.75), 1400)

# Note that we've told Git to IGNORE these data files; see ".gitignore" file

# Confirm that things look good -- it's a spatial raster object,
# global half degree resolution, from January 1901 to December 2022
print(spei)

# Extract our points of interest. terra::extract() will handle making sure
# the coordinates get mapped to the correct grid cell(s) in the data
spei_coords <- terra::extract(spei, coords[2:3])
# 1465 columns! Why? Because that's 122 years * 12 months + "ID" column

# Reshape data into a more manageable form
library(tidyr)
spei_monthly <- pivot_longer(spei_coords, -ID)
spei_monthly <- separate(spei_monthly, name, into = c("spei", "entry"), convert = TRUE)
# The SPEI data don't seem to provide 'time' explicitly in the netcdf, so
# compute it from the entries
spei_monthly$year <- ceiling(spei_monthly$entry / 12) + 1900
spei_monthly$month <- (spei_monthly$entry - 1) %% 12 + 1
spei_monthly$time <- with(spei_monthly, year + (month-1) / 12)
spei_monthly <- dplyr::left_join(spei_monthly, coords, by = "ID")

# Not sure about the best way to plot this
p <- ggplot(spei_monthly, aes(time, value, color = place)) +
  geom_point(size = 0.25) + geom_area() +
  scale_x_discrete(name = "Year", limits = seq(from = 1900, to = 2020, by = 10)) +
  scale_y_discrete(name = "SPEI", limits = seq(from = -4, to = 4, by = 1))
print(p)


