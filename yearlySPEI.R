library(terra)
library(lubridate)
spei <- rast("spei12.nc")[[1:24]]

# Make annual composite by averages per month
annualSPEI <- function(instack, outname){ # instack = raster stack (SPEI), outname = what to output
  # Make list of month groupings
  date_ls <- seq(as.Date('1902/1/15'), by = 'month', length.out = nlyr(instack))
  annual_composite_grp = year(as.Date(date_ls))
  # Take mean, excluding NAs  (should we consider NAs as Zeros?)
  f <- function(v){ # v is a raster
    tapply(v, annual_composite_grp, mean)
  }
  annual_mean_composite = app(instack, f)
  writeRaster(annual_mean_composite, outname, overwrite=TRUE) # makes a .tif file, open with canopy app (nasa)
}

annualSPEI(spei, "speiAnnualTest.tif")
speiAnnual <- rast("speiAnnualTest.tif")
plot(speiAnnual)
plot(spei[[13]])

# next steps
# run for entire time series (1902 - 2022)
# apply extract funtion to annual raster
# have SPEI and SRDB in the same data frame
# questions:
## average Z score?
# https://www.qgis.org/

