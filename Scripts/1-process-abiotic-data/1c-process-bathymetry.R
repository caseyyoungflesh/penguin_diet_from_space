#################
# 1c - Process bathymetry data
#
# Calculate shelf area (area less than 1000m depth)
#################


# load packages -----------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)
library(here)


# set dir -----------------------------------------------------------------

dir <- paste0(here::here(), '/')


# load data ---------------------------------------------------------------

#bath raster
#https://ibcso.org/
#https://doi.org/10.1038/s41597-022-01366-7
#500m resolution in IBCSO Polar Steroegraphic (9354)
br <- terra::rast(paste0(dir, 'Data/L0/Bathymetry_data/IBCSO_v2_bed.tif'))

#colony locations
sites <- read.csv(paste0(dir, 'Data/L0/SiteLocations.csv'))
sl <- sites %>%
  dplyr::select(site_id, mn_easting, mn_northing) %>%
  sf::st_as_sf(coords = c('mn_easting', 'mn_northing')) %>%
  #set crs
  sf::st_set_crs(3031) %>%
  #reproject to 9354
  sf::st_transform(9354) %>%
  #150km buffer
  sf::st_buffer(150000)


# extract -----------------------------------------------------------------

#get area of pixels that are < 1000m depth (shelf area) 0 in km^2
#1k meters is shelf break according to Clarke et al. 1998 Polar Biology
#few min run time for each
tt <- proc.time()
sh_ar <- terra::extract(br, terra::vect(sl),
                         fun = function(x) sum(x > -1000 & x < 0, na.rm = TRUE) * 0.25)
proc.time() - tt

#total number of below sea level pixels
tot_ar <- terra::extract(br, terra::vect(sl),
                         fun = function(x) sum(x < 0, na.rm = TRUE) * 0.25)

#merge with site names
sa_df <- data.frame(site = sl$site_id, 
                    sites,
                    shelf_area = round(sh_ar$elevation, 0), 
                    total_area = round(tot_ar$elevation, 0)) %>%
  dplyr::mutate(frac_shelf = round(shelf_area / total_area, 3)) %>%
  dplyr::rename(easting = mn_easting, northing = mn_northing) %>%
  dplyr::relocate(longitude, latitude, .after = northing)

#save to file
write.csv(sa_df, paste0(dir, 'Data/L1/shelf_area.csv'), row.names = FALSE)

