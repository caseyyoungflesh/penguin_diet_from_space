######################
# 1a - Get centroids for all sites
######################


# load packages -----------------------------------------------------------

library(tidyverse)
library(sf)
library(here)


# set dir -----------------------------------------------------------------

dir <- paste0(here::here(), '/')
run_date <- Sys.Date()
sat_process_date <- '2025-07-16'


# load data ---------------------------------------------------------------

#load Landsat reflectances
guano_ret <- readRDS(paste0(dir, 'Data/L0/Landsat/landsat_reflectance_L457-', sat_process_date, '.rds'))


# get centroids -----------------------------------------------------------

#characteristics
site_char <- guano_ret %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(site_id) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::select(site_id, ccamlr_id, region)

site_cen <- guano_ret %>%
  dplyr::mutate(easting = sf::st_coordinates(geometry)[,1],
              northing = sf::st_coordinates(geometry)[,2]) %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(site_id) %>%
  dplyr::summarize(mn_easting = mean(easting),
                   mn_northing = mean(northing)) %>%
  dplyr::left_join(site_char, by = 'site_id')

site_ll <- sf::st_as_sf(site_cen, 
                        coords = c('mn_easting', 'mn_northing')) %>%
  #set crs
  sf::st_set_crs(3031) %>%
  #transform to ll
  sf::st_transform(4326) %>%
  sf::st_coordinates() %>%
  as.data.frame() %>%
  dplyr::select(X, Y) %>%
  dplyr::rename(longitude = X, latitude = Y)

site2 <- cbind(site_cen, site_ll)


# write to file -----------------------------------------------------------

write.csv(site2, paste0(dir, 'Data/L0/SiteLocations.csv'), row.names = FALSE)
