######################
# 1b - Process SIC data
#
######################


# load packages -----------------------------------------------------------

library(rvest)
library(tidyverse)
library(terra)
library(sf)
library(here)


# set dir -----------------------------------------------------------------

dir <- paste0(here::here(), '/')


# DL SIC data -----------------------------------------------------------------

#dl from NSIDC
#https://nsidc.org/data/user-resources/help-center/how-access-and-download-noaansidc-data
#v4 was used for analysis
# url <- 'https://noaadata.apps.nsidc.org/NOAA/G02202_V4/south/monthly/'
#v5 is current release
url <- 'https://noaadata.apps.nsidc.org/NOAA/G02202_V5/south/monthly/'


page <- rvest::read_html(url)
files <- page %>% 
  rvest::html_nodes("a") %>% 
  rvest::html_attr('href')

for (i in 1:length(files))
{
  #i <- 2
  print(paste0('downloading file ', i, ' of ', length(files)))
  u <- paste(url, files[i], sep = '/')
  download.file(u, paste0(dir, 'Data/L0/SIC/', files[i]), mode = 'wb')
}
  

# load data ---------------------------------------------------------------

#read in ncdfs using terra (stacked - climate data record)
#25km grid
#NSIDC Sea Ice Polar Stereographic South
#EPSG::3412
lf <- list.files(paste0(dir, 'Data/L0/SIC/'), full.names = TRUE)
sic_stack <- terra::rast(lf, subds = 'cdr_seaice_conc_monthly')
terra::crs(sic_stack) <- 'epsg:3412'

#get years/months
dates <- paste0('date_', substring(lf, nchar(lf[1]) - 19, nchar(lf[1]) - 14))
# months <- substring(lf, nchar(lf[1]) - 15, nchar(lf[1]) - 14)

#colony locations
sites <- read.csv(paste0(dir, 'Data/L0/SiteLocations.csv'))
sl_150 <- sites %>%
  dplyr::select(site_id, mn_easting, mn_northing) %>%
  sf::st_as_sf(coords = c('mn_easting', 'mn_northing')) %>%
  #set crs
  sf::st_set_crs(3031) %>%
  #reproject to 3412
  sf::st_transform(3412) %>%
  #150km buffer
  sf::st_buffer(150000)


# extract -----------------------------------------------------------------

colony_sic_150 <- terra::extract(sic_stack, terra::vect(sl_150), 
                          function(x) mean(x, na.rm = TRUE))
colnames(colony_sic_150) <- c('ID', dates)

#merge with site names and switch to long format
sic_df <- data.frame(sites,
                     colony_sic_150[,-1]) %>%
  tidyr::pivot_longer(cols = tidyselect::starts_with('date'),
                      names_to = 'date',
                      values_to = 'sic_150km') %>%
  dplyr::mutate(year = as.numeric(substring(date, 6, 9)),
                month = as.numeric(substring(date, 10, 11))) %>%
  dplyr::select(-date) %>%
  dplyr::rename(easting = mn_easting, 
                northing = mn_northing) %>%
  dplyr::relocate(sic_150km, .after = month) %>%
  dplyr::relocate(longitude, latitude, .after = northing)

#save to file
write.csv(sic_df, paste0(dir, 'Data/L1/colony_sic.csv'), row.names = FALSE)
