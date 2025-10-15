######################
# 3 - Predict diet values based on PLSR model and spectra
#
# save out diet predictions (into L3)
######################


# load packages -----------------------------------------------------------

library(tidyverse)
library(sf)
library(here)


# set dir -----------------------------------------------------------------

dir <- paste0(here::here(), '/')
run_date <- Sys.Date()
sat_process_date <- '2025-07-16'
plsr_run_date <- '2024-12-06'


# load data ---------------------------------------------------------------

#plsr ensemble model coefs
L4_coef_en <- readRDS(paste0(dir, 'Data/L2/', plsr_run_date, '/L4_coef_en.rds'))
L5_coef_en <- readRDS(paste0(dir, 'Data/L2/', plsr_run_date, '/L5_coef_en.rds'))
L7_coef_en <- readRDS(paste0(dir, 'Data/L2/', plsr_run_date, '/L7_coef_en.rds'))

#load Landsat reflectances
guano_ret <- readRDS(paste0(dir, 'Data/L0/Landsat/landsat_reflectance_L457-', sat_process_date, '.rds'))

#all good pixels
unique(guano_ret$bad)

#all guano
unique(guano_ret$guano_bin)


# predict dN --------------------------------------------------------------

#predict dN function
#filter for relevant fields

pred_dN_fun <- function(data, coef)
{
  #filter for relevant fields
  LX_ret_bands <- dplyr::select(data, blue, green, red, nir, swir1, swir2) %>%
    sf::st_drop_geometry()
  
  #predict
  #to predict (next script)
  LX_pred_en_ch <- matrix(NA,
                          nrow = NROW(coef),
                          ncol = NROW(LX_ret_bands))
  for (i in 1:NROW(coef))
  {
    #i <- 1
    LX_pred_en_ch[i,] <- as.matrix(cbind(1, LX_ret_bands)) %*% 
      as.numeric(coef[i,]) 
  }

  return(LX_pred_en_ch)
}

#run function for each sensor
guano_ret4 <- dplyr::filter(guano_ret, satellite == 4)
guano_ret5 <- dplyr::filter(guano_ret, satellite == 5)
guano_ret7 <- dplyr::filter(guano_ret, satellite == 7)

L4_pred_ch <- pred_dN_fun(data = guano_ret4, 
                       coef = L4_coef_en)
L5_pred_ch <- pred_dN_fun(data = guano_ret5, 
                       coef = L5_coef_en)
L7_pred_ch <- pred_dN_fun(data = guano_ret7, 
                       coef = L7_coef_en)


# combine data ------------------------------------------------------------

#add dN
guano_ret4$dN_mn <- apply(L4_pred_ch, 2, mean)
guano_ret5$dN_mn <- apply(L5_pred_ch, 2, mean)
guano_ret7$dN_mn <- apply(L7_pred_ch, 2, mean)

guano_ret4$dN_sd <- apply(L4_pred_ch, 2, sd)
guano_ret5$dN_sd <- apply(L5_pred_ch, 2, sd)
guano_ret7$dN_sd <- apply(L7_pred_ch, 2, sd)

gr_comb <- rbind(guano_ret4, guano_ret5, guano_ret7)


# add easting/northing ----------------------------------------------------

gr_comb2 <- gr_comb %>%
  dplyr::mutate(easting = sf::st_coordinates(geometry)[,1],
                northing = sf::st_coordinates(geometry)[,2]) %>%
  sf::st_drop_geometry() %>%
  dplyr::select(-bad, -guano_bin)


# write out ---------------------------------------------------------------

saveRDS(gr_comb2, paste0(dir, 'Data/L3/diet_pred-', run_date, '.rds'))

