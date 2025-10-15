######################
# 4 - diet est
#
######################


# Load packages -----------------------------------------------------------

library(MCMCvis)
library(tidyverse)
library(cmdstanr)
library(here)


# set dir -----------------------------------------------------------------

dir <- paste0(here::here(), '/')
run_date <- Sys.Date()
diet_pred_run_date <- '2025-07-16'


# read in data -------------------------------------------------------------

#diet prediction
dpp1 <- readRDS(paste0(dir, 'Data/L3/diet_pred-', diet_pred_run_date, '.rds')) %>%
  dplyr::filter(!is.na(site_id)) %>%
  dplyr::mutate(date = paste0(year, '-', month, '-', day)) %>%
  #num pixels each site, year, doy
  dplyr::group_by(site_id, season, month, day) %>%
  dplyr::mutate(n_pix = n()) %>%
  dplyr::ungroup() %>%
  #avg num pix each site
  dplyr::group_by(site_id) %>%
  dplyr::mutate(avg_n_pix = mean(n_pix)) %>%
  dplyr::ungroup() %>%
  #only colonies with at least 2 pixels on average
  dplyr::filter(avg_n_pix >= 2)
  
dp <- dpp1


# read in lat/lon -------------------------------------------------

sites <- read.csv(paste0(dir, 'Data/L0/SiteLocations.csv'))


# prediction uncertainty ----------------------------------------------------------

#r2 = 1 - (resid var / total var)
#resid var / total var = (1 - r2)
#(1 - r2) * total var = resid var

#r2 from model prediction - 2-PLSR.R
tau <- sqrt((1 - 0.45) * var(dp$dN_mn))


# get DOY ------------------------------------------------------------

#days from Nov 1 for each 'season'
#doy for each day
tday <- lubridate::yday(lubridate::ymd(dp$date))
#doy for Nov 1 for associated season minus 1
n1_doy <- (lubridate::yday(lubridate::ymd(paste0(dp$season, 
                                                 '-11-01'))) - 1)
#number days between Nov 1 and end of year
n1_delta <- ((lubridate::yday(lubridate::ymd(paste0(dp$season, 
                                                    '-12-31'))) - 1) - n1_doy)
#year and season the same
gt180 <- which(tday > 180)
#year one greater than season
lt180 <- which(tday < 180)

#fill blank vector
doy <- rep(NA, NROW(dp))
doy[gt180] <- tday[gt180] - n1_doy[gt180]
doy[lt180] <- tday[lt180] + n1_delta[lt180]

#add to df
dp$doy <- doy


# indexing ----------------------------------------------------------------

dp2 <- dplyr::arrange(dp, site_id, season, doy, cell) %>%
  dplyr::group_by(site_id, season, doy) %>%
  dplyr::mutate(syd_idx = dplyr::cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(site_id, season) %>%
  dplyr::mutate(sy_idx = dplyr::cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(site_id) %>%
  dplyr::mutate(s_idx = dplyr::cur_group_id()) %>%
  dplyr::ungroup()

d_sy <- dplyr::select(dp2, site_id, ccamlr_id, region,
                      season, doy, n_pix, avg_n_pix,
                      syd_idx, sy_idx, s_idx) %>%
  dplyr::distinct() %>%
  dplyr::arrange(syd_idx)


# explore data ------------------------------------------------------------

#Nov - Feb
unique(dp2$month)

#years
range(dp2$year)
range(dp2$season)

#unique scenes - 1727 (some scenes multiple colonies)
length(unique(dp2$Landsat.Product.Identifier.L1))

#unique colonies - 119
length(unique(dp2$site_id))

#unique site/days - 3938
NROW(dplyr::distinct(dp2, site_id, year, month, day))

#number of pixels - 179,346
NROW(dp2)

#number of unique days per season for each colony
# 1-21
nud <- dplyr::distinct(dp2, site_id, season, doy) %>%
  dplyr::group_by(site_id, season) %>%
  dplyr::summarize(N = n()) %>%
  dplyr::ungroup()
range(nud$N)

#number of total captures for each colony
# 1-128
ntc <- dplyr::group_by(nud, site_id) %>%
  dplyr::summarize(N_total = sum(N))
range(ntc$N_total)


# fit model ---------------------------------------------------------------

DATA_est <- list(N = NROW(dp2),
             Nsyd = length(unique(dp2$syd_idx)),
             Ns = length(unique(dp2$s_idx)),
             y = dp2$dN_mn,
             sd_y = dp2$dN_sd,
             tau = tau,
             syd_idx = dp2$syd_idx,
             s_idx2 = d_sy$s_idx, #s_idx for each site/year/DOY
             pro_data = dp2)

options(mc.cores = parallel::detectCores())

# DELTA <- 0.92
# TREE_DEPTH <- 10
# STEP_SIZE <- 0.005
CHAINS <- 4
ITER <- 9000

#compile model
est_mod <- cmdstanr::cmdstan_model(paste0(dir, 
                                      'Scripts/model_files/4-diet-est.stan'),
                               force_recompile = TRUE)

print("STARTING TO FIT MODEL")

#sample
est_fit <- est_mod$sample(
  data = DATA_est,
  chains = CHAINS,
  iter_sampling = ITER / 2,
  iter_warmup = ITER / 2,
  parallel_chains = CHAINS,
  refresh = 500)
# step_size = STEP_SIZE
# max_treedepth = TREE_DEPTH
# adapt_delta = DELTA


# subset model object -----------------------------------------------------

# function to subset cmdstanr object
# model_fit is the fit object from `cmdstanr`
# `target_params` argument ignores square brackets (i.e,. []).
# `target_params` = 'alpha' will return 'alpha[1]', 'alpha[2]', etc.
# `thin` argument saves every X iteration
subset_cmdstanr <- function(model_fit, target_params, thin = 1)
{
  #get list of temp cmdstan files
  cmdstan_files <- model_fit$output_files()
  
  #get param names
  all_params <- model_fit$metadata()$model_params
  
  #remove []
  all_params_ISB <- vapply(strsplit(all_params,
                                    split = "[", fixed = TRUE), 
                           `[`, 1, FUN.VALUE=character(1))
  
  #get idx for target params
  f_ind <- which(all_params_ISB %in% target_params)
  
  #add 6 to all indices except 1, if it exists to account for:
  #accept_stat__,stepsize__,treedepth__,n_leapfrog__,divergent__,energy__
  #which are in csv but not all_params
  #add 1 to beginning of vec if not there
  f_ind2 <- f_ind
  if (sum(f_ind2 == 1))
  {
    f_ind2[-1] <- f_ind2[-1] + 6
  } else {
    f_ind2 <- c(1, f_ind2 + 6)
  }
  
  #add 2-7 (metadata)
  f_ind2 <- append(f_ind2, 2:7, after = 1)
  
  #when feeding >0.5 million args (columns) to cut, there are issues when passing
  #too many indices. Get areas where range ('-') can be used
  #get diff between indices
  del <- diff(f_ind2)
  #add max indx to vec
  del_idx <- c(which(del > 1), length(f_ind2))
  
  #start with first var
  fv <- paste0('-f', f_ind2[1])
  #if indices are a series, inset - max idx
  if (sum(del > 1) == 0)
  {
    fv <- paste0(fv, '-', tail(f_ind2, 1))
  } else {
    while (length(del) > 0)
    {
      #find first diff > 1
      if (sum(del > 1) != 0)
      {
        md <- min(which(del > 1))
        
        # add to call
        if (md > 1)
        {
          fv <- paste0(fv, '-', f_ind2[md], ',', f_ind2[md+1])
        } else {
          fv <- paste0(fv, ',', f_ind2[md+1])
        }
        
        # reorg del and f_ind2
        del <- del[-c(1:md)]
        f_ind2 <- f_ind2[-c(1:md)]
      } else {
        #if series to end
        fv <- paste0(fv, '-', tail(f_ind2, 1))
        del <- del[-c(1:length(del))]
      }
    } 
  }
  
  #function to subset draws in csv to df
  awk_fun <- function(file, IDX, thin)
  {
    print(paste0('subsetting ', file))
    #remove header
    p1_awk_call <- paste0("awk -F: '/^[^#]/ {print}' ", file, 
                          # select params of interest
                          " | cut -d \",\" ", IDX) 
    if (thin == 1)
    {
      awk_call <- paste0(p1_awk_call, 
                         " > ",
                         strsplit(file, '.csv')[[1]], '-draws-subset.csv')
    } else {
      #thin after keeping first line
      #https://unix.stackexchange.com/questions/648113/how-to-skip-every-three-lines-using-awk
      awk_call <- paste0(p1_awk_call, 
                         " | awk 'NR%", thin, "==1' > ", 
                         strsplit(file, '.csv')[[1]], '-draws-subset.csv')
    }
    system(awk_call)
    
    #get header only and write to file
    call2 <- paste0("head -n 47 ", file, " > ", 
                    strsplit(file, '.csv')[[1]], '-header.csv')
    system(call2)
    
    #get elapsed time
    call3 <- paste0("tail -n 5 ", file, " > ", 
                    strsplit(file, '.csv')[[1]], '-time.csv')
    system(call3)
    
    print(paste0('combining and writing to file'))
    #combine and write to file
    awk_call4 <- paste0('cat ', 
                        strsplit(file, '.csv')[[1]], '-header.csv ', 
                        strsplit(file, '.csv')[[1]], '-draws-subset.csv ',
                        strsplit(file, '.csv')[[1]], '-time.csv > ',
                        strsplit(file, '.csv')[[1]], '-subset-comb.csv')
    system(awk_call4)
  }
  
  #run awk fun and put all chains into list
  invisible(lapply(cmdstan_files, FUN = function(x) awk_fun(file = x, 
                                                            IDX = fv,
                                                            thin = thin)))
  
  #list files in dir
  lf <- list.files(dirname(cmdstan_files)[1], full.names = TRUE)
  #get names of comb files
  new_files <- grep('-subset-comb.csv', lf, value = TRUE)
  #make sure only files pertaining to this run
  tt <- strsplit(basename(cmdstan_files), '.csv')[[1]]
  new_files2 <- grep(substring(tt, nchar(tt)-5, nchar(tt)), new_files, value = TRUE)
  
  print(paste0('creating new cmdstanr object'))
  model_fit2 <- cmdstanr::as_cmdstan_fit(new_files2)
  
  return(model_fit2)
}

pp <- c('lp__', 'z', 'w', 'sigma',
        'mu_w', 'sigma_w', 'mm', 'sm')

est_fit2 <- subset_cmdstanr(model_fit = est_fit, 
                        target_params = pp,
                        thin = 2)


# save summary space ------------------------------------------------------------

#save out summary, model fit, data
MCMCvis::MCMCdiag(est_fit2,
                  round = 4,
                  file_name = paste0('diet-est-results-', run_date),
                  dir = paste0(dir, 'Results'),
                  mkdir = paste0('diet-est-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('diet-est-fit-', run_date),
                  add_obj = list(DATA_est),
                  add_obj_names = paste0('diet-est-data-', run_date),
                  cp_file = c(paste0(dir, 'Scripts/Model_files/4-diet-est.stan'),
                              paste0(dir, 'Scripts/4-diet-est.R')),
                  cp_file_names = c(paste0('4-diet-est-', run_date, '.stan'),
                                    paste0('4-diet-est-', run_date, '.R')))

# result_date <- '2025-05-29'
# est_fit2 <- readRDS(paste0(dir, '/Results/diet-est-', result_date,
#                       '/diet-est-fit-', result_date, '.rds'))
# dp2 <- readRDS(paste0(dir, '/Results/diet-est-', result_date,
#                       '/diet-est-data-', result_date, '.rds'))$pro_data
# library(shinystan)
# shinystan::launch_shinystan(est_fit2)


# dietary index -----------------------------------------------------------

#posterior ch for w
w_ch <- MCMCvis::MCMCchains(est_fit2, params = 'w')

#min-max transform
#x - min / max - min
#instead of min/max use 0.1 and 99.9 quantiles
LQ <- 0.001
UQ <- 0.999
# mm <- max(w_ch) - min(w_ch)
# comb_DI_ch <- (w_ch - min(w_ch)) / mm
qs <- quantile(w_ch, probs = c(LQ, UQ))
mm <- diff(qs)
comb_DI_ch <- (w_ch - qs[1]) / mm

# ~ 0.2 percent of posterior iterations are outside of 0, 100 range
(sum(comb_DI_ch > 1 | comb_DI_ch < 0) / length(comb_DI_ch)) * 100

d_sy$DI <- apply(comb_DI_ch, 2, mean) * 100
d_sy$sd_DI <- apply(comb_DI_ch, 2, sd) * 100
d_sy$w <- apply(w_ch, 2, mean)
d_sy$sd_w <- apply(w_ch, 2, sd)

# plot(d_sy$w, d_sy$DI)
# range(d_sy$DI)

saveRDS(d_sy, paste0(dir, 'Results/diet-est-', run_date, '/DI_df.rds'))


# Sat availability plot by site -----------------------------------------------

sat_images_tplt <- dplyr::group_by(dp2, site_id, season, doy) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(site_id, season) %>%
  #number images each site/year
  dplyr::count() %>%
  dplyr::ungroup() %>%
  #fill in missing years
  dplyr::right_join(tidyr::expand(dp2, site_id, season = min(dp2$season):max(dp2$season)), 
                    by = c('site_id', 'season')) %>%
  dplyr::arrange(site_id, season) %>%
  #fill in missing values with 0
  tidyr::replace_na(list(n = 0)) %>%
  #join with lat/lon
  dplyr::left_join(dplyr::distinct(dplyr::select(sites, site_id, longitude, latitude)), 
                   by = c('site_id')) %>%
  dplyr::arrange(site_id, season) %>%
  #negative lats, abs and add 180
  dplyr::mutate(lsrt = if_else(longitude < 0, longitude + 360, longitude)) %>%
  dplyr::arrange(lsrt, season)

#range vals
range(sat_images_tplt$n)

#reorder factors
sat_images_tplt$site_id <- factor(sat_images_tplt$site_id, 
                                  levels = unique(sat_images_tplt$site_id))

sat_plt <- ggplot(data = sat_images_tplt, 
                  aes(x = season, y = site_id)) + 
  geom_tile(aes(fill = cut(n, breaks = c(-1, 0, 2, 4, 5, 100),
                           labels = c('0', '1-2', '3-4', '4-5',
                                      '6-21'))),
            color = 'black', linewidth = 0.5) +
  scale_fill_manual(values = colorRampPalette(c("white","red"))(6), 
                    na.value = "darkgrey", name = "# satellite images") +
  theme_minimal()

ggsave(sat_plt, filename = paste0(dir, 'Results/diet-est-', run_date, '/sat_plt.pdf'),
       width = 5, height = 10)

# dplyr::left_join(dplyr::distinct(sat_images_tplt, site_id),
#                  dplyr::distinct(dp2, site_id, ccamlr_id, region), by = 'site_id') %>%
#   print(n = 500)
# dplyr::select(sat_images_tplt, site_id, latitude, longitude, lsrt) %>%
#   dplyr::distinct() %>%
#   dplyr::arrange(desc(lsrt)) %>%
#   print(n = 500)


# Sat availability plot ---------------------------------------------------

#colors 
#85d0bc
#f7bc3b
#ec8b26

#num scenes over time
sat_images <- dplyr::group_by(dp2, Landsat.Product.Identifier.L1) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  ggplot(aes(season, fill = factor(satellite))) +
  geom_bar() +
  scale_fill_manual(values=c("#85d0bc", "#f7bc3b", "#ec8b26")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

ggsave(sat_images, filename = paste0(dir, 'Results/diet-est-', run_date, '/sat_images.pdf'),
       width = 6, height = 5)
