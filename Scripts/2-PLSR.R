######################
# 2 - Convolve guano spectra into Landsat bands and run PLSR
# 
######################


# load packages -----------------------------------------------------------

library(tidyverse)
library(pls)
library(here)


# set/create dir -----------------------------------------------------------------

dir <- paste0(here::here(), '/')
run_date <- Sys.Date()
#create dir
dir.create(paste0(dir, 'Data/L2/', run_date))


# load data --------------------------------------------------------

#load SIA/spectra data
SIA_spectra <- read.csv(paste0(dir, 'Data/L1/SIA_spectra.csv'), 
                        header = TRUE, stringsAsFactors = FALSE)

#number of sampling locations
length(unique(SIA_spectra$Site))

waves <- 350:2500
wave_idx <- grep('Wave', colnames(SIA_spectra))
spec_mat <- as.matrix(SIA_spectra[, wave_idx])
SIA_data <- SIA_spectra[, -wave_idx]

#read in Relative Spectral Responses
#spectral response curves: https://landsat.usgs.gov/spectral-characteristics-viewer
#Landsat 4
L4_blue <- readxl::read_excel(paste0(dir, 
                                     'Data/L0/RSR/L4_TM_RSR.xlsx'), 
                              sheet = 'Blue-L4 TM')
L4_green <- readxl::read_excel(paste0(dir, 
                                      'Data/L0/RSR/L4_TM_RSR.xlsx'), 
                               sheet = 'Green-L4 TM')
L4_red <- readxl::read_excel(paste0(dir, 
                                    'Data/L0/RSR/L4_TM_RSR.xlsx'), 
                             sheet = 'Red-L4 TM')
L4_nir <- readxl::read_excel(paste0(dir, 
                                    'Data/L0/RSR/L4_TM_RSR.xlsx'), 
                             sheet = 'NIR-L4 TM')
L4_swir5 <- readxl::read_excel(paste0(dir, 
                                      'Data/L0/RSR/L4_TM_RSR.xlsx'), 
                               sheet = 'SWIR(5)-L4 TM')
L4_swir7 <- readxl::read_excel(paste0(dir, 
                                      'Data/L0/RSR/L4_TM_RSR.xlsx'), 
                               sheet = 'SWIR(7)-L4 TM')
#Landsat 5
L5_blue <- readxl::read_excel(paste0(dir, 
                                     'Data/L0/RSR/L5_TM_RSR.xlsx'), 
                              sheet = 'Blue-L5 TM')
L5_green <- readxl::read_excel(paste0(dir, 
                                      'Data/L0/RSR/L5_TM_RSR.xlsx'), 
                               sheet = 'Green-L5 TM')
L5_red <- readxl::read_excel(paste0(dir, 
                                    'Data/L0/RSR/L5_TM_RSR.xlsx'), 
                             sheet = 'Red-L5 TM')
L5_nir <- readxl::read_excel(paste0(dir, 
                                    'Data/L0/RSR/L5_TM_RSR.xlsx'), 
                             sheet = 'NIR-L5 TM')
L5_swir5 <- readxl::read_excel(paste0(dir, 
                                      'Data/L0/RSR/L5_TM_RSR.xlsx'), 
                               sheet = 'SWIR(5)-L5 TM')
L5_swir7 <- readxl::read_excel(paste0(dir, 
                                      'Data/L0/RSR/L5_TM_RSR.xlsx'), 
                               sheet = 'SWIR(7)-L5 TM')
#Landsat 7
L7_blue <- readxl::read_excel(paste0(dir, 
                                     'Data/L0/RSR/L7_ETM_RSR.xlsx'), 
                              sheet = 'Blue-L7')
L7_green <- readxl::read_excel(paste0(dir, 
                                      'Data/L0/RSR/L7_ETM_RSR.xlsx'), 
                               sheet = 'Green-L7')
L7_red <- readxl::read_excel(paste0(dir, 
                                    'Data/L0/RSR/L7_ETM_RSR.xlsx'), 
                             sheet = 'Red-L7')
L7_nir <- readxl::read_excel(paste0(dir, 
                                    'Data/L0/RSR/L7_ETM_RSR.xlsx'), 
                             sheet = 'NIR-L7')
L7_swir5 <- readxl::read_excel(paste0(dir, 
                                      'Data/L0/RSR/L7_ETM_RSR.xlsx'),
                               sheet = 'SWIR(5)-L7')
L7_swir7 <- readxl::read_excel(paste0(dir, 
                                      'Data/L0/RSR/L7_ETM_RSR.xlsx'),
                               sheet = 'SWIR(7)-L7')


# process RSR -------------------------------------------------------------

#Landsat 4
L4_rsr <- dplyr::left_join(data.frame(Wavelength = 350:2500), 
                           L4_blue, by = 'Wavelength') %>%
  dplyr::left_join(L4_green, by = 'Wavelength') %>%
  dplyr::left_join(L4_red, by = 'Wavelength') %>%
  dplyr::left_join(L4_nir, by = 'Wavelength') %>%
  dplyr::left_join(L4_swir5, by = 'Wavelength') %>%
  dplyr::left_join(L4_swir7, by = 'Wavelength') %>%
  replace(is.na(.), 0)
colnames(L4_rsr) <- c('Wavelength', 'Band_1', 'Band_2', 'Band_3',
                      'Band_4', 'Band_5', 'Band_7')
#Landsat 5
L5_rsr <- dplyr::left_join(data.frame(Wavelength = 350:2500), 
                           L5_blue, by = 'Wavelength') %>%
  dplyr::left_join(L5_green, by = 'Wavelength') %>%
  dplyr::left_join(L5_red, by = 'Wavelength') %>%
  dplyr::left_join(L5_nir, by = 'Wavelength') %>%
  dplyr::left_join(L5_swir5, by = 'Wavelength') %>%
  dplyr::left_join(L5_swir7, by = 'Wavelength') %>%
  replace(is.na(.), 0)
colnames(L5_rsr) <- c('Wavelength', 'Band_1', 'Band_2', 'Band_3',
                      'Band_4', 'Band_5', 'Band_7')
#Landsat 7
L7_rsr <- dplyr::left_join(data.frame(Wavelength = 350:2500), 
                           L7_blue, by = 'Wavelength') %>%
  dplyr::left_join(L7_green, by = 'Wavelength') %>%
  dplyr::left_join(L7_red, by = 'Wavelength') %>%
  dplyr::left_join(L7_nir, by = 'Wavelength') %>%
  dplyr::left_join(L7_swir5, by = 'Wavelength') %>%
  dplyr::left_join(L7_swir7, by = 'Wavelength') %>%
  replace(is.na(.), 0)
colnames(L7_rsr) <- c('Wavelength', 'Band_1', 'Band_2', 'Band_3',
                      'Band_4', 'Band_5', 'Band_7')

# #plot RSR
# tidyr::pivot_longer(L4_rsr,
#                     cols = starts_with('Band'),
#                     names_to = 'Band',
#                     values_to = 'Value') %>%
#   ggplot(aes(Wavelength, Value, color = Band)) +
#   geom_line(linewidth = 1.1) +
#   xlim(c(350, 2500)) +
#   theme_bw()


# plot guano spectra ------------------------------------------------------

#get mean and sd reflectance across all samples
spec_mn <- apply(spec_mat, 2, mean)
spec_sd <- apply(spec_mat, 2, sd)

#wavelengths
wave_vals <- as.numeric(sapply(strsplit(colnames(spec_mat), '_'), `[`, 2))

#df to plt
tplt <- data.frame(wave = wave_vals,
                   spec_mn = spec_mn,
                   spec_sd = spec_sd)

gp <- ggplot(data = tplt) +
  # geom_line(aes(wave, spec_mn - spec_sd),
  #           linetype = 'dashed',
  #           color = 'grey',
  #           linewidth = 1.1) +
  # geom_line(aes(wave, spec_mn + spec_sd),
  #           linetype = 'dashed',
  #           color = 'grey',
  #           linewidth = 1.1) +
  geom_ribbon(aes(wave, 
                  ymin = spec_mn - spec_sd,
              ymax = spec_mn + spec_sd),
              alpha = 0.4) +
  geom_line(aes(wave, spec_mn),
            linewidth = 1.1) +
  # xlim(c(350, 2500)) +
  ylab('Reflectance') +
  xlab('Wavelength') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

ggsave(gp, filename = paste0(dir, 'Data/L2/', run_date, '/guano_spectra.pdf'),
       width = 6, height = 5)


#spectra with L7 bands
rng_b1 <- range(dplyr::filter(L7_rsr, Band_1 > 0)$Wavelength)
rng_b2 <- range(dplyr::filter(L7_rsr, Band_2 > 0)$Wavelength)
rng_b3 <- range(dplyr::filter(L7_rsr, Band_3 > 0)$Wavelength)
rng_b4 <- range(dplyr::filter(L7_rsr, Band_4 > 0)$Wavelength)
rng_b5 <- range(dplyr::filter(L7_rsr, Band_5 > 0)$Wavelength)
rng_b7 <- range(dplyr::filter(L7_rsr, Band_7 > 0)$Wavelength)

ymin <- 0
ymax <- 0.6
alpha <- 0.01

gp2 <- ggplot(data = tplt) +
  geom_rect(xmin = rng_b1[1],
            xmax = rng_b1[2],
            ymin = ymin,
            ymax = ymax,
            fill = 'lightblue',
            color = 'black',
            alpha = alpha) +
  geom_rect(xmin = rng_b2[1],
            xmax = rng_b2[2],
            ymin = ymin,
            ymax = ymax,
            fill = 'lightblue',
            color = 'black',
            alpha = alpha) +
  geom_rect(xmin = rng_b3[1],
            xmax = rng_b3[2],
            ymin = ymin,
            ymax = ymax,
            fill = 'lightblue',
            color = 'black',
            alpha = alpha) +
  geom_rect(xmin = rng_b4[1],
            xmax = rng_b4[2],
            ymin = ymin,
            ymax = ymax,
            fill = 'lightblue',
            color = 'black',
            alpha = alpha) +
  geom_rect(xmin = rng_b5[1],
            xmax = rng_b5[2],
            ymin = ymin,
            ymax = ymax,
            fill = 'lightblue',
            color = 'black',
            alpha = alpha) +
  geom_rect(xmin = rng_b7[1],
            xmax = rng_b7[2],
            ymin = ymin,
            ymax = ymax,
            fill = 'lightblue',
            color = 'black',
            alpha = alpha) +
  geom_ribbon(aes(wave, 
                  ymin = spec_mn - spec_sd,
                  ymax = spec_mn + spec_sd),
              alpha = 0.4) +
  geom_line(aes(wave, spec_mn),
            linewidth = 1.1) +
  # xlim(c(350, 2500)) +
  ylab('Reflectance') +
  xlab('Wavelength') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

ggsave(gp2, filename = paste0(dir, 'Data/L2/', run_date, '/guano_spectra_L7_bands.pdf'),
       width = 6, height = 5)


# convolve spectra --------------------------------------------------------

spec_mat_n <- spec_mat

#multiply target spectra by RSR for band
#sum over scaled values and divide by sum of RSR for band
#target is vector of normalized reflectance values
c_fun <- function(target, RSR)
{
  return(colSums(target * RSR[,-1]) / colSums(RSR[,-1]))
}

#apply fun to all samples - samples in rows, bands in cols
L4_spectra <- t(apply(spec_mat_n, 1, 
                    FUN = function(x) c_fun(target = x, RSR = L4_rsr)))
L5_spectra <- t(apply(spec_mat_n, 1, 
                    FUN = function(x) c_fun(target = x, RSR = L5_rsr)))
L7_spectra <- t(apply(spec_mat_n, 1, 
                    FUN = function(x) c_fun(target = x, RSR = L7_rsr)))

#merge convolved spectra with sample info
L4_sp_data <- cbind(SIA_data, L4_spectra)
L5_sp_data <- cbind(SIA_data, L5_spectra)
L7_sp_data <- cbind(SIA_data, L7_spectra)


# PLSR settings --------------------------------------------------------------

pls.options(plsralg = "oscorespls")


# VIP fun ---------------------------------------------------------------------

VIP <- function(object)
{
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, '*')
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}


# Split into train/val and run PLSR iteratively -------------------------------

plsr_fun <- function(data, PROP_CV = 0.65, ITER = 250)
{
  # CV model validation
  resamples <- ITER
  output_CV_stats <- data.frame(Sample = rep(NA, resamples * 7),
                                Comps = NA,
                                R2_cal = NA, 
                                R2_val = NA, 
                                RMSEP_cal = NA,
                                RMSEP_val = NA,
                                PERC_RMSEP_cal = NA,
                                PERC_RMSEP_val = NA,
                                Bias_val = NA,
                                SD_resid_val = NA)
  
  output_CV_coefs <- data.frame(Sample = rep(NA, resamples * 6),
                                Comps = NA,
                                Int = NA,
                                B1 = NA,
                                B2 = NA,
                                B3 = NA,
                                B4 = NA,
                                B5 = NA,
                                B7 = NA)
  
  output_CV_vips <- data.frame(Sample = rep(NA, resamples * 6),
                               Comps = NA,
                               B1 = NA,
                               B2 = NA,
                               B3 = NA,
                               B4 = NA,
                               B5 = NA,
                               B7 = NA)
  
  set.seed(1)
  #proportion split (train)
  prop_cv <- PROP_CV
  coef_vip_counter <- 1
  stats_counter <- 1
  data$nn <- 1:NROW(data)
  for (i in 1:resamples)
  {
    #i <- 1
    print(paste0('iteration: ', i, ' of ', resamples))
    
    #stratified sampling
    #divide into 5 bins based on quantiles
    bins <- quantile(data$d_15N, probs = c(0, 0.2, 0.4, 0.6, 0.8))
    #assign bin #
    data$bin <- findInterval(data$d_15N, bins)
    #draw train/test from each bin
    train_rows <- rep(NA, floor(prop_cv * nrow(data)))
    counter <- 1
    for (j in 1:length(bins))
    {
      #j <- 1
      #nn is row number in data
      fd <- dplyr::filter(data, bin == j)$nn
      sr <- sample(fd, size = floor(prop_cv * length(fd)))
      train_rows[counter:(counter + length(sr) - 1)] <- sr
      
      counter <- counter + length(sr)
    }
    train_rows2 <- train_rows[!is.na(train_rows)]
    
    #create train/test
    cal_data_cv <- data[train_rows2,]
    val_data_cv <- data[-train_rows2,]
    dimsCal <- dim(cal_data_cv)
    dimsVal <- dim(val_data_cv)
    
    ### Build PLSR model with training data
    pls_cv <- pls::plsr(d_15N ~ Spectra, 
                        scale = FALSE,
                        ncomp = 6, 
                        validation = "none", 
                        data = cal_data_cv, 
                        trace = FALSE)
    
    # predict training data
    pred_cal_data <- predict(pls_cv,
                             newdata = cal_data_cv$Spectra,
                             ncomp = 1:6,
                             type = "response")[,1,]
    
    #predict val data
    pred_val_data <- predict(pls_cv, 
                             newdata = val_data_cv$Spectra,
                             ncomp = 1:6, 
                             type = "response")[,1,]
    
    #df indices to fill
    coef_vip_idx <- coef_vip_counter:(coef_vip_counter + 5)
    stats_idx <- stats_counter:(stats_counter + 6)
    
    #fill coefs df
    coef <- coef(pls_cv, 
                 ncomp = 1:6, 
                 intercept = TRUE)[,1,]
    
    output_CV_coefs$Sample[coef_vip_idx] <- i
    output_CV_coefs$Comps[coef_vip_idx] <- 1:6
    output_CV_coefs$Int[coef_vip_idx] <- coef[1,]
    output_CV_coefs$B1[coef_vip_idx] <- coef[2,]
    output_CV_coefs$B2[coef_vip_idx] <- coef[3,]
    output_CV_coefs$B3[coef_vip_idx] <- coef[4,]
    output_CV_coefs$B4[coef_vip_idx] <- coef[5,]
    output_CV_coefs$B5[coef_vip_idx] <- coef[6,]
    output_CV_coefs$B7[coef_vip_idx] <- coef[7,]
    
    #fill VIPS df
    vips <- VIP(pls_cv)
    
    output_CV_vips$Sample[coef_vip_idx] <- i
    output_CV_vips$Comps[coef_vip_idx] <- 1:6
    output_CV_vips$B1[coef_vip_idx] <- vips[1,]
    output_CV_vips$B2[coef_vip_idx] <- vips[2,]
    output_CV_vips$B3[coef_vip_idx] <- vips[3,]
    output_CV_vips$B4[coef_vip_idx] <- vips[4,]
    output_CV_vips$B5[coef_vip_idx] <- vips[5,]
    output_CV_vips$B7[coef_vip_idx] <- vips[6,]
    
    coef_vip_counter <- coef_vip_counter + 6
    
    #fill stats df
    cal_residuals <- cal_data_cv$d_15N - pred_cal_data
    val_residuals <- val_data_cv$d_15N - pred_val_data
    
    output_CV_stats$Sample[stats_idx] <- i
    output_CV_stats$Comps[stats_idx] <- 0:6
    output_CV_stats$R2_cal[stats_idx] <- pls::R2(pls_cv, newdata = cal_data_cv)$val[,,(0:6) + 1]
    output_CV_stats$R2_val[stats_idx] <- pls::R2(pls_cv, newdata = val_data_cv)$val[,,(0:6) + 1]
    RMSEP_cal <- pls::RMSEP(pls_cv, newdata = cal_data_cv)$val[,,(0:6) + 1]
    RMSEP_val <- pls::RMSEP(pls_cv, newdata = val_data_cv)$val[,,(0:6) + 1]
    output_CV_stats$RMSEP_cal[stats_idx] <- RMSEP_cal
    output_CV_stats$RMSEP_val[stats_idx] <- RMSEP_val
    output_CV_stats$PERC_RMSEP_cal[stats_idx] <- (RMSEP_cal/(max(cal_data_cv$d_15N) - 
                                                               min(cal_data_cv$d_15N)))*100
    output_CV_stats$PERC_RMSEP_val[stats_idx] <- (RMSEP_val/(max(val_data_cv$d_15N) - 
                                                               min(val_data_cv$d_15N)))*100
    #NA for intercept
    output_CV_stats$Bias_val[stats_idx] <- c(NA, apply(pred_val_data, 2, mean) - 
                                               mean(val_data_cv$d_15N))
    output_CV_stats$SD_resid_val[stats_idx] <- c(NA, apply(val_residuals, 2, sd))
    #advance counter
    stats_counter <- stats_counter + 7
    
    # Remove temp objects - FIX
    rm(cal_data_cv, val_data_cv,
       pls_cv, 
       pred_cal_data, pred_val_data,
       coef_vip_idx, stats_idx,
       coef, vips,
       cal_residuals, val_residuals,
       RMSEP_cal, RMSEP_val)
  }
  
  return(list(output_CV_stats = output_CV_stats,
              output_CV_coefs = output_CV_coefs,
              output_CV_vips = output_CV_vips))
}


# make RMSEP, R2, Coef, VIP plots -----------------------------------------

plt_fun <- function(plsr_run, plot_dir, name)
{
  #RMSEP
  dd <- dplyr::group_by(plsr_run$output_CV_stats, Comps) %>%
    dplyr::summarize(mn_RMSEP_cal = mean(RMSEP_cal),
                     sd_RMSEP_cal = sd(RMSEP_cal),
                     mn_RMSEP_val = mean(RMSEP_val),
                     sd_RMSEP_val = sd(RMSEP_val),
                     mn_R2_cal = mean(R2_cal),
                     sd_R2_cal = sd(R2_cal),
                     mn_R2_val = mean(R2_val),
                     sd_R2_val = sd(R2_val))
  
  #coefs
  dd_coefs <- dplyr::group_by(plsr_run$output_CV_coefs, Comps) %>%
    dplyr::summarize(mn_Int = mean(Int),
                     sd_Int = sd(Int),
                     mn_B1 = mean(B1),
                     sd_B1 = sd(B1),
                     mn_B2 = mean(B2),
                     sd_B2 = sd(B2),
                     mn_B3 = mean(B3),
                     sd_B3 = sd(B3),
                     mn_B4 = mean(B4),
                     sd_B4 = sd(B4),
                     mn_B5 = mean(B5),
                     sd_B5 = sd(B5),
                     mn_B7 = mean(B7),
                     sd_B7 = sd(B7)) %>%
    dplyr::filter(Comps == 5)
  
  dd_coefs_mn <- data.frame(pos = 1:6,
                            coefs = t(dplyr::select(dd_coefs, 
                                                    tidyselect::starts_with('mn'),
                                                    -mn_Int)))
  dd_coefs_ci <- data.frame(pos = 1:6, 
                            lci = dd_coefs_mn$coefs - 
                              t(dplyr::select(dd_coefs, 
                                              tidyselect::starts_with('sd'),
                                              -sd_Int))[,1],
                            uci = dd_coefs_mn$coefs + 
                              t(dplyr::select(dd_coefs, 
                                              tidyselect::starts_with('sd'),
                                              -sd_Int))[,1])
  
  #VIP
  dd_vips <- dplyr::group_by(plsr_run$output_CV_vips, Comps) %>%
    dplyr::summarize(mn_B1 = mean(B1),
                     sd_B1 = sd(B1),
                     mn_B2 = mean(B2),
                     sd_B2 = sd(B2),
                     mn_B3 = mean(B3),
                     sd_B3 = sd(B3),
                     mn_B4 = mean(B4),
                     sd_B4 = sd(B4),
                     mn_B5 = mean(B5),
                     sd_B5 = sd(B5),
                     mn_B7 = mean(B7),
                     sd_B7 = sd(B7)) %>%
    dplyr::filter(Comps == 5)
  
  dd_vips_mn <- data.frame(pos = 1:6,
                           vips = t(dplyr::select(dd_vips, 
                                                  tidyselect::starts_with('mn'))))
  dd_vips_ci <- data.frame(pos = 1:6, 
                           lci = dd_vips_mn$vips - 
                             t(dplyr::select(dd_vips, 
                                             tidyselect::starts_with('sd')))[,1],
                           uci = dd_vips_mn$vips + 
                             t(dplyr::select(dd_vips, 
                                             tidyselect::starts_with('sd')))[,1])
  
  #plot RMSEP
  rmsep_plt <- ggplot() +
    geom_line(data = dd, 
              aes(Comps, mn_RMSEP_cal), 
              color = 'black',
              linewidth = 1.2) +
    geom_line(data = dd, 
              aes(Comps, mn_RMSEP_val), 
              color = 'red',
              linewidth = 1.2,
              linetype = 'dashed') +
    ggtitle(name) +
    xlab('# Comps') +
    ylab('Model RMSEP') +
    scale_x_continuous(breaks = 0:6) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(color = 'black'))
  
  #plot R2
  r2_plt <- ggplot() +
    geom_line(data = dd, 
              aes(Comps, mn_R2_cal), 
              color = 'black',
              linewidth = 1.2) +
    geom_line(data = dd, 
              aes(Comps, mn_R2_val), 
              color = 'red',
              linewidth = 1.2,
              linetype = 'dashed') +
    ggtitle(name) +
    xlab('# Comps') +
    ylab('Model R2') +
    scale_y_continuous(breaks = seq(0.0, 0.6, by = 0.1)) +
    scale_x_continuous(breaks = 0:6) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(color = 'black'))
  
  #plot coefs
  coef_plt <- ggplot() +
    geom_hline(yintercept = 0, 
               linetype = 'dashed',
               linewidth = 1.1,
               alpha = 0.5) +
    geom_ribbon(data = dd_coefs_ci,
                aes(x = pos,
                    ymin = lci,
                    ymax = uci),
                alpha = 0.2) +
    geom_line(data = dd_coefs_mn, 
              aes(pos, coefs), 
              color = 'black',
              linewidth = 1.2) +
    ggtitle(name) +
    xlab('Band #') +
    ylab('Model Coefs') +
    scale_x_continuous(breaks = 1:6) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = 'black'))
  
  #plot VIP
  vip_plt <- ggplot() +
    geom_line(data = dd_vips_mn, 
              aes(pos, vips), 
              color = 'black',
              linewidth = 1.2) +
    ggtitle(name) +
    xlab('Band #') +
    ylab('Model VIP') +
    scale_x_continuous(breaks = 1:6) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = 'black'))
  
  #save out plots
  ggsave(filename = paste0(plot_dir, name, '_RMSEP_plot.pdf'),
                           rmsep_plt)
  ggsave(filename = paste0(plot_dir, name, '_R2_plot.pdf'),
                           r2_plt)
  ggsave(filename = paste0(plot_dir, name, '_Coef_plot.pdf'),
                           coef_plt)
  ggsave(filename = paste0(plot_dir, name, '_VIP_plot.pdf'),
                           vip_plt)
}


# ensemble coefs ----------------------------------------------------------

#make ensemble predication plot, return ensemble coefs
en_fun <- function(plsr_run, plsr_data, COMPS = 5, name, plot_dir)
{
  #ensemble coefs
  LX_cv_coefs_cX <- dplyr::filter(plsr_run$output_CV_coefs, Comps == COMPS)%>%
    dplyr::select(-Sample, - Comps)
  LX_cv_coefs_cX_mn <- dplyr::summarize_all(LX_cv_coefs_cX, mean)
  LX_en_pred <- cbind(1, plsr_data$Spectra) %*% as.numeric(LX_cv_coefs_cX_mn)
  
  rng_vals <- range(c(LX_en_pred, plsr_data$d_15N))
  # plot(LX_en_pred, 
  #      LX_plsr_data$d_15N,
  #      xlab = 'Predicted d_15N',
  #      ylab = 'Observed d_15N',
  #      pch = 19,
  #      col = 'grey',
  #      main = name,
  #      xlim = rng_vals,
  #      ylim = rng_vals)
  # abline(0, 1, lty = 2, lwd = 2)
  
  #iterating over each model
  LX_en_pred_ch <- matrix(NA,
                          nrow = NROW(LX_cv_coefs_cX),
                          ncol = NROW(plsr_data))
  for (i in 1:NROW(LX_cv_coefs_cX))
  {
    #i <- 1
    LX_en_pred_ch[i,] <- cbind(1, plsr_data$Spectra) %*%
      as.numeric(LX_cv_coefs_cX[i,])
  }
  
  LX_en_pred_unc <- data.frame(obs = plsr_data$d_15N,
                               pred_mn = apply(LX_en_pred_ch, 2, mean),
                               pred_sd = apply(LX_en_pred_ch, 2, sd)) %>%
    dplyr::mutate(pred_LCI = pred_mn - pred_sd,
                  pred_UCI = pred_mn + pred_sd)
  
  pred_plt <- ggplot(data = LX_en_pred_unc) +
    geom_point(aes(pred_mn, obs),
               alpha = 0.6,
               size = 2) +
    geom_errorbarh(aes(xmin = pred_LCI,
                       xmax = pred_UCI,
                       y = obs),
                   alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1,
                linetype = 'dashed',
                alpha = 0.6,
                linewidth = 1.1,
                col = 'red') +
    xlim(rng_vals) +
    ylim(rng_vals) +
    xlab('Predicted d_15N') +
    ylab('Observed d_15N') +
    ggtitle(name) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = 'black'))
  
  ggsave(filename = paste0(plot_dir, name, '_Pred_plot.pdf'),
         pred_plt)
  
  return(LX_cv_coefs_cX)
}


# Landsat 4 ---------------------------------------------------------------

#format data
L4_plsr_data <- data.frame(L4_sp_data) %>%
  dplyr::mutate(Spectra = I(as.matrix(L4_sp_data[,grep('Band', colnames(L4_sp_data))]))) %>%
  dplyr::select(-starts_with('Band'))

#run PLSR
L4_run <- plsr_fun(L4_plsr_data, PROP_CV = 0.65, ITER = 1000)

#test for # comps - 5 is best
t.test(dplyr::filter(L4_run$output_CV_stats, Comps == 5)$RMSEP_val,
       dplyr::filter(L4_run$output_CV_stats, Comps == 6)$RMSEP_val)
t.test(dplyr::filter(L4_run$output_CV_stats, Comps == 5)$R2_val,
       dplyr::filter(L4_run$output_CV_stats, Comps == 6)$R2_val)

#plots
plt_fun(plsr_run = L4_run, 
        plot_dir = paste0(dir, 'Data/L2/', run_date, '/'),
        name = 'L4')

#L4 summary
#r2
L4_cv_stats_c5 <- dplyr::filter(L4_run$output_CV_stats, Comps == 5)
hist(L4_cv_stats_c5$R2_val)
median(L4_cv_stats_c5$R2_val)
IQR(L4_cv_stats_c5$R2_val)

#bias
hist(L4_cv_stats_c5$Bias_val)

#RMSEP
hist(L4_cv_stats_c5$RMSEP_val)

#sd resid
hist(L4_cv_stats_c5$SD_resid_val)
median(L4_cv_stats_c5$SD_resid_val)

#ensemble pred/obs plot and ensemble predictions
L4_coef_en <- en_fun(plsr_run = L4_run, 
                     plsr_data = L4_plsr_data, 
                     COMPS = 5, 
                     name = 'L4', 
                     plot_dir = paste0(dir, 'Data/L2/', run_date, '/'))
                     

# Landsat 5 ---------------------------------------------------------------

#format data
L5_plsr_data <- data.frame(L5_sp_data) %>%
  dplyr::mutate(Spectra = I(as.matrix(L5_sp_data[,grep('Band', colnames(L5_sp_data))]))) %>%
  dplyr::select(-starts_with('Band'))

#run PLSR
L5_run <- plsr_fun(L5_plsr_data, PROP_CV = 0.65, ITER = 1000)

#test for # comps - 5 is best
t.test(dplyr::filter(L5_run$output_CV_stats, Comps == 5)$RMSEP_val,
       dplyr::filter(L5_run$output_CV_stats, Comps == 6)$RMSEP_val)
t.test(dplyr::filter(L5_run$output_CV_stats, Comps == 5)$R2_val,
       dplyr::filter(L5_run$output_CV_stats, Comps == 6)$R2_val)

#plots
plt_fun(plsr_run = L5_run, 
        plot_dir = paste0(dir, 'Data/L2/', run_date, '/'),
        name = 'L5')

#L5 summary
#r2
L5_cv_stats_c5 <- dplyr::filter(L5_run$output_CV_stats, Comps == 5)
hist(L5_cv_stats_c5$R2_val)
median(L5_cv_stats_c5$R2_val)
IQR(L5_cv_stats_c5$R2_val)

#bias
hist(L5_cv_stats_c5$Bias_val)

#RMSEP
hist(L5_cv_stats_c5$RMSEP_val)

#sd resid
hist(L5_cv_stats_c5$SD_resid_val)
median(L5_cv_stats_c5$SD_resid_val)

#ensemble pred/obs plot and ensemble predictions
L5_coef_en <- en_fun(plsr_run = L5_run, 
                     plsr_data = L5_plsr_data, 
                     COMPS = 5, 
                     name = 'L5', 
                     plot_dir = paste0(dir, 'Data/L2/', run_date, '/'))


# Landsat 7 ---------------------------------------------------------------

#format data
L7_plsr_data <- data.frame(L7_sp_data) %>%
  dplyr::mutate(Spectra = I(as.matrix(L7_sp_data[,grep('Band', colnames(L7_sp_data))]))) %>%
  dplyr::select(-starts_with('Band'))

#run PLSR
L7_run <- plsr_fun(L7_plsr_data, PROP_CV = 0.65, ITER = 1000)

#test for # comps - 5 is best
t.test(dplyr::filter(L7_run$output_CV_stats, Comps == 5)$RMSEP_val,
       dplyr::filter(L7_run$output_CV_stats, Comps == 6)$RMSEP_val)
t.test(dplyr::filter(L7_run$output_CV_stats, Comps == 5)$R2_val,
       dplyr::filter(L7_run$output_CV_stats, Comps == 6)$R2_val)

#plots
plt_fun(plsr_run = L7_run, 
        plot_dir = paste0(dir, 'Data/L2/', run_date, '/'),
        name = 'L7')

#L7 summary
#r2
L7_cv_stats_c5 <- dplyr::filter(L7_run$output_CV_stats, Comps == 5)
hist(L7_cv_stats_c5$R2_val)
median(L7_cv_stats_c5$R2_val)
IQR(L7_cv_stats_c5$R2_val)

#bias
hist(L7_cv_stats_c5$Bias_val)

#RMSEP
hist(L7_cv_stats_c5$RMSEP_val)

#sd resid
hist(L7_cv_stats_c5$SD_resid_val)
median(L7_cv_stats_c5$SD_resid_val)

#ensemble pred/obs plot and ensemble predictions
L7_coef_en <- en_fun(plsr_run = L7_run, 
                     plsr_data = L7_plsr_data, 
                     COMPS = 5, 
                     name = 'L7', 
                     plot_dir = paste0(dir, 'Data/L2/', run_date, '/'))
                     

# model summaries ---------------------------------------------------------

#sd of residuals
median(L4_cv_stats_c5$SD_resid_val)
median(L5_cv_stats_c5$SD_resid_val)
median(L7_cv_stats_c5$SD_resid_val)


# save out model results --------------------------------------------------

saveRDS(L4_coef_en, paste0(dir, 'Data/L2/', run_date, '/L4_coef_en.rds'))
saveRDS(L5_coef_en, paste0(dir, 'Data/L2/', run_date, '/L5_coef_en.rds'))
saveRDS(L7_coef_en, paste0(dir, 'Data/L2/', run_date, '/L7_coef_en.rds'))


# pred plot just L7 -------------------------------------------------------

#LX_cv_coefs_cX
L7_coef_en

# plsr_run = L7_run, 
# plsr_data = L7_plsr_data, 
# COMPS = 5, 
# name = 'L7', 
# plot_dir = paste0(dir, 'Data/L2/', run_date, '/')

L7_coefs_en_mn <- dplyr::summarize_all(L7_coef_en, mean)
L7_en_pred <- cbind(1, L7_plsr_data$Spectra) %*% as.numeric(L7_coefs_en_mn)
rng_vals <- range(c(L7_en_pred, L7_plsr_data$d_15N))

#iterating over each model
L7_en_pred_ch <- matrix(NA,
                        nrow = NROW(L7_coef_en),
                        ncol = NROW(L7_plsr_data))
for (i in 1:NROW(L7_coef_en))
{
  #i <- 1
  L7_en_pred_ch[i,] <- cbind(1, L7_plsr_data$Spectra) %*%
    as.numeric(L7_coef_en[i,])
}

L7_en_pred_unc <- data.frame(obs = L7_plsr_data$d_15N,
                             pred_mn = apply(L7_en_pred_ch, 2, mean),
                             pred_sd = apply(L7_en_pred_ch, 2, sd)) %>%
  dplyr::mutate(pred_LCI = pred_mn - pred_sd,
                pred_UCI = pred_mn + pred_sd)

L7_pred_plt <- ggplot(data = L7_en_pred_unc) +
  geom_abline(intercept = 0, slope = 1,
              linetype = 'dashed',
              alpha = 0.7,
              linewidth = 1.1,
              color = 'red') +
  geom_errorbarh(aes(xmin = pred_LCI,
                     xmax = pred_UCI,
                     y = obs),
                 alpha = 0.6) +
  geom_point(aes(pred_mn, obs),
             alpha = 0.7,
             size = 2) +
  xlim(rng_vals) +
  ylim(rng_vals) +
  xlab('Predicted d_15N') +
  ylab('Observed d_15N') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

ggsave(L7_pred_plt, filename = paste0(dir, 'Data/L2/', run_date, '/L7_pred_pub.pdf'),
       width = 6, height = 5)
