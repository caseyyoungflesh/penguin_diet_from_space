######################
# 5 - model diet ~ cov
#
# diet ~ DOY (within years)
# diet ~ SIC (across years)
# diet ~ bathymetry and SIC (across space)
######################


# Load packages -----------------------------------------------------------

library(MCMCvis)
library(tidyverse)
library(cmdstanr)
library(here)


# set dir -----------------------------------------------------------------

dir <- paste0(here::here(), '/')
run_date <- Sys.Date()
# run_date <- '2024-11-21'
# diet_est_run_date <- '2024-12-06'
diet_est_run_date <- '2025-07-16'


# read in data -------------------------------------------------------------

#output from 4-diet-est.R
dd <- readRDS(paste0(dir, 'Results/diet-est-', diet_est_run_date, '/DI_df.rds'))


# data stats --------------------------------------------------------------

#N captures each site
range(dplyr::count(dd, site_id)$n)
#N captures each site/year
range(dplyr::count(dd, site_id, season)$n)
#N years each site
dplyr::count(dd, site_id, season) %>%
  dplyr::count(site_id) %>%
  print(n = 100)


# SIC ---------------------------------------------------------------------

sic <- read.csv(paste0(dir, 'Data/L1/colony_sic.csv')) %>%
  dplyr::select(site_id, year, month, 
                site_easting = easting, 
                site_northing = northing, 
                longitude,
                latitude,
                sic_150km) %>%
  dplyr::filter(!is.na(site_id)) %>%
  dplyr::filter(site_id %in% unique(dd$site_id))

#SIC available for each site, year, month
#aggregate to site
#avg SIC Nov - Feb
#avg across years

#create season in sic
sic2 <- dplyr::filter(sic, 
                      month %in% c(11, 12, 1, 2))

#jan/feb idx
jf_idx <- which(sic2$month %in% c(1, 2))
#nov/dec idx
nd_idx <- which(sic2$month %in% c(11, 12))
sic2$season <- NA
sic2$season[jf_idx] <- sic2$year[jf_idx] - 1
sic2$season[nd_idx] <- sic2$year[nd_idx]

#sic each site/year
sic2b <- dplyr::group_by(sic2, site_id, season) %>%
  #average across months for each year (season)
  dplyr::summarize(sic_year = mean(sic_150km, na.rm = TRUE)) %>%
  dplyr::ungroup()

sic3 <- dplyr::left_join(sic2b, dplyr::distinct(dplyr::select(dd, site_id, 
                                                 season, s_idx, sy_idx)), 
                   by = c('site_id', 'season')) %>%
  dplyr::filter(!is.na(sy_idx)) %>%
  dplyr::group_by(site_id) %>%
  dplyr::mutate(sic_year_sc = scale(sic_year, scale = FALSE)[,1],
                mn_sic = mean(sic_year)) %>%
  dplyr::ungroup()

#sic each site
sic4 <- dplyr::group_by(sic2, site_id, season) %>%
  #average across months for each year (season)
  dplyr::summarize(sic_year = mean(sic_150km, na.rm = TRUE)) %>%
  dplyr::ungroup()  %>%
  dplyr::group_by(site_id) %>%
  dplyr::summarize(sic_site = mean(sic_year, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::distinct(dplyr::select(dd, site_id, s_idx)), 
                   by = 'site_id')

#add coord to df
site_coords <- dplyr::select(sic, site_id, 
                             site_easting, site_northing,
                             longitude, latitude) %>%
  dplyr::distinct()


# bath --------------------------------------------------------------------

bath <- read.csv(paste0(dir, 'Data/L1/shelf_area.csv')) %>%
  dplyr::select(site_id, shelf_area, frac_shelf) %>%
  dplyr::filter(!is.na(site_id)) %>%
  dplyr::filter(site_id %in% unique(dd$site_id))

bath2 <- dplyr::left_join(bath,
                          dplyr::distinct(dplyr::select(dd, site_id, s_idx)), 
                          by = 'site_id') %>%
  dplyr::mutate(shelf_area_sc = scale(shelf_area, scale = TRUE)[,1])


# DOY processing ------------------------------------------------------

#DOY
#120 day span
range(dd$doy)
#center DOY on day 60 (Dec 30)
dd$doy_centered <- dd$doy - 60

#make polynomial
poly_doy <- poly(dd$doy_centered, 2, raw = FALSE)


# join with site coords ---------------------------------------------------

dd2 <- dplyr::left_join(dd, site_coords, by = 'site_id')

#s index for each site, year
d_s <- dplyr::select(dd2, site_id, ccamlr_id, region,
                     season, avg_n_pix,
                     sy_idx, s_idx) %>%
  dplyr::distinct() %>%
  dplyr::arrange(sy_idx) %>%
  dplyr::left_join(site_coords, by = 'site_id')

#number sites
length(unique(dd2$site_id))

#number doy/site/years
NROW(dd2)


# fit model ---------------------------------------------------------------

DATA_model <- list(Nsyd = length(unique(dd2$syd_idx)),
                 Nsy = length(unique(dd2$sy_idx)),
                 Ns = length(unique(dd2$s_idx)),
                 DI = dd2$DI,
                 sd_DI = dd2$sd_DI,
                 sy_idx = dd2$sy_idx,
                 s_idx = d_s$s_idx, #s_idx for each site/year,
                 s_idx2 = dd2$s_idx, #s_idx for each site/year/DOY
                 doy1 = poly_doy[,1],
                 doy2 = poly_doy[,2],
                 SIC_year = sic3$sic_year_sc,
                 SIC_site = sic4$sic_site,
                 bath = bath2$shelf_area_sc,
                 pro_data = dd2)

options(mc.cores = parallel::detectCores())

# DELTA <- 0.92
# TREE_DEPTH <- 10
# STEP_SIZE <- 0.005
CHAINS <- 4
ITER <- 3000

#compile model
mod <- cmdstanr::cmdstan_model(paste0(dir, 
                                          'Scripts/model_files/5-diet-model.stan'),
                                   force_recompile = TRUE)

print("STARTING TO FIT MODEL")

#sample
mod_fit <- mod$sample(
  data = DATA_model,
  chains = CHAINS,
  iter_sampling = ITER / 2,
  iter_warmup = ITER / 2,
  parallel_chains = CHAINS,
  refresh = 500)
# step_size = STEP_SIZE
# max_treedepth = TREE_DEPTH
# adapt_delta = DELTA


# PPC sim --------------------------------------------------------

#returns list with y, cn_id, sp, y_rep, and mu
PPC_fun <- function(Nsim = 100,
                    seed = 1)
{
  print('Extracting posteriors')
  # y2 ~ normal(z, tau)
  v_ch <- MCMCvis::MCMCchains(mod_fit, params = 'v')

  #random iterations to draw from posteriors
  idx <- sample(1:NROW(v_ch), Nsim)

  print('Simulating data')
  DI_rep <- matrix(NA, nrow = Nsim, ncol = DATA_model$Nsyd)

  set.seed(seed)
  for (i in 1:length(idx))
  {
    print(paste0('iteration: ', i, ' of ', length(idx)))
    eps <- rnorm(DATA_model$Nsyd, 0, DATA_model$sd_DI)
    DI_rep[i,] <- v_ch[i,] + eps
  }

  return(list(DI = DATA_model$DI,
              DI_rep = DI_rep))
}

#run function
ppc_sim <- PPC_fun(Nsim = 100)


# save summary space ------------------------------------------------------------

#save out summary, model fit, data
MCMCvis::MCMCdiag(mod_fit,
                  round = 4,
                  file_name = paste0('diet-model-results-', run_date),
                  dir = paste0(dir, 'Results'),
                  mkdir = paste0('diet-model-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('diet-model-fit-', run_date),
                  add_obj = list(DATA_model, ppc_sim),
                  add_obj_names = c(paste0('diet-model-data-', run_date),
                                    paste0('ppc-sim-', run_date)),
                  cp_file = c(paste0(dir, 'Scripts/Model_files/5-diet-model.stan'),
                              paste0(dir, 'Scripts/5-diet-model.R')),
                  cp_file_names = c(paste0('5-diet-model-', run_date, '.stan'),
                                    paste0('5-diet-model-', run_date, '.R')))

fig_dir <- paste0(dir, 'Results/diet-model-', run_date, '/')

# result_date <- '2025-07-17'
# mod_fit <- readRDS(paste0(dir, '/Results/diet-model-', result_date,
#                       '/diet-model-fit-', result_date, '.rds'))
# dd2 <- readRDS(paste0(dir, '/Results/diet-model-', result_date,
#                       '/diet-model-data-', result_date, '.rds'))$pro_data


# variables ---------------------------------------------------------------

# v = latent DI for each site/year/DOY
# alpha = diet for each site/year (at Dec 31, scaled day 0 AKA day 60)
# gamma = diet for each site (at mean site SIC)
# kappa = diet across all sites (grand mean)

# beta1 = linear effect of DOY for each site
# beta2 = quadratic effect of DOY for each site
# theta = effect of SIC on diet over time (across years) for each site
# zeta1 = effect of SIC on diet over space
# zeta2 = effect of bathymetry (shelf area) on diet over space


# PPC ---------------------------------------------------------------------

#posterior predictive check
y_val <- DATA_model$DI

pdf(paste0(fig_dir, 'PPC.pdf'), height = 5, width = 5)
plot(density(ppc_sim$DI), lwd = 2, main = 'PPC', xlab = 'Value', ylim = c(0, 0.02))
for (i in 1:NROW(ppc_sim$DI_rep))
{
  lines(density(ppc_sim$DI_rep[i,]), col = rgb(1,0,0,0.1))
}
dev.off()


# cat plots ---------------------------------------------------------------

pdf(paste0(fig_dir, 'DOY_linear-', run_date, '.pdf'),
    height = 10, width = 5)
MCMCvis::MCMCplot(mod_fit,
                  params = 'beta1',
                  labels = unique(d_s$site_id),
                  sz_labels = 0.25,
                  ci = c(89, 89),
                  sz_thick = 3,
                  sz_thin = 3,
                  main = 'DOY - Linear effect',
                  guide_lines = TRUE)
dev.off()

pdf(paste0(fig_dir, 'DOY_quadratic-', run_date, '.pdf'),
    height = 10, width = 5)
MCMCvis::MCMCplot(mod_fit,
                  params = 'beta2',
                  labels = unique(d_s$site_id),
                  sz_labels = 0.25,
                  ci = c(89, 89),
                  sz_thick = 3,
                  sz_thin = 3,
                  main = 'DOY - Quadratic effect',
                  guide_lines = TRUE)
dev.off()

pdf(paste0(fig_dir, 'SIC_effect-', run_date, '.pdf'),
    height = 10, width = 5)
MCMCvis::MCMCplot(mod_fit,
                  params = 'theta',
                  labels = unique(d_s$site_id),
                  sz_labels = 0.25,
                  ci = c(89, 89),
                  sz_thick = 3,
                  sz_thin = 3,
                  main = 'SIC effect',
                  guide_lines = TRUE)
dev.off()


# diet ~ DOY individual sites ---------------------------------------------------

# gamma_j + beta1_j * doy1 + beta2_j * doy2

gamma_mn <- MCMCvis::MCMCpstr(mod_fit, params = 'gamma')[[1]]
beta1_mn <- MCMCvis::MCMCpstr(mod_fit, params = 'beta1')[[1]]
beta2_mn <- MCMCvis::MCMCpstr(mod_fit, params = 'beta2')[[1]]

coefs_poly <- attr(poly_doy, 'coefs')

#simulate data for overall fit
out_new_doy <- array(NA, c(length(gamma_mn), length(-60:60)))
new_doy <- -60:60
poly_new_doy <- poly(new_doy, 2, raw = FALSE, coefs = coefs_poly)

#get pred values
for (i in 1:length(gamma_mn))
{
  #i <- 1
  out_new_doy[i,] <- gamma_mn[i] +
    (beta1_mn[i] * poly_new_doy[,1]) +
    (beta2_mn[i] * poly_new_doy[,2])
}

#reshape before plotting
tplt <- data.frame(site_id = rep(unique(d_s$site_id), length(-60:60)), 
                   DOY = rep(-60:60, each = length(unique(d_s$site_id))),
                   pred = reshape2::melt(out_new_doy)[,3])#,
                   # gamma_mn = rep(gamma_mn, length(-60:60)))

#plot 
p1 <- ggplot(tplt, aes(DOY, pred, group = factor(site_id))) +
# p1 <- ggplot(tplt, aes(DOY, pred, group = factor(site_id), color = gamma_mn)) +
  geom_line(alpha = 0.3, linewidth = 1.5) +
  #geom_point(data = data_tplt, aes(sc_jg, mn_p), alpha = 0.1) +
  theme_bw() +
  theme(legend.position = 'none') +
  ylab('DI') +
  xlab('DOY') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

ggsave(p1, filename = paste0(fig_dir, 'DOY_ind-', run_date, '.pdf'),
       width = 5, height = 5)


# diet ~ DOY cross-sites ------------------------------------------------------

# kappa + mu_beta1 * doy1 + mu_beta2 * doy2

kappa_ch <- MCMCvis::MCMCchains(mod_fit, params = 'kappa')
mu_beta1_ch <- MCMCvis::MCMCchains(mod_fit, params = 'mu_beta1')
mu_beta2_ch <- MCMCvis::MCMCchains(mod_fit, params = 'mu_beta2')

coefs_poly <- attr(poly_doy, 'coefs')

#simulate data for overall fit
out_new_doy2 <- array(NA, c(NROW(kappa_ch), length(-60:60)))
new_doy <- -60:60
poly_new_doy <- poly(new_doy, 2, raw = FALSE, coefs = coefs_poly)

#get pred values
for (i in 1:NROW(kappa_ch))
{
  print(paste0('processing iter: ', i, ' of ', NROW(kappa_ch)))
  out_new_doy2[i,] <- kappa_ch[i,1] +
    (mu_beta1_ch[i,1] * poly_new_doy[,1]) +
    (mu_beta2_ch[i,1] * poly_new_doy[,2])
}

#reshape before plotting
tplt2 <- data.frame(DOY = -60:60,
                    pred_mn = apply(out_new_doy2, 2, mean),
                    pred_LCI = apply(out_new_doy2, 2, 
                                     function(x) quantile(x, probs = 0.055)),
                    pred_UCI = apply(out_new_doy2, 2, 
                                     function(x) quantile(x, probs = 0.945)))

#plot 
p2 <- ggplot(tplt2) +
  #p1 <- ggplot(pred_tplt, aes(async_index, pro_index_mean, color = factor(sci_name))) +
  geom_ribbon(aes(x = DOY, ymin = pred_LCI, ymax = pred_UCI),
              alpha = 0.35) +
  geom_line(aes(DOY, pred_mn), 
            alpha = 0.8, 
            linewidth = 1.5) +
  #geom_point(data = data_tplt, aes(sc_jg, mn_p), alpha = 0.1) +
  theme_bw() +
  theme(legend.position = 'none') +
  ylab('DI') +
  xlab('DOY') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

ggsave(p2, filename = paste0(fig_dir, 'DOY_cross-', run_date, '.pdf'),
       width = 5, height = 5)



# diet ~ DOY cross and individual ----------------------------------------------

p3 <- ggplot() +
  #p1 <- ggplot(pred_tplt, aes(async_index, pro_index_mean, color = factor(sci_name))) +
  geom_line(data = tplt, 
            aes(DOY, pred, group = factor(site_id)),
            alpha = 0.15, linewidth = 1.1) +
  geom_line(data = tplt2, 
            aes(DOY, pred_mn), 
            alpha = 0.8, 
            linewidth = 2.2) +
  theme_bw() +
  theme(legend.position = 'none') +
  ylab('DI') +
  xlab('DOY') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18)) +
  ylim(c(0, 101))

ggsave(p3, filename = paste0(fig_dir, 'DOY_cross_ind-', run_date, '.pdf'),
       width = 5, height = 5)


# diet ~ SIC (over time) individual sites --------------------------------------

# gamma_j + theta_j * SIC_year

gamma_mn <- MCMCvis::MCMCpstr(mod_fit, params = 'gamma')[[1]]
theta_mn <- MCMCvis::MCMCpstr(mod_fit, params = 'theta')[[1]]

#simulate data for overall fit
out_new_sic <- array(NA, c(length(gamma_mn), 50))
out_ss_sic <- array(NA, c(length(gamma_mn), 50))
act_ss_sic <- array(NA, c(length(gamma_mn), 50))
new_sic <- seq(min(sic3$sic_year_sc), max(sic3$sic_year_sc), length = 50)
act_sic <- seq(min(sic3$sic_year_sc + sic3$mn_sic), 
               max(sic3$sic_year_sc + sic3$mn_sic), length = 50) 

sid <- unique(sic3$site_id)
#get pred values
for (i in 1:length(gamma_mn))
{
  #i <- 1
  #entire range SIC
  out_new_sic[i,] <- gamma_mn[i] + theta_mn[i] * new_sic
  
  #site-specific range
  tsic <- dplyr::filter(sic3, site_id == sid[i])
  #sic
  act_ss_sic[i,] <- c(tsic$sic_year_sc + tsic$mn_sic, rep(NA, 50 - length(tsic$sic_year_sc)))
  #DI
  out_ss_sic[i,] <- c(gamma_mn[i] + theta_mn[i] * tsic$sic_year_sc, rep(NA, 50 - length(tsic$sic_year_sc)))
}


#reshape before plotting
tplt_sic <- data.frame(site_id = rep(unique(d_s$site_id), length(new_sic)),
                       SIC = rep(act_sic, each = length(unique(d_s$site_id))),
                       pred = reshape2::melt(out_new_sic)[,3])

# #plot 
p4 <- ggplot(tplt_sic, aes(SIC, pred, group = factor(site_id))) +
  #p1 <- ggplot(pred_tplt, aes(async_index, pro_index_mean, color = factor(sci_name))) +
  geom_line(alpha = 0.3, linewidth = 1.5) +
  #geom_point(data = data_tplt, aes(sc_jg, mn_p), alpha = 0.1) +
  theme_bw() +
  theme(legend.position = 'none') +
  ylab('DI') +
  xlab('SIC') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

ggsave(p4, filename = paste0(fig_dir, 'SIC_ind-', run_date, '.pdf'),
       width = 5, height = 5)


# diet ~ SIC (over time) cross sites --------------------------------------

# kappa + mu_theta * SIC_year

kappa_ch <- MCMCvis::MCMCchains(mod_fit, params = 'kappa')
mu_theta_ch <- MCMCvis::MCMCchains(mod_fit, params = 'mu_theta')

#simulate data for overall fit
out_new_sic2 <- array(NA, c(length(kappa_ch), 50))

#get pred values
for (i in 1:NROW(kappa_ch))
{
  print(paste0('processing iter: ', i, ' of ', NROW(kappa_ch)))
  out_new_sic2[i,] <- kappa_ch[i,1] +
    mu_theta_ch[i,1] * new_sic
}

#reshape before plotting
tplt_sic2 <- data.frame(SIC = act_sic,
                        pred_mn = apply(out_new_sic2, 2, mean),
                        pred_LCI = apply(out_new_sic2, 2, 
                                         function(x) quantile(x, probs = 0.055)),
                        pred_UCI = apply(out_new_sic2, 2, 
                                         function(x) quantile(x, probs = 0.945)))

#plot
p5 <- ggplot(tplt_sic2) +
  #p1 <- ggplot(pred_tplt, aes(async_index, pro_index_mean, color = factor(sci_name))) +
  geom_ribbon(aes(x = SIC, ymin = pred_LCI, ymax = pred_UCI),
              alpha = 0.35) +
  geom_line(aes(SIC, pred_mn),
            alpha = 0.8,
            linewidth = 1.5) +
  #geom_point(data = data_tplt, aes(sc_jg, mn_p), alpha = 0.1) +
  theme_bw() +
  theme(legend.position = 'none') +
  ylab('DI') +
  xlab('SIC') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

ggsave(p5, filename = paste0(fig_dir, 'SIC_cross-', run_date, '.pdf'),
       width = 5, height = 5)


# diet ~ SIC (over time) cross and individual ----------------------------------

p6 <- ggplot() +
  geom_line(data = tplt_sic, 
            aes(SIC, pred, group = factor(site_id)),
            alpha = 0.1, linewidth = 1.1) +
  # geom_ribbon(data = tplt_sic2,
  #             aes(x = SIC, ymin = pred_LCI, ymax = pred_UCI),
  #             alpha = 0.35) +
  geom_line(data = tplt_sic2, 
            aes(SIC, pred_mn),
            alpha = 0.8,
            linewidth = 2.2) +
  #geom_point(data = data_tplt, aes(sc_jg, mn_p), alpha = 0.1) +
  theme_bw() +
  theme(legend.position = 'none') +
  ylab('DI') +
  xlab('SIC') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18)) +
  ylim(c(0, 101))

ggsave(p6, filename = paste0(fig_dir, 'SIC_ind_cross-', run_date, '.pdf'),
       width = 5, height = 5)


# diet variation over time ------------------------------------------------

d_s2 <- d_s
d_s2$alpha_mn <- MCMCvis::MCMCpstr(mod_fit, params = 'alpha')[[1]]

# ggplot(d_s2, aes(season, alpha_mn, color = site_id)) +
#   geom_point() +
#   # geom_line(stat = 'smooth',
#   #           method = lm,
#   #           linewidth = 1.2,
#   #           alpha = 0.5) +
#   geom_line(linewidth = 1.2,
#             alpha = 0.5) +
#   theme(legend.position = 'none')


# diet ~ SIC (over time) with data --------------------------------------

out_ss_sic_rs <- reshape2::melt(out_ss_sic)
colnames(out_ss_sic_rs) <- c('site_n', 'year_slot', 'pred_DI')

act_ss_sic_rs <- reshape2::melt(act_ss_sic)
colnames(act_ss_sic_rs) <- c('site_n', 'year_slot', 'SIC')

out_ss_mrg <- dplyr::left_join(out_ss_sic_rs, act_ss_sic_rs, 
                               by = c('site_n', 'year_slot')) %>%
  dplyr::left_join(data.frame(site_id = unique(d_s$site_id), 
                              site_n = 1:length(unique(d_s$site_id))),
                   by = 'site_n') %>%
  dplyr::filter(!is.na(pred_DI)) %>%
  dplyr::arrange(site_id, year_slot)

out_ss_mrg$alpha_mn <- d_s2$alpha_mn


#reshape before plotting
ss_dl <- ggplot(data = out_ss_mrg, aes(SIC, alpha_mn, col = site_id)) +
  geom_point(alpha = 0.5) +
  geom_line(stat = 'smooth',
            method = lm,
            linewidth = 1.1,
            alpha = 0.6) +
  theme_bw() +
  theme(legend.position = 'none') +
  ylab('DI') +
  xlab('SIC') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

ggsave(ss_dl, filename = paste0(fig_dir, 'SIC_data_line-', run_date, '.pdf'),
       width = 5, height = 5)


# diet variation across space ---------------------------------------------

d_o <- dplyr::select(d_s, -sy_idx, -season) %>%
  dplyr::distinct()

#posteriors
d_o$gamma_mn <- MCMCvis::MCMCpstr(mod_fit, params = 'gamma')[[1]]
d_o$gamma_sd <- MCMCvis::MCMCpstr(mod_fit, params = 'gamma', fun = sd)[[1]]
d_o$theta_mn <- MCMCvis::MCMCpstr(mod_fit, params = 'theta')[[1]]
d_o$theta_sd <- MCMCvis::MCMCpstr(mod_fit, params = 'theta', fun = sd)[[1]]
d_o$beta1_mn <- MCMCvis::MCMCpstr(mod_fit, params = 'beta1')[[1]]
d_o$beta2_mn <- MCMCvis::MCMCpstr(mod_fit, params = 'beta2')[[1]]

#read in shp file
Ant <- sf::st_read(paste0(dir, 'Data/L0/Coastline_low_res_polygon/Coastline_low_res_polygon.shp'))

#plot data
AC <- ggplot() + 
  theme_void() +
  geom_sf(data = Ant, 
          fill = 'grey',
          color = 'black')

# BREAKS <- c(0, 30, 35, 40, 45, 50, 55, 100)
BREAKS <- quantile(d_o$gamma_mn, probs = seq(0, 1, by = 0.2))
BREAKS[1] <- 0
BREAKS[length(BREAKS)] <- 100

pmap <- AC +
  geom_point(data = d_o,
             inherit.aes = FALSE,
             size = 8,
             alpha = 0.5,
             aes(site_easting, site_northing, 
                 color = cut(gamma_mn, breaks = BREAKS))) +
  scale_color_viridis_d('Diet') +
  geom_point(data = d_o,
             inherit.aes = FALSE,
             size = 8,
             shape = 21,
             alpha = 0.4,
             stroke = 1,
             color = 'black',
             aes(site_easting, site_northing))

ggsave(pmap, filename = paste0(fig_dir, 'diet_map-', run_date, '.pdf'),
       width = 5, height = 5)

#just continent and sites
smap <- ggplot() + 
  theme_void() +
  geom_sf(data = Ant, 
          fill = 'black',
          color = 'black') +
  geom_point(data = d_o,
             inherit.aes = FALSE,
             size = 8,
             alpha = 1,
             aes(site_easting, site_northing),
             color = '#2596be')

ggsave(smap, filename = paste0(fig_dir, 'site_map-', run_date, '.pdf'),
       width = 5, height = 5)

  
# diet ~ SIC space -------------------------------------------------

#kappa + zeta1 * SIC_site + zeta2 * bath
kappa_ch <- MCMCvis::MCMCchains(mod_fit, params = 'kappa')
zeta1_ch <- MCMCvis::MCMCchains(mod_fit, params = 'zeta1')
zeta2_ch <- MCMCvis::MCMCchains(mod_fit, params = 'zeta2')

#add site avg SIC
d_o$sic_site <- sic4$sic_site
#add site avg Bath
d_o$shelf_area <- bath2$shelf_area #used in model (scaled version)
d_o$shelf_area_sc <- bath2$shelf_area_sc

sic_sim <- seq(min(d_o$sic_site), max(d_o$sic_site), length.out = 50)
mu_sic <- array(NA, c(NROW(kappa_ch), 50))
for (i in 1:NROW(kappa_ch))
{
  mu_sic[i,] <- kappa_ch[i,] + zeta1_ch[i,] * sic_sim
}

#model fit SIC
tplt7 <- data.frame(np = sic_sim,
                    mu_sic_LCI = apply(mu_sic, 2, function(x) quantile(x, probs = 0.055)),
                    mu_sic_UCI = apply(mu_sic, 2, function(x) quantile(x, probs = 0.945)),
                    mu_sic_mean = apply(mu_sic, 2, mean))

#SIC plot
p7 <- ggplot() +
  geom_errorbar(data = d_o,
                aes(x = sic_site,
                    ymin = gamma_mn - gamma_sd,
                    ymax = gamma_mn + gamma_sd),
                alpha = 0.4) +
  geom_point(data = d_o,
             aes(x = sic_site,
                 y = gamma_mn),
             alpha = 0.4,
             size = 2.5) +
  geom_ribbon(data = tplt7,
              aes(x = np, ymin = mu_sic_LCI, ymax = mu_sic_UCI),
              alpha = 0.2) +
  geom_line(data = tplt7,
            aes(x = np, mu_sic_mean),
            alpha = 0.8,
            linewidth = 1.5) +
  theme_bw() +
  theme(legend.position = 'none') +
  ylab('DI') +
  xlab('SIC') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

ggsave(p7, filename = paste0(fig_dir, 'diet-SIC-', run_date, '.pdf'),
       width = 5, height = 5)


# diet ~ bath -------------------------------------------------

bath_sim <- seq(min(d_o$shelf_area_sc), max(d_o$shelf_area_sc), length.out = 50)
mu_bath <- array(NA, c(NROW(kappa_ch), 50))
#at avg SIC
for (i in 1:NROW(kappa_ch))
{
  mu_bath[i,] <- kappa_ch[i,] + + zeta1_ch[i] * mean(d_o$sic_site) + zeta2_ch[i,] * bath_sim
}

#model fit bath
#convert to actual shelf area
tplt8 <- data.frame(np = (bath_sim * sd(d_o$shelf_area)) + mean(d_o$shelf_area),
                    mu_bath_LCI = apply(mu_bath, 2, function(x) quantile(x, probs = 0.055)),
                    mu_bath_UCI = apply(mu_bath, 2, function(x) quantile(x, probs = 0.945)),
                    mu_bath_mean = apply(mu_bath, 2, mean))

#Bath plot
p8 <- ggplot() +
  geom_errorbar(data = d_o,
                aes(x = shelf_area,
                    ymin = gamma_mn - gamma_sd,
                    ymax = gamma_mn + gamma_sd),
                alpha = 0.4) +
  geom_point(data = d_o,
             aes(x = shelf_area,
                 y = gamma_mn),
             alpha = 0.4,
             size = 2.5) +
  geom_ribbon(data = tplt8,
              aes(x = np, ymin = mu_bath_LCI, ymax = mu_bath_UCI),
              alpha = 0.2) +
  geom_line(data = tplt8,
            aes(x = np, mu_bath_mean),
            alpha = 0.8,
            linewidth = 1.5) +
  theme_bw() +
  theme(legend.position = 'none') +
  ylab('DI') +
  xlab('Shelf Area') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

ggsave(p8, filename = paste0(fig_dir, 'diet-bath-', run_date, '.pdf'),
       width = 5, height = 5)


# site/year ---------------------------------------------------------------

d_s2_plt <- dplyr::right_join(d_s2, tidyr::expand(d_s2, site_id, 
                                               season = min(dd2$season):max(dd2$season)), 
                            by = c('site_id', 'season')) %>%
  dplyr::arrange(site_id, season) %>%
  #negative lats, abs and add 180
  dplyr::mutate(lsrt = if_else(longitude < 0, longitude + 360, longitude)) %>%
  dplyr::arrange(lsrt, season)

#reorder factors
d_s2_plt$site_id <- factor(d_s2_plt$site_id, levels = unique(d_s2_plt$site_id))

#AP at the top, then counterclockwise around continent
sy_plt <- ggplot() +
  geom_tile(data = d_s2_plt, aes(x = season, y = site_id, fill = alpha_mn),
            color = 'black', linewidth = 0.5) +
  scale_fill_viridis_c('Diet',
                       limits = range(d_s2_plt$alpha_mn),
                       na.value = 'white') +
  theme_minimal()

ggsave(sy_plt, filename = paste0(fig_dir, 'site-year-diet-', run_date, '.pdf'),
       width = 4, height = 10)

# dplyr::select(d_s2_plt, site_id, ccamlr_id, region, latitude, longitude, lsrt) %>%
#   dplyr::filter(!is.na(ccamlr_id)) %>%
#   dplyr::distinct() %>%
#   dplyr::arrange(desc(lsrt)) %>%
#   print(n = 500)


# SIC ~ time --------------------------------------------------------------

sic_time <- ggplot() +
  annotate('rect', 
           xmin = 1985, xmax = 2012, ymin = 0, ymax = 1,
           alpha = 0.3) +
  geom_line(data = sic2b, aes(season, sic_year, group = site_id),
            stat = 'smooth',
            # method = lm,
            linewidth = 1,
            alpha = 0.5) +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(legend.position = 'none') +
  ylab('SIC') +
  xlab('Year') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

ggsave(sic_time, filename = paste0(fig_dir, 'sic-year-gam-', run_date, '.pdf'),
       width = 5, height = 5)


# percent var -------------------------------------------------------------

#calculate how much variation over range of SIC each location
di_rng <- apply(out_ss_sic, 1, function(x) diff(range(x, na.rm = TRUE)))
#calculate variance at each location
di_var <- apply(out_ss_sic, 1, function(x) var(x, na.rm = TRUE))
#remote site with one year of data
di_rng2 <- di_rng[di_rng > 0]
di_var2 <- di_var[!is.na(di_var)]
hist(di_rng2)
mean(di_rng2)
hist(di_var2, breaks = 20)
mean(di_var2)

#variance across sites
vs <- var(d_o$gamma_mn)
frac_v <- di_var2 / vs
hist(frac_v)
mean(frac_v)
median(frac_v)
rngs <- diff(range(d_o$gamma_mn))
frac_r <- di_rng2 / rngs
hist(frac_r)
mean(frac_r)
median(frac_r)

sid <- unique(d_s2$site_id)
pv <- rep(NA, length(sid))
for (i in 1:length(sid))
{
  #i <- 1
  tt2 <- dplyr::filter(d_s2, site_id == sid[i])
  
  #var fit / (var fit + var res)
  pv[i] <- di_var[i] / (di_var[i] + var(tt2$alpha_mn - out_ss_sic[i,1:NROW(tt2)]))
}

hist(pv, breaks = 30)
mean(pv, na.rm = TRUE)
median(pv, na.rm = TRUE)

