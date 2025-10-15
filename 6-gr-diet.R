######################
# 6 - growth rate ~ diet
######################


# Load packages -----------------------------------------------------------

# devtools::install_github('trafficonese/leaflet.esri')
# devtools::install_github('CCheCastaldo/mapppdr')

library(MCMCvis)
library(tidyverse)
library(cmdstanr)
library(mapppdr)
library(here)


# set dir -----------------------------------------------------------------

dir <- paste0(here::here(), '/')
run_date <- Sys.Date()
diet_model_run_date <- '2025-07-17'
CCC_run_date <- '2025-07-16'


# read in data -------------------------------------------------------------

#output from 5-diet-model.R
p_fit <- readRDS(paste0(dir, '/Results/diet-model-', diet_model_run_date,
                      '/diet-model-fit-', diet_model_run_date, '.rds'))
p_data <- readRDS(paste0(dir, '/Results/diet-model-', diet_model_run_date,
                      '/diet-model-data-', diet_model_run_date, '.rds'))$pro_data

#from Che-Castaldo et al. 2017 Nature Communications
nsl <- readRDS(paste0(dir, 'Data/L2/Che-Castaldo-et-al-Nat-Comms-', CCC_run_date, '.rds'))


# extract and process -----------------------------------------------------

#avg diet
gamma <- MCMCvis::MCMCpstr(p_fit, params = 'gamma')[[1]]
gamma_sd <- MCMCvis::MCMCpstr(p_fit, params = 'gamma', fun = sd)[[1]]

dd8 <- dplyr::select(p_data, site_id, ccamlr_id, region, s_idx,
                     avg_n_pix, site_easting, site_northing) %>%
  dplyr::distinct() %>%
  #merge with full name and lat lon
  dplyr::left_join(mapppdr::sites,
                   by = c('site_id', 'ccamlr_id', 'region')) %>%
  #merge with growth rate
  dplyr::left_join(nsl, by = 'site_id')

pg_na <- which(is.na(dd8$pg))

#9 sites with no growth info
length(pg_na)

#remove from data
dd8_2 <- dd8[-pg_na,]
dd8_2$gamma <- gamma[-pg_na]
dd8_2$gamma_sd <- gamma_sd[-pg_na]


# fit model ---------------------------------------------------------------

DATA <- list(N = NROW(dd8_2),
             pg = dd8_2$pg,
             sd_pg = dd8_2$pg_sd,
             gamma = dd8_2$gamma,
             sd_gamma = dd8_2$gamma_sd,
             pro_data = dd8_2)

options(mc.cores = parallel::detectCores())

# DELTA <- 0.92
# TREE_DEPTH <- 10
# STEP_SIZE <- 0.005
CHAINS <- 4
ITER <- 9000

#compile model
mod <- cmdstanr::cmdstan_model(paste0(dir, 
                                      'Scripts/model_files/6-gr-diet.stan'),
                               force_recompile = TRUE)

print("STARTING TO FIT MODEL")

#sample
mod_fit <- mod$sample(
  data = DATA,
  chains = CHAINS,
  iter_sampling = ITER / 2,
  iter_warmup = ITER / 2,
  parallel_chains = CHAINS,
  refresh = 500)
# step_size = STEP_SIZE
# max_treedepth = TREE_DEPTH
# adapt_delta = DELTA


# save summary space ------------------------------------------------------------

#save out summary, model fit, data
MCMCvis::MCMCdiag(mod_fit,
                  round = 4,
                  file_name = paste0('gr-diet-results-', run_date),
                  dir = paste0(dir, 'Results'),
                  mkdir = paste0('gr-diet-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('gr-diet-fit-', run_date),
                  add_obj = list(DATA),
                  add_obj_names = paste0('gr-diet-data-', run_date),
                  cp_file = c(paste0(dir, 'Scripts/Model_files/8-gr-diet.stan'),
                              paste0(dir, 'Scripts/8-gr-diet.R')),
                  cp_file_names = c(paste0('6-gr-diet-', run_date, '.stan'),
                                    paste0('6-gr-diet-', run_date, '.R')))

fig_dir <- paste0(dir, 'Results/gr-diet-', run_date, '/')

# result_date <- '2025-07-17'
# fig_dir <- paste0(dir, 'Results/gr-diet-', result_date, '/')
# mod_fit <- readRDS(paste0(dir, '/Results/gr-diet-', result_date,
#                       '/gr-diet-fit-', result_date, '.rds'))
# dd2 <- readRDS(paste0(dir, '/Results/gr-diet-pc-', result_date,
#                       '/gr-diet-data-', result_date, '.rds'))$pro_data


# PPC ---------------------------------------------------------------------

#posterior predictive check - response
y_val <- DATA$pg
ppc_sim <- MCMCvis::MCMCchains(mod_fit, params = 'pg_rep')

pdf(paste0(fig_dir, 'PPC.pdf'), height = 5, width = 5)
plot(density(y_val), lwd = 2, main = 'PPC', xlab = 'Value', 
     xlim = c(0.89, 1.12), ylim = c(0, 15))
for (i in 1:100)
{
  lines(density(ppc_sim[i,]), col = rgb(1,0,0,0.1))
}
dev.off()

#posterior predictive check - predictor
x_val <- DATA$gamma
ppc_sim_cov <- MCMCvis::MCMCchains(mod_fit, params = 'gamma_rep')

pdf(paste0(fig_dir, 'PPC_cov.pdf'), height = 5, width = 5)
plot(density(x_val), lwd = 2, main = 'PPC', xlab = 'Value', 
     ylim = c(0, 0.025))
for (i in 1:100)
{
  lines(density(ppc_sim_cov[i,]), col = rgb(1,0,0,0.1))
}
dev.off()


# gr ~ diet ------------------------------------------------------

# alpha + beta * kappa
mu_pg_ch <- MCMCvis::MCMCchains(mod_fit, params = 'mu_pg')
kappa_ch <- MCMCvis::MCMCchains(mod_fit, params = 'kappa')
alpha_ch <- MCMCvis::MCMCchains(mod_fit, params = 'alpha')
beta_ch <- MCMCvis::MCMCchains(mod_fit, params = 'beta')

tplt <- data.frame(gamma_mn = DATA$gamma,
                   gamma_sd = DATA$sd_gamma,
                   # DI_mn = DATA$pc + mean(eta_mn),
                   # DI_sd = DATA$sd_pc,
                   pg_mn = DATA$pg,
                   pg_sd = DATA$sd_pg)

rng_gamma <- range(c(tplt$gamma_mn - tplt$gamma_sd),
                   c(tplt$gamma_mn + tplt$gamma_sd))
np <- seq(rng_gamma[1], rng_gamma[2], length.out = 50)

mu_pred <- array(NA, c(NROW(alpha_ch), 50))
for (i in 1:NROW(alpha_ch))
{
  mu_pred[i,] <- alpha_ch[i,] + beta_ch[i,] * np
}

#model fit
tplt2 <- data.frame(np = np,
                   mu_pred_LCI = apply(mu_pred, 2, function(x) quantile(x, probs = 0.055)),
                   mu_pred_UCI = apply(mu_pred, 2, function(x) quantile(x, probs = 0.945)),
                   mu_pred_mean = apply(mu_pred, 2, mean))

#plot
p1 <- ggplot() +
  geom_errorbar(data = tplt,
                aes(x = gamma_mn,
                    ymin = pg_mn - pg_sd,
                    ymax = pg_mn + pg_sd),
                alpha = 0.4) +
  geom_errorbarh(data = tplt,
                 aes(y = pg_mn,
                     xmin = gamma_mn - gamma_sd,
                     xmax = gamma_mn + gamma_sd),
                 alpha = 0.4) +
  geom_point(data = tplt,
             aes(x = gamma_mn,
                 y = pg_mn),
             alpha = 0.4) +
  geom_ribbon(data = tplt2,
              aes(x = np, ymin = mu_pred_LCI, ymax = mu_pred_UCI),
              alpha = 0.2) +
  geom_line(data = tplt2,
            aes(x = np, mu_pred_mean),
            alpha = 0.8,
            linewidth = 1.5) +
  theme_bw() +
  theme(legend.position = 'none') +
  ylab('Pop growth rate') +
  xlab('Mean diet') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linewidth = 2),
        axis.ticks = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

ggsave(p1, filename = paste0(fig_dir, 'gr-diet2-', run_date, '.pdf'),
       width = 6, height = 5)
