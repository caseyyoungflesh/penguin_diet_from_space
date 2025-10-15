
data {
  int<lower=0> Nsyd;    // number of site, years, DOY
  int<lower=0> Nsy;     // number of site, years
  int<lower=0> Ns;      // number of sites
  vector[Nsyd] DI;      // posterior mean DI
  vector[Nsyd] sd_DI;   // associated uncertainty
  array[Nsyd] int<lower=0, upper=Nsy> sy_idx;    // site, year index for each site, year, DOY
  array[Nsy] int<lower=0, upper=Ns> s_idx;       // site index for each site, year
  array[Nsyd] int<lower=0, upper=Ns> s_idx2;     // site index for each site, year, DOY
  vector[Nsyd] doy1;    // DOY - centered at particular day (mid-point of season) - poly 1 (linear)
  vector[Nsyd] doy2;    // DOY - centered at particular day (mid-point of season) - poly 2 (quadratic)
  vector[Nsy] SIC_year; // SIC for each site
  vector[Ns] SIC_site;  // SIC for each site
  vector[Ns] bath;      // bathymetry for each site
}

parameters {
  vector[Nsyd] v_raw;
  // vector[Nsy] alpha;    // diet for site, year (at midpoint of season)
  vector[Nsy] alpha_raw;
  real<lower=0> sigma_v;
  vector[Ns] gamma;        // diet for site (at midpoint of time series)
  // vector[Ns] theta;     // effect of cov on diet
  vector[Ns] theta_raw;
  real<lower=0> sigma_alpha;
  real mu_theta;
  real<lower=0> sigma_theta;
  real mu_beta1;
  real mu_beta2;
  vector<lower=0>[2] sigma_bb;
  matrix[2, Ns] z_bb;                     // z-scores
  cholesky_factor_corr[2] L_Rho_bb;       // cholesky factor of corr matrix
  real kappa;
  real zeta1;
  real zeta2;
  real<lower=0> sigma_gamma;
}

transformed parameters {
  // vector[Ns] theta;
  vector[Nsy] alpha;
  vector[Ns] beta1;
  vector[Ns] beta2;
  vector[Ns] theta;                      // effect of cov on diet
  matrix[Ns, 2] bb;
  matrix[2, 2] Rho_bb;                   // corr matrix
  vector[Nsyd] v;
  
  // implies theta ~ normal(mu_theta, sigma_theta);
  theta = theta_raw * sigma_theta + mu_theta;
  
  // implies alpha ~ normal(gamma[s_idx] + theta[s_idx] .* SIC, sigma_alpha);
  alpha = alpha_raw * sigma_alpha + gamma[s_idx] + theta[s_idx] .* SIC_year;
  
  // cholesky factor of covariance matrix multiplied by z score
  // implies gt ~ MVN(0, sigma)
  bb = (diag_pre_multiply(sigma_bb, L_Rho_bb) * z_bb)';
  // implies Rho = L_Rho * L_Rho';
  Rho_bb = multiply_lower_tri_self_transpose(L_Rho_bb);
  
  beta1 = mu_beta1 + bb[,1];
  beta2 = mu_beta2 + bb[,2];
  
  // implies v ~ N(alpha[sy_idx] + beta1[s_idx2] .* doy1 + beta2[s_idx2] .* doy2, sigma_v);
  v = v_raw * sigma_v + alpha[sy_idx] + beta1[s_idx2] .* doy1 + beta2[s_idx2] .* doy2;
}

model {
  //priors
  sigma_v ~ normal(0, 5);
  sigma_alpha ~ normal(0, 5);
  mu_theta ~ normal(0, 20);
  sigma_theta ~ normal(0, 15);
  mu_beta1 ~ normal(0, 500);
  mu_beta2 ~ normal(0, 500);
  sigma_bb ~ normal(0, 200);
  kappa ~ normal(30, 20);
  zeta1 ~ normal(0, 20);
  zeta2 ~ normal(0, 10);
  sigma_gamma ~ normal(0, 15);
  
  // for non-centered
  // theta_raw ~ std_normal();
  alpha_raw ~ std_normal();
  v_raw ~ std_normal();
  theta_raw ~ std_normal();
  
  to_vector(z_bb) ~ std_normal();
  L_Rho_bb ~ lkj_corr_cholesky(2);
  
  // likelihood
  // see above for non-centered
  // theta ~ normal(mu_theta, sigma_theta);
  
  // gamma_j ~ (kappa + zeta1 * SIC + zeta2 * bathymetry) - how diet varies around continent
  gamma ~ normal(kappa + zeta1 * SIC_site + zeta2 * bath, sigma_gamma);

  // see above for non-centered
  // alpha_jm ~ N(gamma_j + theta_j * year, sigma_alpha) - how diet changes over time
  // alpha ~ normal(gamma[s_idx] + theta[s_idx] .* year, sigma_alpha);
  
  // see above for non-centered
  // v_jkm ~ N(alpha_jm + beta1_j * DOY_jkm + beta2 * DOY_jkm^2, sigma) - trend within year; alpha = diet for site/year
  // v ~ normal(alpha[sy_idx] + beta[sy_idx] .* doy, sigma_v);
  
  DI ~ normal(v, sd_DI);
}

