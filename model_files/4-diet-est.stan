
data {
  int<lower=0> N;
  int<lower=0> Nsyd;    // number of site, years, DOY
  int<lower=0> Ns;      // number of sites
  vector[N] y;          // predicted diet
  vector[N] sd_y;       // associated coefficient uncertainty
  real tau;             // model prediction uncertainty (from r2)
  array[N] int<lower=0, upper=Nsyd> syd_idx;   // site, year, DOY index for each site, year, DOY, cell
  array[Nsyd] int<lower=0, upper=Ns> s_idx2;   // site index for each site, year, DOY
}

parameters {
  // vector[N] y2;      // response accounting for model coefficient uncertainty
  vector[N] y2_raw;
  // vector[N] z;       // response accounting for model coefficient uncertainty and model pred uncertainty
  vector[N] z_raw;
  // vector[Nsyd] w;    // response aggregated to site/year/DOY
  vector[Nsyd] w_raw;
  vector[Ns] mu_w;
  real<lower=0> sigma;
  real<lower=0> sigma_w;
  real mm;
  real<lower=0> sm;
}

transformed parameters {
  vector[Nsyd] w;
  vector[N] z;
  vector[N] y2;
  
  // implies w ~ N(m_w[s_idx2], sigma_w);
  w = w_raw * sigma_w + mu_w[s_idx2];
  
  // implies z ~ N(w[syd_idx], sigma);
  z = z_raw * sigma + w[syd_idx];
  
  // implies y2 ~ normal(z, tau);
  y2 = y2_raw * tau + z;
}

model {
  //priors
  sigma ~ std_normal();
  sigma_w ~ normal(0, 2);
  mm ~ normal(30, 10);
  sm ~ normal(0, 3);
  
  // for non-centered
  w_raw ~ std_normal();
  z_raw ~ std_normal();
  y2_raw ~ std_normal();
  
  // mu_w_j ~ N(mm, sm);
  mu_w ~ normal(mm, sm);
  
  // see above for non-centered
  // z_ijkm ~ N(w_jkm, sigma) - aggregate to site/DOY/year (i.e., scene)
  // z ~ normal(w[syd_idx], sigma);
  
  // see above for non-centered
  // y2_ijkm ~ N(z_ijkm, tau) - model prediction uncertainty
  // y2 ~ normal(z, tau);

  // y1_ijkm ~ N(y2_ijkm, sd_y_ijkm) - model coefficient uncertainty
  y ~ normal(y2, sd_y);
}

