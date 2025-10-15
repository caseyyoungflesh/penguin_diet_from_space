
data {
  int<lower=0> N;       // number of sites
  vector[N] pg;         // posterior mean pop growth rates
  vector[N] sd_pg;      // associated uncertainty
  vector[N] gamma;      // site-level diet
  vector[N] sd_gamma;   // associated uncertainty
}

parameters {
  vector[N] mu_pg_raw;
  vector[N] kappa_raw;
  real alpha;
  real beta;
  real<lower=0> sigma;
  real mu_kappa;
  real<lower=0> sigma_kappa;
}

transformed parameters {
  vector[N] kappa;
  vector[N] mu_pg;
  
  // implies kappa_i ~ normal(mu_kappa, sigma_kappa);
  kappa = kappa_raw * sigma_kappa + mu_kappa;
  
  // implies mu_pg_i ~ normal(alpha + beta * kappa_i, sigma);
  mu_pg = mu_pg_raw * sigma + alpha + beta * kappa;
}

model {
  //priors
  sigma ~ std_normal();
  alpha ~ std_normal();
  beta ~ std_normal();
  mu_kappa ~ normal(40, 10);
  sigma_kappa ~ normal(10, 20);

  // for non-centered
  kappa_raw ~ std_normal();
  mu_pg_raw ~ std_normal();
  
  // kappa ~ normal(mu_kappa, sigma_kappa);
  
  // observation model covariate
  gamma ~ normal(kappa, sd_gamma);

  // process model
  // mu_pg ~ normal(alpha + beta * kappa, sigma);
  
  // observation model response
  pg ~ normal(mu_pg, sd_pg);
}

generated quantities {
  array[N] real pg_rep;
  array[N] real gamma_rep;
  
  pg_rep = normal_rng(mu_pg, sd_pg);
  gamma_rep = normal_rng(kappa, sd_gamma);
}
