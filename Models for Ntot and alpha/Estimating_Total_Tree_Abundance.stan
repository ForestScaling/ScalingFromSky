data {
  int<lower=1> K;                // Number of bins
  real<lower=0> bin_min[K];      // Lower bounds of DBH for bins
  real<lower=0> bin_max[K];      // Upper bounds of DBH for bins
  real<lower=0> N_obs[K];        // Observed counts per bin (continuous)
  real<lower=0> x_min;           // Minimum DBH threshold
  real<lower=0> x_max;           // Maximum DBH threshold
  int<lower=1> n_alpha_samples;  // Number of pre-sampled alphas
  real alpha_samples[n_alpha_samples]; // Pre-sampled alpha values
  real<lower=0> N_tot_prior_mean; // Mean of the prior for N_tot
  real<lower=0> N_tot_prior_sd;   // SD of the prior for N_tot
  real<lower=0, upper=1> breakpoint_norm; // Normalized breakpoint value
  real<lower=0, upper=1> LAI_norm;        // Normalized LAI value
}

parameters {
  real<lower=0> N_tot;           // Total number of trees with DBH >= x_min
  real<lower=0> sigma;           // Standard deviation for likelihood
}

model {
  vector[K] lambda;
  real log_lik[n_alpha_samples];

  // Compute adjustment factor
  real adjustment_factor = 1 - sqrt(LAI_norm * breakpoint_norm);
  real adjusted_N_tot = N_tot / adjustment_factor;

  // Priors
  N_tot ~ normal(N_tot_prior_mean, N_tot_prior_sd);
  sigma ~ normal(0, 1);

  // Marginalize over alpha uncertainty
  for (j in 1:n_alpha_samples) {
    real alpha_sample = alpha_samples[j];
    real trunc_correction = (x_min^(1 - alpha_sample) - x_max^(1 - alpha_sample));
    for (k in 1:K) {
      lambda[k] = adjusted_N_tot * (bin_min[k]^(1 - alpha_sample) - bin_max[k]^(1 - alpha_sample)) / trunc_correction;
    }
    log_lik[j] = sum(normal_lpdf(N_obs | lambda, sigma));
  }

  target += log_sum_exp(log_lik) - log(n_alpha_samples);
}