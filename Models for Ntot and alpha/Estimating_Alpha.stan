data {
  int<lower=0> N;                  // Number of observations
  real<lower=0> x_min;             // Minimum threshold for the distribution (e.g., 10)
  real<lower=0> trunc_point;       // Observed data starts at this truncation threshold (e.g., 20)
  real<lower=trunc_point> trunc_upper;  // Observed data is capped at this upper threshold (e.g., 50)
  vector<lower=trunc_point>[N] x;  // Observed data, limited to range [trunc_point, trunc_upper]
  real<lower=0, upper=1> LAI_norm; // Normalized LAI score (0 to 1)
  real<lower=0, upper=1> breakpoint_norm; // Normalized breakpoint distance (0 to 1)
  real prior_mean;                 // Mean of the truncated normal prior for alpha
  real prior_sd;                   // Standard deviation of the truncated normal prior for alpha
}

parameters {
  real<lower=0> alpha;             // Shape parameter for Pareto distribution (must be positive)
}

transformed parameters {
  // Adjustment factor: Penalizes the likelihood based on LAI and breakpoint distance
  real adjustment_factor = 1 - sqrt(LAI_norm * breakpoint_norm);
}

model {
  // Truncated normal prior for alpha, with lower bound at 0
  alpha ~ normal(prior_mean, prior_sd) T[0, ]; // Truncated normal with lower bound at 0
  
  // Truncated cumulative probability in the observed range
  real p_trunc = pareto_cdf(trunc_upper | x_min, alpha) - pareto_cdf(trunc_point | x_min, alpha);
  
  // Likelihood: Adjusted for observational biases
  for (n in 1:N) {
    target += adjustment_factor * (pareto_lpdf(x[n] | x_min, alpha) - log(p_trunc));
  }
}