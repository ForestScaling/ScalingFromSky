##########################################
# Estimating Total Tree Abundance from   #
# Remote Sensing Data                   #
##########################################

##########################################
# Purpose:                              #
# We use observed tree size data (specifically, DBH—diameter at breast height) to #
# estimate \(N_{\text{tot}}\) using a **Bayesian approach**. This method allows #
# us to incorporate prior knowledge of \(N_{\text{tot}}\) and the **shape parameter**  #
# (\(\alpha\)) of the tree size distribution (assumed to follow a **Pareto distribution**) #
# into the estimation process. The observed data we work with include the number of trees #
# (\(N_{\text{obs}}\)) within a specific size range, bounded by a cutoff diameter (over #
# which remote sensing data is more accurate) and the maximum diameter.             #
##########################################


# Load necessary libraries
library(rstan)
library(dplyr)
library(forestscaling) #Not on CRAN, see Grady et al. 2024 (reference in
# citation list of main paper)
library(data.table)

# Define the Stan model
stan_totaltrees_model <- "data {
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
}"

# Compile the Stan model
stan_totaltrees_model_compile <- stan_model(model_code = stan_totaltrees_model)

# Load and process Leaf Area Index (LAI) data. This is obtained from Google Earth Engine (MODIS satellite)
Leaf_area_index <- data.frame(
  site = c("ABBY", "BART", "BLAN", "BONA", "CLBJ", "CPER", "DCFS", "DEJU", "DELA", 
           "DSNY", "GRSM", "GUAN", "HARV", "HEAL", "JERC", "JORN", "KONZ", "LAJA", 
           "LENO", "MLBS", "MOAB", "NIWO", "NOGP", "ONAQ", "ORNL", "OSBS", "PUUM", 
           "RMNP", "SCBI", "SERC", "SJER", "SOAP", "SRER", "STEI", "TALL", "TEAK", 
           "TREE", "UKFS", "UNDE", "WREF", "YELL"),
  
  Leaf_area_index = c(40.78947368, 49.46551724, 43.07142857, 34.88461538, 23.25, 
                      9.928571429, 20, 14.81034483, 43.91666667, NA, 50.03125, 
                      33.55, 54.26923077, 15.07692308, 16.80769231, 3, 36.72222222, 
                      NA, 51.5, 46.86111111, 3.777777778, 8.25, 16.66666667, 6.55, 
                      41.67647059, 16.08333333, NA, 12.78571429, 56.36363636, 47.875, 
                      11.73529412, 22.96153846, 3.076923077, 53.04545455, 51.52631579, 
                      13.05263158, 53.22, 46.73809524, 52.33333333, 37.61904762, 11.375)
)

# Remove rows with NA in the Leaf_area_index column
Leaf_area_index <- na.omit(Leaf_area_index)

# Divide the values in the Leaf_area_index column by 10
Leaf_area_index$Leaf_area_index <- Leaf_area_index$Leaf_area_index / 10

# Display the final data frame
print(Leaf_area_index)

# Initialize lists for results
bayesianlist <- list()
ntot_model_list <- list()
scorefix <- fread("remote_sensing_updated_scorefix.csv") %>%
  select(left, top, bottom, right, score, final_height, dbh, plotID, siteID)

# Loop through sites
for (i in sitevector) {
  site_data <- filter(scorefix, siteID == i, dbh >= 10, dbh <= 50)
  ## the data frame "breakpoints" is made of two columns. The site column has the NEON sites, wile
  ## breakpoint has the x value for xlower as calculated in the script for estimating alpha
  breakpoint <- 10^(filter(breakpoints, site == i) %>% pull(breakpoint))
  site_data <- filter(site_data, dbh >= breakpoint)
  
  if (nrow(site_data) < 25) next
  
  # Bin data according to logscale, as in the forestscaling R package. This is better for size-abundance scaling estimations 
  # than histograms
  binned_data <- forestscaling::logbin(site_data$dbh, site_data$score, n = 8) %>%
    rowwise() %>%
    mutate(avg_score = mean(filter(site_data, dbh >= bin_min & dbh < bin_max)$score, na.rm = TRUE)) %>%
    ungroup() %>%
    drop_na()
  
  # Sample alphas
  alpha_samples <- sample(rstan::extract(results_model[[i]], "alpha")$alpha, 1000, replace = FALSE)
  DBH_max <- max(filter(scorefix, siteID == i)$dbh, na.rm = TRUE)
  breakpoint_norm <- (breakpoint - 10) / (min(DBH_max, 50) - 10)
  LAI_norm <- Leaf_area_index %>% filter(siteID == i) %>% pull(Leaf_area_index)
  
  if (is.na(LAI_norm)) next
  
  # Prepare Stan data
  stan_data <- list(
    K = nrow(binned_data),
    bin_min = binned_data$bin_min,
    bin_max = binned_data$bin_max,
    N_obs = binned_data$bin_count,
    x_min = 10,
    x_max = min(DBH_max, 50),
    LAI_norm = LAI_norm,
    breakpoint_norm = max(0, breakpoint_norm),
    n_alpha_samples = length(alpha_samples),
    alpha_samples = alpha_samples,
    N_tot_prior_mean = 1250,
    N_tot_prior_sd = 625
  )
  
  # Fit the Stan model
  bayesianNtotmod <- sampling(stan_totaltrees_model_compile, 
                              stan_data, chains = 4, warmup = 6000, iter = 9000)
  bayesianlist[[length(bayesianlist) + 1]] <- summarize_draws(bayesianNtotmod) %>% 
    mutate(site = i)%>%
    filter(variable == "N_tot")%>%
    mutate(N_tot_real = full_neon_data%>%
             filter(siteID == i&stemDiameter>=10&stemDiameter<=50)%>%
             filter(plotID %in% unique(site_data$plotID))%>%
             nrow())%>%
    mutate(breakpoint = breakpoint)
  print(paste(i))
  
  ntot_model_list[[i]]<-bayesianNtotmod
}
