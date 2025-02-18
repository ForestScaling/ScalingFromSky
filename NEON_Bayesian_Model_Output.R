# ###########################################################
#              Testing Abundance Output for Field Data Replication
# ###########################################################
#
# Purpose:
# This script tests whether the abundance estimates obtained from the posterior distributions
# of the remote sensing model replicate field data. It interpolates remote sensing-based abundance
# estimates and compares them with field data by fitting a mixed-effects model with `brms`.
#
# The script:
# 1. **Reads and combines data**: It loads and combines the abundance data from the "NEON_combined_results_prep_brmmultiple" directory.
# 2. **Data preparation**: The script performs log transformation of abundance and calculates summary statistics (mean and SD).
# 3. **Fits a Bayesian model**: It fits a mixed-effects model to the abundance data using the `brms` package.
# 4. **Saves the model**: Finally, it saves the fitted model for further analysis.

# ###########################################################

# Load necessary libraries
library(dplyr)        # For data manipulation and summary statistics
library(data.table)   # For efficient data manipulation
library(brms)         # For Bayesian mixed-effects models
library(posterior)    # For working with posterior distributions
library(purrr)        # For functional programming (map functions)

# Step 1: Load and combine data
# Read the first 500 files from the "NEON_combined_results_prep_brmmultiple" directory
data <- list.files("NEON_combined_results_prep_brmmultiple", full.names = TRUE)[1:500] %>%
  lapply(fread)  # Use fread from data.table to read each file into a list

# Combine the list of data frames into one large data frame
data <- data %>%
  rbindlist(idcol = TRUE)  # Add a column to identify the source file (idcol = TRUE)

# Step 2: Data transformation and summarization
# Group data by size class, method, and site
data <- data %>%
  group_by(size_class, method, site) %>%
  mutate(abundance = log10(abundance)) %>%  # Log transform the abundance for better model behavior
  summarize(
    abundance_mean = mean(abundance),  # Calculate the mean abundance for each group
    abundance_sd = sd(abundance)      # Calculate the standard deviation of abundance for each group
  ) %>%
  ungroup()  # Remove grouping after summarization

# Step 3: Fit a Bayesian mixed-effects model
# The model predicts log-transformed abundance mean using size class and method as predictors
# It includes a random effect for site, with a method-specific slope for size class
neon_randommodel <- brms::brm(
  formula = abundance_mean | se(abundance_sd, sigma = TRUE) ~ log10(size_class) * method +  
    (1 + log10(size_class) | site),  # Random intercept and slope for size class within each site
  
  data    = data,  # Data for the model
  warmup  = 6000,  # Number of warmup iterations
  iter    = 9000,  # Total number of iterations
  chains  = 4,     # Number of MCMC chains
  cores   = 4,     # Number of cores for parallel computation
  threads = threading(60)  # Number of threads for parallel processing within a chain (adjust based on your machine and whether you need it)
)

# Step 4: Save the fitted model for future analysis
save(neon_randommodel, file = "neon_randommodel.RData")

