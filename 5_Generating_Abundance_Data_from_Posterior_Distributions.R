# ###########################################################
#              Generating Abundance Data from Posterior Distributions
# ###########################################################
#
# Purpose:
# This script is designed to generate posterior distributions for potential tree abundances
# across a range of DBH values for different size classes, based on alpha and Ntot parameters 
# estimated from the posterior distribution of a Bayesian model. 
# 
# The script processes data for multiple forest sites, where both the alpha (shape of the 
# size-abundance distribution) and Ntot (total tree count) have already been estimated 
# from posterior samples for each site.
#
# For each site, the script generates abundance data in two ways:
# 1. **Indirect method**: Using posterior samples for both alpha and Ntot, the script calculates
#    the expected abundance at different DBH size classes. The `generate_abundance_data()` function
#    uses the Pareto distribution (via `VGAM::dpareto`) to calculate density for each size class, 
#    then multiplies this by the total number of trees (Ntot) to estimate abundance.
#
# 2. **Direct method**: Using the real (observed) value of alpha for the site, along with a fixed 
#    value of Ntot (calculated directly from the field data), the script generates abundance 
#    estimates again based on the same DBH size classes. This serves as a comparison to the indirect method.
#
# The results are combined, and each posterior realization is saved for further analysis. The 
# generated data will be used for assessing the reliability and variability of the predicted 
# abundance at different size classes and for different sites.
#
# The posterior samples for each site are drawn from Stan model results, and the final output 
# is stored in `.rds` format for further processing.
#
# ###########################################################

generate_abundance_data <- function(alpha_posterior, ntot_posterior, dbh_values, site_label, method_label, xmin = 10) {
  # Generate all combinations of alpha and ntot
  posterior_combinations <- expand.grid(
    alpha = alpha_posterior,
    ntot = ntot_posterior
  )
  
  # Repeat dbh values for each posterior combination
  dbh_rep <- rep(dbh_values, times = nrow(posterior_combinations))
  
  # Repeat alpha and ntot values for each dbh value
  alpha_rep <- rep(posterior_combinations$alpha, each = length(dbh_values))
  ntot_rep <- rep(posterior_combinations$ntot, each = length(dbh_values))
  
  # Calculate Pareto density for all combinations
  densPareto <- VGAM::dpareto(x = dbh_rep, shape = alpha_rep, scale = xmin)
  
  # Calculate abundance for all combinations
  abundance <- ntot_rep * densPareto
  
  # Create the final data frame
  abundance_data <- data.frame(
    size_class = dbh_rep,
    abundance = abundance,
    alpha = alpha_rep,
    ntot = ntot_rep,
    site = site_label,
    method = method_label
  )
  
  return(abundance_data)
}

# Example DBH size classes
dbh_values <- seq(10, 50, length.out=100)  

# Loop through each site and generate the posterior abundance data
for(site in sitevector){
  
  # Extract posterior samples for alpha and Ntot
  alpha_posterior <- as.vector(rstan::extract(results_model[[site]], "alpha")$alpha)
  ntot_posterior <- as.vector(rstan::extract(ntot_model_list[[site]], "N_tot")$N_tot)
  
  # Generate abundance data for the indirect method (using posterior for alpha and ntot)
  posterior_datasets1 <- lapply(1:12000, function(i) {
    abundance_data_indirect <- generate_abundance_data(
      alpha_posterior = alpha_posterior[i],
      ntot_posterior = ntot_posterior[i],
      dbh_values = dbh_values,
      site_label = site,
      method_label = "Remote"  # Label for indirect method
    )
    
    # Get real alpha and Ntot for the direct method
    alpha_posterior_real <- as.vector(rstan::extract(realalpha[[site]], "alpha")$alpha)
    ntot_posterior_real <- full_neon_data %>%
      filter(siteID == site) %>%
      filter(stemDiameter >= 10 & stemDiameter <= 50) %>%
      nrow()
    
    # Generate abundance data for the direct method (using fixed posterior for ntot_real)
    abundance_data_direct <- generate_abundance_data(
      alpha_posterior = alpha_posterior_real[i],
      ntot_posterior = ntot_posterior_real,  # Fixed ntot for direct method
      dbh_values = dbh_values,
      site_label = site,
      method_label = "Field"  # Label for direct method
    )
    
    # Combine both indirect and direct method data for this posterior realization
    combined_data <- rbind(abundance_data_indirect, abundance_data_direct)
    combined_data$posterior_index <- site  # Add index to track posterior realization
    
    return(combined_data)
  })
  
  # Perform garbage collection to manage memory usage
  gc()
  
  # Save the posterior datasets to an .rds file for future analysis
  saveRDS(posterior_datasets1, file = paste0("posterior_NEONpareto//", site, ".rds"))
  
  # Remove temporary variables to free up memory
  rm(posterior_datasets1)
  gc()
  
  # Print progress message with the current site
  print(paste(site))
}
