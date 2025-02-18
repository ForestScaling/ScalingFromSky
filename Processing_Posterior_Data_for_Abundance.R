# ###########################################################
#                    Processing Posterior Data for Abundance
# ###########################################################
#
# Purpose:
# This script processes the posterior abundance data generated in the previous script. 
# The goal is to combine the individual posterior distributions for different sites 
# into a single dataset for each posterior realization, which can be used for further 
# analysis or model fitting.
#
# The script performs the following steps:
# 1. **Loading posterior data**: It loads the results saved as `.rds` files in the "posterior_NEONpareto" directory.
# 2. **Parallelization setup**: It sets up parallel processing to speed up the combination of large datasets. 
# 3. **Processing each element**: For each element in the posterior data lists, it reads the data from all sites and 
#    combines the corresponding entries into a single data frame.
# 4. **Saving results**: The combined data for each posterior element is saved as a `.csv` file in the "NEON_combined_results_prep_brmmultiple" directory.
# 
# The final output is a set of combined `.csv` files, each containing the processed posterior data 
# for one posterior realization, which can then be used for further analysis or model fitting.
#
# ###########################################################

library(data.table)     # For rbindlist() and fwrite()
library(future)         # For parallelization
library(future.apply)   # Parallelized versions of apply functions

# Step 1: Get file paths for the saved lists
file_paths <- list.files("posterior_NEONpareto", full.names = TRUE)

# Step 2: Determine the number of elements in the lists
# Read the first file to get the number of elements
first_list <- readRDS(file_paths[1])
num_elements <- length(first_list)

# Step 3: Create an output directory
output_dir <- "NEON_combined_results_prep_brmmultiple"
dir.create(output_dir, showWarnings = FALSE)

# Step 4: Set up parallelization
plan(multisession, workers = parallel::detectCores())  # Use all but one core

# Step 5: Function to process a single element
process_element <- function(i) {
  # Initialize an empty list to hold data frames for the current element
  current_element_data <- list()
  
  for (path in file_paths) {
    # Read each list
    list_data <- readRDS(path)
    
    # Extract the current element and add it to the list
    current_element_data <- append(current_element_data, list(list_data[[i]]))
  }
  
  # Combine all data frames for the current element
  combined_data <- rbindlist(current_element_data)
  
  # Save the combined data frame to the output directory
  fwrite(combined_data, file.path(output_dir, paste0("combined_element_", i, ".csv")))
  
  # Clear variables to free up memory
  rm(current_element_data, combined_data, list_data)
  gc()
  
  # Return success message
  paste0("Element ", i, " processed and saved.")
}

# Step 6: Parallel loop over elements
results <- future_lapply(seq_len(num_elements), process_element)

# Step 7: Verify output
cat("Combined data frames saved to:", output_dir, "\n")
