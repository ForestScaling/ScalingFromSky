# -------------------------------------------------------------------------
# Preparing NEON Data for Analysis: Removing Multibole Trees
# -------------------------------------------------------------------------
#
# Purpose:
# The NEON vegetation structure dataset provides valuable information on 
# tree abundance and size distributions. However, one challenge in the 
# dataset is the presence of **multibole trees**, where trees with multiple 
# boles (stems) are recorded as separate entries under different 
# `individualID`s. Additionally, for this analysis, we focus on trees with 
# a **DBH ≥ 10 cm** since smaller trees (DBH < 10 cm) are measured across 
# variable nested plot sizes that are unmapped in space, making their 
# abundance measurements less comparable across plots.
#
# This script prepares the dataset for further analysis by removing multibole 
# entries and excluding trees with DBH < 10 cm, ensuring that only trees with 
# consistent and comparable measurements remain. This is essential for accurate 
# estimation of tree abundance and size distributions.
#
# Objectives:
# 1. Load the NEON dataset and filter it based on selected plots.
# 2. Remove multibole entries from the dataset to ensure each tree is unique.
# 3. Exclude trees with DBH < 10 cm to ensure the data is comparable across plots.
# 4. Prepare the dataset for subsequent analyses by retaining the most recent 
#    measurements for each tree.

# -------------------------------------------------------------------------
# Libraries
# -------------------------------------------------------------------------
library(geoNEON)
library(rlang)
library(data.table)
library(dplyr)
library(neonUtilities)
library(stringr)

# -------------------------------------------------------------------------
# Data Loading
# -------------------------------------------------------------------------
# Load the NEON vegetation structure data
load("stackedNEONtables.Rdata")

# Filter the data to include only the necessary plots
plotd <- remote_sensing_updated %>% distinct(plotID) %>% pull(plotID)

# Subset the NEON data for the selected plots, focusing on DBH >= 10 cm
full_neon_data <- vst$vst_apparentindividual %>%
  filter(plotID %in% plotd, stemDiameter >= 10) %>%  # Exclude DBH < 10 cm
  group_by(siteID, plotID, subplotID) %>%
  distinct(individualID, .keep_all = TRUE) %>%
  ungroup()

# -------------------------------------------------------------------------
# Data Preparation
# -------------------------------------------------------------------------
# Remove multibole trees from the dataset
remove_multibole <- function(df) {
  # Extract the last section and siteID from individualID
  id_split <- str_split_fixed(df$individualID, "\\.", 5)
  last_id_section <- id_split[, 5]
  siteID <- id_split[, 4]
  
  # Identify whether the last section is a bole (contains letters)
  last_digits <- gsub("[^0-9]", "", last_id_section)
  is_bole <- last_id_section != last_digits
  
  # Combine relevant columns into a lookup table
  id_lut <- df %>%
    dplyr::mutate(
      last_digits = last_digits,
      last_id_section = last_id_section,
      siteID = siteID,
      is_bole = is_bole
    )
  
  # Identify IDs with multiple entries
  multiple_ids <- id_lut %>%
    dplyr::count(last_digits, siteID) %>%
    filter(n > 1) %>%
    select(last_digits, siteID)
  
  # Filter out multibole entries
  remove_ids <- id_lut %>%
    inner_join(multiple_ids, by = c("last_digits", "siteID")) %>%
    filter(is_bole) %>%
    pull(individualID)
  
  # Return dataframe without multibole entries
  df %>%
    filter(!individualID %in% remove_ids)
}

# Step 3: Remove multibole entries and retain the most recent measurements
full_neon_data <- full_neon_data %>%
  remove_multibole()
