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

# Download and load NEON vegetation structure data
# vst <- loadByProduct(dpID = "DP1.10098.001", check.size = FALSE)

# Save the data for future use
# save(vst, file = "stackedNEONtables.Rdata")

# Load the NEON vegetation structure data
load("stackedNEONtables.Rdata")

# Filter the data to include only the necessary plots
plotd <- c("ABBY_002", "ABBY_009", "ABBY_004", "ABBY_016", "ABBY_067", "ABBY_065", "ABBY_061", "ABBY_014", "ABBY_006", "ABBY_019", "ABBY_063", "ABBY_069", "ABBY_073", "ABBY_074", "ABBY_064", "ABBY_070", "ABBY_075", "ABBY_068", "ABBY_062", "ABBY_066", "ABBY_025", "ABBY_076", 
                    "BART_015", "BART_013", "BART_023", "BART_010", "BART_016", "BART_075", "BART_011", "BART_007", "BART_001", "BART_019", "BART_002", "BART_024", "BART_003", "BART_027", "BART_047", "BART_071", "BART_034", "BART_072", "BART_032", "BART_046", "BART_050", "BART_041", "BART_040", "BART_042", 
                    "BART_039", "BART_073", "BART_051", "BART_036", "BART_012", "BART_025", "BART_026", "BART_006", "BART_018", "BART_005", "BART_070", "BART_044", "BART_037", "BART_074", "BART_033", "BART_004", "BART_028", "BLAN_001", "BLAN_009", "BLAN_015", "BLAN_044", "BLAN_050", "BLAN_047", "BLAN_048", 
                    "BLAN_063", "BLAN_052", "BLAN_036", "BLAN_032", "BLAN_035", "BLAN_041", "BLAN_011", "BLAN_005", "BLAN_017", "BLAN_007", "BLAN_018", "BLAN_012", "BLAN_013", "BONA_012", "BONA_011", "BONA_018", "BONA_008", "BONA_089", "BONA_078", "BONA_076", "BONA_083", "BONA_093", "BONA_071", "BONA_077", "BONA_080", 
                    "BONA_081", "BONA_086", "BONA_070", "BONA_016", "BONA_075", "BONA_079", "BONA_087", "BONA_088", "BONA_072", "BONA_082", "BONA_073", "BONA_006", "BONA_001", "BONA_026", "BONA_009", "BONA_007", "BONA_021", "BONA_005", "BONA_023", "BONA_084", "BONA_013", "BONA_014", "BONA_003", "BONA_002", "BONA_004", 
                    "CLBJ_032", "CLBJ_043", "CLBJ_040", "CLBJ_051", "CLBJ_042", "CLBJ_039", "CLBJ_047", "CLBJ_002", "CLBJ_007", "CLBJ_008", "CLBJ_036", "CLBJ_015", "CLBJ_028", "CLBJ_005", "CLBJ_011", "CLBJ_003", "CLBJ_021", "CLBJ_049", "CLBJ_018", "CLBJ_016", "CLBJ_009", "CLBJ_022", "CLBJ_017", "CLBJ_014", "CLBJ_033", 
                    "CLBJ_041", "CLBJ_038", "CLBJ_046", "CLBJ_052", "CLBJ_044", "CLBJ_025", "CLBJ_006", "CLBJ_010", "CLBJ_001", "CLBJ_004", "CLBJ_023", "CLBJ_030", "CLBJ_026", "CLBJ_045", "CLBJ_012", "CLBJ_027", "CLBJ_029", "CLBJ_019", "CLBJ_024", "DEJU_061", "DEJU_049", "DEJU_047", "DEJU_048", "DEJU_051", "DEJU_060", 
                    "DEJU_062", "DEJU_056", "DEJU_045", "DEJU_050", "DEJU_057", "DEJU_052", "DEJU_014", "DEJU_006", "DEJU_018", "DEJU_055", "DEJU_063", "DEJU_046", "DEJU_015", "DEJU_002", "DEJU_058", "DEJU_017", "DEJU_009", "DEJU_016", "DEJU_021", "DEJU_023", "DEJU_020", "DEJU_024", "DEJU_010", "DEJU_054", "DEJU_053", 
                    "DEJU_064", "DEJU_059", "DEJU_019", "DELA_013", "DELA_008", "DELA_022", "DELA_023", "DELA_003", "DELA_015", "DELA_004", "DELA_002", "DELA_012", "DELA_021", "DELA_014", "DELA_009", "DELA_001", "DELA_006", "DELA_017", "DELA_041", "DELA_044", "DELA_047", "DELA_048", "DELA_018", "DELA_053", "DELA_052", 
                    "DELA_040", "DELA_056", "DELA_045", "DELA_059", "DELA_055", "DELA_039", "DELA_037", "DELA_051", "DELA_020", "DELA_049", "DELA_038", "DELA_016", "DELA_043", "DELA_050", "DELA_046", "DELA_042", "DELA_054", "DELA_011", "DELA_005", "DSNY_001", "DSNY_008", "DSNY_014", "DSNY_010", "DSNY_007", "DSNY_006", 
                    "DSNY_013", "DSNY_034", "DSNY_018", "DSNY_011", "DSNY_025", "DSNY_003", "DSNY_012", "GUAN_056", "GUAN_042", "GUAN_048", "GUAN_015", "GUAN_059", "GUAN_007", "GUAN_003", "GUAN_050", "GUAN_051", "GUAN_057", "GUAN_053", "GUAN_018", "GUAN_061", "GUAN_060", "GUAN_049", "GUAN_001", "GUAN_004", "GUAN_017", 
                    "GUAN_008", "GUAN_005", "GUAN_006", "GUAN_016", "GUAN_019", "GUAN_044", "GUAN_052", "GUAN_047", "GUAN_055", "GUAN_054", "GUAN_045", "GUAN_043", "GUAN_058", "GUAN_046", "GUAN_011", "GUAN_010", "GUAN_002", "GUAN_009", "GUAN_012", "HARV_026", "HARV_027", "HARV_006", "HARV_016", "HARV_011", "HARV_052", 
                    "HARV_048", "HARV_044", "HARV_046", "HARV_035", "HARV_040", "HARV_047", "HARV_033", "HARV_041", "HARV_043", "HARV_038", "HARV_049", "HARV_037", "HARV_034", "HARV_050", "HARV_015", "HARV_020", "HARV_008", "HARV_005", "HARV_017", "HARV_025", "HARV_021", "HARV_012", "HARV_013", "HARV_023", "HARV_001", 
                    "HARV_010", "HARV_014", "HARV_024", "HARV_051", "HARV_042", "HARV_039", "HARV_002", "HARV_022", "HARV_004", "HARV_036", "HARV_045", "HEAL_068", "HEAL_059", "HEAL_016", "HEAL_026", "HEAL_047", "HEAL_025", "HEAL_014", "HEAL_004", "JERC_016", "JERC_010", "JERC_005", "JERC_004", "JERC_002", "JERC_011", 
                    "JERC_001", "JERC_009", "JERC_013", "JERC_007", "JERC_015", "JERC_014", "JERC_003", "JERC_006", "JERC_008", "LODE_005", "LODE_009", "LODE_017", "LODE_010", "LODE_002", "LODE_011", "LODE_003", "LODE_016", "LODE_006", "LODE_018", "LODE_012", "LODE_004", "LODE_015", "LODE_014", "LODE_013", "LODE_008", 
                    "LODE_001", "MUKR_005", "MUKR_006", "MUKR_004", "MUKR_010", "MUKR_008", "MUKR_009", "MUKR_007", "MUKR_011", "MUKR_002", "MUKR_003", "MUKR_001", "MUKR_012", "NIST_017", "NIST_009", "NIST_019", "NIST_004", "NIST_018", "NIST_002", "NIST_016", "NIST_014", "NIST_007", "NIST_013", "NIST_003", "NIST_012", 
                    "NIST_015", "NIST_005", "NIST_008", "NIST_010", "NIST_011", "NIST_006", "NIST_001", "SERC_051", "SERC_050", "SERC_052", "SERC_056", "SERC_033", "SERC_037", "SERC_038", "SERC_046", "SERC_054", "SERC_047", "SERC_060", "SERC_055", "SERC_045", "SERC_048", "SERC_061", "SERC_039", "SERC_053", "SERC_049", 
                    "SERC_059", "SERC_041", "SERC_040", "SERC_044", "SERC_036", "SERC_034", "SERC_032", "SERC_058", "SERC_035", "SERC_031", "SERC_030", "SERC_028", "SERC_029", "SERC_027", "SERC_026", "SERC_025")


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
