[![DOI](https://zenodo.org/badge/884536370.svg)](https://doi.org/10.5281/zenodo.15593238)

# ScalingFromSky

## Overview


This repository contains code and data for the project *Combining Theory with Remote Sensing to Improve Forest Predictions Across Space*. The goal of this work is to combine remote sensing with size-abundance distribution theory to infer the full forest structure from canopy observations.

### Background and Approach

Remote sensing (RS) provides large-scale measurements of forest canopies but often underrepresents smaller, understory trees. However, forests typically follow a well-established **size-abundance relationship**: tree abundance declines with size, often following a power law. This relationship allows us to infer the number of smaller trees from observations of larger canopy trees.

Our approach leverages this idea:
- **If we observe the canopy**, we can estimate the underlying **size-abundance distribution**.
- **The shape of this distribution** is governed by a parameter **alpha**, which determines how quickly tree abundance decreases with size.
- **The total abundance of trees (N_tot)** is also critical: it scales the size-abundance curve from relative density to actual tree counts.

**Note**:  We focus on using a Pareto Distribution in this Github code. However, the approach will still work if a Weibull is used for the interpolation process instead and larger trees are included, as *x_breakpoint* (see paper) occurs at intermediate tree DBH sizes. There will just be more than one breakpoint from the regression, from which *x_breakpoint* must be selected.

Combining RS canopy data with ecological theory allows us to estimate both **alpha** and **N_tot**. With these two parameters, we can reconstruct the size-abundance curve for a plot, extending beyond what is visible from the canopy alone.

### Key Outputs:
- Estimates of **alpha** (size-abundance slope) from remote sensing data.
- Estimates of **N_tot** (total number of trees in a plot).
- Reconstructed size-abundance distributions spanning from small understory trees to large canopy trees.
## What's in this repository?

### Models for Ntot and alpha
This folder contains the STAN codes for the two models we developed to estimate `alpha` and `N_tot`. It also contains a `ReadMe` file explaining how the models work.

### NEON and ForestGEO Data Preparation and Analysis Pipeline
This set of R scripts processes NEON forest inventory data, remote sensing data, and model outputs to estimate total tree abundance (`N_tot`) and the size-abundance scaling parameter (`alpha`). It also includes applications to ForestGEO sites. The pipeline includes the following scripts:

1. **`1_Preparing_NEON_Data_for_Analysis.R`**  
   Cleans and prepares NEON vegetation structure data for size-abundance analysis, filtering and summarizing tree diameter data (DBH) by site.

2. **`2_Inform_Prior_Random_Forest_Model_Alpha.R`**  
   Fits a random forest model to **FIA (Forest Inventory and Analysis)** data to predict `alpha` based on forest structure and environmental variables. This model can then be used to estimate an informed prior for `alpha` at any site across the country based on site-specific input data.

3. **`3_Estimating_Alpha_from_NEON_Data.R`**  
   Uses a Bayesian model (Stan) to estimate the power-law parameter `alpha` from NEON DBH distributions, accounting for truncation due to canopy occlusion.

4. **`4_Estimating_Total_Tree_Abundance_from_Remote_Sensing_Data.R`**  
   Combines remote sensing estimates with prior information to estimate total tree abundance (`N_tot`) using a Bayesian framework.

5. **`5_Generating_Abundance_Data_from_Posterior_Distributions.R`**  
   Simulates abundance estimates from the posterior distributions of the `N_tot` and `alpha` models, generating a distribution of possible tree abundance values for each site.

6. **`6_Processing_Posterior_Data_for_Abundance.R`**  
   Aggregates and processes the abundance estimates generated from `5_Generating_Abundance_Data_from_Posterior_Distributions.R` across multiple posterior draws.

7. **`7_NEON_Bayesian_Model_Output.R`**  
   Fits a Bayesian hierarchical model to compare estimated tree abundances across different methods (e.g., remote sensing interpolation vs. field inventory) and size classes.

8. **`8_ApplyMethod_HARVForestGeo.R`**  
   Applies the abundance estimation pipeline to the **Harvard Forest ForestGEO** plot, using remote sensing data and inventory data from the site.

9. **`9_ApplyMethod_SERCForestGeo.R`**  
   Applies the abundance estimation pipeline to the **SERC ForestGEO** plot, using remote sensing data and inventory data from the site.
