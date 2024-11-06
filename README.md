# Recovering Parameters in Truncated Pareto Distributions for Forest DBH Analysis
This repository provides an RMarkdown file and accompanying code for estimating parameters in a truncated Pareto (power law) distribution. This project is designed for analyzing size-abundance distributions in forests based on remote sensing of tree diameters at breast height (DBH). Our aim is to accurately recover alpha given only a subset of the distribution—data from the largest observable trees, which are visible through remote sensing.

### Project Overview
Forest DBH distributions often follow a power-law pattern, where larger trees are less frequent than smaller ones. Estimating the alpha parameter of this distribution provides critical insights into forest structure and ecological dynamics. However, in remote sensing data, only larger trees may be directly visible due to occlusion from canopy cover. This RMarkdown document demonstrates methods for estimating the entire distribution using Bayesian modeling in Stan.

### Methodology
The document details the following workflow:

Simulating DBH Data: Using a truncated Pareto distribution with random alpha values, we simulate data to mimic remotely sensed DBH observations.
Dynamic Truncation Approach: Since higher values of alpha lead to fewer trees in the observable size range, we iteratively adjust the lower truncation point to include additional trees when necessary. This approach improves estimation by maintaining a minimum data sample size.

Bayesian Model Fitting: We use Stan to fit a truncated Pareto distribution model to the observed data, handling both lower and upper truncations (e.g., observable DBHs are restricted to a maximum of 50 cm).
Evaluation of Estimation Accuracy: After fitting the model, we calculate the mean absolute error and root mean squared error for alpha estimates and visualize the relationship between true and estimated values.
Key Functions

Dynamic Adjustment of Lower Truncation: To ensure stable estimation, we dynamically adjust the lower truncation threshold (e.g., trunc_point) if the number of observations is insufficient.
Handling Truncations in Stan: The Stan code models both the lower and upper truncation adjustments, accounting for the limited range of DBH values observable via remote sensing.

### Repository Contents
STAN_elegant_approach.Rmd: Main RMarkdown file with code and explanations for each step of the alpha estimation process.

### Running the Analysis
#### Prerequisites
R and RStudio
rstan package
ggplot2 for visualization
data.table for data handling

### Instructions
Clone this repository to your local machine.
Open analysis.Rmd in RStudio.
Run the entire RMarkdown document to simulate data, fit the model, and evaluate results. The script should dynamically adjust the lower truncation point and plot the results.

### Future Directions
1. Testing the method on real world data
2. Include additional method to estimate the total number of trees that should be present in the plot, based on observed tree counts and the estimated alpha value
