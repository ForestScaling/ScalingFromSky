# Recovering Parameters in Truncated Pareto Distributions for Forest DBH Analysis
This repository provides an RMarkdown file and accompanying code for estimating parameters in a truncated Pareto (power law) distribution. The project focuses on analyzing size-abundance distributions in forests based on remote sensing of tree diameters at breast height (DBH), specifically targeting the recovery of the alpha parameter from a subset of the distribution visible via remote sensing.

### Project Overview
Forest DBH distributions often follow a power-law pattern, with larger trees being less frequent than smaller ones. However, remote sensing data typically only captures the largest observable trees due to canopy cover occlusion. This RMarkdown document demonstrates Bayesian modeling techniques in Stan for estimating the entire distribution from this limited subset.

### Repository Contents
STAN_elegant_approach.Rmd: RMarkdown file with code and explanations for each step in the alpha estimation process.
Simulating_Forests.R: Script to simulate forest data for testing the model, using distributions based on evergreen/deciduous species and global biome equations.
Estimating_Total_Tree_Abundance.R: Script for estimating the total number of trees in a plot from observed data alone. This approach aims to extrapolate from visible trees to an estimate for the entire distribution.
simulatedforests: Folder with 500 simulated forests each from Nearctic coniferous forests and Nearctic mixed forests. The alpha used to simulate each forest is in its own column in the csv file, and each tree is marked as either "visible" or "invisible" to an overhead remote sensing platform.

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
