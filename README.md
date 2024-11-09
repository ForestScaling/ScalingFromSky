# Recovering Parameters in Truncated Pareto Distributions for Forest DBH Analysis
This repository provides an RMarkdown file and accompanying code for estimating parameters in a truncated Pareto (power law) distribution. The project focuses on analyzing size-abundance distributions in forests based on remote sensing of tree diameters at breast height (DBH), specifically targeting the recovery of the alpha parameter from a subset of the distribution visible via remote sensing.

### Project Overview
Forest DBH distributions often follow a power-law pattern, with larger trees being less frequent than smaller ones. However, remote sensing data typically only captures the largest observable trees due to canopy cover occlusion. This RMarkdown document demonstrates Bayesian modeling techniques in Stan for estimating the entire distribution from this limited subset.

### Repository Contents
STAN_elegant_approach.Rmd: Main RMarkdown file with code and explanations for each step in the alpha estimation process.
Simulating_Forests.R: Script to simulate forest data for testing the model, using distributions based on evergreen/deciduous species and global biome equations.
Estimating_Total_Tree_Abundance.R: Script for estimating the total number of trees in a plot from observed data alone. This approach aims to extrapolate from visible trees to an estimate for the entire distribution.

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
2. Adjust the forest simulation code. Right now it's based on evergreen and deciduous for one equation and the global biome equation for the other. If we switch to angiosperm and gymnosperm, we can be more specific with both equations.
