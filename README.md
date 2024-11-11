# Recovering Parameters in Truncated Pareto Distributions for Forest DBH Analysis
This repository provides an RMarkdown file and accompanying code for estimating parameters in a truncated Pareto (power law) distribution. The project focuses on analyzing size-abundance distributions in forests based on remote sensing of tree diameters at breast height (DBH), specifically targeting the recovery of the alpha parameter from a subset of the distribution visible via remote sensing.

### Project Overview
Forest DBH distributions often follow a power-law pattern, with larger trees being less frequent than smaller ones. However, remote sensing data typically only captures the largest observable trees due to canopy cover occlusion. This RMarkdown document demonstrates Bayesian modeling techniques in Stan for estimating the entire distribution from this limited subset.

### Repository Contents
STAN_elegant_approach.Rmd: RMarkdown file with code and explanations for each step in the alpha estimation process.

Simulating_Forests.R: Script to simulate forest data for testing the model, using distributions based on evergreen/deciduous species and global biome equations.

![image](https://github.com/user-attachments/assets/7a09657a-17ee-440e-bf7f-1df0bf892c7f)

**Fig 1.** _Example of a small simulated forest, showing the various estimated parameters. Visibility is calculated last (green = invisible, orange = visible), with the green tree invisible to a remote sensing platform due to being covered by a larger, taller tree._

Estimating_Total_Tree_Abundance.R: Script for estimating the total number of trees in a plot from observed data alone. This approach aims to extrapolate from visible trees to an estimate for the entire distribution.

simulatedforests: Folder with 500 simulated forests each from Nearctic coniferous forests and Nearctic mixed forests. The alpha used to simulate each forest is in its own column in the csv file, and each tree is marked as either "visible" or "invisible" to an overhead remote sensing platform.

![image](https://github.com/user-attachments/assets/8699beea-ce4d-4200-888c-c2b328571349)

**Fig. 2.** _Comparison of estimated alpha values for simulated mixed forest data (based on "remote sensing data", derived from code that determines which trees are visible from above) against the true alpha values for the forest._


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
2. Incorporate error in the estimate of alpha into future predictions (uncertainty propagation)
3. Update/spend more time on the estimation of total tree abundance Ntot based on observed tree abundance Nobs and estimated alpha, with uncertainty in alpha incorporated
