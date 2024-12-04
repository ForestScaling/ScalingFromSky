# Recovering Parameters in Truncated Pareto Distributions for Forest DBH Analysis
This repository provides an RMarkdown file and accompanying code for estimating parameters in a truncated Pareto (power law) distribution. The project focuses on analyzing size-abundance distributions in forests based on remote sensing of tree diameters at breast height (DBH), specifically targeting the recovery of the alpha parameter from a subset of the distribution visible via remote sensing.

### Project Overview
Forest DBH distributions often follow a power-law pattern, with larger trees being less frequent than smaller ones. However, remote sensing data typically only captures the largest observable trees due to canopy cover occlusion. This RMarkdown document demonstrates Bayesian modeling techniques in Stan for estimating the entire distribution from this limited subset.

### Repository Contents
_STAN_elegant_approach.Rmd_: RMarkdown file with code and explanations for each step in the alpha estimation process.

_Simulating_Forests.Rmd_: Script to simulate forest data for testing the model, using distributions based on evergreen/deciduous species and global biome equations.

![image](https://github.com/user-attachments/assets/7a09657a-17ee-440e-bf7f-1df0bf892c7f)

**Fig 1.** _Example of a small simulated forest, showing the various estimated parameters. Visibility is calculated last (green = invisible, orange = visible), with the green tree invisible to a remote sensing platform due to being covered by a larger, taller tree._

_Estimating_Total_Tree_Abundance.Rmd_: Script for estimating the total number of trees in a plot from observed data alone. This approach aims to extrapolate from visible trees to an estimate for the entire distribution.

_simulatedforests_: Folder with 500 simulated forests each from Nearctic coniferous forests and Nearctic mixed forests. The alpha used to simulate each forest is in its own column in the csv file, and each tree is marked as either "visible" or "invisible" to an overhead remote sensing platform.

_Real_world_adjustments_stan.Rmd_: Script that adjusts the model created in STAN_elegant_approach.Rmd to estimate alpha from real-world data. As of now, STAN_elegant_approach appears to works well when the canopy is segmented perfectly, but falters significantly when segmentation is so-so. We include a prior based on the maximum height of trees in the plot, which we know is predictive of alpha from Duncanson et al., 2015. The prior is coded to become more informative when uncertainty in the crown segmentation process is high (low crown scores) and less informative when uncertainty is low (high crown scores). It _should_ be that the prior is not considered when the canopy is known perfectly, but this has not yet been tested.   

![image](https://github.com/user-attachments/assets/8699beea-ce4d-4200-888c-c2b328571349)

**Fig. 2.** _Comparison of estimated alpha values for simulated mixed forest data (based on "remote sensing data", derived from code that determines which trees are visible from above) against the true alpha values for the forest._

_Testing Bayesian Method on NEON Vegetation Structure Data (Tree Abundance)_: Script that tests a Bayesian method for estimating tree abundance using NEON vegetation structure data. The NEON dataset provides detailed measurements of tree stem diameters (DBH) and includes information about the sampling design. A key challenge in working with this dataset is the use of nested subplots to measure trees with DBH <10 cm. These nested subplots may or may not be used on a plot-by-plot basis, depending on the density of smaller trees and other local factors. This variation introduces complexity because plots with nested subplots measure smaller trees over different areas, while trees with DBH ≥10 cm are typically measured across the full plot area. To avoid potential biases and simplify the analysis, we will focus exclusively on plots where nested subplots were *not* used, ensuring that all trees, regardless of size, were measured within the full plot area. This script filters the NEON dataset to retain only plots where nested subplots were not used, sets up a prior for tree abundance in each location using the TreeMap dataset (published online), and then tests a Bayesian approach for estimating total tree abundance (`N_tot`).

![image](https://github.com/user-attachments/assets/0a703426-77cf-4fd8-b80c-07c4221e7401)
**Fig. 3.** _A  plot of the TreeMap dataset, cropped to just the extent of the Harvard Forest NEON site. The raster is displaying estimates of trees per acre on a 30x30m scale. In our scripts, we convert this estimate to trees per m<sup>2</sup> and then multiply by the area of the plot in question in m<sup>2</sup>._

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
