
# Stan Models for Tree Density and Truncation Analysis

This folder contains Stan models designed to estimate tree density and shape parameters for Pareto distributions using observational data. The models are specifically designed to handle truncated distributions and account for observational biases. There are two key models:

## 1. Pareto Distribution with Adjustment Factor (Model 1)

### Purpose:
This model estimates the shape parameter (`alpha`) of a truncated Pareto distribution. The model adjusts for observational biases using factors like normalized Leaf Area Index (LAI) and a breakpoint distance, which represent environmental variables affecting tree visibility.

### Inputs:
- `N`: Number of observations.
- `x_min`: Minimum threshold for the Pareto distribution.
- `trunc_point`: Lower truncation threshold (e.g., 20 cm DBH).
- `trunc_upper`: Upper truncation threshold (e.g., 50 cm DBH).
- `x`: Observed tree DBH values, truncated between `trunc_point` and `trunc_upper`.
- `LAI_norm`: Normalized value of Leaf Area Index (0 to 1).
- `breakpoint_norm`: Normalized value of breakpoint distance (0 to 1).
- `prior_mean`: Mean of the truncated normal prior for the shape parameter (`alpha`).
- `prior_sd`: Standard deviation of the truncated normal prior for `alpha`.

### Outputs:
- The estimated shape parameter (`alpha`) of the Pareto distribution.

### Model Details:
1. **Prior for `alpha`:** The model assumes that `alpha` follows a truncated normal distribution with a mean (`prior_mean`) and standard deviation (`prior_sd`), constrained to be greater than 0.
2. **Likelihood:** The likelihood is based on the truncated Pareto cumulative distribution function (CDF) within the range `[trunc_point, trunc_upper]`. The model adjusts the likelihood by a factor derived from the LAI and breakpoint distance, which are assumed to penalize the likelihood based on observational biases.
3. **Adjustment Factor:** An adjustment factor is computed to correct for biases, where a lower LAI and breakpoint distance lead to a stronger penalty in the likelihood.

## 2. Tree Density Estimation with Pre-sampled `alpha` Values (Model 2)

### Purpose:
This model estimates the total number of trees (`N_tot`) with DBH greater than or equal to `x_min` based on observed counts within defined DBH bins. It uses pre-sampled `alpha` values to marginalize over the uncertainty of the shape parameter in the Pareto distribution and adjusts for observational biases.

### Inputs:
- `K`: Number of DBH bins.
- `bin_min[K]`: Lower bounds of DBH for each bin.
- `bin_max[K]`: Upper bounds of DBH for each bin.
- `N_obs[K]`: Observed counts of trees per bin (continuous values).
- `x_min`: Minimum DBH threshold.
- `x_max`: Maximum DBH threshold.
- `n_alpha_samples`: Number of pre-sampled `alpha` values.
- `alpha_samples`: Array of pre-sampled `alpha` values.
- `N_tot_prior_mean`: Mean of the prior for `N_tot`.
- `N_tot_prior_sd`: Standard deviation of the prior for `N_tot`.
- `LAI_norm`: Normalized Leaf Area Index value (0 to 1).
- `breakpoint_norm`: Normalized breakpoint value (0 to 1).

### Outputs:
- The estimated total number of trees with DBH greater than or equal to `x_min` (`N_tot`).

### Model Details:
1. **Prior for `N_tot`:** The model assumes that the total number of trees (`N_tot`) follows a normal prior with mean (`N_tot_prior_mean`) and standard deviation (`N_tot_prior_sd`).
2. **Likelihood:** The likelihood is based on the observed counts of trees (`N_obs`) within predefined DBH bins, modeled using a normal distribution. The model computes an expected count (`lambda[k]`) for each bin using the pre-sampled `alpha` values.
3. **Adjustment Factor:** Similar to the first model, an adjustment factor is used to account for biases based on LAI and breakpoint values. This adjustment factor modifies the estimate of `N_tot`.
4. **Marginalization over `alpha`:** The model marginalizes over the uncertainty of the `alpha` parameter by using pre-sampled `alpha` values. For each sample, the model computes the expected tree counts in each bin and then calculates the log-likelihood of the observed counts.

## How These Models Work Together

- **Model 1** estimates the shape parameter (`alpha`) of the Pareto distribution, which governs the distribution of tree sizes (DBH). This model is designed to account for truncation in the data and biases due to observational factors like LAI and breakpoint distance.
- **Model 2** estimates the total number of trees (`N_tot`) with DBH greater than or equal to `x_min`, using the pre-sampled `alpha` values from Model 1. This model accounts for uncertainties in `alpha` by marginalizing over multiple samples and also adjusts for observational biases.

Together, these models allow for a comprehensive understanding of tree density and distribution, incorporating both environmental factors and sampling biases.

## Running the Models

1. **Pre-sample `alpha` values:** For Model 2, you will need to pre-sample `alpha` values from the posterior of Model 1 or another source.
2. **Data Preparation:** Ensure that your data is formatted correctly, with appropriate truncation points and observed tree counts.
3. **Stan Compilation:** To run the models, compile the Stan models in your R environment using the `rstan` package:

    ```r
    model1 <- stan_model("model1.stan")
    model2 <- stan_model("model2.stan")
    ```

4. **Sampling:** Once the models are compiled, run the sampling process to estimate the posterior distributions of the parameters:

    ```r
    fit1 <- sampling(model1, data = data_model1)
    fit2 <- sampling(model2, data = data_model2)
    ```

## Conclusion

These Stan models provide a framework for estimating tree density and distribution in forest ecosystems, accounting for biases due to sampling and environmental factors. By combining truncation modeling with adjustments for observational conditions, they offer a robust way to estimate tree abundance and Pareto distribution parameters.