
# ScalingFromSky

## Overview


This repository contains code and data for the project *Combining Theory with Remote Sensing to Improve Forest Predictions Across Space*. The goal of this work is to combine remote sensing with size-abundance distribution theory to infer the full forest structure from canopy observations.

### Background and Approach

Remote sensing (RS) provides large-scale measurements of forest canopies but often underrepresents smaller, understory trees. However, forests typically follow a well-established **size-abundance relationship**: tree abundance declines with size, often following a power law. This relationship allows us to infer the number of smaller trees from observations of larger canopy trees.

Our approach leverages this idea:
- **If we observe the canopy**, we can estimate the underlying **size-abundance distribution**.
- **The shape of this distribution** is governed by a parameter **alpha**, which determines how quickly tree abundance decreases with size.
- **The total abundance of trees (N_tot)** is also critical: it scales the size-abundance curve from relative density to actual tree counts.

Combining RS canopy data with ecological theory allows us to estimate both **alpha** and **N_tot**. With these two parameters, we can reconstruct the size-abundance curve for a plot, extending beyond what is visible from the canopy alone.

### Key Outputs:
- Estimates of **alpha** (size-abundance slope) from remote sensing data.
- Estimates of **N_tot** (total number of trees in a plot).
- Reconstructed size-abundance distributions spanning from small understory trees to large canopy trees.

## What's in this repository?