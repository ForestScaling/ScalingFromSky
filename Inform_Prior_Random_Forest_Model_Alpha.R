# ###########################################################
#                Random Forest Model for Predicting Alpha
# ###########################################################
#
# Purpose: 
# This script aims to develop a random forest model to predict the alpha parameter
# for size-abundance distributions in forest plots using existing FIA data. The alpha
# parameter describes the shape of the size-abundance distribution, and each plot has
# already had its alpha estimated. 
# 
# The goal is to use the trained random forest model to predict alpha for future forest
# sites, leveraging features like tree height, location (latitude and longitude), and
# elevation. These predicted alpha values will be used as informed priors in a Bayesian 
# model for estimating alpha across a larger spatial range of forest sites.
#
# Additionally, the script incorporates uncertainty in the predicted alpha values by 
# calculating both the mean and standard deviation of the predictions, which will then
# be incorporated into the priors for future alpha estimation. 
#
# To assess the performance of the model, cross-validation is performed to test how well 
# the random forest generalizes to unseen data. This will provide a robust estimate of 
# model accuracy and help in understanding the variability in alpha predictions.


# Load necessary libraries
library(dplyr)
library(tidyr)
library(data.table)
library(ranger)
library(caTools)
library(spm)

# Load the datasets
dataframe <- fread("C:\\Users\\adam.local\\Downloads\\xmin12.7cm_xmax50cm_recentclip.csv")
subplotonly <- fread("C:\\Users\\adam.local\\Downloads\\slope_xmin12.7cm_xmax50cm_alphatest.csv")

# Prepare the FIA data by joining and transforming the dataset
FIAData <- subplotonly %>%
  inner_join(dataframe) %>%
  separate(geometry, into = c("Lon", "Lat"), sep = "\\|") %>%
  mutate(Lat = as.numeric(Lat), Lon = as.numeric(Lon)) %>%
  mutate(Max_Height_meters = MAX_HEIGHT * 0.3048)

# Train the random forest model to predict alpha based on selected features
rangermodel <- ranger(mean ~ Max_Height_meters + Lon + Lat + ELEV, 
                      num.threads = 11,
                      data = FIAData %>%
                        drop_na(), 
                      keep.inbag = TRUE,
                      write.forest = TRUE, 
                      verbose = TRUE)

# Perform cross-validation using rgcv function
cvranger <- rgcv(trainx = FIAData %>%
                   select(Max_Height_meters, Lon, Lat, ELEV), 
                 trainy = FIAData %>%
                   pull(mean), 
                 verbose = TRUE)

# Predict alpha for future forest sites using the random forest model
neonpredict <- predict(rangermodel, 
                       data = (remote_sensing_updated %>%
                                 distinct(siteID, .keep_all = TRUE) %>%
                                 mutate(Max_Height_meters = final_height) %>%
                                 inner_join(vst$vst_perplotperyear %>%
                                              drop_na(decimalLongitude, decimalLatitude) %>%
                                              distinct(siteID, .keep_all = TRUE) %>%
                                              mutate(Lat = decimalLatitude, Lon = decimalLongitude, ELEV = elevation) %>%
                                              select(Lat, Lon, ELEV, siteID)) %>%
                                 data.frame()), 
                       predict.all = TRUE)

# Assuming 'neonpredict$predictions' is a matrix of predictions from the ranger model
# Calculate mean and standard deviation for each data point (each row)
prior_mean <- apply(neonpredict$predictions, 1, mean)  # Mean for each row
prior_sd <- apply(neonpredict$predictions, 1, sd)      # Standard deviation for each row

# Prepare the final output with predicted means and standard deviations
neonregressionoutput <- (remote_sensing_updated %>%
                           distinct(siteID, .keep_all = TRUE) %>%
                           mutate(Max_Height_meters = final_height) %>%
                           inner_join(vst$vst_perplotperyear %>%
                                        drop_na(decimalLongitude, decimalLatitude) %>%
                                        distinct(siteID, .keep_all = TRUE) %>%
                                        mutate(Lat = decimalLatitude, Lon = decimalLongitude, ELEV = elevation) %>%
                                        select(Lat, Lon, ELEV, siteID))) %>%
  select(siteID) %>%
  cbind(prior_mean) %>%
  cbind(prior_sd)

# The output 'neonregressionoutput' will contain siteID along with the predicted
# prior mean and standard deviation for alpha, which can be used in the Bayesian
# model as informed priors.