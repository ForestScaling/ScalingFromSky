
library(terra)
library(sf)

sercheight<-terra::rast("C:\\Users\\adam.local\\Downloads\\SERCcanopyheight\\NEON_D02_SERC_DP3_364000_4305000_CHM.tif")
SERCheight<-terra::extract(terra::rast("C:\\Users\\adam.local\\Downloads\\SERCcanopyheight\\NEON_D02_SERC_DP3_364000_4305000_CHM.tif"),
                           read_sf("C:\\Users\\adam.local\\Downloads\\SERCForestGEOPredictions\\SERCForestGEOPredictions.shp"), fun = max)%>%
  
  # Rename the extracted height column
  rename(Height = NEON_D02_SERC_DP3_364000_4305000_CHM) %>%
  
  # Combine with 'bencrown'
  cbind(read_sf("C:\\Users\\adam.local\\Downloads\\SERCForestGEOPredictions\\SERCForestGEOPredictions.shp")) %>%
  
  # Filter out heights that are equal to 0, since that's just missing data
  filter(Height != 0)
rename<-dplyr::rename
SERCheight<-terra::extract(terra::rast("sercheight.tif"),
                           read_sf("C:\\Users\\adam.local\\Downloads\\SERCForestGEOPredictions\\SERCForestGEOPredictions.shp"), fun = max)%>%
  
  # Rename the extracted height column
  rename(Height = b1) %>%
  
  # Combine with 'bencrown'
  cbind(read_sf("C:\\Users\\adam.local\\Downloads\\SERCForestGEOPredictions\\SERCForestGEOPredictions.shp")) %>%
  
  # Filter out heights that are equal to 0, since that's just missing data
  filter(Height != 0)

sercelevation<-terra::rast("NEON_elevation//NEON_D02_SERC_DP3_364000_4305000_DTM.tif")

sercelevation<-project(sercelevation, "WGS84")
# Define bounding box serc
g_bbox <- st_bbox(c(xmin = 364544.2, 
                    ymin = 4305426, 
                    xmax = 364952.8, 
                    ymax = 4305832), crs = 32618)

# Convert bounding box to polygon and transform
g_poly <- st_as_sfc(g_bbox) %>%
  st_transform(32618)


# Create grid within the polygon
g_grid <- st_make_grid(g_poly, units::as_units(10000, "m^2"))#cellsize = 150 / sqrt(2))

# Assign unique IDs to the plots
plot_ids <- paste0("plot_", seq_len(length(g_grid)))

# Convert grid to sf object
plots <- st_as_sf(g_grid) %>%
  mutate(plot_id = plot_ids)
mapview::mapview(plots, legend = FALSE)
library(itcSegment)
serc1<-fread("C:\\Users\\adam.local\\Downloads\\SERC_AllMeasurements\\PlotDataReport04-19-2024_417675680.txt")
# Spatial join between HarvFieldData and plots.

SERCfragmented <- st_join(SERCheight%>%
                            st_sf(), plots, largest = TRUE)
SERCfragmented$perimeter<-as.numeric(st_perimeter(SERCfragmented))
SERCfragmented$area<-as.numeric(st_area(SERCfragmented))
SERCfragmented<-SERCfragmented%>%
  dplyr::mutate(Diameter = 0.5*(sqrt(perimeter^2-(8*area)))) %>%
  # mutate(Diameter1 = sqrt(crown_area))
  # Calculate dbh using formula from Jucker et al. 2017
  mutate(dbh = dbh(H=Height, CA = Diameter, biome = 11))
# mapview::mapview(SERCfragmented, legend = FALSE)
serc1<-serc1%>%
  mutate(X =  364544.2+PX, Y = 4305426 + PY)%>%
  st_as_sf(coords = c("X","Y"),crs = 32618)
sercfieldfragmented <- st_join(serc1, plots, largest = TRUE)%>%
  rename(dbh=DBH)
# mapview::mapview(sercfieldfragmented, legend = FALSE)


data<-list()
model <- stan_model(file = "C:\\Users\\adam.local\\Downloads\\density1_simplified.stan")
library(posterior)
sercrealalpha<-list()
for (i in unique(sercfieldfragmented$plot_id)){
  test<-sercfieldfragmented%>%
    data.frame()%>%
    filter(plot_id == i)%>%
    filter(dbh >= 10 & dbh <= 50)
  # Number of observations
  N <- length(test$dbh)
  # Extract diameters
  x <- test$dbh
  
  # Minimum diameter
  x_min <- 10
  
  # Maximum diameter
  x_max <- max(test$dbh)
  
  # Prepare data for Stan model
  stan_dat <- list(N = N, x = x, x_min = x_min)
  
  # Run Stan model
  fit1 <- sampling(model, data = stan_dat, iter = 9000, warmup = 6000, chains = 4)
  data[[length(data)+1]]<-summarize_draws(fit1)%>%
    mutate(plot_id = i,
           x_min = x_min,
           Ntotal = N)
  sercrealalpha[[i]]<-fit1
}


sercfielddata<-data%>%
  rbindlist()%>%
  filter(variable =="alpha")%>%
  mutate(label = "Field")


fieldserc<-list()


SERCfragmented<-vect(SERCfragmented)%>%
  project("WGS84")
tidy_pred_SERC<-list()
predicted<-list()
newdata1<-list()
plotswgs<-st_centroid(plots%>%st_transform("WGS84"))%>%st_coordinates()%>%
  cbind(plots)
library(tidyr)
for (i in na.omit(unique(sercfieldfragmented$plot_id))){
  SERCfragmentedsubset<-subset(SERCfragmented, sercfieldfragmented$plot_id == i)
  # SERCfragmentedsubset <- SERCfragmented
  newdata <- expand_grid(Max_Height_meters = (SERCfragmented%>%
                                                data.frame()%>%
                                                filter(plot_id == i)%>%
                                                summarize(maxht = max(Height)))$maxht* 3.281,
                         ELEV = terra::extract(sercelevation, vect(ext(SERCfragmentedsubset)), fun= "mean"
                         )[1,2])%>%
    
    mutate(
      Lat= (plotswgs%>%
              filter(plot_id == i))$Y,
      Lon = (plotswgs%>%
               filter(plot_id == i))$X)
  
  newdata1[[length(newdata1)+1]]<-newdata%>%
    mutate(
      plot_id = i)
}

library(ranger)
#NEED TO CHANGE RANGER MODEL FOR DIFFERENT XMIN
# rangermodel<-ranger(mean~Max_Height+LON+LAT+elevation, 
#                     num.threads = 11,
#                     data= train%>%
#                       drop_na()%>%
#                       mutate(Max_Height = 0.3048*Max_Height), keep.inbag = TRUE,
#                     write.forest = TRUE, verbose= TRUE)

sercpredicse<-predict(rangermodel, newdata1%>%
                        rbindlist()%>%
                        mutate(ELEV = ifelse(is.na(ELEV)==TRUE, 5.23, ELEV))
                      # rename(Max_Height=Max_Height_meters,
                      #        LON=Lon,LAT=Lat,elevation = ELEV)
                      , 
                      predict.all=TRUE)
# Calculate mean and standard deviation for each data point (each row)
serc_prior_mean <- apply(sercpredicse$predictions, 1, mean) # Mean for each row
serc_prior_sd <- apply(sercpredicse$predictions, 1, sd)      # Standard deviation for each row

sercregressionoutput<-(newdata1%>%
                         rbindlist())%>%
  cbind(serc_prior_mean)%>%
  cbind(serc_prior_sd)

stan_alpha_laibreak<-"data {
  int<lower=0> N;                  // Number of observations
  real<lower=0> x_min;             // Minimum threshold for the distribution (e.g., 3)
  real<lower=0> trunc_point;       // Observed data starts at this truncation threshold (e.g., 20)
  real<lower=trunc_point> trunc_upper;  // Observed data is capped at this upper threshold (e.g., max DBH for site)
  vector<lower=trunc_point>[N] x;  // Observed data, limited to range [trunc_point, trunc_upper]
  real<lower=0, upper=1> LAI_norm; // Normalized LAI score (0 to 1)
  real<lower=0, upper=1> breakpoint_norm; // Normalized breakpoint distance (0 to 1)
  real prior_mean;                 // Mean of the truncated normal prior for alpha
  real prior_sd;                   // Standard deviation of the truncated normal prior for alpha
}

parameters {
  real<lower=0> alpha;             // Shape parameter for Pareto distribution (must be positive)
}

transformed parameters {
  // Adjustment factor: Penalizes the likelihood based on LAI and breakpoint distance
  real adjustment_factor = 1 - sqrt(LAI_norm * breakpoint_norm);
}

model {
  // Truncated normal prior for alpha, with lower bound at 0
  alpha ~ normal(prior_mean, prior_sd) T[0, ]; // Truncated normal with lower bound at 0
  
  // Truncated cumulative probability in the observed range
  real p_trunc = pareto_cdf(trunc_upper | x_min, alpha) - pareto_cdf(trunc_point | x_min, alpha);
  
  // Likelihood: Adjusted for observational biases
  for (n in 1:N) {
    target += adjustment_factor * (pareto_lpdf(x[n] | x_min, alpha) - log(p_trunc));
  }
}
"

stan_alpha_laibreak_model<-stan_model(model_code=stan_alpha_laibreak)



# results_bootstrap <- data.frame(
#   siteID = character(),
#   Alpha = numeric()
# )

library(splus2R)
library(segmented)
# sitevector<-results_bootstrap$site
results_serc<-data.frame(
  variable = character(),
  site = character(),
  mean=numeric(),
  median=numeric(),
  alpha = numeric(),
  sd = numeric(),
  mad=numeric(),
  q5 = numeric(),
  q95 = numeric(),
  rhat = numeric(),
  rhat = numeric(),
  ess_bulk = numeric(),
  ess_tail = numeric(),
  R2_kernel= numeric(),
  breakpoint =numeric()
)
serc_bayesmodel<-list()
for(site in sercregressionoutput$plot_id){#unique(remote_sensing_updated$siteID)[-6]) {
  site_data <- subset(SERCfragmented%>%
                        data.frame()%>%
                        # select(-geometry)%>%
                        filter(dbh<=50),
                      plot_id == site)
  # site_data <- SERCfragmented%>%
  #                       data.frame()%>%
  #                       # select(-geometry)%>%
  #                       filter(dbh<=50)
  if(nrow(site_data%>%filter(dbh>10))<25){
    next
  }else{
    original_data<-site_data%>%
      pull(dbh)
    
    library(dplyr)
    # Kernel Density Estimation for DBH data
    kde <- density(original_data, bw = "SJ")  # Adjust bandwidth as needed
    xmin<-min(original_data)
    # Truncate the KDE results at x = 50
    kde <- list(
      x = kde$x[kde$x <= max(original_data) & kde$x>=xmin],
      y = kde$y[kde$x <=  max(original_data) & kde$x>=xmin]
    )
    
    # Filter the KDE results to only keep values near the observed data (threshold = 1.5)
    filtered_kde <- data.frame(
      x = kde$x,
      y = kde$y
    ) %>%
      rowwise() %>%
      filter(any(abs(x - original_data) <= 0.5)) %>%
      ungroup()
    # Assume original_data and observed_dbhs are defined
    n_bootstrap <- 1000  # Number of bootstrap samples
    xmin <- min(original_data)
    x_values <- filtered_kde$x  # Unique and sorted x-values from data
    
    # Generate bootstrap KDE
    n_bootstrap <- 1000  # Number of bootstrap samples
    
    
    # Use unique values of the original data as x_values
    x_values <- filtered_kde$x # Unique and sorted x-values from data
    
    # Perform bootstrap sampling and KDE estimation
    bootstrap_kdes <- replicate(
      n_bootstrap,
      {
        sample_data <- sample(original_data, size = length(original_data), replace = TRUE)
        density(sample_data)  # Kernel Density Estimation
      },
      simplify = FALSE
    )
    # Interpolate densities onto x_values
    densities1 <- sapply(bootstrap_kdes, function(kde) approx(kde$x, kde$y, xout = x_values)$y)
    densities <- rowMeans(densities1, na.rm = TRUE)
    # Total observed data (number of trees)
    total_trees <- length(original_data)
    
    # Estimated abundances by multiplying the densities by the total number of trees
    # estimated_abundance <- densities * total_trees
    # Log-transform densities and calculate mean and SE
    mean_log_densities <- log10(densities)
    # mean_log_densities <- rowMeans(log_densities, na.rm = TRUE)
    se_log_densities <- apply(densities1, 1, sd, na.rm = TRUE) / sqrt(n_bootstrap) / 
      apply(densities1, 1, sd, na.rm = TRUE) # Standard error
    
    # Combine results into a data frame
    bootstrap_kde_log <- data.frame(
      log_x = log10(x_values),  # Log10-transform x-values
      mean_log_density = mean_log_densities,
      se_log_density = se_log_densities
    )
    
    
    data<-peaks(bootstrap_kde_log$mean_log_density)%>%
      cbind(bootstrap_kde_log)%>%
      filter(.==TRUE)
    
    # data<-data%>%
    #       filter(mean_log_density>0)
    # Apply 10^x to every numeric column
    
    ##ACCOUNT FOR LOCAL MINIMA
    breakpoint = data%>%
      filter(log_x<=quantile(log10(original_data), 0.75))%>%
      # filter(mean_log_density  == max(mean_log_density ))%>%
      filter(log_x==max(log_x))%>%
      pull(log_x)
    
    segmented_model<-lm(mean_log_density~log_x,bootstrap_kde_log)%>%selgmented()
    
    # Extract the estimated breakpoints
    possiblebreakpoints <- segmented_model$psi[, "Est."]
    
    # Add the minimum and maximum x-values as segment boundaries
    segment_boundaries <- c(min(bootstrap_kde_log$log_x), possiblebreakpoints, max(bootstrap_kde_log$log_x))
    
    # Calculate the slope of each segment
    segment_slopes <- slope(segmented_model)$log_x[, 1]
    
    # Create a data frame with segment details
    segments_df <- data.frame(
      left_x = head(segment_boundaries, -1),     # Left boundary of each segment
      right_x = tail(segment_boundaries, -1),   # Right boundary of each segment
      slope = segment_slopes                    # Slope of each segment
    )
    
    if(segments_df[nrow(segments_df),]$slope<=-10 | segments_df[nrow(segments_df),]$slope>0){
      bayesian_data<-site_data%>%
        filter(dbh>=10^breakpoint)%>%
        filter(dbh>=10)%>%
        filter(dbh<=10^segments_df[nrow(segments_df),]$left_x)
      bootstrap_kde_log<-bootstrap_kde_log%>%
        filter(log_x>=breakpoint)%>%
        filter(log_x>=log10(10))%>%
        filter(log_x<=segments_df[nrow(segments_df),]$left_x)
    }else{
      bayesian_data<-site_data%>%
        filter(dbh>=10^breakpoint)%>%
        filter(dbh>=10)
      
      bootstrap_kde_log<-bootstrap_kde_log%>%
        filter(log_x>=breakpoint)%>%
        filter(log_x>=log10(10))
    }
    
    if(nrow(bayesian_data) == 0){
      print(paste(site, "does not meet the requirements to fit a model"))
    }else{
      # Transform slope prior mean from alpha to slope
      site_priors <- sercregressionoutput %>%
        # mutate(slope_prior_mean = -(prior_mean + 1))%>%
        filter(plot_id==site)
      DBH_max=bayesian_data$dbh%>%max()
      breakpoint_norm=(((10^breakpoint) - 10) / (ifelse(DBH_max>50,50,DBH_max) - 10))
      # breakpoint_norm=(((10^breakpoint) - 10) / (DBH_max - 10))
      
      prior_mean<-site_priors$serc_prior_mean
      prior_sd <- site_priors$serc_prior_sd
      if(length((Leaf_area_index%>%rename(siteID = site)%>%dplyr::filter(siteID == "SERC")%>%pull(Leaf_area_index) / 10))==0){
        next
      }else{
        stan_data<-list(
          N=bayesian_data%>%nrow(),
          x_min = 10,
          trunc_point = 10^breakpoint,
          trunc_upper = DBH_max,
          x = bayesian_data$dbh,
          LAI_norm = (Leaf_area_index%>%rename(siteID = site)%>%dplyr::filter(siteID == "SERC")%>%
                        pull(Leaf_area_index) / 10),
          breakpoint_norm=ifelse(breakpoint_norm >0, breakpoint_norm, 0),
          prior_mean = prior_mean,
          prior_sd=prior_sd
        )
        
        bayesian_model<-sampling(stan_alpha_laibreak_model, stan_data,
                                 chains = 4, iter = 9000, warmup=6000)
        
        results_serc<-rbind(results_serc,
                            summarise_draws(bayesian_model)%>%
                              filter(variable == "alpha")%>%
                              # mutate(alpha =  abs(mean)-1)%>%
                              mutate(site = site,breakpoint = 10^breakpoint,
                                     R2_kernel = performance::r2(lm(mean_log_density~log_x, bootstrap_kde_log))$R2))
        serc_bayesmodel[[site]]<-bayesian_model
      }
    }
  }
}

ggplot(results_serc%>%
         filter(R2_kernel >=0.8)%>%
         # filter(alpha >0 &alpha<5)%>%
         rename(Alpha_Mean=mean)%>%
         rename(plot_id = site)%>%
         # separate(plotID, into = c("siteID","number"), remove = FALSE)%>%
         # inner_join(remote_sensing_updated%>%
         #              group_by(siteID)%>%
         #              summarize(score=mean(score),
         #                        max_dbh = max(dbh)))%>%
         # filter(Site != "HEAL" & Site != "DEJU"&
         #          Site != "BONA" & Site != "SRER")%>%
         inner_join(sercfielddata%>%
                      dplyr::select(alpha, plot_id)%>%
                      # inner_join(tausandalphasresults%>%
                      #              rename(alpha=mean)%>%
                      select(mean, plot_id)%>%
                      rename(alpha=mean)), aes(alpha, Alpha_Mean,label = plot_id)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  # labs(title = "Recovering Ntot for NEON plots\n(Model Output)",
  #      x = "True Ntot", y = "Estimated Ntot") +
  geom_text() +# geom_errorbar(aes(ymin = Alpha_Lower_CI, ymax = Alpha_Upper_CI))+
  theme_minimal()+geom_smooth(method = "lm")#+geom_errorbar(aes(ymin=q5, ymax = q95))
# 
# 
# finalresults<-results_bootstrap%>%
#   # filter(R2 >=0.8)%>%
#   filter(site != "ABBY")%>%
#   filter(alpha >0 &alpha<5)%>%
#   rename(Alpha_Mean=mean)%>%
#   rename(siteID = site)%>%
#   # separate(plotID, into = c("siteID","number"), remove = FALSE)%>%
#   # inner_join(remote_sensing_updated%>%
#   #              group_by(siteID)%>%
#   #              summarize(score=mean(score),
#   #                        max_dbh = max(dbh)))%>%
#   # filter(Site != "HEAL" & Site != "DEJU"&
#   #          Site != "BONA" & Site != "SRER")%>%
#   inner_join(alpha_results%>%
#                dplyr::select(alpha, siteID))
# lm(data=finalresults,Alpha_Mean~alpha) %>%summary()


results_serc_ntot<-data.frame(
  variable = character(),
  site = character(),
  mean=numeric(),
  median=numeric(),
  alpha = numeric(),
  sd = numeric(),
  mad=numeric(),
  q5 = numeric(),
  q95 = numeric(),
  rhat = numeric(),
  rhat = numeric(),
  ess_bulk = numeric(),
  ess_tail = numeric(),
  breakpoint =numeric()
)
serc_bayesmodel_ntot<-list()
for(site in results_serc$site){
  site_data <- SERCfragmented%>%
    data.frame()%>%
    filter(plot_id == site)%>%
    filter(dbh<=50&dbh>=(results_serc%>%
                           rename(plot_id = site)%>%
                           filter(plot_id == site)%>%
                           pull(breakpoint)))
  # filter(dbh>=quantile(site_data$dbh, 0.75))
  
  
  binned_data<-forestscaling::logbin(site_data%>%
                                       # filter(dbh>=site_data$dbh%>%quantile(0.75) & dbh<=50)%>%
                                       pull(dbh), n =8)
  breakpoint<-results_serc%>%
    rename(plot_id = site)%>%
    filter(plot_id == site)%>%
    pull(breakpoint)
  # Add the avg_score column by calculating average scores for each bin
  binned_data <- binned_data %>%
    rowwise() %>%
    # mutate(avg_score = mean(
    #   site_data %>%
    #     filter(dbh >= bin_min & dbh < bin_max) %>%
    #     pull(score),
    #   na.rm = TRUE
    # ))
    # %>%
    ungroup()%>%
    tidyr::drop_na()
  
  
  # alpha_samples <- rnorm(n_alpha_samples, mean = taus%>%pull(mean), sd = taus%>%pull(sd))
  alpha_samples <- as.vector(rstan::extract(serc_bayesmodel[[site]], "alpha")$alpha)[1:1500]
  DBH_max=site_data%>%
    pull(dbh)%>%max(na.rm=TRUE)
  # Normalize Breakpoint (relative to x_min and x_max) and add 1
  
  breakpoint_norm=((breakpoint - 10) / (ifelse(DBH_max>50,50,DBH_max) - 10))
  # Example Stan data preparation
  stan_data <- list(
    K = nrow(binned_data),  # Number of bins
    bin_min = binned_data$bin_min,
    bin_max = binned_data$bin_max,
    N_obs = binned_data$bin_count,
    x_min = 10,
    LAI_norm = (4.7875000 / 10),
    breakpoint_norm=ifelse(breakpoint_norm >0, breakpoint_norm, 0),
    x_max=ifelse(DBH_max>50,50,DBH_max),
    # avg_scores=binned_data$avg_score,
    # Prior for alpha
    n_alpha_samples = length(alpha_samples),
    alpha_samples = alpha_samples,
    
    # Prior for Ntot
    N_tot_prior_mean =231,
    N_tot_prior_sd = 200)
  bayesianNtotmod<-sampling(stan_laibreakpoint_noscore_model,
                            stan_data,
                            chains =4,
                            warmup=2000,
                            iter=5000,cores = 4)
  print(paste(site))
  results_serc_ntot<-rbind(results_serc_ntot,
                           summarise_draws(bayesianNtotmod)%>%
                             filter(variable == "N_tot")%>%
                             # mutate(alpha =  abs(mean)-1)%>%
                             mutate(site = site))
  serc_bayesmodel_ntot[[site]]<-bayesianNtotmod
}

results_serc_ntot%>%
  inner_join(sercfielddata%>%
               select(plot_id, Ntotal)%>%
               rename(site=plot_id))

generate_abundance_data <- function(alpha_posterior, ntot_posterior, dbh_values, site_label, method_label, xmin = 10) {
  # Generate all combinations of alpha and ntot
  posterior_combinations <- expand.grid(
    alpha = alpha_posterior,
    ntot = ntot_posterior
  )
  
  # Repeat dbh values for each posterior combination
  dbh_rep <- rep(dbh_values, times = nrow(posterior_combinations))
  
  # Repeat alpha and ntot values for each dbh value
  alpha_rep <- rep(posterior_combinations$alpha, each = length(dbh_values))
  ntot_rep <- rep(posterior_combinations$ntot, each = length(dbh_values))
  
  # Calculate Pareto density for all combinations
  densPareto <- VGAM::dpareto(x = dbh_rep, shape = alpha_rep, scale = xmin)
  
  # Calculate abundance for all combinations
  abundance <- ntot_rep * densPareto
  
  # Create the final data frame
  abundance_data <- data.frame(
    size_class = dbh_rep,
    abundance = abundance,
    alpha = alpha_rep,
    ntot = ntot_rep,
    site = site_label,
    method = method_label
  )
  
  return(abundance_data)
}



for(site in results_serc_ntot$site){
  alpha_posterior <- as.vector(rstan::extract(serc_bayesmodel[[site]], "alpha")$alpha)
  ntot_posterior<-as.vector(rstan::extract(serc_bayesmodel_ntot[[site]], "N_tot")$N_tot)
  
  posterior_datasets1 <- lapply(1:12000, function(i) {
    # Generate abundance data for the indirect method (using posterior for alpha and ntot)
    abundance_data_indirect <- generate_abundance_data(
      alpha_posterior = alpha_posterior[i],
      ntot_posterior = ntot_posterior[i],
      dbh_values = dbh_values,
      site_label = site,
      method_label = "Remote"  # Label for indirect method
    )
    
    alpha_posterior_real<-as.vector(rstan::extract(sercrealalpha[[site]],"alpha")$alpha)
    ntot_posterior_real<-sercfieldfragmented%>%
      filter(plot_id == site)%>%
      filter(dbh>=10&dbh<=50)%>%
      nrow()
    
    # Generate abundance data for the direct method (using fixed posterior for ntot_real)
    abundance_data_direct <- generate_abundance_data(
      alpha_posterior = alpha_posterior_real[i],
      ntot_posterior = ntot_posterior_real,  # Fixed ntot for direct method
      dbh_values = dbh_values,
      site_label = site,
      method_label = "Field"  # Label for direct method
    )
    
    # Combine both indirect and direct method data for this posterior realization
    combined_data <- rbind(abundance_data_indirect, abundance_data_direct)
    combined_data$posterior_index <- site  # Add index to track posterior realization
    
    return(combined_data)
  })
  gc()
  saveRDS(posterior_datasets1, file=paste0("posterior_SERC//",site,".rds"))
  rm(posterior_datasets1)
  gc()
  print(paste(site))
}



final_datasets1<-fread("SERC_posterior_abundancecombos.csv")
final_datasets1<-final_datasets1%>%
  group_by(size_class, site, method)%>%
  mutate(log_abundance = log10(abundance),
         log_sizeclass = log10(size_class))%>%
  summarize(log_abundance_mean = mean(log_abundance),
            log_abundance_sd = sd(log_abundance))%>%
  # summarize(abundance_mean = mean(abundance),
  # abundance_sd=sd(abundance))%>%
  ungroup()



brmmodel<-brm(log_abundance_mean | mi(log_abundance_sd) ~ log_sizeclass + method + log_sizeclass:method,
              data= final_datasets1%>%
                mutate(log_sizeclass = log10(size_class))%>%
                filter(size_class %% 5 == 0), chains = 4,
              iter = 9000, warmup = 6000)
brmmodel
conditional_effects(brmmodel, "log_sizeclass:method")

ce <- conditional_effects(
  brmmodel,"log_sizeclass:method",method = "posterior_predict"
)

ce_data <- ce$log_sizeclass
ce_data$estimate_original <- 10^ce_data$estimate__
ce_data$lower_original <- 10^ce_data$lower__
ce_data$upper_original <- 10^ce_data$upper__

ce_data$sizeclass_original <- 10^ce_data$log_sizeclass

ggplot(ce_data, aes(x = sizeclass_original, y = estimate_original, fill = method)) +
  geom_line(aes(color = method), linewidth=1) +
  geom_ribbon(aes(ymin = lower_original, ymax = upper_original), alpha = 0.2) +
  labs(x = "DBH (cm)", y = "Trees in Plot (Abundance)") +
  scale_x_log10()+scale_y_log10()  # Optionally set a log10 scale for the x-axis


