library(data.table)
library(sf)
library(dplyr)
HarvFieldData<-fread("R:/Adam Eichenwald/HarvFieldData.csv")%>%
  data.frame()
HarvFieldData<-st_as_sf(HarvFieldData%>%
                          separate(geometry,sep = "\\|", 
                                   into =c("x","y")),
                        coords = c("x","y"),crs=32618)

write_sf(HarvFieldData, "HarvFieldData.shp")

(head(HarvFieldData))

bencrownfragmented<-read_sf("R:/Adam Eichenwald/HARVplotsGEOJan25.shp")%>%
  st_set_crs(32618)
# write_sf(bencrown, "R:/Adam Eichenwald/HarvForestGeo.shp")


# Define bounding box
g_bbox <- st_bbox(c(xmin = 731592, 
                    ymin = 4713222, 
                    xmax = 732295, 
                    ymax = 4713725), crs = 32618)

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

# Spatial join between HarvFieldData and plots.
HarvFieldfragmented <- st_join(HarvFieldData, plots)
mapview::mapview(plots, legend = FALSE)

# bencrownfragmented <- st_join(bencrown, plots, largest = TRUE)
# bencrownfragmented<-bencrown
# plot_remove<-st_intersects(plots%>%st_transform("WGS84"),swamp, sparse = FALSE)%>%data.frame()%>%cbind(plots)%>%filter(.==FALSE)


Harv2019<-terra::rast("C:\\Users\\adam.local\\Downloads\\NEON_D01_HARV_DP3_731000_4713000_CHM.tif")
Harv2019_1<-terra::rast("C:\\Users\\adam.local\\Downloads\\harvfield\\NEON_D01_HARV_DP3_732000_4713000_CHM.tif")
mosaicHarv<-mosaic(Harv2019, Harv2019_1)

print(mosaicHarv)


e<-ext(read_sf("C:\\Users\\adam.local\\Documents\\HarvForestGeo.shp")%>%
         st_set_crs(32618))
mosaic2<-terra::crop(mosaicHarv,e)
mosaic2<-terra::rast("Harvheight.tif")
mapview::mapview(mosaic2)


HarvField<-read_sf("HARV_FieldData_ForestGEO_plotcorrected.shp")
HarvField<-HarvField%>%
  drop_na(plot_id)%>%
  filter(dbh>=10&dbh<=50)

harvdata<-list()
harvrealalpha<-list()
for (i in unique(HarvField$plot_id)){
  test<-HarvField%>%
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
  harvdata[[length(harvdata)+1]]<-summarize_draws(fit1)%>%
    mutate(plot_id = i,
           x_min = x_min,
           Ntotal = N)
  harvrealalpha[[i]]<-fit1
}


harvfielddata<-harvdata%>%
  rbindlist()%>%
  filter(variable =="alpha")%>%
  mutate(label = "Field")

library(tidyverse)
library(itcSegment)
bencrownfragmented$perimeter<-as.numeric(st_perimeter(bencrownfragmented))
bencrownfragmented$area<-as.numeric(st_area(bencrownfragmented))
bencrownfragmented<-bencrownfragmented%>%
  dplyr::mutate(Diameter = 0.5*(sqrt(perimeter^2-(8*area)))) %>%
  # mutate(Diameter1 = sqrt(crown_area))
  # Calculate dbh using formula from Jucker et al. 2017
  mutate(dbh = dbh(H=height, CA = Diameter, biome = 12))

bencrownfragmented<-bencrownfragmented%>%
  filter(plot_id != "plot_1"&plot_id != "plot_2"&
           plot_id!= "plot_10"& plot_id != "plot_9")







predicted<-list()
newdata1<-list()
bencrownfragmented<-bencrownfragmented%>%
  st_transform("WGS84")
plotswgs<-st_centroid(plots%>%st_transform("WGS84"))%>%st_coordinates()%>%
  cbind(plots)
for (i in na.omit(unique(bencrownfragmented$plot_id))){
  bencrownfragmentedsubset<-subset(bencrownfragmented, bencrownfragmented$plot_id == i)
  # SERCfragmentedsubset <- SERCfragmented
  newdata <- expand_grid(Max_Height_meters = (bencrownfragmented%>%
                                                data.frame()%>%
                                                filter(plot_id == i)%>%
                                                summarize(maxht = max(height)))$maxht* 3.281,
                         ELEV = 1141.73)%>%
    
    mutate(
      Lat= (plotswgs%>%
              filter(plot_id == i))$Y,
      Lon = (plotswgs%>%
               filter(plot_id == i))$X)
  
  newdata1[[length(newdata1)+1]]<-newdata%>%
    mutate(
      plot_id = i)
}

newdata1

library(ranger)

harvpredicse<-predict(rangermodel, newdata1%>%
                        rbindlist(),
                      # rename(Max_Height=Max_Height_meters,
                      #        LON=Lon,LAT=Lat,elevation = ELEV)
                      predict.all=TRUE)
# Calculate mean and standard deviation for each data point (each row)
harv_prior_mean <- apply(harvpredicse$predictions, 1, mean) # Mean for each row
harv_prior_sd <- apply(harvpredicse$predictions, 1, sd)      # Standard deviation for each row

harvregressionoutput<-(newdata1%>%
                         rbindlist())%>%
  cbind(harv_prior_mean)%>%
  cbind(harv_prior_sd)

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
# sitevector<-results_bootstrap$site
results_harv<-data.frame(
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
harv_bayesmodel<-list()
library(splus2R)
library(segmented)
for(site in unique(bencrownfragmented$plot_id)){#unique(remote_sensing_updated$siteID)[-6]) {
  site_data <- subset(bencrownfragmented%>%
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
      site_priors <- harvregressionoutput %>%
        # mutate(slope_prior_mean = -(prior_mean + 1))%>%
        filter(plot_id==site)
      DBH_max=bayesian_data$dbh%>%max()
      breakpoint_norm=(((10^breakpoint) - 10) / (ifelse(DBH_max>50,50,DBH_max) - 10))
      # breakpoint_norm=(((10^breakpoint) - 10) / (DBH_max - 10))
      
      prior_mean<-site_priors$harv_prior_mean
      prior_sd <- site_priors$harv_prior_sd
      if(length((Leaf_area_index%>%rename(siteID = site)%>%dplyr::filter(siteID == "SERC")%>%pull(Leaf_area_index) / 10))==0){
        next
      }else{
        stan_data<-list(
          N=bayesian_data%>%nrow(),
          x_min = 10,
          trunc_point = 10^breakpoint,
          trunc_upper = DBH_max,
          x = bayesian_data$dbh,
          LAI_norm = (Leaf_area_index%>%rename(siteID = site)%>%dplyr::filter(siteID == "HARV")%>%
                        pull(Leaf_area_index) / 10),
          breakpoint_norm=ifelse(breakpoint_norm >0, breakpoint_norm, 0),
          prior_mean = prior_mean,
          prior_sd=prior_sd
        )
        
        bayesian_model<-sampling(stan_alpha_laibreak_model, stan_data,
                                 chains = 4, iter = 9000, warmup=6000)
        
        results_harv<-rbind(results_harv,
                            summarise_draws(bayesian_model)%>%
                              filter(variable == "alpha")%>%
                              # mutate(alpha =  abs(mean)-1)%>%
                              mutate(site = site,breakpoint = 10^breakpoint,
                                     R2_kernel = performance::r2(lm(mean_log_density~log_x, bootstrap_kde_log))$R2))
        harv_bayesmodel[[site]]<-bayesian_model
      }
    }
  }
}





results_harv_ntot<-data.frame(
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
harv_bayesmodel_ntot<-list()
for(site in (HarvField%>%
             filter(plot_id != "plot_1"&plot_id != "plot_2"&
                    plot_id!= "plot_10"& plot_id != "plot_9"))$plot_id%>%unique()){
  site_data <- bencrownfragmented%>%
    data.frame()%>%
    filter(plot_id == site)%>%
    filter(dbh<=50&dbh>=(results_harv%>%
                           rename(plot_id = site)%>%
                           filter(plot_id == site)%>%
                           pull(breakpoint)))
  # filter(dbh>=quantile(site_data$dbh, 0.75))
  
  
  binned_data<-forestscaling::logbin(site_data%>%
                                       # filter(dbh>=site_data$dbh%>%quantile(0.75) & dbh<=50)%>%
                                       pull(dbh), n =8)
  breakpoint<-results_harv%>%
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
  alpha_samples <- as.vector(rstan::extract(harv_bayesmodel[[site]], "alpha")$alpha)[1:1500]
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
    LAI_norm = (5.4269231 / 10),
    breakpoint_norm=ifelse(breakpoint_norm >0, breakpoint_norm, 0),
    x_max=ifelse(DBH_max>50,50,DBH_max),
    # avg_scores=binned_data$avg_score,
    # Prior for alpha
    n_alpha_samples = length(alpha_samples),
    alpha_samples = alpha_samples,
    
    # Prior for Ntot
    N_tot_prior_mean =630,
    N_tot_prior_sd = 200)
  bayesianNtotmod<-sampling(stan_laibreakpoint_noscore_model,
                            stan_data,
                            chains =4,
                            warmup=2000,
                            iter=5000,cores = 4)
  print(paste(site))
  results_harv_ntot<-rbind(results_harv_ntot,
                           summarise_draws(bayesianNtotmod)%>%
                             filter(variable == "N_tot")%>%
                             # mutate(alpha =  abs(mean)-1)%>%
                             mutate(site = site))
  harv_bayesmodel_ntot[[site]]<-bayesianNtotmod
}

results_harv_ntot%>%
  inner_join(harvfielddata%>%
               select(plot_id, Ntotal)%>%
               rename(site=plot_id))
dbh_values <- seq(10, 50, length.out=100)  # Example DBH size classes

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



for(site in (HarvField%>%
             filter(plot_id != "plot_1"&plot_id != "plot_2"&
                    plot_id!= "plot_10"& plot_id != "plot_9"))$plot_id%>%unique()){
  alpha_posterior <- as.vector(rstan::extract(harv_bayesmodel[[site]], "alpha")$alpha)
  ntot_posterior<-as.vector(rstan::extract(harv_bayesmodel_ntot[[site]], "N_tot")$N_tot)
  
  posterior_datasets1 <- lapply(1:12000, function(i) {
    # Generate abundance data for the indirect method (using posterior for alpha and ntot)
    abundance_data_indirect <- generate_abundance_data(
      alpha_posterior = alpha_posterior[i],
      ntot_posterior = ntot_posterior[i],
      dbh_values = dbh_values,
      site_label = site,
      method_label = "Remote"  # Label for indirect method
    )
    
    alpha_posterior_real<-as.vector(rstan::extract(harvrealalpha[[site]],"alpha")$alpha)
    ntot_posterior_real<-HarvField%>%
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
  saveRDS(posterior_datasets1, file=paste0("posterior_HARV_allometryfix//",site,".rds"))
  rm(posterior_datasets1)
  gc()
  print(paste(site))
}






final_datasets1<-fread("HARV_posterior_abundancecombos_allometryfix.csv")
final_datasets1<-final_datasets1%>%
  group_by(size_class, site, method)%>%
  mutate(log_abundance = log10(abundance),
         log_sizeclass = log10(size_class))%>%
  # summarize(log_abundance_mean = mean(log_abundance),
  #           log_abundance_sd = sd(log_abundance))%>%
  summarize(abundance_mean = mean(abundance),
            abundance_sd=sd(abundance))%>%
  ungroup()



brmmodel<-brm(log_abundance_mean | mi(log_abundance_sd) ~ log_sizeclass + method + log_sizeclass:method + 
                (1|site),
              data= final_datasets1%>%
                mutate(log_sizeclass = log10(size_class)), chains = 4, cores = 4,
              iter = 9000, warmup = 6000)

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
  labs(x = "DBH (cm)", y = "Trees in Plot (Abundance)") #+
# scale_x_log10()+scale_y_log10()  # Optionally set a log10 scale for the x-axis

library(bayestestR)

sexit(brmmodel)

test<-read_sf("HarvFieldData_heightinclude.shp")

ggplot(test%>%
         data.frame()%>%
         select(RASTERVALU, dbh)%>%
         rename(height=RASTERVALU)%>%
         mutate(label = "Field")%>%
         rbind(bencrownfragmented%>%
                 data.frame()%>%
                 select(height,dbh)%>%
                 mutate(label = "Remote")), aes(height, dbh))+geom_point()+facet_wrap("label")+
  scale_y_log10()



library(dplyr)
library(rstan)

# Prepare the dataset
final_datasets1 <- final_datasets1 %>%
  mutate(
    log_sizeclass = log10(size_class),
    site_id = as.integer(factor(site)), # Convert site factor to integer IDs
    method_id = as.integer(factor(method)) # Convert method factor to integer IDs
  )

# Response variable
Y <- final_datasets1$abundance_mean

# Measurement error (NA-safe handling)
noise <- ifelse(is.na(final_datasets1$abundance_sd), 0, final_datasets1$abundance_sd)

# Non-missing values for measurement error
Jme <- which(!is.na(final_datasets1$abundance_sd))
Nme <- length(Jme)

# Fixed effect predictors (including interaction)
X <- model.matrix(~  size_class*method, data = final_datasets1)

# Grouping structure for site random effect
J_1 <- final_datasets1$site
N_1 <- length(unique(J_1))

# Random intercept predictor (usually just a column of ones)
Z_1_1 <- rep(1, nrow(final_datasets1))

# Define the stan data list
stan_data <- list(
  N = nrow(final_datasets1),
  Y = Y,
  noise = noise,
  Nme = Nme,
  Jme = Jme,
  K = ncol(X),
  X = X,
  Kc = ncol(X) - 1,  # Number of predictors after centering
  N_1 = N_1,
  M_1 = 1,
  J_1 = J_1,
  Z_1_1 = Z_1_1,
  prior_only = 0,
  x_min = 10
)

stan_pareto_data_test<-
  "data {
  int<lower=1> N;             // Number of observations
  vector[N] dbh;              // Predictor (dbh), must be ≥ 10
  real<lower=0> abundance[N];  // Response variable (abundance)
  int<lower=1, upper=2> method[N];  // Method as a categorical index (1 or 2)
}

parameters {
  real intercept[2];  // Separate intercepts for method 1 and 2
  real beta_dbh[2];   // Separate slopes for dbh for method 1 and 2
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] mu;
  
  for (i in 1:N) {
    mu[i] = exp(intercept[method[i]] + beta_dbh[method[i]] * log(dbh[i]));
  }
}

model {
  // Priors
  intercept ~ normal(0, 5);
  beta_dbh ~ normal(0, 5);
  sigma ~ normal(0, 2);

  // Likelihood
  abundance ~ lognormal(log(mu), sigma);  
}

"

stan_pareto_data_test_compile<-stan_model(model_code = stan_pareto_data_test)
# Assuming your data is ready:
final_datasets2<-final_datasets1%>%
  # filter(method=="Remote")%>%
  mutate(dbh=size_class)#%>%
# mutate(method = ifelse(method == "Field",1,2))

stan_data <- list(
  N = nrow(final_datasets2),
  abundance = final_datasets2$abundance_mean,
  dbh = final_datasets2$size_class,
  method = final_datasets2$method
)

# Fit the model
fit <- sampling(
  stan_pareto_data_test_compile,  # path to the Stan code
  data = stan_data,
  chains = 4,
  cores = 4,
  iter = 5000,
  warmup = 3000
)

fit_power <- brms::brm(
  formula = brms::bf(abundance_mean ~ c * dbh^(-d),  
                     c + d ~ 1, 
                     nl=TRUE),
  data    = final_datasets2,  # Replace with your actual dataset
  # prior   = prior(lognormal(0,0.5), nlpar = "c", lb = 0) + 
  #   prior(lognormal(0,0.5), nlpar = "d", lb = 0),
  family  = gaussian(),  # Use lognormal if abundance is skewed, otherwise use gaussian()
  # control = list(adapt_delta = 0.99),
  chains = 3, cores = 3
)
fit_power <- brms::brm(
  formula = brms::bf(abundance_mean ~ c * dbh^(-d),  
                     c ~ method,  # Method affects the scaling factor
                     d ~ method,  # Method affects the power exponent
                     nl=TRUE),
  data    = final_datasets2,
  prior   = prior(lognormal(0,0.5), nlpar = "c", lb = 0) + 
    prior(lognormal(0,0.5), nlpar = "d", lb = 0),
  family  = lognormal(),
  chains  = 4, cores = 4, warmup=6000, iter = 9000
)
full_values <-distinct(final_datasets1, size_class)%>%pull(size_class)
downsampled_values <- full_values[seq(1, length(full_values), length.out = 15)]
final_datasets1<-fread("HARV_posterior_abundancecombos_allometryfix.csv")

# Create a list of datasets with random samples for each group
result_list <- purrr::map(1:20, ~ {  # Change 5 to the number of datasets you want
  final_datasets1 %>%
    distinct(alpha, method, site) %>%
    group_by(method, site) %>%
    mutate(random_value = sample(alpha, size = 1)) %>%
    filter(alpha == random_value) %>%
    select(-random_value)%>%
    inner_join(final_datasets1)%>%
    filter(size_class %in% downsampled_values)# Remove the random_value column
})

# View the resulting list
print(result_list)


fit_power_multiple <- brms::brm_multiple(
  formula = brms::bf(abundance ~ c * size_class^(-d),  
                     c ~ method,  # Method affects the scaling factor
                     d ~ method,  # Method affects the power exponent
                     nl=TRUE),
  data    = result_list,
  prior   = prior(lognormal(0,0.5), nlpar = "c", lb = 0) + 
    prior(lognormal(0,0.5), nlpar = "d", lb = 0),
  family  = lognormal(),
  chains  = 4,warmup=1000, iter = 4000
)

library(future)
plan(multisession, workers = 11)
# Get the time before running the function
start_time <- Sys.time()
fit_power_multiple_update<-update(fit_power_multiple, newdata=result_list[1:11])

# Get the time before running the function
end_time <- Sys.time()
plan(sequential)


# Calculate the elapsed time
elapsed_time <- end_time - start_time

# Print the result
print(elapsed_time)

# View the summary
summary(fit_power)

load("C:/Users/adam.local/Downloads/HARV_brmmultiple_3000.RData")

ceffplot<-plot(conditional_effects(fit_power_multiple_update, "size_class:method"), plot = FALSE)[[1]]


ceffplot<-plot(conditional_effects(fit_power_multiple_update, "size_class:method", ndraws=800),
               plot = FALSE)[[1]]
ceffplot+xlab("DBH (cm)")+ylab("Tree Abundance")+theme_minimal()+
  scale_x_log10()+scale_y_log10()+labs(fill = "Measurement Method",
                                       color = "Measurement Method")+
  labs(title = "Smithsonian Environmental Research Center")
