
stan_alpha_laibreak<-"data {
  int<lower=0> N;                  // Number of observations
  real<lower=0> x_min;             // Minimum threshold for the distribution (e.g., 10)
  real<lower=0> trunc_point;       // Observed data starts at this truncation threshold (e.g., 20)
  real<lower=trunc_point> trunc_upper;  // Observed data is capped at this upper threshold (e.g., 50)
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

sitevector<-c("ABBY", "BART", "BLAN", "BONA", "CLBJ", "DELA", "GUAN", "HARV", "JERC", 
              "KONZ", "LENO", "MLBS", "NIWO", "OSBS", "RMNP", "SCBI", "SERC", "SJER", 
              "SOAP", "TALL", "TEAK", "WREF", "YELL")

results_bootstrap<-data.frame(
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
  R2_kernel= numeric()
)
results_model<-list()
for(site in sitevector){
  site_data <- subset(remote_sensing_updated_sitefix%>%
                        data.frame()%>%
                        select(-geometry)%>%
                        filter(dbh<=50),
                      siteID == site)

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
      filter(any(abs(x - observed_dbhs) <= 0.5)) %>%
      ungroup()
    # Assume original_data and observed_dbhs are defined
    n_bootstrap <- 1000  # Number of bootstrap samples
    xmin <- min(observed_dbhs)
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

    breakpoint = data%>%
      filter(log_x<=quantile(log10(original_data), 0.75))%>%
      # filter(mean_log_density  == max(mean_log_density ))%>%
      filter(log_x==max(log_x))%>%
      pull(log_x)
    
    segmented_model<-lm(mean_log_density~log_x,bootstrap_kde_log%>%
                          filter(log_x>=breakpoint))%>%selgmented()
    
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
      site_priors <- neonregressionoutput %>%
        filter(siteID==site)
      DBH_max=bayesian_data$dbh%>%max()
      breakpoint_norm=(((10^breakpoint) - 10) / (ifelse(DBH_max>50,50,DBH_max) - 10))
      prior_mean<-site_priors$prior_mean
      prior_sd <- site_priors$prior_sd
      if(length((Leaf_area_index%>%rename(siteID = site)%>%dplyr::filter(siteID == site)%>%pull(Leaf_area_index) / 10))==0){
        next
      }else{
        stan_data<-list(
          N=bayesian_data%>%nrow(),
          x_min = 10,
          trunc_point = 10^breakpoint,
          trunc_upper = DBH_max,
          x = bayesian_data$dbh,
          LAI_norm = (Leaf_area_index%>%rename(siteID = site)%>%dplyr::filter(siteID == site)%>%
                        pull(Leaf_area_index) / 10),
          breakpoint_norm=ifelse(breakpoint_norm >0, breakpoint_norm, 0),
          prior_mean = prior_mean,
          prior_sd=prior_sd
        )
        
        bayesian_model<-sampling(stan_alpha_laibreak_model, stan_data,
                                 chains = 4, iter = 9000, warmup=6000)
    
        results_bootstrap<-rbind(results_bootstrap,
                                 summarise_draws(bayesian_model)%>%
                                   filter(variable == "alpha")%>%
                                   # mutate(alpha =  abs(mean)-1)%>%
                                   mutate(site = site,
                                          R2_kernel = performance::r2(lm(mean_log_density~log_x, bootstrap_kde_log))$R2))
        print(paste(site))
        
        results_model[[site]]<-bayesian_model
      }
    }
  }
}

alphaplot<-ggplot(results_bootstrap%>%
                    filter(R2_kernel >=0.8)%>%
                    rename(Alpha_Mean=mean)%>%
                    rename(siteID = site)%>%
                    inner_join(alpha_results%>%
                                 dplyr::select(alpha, siteID)%>%
                                 select(alpha, siteID))%>%
                    inner_join(site_eval%>%
                                 mutate(slope_diff = abs(slope - 1))), aes(alpha, Alpha_Mean,color=slope_diff,
                                                                           label = siteID)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(x = expression(paste("True ", alpha)), y = expression(paste("Remote Sensing Estimated ", alpha))) +
  geom_text(nudge_y =0.1) +
  theme_minimal()+geom_smooth(method = "lm")
alphaplot

ggsave("alphaplot.png",device=cairo_pdf())

finalresults<-results_bootstrap%>%
  filter(R2_kernel >=0.8)%>%
  filter(alpha >0 &alpha<5)%>%
  rename(Alpha_Mean=mean)%>%
  rename(siteID = site)%>%
  inner_join(alpha_results%>%
               dplyr::select(alpha, siteID))
lm(data=finalresults,Alpha_Mean~alpha) %>%summary()