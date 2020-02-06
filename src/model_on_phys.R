##
## effects SPECIES AND NITROGEN on physiology
##
###########################################################################################
###########################################################################################
# Load libraries
library(rstan) # Bayesian model compiler and sampler
options(mc.cores = parallel::detectCores()) # option to make Stan parallelize
library(bayesplot) # Plots of mcmc output


###########################################################################################


phys_data <- prep_trait_data()

compile_and_fit_phys_mixture_model <- 
  function(df, responsevar, iter = 100000){
    ## Bayesian mixture model for bimodal distributions
    
    #Stan model
    phys_model_mixture <- "
      data {
      int<lower=1> K; // number of peaks
      int<lower=1> N; // list of data
      real phys[N]; // observations
      real spp[N]; // covariates
      real Nit[N];
      }
      parameters {
      real intercept[K]; // locations of mixture components
      real<lower=0> sigma[K];  // scales of mixture components
      real B1; // regression coefficient
      real B2;
      real B3;
      }
      model {
      real ps[K]; // temp for log component densities
      intercept ~ normal(0, 30);
      B1 ~ normal(0,100);
      B2 ~ normal(0,100);
      B3 ~ normal(0,100);
      for (n in 1:N) {
      for (k in 1:K) {
      ps[k] = normal_lpdf(phys[n] | intercept[k] + spp[n]*B1 + Nit[n]*B2 + spp[n]*Nit[n]*B3 , sigma[k]);
      }
      target += log_sum_exp(ps);
      }
      }
      generated quantities {
      real draws1[N];
      vector[1] BG_val;
      vector[1] SC_val;
      for (n in 1:N) {
      for (k in 1:K) {
      draws1[n] = normal_rng( intercept[k] + spp[n]*B1 + Nit[n]*B2 + spp[n]*Nit[n]*B3, sigma[k]);
      }
      }
      for (k in 1:K) {
      BG_val[1] = intercept[k];
      SC_val[1] = intercept[k] + B1;
      }
      }
    "
    
    # Compile model
    comp.mixture <- stan_model(model_code = phys_model_mixture, 
                               model_name = 'phys.model.mixture')
    
    # Select data
    modeldat <- list(
      'N' = nrow(df),
      'K' = 2, # Number of peaks in multimodal distribution
      'phys' = responsevar,
      'spp' = as.numeric(df$spp) - 1,  # Want zeros and ones for species
      'Nit' = df$nitrogen
    )
    # Perform MCMC sampling
    fit <-
      sampling(
        comp.mixture,
        data = modeldat,
        iter = iter,
        warmup = iter / 2,
        thin = 1,
        chains = 2
      )
    return(fit)
}


compile_and_fit_phys_normal_model <-
  function(df, responsevar, iter = 100000){
    ## Bayesian model for normal distributions
    
    # Stan model
    phys_model_normal <- "
      data {
      int<lower=0> N; //list of data
      vector[N] spp;
      vector[N] Nit;
      vector[N] phys;
      }
      parameters {
      real B0;
      real B1;
      real B2;
      real B3;
      real<lower=0> sigma;
      }
      transformed parameters{
      vector[N] mu;
      for (n in 1:N)
      mu[n] = B0 + spp[n]*B1 + Nit[n]*B2 + spp[n]*Nit[n]*B3;
      }
      model {
      sigma ~ cauchy(0, 10);
      B0 ~ normal(0,10000);
      B1 ~ normal(0,10000);
      B2 ~ normal(0,10000);
      B3 ~ normal(0,10000);
      phys ~ normal(mu,sigma);
      }
      generated quantities{
      vector[N] draws1;
      vector[1] BG_val;
      vector[1] SC_val;
      for(n in 1:N)
      draws1[n] = normal_rng(mu[n],sigma); //posterior draws
      BG_val[1] = B0;
      SC_val[1] = B0 + B1;
      }
    "
    
    # Compile model
    comp.normal <- stan_model(model_code = phys_model_normal, 
                               model_name = 'phys.model.normal')
    
    # Select data
    modeldat <- list(
      'N' = nrow(df),
      'phys' = responsevar,
      'spp' = as.numeric(df$spp) - 1,  # Want zeros and ones for species
      'Nit' = df$nitrogen
    )
    # Perform MCMC sampling
    fit <-
      sampling(
        comp.normal,
        data = modeldat,
        iter = iter,
        warmup = iter / 2,
        thin = 1,
        chains = 2
      )
    return(fit)
}


compile_and_fit_phys_gamma_model <- 
  function(df, responsevar, iter = 100000){
    ## Bayesian model for gamma distributions
    
    # Stan model
    phys_model_gamma <- "
      data {
      int<lower=0> N; //list of data
      vector[N] spp;
      vector[N] Nit;
      vector[N] phys;
      }
      parameters {
      real B0;
      real B1;
      real B2;
      real B3;
      real<lower=0> phi;
      }
      transformed parameters{
      vector[N] mu;
      vector[N] alpha;
      vector[N] beta;
      for (n in 1:N)
      mu[n] = exp( B0 + spp[n]*B1 + Nit[n]*B2 + spp[n]*Nit[n]*B3 );
      alpha = mu .* mu / phi;
      beta = mu / phi;
      }
      model {
      phi ~ cauchy(0, 10);
      B0 ~ normal(0,1000);
      B1 ~ normal(0,1000);
      B2 ~ normal(0,1000);
      B3 ~ normal(0,1000);
      phys ~ gamma(alpha,beta);
      }
      generated quantities{
      vector[N] draws1;
      vector[1] BG_val;
      vector[1] SC_val;
      for(n in 1:N)
      draws1[n] = gamma_rng(alpha[n],beta[n]); //posterior draws
      BG_val[1] = exp( B0 );
      SC_val[1] = exp (B0 + B1 );
      }
    "
    
    # Compile model
    comp.gamma <- stan_model(model_code = phys_model_gamma, 
                              model_name = 'phys.model.gamma')
    
    # Select data
    modeldat <- list(
      'N' = nrow(df),
      'phys' = responsevar,
      'spp' = as.numeric(df$spp) - 1,  # Want zeros and ones for species
      'Nit' = df$nitrogen
    )
    # Perform MCMC sampling
    fit <-
      sampling(
        comp.gamma,
        data = modeldat,
        iter = iter,
        warmup = iter / 2,
        thin = 1,
        chains = 2
      )
    return(fit)
     
}


write_model_statistics_mixture <- function(fit, name) {
  ## gather 95% CIs
  print("Calculating model statistics..")
  model_stats <-
    cbind(
      summary(fit)$summary[, 1, drop = F],
      summary(fit)$summary[, 4, drop = F],
      summary(fit)$summary[, 8, drop = F],
      summary(fit)$summary[, 10, drop = F]
    )
  write.csv(model_stats,
            file = paste("output/stanmodel_results_", name, "_stats.csv", sep = ""))
  
  ## Pr values
  print("Calculating and writing Pr values..")
  # Retain only stats for parameters
  Pr_vals <-
    rbind(model_stats["intercept[1]", ],
          model_stats["intercept[2]", ],
          model_stats["B1", ],
          model_stats["B2", ],
          model_stats["B3", ])
  Pr_vals <- Pr_vals[,-4] # Drop rhat for calculating Pr
  Pr_calc <- as.data.frame(NULL)
  # Iterate through parameters
  for (r in 1:nrow(Pr_vals)) {
    max_ <- max(Pr_vals[r,])
    min_ <- min(Pr_vals[r,])
    # Calculate the proportion of overlap on zero
    if (max_ < 0 & min_ < 0) {
      pr. <- 1
    } else if (max_ > 0 & min_ > 0) {
      pr. <- 1
    } else if (abs(max_) > abs(min_)) {
      pr. <- abs(max_) / (abs(min_) + abs(max_))
    } else if (abs(max_) < abs(min_)) {
      pr. <- abs(min_) / (abs(min_) + abs(max_))
    } else {
      pr. <- "NA"
    }
    Pr_calc[1, r] <- pr.
  }
  c("intercept[1]", "intercept[2]", "B1", "B2", "B3")
  write.csv(Pr_calc,
            file = paste("output/stanmodel_results_", name, "_pr.csv", sep = ""))
}


model_phys <- function() {
  # Fit models on three principal components:
  # Load data
  dat <- prep_trait_data()
  
  # Photosynthetic rate - data are normal/mixture
  fit <-
    compile_and_fit_phys_mixture_model(dat, responsevar = dat$photo)
  plot_main_effects(fit, name = "photo")
  write_model_statistics_mixture(fit, name = "photo")
  plot_posterior_checks(fit, responsevar = dat$photo, name = "photo")
  
  # Stomatal conductance - data are gamma
  fit <-
    compile_and_fit_phys_gamma_model(dat, responsevar = dat$cond + 1)
  plot_main_effects(fit, name = "cond")
  write_model_statistics(fit, name = "cond")
  plot_posterior_checks(fit, responsevar = dat$photo + 1, name = "cond")
  
  # Intercellular CO2 - data are normal/mixture
  fit <-
    compile_and_fit_phys_mixture_model(dat, responsevar = dat$Ci)
  plot_main_effects(fit, name = "Ci")
  write_model_statistics_mixture(fit, name = "Ci")
  plot_posterior_checks(fit, responsevar = dat$Ci, name = "Ci")
}





# 
# 
# 
# 
# ## Data are normal/mixture
# phys.model.mixture(
#   plot_title = "Photosynthetic rate",
#   responsevar = phys.data$photo,
#   outfile1 = "Stanmodel_phys_results_PHOTO_stats.csv",
#   outfile2 = "Stanmodel_phys_results_PHOTO_Pr.csv",
#   normalityplot1 = "Stanmodel_phys_results_PHOTO_posterior_check1.pdf",
#   normalityplot2 = "Stanmodel_phys_results_PHOTO_posterior_check2.pdf"
# )
# ## Data are gamma
# phys.model.gamma(
#   plot_title = "Stomatal conductance",
#   responsevar = phys.data$cond + 1,
#   ## Transformation of conductance necessary, won't converge otherwise
#   outfile1 = "Stanmodel_phys_results_COND_stats.csv",
#   outfile2 = "Stanmodel_phys_results_COND_Pr.csv",
#   normalityplot1 = "Stanmodel_phys_results_COND_posterior_check1.pdf",
#   normalityplot2 = "Stanmodel_phys_results_COND_posterior_check2.pdf"
# )
# ## Data are normal/mixture
# phys.model.mixture(
#   plot_title = "Intercellular CO2",
#   responsevar = phys.data$Ci,
#   outfile1 = "Stanmodel_phys_results_Ci_stats.csv",
#   outfile2 = "Stanmodel_phys_results_Ci_Pr.csv",
#   normalityplot1 = "Stanmodel_phys_results_Ci_posterior_check1.pdf",
#   normalityplot2 = "Stanmodel_phys_results_Ci_posterior_check2.pdf"
# )
# ## Data are gamma
# phys.model.gamma(
#   plot_title = "Transpiration rate",
#   responsevar = phys.data$Trmmol,
#   outfile1 = "Stanmodel_phys_results_Trmmol_stats.csv",
#   outfile2 = "Stanmodel_phys_results_Trmmol_Pr.csv",
#   normalityplot1 = "Stanmodel_phys_results_Trmmol_posterior_check1.pdf",
#   normalityplot2 = "Stanmodel_phys_results_Trmmol_posterior_check2.pdf"
# )
# # Data are normal
# phys.model.normal(
#   plot_title = "Vapor pressure deficit",
#   responsevar = phys.data$VpdL,
#   outfile1 = "Stanmodel_phys_results_VpdL_stats.csv",
#   outfile2 = "Stanmodel_phys_results_VpdL_Pr.csv",
#   normalityplot1 = "Stanmodel_phys_results_VpdL_posterior_check1.pdf",
#   normalityplot2 = "Stanmodel_phys_results_VpdL_posterior_check2.pdf"
# )
# ## Data are gamma
# phys.model.gamma(
#   plot_title = "Water use efficiency",
#   responsevar = phys.data$iWUE,
#   outfile1 = "Stanmodel_phys_results_iWUE_stats.csv",
#   outfile2 = "Stanmodel_phys_results_iWUE_Pr.csv",
#   normalityplot1 = "Stanmodel_phys_results_iWUE_posterior_check1.pdf",
#   normalityplot2 = "Stanmodel_phys_results_iWUE_posterior_check2.pdf"
# )