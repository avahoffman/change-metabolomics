###########################################################################################
##
## effects SPECIES AND NITROGEN on physiology
## Keep in mind, recommended lawn application is about 5g/m2 of Nitrogen!
##
###########################################################################################
###########################################################################################
# Load libraries
library(rstan) # Bayesian model compiler and sampler
options(mc.cores = parallel::detectCores()) # option to make Stan parallelize
rstan_options(auto_write = TRUE)
library(bayesplot) # Plots of mcmc output


###########################################################################################


compile_and_fit_phys_mixture_model <-
  function(df, responsevar, iter = 100000) {
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
      'K' = 2,
      # Number of peaks in multimodal distribution
      'phys' = responsevar,
      'spp' = as.numeric(df$spp) - 1,
      # Want zeros and ones for species
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
  function(df, responsevar, iter = 100000) {
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
      'spp' = as.numeric(df$spp) - 1,
      # Want zeros and ones for species
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
  function(df, responsevar, iter = 100000) {
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
      'spp' = as.numeric(df$spp) - 1,
      # Want zeros and ones for species
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


model_phys <- function() {
  # Fit models on three principal components:
  # Load data
  dat <- prep_trait_data()
  
  # Photosynthetic rate - data are normal/mixture
  fit <-
    compile_and_fit_phys_mixture_model(dat,
                                       responsevar = dat$photo,
                                       iter = iter)
  plot_main_effects(fit, name = "photo")
  write_model_statistics(fit, name = "photo", mixture = T)
  plot_posterior_checks(fit,
                        responsevar = dat$photo,
                        name = "photo")
  
  # Stomatal conductance - data are gamma
  fit <-
    compile_and_fit_phys_gamma_model(dat,
                                     responsevar = dat$cond + 1,
                                     iter = iter)
  plot_main_effects(fit, name = "cond")
  write_model_statistics(fit, name = "cond")
  plot_posterior_checks(fit,
                        responsevar = dat$cond + 1,
                        name = "cond")
  
  # Intercellular CO2 - data are normal/mixture
  fit <-
    compile_and_fit_phys_mixture_model(dat,
                                       responsevar = dat$Ci,
                                       iter = iter)
  plot_main_effects(fit, name = "Ci")
  write_model_statistics(fit, name = "Ci", mixture = T)
  plot_posterior_checks(fit,
                        responsevar = dat$Ci,
                        name = "Ci")
  
  # Transpiration rate - data are gamma
  fit <-
    compile_and_fit_phys_gamma_model(dat,
                                     responsevar = dat$Trmmol,
                                     iter = iter)
  plot_main_effects(fit, name = "Trmmol")
  write_model_statistics(fit, name = "Trmmol")
  plot_posterior_checks(fit,
                        responsevar = dat$Trmmol,
                        name = "Trmmol")
  
  # Vapor pressure deficit - data are normal
  fit <-
    compile_and_fit_phys_normal_model(dat,
                                      responsevar = dat$VpdL,
                                      iter = iter)
  plot_main_effects(fit, name = "VpdL")
  write_model_statistics(fit, name = "VpdL")
  plot_posterior_checks(fit,
                        responsevar = dat$VpdL,
                        name = "VpdL")
  
  # Water use efficiency - data are gamma
  fit <-
    compile_and_fit_phys_gamma_model(dat,
                                     responsevar = dat$iWUE,
                                     iter = iter)
  plot_main_effects(fit, name = "iWUE")
  write_model_statistics(fit, name = "iWUE")
  plot_posterior_checks(fit,
                        responsevar = dat$iWUE,
                        name = "iWUE")
}
