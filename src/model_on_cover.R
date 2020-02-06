##
## effects of NITROGEN on community
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


compile_and_fit_cover_mixture_model <-
  function(df, responsevar, iter = 100000) {
    ## Bayesian mixture model for bimodal distributions
    
    # Stan model
    community_model_mixture <- "
      data {
      int<lower=1> K; // number of peaks
      int<lower=1> N; // list of data
      real comm[N]; // observations
      real Nit[N];
      }
      parameters {
      real intercept[K]; // locations of mixture components
      real<lower=0> sigma[K];  // scales of mixture components
      real B2;
      }
      model {
      real ps[K]; // temp for log component densities
      intercept ~ normal(0, 30);
      B2 ~ normal(0,100);
      for (n in 1:N) {
      for (k in 1:K) {
      ps[k] = normal_lpdf(comm[n] | intercept[k] + Nit[n]*B2 , sigma[k]);
      }
      target += log_sum_exp(ps);
      }
      }
      generated quantities {
      real draws1[N];
      for (n in 1:N) {
      for (k in 1:K) {
      draws1[n] = normal_rng( intercept[k] + Nit[n]*B2 , sigma[k]);
      }
      }
      }
    "
    
    # Compile model
    comp.mixture <-
      stan_model(model_code = community_model_mixture,
                 model_name = 'community.model.mixture')
    
    # Select data
    modeldat <- list(
      'N' = nrow(df),
      'K' = 2,
      # Number of peaks in multimodal distribution
      'comm' = responsevar,
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


compile_and_fit_cover_gamma_model <-
  function(df, responsevar, iter = 100000) {
    ## Bayesian model for gamma distributions
    
    # Stan model
    community_model_gamma <- "
      data {
      int<lower=0> N; //list of data
      vector[N] Nit;
      vector[N] comm;
      }
      parameters {
      real B0;
      real B2;
      real<lower=0> phi;
      }
      transformed parameters{
      vector[N] mu;
      vector[N] alpha;
      vector[N] beta;
      for (n in 1:N)
      mu[n] = exp( B0 + Nit[n]*B2 );
      alpha = mu .* mu / phi;
      beta = mu / phi;
      }
      model {
      phi ~ cauchy(0, 10);
      B0 ~ normal(0,1000);
      B2 ~ normal(0,1000);
      comm ~ gamma(alpha,beta);
      }
      generated quantities{
      vector[N] draws1;
      for(n in 1:N)
      draws1[n] = gamma_rng(alpha[n],beta[n]); //posterior draws
      }
    "
    
    # Compile model
    comp.gamma <- stan_model(model_code = community_model_gamma,
                             model_name = 'community.model.gamma')
    
    # Select data
    modeldat <- list('N' = nrow(df),
                     'comm' = responsevar,
                     'Nit' = df$nitrogen)
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


model_cover <- function() {
  # Load data
  dat <-
    na.omit(read.csv("data/SpecAbund_community_design.csv"))
  
  # BOGR cover
  fit <-
    compile_and_fit_cover_gamma_model(dat,
                                      responsevar = dat$BOGR_cover, 
                                      iter = iter)
  plot_main_effects(fit, name = "BOGR_cover", interaction = F)
  write_model_statistics(fit,
                         name = "BOGR_cover",
                         mixture = F,
                         interaction = F)
  plot_posterior_checks(fit, 
                        responsevar = dat$BOGR_cover, 
                        name = "BOGR_cover")
  
  # SPCO cover
  fit <-
    compile_and_fit_cover_gamma_model(dat,
                                      responsevar = dat$SPCO_cover, 
                                      iter = iter)
  plot_main_effects(fit, name = "SPCO_cover", interaction = F)
  write_model_statistics(fit,
                         name = "SPCO_cover",
                         mixture = F,
                         interaction = F)
  plot_posterior_checks(fit, 
                        responsevar = dat$SPCO_cover, 
                        name = "SPCO_cover")
  
  # ELEL cover
  fit <-
    compile_and_fit_cover_gamma_model(dat,
                                      responsevar = dat$ELEL_.cover, 
                                      iter = iter)
  plot_main_effects(fit, name = "ELEL_.cover", interaction = F)
  write_model_statistics(fit,
                         name = "ELEL_.cover",
                         mixture = F,
                         interaction = F)
  plot_posterior_checks(fit, 
                        responsevar = dat$ELEL_.cover, 
                        name = "ELEL_.cover")
  
  # Grass biomass
  fit <-
    compile_and_fit_cover_gamma_model(dat,
                                      responsevar = dat$grass_bio, 
                                      iter = iter)
  plot_main_effects(fit, name = "grass_bio", interaction = F)
  write_model_statistics(fit,
                         name = "grass_bio",
                         mixture = F,
                         interaction = F)
  plot_posterior_checks(fit, 
                        responsevar = dat$grass_bio, 
                        name = "grass_bio")
  
  # Forb biomass
  fit <-
    compile_and_fit_cover_gamma_model(dat,
                                      responsevar = dat$forb_bio, 
                                      iter = iter)
  plot_main_effects(fit, name = "forb_bio", interaction = F)
  write_model_statistics(fit,
                         name = "forb_bio",
                         mixture = F,
                         interaction = F)
  plot_posterior_checks(fit, 
                        responsevar = dat$forb_bio, 
                        name = "forb_bio")
  
  # Total biomass
  fit <-
    compile_and_fit_cover_gamma_model(dat,
                                      responsevar = dat$total_bio, 
                                      iter = iter)
  plot_main_effects(fit, name = "total_bio", interaction = F)
  write_model_statistics(fit,
                         name = "total_bio",
                         mixture = F,
                         interaction = F)
  plot_posterior_checks(fit, 
                        responsevar = dat$total_bio, 
                        name = "total_bio")
}
