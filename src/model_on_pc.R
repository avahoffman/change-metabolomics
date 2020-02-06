##
## effects SPECIES AND NITROGEN on metabolome PCs
## Keep in mind, recommended lawn application is about 5g/m2 of Nitrogen!
##
###########################################################################################
###########################################################################################
# Load libraries
library(rstan) # Bayesian model compiler and sampler
options(mc.cores = parallel::detectCores()) # option to make Stan parallelize
library(bayesplot) # Plots of mcmc output


###########################################################################################


check_pc_distributions <- function() {
  # Check histograms to determine how data is distributed
  metab_data <-
    read.csv(file = "output/pca_scores_predictors.csv", header = T)
  hist(metab_data$PC1)
  hist(metab_data$PC2)
  hist(metab_data$PC3)
}


compile_and_fit_pc_normal_model <-
  function(df, responsevar, iter = 100000) {
    # Stan model
    pc_model_normal <- "
    data {
    int<lower=0> N; //list of data
    vector[N] spp;
    vector[N] Nit;
    vector[N] PC;
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
    PC ~ normal(mu,sigma);
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
    comp.normal <- stan_model(model_code = pc_model_normal,
                              model_name = 'pc.model.normal')
    # Select data
    modeldat <-
      list(
        'N' = nrow(df),
        'PC' = responsevar,
        'spp' = as.numeric(df$spp) - 1, # Want zeros and ones for species
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


model_first_three_pcs <- function() {
  # Fit models on three principal components:
  # Load data
  dat <-
    read.csv(file = "output/pca_scores_predictors.csv", header = T)
  
  # PC1
  fit <-
    compile_and_fit_pc_normal_model(dat, responsevar = dat$PC1)
  plot_main_effects(fit, name = "PC1")
  write_model_statistics(fit, name = "PC1")
  plot_posterior_checks(fit, responsevar = dat$PC1, name = "PC1")
  
  # PC2
  fit <-
    compile_and_fit_pc_normal_model(dat, responsevar = dat$PC2)
  plot_main_effects(fit, name = "PC2")
  write_model_statistics(fit, name = "PC2")
  plot_posterior_checks(fit, responsevar = dat$PC2, name = "PC2")
  
  # PC3
  fit <-
    compile_and_fit_pc_normal_model(dat, responsevar = dat$PC3)
  plot_main_effects(fit, name = "PC3")
  write_model_statistics(fit, name = "PC3")
  plot_posterior_checks(fit, responsevar = dat$PC3, name = "PC3")
  
}
