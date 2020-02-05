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


compile_and_fit_normal_model <-
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
        'spp' = as.numeric(df$spp) - 1,
        ## want zeros and ones for species
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


plot_main_effects <-  function(fit, name) {
  # What do the main effects look like?
  print("Plotting main effects..")
  # Select posterior draws
  posterior <- as.array(fit)
  dimnames(posterior)$parameters[2] <- "Species effect"
  dimnames(posterior)$parameters[3] <- "Nitrogen effect"
  dimnames(posterior)$parameters[4] <- "Interaction effect"
  
  # Plot species, nitrogen, and interaction effect
  param_plot <-
    mcmc_dens(
      posterior,
      pars = c("Species effect",
               "Nitrogen effect",
               "Interaction effect"),
      facet_args = list(ncol = 1)
    ) +
    geom_vline(xintercept = 0, lty = 3) +
    theme_minimal(base_size = 20) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle(name)
  # Save plot
  ggsave(param_plot,
         file = paste("output/stanmodel_results_",
                      name,
                      "_maineffects.pdf",
                      sep = ""))
}


write_model_statistics <- function(fit, name) {
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
    rbind(model_stats["B0",], 
          model_stats["B1",], 
          model_stats["B2",], 
          model_stats["B3",])
  Pr_vals <- Pr_vals[,-4] # Drop rhat for calculating Pr
  Pr_calc <- as.data.frame(NULL)
  # Iterate through parameters
  for (r in 1:nrow(Pr_vals)) {
    max_ <- max(Pr_vals[r,])
    min_ <- min(Pr_vals[r,])
    # Calculate the proportion of overlap on zero
    if (max_ < 0 & min_ < 0) {
      pr. <- 0
    } else if (max_ > 0 & min_ > 0) {
      pr. <- 0
    } else if (abs(max_) > abs(min_)) {
      pr. <- abs(min_) / (abs(min_) + abs(max_))
    } else if (abs(max_) < abs(min_)) {
      pr. <- abs(max_) / (abs(min_) + abs(max_))
    } else {
      pr. <- "NA"
    }
    Pr_calc[1, r] <- pr.
  }
  colnames(Pr_calc) <- c("B0", "B1", "B2", "B3")
  write.csv(Pr_calc,
            file = paste("output/stanmodel_results_", name, "_pr.csv", sep = ""))
}


plot_posterior_checks <- function(fit, responsevar, name) {
  # Posterior predictive checks
  print("Plotting posterior predictive checks")
  list_of_draws <- extract(fit)
  yrep <- list_of_draws$draws1
  # Two ways to plot actual and modeled data
  np1 <- ppc_dens_overlay(responsevar, yrep[1:400,]) +
    theme_minimal(base_size = 20) +
    xlab(name)
  ggsave(np1,
         file = paste(
           "output/stanmodel_results_",
           name,
           "_predictivecheck1.pdf",
           sep = ""
         ))
  pdf(file = paste(
    "output/stanmodel_results_",
    name,
    "_predictivecheck2.pdf",
    sep = ""
  ))
  hist(responsevar,
       prob = T,
       breaks = 10,
       main = name)
  lines(density(list_of_draws$draws1), col = "red")
  dev.off()
}


model_first_three_pcs <- function() {
  # Fit models on three principal components:
  # Load data
  dat <-
    read.csv(file = "output/pca_scores_predictors.csv", header = T)
  
  # PC1
  fit <-
    compile_and_fit_normal_model(dat, responsevar = dat$PC1)
  plot_main_effects(fit, name = "PC1")
  write_model_statistics(fit, name = "PC1")
  plot_posterior_checks(fit, responsevar = dat$PC1, name = "PC1")
  
  # PC2
  fit <-
    compile_and_fit_normal_model(dat, responsevar = dat$PC2)
  plot_main_effects(fit, name = "PC1")
  write_model_statistics(fit, name = "PC2")
  plot_posterior_checks(fit, responsevar = dat$PC2, name = "PC2")
  
  # PC3
  fit <-
    compile_and_fit_normal_model(dat, responsevar = dat$PC3)
  plot_main_effects(fit, name = "PC1")
  write_model_statistics(fit, name = "PC3")
  plot_posterior_checks(fit, responsevar = dat$PC3, name = "PC3")
  
}
