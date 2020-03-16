# Utility functions for project
###########################################################################################

# Data processing

prep_trait_data <- function() {
  phys_data <-
    read.csv("data/clean_SpecAbund_design.csv") %>%
    filter(!(is.na(photo))) %>%
    filter(!(X == "A2_Bogr")) %>% # all zeros, not helpful datapoint
    mutate(iWUE = photo / cond) # generate new variable for water use efficiency
  return(phys_data)
}


# Stan modeling / plotting

plot_main_effects <-  function(fit, name, interaction = T) {
  # What do the main effects look like?
  print("Plotting main effects..")
  # Select posterior draws
  posterior <- as.array(fit)
  dimnames(posterior)$parameters <-
    gsub("B2", "Nitrogen effect", dimnames(posterior)$parameters)
  if (interaction) {
    dimnames(posterior)$parameters <-
      gsub("B1", "Species effect", dimnames(posterior)$parameters)
    dimnames(posterior)$parameters <-
      gsub("B3", "Interaction effect", dimnames(posterior)$parameters)
  }
  
  # Plot species, nitrogen, and interaction effect
  if (interaction){
    pars = c("Species effect",
             "Nitrogen effect",
             "Interaction effect")
  } else {
    pars = c("Nitrogen effect")
  }
  pars <- 
  param_plot <-
    mcmc_dens(
      posterior,
      pars = pars,
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


write_model_statistics <-
  function(fit,
           name,
           mixture = F,
           interaction = T) {
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
    if (mixture & interaction) {
      Pr_vals <-
        rbind(model_stats["intercept[1]", ],
              model_stats["intercept[2]", ],
              model_stats["B1", ],
              model_stats["B2", ],
              model_stats["B3", ])
    } else if (!(interaction)) {
      if (mixture) {
        Pr_vals <-
          rbind(model_stats["intercept[1]",],
                model_stats["intercept[2]",],
                model_stats["B2",])
      } else {
        Pr_vals <- rbind(model_stats["B0",],
                         model_stats["B2",])
      }
    } else {
      Pr_vals <-
        rbind(model_stats["B0",],
              model_stats["B1",],
              model_stats["B2",],
              model_stats["B3",])
    }
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
    if (mixture & interaction) {
      colnames(Pr_calc) <-
        c("intercept[1]", "intercept[2]", "B1", "B2", "B3")
    } else if (!(interaction)) {
      if (mixture) {
        colnames(Pr_calc) <- c("intercept[1]", "intercept[2]", "B2")
      } else {
        colnames(Pr_calc) <- c("B0", "B2")
      }
    } else {
      colnames(Pr_calc) <- c("B0", "B1", "B2", "B3")
    }
    write.csv(Pr_calc,
              file = paste("output/stanmodel_results_", name, "_pr.csv", sep = ""))
  }


plot_posterior_checks <- function(fit, responsevar, name) {
  # Posterior predictive checks
  print("Plotting posterior predictive checks")
  list_of_draws <- rstan::extract(fit)
  yrep <- list_of_draws$draws1
  # Two ways to plot actual and modeled data
  np1 <- ppc_dens_overlay(responsevar, yrep[1:400, ]) +
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


# Plotting

theme_sigmaplot <-
  function(xticks = TRUE,
           ticklen = -0.25) {
    # This function adds Sigma-plot like theme elements to a ggplot object.
    # Use as an additional arg, eg:
    # ggplot() + theme_sigmaplot()
    
    sigmaplot <-
      theme(
        panel.background = element_blank(),
        panel.border = element_rect(size = 1, fill = NA),
        legend.key = element_rect(fill = NA),
        axis.ticks.length.y = unit(ticklen, "cm"),
        axis.ticks.length.y.right = unit(ticklen, "cm"),
        axis.ticks.length.x = unit(ticklen, "cm"),
        axis.text.x = element_text(
          color = "black",
          margin = margin(
            t = 10,
            r = 0,
            b = 0,
            l = 0,
            unit = "pt"
          )
        ),
        axis.text.y = element_text(
          hjust = 1,
          color = "black",
          margin = margin(
            t = 0,
            r = 10,
            b = 0,
            l = 0,
            unit = "pt"
          )
        ),
      )
    if (!xticks) {
      sigmaplot <- sigmaplot +
        theme(axis.ticks.x = element_blank())
    }
    return(sigmaplot)
  }


legend_custom <-
  function() {
    obj <- theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = 0.5,
      legend.title = element_blank()
    )
    
    return(obj)
  }
