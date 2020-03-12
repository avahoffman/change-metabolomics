##
## Make principal components
##
###########################################################################################
###########################################################################################
# Load libraries
library(dplyr)
library(ggplot2)
library(cowplot)

###########################################################################################


run_pca <- function() {
  metab_design <-
    read.csv(file = "data/clean_SpecAbund_design.csv", header = T)
  metab_data <-
    as.matrix(read.csv(
      file = "data/clean_SpecAbund_metabolites.csv",
      header = T,
      row.names = 1
    ))
  
  metabolites <- t(metab_data)
  prcomp_analysis <-
    prcomp(metabolites, center = T, scale = T)
  
  # Save loadings
  write.csv(prcomp_analysis$rotation,
            file = "output/pca_loadings.csv")
  
  # Save summary of variance explained
  write.csv(summary(prcomp_analysis)$importance,
            file = "output/pca_proportionofvariance.csv")
  
  # Save scores
  write.csv(prcomp_analysis$x,
            file = "output/pca_scores.csv")
  
  combined_data <-
    full_join(
      metab_design %>%
        rename(sample_name = X) %>%
        mutate(sample_name = as.character(sample_name)) %>%
        select(-(X.1)),
      as.data.frame(prcomp_analysis$x) %>%
        mutate(sample_name = rownames(as.data.frame(
          prcomp_analysis$x
        ))),
      by = "sample_name"
    )
  
  write.csv(combined_data,
            file = "output/pca_scores_predictors.csv")
  
  return(combined_data)
}


## plot of PC 1

plot_pca <- function(combined_data, filename = NA) {
  gg <-
    ggplot() +
    geom_point(data = combined_data %>%
                 filter(!(is.na(
                   get(pca_component)
                 ))),
               aes(
                 x = nitrogen,
                 y = get(pca_component),
                 col = as.factor(spp),
                 shape = as.factor(spp)
               )) +
    geom_jitter(size = 3,
                width = 0.2,
                height = 0) +
    scale_shape_manual(values = c(bogr_shape, spco_shape),
                       guide = F) +
    theme_cowplot() +
    scale_color_manual(
      values = c(bogr_color,
                 spco_color),
      labels = c(B.gracilis, S.coccinea)
    ) +
    xlab(expression(paste("N addition (g ", m ^ {
      -2
    }, ')'))) +
    ylab(pca_yaxis) +
    guides(colour = guide_legend(override.aes = list(shape = c(
      bogr_shape, spco_shape
    )))) +
    legend_custom() +
    theme(legend.text.align = 0)
  
  gg
  
  if (!(is.na(filename))) {
    ggsave(file = filename,
           height = 3.5,
           width = 3.5)
  }
  
  return(gg)
}
