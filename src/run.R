# Execution script for project pipeline.

###########################################################################################
# Set working directory for the repository (should be the git repo):
wd <-
  "/Users/hoffman ava/change-metabolomics/"

setwd(wd)

# General functions and configuration
source("src/config.R")
source("src/utils.R")

# Specific functions
source("src/clean.R")
source("src/pca.R")

###########################################################################################
# Run pipeline

clean_design_and_metabolites()
plot_pca(run_pca(),
         filename = "figures/component_v_N.pdf")
