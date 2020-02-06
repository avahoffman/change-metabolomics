# Execution script for project pipeline.

###########################################################################################
# Set working directory for the repository (should be the git repo):
wd <-
   "/Users/hoffman ava/change-metabolomics/"
# "/Users/avahoffman/Dropbox/Research/SGS_Metabolomics/Change_2016-17_Nate/change-metabolomics/"

setwd(wd)

# General functions and configuration
source("src/config.R")
source("src/utils.R")

# Specific functions
source("src/clean.R")
source("src/pca.R")
source("src/plot_traits.R")
source("src/plot_cover.R")
source("src/model_on_pc.R")
source("src/model_on_phys.R")

###########################################################################################
# Run pipeline

# Clean data
clean_design_and_metabolites()

# Run PCA and plot
plot_pca(run_pca(), filename = "figures/component_v_N.pdf")

# Bayes model for PCs
model_first_three_pcs()

# Bayes models for phys
model_phys()

# Phys plots against N
arrange_phys_plots()

# Cover against N
plot_nitrogen_and_cover(filename = "figures/cover_v_N.pdf")
