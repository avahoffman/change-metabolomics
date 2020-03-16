###########################################################################################
# Execution script for project pipeline

###########################################################################################
# Set working directory for the repository (should be the git repo)
wd <-
   "/Users/avahoffman/Dropbox/Research/SGS_Metabolomics/Change_2016-17_Nate/change-metabolomics/"
setwd(wd)

###########################################################################################
# Load sources
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
source("src/model_on_cover.R")
source("src/metab_modules.R")
source("src/sem.R")

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

# Bayes models for cover
model_cover()

# Phys plots against N
arrange_phys_plots()

# Cover against N
plot_nitrogen_and_cover(filename = "figures/cover_v_N.png")

# Metabolite network/module analysis
make_module_heatmap()

# SEM
plot_sem()
