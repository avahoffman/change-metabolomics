###########################################################################################
# Execution script for project pipeline

###########################################################################################
# Set working directory for the repository (should be the git repo)
wd <-
  here() # Replace with manual git repo path if desired, e.g. 
         # "/Users/me/change-metabolomics/"

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
source("src/model_on_pc.R")
source("src/model_on_phys.R")
source("src/metab_modules.R")
source("src/MIC_plots.R")

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

# Metabolite network/module analysis
make_module_heatmap()

# Nitrogen MIC curves
gather_nit_plots(filename = "figures/MIC.pdf", span_ = 0.75)

# Make a nice combined figure for pub
png("figures/PCA_MIC.png",
    height = 400,
    width = 800)
grid.arrange(
  plot_grid(
    plot_pca(run_pca()),
    nrow = 1,
    labels = c("a"),
    label_size = 18
  ),
  gather_nit_plots(span_ = 0.75)[[1]],
  gather_nit_plots(span_ = 0.75)[[2]],
  nrow = 2,
  layout_matrix =
    rbind(c(1, 2),
          c(1, 3)),
  heights = c(100,11),
  widths = c(1,2)
)
dev.off()
