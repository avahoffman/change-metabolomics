# Config file for project

###########################################################################################

# Files

design_file <- "data/SpecAbund_design.csv"
clean_design_file <- "data/clean_SpecAbund_design.csv"
metab_file <- "data/SpecAbund_metabolites.csv"
clean_metab_file <- "data/clean_SpecAbund_metabolites.csv"

# PCA plotting

# Which PC component should be used? Choose PC1, PC2, or PC3
pca_component <- "PC1"

# What should the y axis label be?
pca_yaxis <- "PC1 (40%)"

###########################################################################################

# MCMC sampling

# How many iterations to use? The functions that use this variable will discard the first 
# half of iter as warmup.
iter <- 100000

###########################################################################################

# Plotting

# Follow viridis plotting convention, but these can be any hex color.
bogr_color <- "#73D055FF"
spco_color <- "#482677FF"
elel_color <- "#287D8EFF"

# Point shapes for different species
bogr_shape <- 16
spco_shape <- 8
elel_shape <- 3

# Species names - ensures these are italicized throughout
B.gracilis <- expression(italic("B. gracilis"))
E.elymoides <- expression(italic("E. elymoides"))
S.coccinea <- expression(italic("S. coccinea"))

# Labels for physiology plot y axes
photo_ylab <-
  expression(atop("Net photsynthetic rate",
          "("~mu~"mol "~CO[2] ~ m ^ -2~s ^ -1~ ")"))
cond_ylab <- 
  expression(atop("Stomatal conductance", 
                  "(mol "~H[2]~O~ m^-2~s^-1~")"))
ci_ylab <- 
  expression(atop("Intercellular "~CO[2],
                  "("~mu~"mol "~CO[2]~ mol^-1~")"))
iWUE_ylab <- 
  expression(atop("Water use efficiency",
                  "("~mu~"mol"~CO[2]~" / "~"mol"~ H[2]~O~m^-2~s^-1~")"))

# Labels for the heatmap (WGCNA modules)
heatmap_plotnames <- 
  c(
    "nitrogen",
    "species",
    "photo. rate",
    "conductance",
    "intercell. CO2",
    "water use efficiency"
  )

# Color palette for the heatmap (currently two colors, with white in between)
heatmap_colorpal <- 
  colorRampPalette(c("#38598CFF","white","#92D741FF"))

# Metabolites from MIC also in module 16 - Currently this needs to be determined manually
# by looking at metabolites with significant MIC output and metabolites found in module 16.
# TODO: Make this automatic
overlap_metab <- as.factor(c("V224","V713"))
