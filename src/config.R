# Config file for project

###########################################################################################

# PCA plotting
pca_component <- "PC1"
pca_yaxis <- "PC1 (40%)"

###########################################################################################

# MCMC sampling
iter <- 100000

###########################################################################################

# Plotting

# Follow viridis plotting convention
bogr_color <- "#73D055FF"
spco_color <- "#482677FF"
elel_color <- "#287D8EFF"

# Point shapes for different species
bogr_shape <- 16
spco_shape <- 8
elel_shape <- 3

# Species names
B.gracilis <- expression(italic("B. gracilis"))
E.elymoides <- expression(italic("E. elymoides"))
S.coccinea <- expression(italic("S. coccinea"))

# Labels for physiology plot
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

# Color palette for the SEM plot
sem_pos_color <- "#92D741FF"
sem_neg_color <- "#38598CFF"