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

bogr_shape <- 16
spco_shape <- 8
elel_shape <- 3

B.gracilis <- expression(italic("B. gracilis"))
E.elymoides <- expression(italic("E. elymoides"))
S.coccinea <- expression(italic("S. coccinea"))

photo_ylab <- expression(paste("Net photsynthetic rate (", mu,"mol ",CO[2]~ m^-2, s^-1,")"))
cond_ylab <- expression(paste("Stomatal conductance (mol ", H[2]~O~ m^-2, s^-1,")"))
ci_ylab <- expression(paste("Intercellular ", CO[2]," (", mu,"mol ",CO[2]~ mol^-1,")"))
iWUE_ylab <- expression(paste("Water use efficiency (", mu,"mol ",CO[2]," / ","mol ", H[2]~O~ m^-2, s^-1,")"))

heatmap_plotnames <- 
  c(
    "nitrogen",
    "species",
    "BOGR cover",
    "SPCO cover",
    "ELEL cover",
    "total biomass",
    "grass biomass",
    "forb biomass",
    "photo. rate",
    "conductance",
    "intercell. CO2",
    "water use efficiency"
  )

heatmap_colorpal <- 
  colorRampPalette(c("#38598CFF","white","#92D741FF"))

sem_pos_color <- "#92D741FF"
sem_neg_color <- "#38598CFF"