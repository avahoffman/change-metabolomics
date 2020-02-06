# Config file for project

###########################################################################################

# PCA plotting
pca_component <- "PC1"
pca_yaxis <- "PC1 (59%)"

###########################################################################################

# MCMC sampling
iter <- 1000

###########################################################################################

# Plotting

bogr_color <- "#1b9e77"
spco_color <- "#d95f02"
elel_color <- "#7570b3"

bogr_shape <- 16
spco_shape <- 8
elel_shape <- 3

B.gracilis <- expression(italic("B. gracilis"))
E.elymoides <- expression(italic("E. elymoides"))
S.coccinea <- expression(italic("S. coccinea"))
# B.gracilis <- expression(italic("Grass species       "))
# S.coccinea <- expression(italic("Wildflower species"))

photo_ylab <- expression(paste("Net photsynthetic rate (", mu,"mol ",CO[2]~ m^-2, s^-1,")"))
cond_ylab <- expression(paste("Stomatal conductance (mol ", H[2]~O~ m^-2, s^-1,")"))
ci_ylab <- expression(paste("Intercellular ", CO[2]," (", mu,"mol ",CO[2]~ mol^-1,")"))
iWUE_ylab <- expression(paste("Water use efficiency (", mu,"mol ",CO[2]," / ","mmol ",H[2]~O,")"))