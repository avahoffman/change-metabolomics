###########################################################################################
##
## clustering using WGCNA
##
###########################################################################################
###########################################################################################
# Load libraries
## Several dependencies required for WCGNA prior to install
# install.packages("nlme")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("impute","preprocessCore", "GO.db", "AnnotationDbi"))
# install.packages("WGCNA")
# Note that there are some default compiler changes in Mac OS Catalina. These can be worked
# around using R from terminal and xcode command line tools.
library(WGCNA)

###########################################################################################


wgcna_data <- function() {
  dat <-
    read.csv("data/clean_SpecAbund_metabolites.csv")
  phys_dat <-
    read.csv("data/phys_metabs_for_modules.csv")
  
  metabs <- as.data.frame(t(dat[,-c(1)]))
  
  ## transpose (so all metabolites as columns and samples as rows). also get rid of any extra data (just want samples)
  names(metabs) <- dat$X
  rownames(metabs) = names(dat)[-c(1)]
  # gsg = goodSamplesGenes(metabs, verbose = 3)
  # n_metabs = ncol(metabs)
  # n_samples = nrow(metabs)
  ## remove columns that hold information we do not need.
  all_phys = phys_dat[,-c(1, 4, 6)]
  samples = rownames(metabs)
  phys_rows = match(samples, all_phys$Samples)
  dat_phys = all_phys[phys_rows,-1]
  rownames(dat_phys) = all_phys[phys_rows, 1]
  collectGarbage()
  
  return(list(metabs, dat_phys))
}


make_sample_dendro_heatmap <-
  function(metabs, dat_phys) {
    ## Re-cluster samples
    sample_tree <- hclust(dist(metabs), method = "average")
    
    ## ensure no predictor vars are text factors, need to convert to 0s and 1s
    phys_colors <- numbers2colors(dat_phys, signed = FALSE)
    
    # Plot the sample dendrogram and the colors underneath.
    pdf(file = "figures/metab_dendro_with_predictors.pdf",
        width = 10,
        height = 5)
    plotDendroAndColors(sample_tree,
                        phys_colors,
                        groupLabels = names(dat_phys),
                        main = "Sample dendrogram and trait heatmap")
    dev.off()
  }


determine_power_for_modules <-
  function(metabs) {
    # Determine what power to use for the module analysis
    
    powers <-
      c(c(1:10),
        seq(from = 2, to = 15, by = 2))
    # Call the network topology analysis function (may need to use base R)
    sft <-
      pickSoftThreshold(metabs,
                        powerVector = powers,
                        verbose = 5)
    
    # Plot the results:
    # Plot the sample dendrogram and the colors underneath.
    pdf(file = "figures/power_threshold_modules.pdf",
        width = 8,
        height = 5)
    par(mfrow = c(1, 2)) # Plots beside each other
    cex <- 0.7 # Text size
    
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(
      sft$fitIndices[, 1],-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
      xlab = "Soft Threshold (power)",
      ylab = "Scale Free Topology Model Fit,signed R^2",
      type = "n",
      main = paste("Scale independence")
    )
    
    text(
      sft$fitIndices[, 1],
      -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
      labels = powers,
      cex = cex,
      col = "red"
    )
    
    # this line corresponds to using an R^2 cut-off of h
    abline(h = 0.90, col = "red")
    # Mean connectivity as a function of the soft-thresholding power
    plot(
      sft$fitIndices[, 1],
      sft$fitIndices[, 5],
      xlab = "Soft Threshold (power)",
      ylab = "Mean Connectivity",
      type = "n",
      main = paste("Mean connectivity")
    )
    text(
      sft$fitIndices[, 1],
      sft$fitIndices[, 5],
      labels = powers,
      cex = cex,
      col = "red"
    )
    dev.off()
  }


make_module_network <-
  function(metabs, makeplot = F) {
    net = blockwiseModules(
      metabs,
      power = 8,
      TOMType = "unsigned",
      minModuleSize = 20,
      reassignThreshold = 0,
      mergeCutHeight = 0.25,
      numericLabels = TRUE,
      pamRespectsDendro = FALSE,
      saveTOMs = FALSE,
      verbose = 3
    )
    # Shows module rank and size (num. genes) - module "0" is genes outside of modules
    print(table(net$colors))
    if (makeplot) {
      # open a graphics window
      sizeGrWindow(12, 9)
      # Convert labels to colors for plotting
      mergedColors <-  labels2colors(net$colors)
      # Plot the dendrogram and the module colors underneath
      pdf(file = "Dendrogram_modules.pdf",
          width = 9,
          height = 7)
      plotDendroAndColors(
        net$dendrograms[[1]],
        mergedColors[net$blockGenes[[1]]],
        "Module colors",
        dendroLabels = FALSE,
        hang = 0.03,
        addGuide = TRUE,
        guideHang = 0.05
      )
      dev.off()
    }
    
    return(net)
  }


make_module_heatmap <-
  function(metabs, dat_phys, net) {
    # Use functions above
    metabs <- wgcna_data()[[1]]
    dat_phys <- wgcna_data()[[2]]
    net <- make_module_network(wgcna_data()[[1]])
    
    # Save raw module assignments for each metabolite
    write.csv(cbind(
      labels2colors(net$colors),
      read.csv("data/clean_SpecAbund_metabolites.csv")
    ),
    file = "output/module_assignments_raw.csv")
    
    MEs <-
      orderMEs(moduleEigengenes(metabs, labels2colors(net$colors))$eigengenes)
    # Gathere correlation and p values
    moduleTraitCor <- cor(MEs, dat_phys, use = "p")
    moduleTraitPvalue <-
      p.adjust(corPvalueStudent(moduleTraitCor, nrow(metabs)), "BH")
    # Will display correlations and their p-values
    textMatrix <-  paste(signif(moduleTraitCor, 2),
                         "\n(",
                         signif(moduleTraitPvalue, 1),
                         ")",
                         sep = "")
    
    # Make heatmap plot
    sizeGrWindow(10, 6)
    par(mar = c(6, 8.5, 3, 3))
    ## Display the correlation values within a heatmap plot
    pdf(file = "figures/module_phys_heatmap.pdf",
        width = 8,
        height = 6)
    labeledHeatmap(
      Matrix = moduleTraitCor,
      xLabels = heatmap_plotnames,
      yLabels = c(1:16),
      ##change this back to names(MEs) for color naming convention instead
      ySymbols = names(MEs),
      colorLabels = F,
      colors = heatmap_colorpal(50),
      textMatrix = textMatrix,
      setStdMargins = T,
      cex.text = 0.5,
      #main = paste("Module-trait relationships"),
      zlim = c(-1, 1)
    )
    dev.off()
  }
