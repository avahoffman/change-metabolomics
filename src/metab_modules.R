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
library(WGCNA)

###########################################################################################


wgcna_data <- function(){
  dat <- 
    read.csv("data/clean_SpecAbund_metabolites.csv")
  phys_dat <-  
    read.csv("data/phys_metabs_for_modules.csv")
  
  metabs <- as.data.frame(t(dat[, -c(1)]))
  
  ## transpose (so all metabolites as columns and samples as rows). also get rid of any extra data (just want samples)
  names(metabs) <- dat$X
  rownames(metabs) = names(dat)[-c(1)]
  # gsg = goodSamplesGenes(metabs, verbose = 3)
  # n_metabs = ncol(metabs)
  # n_samples = nrow(metabs)
  ## remove columns that hold information we do not need.
  all_phys = phys_dat[, -c(1, 4, 6)]
  samples = rownames(metabs)
  phys_rows = match(samples, all_phys$Samples)
  dat_phys = all_phys[phys_rows, -1]
  rownames(dat_phys) = all_phys[phys_rows, 1]
  collectGarbage()

  return(list(metabs, samples, dat_phys))
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


make_sample_dendro_heatmap(wgcna_data()[[1]], wgcna_data()[[3]])

