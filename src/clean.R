###########################################################################################
# Data cleanup
# 
# want to control for the method, because of batch and operator effects. One batch was
# 80% methanol first, and the second day was 80% water first. Also want to remove any out-
# liers.
#
###########################################################################################
# Load libraries
library(dplyr)
## note: sva is a bioconductor package, and there may be other dependencies
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("sva"))
library(sva) ## batch correction on raw data

###########################################################################################


clean_design_and_metabolites <- function() {
  
  # Read and clean design
  if (!(file.exists(clean_design_file))) {
    
    # Load design
    metab_design <-
      read.csv(file = design_file, header = T) %>%
      mutate(X = as.character(X)) %>%
      
      # The below sample is an outlier in preliminary pca analysis..
      # removed for all subsequent work
      filter(!(X == "A3_Spco")) %>%
      
      # The value C28_BOGR *I think* should actually be D28
      mutate(block = replace(block, X == "C28_Bogr", "D")) %>%
      mutate(X = replace(X, X == "C28_Bogr", "D28_Bogr"))
    
    # Write clean data
    write.csv(metab_design, file = clean_design_file)
  }
  
  # Read and clean metabolites
  if (!(file.exists(clean_metab_file))) {
    
    # Load metabolites
    metab_data <-
      read.csv(file = metab_file,
               header = T,
               row.names = 1) %>%
      
      # the below sample is an outlier in preliminary pca analysis.. 
      # removed for all subsequent work
      select(-A3_Spco)
    
    # Log transform
    metab_data <-
      as.matrix(log(metab_data), 10)
    
    # Load design
    metab_design <-
      read.csv(file = clean_design_file, header = T)
    
    # Independent vars
    pheno <- metab_design[, 1:5]
    
    # Batch effect column
    batch <- metab_design$batch
    
    # Determine how the design should be
    modcombat <- model.matrix( ~ as.factor(spp) * nitrogen, data = pheno)
    
    # Correct for batch
    combat_metabolite_data <- ComBat(
      dat = metab_data,
      batch = batch,
      mod = modcombat,
      par.prior = FALSE,
      prior.plots = FALSE,
      ref.batch = "normal"
    )
    
    # Rewrite data
    write.csv(combat_metabolite_data, file = clean_metab_file)
  }
  
}