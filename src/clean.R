##
## DATA CLEANUP
##
## want to control for the method, because of batch and operator effects. One batch was
## 80% methanol first, and the second day was 80% water first. Also want to remove any out-
## liers.
##
###########################################################################################
###########################################################################################
# Load libraries
library(dplyr)
## note: sva is a bioconductor package, and there may be other dependencies
# source("https://bioconductor.org/biocLite.R")
# biocLite("sva")
library(sva) ## batch correction on raw data


###########################################################################################


clean_design_and_metabolites <- function() {
  # Read and clean design
  if (!(file.exists("data/clean_SpecAbund_design.csv"))) {
    metab_design <-
      read.csv(file = "data/SpecAbund_design.csv", header = T) %>%
      mutate(X = as.character(X)) %>%
      # the below sample is an outlier in preliminary pca analysis.. removed for all subsequent work
      filter(!(X == "A3_Spco")) %>%
      # The value C28_BOGR *I think* should actually be D28
      mutate(block = replace(block, X == "C28_Bogr", "D")) %>%
      mutate(X = replace(X, X == "C28_Bogr", "D28_Bogr"))
    write.csv(metab_design, file = "data/clean_SpecAbund_design.csv")
  }
  
  # Read and clean metabolites
  if (!(file.exists("data/clean_SpecAbund_metabolites.csv"))) {
    # Load metabolites
    metab_data <-
      read.csv(file = "data/SpecAbund_metabolites.csv",
               header = T,
               row.names = 1) %>%
      # the below sample is an outlier in preliminary pca analysis.. removed for all subsequent work
      select(-A3_Spco)
    
    # Log transform
    metab_data <-
      as.matrix(log(metab_data), 10)
    
    # Load design
    metab_design <-
      read.csv(file = "data/clean_SpecAbund_design.csv", header = T)
    # Independent vars
    pheno <- metab_design[, 1:5]
    # Batch effect column
    batch <- metab_design$batch
    
    # mod <-
    #   model.matrix(~ as.factor(spp) * as.factor(nitrogen) + as.factor(block), data =
    #                  pheno)
    # mod0 <- model.matrix(~ 1, data = pheno)
    modcombat <- model.matrix( ~ as.factor(spp) * nitrogen, data = pheno)
    combat_metabolite_data <- ComBat(
      dat = metab_data,
      batch = batch,
      mod = modcombat,
      par.prior = FALSE,
      prior.plots = FALSE,
      ref.batch = "normal"
    )
    write.csv(combat_metabolite_data, file = "data/clean_SpecAbund_metabolites.csv")
  }
  
}