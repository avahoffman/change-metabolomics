###########################################################################################
# Execution script for analysis run on the computing cluster

###########################################################################################
# Load sources
# General functions and configuration
source("MIC_score.R")

# The following should be in the executing directory:
# This file
# "data/clean_SpecAbund_metabolites.csv"
# "data/phys_metabs_for_modules.csv"
# MIC_score.R
# MIC_exe.sh

###########################################################################################
# Run 

run_metab_data_MIC(spp = 0)
run_metab_data_MIC(spp = 1)
