#!/bin/bash


#SBATCH --time=72:00:00
#SBATCH --job-name=MIC_scor
#SBATCH --mail-type=ALL
#SBATCH --mail-user=avamariehoffman@gmail.com
#SBATCH --error=error.%j.out
#SBATCH --output=out.%j.out
#SBATCH --qos=normal
#SBATCH --partition=smem
#SBATCH --ntasks=48


# Written by:    Ava Hoffman
# Date:              30 Sept 2020
# Purpose:          Calculate MIC scores for relationships between metabolites and something else

# purge all existing modules
module purge

# load any modules needed to run your program
module load R

# Set paths

export R_LIBS=/home/cdmckee@colostate.edu/R/x86_64-pc-linux-gnu-library/3.5
export PATH=$PATH:/home/cdmckee@colostate.edu/R/x86_64-pc-linux-gnu-library/3.5

#############

Rscript run_summit.R

#############

echo "== End of Job =="