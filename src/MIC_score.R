###########################################################################################
##
## Calculate MIC scores
## Adopted from Mead et al 2019 (Mol Ecol) - thank you!!!
##
###########################################################################################
###########################################################################################
# Load libraries
library(minerva)
library(foreach)
library(doParallel)

# Note that 24 is a lot of nodes! This is meant to be run on a high memory cluster.
cl <- makeCluster(24)
registerDoParallel(cl)

###########################################################################################


run_metab_data_MIC <-
  function(test_ = F,
           spp = 0,
           x_ = "nitrogen") {
    # This function cycles through metabolites and calculates the MIC against a variable 
    # of interest.
    # test_: bool whether this should be a test run, with more reasonable criteria to do a
    # debugging run on a laptop
    # spp: 0 is bogr, 1 is spco
    # x_: the x variable of interest across which the metabolite change will be calculated
    # 
    # output: writes a .csv containing pvalues and MIC scores for each metabolite for a 
    # given x_
    
    # Read in metabolite data
    dat <-
      read.csv("data/clean_SpecAbund_metabolites.csv")
    # Transpose data
    metabs <- as.data.frame(t(dat[, -c(1)]))
    
    # Read in phys / nitrogen data
    phys_dat <-
      read.csv("data/phys_metabs_for_modules.csv")
    
    # Set bootstraps and number of metabolites. If not testing, the loop should do all 
    # of them
    if (!(test_)) {
      # Full
      n_metabs_ <- ncol(metabs)
      nBootstraps = 5000
    } else {
      # Test
      n_metabs_ <- 3
      nBootstraps = 5
    }
    
    # Parallelize the tests across metabolites. Results will be row bound.
    # Full version should be 1:(number of metabolites)
    spp_frame <-
      foreach(i = 1:n_metabs_, .combine = rbind) %dopar% {
        library(foreach)
        library(doParallel)
        library(minerva)
        
        # Function to test for difference in null and actual MIC
        pcalc <-
          function(nullDist,
                   observedStatistic,
                   sided = "two") {
            switch (sided, two = {
              # two-sided test: actual value is different than null
              # (either bigger or smaller)
              # calculate for most likely scenario, then multiply by 2
              # to account for using both tails
              pValue <-
                (min(
                  sum(nullDist > observedStatistic),
                  sum(nullDist < observedStatistic)
                ) / length(nullDist)) * 2
            },
            one = {
              # one-sided test: actual is bigger or smaller than null
              # get null observations that are less than actual
              # return pvalue for probability that observed is
              # greater than null
              # (would have to change for prob observed is less than null)
              pValue <-
                sum(nullDist > observedStatistic) / length(nullDist)
            },
            
            {
              stop("Unknown sidedness")
            })
            return(pValue)
          }
        
        # Specific metabolite ID
        metab_id <- names(metabs[i])
        
        # Y variable is the metabolite in question
        yvar <- metabs[phys_dat$spp == spp, i]
        
        # X var depends on function argument
        if (x_ == "nitrogen") {
          xvar <- phys_dat$nitrogen[phys_dat$spp == spp]
        } else if (x_ == "photo") {
          xvar <- phys_dat$photo[phys_dat$spp == spp]
        } else if (x_ == "cond") {
          xvar <- phys_dat$cond[phys_dat$spp == spp]
        } else if (x_ == "Ci") {
          xvar <- phys_dat$Ci[phys_dat$spp == spp]
        } else if (x_ == "iWUE") {
          xvar <- phys_dat$iWUE[phys_dat$spp == spp]
        }
        
        # Make sure that NAs are excluded, since they will interfere with the 
        # test
        yvar_no_na <- yvar[!(is.na(xvar))]
        xvar_no_na <- xvar[!(is.na(xvar))]
        
        # For each metabolite, you want to do null MIC draws to get a null
        # distribution. This is what will be tested against the actual MIC
        MIC.null <- vector()
        for (b in 1:nBootstraps) {
          nsamples <- length(xvar_no_na)
          resampledX <-
            xvar_no_na[sample(1:nsamples, replace = F)]
          mic <- mine(resampledX, yvar_no_na, n.cores = 2)
          mic$MIC
          # Append
          MIC.null[b] <- mic$MIC
        }
        
        # Calculate actual MIC
        MIC.actual <- mine(xvar_no_na, yvar_no_na)$MIC
        
        # Test for a difference between Null and Actual
        p.value = pcalc(
          nullDist = MIC.null,
          observedStatistic = MIC.actual,
          sided = "one"
        )
        
        # Assemble results in a dataframe
        out <-
          data.frame(
            metab_id = rep(metab_id, length(xvar_no_na)),
            metab = yvar_no_na,
            xvar_no_na,
            spp = spp,
            pval = rep(p.value, length(xvar_no_na)),
            mic = rep(MIC.actual, length(xvar_no_na))
          )
        
        # Have to explicitly name this column to the variable
        colnames(out)[3] <- x_
        
        # Output the dataframe so that results can be appended onto one another
        out
      }
    
    # Write results in a .csv that depends on the x_ in question
    if (x_ == "nitrogen") {
      write.csv(spp_frame,
                file = paste(spp, "_metabolomic_MIC_output_nitrogen.csv", sep = ""))
    } else if (x_ == "photo") {
      write.csv(spp_frame,
                file = paste(spp, "_metabolomic_MIC_output_photo.csv", sep = ""))
    } else if (x_ == "cond") {
      write.csv(spp_frame,
                file = paste(spp, "_metabolomic_MIC_output_cond.csv", sep = ""))
    } else if (x_ == "Ci") {
      write.csv(spp_frame,
                file = paste(spp, "_metabolomic_MIC_output_Ci.csv", sep = ""))
    } else if (x_ == "iWUE") {
      write.csv(spp_frame,
                file = paste(spp, "_metabolomic_MIC_output_iWUE.csv", sep = ""))
    }
  }
