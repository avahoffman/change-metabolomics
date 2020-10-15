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

cl <- makeCluster(24)
registerDoParallel(cl)

###########################################################################################


run_metab_data_MIC <-
  function(test_ = F,
           spp = 0,
           x_ = "nitrogen") {
    dat <-
      read.csv("data/clean_SpecAbund_metabolites.csv")
    metabs <- as.data.frame(t(dat[, -c(1)]))
    
    phys_dat <-
      read.csv("data/phys_metabs_for_modules.csv")
    
    if (!(test_)) {
      n_metabs_ <- ncol(metabs)
      nBootstraps = 5000
    } else {
      n_metabs_ <- 3
      nBootstraps = 5
    }
    
    spp_frame <-
      foreach(i = 1:n_metabs_, .combine = rbind) %dopar% {
        #ncol(metabs)
        library(foreach)
        library(doParallel)
        library(minerva)
        
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
        
        metab_id <- names(metabs[i])
        print(metab_id)
        yvar <- metabs[phys_dat$spp == spp, i]
        if (x_ == "nitrogen") {
          xvar <- phys_dat$nitrogen[phys_dat$spp == spp]
        } else if (x_ == "photo") {
          xvar <- phys_dat$photo[phys_dat$spp == spp]
        }
        
        yvar_no_na <- yvar[!(is.na(xvar))]
        xvar_no_na <- xvar[!(is.na(xvar))]
        
        MIC.null <- vector()
        for (b in 1:nBootstraps) {
          nsamples <- length(xvar_no_na)
          resampledX <-
            xvar_no_na[sample(1:nsamples, replace = F)]
          mic <- mine(resampledX, yvar_no_na , n.cores = 2)
          mic$MIC
          MIC.null[b] <- mic$MIC
        }
        
        MIC.actual <- mine(xvar_no_na, yvar_no_na)$MIC
        p.value = pcalc(
          nullDist = MIC.null,
          observedStatistic = MIC.actual,
          sided = "one"
        )
        
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
        
        out
      }
    
    if (x_ == "nitrogen") {
      write.csv(spp_frame,
                file = paste(spp, "_metabolomic_MIC_output_nitrogen.csv", sep = ""))
    } else if (x_ == "photo") {
      write.csv(spp_frame,
                file = paste(spp, "_metabolomic_MIC_output_photo.csv", sep = ""))
    }
  }
