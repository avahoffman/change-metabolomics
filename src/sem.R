###########################################################################################
##
## Path analysis / CFA / Structural equation modelling
##
###########################################################################################
###########################################################################################
# Load libraries

## working model
setwd(wd)
library(dplyr)
library(lavaan) # Requires gfortran on machine
library(semPlot)
# OpenMx dependency may need to be installed off of  old binary off CRAN website, for some
# reason new edition doesn't compile on MacOS
# https://cran.r-project.org/bin/macosx/el-capitan/contrib/3.4/OpenMx_2.7.10.tgz

# may not want to use chisq as a measure of model fit..
# https://www.researchgate.net/post/Is_it_necessary_that_in_model_fit_my_Chi-square_valuep-Value_must_be_non-significant_in_structure_equation_modeling_AMOS
# http://davidakenny.net/cm/fit.htm


###########################################################################################


get_sem_data <-
  function() {
    # Load data
    corr_data <-
      read.csv(file = "sem/all_data_for_corr.csv", header = T) %>%
      mutate(spp, spp = as.character(spp)) %>%
      mutate(spp = replace(spp, spp == "Bogr", 1)) %>%
      mutate(spp = replace(spp, spp == "Spco", 2))
    
    # Scale data
    corr_data[, 7:20] <- scale(corr_data[, 7:20])
    
    # Replace NAs with means
    corr_data <-
      corr_data %>%
      mutate(photo = replace(photo, is.na(photo), mean(corr_data$photo, na.rm =
                                                         T))) %>%
      mutate(cond = replace(cond, is.na(cond), mean(corr_data$cond, na.rm =
                                                      T))) %>%
      mutate(Ci = replace(Ci, is.na(Ci), mean(corr_data$Ci, na.rm = T))) %>%
      mutate(iWUE = replace(iWUE, is.na(iWUE), mean(corr_data$iWUE, na.rm =
                                                      T))) %>%
      mutate(PC1 = replace(PC1, is.na(PC1), mean(corr_data$PC1, na.rm =
                                                   T))) %>%
      mutate(PC2 = replace(PC2, is.na(PC2), mean(corr_data$PC2, na.rm =
                                                   T))) %>%
      mutate(M10 = replace(M10, is.na(M10), mean(corr_data$M10, na.rm =
                                                   T))) %>%
      mutate(M16 = replace(M16, is.na(M16), mean(corr_data$M16, na.rm =
                                                   T)))
    
    # View correlations
    print(cor(corr_data[,c( 7:18,21,22)]))
    
    return(corr_data)
  }


fit_sem_model <-
  function(corr_data, bogr = T) {
    ## Fit model for B. gracilis coveronly
    if (bogr) {
      corr_data_BOGR <-
        corr_data %>%
        filter(spp == 1)
      
      # SEM model
      model <- '
        # measurement model
        traits =~ BOGR_cover
        metabolome =~ traits

        # regressions

        ELEL_cover ~ nitrogen
        BOGR_cover ~ ELEL_cover + nitrogen
         M16 + M10 + photo + iWUE ~ nitrogen
        metabolome ~ M16 + M10
        traits ~ photo + iWUE

        ##covariance
        M16~~M10
        photo~~iWUE
      '
      
      # Fit model
      fit <- sem(model, data = corr_data_BOGR)
      fitMeasures(fit, c("nfi", "ifi"))
      
      # Predict latent variables
      lavPredict(fit)
      
      # All fit parameters
      parameterestimates(fit)
      sink("sem/sem_bogr_cover_summary.txt")
      summary(
        fit,
        fit.measures = TRUE,
        standardized = T,
        rsquare = T
      )
      fitMeasures(fit)
      sink()
      
      return(fit)
    } else if (!(bogr)) {
      corr_data_SPCO <-
        corr_data %>%
        filter(spp == 2)
      
      # SEM model
      model <- '
        # measurement model
        traits =~ SPCO_cover
        metabolome =~ traits

        # regressions

        ELEL_cover ~ nitrogen
        SPCO_cover ~ nitrogen + ELEL_cover
        M16 + M10 + photo + iWUE ~ nitrogen
        metabolome ~ M16 + M10
        traits ~ photo + iWUE

        ##covariance
        M16~~M10
        photo~~iWUE
      '
      
      # Fit model
      fit <- sem(model, data = corr_data_SPCO)
      fitMeasures(fit, c("nfi", "ifi"))
      
      # Predict latent variables
      lavPredict(fit)
      
      # All fit parameters
      parameterestimates(fit)
      sink("sem/sem_spco_cover_summary.txt")
      summary(
        fit,
        fit.measures = TRUE,
        standardized = T,
        rsquare = T
      )
      fitMeasures(fit)
      sink()
    }
  }


setup_plot_params <-
  function(fit) {
    ## custom labels
    # labs <- unlist(fit@pta$vnames$ov)
    # labs
    # labs[1] <- "BOGR\ncover"
    # labs[6] <- "ELEL\ncover"
    # labs[5] <- "M10"
    # labs[4] <- "M16"
    # labs[3] <- "photo"
    # labs[2] <- "CO2"
    # labs[7] <- "N"
    
    # Add custom  colors, will cycle through and color purple if neg and significant,
    # green if pos and significant, grey if NS
    params <-
      as.data.frame(parameterEstimates(fit, standardized = T))
    useparams <- data.frame()
    params$col <- rep(sem_neg_color, nrow(params))
    params$pvalue <- round(params$pvalue, 3)
    params$std.all <- round(params$std.all, 2)
    
    for (i in 1:nrow(params)) {
      if (!(params$lhs[i] == params$rhs[i])) {
        useparams <- rbind(useparams, params[i, ])
      }
    }
    for (i in 1:nrow(useparams)) {
      if (useparams$est[i] > 0) {
        useparams$col[i] <- sem_pos_color
      }
      if (!(is.na(useparams$pvalue[i]))) {
        if (useparams$pvalue[i] > 0.05) {
          useparams$col[i] <- "light grey"
          useparams$std.all[i] <- ""
        }
      }
      
    }
    useparams$col[nrow(useparams)] <- "light grey"
    useparams$std.all[nrow(useparams)] <- ""
    
    return(useparams)
  }


plot_sem <-
  function() {
    fit <- fit_sem_model(get_sem_data())
    useparams <- setup_plot_params(fit)
    
    labs <- c(
      "BOGR\ncover",
      "iWUE",
      "photo",
      "M16",
      "M10",
      "ELEL\ncover",
      "N"
    )
    
    ## weighted diagram
    pdf(file = "sem/estimates_diagram_BOGR_cover.pdf",
        width = 4,
        height = 2)
    graph <- semPaths(
      fit,
      layout = "tree2",
      nCharNodes = 0,
      what = "std",
      nodeLabels = c(labs, "traits", "metabolome"),
      edgeLabels = useparams$std.all,
      fixedStyle = c(1),
      ## all solid lines, dashed line indicates fixed parameter estimates
      residuals = F,
      edge.color = useparams$col,
      weighted = F
      #curve = 1.2
      #edge.label.position = 0.35
    )
    dev.off()
    
    ## Unweightd diagram
    pdf(file = "sem/unweighted_diagram_BOGR_cover.pdf",
        width = 4,
        height = 2)
    graph <- semPaths(
      fit,
      layout = "tree2",
      nCharNodes = 0,
      what = "path",
      nodeLabels = c(labs, "traits", "metabolome"),
      fixedStyle = c(1),
      ## all solid lines, dashed line indicates fixed parameter estimates
      residuals = F,
      weighted = F
    )
    dev.off()
    fitMeasures(fit, c("nnfi", "cfi"))
    
  }