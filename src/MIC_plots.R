###########################################################################################
##
## Plot Physiological traits
##
###########################################################################################
###########################################################################################
# Load libraries
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)
library(patchwork)

###########################################################################################


make_mic_plots_nit <-
  function(dat_bogr, dat_spco, span_ = 1) {
    df_final <-
      
      # Note on p value cutoff:
      # if p cutoff is < 0.001, that means for every 2570 metabs,
      # we expect 2.6 to have a significant relationship by chance
      #
      # 2570 total metabolites in our dataset
      rbind(read.csv(dat_bogr), read.csv(dat_spco)) %>%
      filter(pval < 0.001)
    
    df_unique <- 
      df_final %>%
      select(metab_id, spp) %>%
      distinct() %>%
      arrange(metab_id)
    
    df_grouped <-
      df_final %>%
      group_by(metab_id, nit, spp, pval, mic) %>%
      summarize(mean_metab = mean(metab, na.rm = T)) %>%
      filter(nit == 0 | nit == 30)
    
    # Wide format to do subtraction
    df_grouped_cast <-
      dcast(df_grouped, metab_id + spp + pval + mic ~ nit,
            value.var = "mean_metab") %>%
      mutate(diff = ifelse(`30` - `0` > 0, "increase", "decrease")) %>%
      select(-c(spp, pval, mic, `0`, `30`))
    
    df_final <-
      df_final %>%
      left_join(df_grouped_cast, by = "metab_id")
    
    plotList <- list()
    
    # Leaving this relic here for now. Very undecided how I want the plots to look
    #for (spp_ in c(0, 1)) {
      for (diff_ in c("increase", "decrease")) {
        p_ <-
          ggplot() +
          stat_smooth(
            #data = df_final[(df_final$diff == diff_ & df_final$spp == 0),],
            data = df_final[(df_final$diff == diff_),],
            aes(nit,
                metab,
                color = as.factor(spp),
                group = metab_id),
            span = span_, # lower numbers are less smoothed
            #color = ifelse(spp_ == 0, bogr_color, spco_color),
            se = F
          ) +
          scale_color_manual(values = c(bogr_color, spco_color)) +
          # stat_smooth(
          #   data = df_final[(df_final$diff == diff_ & df_final$spp == 1),],
          #   aes(nit,
          #       metab,
          #       group = metab_id),
          #   span = span_, # lower numbers are less smoothed
          #   #color = ifelse(spp_ == 0, bogr_color, spco_color),
          #   color = spco_color,
          #   se = F
          # ) +
          theme_cowplot() +
          xlab(expression(paste("N addition (g ", m ^ {
            -2
          }, ')'))) +
          ylab("Metabolite abundance")
          
        ind_ <- 
          #paste(diff_, spp_)
          paste(diff_)
        
        plotList[[ind_]] = p_
        
      }
   # }
    
    return(plotList)
  }


make_mic_plots_phys <-
  function(dat_bogr, dat_spco) {
    df_final <-
      
      # Note on p value cutoff:
      # if p cutoff is < 0.001, that means for every 2570 metabs,
      # we expect 2.6 to have a significant relationship by chance
      #
      # 2570 total metabolites in our dataset
      rbind(read.csv(dat_bogr), read.csv(dat_spco)) %>%
      filter(pval < 0.005)
    
    response_var_part <-
      gsub("temp/0_metabolomic_MIC_output_", "", dat_bogr)
    response_var <-
      gsub(".csv", "", response_var_part)
    
    plotList <- list()
    
    for (spp_ in c(0, 1)) {
      p_ <-
        ggplot() +
        stat_smooth(
          data = df_final[(df_final$spp == spp_),],
          aes(get(response_var),
              metab,
              group = metab_id),
          span = 1,
          color = ifelse(spp_ == 0, bogr_color, spco_color),
          se = F
        ) +
        xlab(response_var)
      
      plotList[[spp_ + 1]] <- p_
      
    }
    
    return(plotList)
  }


gather_nit_plots <- 
  function(filename = NA, span_ = 1){
    nit_plots <-
      make_mic_plots_nit(dat_bogr = "temp/0_metabolomic_MIC_output_nitrogen.csv",
                         dat_spco = "temp/1_metabolomic_MIC_output_nitrogen.csv",
                         span_ = span_)
    
    main <-
      plot_grid(
        nit_plots[[1]],
        nit_plots[[2]],
        nrow = 1,
        align = "h",
        labels = c("b","c"),
        label_size = 18,
        axis = "l"
      )

    gg 
    
    if (!(is.na(filename))) {
      ggsave(file = filename,
             height = 3,
             width = 4.5)
    }
    
    return(gg)
  }


gather_phys_plots <- 
  function(filename = NA){
    cond_plots <-
      make_mic_plots_phys(dat_bogr = "temp/0_metabolomic_MIC_output_cond.csv",
                          dat_spco = "temp/1_metabolomic_MIC_output_cond.csv")
    photo_plots <-
      make_mic_plots_phys(dat_bogr = "temp/0_metabolomic_MIC_output_photo.csv",
                          dat_spco = "temp/1_metabolomic_MIC_output_photo.csv")
    Ci_plots <-
      make_mic_plots_phys(dat_bogr = "temp/0_metabolomic_MIC_output_Ci.csv",
                          dat_spco = "temp/1_metabolomic_MIC_output_Ci.csv")
    iWUE_plots <-
      make_mic_plots_phys(dat_bogr = "temp/0_metabolomic_MIC_output_iWUE.csv",
                          dat_spco = "temp/1_metabolomic_MIC_output_iWUE.csv")
    
    gg <- 
      photo_plots[[1]] +
      photo_plots[[2]] +
      cond_plots[[1]] +
      cond_plots[[2]] +
      Ci_plots[[1]] +
      Ci_plots[[2]] +
      iWUE_plots[[1]] +
      iWUE_plots[[2]] +
      patchwork::plot_layout(ncol = 2, nrow = 4)
    
    gg
    
    if (!(is.na(filename))) {
      ggsave(file = filename,
             height = 3,
             width = 4.5)
    }
    
    return(gg)
  }

gather_nit_plots(filename = "figures/MIC.pdf", span_ = 0.75)

plot_pca(run_pca()) + guides(colour = guide_legend(override.aes = list(linetype = c(1,1)))) +
  legend_custom()
