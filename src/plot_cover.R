###########################################################################################
##
## Plot nitrogen addition vs. cover
##
###########################################################################################
# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)


###########################################################################################


plot_nitrogen_and_cover <-
  function(filename = NA) {
    community_data <-
      na.omit(read.csv("data/SpecAbund_community_design.csv")) %>%
      select(c(nitrogen, BOGR_cover, SPCO_cover, `ELEL_.cover`)) %>%
      gather(spp, cover,-nitrogen)
    
    gg <- 
      ggplot(data = community_data,
           aes(
             x = nitrogen,
             y = cover,
             shape = as.factor(spp),
             col = as.factor(spp)
           )) +
      theme_sigmaplot() +
      geom_smooth(
        data = community_data %>%
          filter(spp == "BOGR_cover"),
        method = lm,
        se = F,
        color = bogr_color
      ) +
      geom_smooth(
        data = community_data %>%
          filter(spp == "ELEL_.cover"),
        method = lm,
        se = F,
        color = elel_color
      ) +
      geom_jitter(size = 3,
                  width = 0.2,
                  height = 0) +
      legend_custom() +
      scale_shape_manual(values = c(bogr_shape,
                                    elel_shape,
                                    spco_shape),
                         guide = F) +
      scale_color_manual(
        values = c(bogr_color,
                   elel_color,
                   spco_color),
        labels = c(B.gracilis,
                   E.elymoides,
                   S.coccinea)
      ) +
      xlab(expression(paste("N addition (g ", m ^ {
        -2
      }, ')'))) +
      ylab("% cover") +
      guides(colour = guide_legend(override.aes = list(
        shape = c(bogr_shape,
                  elel_shape,
                  spco_shape)
      )))
    
    gg
    
    if (!(is.na(filename))) {
      ggsave(file = filename,
             height = 4.5,
             width = 4)
    }
    
    return(gg)
    
  }
