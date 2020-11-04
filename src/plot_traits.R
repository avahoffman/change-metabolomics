###########################################################################################
##
## Plot Physiological traits
##
###########################################################################################
###########################################################################################
# Load libraries
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)

###########################################################################################

phys_plot <-
  function(phys_data,
           trait,
           ylab,
           filename = NA,
           fit_line = F) {
    gg <-
      ggplot(data = phys_data,
             aes(
               x = nitrogen,
               y = get(trait),
               col = as.factor(spp),
               shape = as.factor(spp)
             )) +
      theme_cowplot() +
      ylab(ylab) +
      geom_jitter(size = 3,
                  width = 0.2,
                  height = 0) +
      scale_shape_manual(values = c(bogr_shape,
                                    spco_shape),
                         guide = F) +
      scale_color_manual(
        values = c(bogr_color,
                   spco_color),
        labels = c(B.gracilis,
                   S.coccinea)
      ) +
      xlab(expression(paste("N addition (g ", m ^ {
        -2
      }, ')'))) +
      guides(colour = guide_legend(override.aes = list(shape = c(
        bogr_shape,
        spco_shape
      )))) +
      legend_custom() +
      theme(legend.text.align = 0)
    
    if (fit_line) {
      gg <-
        gg +
        geom_smooth(
          data = subset(phys_data,
                        phys_data$spp == "Bogr"),
          method = lm,
          se = F,
          color = bogr_color
        ) +
        geom_smooth(
          data = subset(phys_data,
                        phys_data$spp == "Spco"),
          method = lm,
          se = F,
          color = spco_color
        )
    }
    
    gg
    
    if (!(is.na(filename))) {
      ggsave(file = filename,
             height = 3,
             width = 4.5)
    }
    
    return(gg)
  }


g_legend <- function(a.gplot) {
  tmp <-
    ggplot_gtable(ggplot_build(a.gplot))
  leg <-
    which(sapply(tmp$grobs,
                 function(x)
                   x$name) == "guide-box")
  legend <-
    tmp$grobs[[leg]]
  
  return(legend)
}


arrange_phys_plots <- function() {
  p1 <-
    phys_plot(
      prep_trait_data(),
      trait = "photo",
      ylab = photo_ylab,
      fit_line = T
    )
  leg <- g_legend(p1)
  p1 <-
    p1 +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    )
  p2 <-
    phys_plot(
      prep_trait_data(),
      trait = "cond",
      ylab = cond_ylab
    ) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    )
  p3 <-
    phys_plot(
      prep_trait_data(),
      trait = "Ci",
      ylab = ci_ylab,
      fit_line = T
    ) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    )
  p4 <-
    phys_plot(
      prep_trait_data(),
      trait = "iWUE",
      ylab = iWUE_ylab,
      fit_line = T
    ) +
    theme(legend.position = "none")
  
  main <-
    plot_grid(
      p1,
      p2,
      p3,
      p4,
      nrow = 4,
      align = "v",
      labels = "auto",
      label_size = 18,
      axis = "l",
      rel_heights = c(1,1,1,1.2)
    )
  
  pdf("figures/physiological_traits.pdf",
      height = 12,
      width = 4)
  grid.arrange(main,
               leg,
               heights = c(100, 7),
               nrow = 2)
  dev.off()
  
}
