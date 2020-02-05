###########################################################################################
##
## Plot Physiological traits
##
###########################################################################################
# Load libraries
library(dplyr)
library(ggplot2)
library(gridExtra)


###########################################################################################


prep_trait_data <- function(){
  phys_data <- 
    read.csv("data/clean_SpecAbund_design.csv") %>%
    filter(!(is.na(photo))) %>%
    filter(!(X == "A2_Bogr")) %>% # all zeros, not helpful datapoint
    mutate(iWUE = photo / Trmmol) # generate new variable for water use efficiency
  return(phys_data)
}


phys_plot <-
  function(phys_data,
           trait,
           ylab,
           grid_label,
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
      theme_sigmaplot() +
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
      labs(colour = "Species") +
      annotate(
        "text",
        Inf,
        Inf,
        label = grid_label,
        hjust = 1,
        vjust = 1,
        size = 10
      ) +
      guides(colour = guide_legend(override.aes = list(shape = c(
        bogr_shape,
        spco_shape
      )))) +
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

phys_plot(prep_trait_data(), trait = "photo", ylab = photo_ylab, grid_label = "(a)", fit_line = T)
phys_plot(prep_trait_data(), trait = "cond", ylab = cond_ylab, grid_label = "(b)")
phys_plot(prep_trait_data(), trait = "Ci", ylab = ci_ylab, grid_label = "(c)")
phys_plot(prep_trait_data(), trait = "iWUE", ylab = iWUE_ylab, grid_label = "(d)")


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


arrange_phys_plots <- function(){
  # TODO: Needs some grid arrange cleanup
  p1 <- 
    phys_plot(prep_trait_data(), 
              trait = "photo", 
              ylab = photo_ylab, 
              grid_label = "(a)", 
              fit_line = T)
  leg <- g_legend(p1)
  p1 <- 
    p1 + 
    theme(legend.position="none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  p2 <- 
    phys_plot(prep_trait_data(), 
              trait = "cond", 
              ylab = cond_ylab, 
              grid_label = "(b)") +
    theme(legend.position="none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  p3 <- 
    phys_plot(prep_trait_data(), 
              trait = "Ci", 
              ylab = ci_ylab, 
              grid_label = "(c)") +
    theme(legend.position = "none")
  p4 <- 
    phys_plot(prep_trait_data(), 
              trait = "iWUE", 
              ylab = iWUE_ylab, 
              grid_label = "(d)") +
    theme(legend.position = "none")
  
  pdf("figures/physiological_traits.pdf",height=8,width=9)
  grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow = 2), 
               leg, widths = c(9,1), ncol = 2)
  dev.off()
  
}
