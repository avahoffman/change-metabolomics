library(dplyr)
library(ggplot2)
library(reshape2)


make_mic_plots <-
  function(dat) {
    df_final <-
      dat %>%
      filter(pval < 0.0003)
    
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
    
    p1 <-
      ggplot() +
      stat_smooth(
        data = df_final[(df_final$diff == "increase"), ],
        aes(nit,
            metab,
            group = metab_id),
        span = 1,
        color = "grey",
        se = F
      )
    
    p2 <-
      ggplot() +
      stat_smooth(
        data = df_final[(df_final$diff == "decrease"), ],
        aes(nit,
            metab,
            group = metab_id),
        span = 1,
        color = "grey",
        se = F
      )
    
    return(list(p1, p2))
  }

spco_plots <- make_mic_plots(read.csv("temp/1_metabolomic_MIC_output_nitrogen.csv"))
bogr_plots <- make_mic_plots(read.csv("temp/0_metabolomic_MIC_output_nitrogen.csv"))

spco_plots[[1]] + 
  spco_plots[[2]] + 
  bogr_plots[[1]] +
  bogr_plots[[2]] +
  patchwork::plot_layout(ncol = 2, nrow = 2)
