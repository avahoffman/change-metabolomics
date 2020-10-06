library(dplyr)
library(ggplot2)
library(reshape2)


dat <-
  read.csv("temp/1_metabolomic_MIC_output_nitrogen.csv")

df_final <- 
  dat %>%
  filter(mic > 0.6)

df_grouped <- 
  df_final %>%
  group_by(metab_id, nit, spp, pval, mic) %>%
  summarize(mean_metab = mean(metab, na.rm = T)) %>%
  filter(nit == 0 | nit == 30)

# Wide format to do subtraction
df_grouped_cast <- 
  dcast(df_grouped, metab_id + spp + pval + mic ~ nit,
      value.var = "mean_metab") %>%
  mutate(diff = ifelse(`30`-`0` > 0, "increase", "decrease")) %>%
  select(-c(spp, pval, mic, `0`, `30`))

df_final <- 
  df_final %>%
  left_join(df_grouped_cast, by = "metab_id")

p1 <- 
  ggplot() +
  stat_smooth(data = df_final[(df_final$diff == "increase"),],
              aes(nit,
                  metab,
                  group = metab_id),
              span = 1,
              color = "grey",
              se = F) 

p2 <- 
  ggplot() +
  stat_smooth(data = df_final[(df_final$diff == "decrease"),],
              aes(nit,
                  metab,
                  group = metab_id),
              span = 1,
              color = "grey",
              se = F) 


p1 + p2 + p2 + patchwork::plot_layout(ncol = 2)
