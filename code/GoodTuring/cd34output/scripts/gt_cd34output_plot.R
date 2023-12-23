main_path = "~/Dropbox/shared/HSC dynamics/"
ratio_df = read.csv(paste0(main_path, "cd34output/datasets/gt_abundance_ratio.csv")) %>% 
  tidyr::pivot_longer(cols=c("gt_ratio"),
                      names_to="method", values_to="ratio")

pred_ratio_df = read.csv(paste0(main_path, "regression_datasets/PRED_gtratio_cd34shared.csv")) %>% 
  dplyr::rename(ratio=pred, TimepointMonths=Timepoint) %>% 
  dplyr::filter(final_fit)

color_pal = ratio_df$colorcode %>% unique() %>%
  setNames(ratio_df$CellType %>% unique())

plot_byage = ratio_df %>% 
  dplyr::filter(ClinicalStudy!="XSCID") %>%
  ggplot() +
  stat_smooth(aes(x=TimepointMonths, y=ratio*100,
                  color=CellType, fill=CellType),
              formula=y~log(x),
              method="loess",
              se=T, level=.7) +
  
  scale_color_manual(values=color_pal) +
  scale_fill_manual(values=color_pal) +
  
  scale_x_continuous(breaks = seq(0, max(ratio_df$TimepointMonths, na.rm = T), 12)) +
  
  ggh4x::facet_nested(factor(AgeGroup, levels=c("0-2", "2-15", "30+")) ~
                        factor(ClinicalStudy, levels=c("MLD","WAS","BTHAL")),
                      scales="free_x", space="free_x") +
  
  theme_bw() +
  theme(strip.text=element_text(face="bold", size=16)) +
  theme(legend.direction="horizontal", legend.position="bottom", legend.box="horizontal") +
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16),
        axis.title=element_text(size=16), plot.title=element_text(size=20)) +
  labs(x="FollowUp months after gene therapy",
       y="Perc. of IS shared with CD34 BM",
       colour="CellType", fill="CellType")


plot_noage = ratio_df %>% 
  dplyr::filter(TimepointMonths <= 60) %>% 
  ggplot() +
  stat_smooth(aes(x=TimepointMonths, y=ratio*100,
                  color=CellType, fill=CellType),
              formula=y~log(x),
              method="loess",
              se=T, level=.7) +
  
  scale_color_manual(values=color_pal) +
  scale_fill_manual(values=color_pal) +
  
  scale_x_continuous(breaks = seq(0, max(ratio_df$TimepointMonths, na.rm = T), 12)) +
  
  ggh4x::facet_nested(. ~
                        factor(ClinicalStudy, levels=c("MLD","WAS","BTHAL","XSCID")),
                      scales="free", space="free_x") +
  
  theme_bw() +
  theme(strip.text=element_text(face="bold", size=16)) +
  theme(legend.direction="horizontal", legend.position="bottom", legend.box="horizontal") +
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16),
        axis.title=element_text(size=16), plot.title=element_text(size=20)) +
  labs(x="FollowUp months after gene therapy",
       y="Perc. of IS shared with CD34 BM",
       colour="CellType", fill="CellType")
  

