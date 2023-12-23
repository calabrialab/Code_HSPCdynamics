main_path = "~/Dropbox/shared/HSC dynamics/"
setwd(main_path)
save_path = paste0(main_path, "diversity/")
require(tidyverse)
require(ggplot2)

# get the color palette
color_pal = openxlsx::read.xlsx(paste0(main_path, "diversity/data/202112.Diversity.ByMarker.xlsx")) %>%
  dplyr::filter(PCRMethod == "SLiM" | ClinicalTrial == "WAS") %>% 
  dplyr::select(colorcode, CellType) %>% unique()
color_pal = color_pal$colorcode %>% setNames(color_pal$CellType)
unnested = readRDS(paste0(save_path, "datasets/PRED_202112.Diversity.ByMarker.Hindex.Rds"))

## Plot regression over time ####

## Hindex
unnested %>%
  dplyr::filter(final_fit) %>%

  ggplot() +
  # geom_smooth(data=subs.tab,
  #             aes(x=TimePoint, y=Hindex, fill=CellType, color=CellType),
  #             alpha=.3, level=.75, linetype="solid", linewidth=0.7) +
  
  geom_smooth(aes(x=Timepoint, y=pred_Hindex, fill=CellType, color=CellType),
              alpha=.6, level=.75, linetype="solid", linewidth=0.7) +
  
  scale_fill_manual(values=color_pal) +
  scale_color_manual(values=color_pal) +
  
  ylab("Hindex") +
  
  ggh4x::facet_nested(Tissue ~ factor(ClinicalStudy, levels=c("MLD","WAS","BTHAL")),
                      scales="free_x", space="free") +
  scale_x_continuous(breaks = seq(0, max(unnested$Timepoint, na.rm = T), 12)) +
  theme_bw() +
  theme(strip.text=element_text(face="bold", size=16)) +
  theme(legend.direction="horizontal", legend.position="bottom", legend.box="horizontal") +
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16),
        axis.title=element_text(size=16), plot.title=element_text(size=20)) +
  labs(x="Months after gene therapy",
       y="Shannon diversity index (H')", colour="Cell Marker", fill="Cell Marker")

ggsave(paste0(save_path, "plots/PRED_hindex.pdf"), height=8, width=12)
ggsave(paste0(save_path, "plots/PRED_hindex.png"), height=8, width=12, dpi=600)




## Hindex scaled
unnested %>%
  dplyr::filter(final_fit) %>%

  ggplot() +
  geom_smooth(aes(x=Timepoint, y=pred_h_zscaled, fill=CellType, color=CellType),
              alpha=.6, level=.75, linetype="solid", linewidth=0.7) +
  
  scale_fill_manual(values=color_pal) +
  scale_color_manual(values=color_pal) +
  
  ylab("Hindex") +
  
  ggh4x::facet_nested(Tissue ~ factor(ClinicalStudy, levels=c("MLD","WAS","BTHAL")),
                      scales="free_x", space="free") +
  scale_x_continuous(breaks = seq(0, max(unnested$Timepoint, na.rm = T), 12)) +
  theme_bw() +
  theme(strip.text=element_text(face="bold", size=16)) +
  theme(legend.direction="horizontal", legend.position="bottom", legend.box="horizontal") +
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16),
        axis.title=element_text(size=16), plot.title=element_text(size=20)) +
  labs(x="Months after gene therapy",
       y="Shannon diversity index (H')", colour="Cell Marker", fill="Cell Marker")


ggsave(paste0(save_path, "plots/PRED_hindex_scaled.pdf"), height=8, width=12)
ggsave(paste0(save_path, "plots/PRED_hindex_scaled.png"), height=8, width=12, dpi=600)







## Slice diversity plots ####
library(ggpubr); library(rstatix)

study_list = c("MLD", "WAS", "BTHAL")
celltype_list = c("CD34", "Erythroid", "Myeloid", "B", "T")
faceting_order = c("MLD", "WAS", "BTHAL")
markerlist = c("CD13", "CD14", "CD15", "CD19", "CD3", "CD34", "GLY", "CD36", "GLYA")
trials_colors = c("darkblue", "forestgreen", "firebrick") %>% setNames(study_list)


patient_h_scaled_summary_filtered = read.csv(
  paste0(save_path, "datasets/pred.diversity.zscore.patient_h_scaled_summary_filtered.csv"),
  row.names=1) %>% 
  dplyr::mutate(final_fit=fitname=="fit2")


my_comparisons = list(
  c("MLD", "WAS"),
  c("MLD", "BTHAL"),
  c("WAS", "BTHAL")
)
perc_spacer_stats = -3
perc_starting_max_stats = 18

p_h_scaled_summary_stable_min3_stats = patient_h_scaled_summary_filtered %>% 
  
  dplyr::filter(final_fit) %>% 
  dplyr::mutate(CellType=factor(CellType, levels=celltype_list),
                ClinicalStudy=factor(ClinicalStudy, levels=study_list)) %>% 
  
  ggplot(aes(x=ClinicalStudy, y=mean, group=ClinicalStudy), na.rm=T, se=TRUE) +
  geom_boxplot(colour="lightblue", fill="lightblue", alpha=.3, outlier.size=0) +
  geom_jitter(aes(color=ClinicalStudy), alpha=0.6, size=4) +
  facet_grid(Tissue ~ CellType) +
  scale_color_manual(values=trials_colors) +
  scale_y_continuous(limits=c(-2, 2.3)) + 
  theme_bw() +
  stat_compare_means(comparisons=my_comparisons, 
                     label.y=seq(from=perc_starting_max_stats, 
                                 to=(perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), 
                                 by=perc_spacer_stats)*10^-1, 
                     p.adjust.method="fdr", label="p.signif") +
  
  theme(strip.text=element_text(face="bold", size=16)) +
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16), axis.title=element_text(size=16), plot.title=element_text(size=20)) +
  labs(title=paste0("Clonal Diversity"), 
       subtitle="Filters: samples with homogeneous PCR technology within the same trial, PB and BM separated, markers separated. Months under analysis: 24-60.\nNormalization: Z-score by patient, tissue, time point.",
       x="Studies", 
       y="Shannon diversity index (Z-score)", 
       colour="CellMarker", fill="CellMarker")

ggsave(plot=p_h_scaled_summary_stable_min3_stats,
       filename=paste0(save_path, "plots/PRED.diversity.Zscore.stats_across_trials.pdf"), height=8, width=12)
ggsave(plot=p_h_scaled_summary_stable_min3_stats,
       filename=paste0(save_path, "plots/PRED.diversity.Zscore.stats_across_trials.png"), height=8, width=12, dpi=600)




## within trials
my_comparisons = list(
  c("Erythroid", "Myeloid"),
  c("Erythroid", "B"),
  c("Erythroid", "T"),
  c("Myeloid", "B"),
  c("Myeloid", "T"),
  c("B", "T")
)

perc_spacer_stats = -3
perc_starting_max_stats = 26
p_h_scaled_summary_stable_min3_stats_lineage = patient_h_scaled_summary_filtered %>% 
  
  dplyr::filter(final_fit) %>% 
  dplyr::mutate(CellType=factor(CellType, levels=celltype_list)) %>% 
  
  ggplot(aes(x=CellType, y=mean, group=CellType), na.rm=T, se=TRUE) +
  
  geom_boxplot(colour="lightblue", fill="lightblue", alpha=.3, outlier.size=0) +
  geom_jitter(aes(color=CellType), alpha=0.6, size=4) +
  facet_grid(Tissue ~ factor(ClinicalStudy, levels=faceting_order)) +
  theme_bw() +
  scale_y_continuous(limits=c(-2, 3)) +
  stat_compare_means(comparisons=my_comparisons, 
                     label.y=seq(from=perc_starting_max_stats, 
                                 to=(perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), 
                                 by=perc_spacer_stats)*10^-1, 
                     p.adjust.method="fdr", label="p.signif") +
  theme(strip.text=element_text(face="bold", size=16)) +
  theme(axis.text.x=element_text(size=16, angle=30), axis.text.y=element_text(size=16), axis.title=element_text(size=16), plot.title=element_text(size=20)) +
  labs(title=paste0("Clonal Diversity"), 
       subtitle="Filters: samples with homogeneous PCR technology within the same trial, PB and BM separated, markers separated. Months under analysis: 24-60.\nNormalization: Z-score by patient, tissue, time point.",
       x="Studies", 
       y="Shannon diversity index (Z-score)", 
       colour="CellMarker", fill="CellMarker")

ggsave(plot=p_h_scaled_summary_stable_min3_stats_lineage,
       filename=paste0(save_path, "plots/PRED.diversity.Zscore.stats_within_trials.pdf"), height=8, width=12)
ggsave(plot=p_h_scaled_summary_stable_min3_stats_lineage,
       filename=paste0(save_path, "plots/PRED.diversity.Zscore.stats_within_trials.png"), height=8, width=12, dpi=600)





