## Entropy
ProjectID <- "CTC"
dest_dir <- "analyses/202112/"
vector_size <- 6219
analysis_folder_date <- "202112"
color_schema <- c("navyblue", "orange")
study_list <- c("MLD", "WAS", "BTHAL")
celltype_list <- c("CD34", "Erythroid", "Myeloid", "B", "T")
faceting_order <- c("MLD", "WAS", "BTHAL")
trials_colors <- c("darkblue", "forestgreen", "firebrick")
id_cols <- c("chr", "integration_locus", "strand", "GeneName", "GeneStrand")
markerlist = c("CD13", "CD14", "CD15", "CD19", "CD3", "CD34", "GLY", "CD36", "GLYA")

all_results_h = read.csv("202112.Diversity.ByMarker.Hindex.csv")


all_results_h_slice <- all_results_h[which(all_results_h$CellMarker %in% markerlist 
                                           & all_results_h$Tissue %in% c("PB", "BM") 
                                           & !is.na(all_results_h$Hindex)
                                           # & all_results_h$nIS >= 10
                                           # & all_results_h$`ng DNA corrected_sum` >= 10
                                           & all_results_h$TimePoint > 0),]

all_results_h_slice_scaled <- all_results_h_slice %>% group_by(SubjectID, Tissue, TimePoint) %>% mutate(h_zscaled=scale(Hindex))
all_results_h_slice_scaled <- all_results_h_slice_scaled %>% group_by(SubjectID, Tissue, TimePoint) %>% mutate(h_normalized_zscaled=scale(Hindex_normalized))
all_results_h_slice_scaled_homo <- all_results_h_slice_scaled[which(all_results_h_slice_scaled$PCRMethod == "SLiM" | all_results_h_slice_scaled$Study == "WAS"),]

# write results in a file
write.xlsx(x = all_results_h_slice_scaled, file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.xlsx", sep = ""))
write.xlsx(x = all_results_h_slice_ngdna_scaled, file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.ngdna_notNA.xlsx", sep = ""))


### zscore stats
all_results_h_slice_scaled_homo_stable <- all_results_h_slice_scaled_homo[which(all_results_h_slice_scaled_homo$TimePoint >= 24 & all_results_h_slice_scaled_homo$TimePoint <= 62),]
all_results_h_slice_scaled_homo_stable <- all_results_h_slice_scaled_homo_stable %>% group_by(SubjectID, Tissue, TimePoint) %>% mutate(n_obs_stt_overtime = n())
all_results_h_slice_scaled_homo_stable_min3 <- all_results_h_slice_scaled_homo_stable %>% group_by(SubjectID, Tissue, TimePoint) %>% mutate(n_obs_stt_overtime = n()) %>% filter(n_obs_stt_overtime >= 3)


# summary stats
patient_h_scaled_summary <- all_results_h_slice_scaled_homo_stable_min3 %>%
  filter(n_obs_stt_overtime >= 3) %>%
  group_by(Study, SubjectID, Tissue, CellMarker, CellType) %>%
  # get_summary_stats(PercNIS_onOverallSharedCD34BM_zscaled, type = "mean_sd")
  get_summary_stats(h_zscaled, type = "full")

patient_h_scaled_summary_filtered <- patient_h_scaled_summary %>% filter(n >= 2)
patient_h_scaled_summary_filtered$CellLineage <- factor(patient_h_scaled_summary_filtered$CellType, levels = celltype_list)
patient_h_scaled_summary_filtered$Study <- factor(patient_h_scaled_summary_filtered$Study, levels = study_list)
# patient_h_scaled_summary_filtered$AgeGroup_Sorted <- factor(patient_h_scaled_summary_filtered$AgeGroup, levels = age_list)

write.xlsx(x = patient_h_scaled_summary_filtered, file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".", paste(study_list, collapse = "_"), ".diversity.zscore.24-60m.patient_h_scaled_summary_filtered.xlsx", sep = ""))





# do plots
my_comparisons <- list(
  c("MLD", "WAS"),
  c("MLD", "BTHAL"),
  c("WAS", "BTHAL")
)
perc_spacer_stats <- -3
perc_starting_max_stats <- 18
p_h_scaled_summary_stable_min3_stats <- ggplot(data = patient_h_scaled_summary_filtered,
                                                # aes(x = Study, y = mean, colour = CellMarker, fill = CellMarker, group = Study, shape = Tissue ), na.rm = T, se = TRUE) +
                                                aes(x = Study, y = mean, group = Study), na.rm = T, se = TRUE) +
  # scale_color_manual(values = trials_colors) +
  # scale_fill_manual(values = trials_colors) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  # geom_violin(colour = "darkblue", fill = "lightblue", alpha = 0.6) +
  # geom_line(size=3, alpha = .7) +
  # facet_wrap(. ~ CellLineage, ncol = 4) +
  facet_grid(Tissue ~ CellLineage) +
  scale_y_continuous(limits = c(-2, 2.3)) + 
  theme_bw() +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(30, 35, 40), p.adjust.method = "fdr", label = "p.signif") +
  stat_compare_means(comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Clonal Diversity"), 
       subtitle = "Filters: samples with homogeneous PCR technology within the same trial, PB and BM separated, markers separated. Months under analysis: 24-60.\nNormalization: Z-score by patient, tissue, time point.",
       x = "Studies", 
       y = "Shannon diversity index (Z-score)", 
       colour = "CellMarker", fill = "CellMarker")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.HomoPCR.Zscore.statistics_across_trials.pdf", sep = ""), height=8, width=12)
plot(p_h_scaled_summary_stable_min3_stats)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.HomoPCR.Zscore.statistics_across_trials.png", sep = ""), height=8, width=12, units = "in", res = 300)
plot(p_h_scaled_summary_stable_min3_stats)
dev.off()

p_h_scaled_summary_stable_min3_stats_newcolors <- ggplot(data = patient_h_scaled_summary_filtered,
                                               aes(x = Study, y = mean, colour = Study, fill = Study, group = Study), na.rm = T, se = TRUE) +
                                               #aes(x = Study, y = mean, group = Study), na.rm = T, se = TRUE) +
  scale_color_manual(values = trials_colors) +
  scale_fill_manual(values = trials_colors) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  # geom_violin(colour = "darkblue", fill = "lightblue", alpha = 0.6) +
  # geom_line(size=3, alpha = .7) +
  # facet_wrap(. ~ CellLineage, ncol = 4) +
  facet_grid(Tissue ~ CellLineage) +
  scale_y_continuous(limits = c(-2, 2.3)) + 
  theme_bw() +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(30, 35, 40), p.adjust.method = "fdr", label = "p.signif") +
  stat_compare_means(comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Clonal Diversity"), 
       subtitle = "Filters: samples with homogeneous PCR technology within the same trial, PB and BM separated, markers separated. Months under analysis: 24-60.\nNormalization: Z-score by patient, tissue, time point.",
       x = "Studies", 
       y = "Shannon diversity index (Z-score)", 
       colour = "Study", fill = "Study")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.HomoPCR.Zscore.statistics_across_trials.colored.pdf", sep = ""), height=8, width=12)
plot(p_h_scaled_summary_stable_min3_stats_newcolors)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.HomoPCR.Zscore.statistics_across_trials.colored.png", sep = ""), height=8, width=12, units = "in", res = 300)
plot(p_h_scaled_summary_stable_min3_stats_newcolors)
dev.off()


## find among lineages
my_comparisons <- list(
  c("Erythroid", "Myeloid"),
  c("Erythroid", "B"),
  c("Erythroid", "T"),
  c("Myeloid", "B"),
  c("Myeloid", "T"),
  c("B", "T")
)

perc_spacer_stats <- -3
perc_starting_max_stats <- 26
p_h_scaled_summary_stable_min3_stats_lineage <- ggplot(data = patient_h_scaled_summary_filtered,
                                                        # aes(x = CellLineage, y = mean, colour = CellMarker, fill = CellMarker, group = CellLineage, shape = Tissue), na.rm = T, se = TRUE) +
                                                        aes(x = CellLineage, y = mean, group = CellLineage), na.rm = T, se = TRUE) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  # facet_wrap(. ~ Study, nrow = 1) +
  facet_grid(Tissue ~ Study) +
  theme_bw() +
  scale_y_continuous(limits = c(-2, 3)) +
  stat_compare_means(comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16, angle = 30), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Clonal Diversity"), 
       subtitle = "Filters: samples with homogeneous PCR technology within the same trial, PB and BM separated, markers separated. Months under analysis: 24-60.\nNormalization: Z-score by patient, tissue, time point.",
       x = "Studies", 
       y = "Shannon diversity index (Z-score)", 
       colour = "CellMarker", fill = "CellMarker")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.HomoPCR.Zscore.statistics_within_trials.pdf", sep = ""), height=8, width=12)
plot(p_h_scaled_summary_stable_min3_stats_lineage)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.HomoPCR.Zscore.statistics_within_trials.png", sep = ""), height=8, width=12, units = "in", res = 300)
plot(p_h_scaled_summary_stable_min3_stats_lineage)
dev.off()

p_h_scaled_summary_stable_min3_stats_lineage_newcolors <- ggplot(data = patient_h_scaled_summary_filtered,
                                                       aes(x = CellLineage, y = mean, colour = CellType, fill = CellType, group = CellLineage), na.rm = T, se = TRUE) +
                                                       # aes(x = CellLineage, y = mean, group = CellLineage), na.rm = T, se = TRUE) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # facet_wrap(. ~ Study, nrow = 1) +
  facet_grid(Tissue ~ Study) +
  theme_bw() +
  scale_y_continuous(limits = c(-2, 3)) +
  stat_compare_means(comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16, angle = 30), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Clonal Diversity"), 
       subtitle = "Filters: samples with homogeneous PCR technology within the same trial, PB and BM separated, markers separated. Months under analysis: 24-60.\nNormalization: Z-score by patient, tissue, time point.",
       x = "Studies", 
       y = "Shannon diversity index (Z-score)", 
       colour = "Lineage", fill = "Lineage")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.HomoPCR.Zscore.statistics_within_trials.colors.pdf", sep = ""), height=8, width=12)
plot(p_h_scaled_summary_stable_min3_stats_lineage_newcolors)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.HomoPCR.Zscore.statistics_within_trials.colors.png", sep = ""), height=8, width=12, units = "in", res = 300)
plot(p_h_scaled_summary_stable_min3_stats_lineage_newcolors)
dev.off()

