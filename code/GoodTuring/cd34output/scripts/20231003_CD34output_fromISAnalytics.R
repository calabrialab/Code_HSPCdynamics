###############################################################
## CD34 output -> ISAnalytics
###############################################################
library(openxlsx)

###############################################################
## functions
###############################################################
source("/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/TIGET/Workbench/isatk/script/R/isa_utils_functions.R")
id_cols <- c("chr", "integration_locus", "strand", "GeneName", "GeneStrand")

###############################################################
### GLOBAL SETTINGS AND FUNCTIONS
###############################################################
# ProjectID <- "MLD"
analysis_folder_date <- "202112"
results_folder_name <- "cd34output/"
this_run_foldername <- "" # "MyBTEry34/" "MyBTEry34_ISgt6/" "MyBTEry34_BMonly/" "MyBTEry34_PBonly/"
allpatients_shared34_stats <- NULL
patients_to_exclude_from_cd34output_analysis <- c("WASCUP03", "WASCUP05")
followup_limit_crosstrial <- 60
study_list <- c("MLD", "WAS", "BTHAL")
celltype_list <- c("CD34", "Erythroid", "Myeloid", "B", "T")
faceting_order <- c("MLD", "WAS", "BTHAL")
trials_colors <- c("darkblue", "forestgreen", "firebrick")
plot_footnote_details <- "All patients and all cell markers"
plot_filename_prefix <- paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.Rev2310.", sep = "")
plot_title <- paste0("CD34 BM outout over time")
data_filename <- paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.Rev2310.data.tsv", sep = "")


# metadata of cell markers and patients
blood_lineages_file <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial MLD/blood_lineages_update.xlsx"
blood_lineages <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "blood_lineages")
blood_lineages_colors <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "colors")
blood_lineages <- merge(x = blood_lineages, y = blood_lineages_colors, all.x = T, by = c("CellType"))
rownames(blood_lineages) <- blood_lineages$CellMarker
scale_color_manual_colors_sortedbyname <- blood_lineages_colors[order(blood_lineages_colors$CellType), "colorcode"]
scale_color_manual_colors_sortedbyname <- colors_lineages[1,] # select only the ones here needed
patients_summary <- read.xlsx(xlsxFile = "source/Clinical/Patient Summary.xlsx", sheet = "summary")

# folder management
dir.create(file.path("analyses/", results_folder_name), showWarnings = FALSE)
dir.create(file.path(paste0("analyses/",results_folder_name), analysis_folder_date), showWarnings = FALSE)
dir.create(file.path(paste0("analyses/", results_folder_name, analysis_folder_date), this_run_foldername), showWarnings = FALSE)
dest_folder <- paste0("analyses/", results_folder_name, analysis_folder_date, "/", this_run_foldername)
root_folder <- paste0("analyses/", results_folder_name)
mld_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial MLD/"
bthal_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial bThalassemia/"
was_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial WAS/"

# collect data from each trial
# ---------------- MLD ----------------
mld_shared34_stats <- read.csv(file = paste(mld_wd, "analyses/12.sharing/202112/20220121_MLD_CD34-Output_SuperGroup.tsv.gz", sep = ""), 
                               header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NA", ''))
mld_shared34_stats_details <- as.data.frame( t(as.data.frame(lapply(strsplit(as.character(mld_shared34_stats$g2), '_', fixed = T), function(x) {c(x)}))) )
names(mld_shared34_stats_details) <- c("SubjectID", "CellMarker", "Tissue", "TimepointMonths")
mld_shared34_stats_details$FollowUp <- as.numeric(as.character(mld_shared34_stats_details$TimepointMonths))
mld_shared34_stats <- cbind(mld_shared34_stats, mld_shared34_stats_details)
mld_shared34_stats$StudyID <- "MLD"

# ---------------- BTHAL ----------------
bthal_shared34_stats <- read.csv(file = paste(bthal_wd, "analyses/13.sharing/202112/20220121_BTHAL_CD34-Output_AggregationGroup.tsv.gz", sep = ""), 
                               header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NA", ''))
bthal_shared34_stats_details <- as.data.frame( t(as.data.frame(lapply(strsplit(as.character(bthal_shared34_stats$g2), '_', fixed = T), function(x) {c(x)}))) )
names(bthal_shared34_stats_details) <- c("SubjectID", "CellMarker", "Tissue", "TimepointMonths")
bthal_shared34_stats_details$FollowUp <- as.numeric(as.character(bthal_shared34_stats_details$TimepointMonths))
bthal_shared34_stats <- cbind(bthal_shared34_stats, bthal_shared34_stats_details)
bthal_shared34_stats$StudyID <- "BTHAL"

# ---------------- WAS ----------------
was_shared34_stats <- read.csv(file = paste(was_wd, "analyses/13.sharing/202112/20220121_WAS_CD34-Output_by-SuperGroup.tsv.gz", sep = ""), 
                                 header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NA", ''))
was_shared34_stats_details <- as.data.frame( t(as.data.frame(lapply(strsplit(as.character(was_shared34_stats$g2), '_', fixed = T), function(x) {c(x)}))) )
names(was_shared34_stats_details) <- c("SubjectID", "CellMarker", "Tissue", "TimepointMonths")
was_shared34_stats_details$FollowUp <- as.numeric(as.character(was_shared34_stats_details$TimepointMonths))
was_shared34_stats <- cbind(was_shared34_stats, was_shared34_stats_details)
was_shared34_stats$StudyID <- "WAS"


### ============= combine data ============= 

common_labels <- Reduce(intersect, list(colnames(mld_shared34_stats), colnames(bthal_shared34_stats), colnames(was_shared34_stats)))
allpatients_shared34_stats <- do.call("rbind", list(mld_shared34_stats[common_labels], was_shared34_stats[common_labels], bthal_shared34_stats[common_labels]))
allpatients_shared34_stats <- merge(x = allpatients_shared34_stats, y = patients_summary, by = c("SubjectID"), all.x = T)
# add age group
allpatients_shared34_stats <- allpatients_shared34_stats %>%
  mutate(AgeGroup = case_when(
    TreatmentAge.years <= 2 ~"0-2",
    TreatmentAge.years <= 15 ~"2-15",
    TreatmentAge.years > 15 ~"30+"
  ))
age_list <- c("0-2", "2-15", "30+")
allpatients_shared34_stats$AgeGroup_Sorted <- factor(allpatients_shared34_stats$AgeGroup, levels = age_list)
# refactor study id for plot order
allpatients_shared34_stats$Study <- factor(allpatients_shared34_stats$StudyID, levels = study_list)
# add cell marker info
allpatients_shared34_stats <- merge(x = allpatients_shared34_stats, y = blood_lineages, by = c("CellMarker"), all.x = T)
allpatients_shared34_stats <- allpatients_shared34_stats[which(allpatients_shared34_stats$CellType %in% colnames(scale_color_manual_colors_sortedbyname)),]
# scale zscore
allpatients_shared34_stats <- allpatients_shared34_stats %>% group_by(SubjectID, Tissue, TimepointMonths) %>% mutate(on_g1_zscaled=scale(on_g1))
allpatients_shared34_stats[is.na(allpatients_shared34_stats)] <- NA
# remaning with more clear col names
allpatients_shared34_stats$PercNIS_onOverallSharedCD34BM <- allpatients_shared34_stats$on_g1
allpatients_shared34_stats$PercNIS_onOverallSharedCD34BM_zscaled <- allpatients_shared34_stats$on_g1_zscaled
# write data file
write.table(x = allpatients_shared34_stats, file = data_filename, sep = "\t", quote = FALSE, col.names = TRUE, row.names = T, na = "")

# plot all data

# general data
p_lines_mbte <- ggplot(data = allpatients_shared34_stats, aes(x = FollowUp, y = on_g1, color = CellType, fill = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  stat_smooth(formula = y ~ log(x), level = 0.4, size = 2) +
  facet_grid(. ~ Study, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats$FollowUp, na.rm = T), 12), limits = c(0, followup_limit_crosstrial) ) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34.")


p_lines_mbte_identicalFU <- ggplot(data = allpatients_shared34_stats, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  facet_grid(. ~ Study, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats$FollowUp, na.rm = T), 12), limits = c(0, followup_limit_crosstrial) ) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34.")


p_lines_mbte_zscore <- ggplot(data = allpatients_shared34_stats, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM_zscaled, color = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  facet_grid(. ~ Study, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats$FollowUp, na.rm = T), 12), limits = c(0, followup_limit_crosstrial) ) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), 
       x = "FollowUp months after GT", y = "Perc. of shared CD34 IS (Z-score)", 
       subtitle = "All clinical trials included, all patients included. Percentages are scaled (Z-score) by Patient, Tissue and Follow-up.")

p_lines_mbte_zscore_agegr <- ggplot(data = allpatients_shared34_stats, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM_zscaled, color = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  facet_grid(AgeGroup_Sorted ~ Study, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats$FollowUp, na.rm = T), 12), limits = c(0, followup_limit_crosstrial) ) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), 
       x = "FollowUp months after GT", y = "Perc. of shared CD34 IS (Z-score)", 
       subtitle = "All clinical trials included, all patients included. Percentages are scaled (Z-score) by Patient, Tissue and Follow-up.")

p_lines_mbte_agegr_nopoints_identicalFU <- ggplot(data = allpatients_shared34_stats, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  facet_grid(AgeGroup_Sorted ~ Study, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats$FollowUp, na.rm = T), 12), limits = c(0, followup_limit_crosstrial) ) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34.")

p_lines_mbte_agegr_nopoints_identicalFU_freey <- ggplot(data = allpatients_shared34_stats, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  facet_grid(AgeGroup_Sorted ~ Study, scales = "free", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats$FollowUp, na.rm = T), 12), limits = c(0, followup_limit_crosstrial) ) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34.")


pdf(file = paste0(plot_filename_prefix, "percentage.pdf", sep = ""), height=4.5, width=10)
print(p_lines_mbte_identicalFU)
dev.off()
png(file = paste0(plot_filename_prefix, "percentage.wrap_trials.png", sep = ""), height=4.5, width=10, units = "in", res = 300)
print(p_lines_mbte_identicalFU)
dev.off()

pdf(file = paste0(plot_filename_prefix, "zscore.wrap_trials.pdf", sep = ""), height=4.5, width=10)
print(p_lines_mbte_zscore)
dev.off()
png(file = paste0(plot_filename_prefix, "zscore.wrap_trials.png", sep = ""), height=4.5, width=10, units = "in", res = 300)
print(p_lines_mbte_zscore)
dev.off()

pdf(file = paste0(plot_filename_prefix, "percentage.wrap_trials_age.pdf", sep = ""), height=7, width=10)
print(p_lines_mbte_agegr_nopoints_identicalFU)
dev.off()
png(file = paste0(plot_filename_prefix, "percentage.wrap_trials_age.png", sep = ""), height=7, width=10, units = "in", res = 300)
print(p_lines_mbte_agegr_nopoints_identicalFU)
dev.off()

pdf(file = paste0(plot_filename_prefix, "zscore.wrap_trials_age.pdf", sep = ""), height=7, width=10)
print(p_lines_mbte_zscore_agegr)
dev.off()
png(file = paste0(plot_filename_prefix, "zscore.wrap_trials_age.png", sep = ""), height=7, width=10, units = "in", res = 300)
print(p_lines_mbte_zscore_agegr)
dev.off()

