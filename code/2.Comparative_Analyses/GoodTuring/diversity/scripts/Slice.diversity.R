####################################################
####  R code for Trial comparison               ####
####################################################

# NOTE: NELLA PRIMA PARTE TROVERETE DEL CODICE DI UTILITA', CHE VI PORTANO A VEDERE COME SI ARRIVA AL DATO (PRESO IN ORIGINE DEI GZ CONDIVISI). QUI SOTTO TROVERETE IL PUNTO DA CUI PARTIRE (CON UN COMMENTO IN MAIUSCOLO COME QUESTO), IMPORTANTO IL FILE EXCEL E POI FACENDONE I PLOT SOTTOSTANTI. VI HO TAGLIATO IL RESTO.

# library(forcats)
# library(annotatr)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(gplots)
# library(xlsx)
library(ggplot2)
library(scales) 
library(splines)
library(gridExtra)
library(stringr)
library(sqldf)
library(plyr)
library(stringr)
library(psych)
library(reshape2)
library(dplyr)
# library(rGREAT)
library(trackViewer)
library(GenomicAlignments)
library(openxlsx)
library(Hmisc)
library(openxlsx)
library(circlize)
library(ggpubr)
library(rstatix)
library(ggbreak)

###############################################################
## globals
###############################################################
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
blood_lineages_file <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial MLD/blood_lineages_update.xlsx"
adults_to_exclude <- c("BTHAL001", "BTHAL003", "BTHAL004")
colors_lineages <- data.frame("B" = "dodgerblue3",
                              ## "CD34" = "green4",
                              "Erythroid" = "red3",
                              "Myeloid" = "orange",
                              "T" = "deepskyblue")




###############################################################
## Diversity (rescaled)
###############################################################
markerlist <- c("CD13", "CD14", "CD15", "CD19", "CD3", "CD34", "GLY", "CD36", "GLYA")
h_base_folder <- paste0("analyses/diversity/")
dir.create(file.path(getwd(), h_base_folder), showWarnings = FALSE)
dir.create(file.path(h_base_folder, analysis_folder_date), showWarnings = FALSE)

mld_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial MLD/"
bthal_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial bThalassemia/"
was_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial WAS/"

# import all data
mld_h <- read.csv(file = paste0(mld_wd, "analyses/01.stats/", analysis_folder_date, "/20211202_MLD_descriptive-stats_SubjectID_CellMarker_Tissue_TimepointMonths.tsv.gz"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
mld_h$ClinicalTrial <- "MLD"
mld_h$Hindex <- ifelse(!is.na(mld_h$fragmentEstimate_sum_shannon), mld_h$fragmentEstimate_sum_shannon, mld_h$seqCount_sum_shannon)
colname_to_convertto_cellmarker <- "SuperGroup"
if (colname_to_convertto_cellmarker %in% colnames(mld_h) & !("CellMarker" %in% colnames(mld_h))) {
  mld_h$CellMarker <- mld_h[,colname_to_convertto_cellmarker]
}

was_h <- read.csv(file = paste0(was_wd, "analyses/01.stats/", analysis_folder_date, "/20211210_WAS_descriptive-stats_SubjectID_CellMarker_Tissue_TimepointMonths_LAM-only.tsv.gz"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
was_h$ClinicalTrial <- "WAS"
was_h$Hindex <- was_h$seqCount_sum_shannon # seqCount_sum_shannon seqCount_sum_simpson seqCount_sum_invsimpson
colname_to_convertto_cellmarker <- "SuperGroup"
if (colname_to_convertto_cellmarker %in% colnames(was_h) & !("CellMarker" %in% colnames(was_h))) {
  was_h$CellMarker <- was_h[,colname_to_convertto_cellmarker]
}

bthal_h <- read.csv(file = paste0(bthal_wd, "analyses/00.stats/", analysis_folder_date, "/20211206_BTHAL_descriptive-stats_SubjectID_AggregationGroup_Tissue_TimepointMonths.tsv.gz"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
bthal_h$ClinicalTrial <- "BTHAL"
bthal_h$Hindex <- bthal_h$fragmentEstimate_sum_shannon
colname_to_convertto_cellmarker <- "AggregationGroup"
if (colname_to_convertto_cellmarker %in% colnames(bthal_h) & !("CellMarker" %in% colnames(bthal_h))) {
  bthal_h$CellMarker <- bthal_h[,colname_to_convertto_cellmarker]
}

# merge
common_labels <- Reduce(intersect, list(colnames(was_h), colnames(mld_h), colnames(bthal_h)))
all_results_h <- Reduce(rbind, list(mld_h[common_labels], was_h[common_labels], bthal_h[common_labels]))

all_results_h$Study <- factor(all_results_h$ClinicalTrial, levels = study_list)
all_results_h[is.na(all_results_h)] <- NA
# all_results_h$Hindex <- all_results_h$fragmentEstimate_sum_shannon
# all_results_h$gdf_names <- rownames(all_results_h)
all_results_h$TimePoint <- all_results_h$TimepointMonths

# normalize as Nat Medicine 21 LiBIS
lam_multiplication_factor <- 2.136
all_results_h$Hindex_normalized <- ifelse( (!is.na(all_results_h$`ng DNA corrected_sum`) & (all_results_h$`ng DNA corrected_sum`) > 1),
                                        ifelse(all_results_h$PCRMethod %in% c("LAM", "LAM-PCR", "LAM-PCR|SLiM", "SLiM|LAM-PCR", "LAM|SLIM", "SLiM|LAM"), 
                                               all_results_h$Hindex * lam_multiplication_factor / log(all_results_h$`ng DNA corrected_sum`), 
                                               all_results_h$Hindex / log(all_results_h$`ng DNA corrected_sum`)
                                            ),
                                            NA
                                          )

# use blood lineages\
# add info of lineages
library(openxlsx)
blood_lineages <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "blood_lineages")
# blood_lineages_colors <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "colors_noery")
blood_lineages_colors <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "colors")
blood_lineages <- merge(x = blood_lineages, y = blood_lineages_colors, all.x = T, by = c("CellType"))
rownames(blood_lineages) <- blood_lineages$CellMarker
all_results_h <- merge(x = all_results_h, y = blood_lineages, by = c("CellMarker"), all.x = T)
# rownames(all_results_h) <- as.character(all_results_h$gdf_names)
all_results_h$CellLineage <- factor(all_results_h$CellType, levels = celltype_list)
all_results_h$Jindex <- all_results_h$Hindex / (log(all_results_h$nIS))

scale_color_manual_colors_sortedbyname <- blood_lineages_colors[order(blood_lineages_colors$CellType), "colorcode"]

# do slice
all_results_h_slice_ngdna <- all_results_h[which(all_results_h$CellMarker %in% markerlist 
                                           & all_results_h$Tissue %in% c("PB", "BM") 
                                           & !is.na(all_results_h$Hindex)
                                           & all_results_h$nIS >= 10
                                           & all_results_h$`ng DNA corrected_sum` >= 10
                                           & all_results_h$TimePoint > 0),]
all_results_h_slice_ngdna_scaled <- all_results_h_slice_ngdna %>% group_by(SubjectID, Tissue, TimePoint) %>% mutate(h_zscaled=scale(Hindex))
all_results_h_slice_ngdna_scaled <- all_results_h_slice_ngdna_scaled %>% group_by(SubjectID, Tissue, TimePoint) %>% mutate(h_normalized_zscaled=scale(Hindex_normalized))

all_results_h_slice <- all_results_h[which(all_results_h$CellMarker %in% markerlist 
                                           & all_results_h$Tissue %in% c("PB", "BM") 
                                           & !is.na(all_results_h$Hindex)
                                           & all_results_h$nIS >= 10
                                           # & all_results_h$`ng DNA corrected_sum` >= 10
                                           & all_results_h$TimePoint > 0),]

all_results_h_slice_scaled <- all_results_h_slice %>% group_by(SubjectID, Tissue, TimePoint) %>% mutate(h_zscaled=scale(Hindex))
all_results_h_slice_scaled <- all_results_h_slice_scaled %>% group_by(SubjectID, Tissue, TimePoint) %>% mutate(h_normalized_zscaled=scale(Hindex_normalized))
all_results_h_slice_scaled_homo <- all_results_h_slice_scaled[which(all_results_h_slice_scaled$PCRMethod == "SLiM" | all_results_h_slice_scaled$ClinicalTrial == "WAS"),]


##### -------- QUESTO FILE QUI SOTTO E' PRESENTE NELLA CARTELLA CONDIVISA QUINDI POTETE PARTIRE DA QUESTO ----------- #####

# write results in a file
write.xlsx(x = all_results_h_slice_scaled, file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.xlsx", sep = ""))
write.xlsx(x = all_results_h_slice_ngdna_scaled, file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.ngdna_notNA.xlsx", sep = ""))


plot_hindex_homogeneouspcr <-
  ggplot(data = all_results_h_slice[which(all_results_h_slice$PCRMethod == "SLiM" | all_results_h_slice$ClinicalTrial == "WAS"),], aes(x = TimePoint, y = Hindex, fill = CellType, color = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(alpha = 0.2, size = 2) +
  # geom_jitter(alpha = 0.8, size = 4) +
  # geom_line(size=3, alpha = .7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  scale_x_continuous(breaks = seq(0, max(all_results_h_slice_ngdna$TimePoint, na.rm = T), 12)) +
  facet_grid( Tissue ~ Study, scales = "free_x", space = "free") +
  theme_bw() +
  # scale_y_log10() + 
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  # facet_grid(ProjectID ~ ., scales = "free_y", space = "free") +
  # facet_grid(ProjectID ~ ., scales = "free_y") +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Diversity index (Shannon)"), subtitle = "Patients with homogeneous PCR methods over time (no mixed patients).",
       x = "Months after gene therapy", 
       y = "Shannon diversity index (H')", 
       colour = "CellMarker", fill = "CellMarker")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.HomoPCR.pdf", sep = ""), height=8, width=12)
plot(plot_hindex_homogeneouspcr)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.HomoPCR.png", sep = ""), height=8, width=12, units = "in", res = 300)
plot(plot_hindex_homogeneouspcr)
dev.off()


plot_hindex_homogeneouspcr_zscaled <-
  ggplot(data = all_results_h_slice_scaled_homo, 
         aes(x = TimePoint, y = h_zscaled, fill = CellType, color = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(alpha = 0.2, size = 2) +
  # geom_jitter(alpha = 0.8, size = 4) +
  # geom_line(size=3, alpha = .7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  scale_x_continuous(breaks = seq(0, max(all_results_h_slice_scaled$TimePoint, na.rm = T), 12)) +
  facet_grid( Tissue ~ Study, scales = "free_x", space = "free") +
  theme_bw() +
  # scale_y_log10() + 
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  # facet_grid(ProjectID ~ ., scales = "free_y", space = "free") +
  # facet_grid(ProjectID ~ ., scales = "free_y") +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Diversity index (Shannon)"), subtitle = "Patients with homogeneous PCR methods over time (no mixed patients). Zscore normalization by Patient, Tissue, TimePoint.",
       x = "Months after gene therapy", 
       y = "Shannon diversity index (H')", 
       colour = "CellMarker", fill = "CellMarker")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.HomoPCR.Zscore.pdf", sep = ""), height=8, width=12)
plot(plot_hindex_homogeneouspcr_zscaled)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.HomoPCR.Zscore.png", sep = ""), height=8, width=12, units = "in", res = 300)
plot(plot_hindex_homogeneouspcr_zscaled)
dev.off()
