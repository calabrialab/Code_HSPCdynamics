####################################################
####  R code for Trial comparison               ####
####################################################
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




source("/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/TIGET/Workbench/isatk/script/R/isa_utils_functions.R")

# acquire DB
onco_db <- read.csv("/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/TIGET/Workbench/isatk/script/R/publicdb/201806_uniprot-Proto-oncogene.tsv", header=TRUE, fill=T, sep='\t', check.names = FALSE)
tumsup_db <- read.csv("/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/TIGET/Workbench/isatk/script/R/publicdb/201806_uniprot-Tumor-suppressor.tsv", header=TRUE, fill=T, sep='\t', check.names = FALSE)
ccgd_db <- read.csv("/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/TIGET/Workbench/isatk/script/R/publicdb/ccgd/ccgd_export_2020-02-26.human.csv", header=TRUE, fill=T, sep=',', check.names = FALSE)
cancermine_db <- read.csv("/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/TIGET/Workbench/isatk/script/R/publicdb/cancermine/cancermine_collated.201910.tsv", header=TRUE, fill=T, sep='\t', check.names = FALSE)

# test file of IS
oncots_df_touse <- loadOncoTSgenes(species = "human")
oncots_df_touse[oncots_df_touse %in% c("<NA>", "NA", "")] <- NA
oncots_df_touse[is.na(oncots_df_touse)] <- NA
clinical_relevant_suspicious_genes <- data.frame("GeneName" = c("DNMT3A", "TET2", "ASXL1", "JAK2", "CBL", "TP53"), 
                                                 "ClinicalRelevance" = TRUE,
                                                 "DOIReference" = "https://doi.org/10.1182/blood-2018-01-829937")
rownames(clinical_relevant_suspicious_genes) <- clinical_relevant_suspicious_genes$GeneName

known_clinical_oncogenes <- data.frame("GeneName" = c("MECOM", "CCND2", "TAL1", "LMO2", "HMGA2"),
                                       "KnownClonalExpension" = TRUE)
rownames(known_clinical_oncogenes) <- known_clinical_oncogenes$GeneName

cancermine_db$GeneName <- cancermine_db$gene_normalized
cancermine_db_cast <- dcast(data = cancermine_db, GeneName ~ role, value.var = "citation_count", fun.aggregate = mean)
names(cancermine_db_cast) <- c("GeneName", "CancerMine_Driver", "CancerMine_Oncogene", "CancerMine_Tumor_Suppressor")
# rownames(cancermine_db) <- cancermine_db$GeneName
# patients_summary_byMarkerTissueFU <- dcast(data = stats, SubjectID + CellMarker + Tissue ~ FollowUp, value.var = "NumIS")

# collect totals
total_oncots <- sqldf("select count(GeneName) as ONCOTS_N_Genes, count(distinct TumorSuppressor) as N_Genes_TumorSuppressor, count(distinct OncoGene) as N_Genes_OncoGene from oncots_df_touse where 1")

total_cancermine <- sqldf("select cm.role as CANCERMINE_Role, count(gene_normalized) as CANCERMINE_N_Genes from cancermine_db as cm where 1 group by cm.role")
total_cancermine_cast <- dcast(data = total_cancermine, "from_CGC" ~ CANCERMINE_Role, value.var = "CANCERMINE_N_Genes", fun.aggregate = mean)
names(total_cancermine_cast) <- c(paste0("FoundIn_CancerMine_", colnames(total_cancermine_cast)))

total_ccgd <- sqldf("select count(HumanName) as CCGD_N_Genes from ccgd_db where HumanName not in ('', 'N/A') and Rank not in ('', 'N/A')")

total_ccgd_and_cosmic <- sqldf("select COSMIC, count(HumanName) as CCGD_N_Genes_COSMIC from ccgd_db where HumanName not in ('', 'N/A') and Rank not in ('', 'N/A') group by COSMIC")
total_ccgd_and_cosmic_cast <- dcast(data = total_ccgd_and_cosmic, "from_CGC" ~ COSMIC, value.var = "CCGD_N_Genes_COSMIC", fun.aggregate = mean)
names(total_ccgd_and_cosmic_cast) <- c(paste0("FoundIn_CCGD_COSMIC_", colnames(total_ccgd_and_cosmic_cast)))

total_ccgd_and_cgc <- sqldf("select CGC, count(HumanName) as CCGD_N_Genes_CGC from ccgd_db where HumanName not in ('', 'N/A') and Rank not in ('', 'N/A') group by CGC")
total_ccgd_and_cgc_cast <- dcast(data = total_ccgd_and_cgc, "from_CGC" ~ CGC, value.var = "CCGD_N_Genes_CGC", fun.aggregate = mean)
names(total_ccgd_and_cgc_cast) <- c(paste0("FoundIn_CCGD_CGC_", colnames(total_ccgd_and_cgc_cast)))


totals_from_db <- cbind(total_ccgd, total_ccgd_and_cosmic_cast, total_ccgd_and_cgc_cast, total_oncots, total_cancermine_cast)
###############################################################
## main
###############################################################


##### =============================================================== #####
##### ------------- Metadata ------------ @ ------------------------- #####
##### =============================================================== #####

# association file
metadataf <- paste0("source/AF/20211129.AF.MLD-WAS-BTHAL-MPSIH.tsv", sep = "") 
label_metadata <- read.csv(metadataf, header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", "ND", ""))
rownames(label_metadata) <- label_metadata$CompleteAmplificationID

patients_summary <- read.xlsx(xlsxFile = "source/Clinical/Patient Summary.xlsx", sheet = "summary")

###############################################################
## DATA
###############################################################

mld_cis_avg <- read.csv(file = "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial MLD/analyses/06.cis/202004/202004.MLD.AllPatients.CIS.results_invivo.avgByPatient.tsv", header=T, fill=T, check.names = FALSE, sep = '\t')
top_hit_genes <- mld_cis_avg[order(-mld_cis_avg$avg_neg_zscore_minus_log2_integration_freq_withtolerance), c("GeneName")]
top_50_hit_genes <- top_hit_genes[1:50]
top_20_hit_genes <- top_hit_genes[1:20]

###############################################################
## CIRCOS
###############################################################
# alm_reads_for_rearrangements_withstatsbyread

human_cytoband <- read.cytoband(species = "hg19")$df
vector_cytoband <- read.csv(file = "source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.cytoband", header=F, fill=T, check.names = FALSE, sep = '\t')
cytoband_rbind <- rbind(human_cytoband)
cytoband <- read.cytoband(cytoband_rbind)
# circos.initializeWithIdeogram(cytoband_df)

extend_chromosomes = function(bed, chromosome, prefix = "zoom_") {
  zoom_bed = bed[bed[[1]] %in% chromosome, , drop = FALSE]
  zoom_bed[[1]] = paste0(prefix, zoom_bed[[1]])
  rbind(bed, zoom_bed)
}
# circos.clear()


cytoband_df = cytoband$df
chromosome = cytoband$chromosome

xrange = cytoband$chr.len
# normal_chr_index = c(1:22,24:25)
# zoomed_chr_index = 23


bed_colnames <- c("chr","start","end","value1", "score", "strand")
bed_was1 <- read.csv(file = "source/WAS/bed/WAS1001.SC-FE.sampleFiltered.noDup.bed", header=F, fill=T, check.names = FALSE, sep = '\t')
names(bed_was1) <- bed_colnames
bed_mld1 <- read.csv(file = "source/MLD/bed/IS_matrix_classic_strand_specific_method_mld_202004_MLD01.no0.annotated.bed", header=F, fill=T, check.names = FALSE, sep = '\t')
names(bed_mld1) <- bed_colnames
bed_bthal1 <- read.csv(file = "source/BTHAL/bed/BTHAL001.collDate.filt.gdf.byPMTF.SS.invivo.bed", header=F, fill=T, check.names = FALSE, sep = '\t')
names(bed_bthal1) <- bed_colnames

bed_hg19_genes_fullbody <- read.csv(file = "/Users/calabria.andrea/Dropbox (HSR Global)/TIGET/Workspace/Andrea/RefGene/hg19.ucsc.gemini.txGenes.bed", header=F, fill=T, check.names = FALSE, sep = '\t')
names(bed_hg19_genes_fullbody) <- bed_colnames
bed_hg19_genes <- read.csv(file = "/Users/calabria.andrea/Dropbox (HSR Global)/TIGET/Workspace/Andrea/RefGene/ucsc.hg19.knowngene.2021.bed", header=F, fill=T, check.names = FALSE, sep = '\t')

bed_random <- generateRandomBed(nr = 5000, fun = function(k) sample(letters, k, replace = TRUE))

bed_top_genes <- bed_hg19_genes_fullbody[which(bed_hg19_genes_fullbody$value1 %in% top_50_hit_genes),]
bed_top_genes <- bed_hg19_genes_fullbody[which(bed_hg19_genes_fullbody$value1 %in% top_20_hit_genes),]
# names(bed_chimera_hg19_annot) <- bed_col_names
# bed_chimera_hg19_annot_samples <- as.data.frame( t(as.data.frame(lapply(strsplit(as.character(bed_chimera_hg19_annot$ISid_hg19), '@', fixed = T), function(x) {c(x)}))) )
# names(bed_chimera_hg19_annot_samples) <- c("sample_label", "name")
# bed_chimera_hg19_annot <- cbind(bed_chimera_hg19_annot, bed_chimera_hg19_annot_samples)

# circos.clear()

# circos.par("start.degree" = 90)
# circos.initializeWithIdeogram(species = "hg19")
# circos.genomicIdeogram()
# circos.genomicDensity(list(bed_mld1, bed_was1, bed_bthal1), col = c("navyblue", "firebrick", "forestgreen"))
# circos.genomicDensity(bed_2, col = c("#FF000080"))

theight <- 0.1
col_fun = colorRamp2(c(-1, 0, 1), c("white", "yellow", "red"))

png(file = paste0("analyses/", analysis_folder_date, "/circos/", analysis_folder_date, ".genomic.density.withGenes.png", sep = ""), height=8, width=8, units = "in", res = 300)
circos.par("start.degree" = 90)
circos.initializeWithIdeogram(species = "hg19")
# circos.genomicIdeogram()
# circos.genomicDensity(bed_mld1, bed_was1, bed_bthal1), col = c("navyblue", "firebrick", "forestgreen"))
# circos.genomicHeatmap(genomicDensity(bed_hg19_genes), col = col_fun)
# circos.genomicLabels(bed, labels.column = 4)
circos.genomicDensity(bed_hg19_genes, col = c("gray70"), track.height = 0.1)
circos.genomicDensity(bed_mld1, col = c("navyblue"), track.height = theight)
circos.genomicDensity(bed_was1, col = c("firebrick"), track.height = theight)
circos.genomicDensity(bed_bthal1, col = c("forestgreen"), track.height = theight)
circos.genomicLabels(bed_top_genes[1:20,1:4], labels.column = 4, side = "inside")

# circos.initializeWithIdeogram(extend_chromosomes(cytoband_df, c("chrV")), chromosome.index = c("chr1", "chrV"),
#                               sector.width = sector.width)
# circos.genomicLink(bed_1, bed_2, col = rand_color(nrow(bed_1), transparency = 0.8), border = NA)
dev.off()

circos.clear()
pdf(file = paste0("analyses/", analysis_folder_date, "/circos/", analysis_folder_date, ".genomic.density.withGenes.pdf", sep = ""), height=8, width=8)
circos.par("start.degree" = 90)
circos.initializeWithIdeogram(species = "hg19")
# circos.genomicIdeogram()
# circos.genomicDensity(bed_mld1, bed_was1, bed_bthal1), col = c("navyblue", "firebrick", "forestgreen"))
circos.genomicDensity(bed_hg19_genes, col = c("gray70"), track.height = 0.08)
circos.genomicDensity(bed_mld1, col = c("navyblue"), track.height = theight)
circos.genomicDensity(bed_was1, col = c("firebrick"), track.height = theight)
circos.genomicDensity(bed_bthal1, col = c("forestgreen"), track.height = theight)
circos.genomicLabels(bed_top_genes[1:20,1:4], labels.column = 4, side = "inside")
dev.off()

circos.clear()
pdf(file = paste0("analyses/circos/", analysis_folder_date, analysis_folder_date, ".genomic.density.withoutGenes.newcolors.pdf", sep = ""), height=8, width=8)
circos.par("start.degree" = 90)
circos.initializeWithIdeogram(species = "hg19")
# circos.genomicIdeogram()
# circos.genomicDensity(bed_mld1, bed_was1, bed_bthal1), col = c("navyblue", "firebrick", "forestgreen"))
# circos.genomicDensity(bed_hg19_genes, col = c("gray70"), track.height = 0.08)
circos.genomicDensity(bed_mld1, col = trials_colors[1], track.height = theight)
circos.genomicDensity(bed_was1, col = trials_colors[2], track.height = theight)
circos.genomicDensity(bed_bthal1, col = trials_colors[3], track.height = theight)
circos.genomicLabels(bed_top_genes[1:20,1:4], labels.column = 4, side = "inside")
dev.off()

circos.clear()
pdf(file = paste0("analyses/", analysis_folder_date, "/circos/", analysis_folder_date, ".genomic.density.withoutGenes.2.pdf", sep = ""), height=8, width=8)
circos.par("start.degree" = 90)
circos.initializeWithIdeogram(species = "hg19")
# circos.genomicIdeogram()
# circos.genomicDensity(bed_mld1, bed_was1, bed_bthal1), col = c("navyblue", "firebrick", "forestgreen"))
# circos.genomicDensity(bed_hg19_genes, col = c("gray70"), track.height = 0.08)
circos.genomicLabels(bed_top_genes[1:20,1:4], labels.column = 4, side = "outside")
circos.genomicDensity(bed_mld1, col = c("navyblue"), track.height = theight)
circos.genomicDensity(bed_was1, col = c("firebrick"), track.height = theight)
circos.genomicDensity(bed_bthal1, col = c("forestgreen"), track.height = theight)
dev.off()

###############################################################
## UMAP
###############################################################





###############################################################
## CIS
###############################################################

cis_base_folder <- paste0("analyses/cis/")
dir.create(file.path(cis_base_folder, analysis_folder_date), showWarnings = FALSE)
# dir.create(file.path(paste0("analyses/", results_folder_name), analysis_folder_date), showWarnings = FALSE)


# import all data
was_cis <- read.csv(file = "source/WAS/cis/202103.WAS.CIS.results_invivo.tsv.gz", header=T, fill=T, check.names = FALSE, sep = '\t')
was_cis$neg_zscore_minus_log2_integration_freq_withtolerance <- was_cis$neg_zscore_minus_log2_int_freq_tolerance
mld_cis <- read.csv(file = "source/MLD/cis/202004.MLD.CIS.results_invivo.tsv.gz", header=T, fill=T, check.names = FALSE, sep = '\t')
bthal_cis <- read.csv(file = "source/BTHAL/cis/202003.BTHAL.CIS.results_invivo.tsv.gz", header=T, fill=T, check.names = FALSE, sep = '\t')

mld_cis <- acquire_CIS_infolder(folder = "source/MLD/cis/202110/", prefix = "MLD", suffix = ".tsv.gz", study = "MLD")
was_cis <- acquire_CIS_infolder(folder = "source/WAS/cis/202110/", prefix = "WAS", suffix = ".tsv.gz", study = "WAS")
bthal_cis <- acquire_CIS_infolder(folder = "source/BTHAL/cis/202110/", prefix = "BTHAL", suffix = ".tsv.gz", study = "BTHAL")

common_labels <- Reduce(intersect, list(colnames(was_cis), colnames(mld_cis), colnames(bthal_cis)))
all_results_cis <- Reduce(rbind, list(was_cis[common_labels], mld_cis[common_labels], bthal_cis[common_labels]))

# add missing annotations
all_results_cis <- merge(x = all_results_cis, y = oncots_df_touse, by = c("GeneName"), all.x = T)
all_results_cis <- merge(x = all_results_cis, y = clinical_relevant_suspicious_genes, by = c("GeneName"), all.x = T)
all_results_cis <- merge(x = all_results_cis, y = known_clinical_oncogenes, by = c("GeneName"), all.x = T)

all_results_cis$minus_log_p_fdr <- -log(all_results_cis$tdist_fdr, base = 10)
significance_threshold_minus_log_p <- -log(0.05, base = 10)
all_results_cis$KnownGeneClass <- ifelse(!is.na(all_results_cis$Onco1_TS2), (ifelse((all_results_cis$Onco1_TS2) == 1, "OncoGene", "TumSuppressor")), "Other")
all_results_cis$CriticalForInsMut <- ifelse(!is.na(all_results_cis$KnownClonalExpension), "True", "False")
all_results_cis[is.na(all_results_cis)] <- NA
all_results_cis$positive_outlier_and_significant <- ifelse((!is.na(all_results_cis$tdist_fdr) & all_results_cis$tdist_fdr < 0.05), TRUE, FALSE)


# configs
annotation_threshold_ontots <- 0.5
annotation_threshold_ontots_log <- -log(annotation_threshold_ontots, base = 10)
annotation_cols_to_get <- c("GeneName", "OncoGene", "TumorSuppressor", "Onco1_TS2", "ClinicalRelevance", "DOIReference", "KnownGeneClass", "KnownClonalExpension", "CriticalForInsMut")
annotation_cols_to_get_asSQLstring <- paste(annotation_cols_to_get, collapse = ", ")
annotation_cols_to_get_onheatmap <- c("KnownGeneClass", "ClinicalRelevance", "CriticalForInsMut")

# write file in vivo
write.table(x = all_results_cis, file = paste(cis_base_folder, analysis_folder_date, "/", analysis_folder_date, ".CIS.results_invivo.tsv", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = T, na = '')

# all_results_cis_onlypositive <- all_results_cis[which(!is.na(all_results_cis$tdist_positive_and_correctedEM) & all_results_cis$tdist_positive_and_correctedEM > 0 & all_results_cis$n_IS_perGene > 3),]
# write.table(x = all_results_cis_onlypositive, file = paste(cis_base_folder, analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".CIS.results_invivo.onlyPositive.tsv", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = T)


# plot_cis_alldetails <- ggplot(data = all_results_cis, aes(y = minus_log_p, x = neg_zscore_minus_log2_int_freq_tolerance, color = KnownGeneClass), na.rm = T, se = TRUE) + 
#   geom_point(aes(size = average_TxLen), alpha = .5) +
#   # geom_hex(bins = 60) +
#   # geom_hex() +
#   geom_hline(yintercept = significance_threshold_minus_log_p, color='coral', size=1, show.legend = T, linetype="dotted") +
#   scale_y_continuous(limits = c(0, max(c((significance_threshold_minus_log_p + 0.5), max(all_results_cis$minus_log_p) )) ) ) + 
#   facet_wrap( ~ SubjectID, ncol = 4) +
#   geom_label_repel(
#     data = subset(all_results_cis, !is.na(Onco1_TS2) & minus_log_p > annotation_threshold_ontots_log),
#     # data = subset(all_results_cis, !is.na(Onco1_TS2) & minus_log2_integration_freq_withtolerance > 2.8),
#     aes(label = GeneName, fill = KnownGeneClass, size = geneIS_frequency_byHitIS), 
#     color = 'white',
#     # size = 3.5, 
#     segment.color = 'black') +
#   geom_label_repel(
#     # data = subset(all_results_cis, tdist_positive_and_correctedEM < 0.05),
#     data = subset(all_results_cis, tdist_positive_and_corrected < 0.05),
#     aes(label = GeneName, fill = KnownGeneClass), 
#     # size = 5,
#     box.padding = unit(0.35, "lines"),
#     point.padding = unit(0.3, "lines"),
#     color = 'white',
#     segment.color = 'black'
#   ) + 
#   geom_label_repel(
#     # data = subset(slice_all_results_cis_consecutive, tdist_positive_and_correctedEM < 0.05),
#     data = subset(all_results_cis, ClinicalRelevance == TRUE ),
#     aes(label = GeneName), 
#     # size = 5,
#     box.padding = unit(0.35, "lines"),
#     point.padding = unit(0.3, "lines"),
#     color = 'orange',
#     segment.color = 'black'
#   ) + 
#   geom_label_repel(
#     # data = subset(all_results_cis_avgByPatient, tdist_positive_and_correctedEM < 0.05),
#     data = subset(all_results_cis, !is.na(KnownClonalExpension) ),
#     aes(label = GeneName), 
#     # size = 5,
#     box.padding = unit(0.35, "lines"),
#     point.padding = unit(0.3, "lines"),
#     color = 'firebrick',
#     segment.color = 'firebrick'
#   ) + 
#   theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0), strip.text.x = element_text(size = 16, colour = "blue", angle = 0)) +
#   theme(strip.text = element_text(face="bold", size=16)) +
#   # theme(strip.background = element_rect(fill="darkblue", colour="white", size=1)) +
#   theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
#   labs(list(title = paste(ProjectID, "- Volcano plot of IS gene frequency and CIS results"), 
#             y = "P-value Grubbs test (-log(p) base 10), Bonferroni correction", 
#             x = "Integration frequency (log2)", 
#             color = "Onco TumSupp Genes", 
#             size = "Avg Transcr. Len", 
#             subtitle = paste0("Significance threshold for annotation labeling: P-value < 0.05 (Bonferroni adjusted; -log = ", (round(-log(0.05, base = 10), 3)), ").\nOnco/TS genes source: UniProt (other genes labeled as 'Other'). Annotated if P-value > ", round(annotation_threshold_ontots_log, 3) )) ) 
# 
# pdf(file = paste(cis_base_folder, analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".AllPatients.CIS.results_invivo.pdf", sep = ""), height=20, width=16)
# print(plot_cis_alldetails)
# dev.off()
# png(file = paste(cis_base_folder, analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".AllPatients.CIS.results_invivo.png", sep = ""), height=20, width=16, units = "in", res = 300)
# print(plot_cis_alldetails)
# dev.off()

# all_results_cis_slice <- all_results_cis[which(all_results_cis$n_IS_perGene > 2 & all_results_cis$average_TxLen > 150 & !(all_results_cis$GeneName %in% c("WAS", "ARSA")) & all_results_cis$SubjectID %in% patients_for_pivotal_report),]
all_results_cis_slice <- all_results_cis[which(all_results_cis$n_IS_perGene > 2 & all_results_cis$average_TxLen > 150 & !(all_results_cis$GeneName %in% c("WAS", "ARSA", "OPTC")) ),]
# all_results_cis_slice$ClinicalTrial <- as.character(substr(as.character(all_results_cis_slice$SubjectID), 1, 3))
# all_results_cis_slice$ClinicalTrial <- sub("BTH", "BTHAL", all_results_cis_slice$ClinicalTrial)
all_results_cis_slice$ClinicalTrial <- all_results_cis_slice$ProjectID

# # all_results_cis_slice <- all_results_cis[which(all_results_cis$n_IS_perGene > 2 & all_results_cis$average_TxLen > 150 & all_results_cis$SubjectID %in% patients_for_pivotal_report),]
# plot_cis_fdr <- ggplot(data = all_results_cis_slice, aes(y = minus_log_p_fdr, x = neg_zscore_minus_log2_int_freq_tolerance, color = KnownGeneClass, fill = KnownGeneClass), na.rm = T, se = TRUE) +
#   # geom_point(aes(size = average_TxLen), alpha = .5) +
#   # geom_hex(bins = 60) +
#   # geom_hex() +
#   geom_point(alpha = .5) +
#   # scale_size(name = "-log(p-value)") +
#   # geom_smooth(method=lm, se=T) +
#   # geom_text(aes(label=ifelse(tophit == TRUE, as.character(GeneName),'')), hjust=0, vjust=0) +
#   # geom_text(aes(label=ifelse(scary == TRUE, as.character(GeneName),'')), hjust=0, vjust=0) + 
#   # geom_density_2d(binwidth = 0.5) +
#   geom_hline(yintercept = significance_threshold_minus_log_p, color='coral', size=1, show.legend = T, linetype="dotted") +
#   scale_y_continuous(limits = c(0, max(c((significance_threshold_minus_log_p + 0.5), max(all_results_cis_slice$minus_log_p_fdr) )) ) ) + 
#   scale_x_continuous(breaks = seq(-4, 4,2)) +
#   # facet_grid( . ~ SubjectID ) +
#   facet_wrap( ~ SubjectID, ncol = 4) +
#   # geom_bin2d(bins = 50) +
#   # geom_label_repel(
#   #   data = subset(all_results_cis, !is.na(Onco1_TS2) & minus_log_p_fdr > annotation_threshold_ontots_log),
#   #   # data = subset(all_results_cis, !is.na(Onco1_TS2) & minus_log2_integration_freq_withtolerance > 2.8),
#   #   aes(label = GeneName, fill = KnownGeneClass), 
#   #   color = 'white',
#   #   # size = 3.5, 
#   #   segment.color = 'black') +
#   geom_label_repel(
#     # data = subset(all_results_cis, tdist_positive_and_correctedEM < 0.05),
#     data = subset(all_results_cis_slice, tdist_fdr < 0.05),
#     aes(label = GeneName), 
#     # size = 5,
#     box.padding = unit(0.35, "lines"),
#     point.padding = unit(0.3, "lines"),
#     color = 'white',
#     segment.color = 'black',
#     max.overlaps = Inf
#   ) + 
#   # geom_label_repel(
#   #   # data = subset(slice_all_results_cis_consecutive, tdist_positive_and_correctedEM < 0.05),
#   #   data = subset(all_results_cis, ClinicalRelevance == TRUE ),
#   #   aes(label = GeneName), 
#   #   # size = 5,
#   #   box.padding = unit(0.35, "lines"),
#   #   point.padding = unit(0.3, "lines"),
#   #   color = 'orange',
#   #   segment.color = 'black'
#   # ) + 
#   # geom_label_repel(
# #   # data = subset(all_results_cis_avgByPatient, tdist_positive_and_correctedEM < 0.05),
# #   data = subset(all_results_cis, !is.na(KnownClonalExpension) ),
# #   aes(label = GeneName), 
# #   # size = 5,
# #   box.padding = unit(0.35, "lines"),
# #   point.padding = unit(0.3, "lines"),
# #   color = 'firebrick',
# #   segment.color = 'firebrick'
# # ) + 
# theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 270), strip.text.x = element_text(size = 16, colour = "blue", angle = 0)) +
#   theme(strip.text = element_text(face="bold", size=16)) +
#   # theme(strip.background = element_rect(fill="darkblue", colour="white", size=1)) +
#   theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
#   labs(title = paste(ProjectID, "- Volcano plot of IS gene frequency and CIS results"), 
#        y = "P-value Grubbs test (-log(p) base 10), FDR correction", 
#        x = "Integration frequency (log2)", 
#        size = "Avg Transcr. Len",
#        color = "Onco TumSupp Genes", 
#        subtitle = paste0("Significance threshold for annotation labeling: P-value < 0.05 (FDR adjusted; -log = ", (round(-log(0.05, base = 10), 3)), ").\nOnco/TS genes source: UniProt (other genes labeled as 'Other'). Annotated if P-value > ", round(annotation_threshold_ontots_log, 3) )) 
# 
# # pdf(file = paste(cis_base_folder, analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".AllPatients.CIS.results_invivo.LiBIS.FDR.pdf", sep = ""), height=16, width=24)
# # print(plot_cis_fdr)
# # dev.off()
# png(file = paste(cis_base_folder, analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".AllPatients.CIS.results.MNC_Whole.FDR.noWAS.png", sep = ""), height=9, width=16, units = "in", res = 300)
# print(plot_cis_fdr)
# dev.off()




# now try to acquire all patients together
# all_results_cis_patient_pvalcorrectedbonferroni <- dcast(data = all_results_cis, GeneName ~ SubjectID, value.var = "tdist_positive_and_corrected", fun.aggregate = mean)
# rownames(all_results_cis_patient_pvalcorrectedbonferroni) <- all_results_cis_patient_pvalcorrectedbonferroni$GeneName
# all_results_cis_patient_pvalcorrectedfdr <- dcast(data = all_results_cis, GeneName ~ SubjectID, value.var = "tdist_fdr", fun.aggregate = mean)
# rownames(all_results_cis_patient_pvalcorrectedfdr) <- all_results_cis_patient_pvalcorrectedfdr$GeneName
# all_results_cis_patient_pvalcorrectedfdr[is.na(all_results_cis_patient_pvalcorrectedfdr)] <- NA
all_results_cis_avgByPatient <- sqldf(paste0("select ", annotation_cols_to_get_asSQLstring, ", ClinicalTrial, count(SubjectID) as ObservedInNPatients, 
                                             avg(tdist_fdr) as avg_tdist_fdr, stdev(tdist_fdr) as stdev_tdist_fdr, max(tdist_fdr) as max_tdist_fdr, min(tdist_fdr) as min_tdist_fdr,
                                             avg(minus_log_p_fdr) as avg_minus_log_p_fdr, stdev(minus_log_p_fdr) as stdev_minus_log_p_fdr, max(minus_log_p_fdr) as max_minus_log_p_fdr,
                                             avg(geneIS_frequency_byHitIS) as avg_geneIS_frequency_byHitIS, stdev(geneIS_frequency_byHitIS) as stdev_geneIS_frequency_byHitIS,
                                             avg(neg_zscore_minus_log2_int_freq_tolerance) as neg_zscore_minus_log2_int_freq_tolerance, avg(neg_zscore_minus_log2_int_freq_tolerance) as avg_neg_zscore_minus_log2_int_freq_tolerance, stdev(neg_zscore_minus_log2_int_freq_tolerance) as stdev_neg_zscore_minus_log2_int_freq_tolerance,
                                             avg(tdist2t) as avg_tdist2t, stdev(tdist2t) as stdev_tdist2t
                                             from all_results_cis_slice 
                                             where average_TxLen > 150
                                             group by ClinicalTrial, GeneName") )

# all_results_cis_avgByPatient <- sqldf(paste0("select ", annotation_cols_to_get_asSQLstring, ", count(SubjectID) as ObservedInNPatients, ClinicalTrial,
#                                           avg(tdist_fdr) as avg_tdist_fdr, stdev(tdist_fdr) as stdev_tdist_fdr,
#                                           avg(minus_log_p_fdr) as avg_minus_log_p_fdr, stdev(minus_log_p_fdr) as stdev_minus_log_p_fdr,
#                                           avg(geneIS_frequency_byHitIS) as avg_geneIS_frequency_byHitIS, stdev(geneIS_frequency_byHitIS) as stdev_geneIS_frequency_byHitIS,
#                                           avg(neg_zscore_minus_log2_integration_freq_withtolerance) as avg_neg_zscore_minus_log2_integration_freq_withtolerance, stdev(neg_zscore_minus_log2_integration_freq_withtolerance) as stdev_neg_zscore_minus_log2_integration_freq_withtolerance,
#                                           avg(tdist2t) as avg_tdist2t, stdev(tdist2t) as stdev_tdist2t
#                                         from all_results_cis_slice 
#                                         where 1 
#                                         group by ClinicalTrial, GeneName") )


all_results_cis_avgByPatient$positive_outlier <- ifelse(all_results_cis_avgByPatient$neg_zscore_minus_log2_int_freq_tolerance > 0, TRUE, FALSE)
# all_results_cis_avgByPatient$positive_outlier <- ifelse(all_results_cis_avgByPatient$avg_neg_zscore_minus_log2_integration_freq_withtolerance > 0, TRUE, FALSE)
# all_results_cis_avgByPatient$positive_outlier_and_significant <- ifelse((all_results_cis_avgByPatient$positive_outlier & (all_results_cis_avgByPatient$avg_tdist_fdr - all_results_cis_avgByPatient$stdev_tdist_fdr) < 0.05), TRUE, FALSE)
# all_results_cis_avgByPatient$positive_outlier_and_significant <- ifelse(is.na(all_results_cis_avgByPatient$stdev_tdist_fdr),
#                                                                         ifelse((all_results_cis_avgByPatient$positive_outlier & all_results_cis_avgByPatient$avg_tdist_fdr < 0.05), TRUE, FALSE),
#                                                                         ifelse((all_results_cis_avgByPatient$positive_outlier & (all_results_cis_avgByPatient$avg_tdist_fdr - all_results_cis_avgByPatient$stdev_tdist_fdr) < 0.05), TRUE, FALSE) )
all_results_cis_avgByPatient$positive_outlier_and_significant <- ifelse( (all_results_cis_avgByPatient$min_tdist_fdr < 0.05), 
                                                                         TRUE, 
                                                                         FALSE)

label_top_n_elements <- 20
all_results_cis_avgByPatient$NTopTargeted <- ifelse(all_results_cis_avgByPatient$avg_neg_zscore_minus_log2_int_freq_tolerance >= sort(all_results_cis_avgByPatient$avg_neg_zscore_minus_log2_int_freq_tolerance, decreasing = T)[label_top_n_elements], TRUE, FALSE)

write.table(x = all_results_cis_avgByPatient, file = paste(cis_base_folder, analysis_folder_date, "/", analysis_folder_date, ".CIS.results_invivo.avg_by_trial.tsv", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = T, na = '')

# all_results_cis_avgByPatient <- read.csv(file = paste(cis_base_folder, analysis_folder_date, "/", analysis_folder_date, ".CIS.results_invivo.avg_by_trial.tsv", sep = ""), header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", ""))


all_results_cis_avgByPatient_slice <- all_results_cis_avgByPatient[which(!(all_results_cis_avgByPatient$GeneName %in% c("WAS", "ARSA", "OPTC"))),]
all_results_cis_avgByPatient_slice$Study <- factor(all_results_cis_avgByPatient_slice$ClinicalTrial, levels = study_list)
write.table(x = all_results_cis_avgByPatient_slice, file = paste(cis_base_folder, analysis_folder_date, "/", analysis_folder_date, ".CIS.results_invivo.avg_by_trial.tsv", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = T, na = '')

plot_cis_avg2 <- ggplot(data = all_results_cis_avgByPatient_slice, aes(y = avg_minus_log_p_fdr, x = avg_neg_zscore_minus_log2_int_freq_tolerance, color = positive_outlier_and_significant), na.rm = T, se = TRUE) +
  geom_point(size = 3, alpha = .5) +
  scale_color_manual(values = c("gray", "violet")) +
  scale_size(name = "Stdev P-value") + 
  facet_grid( . ~ Study) +
  geom_hline(yintercept = significance_threshold_minus_log_p, color='coral', size=1, show.legend = T, linetype="dotted") +
  # scale_x_continuous(breaks = seq(-4, 4,2)) +
  # geom_label_repel(
  #   # data = subset(all_results_cis_avgByPatient, tdist_positive_and_correctedEM < 0.05),
  #   data = subset(all_results_cis_avgByPatient, ClinicalRelevance == TRUE | !is.na(KnownClonalExpension) ),
  #   aes(label = GeneName), 
  #   # size = 5,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines"),
  #   color = 'orange',
  #   segment.color = 'black'
  # ) + 
  geom_label_repel(
    data = subset(all_results_cis_avgByPatient_slice, !is.na(Onco1_TS2) & avg_minus_log_p_fdr > significance_threshold_minus_log_p),
    # data = subset(all_results_cis_avgByPatient, !is.na(Onco1_TS2) & minus_log2_integration_freq_withtolerance > 2.8),
    aes(label = GeneName, fill = KnownGeneClass), 
    color = 'white',
    # size = 3.5, 
    segment.color = 'black') +
  geom_label_repel(
    # data = subset(all_results_cis_avgByPatient, tdist_positive_and_correctedEM < 0.05),
    # data = subset(all_results_cis_avgByPatient, ((avg_tdist_fdr - stdev_tdist_fdr) < 0.05 & positive_outlier_and_significant == TRUE) ),
    data = subset(all_results_cis_avgByPatient_slice, (positive_outlier == TRUE & positive_outlier_and_significant == TRUE) ),
    aes(label = GeneName), 
    # aes(label = GeneName), 
    # size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    color = 'black',
    segment.color = 'black'
  ) + 
  # geom_label_repel(
  #   # data = subset(all_results_cis_avgByPatient, tdist_positive_and_correctedEM < 0.05),
  #   data = subset(all_results_cis_avgByPatient, !is.na(KnownClonalExpension) ),
  #   aes(label = GeneName), 
  #   # size = 5,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines"),
  #   color = 'firebrick',
  #   segment.color = 'firebrick'
  # ) + 
  # geom_label_repel(
  #   # data = subset(all_results_cis_avgByPatient, tdist_positive_and_correctedEM < 0.05),
  #   data = subset(all_results_cis_avgByPatient, NTopTargeted == TRUE ),
  #   aes(label = GeneName, fill = KnownGeneClass), 
  #   # size = 5,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines"),
  #   color = 'white',
  #   segment.color = 'black'
  # ) + 
  # geom_text_repel(data=filter(BTHAL01_bthal01_mergegenename, tophit==TRUE, aes(label=paste(as.character(GeneName), sep = ' '))) ) +
  theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0), strip.text.x = element_text(size = 16, colour = "blue", angle = 0)) +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(strip.background = element_rect(fill="darkblue", colour="white", size=1)) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste("Volcano plot of IS gene frequency and CIS"), 
       y = "P-value Grubbs test (-log10(p)), FDR correction", 
       x = "Integration frequency (log2)", 
       color = "Significant", 
       subtitle = paste0("Average by gene accross patients. Threshold and labeling annotation of genes if P-value < 0.05 (FDR adjusted, -log(p)=", round(significance_threshold_minus_log_p, 3), ").\nOnco/TS genes source: UniProt (other genes labeled as 'Other').")) 


png(file = paste(cis_base_folder, analysis_folder_date, "/", analysis_folder_date, ".CIS.results_invivo.avg_by_trial.figure1.png", sep = ""), height=6, width=13, units = "in", res = 300)
print(plot_cis_avg2)
dev.off()
pdf(file = paste(cis_base_folder, analysis_folder_date, "/", analysis_folder_date, ".CIS.results_invivo.avg_by_trial.figure1.pdf", sep = ""), height=6, width=13)
print(plot_cis_avg2)
dev.off()


plot_cis_avg3 <- ggplot(data = all_results_cis_avgByPatient_slice, aes(y = avg_minus_log_p_fdr, x = avg_neg_zscore_minus_log2_int_freq_tolerance, 
                                                                       color = KnownGeneClass, fill = KnownGeneClass), na.rm = T, se = TRUE) +
  geom_point(size = 3, alpha = .5) +
  # scale_color_manual(values = c("red", "green", "blue")) +
  scale_size(name = "Stdev P-value") + 
  facet_grid( . ~ Study) +
  geom_hline(yintercept = significance_threshold_minus_log_p, color='darkorange', size=1, show.legend = T, linetype="dashed") +
  scale_x_continuous(breaks = seq(-4, 4,2)) +
  geom_label_repel(
    # data = subset(all_results_cis_avgByPatient, tdist_positive_and_correctedEM < 0.05),
    # data = subset(all_results_cis_avgByPatient, ((avg_tdist_fdr - stdev_tdist_fdr) < 0.05 & positive_outlier_and_significant == TRUE) ),
    data = subset(all_results_cis_avgByPatient_slice, (positive_outlier == TRUE & positive_outlier_and_significant == TRUE) ),
    aes(label = GeneName, fill = KnownGeneClass), 
    # aes(label = GeneName), 
    # size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    color = 'white',
    segment.color = 'black'
  ) + 
  theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0), strip.text.x = element_text(size = 16, colour = "blue", angle = 0)) +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(strip.background = element_rect(fill="darkblue", colour="white", size=1)) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste("Volcano plot of IS gene frequency and CIS"), 
       y = "P-value Grubbs test (-log10(p))", 
       x = "Integration frequency (log2)", 
       color = "Gene class", 
       subtitle = paste0("Average by gene accross patients. Threshold and labeling annotation of genes if P-value < 0.05 (FDR adjusted, -log(p)=", round(significance_threshold_minus_log_p, 3), ").\nOnco/TS genes source: UniProt (other genes labeled as 'Other').")) 

png(file = paste(cis_base_folder, analysis_folder_date, "/", analysis_folder_date, ".CIS.results_invivo.avg_by_trial.figure2.png", sep = ""), height=6, width=14, units = "in", res = 300)
print(plot_cis_avg3)
dev.off()
pdf(file = paste(cis_base_folder, analysis_folder_date, "/", analysis_folder_date, ".CIS.results_invivo.avg_by_trial.figure2.pdf", sep = ""), height=6, width=14)
print(plot_cis_avg3)
dev.off()


# prova plot triangolo per geni
# all_results_cis_slice_castfortriplot <- dcast(data = all_results_cis_slice, GeneName ~ ClinicalTrial, value.var = "integration_frequency_withtolerance", fun.aggregate = mean)
all_results_cis_slice_castfortriplot <- dcast(data = all_results_cis_slice, GeneName ~ ClinicalTrial, value.var = "neg_zscore_minus_log2_integration_freq_withtolerance", fun.aggregate = mean)
all_results_cis_slice_castfortriplot$sharing <- apply(all_results_cis_slice_castfortriplot[2:ncol(all_results_cis_slice_castfortriplot)], 1, function(x) {length(x[!is.na(x)])})
genes_observed_in_all_trials <- as.character(all_results_cis_slice_castfortriplot[which(all_results_cis_slice_castfortriplot$sharing == 3), c("GeneName")])
all_results_cis_slice_castfortriplot_allshared <- all_results_cis_slice_castfortriplot[which(all_results_cis_slice_castfortriplot$sharing == 3), c(2:4)]
all_results_cis_slice_castfortriplot_allshared_abs <- abs(all_results_cis_slice_castfortriplot_allshared)

# # non funziona
# ggtern(data=all_results_cis_slice_castfortriplot_allshared_abs, aes(x=MLD, y=WAS, z=BTHAL)) +
#   geom_point() + geom_density_2d(color = "red")

write.table(all_results_cis_slice_castfortriplot_allshared, file = paste(cis_base_folder, analysis_folder_date, "/", analysis_folder_date, ".targetedGenes.zscore.results_invivo.avg_by_trial.csv", sep = ""), sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE, na = '')

# all_results_cis_avgByPatient_slice_castfortriplot <- dcast(data = all_results_cis_avgByPatient_slice, GeneName ~ Study, value.var = "avg_neg_zscore_minus_log2_int_freq_tolerance", fun.aggregate = mean)
# all_results_cis_avgByPatient_slice_castfortriplot$sharing <- apply(all_results_cis_avgByPatient_slice_castfortriplot[2:ncol(all_results_cis_avgByPatient_slice_castfortriplot)], 1, function(x) {length(x[!is.na(x)])})

# try analysis of variance by gene (min 3 patients)
all_results_cis_slice_genecount_onpatients_bytrial <- sqldf(paste0("select ClinicalTrial, GeneName, count(*) as NPatients from all_results_cis_slice where GeneName in ('", paste(genes_observed_in_all_trials, collapse = "','", sep = ''), "') group by ClinicalTrial, GeneName"))
# all_results_cis_slice_genecount_onpatients_bytrial_min3pt <- all_results_cis_slice_genecount_onpatients_bytrial[which(all_results_cis_slice_genecount_onpatients_bytrial$NPatients >= 3), "GeneName"]
all_results_cis_slice_genecount_onpatients_bytrial_cast <- dcast(data = all_results_cis_slice_genecount_onpatients_bytrial, GeneName ~ ClinicalTrial, value.var = "NPatients", fun.aggregate = mean)
all_results_cis_slice_genecount_onpatients_bytrial_cast$minPatients <- apply(all_results_cis_slice_genecount_onpatients_bytrial_cast[2:4], 1, min)
all_results_cis_slice_genecount_onpatients_bytrial_cast_slice <- all_results_cis_slice_genecount_onpatients_bytrial_cast[which(all_results_cis_slice_genecount_onpatients_bytrial_cast$minPatients >= 3 & all_results_cis_slice_genecount_onpatients_bytrial_cast$minPatients <= 9),]
# all_results_cis_slice_genecount_onpatients_bytrial_cast_slice_min4 <- all_results_cis_slice_genecount_onpatients_bytrial_cast[which(all_results_cis_slice_genecount_onpatients_bytrial_cast$minPatients >= 4 & all_results_cis_slice_genecount_onpatients_bytrial_cast$minPatients <= 9),]
# all_results_cis_slice_genecount_onpatients_bytrial_cast_slice_min5 <- all_results_cis_slice_genecount_onpatients_bytrial_cast[which(all_results_cis_slice_genecount_onpatients_bytrial_cast$minPatients >= 5 & all_results_cis_slice_genecount_onpatients_bytrial_cast$minPatients <= 9),]
# all_results_cis_slice_genecount_onpatients_bytrial_cast_slice_min2 <- all_results_cis_slice_genecount_onpatients_bytrial_cast[which(all_results_cis_slice_genecount_onpatients_bytrial_cast$minPatients >= 2 & all_results_cis_slice_genecount_onpatients_bytrial_cast$minPatients <= 9),]
# all_results_cis_slice_genecount_onpatients_bytrial_cast_slice_min8 <- all_results_cis_slice_genecount_onpatients_bytrial_cast[which(all_results_cis_slice_genecount_onpatients_bytrial_cast$minPatients >= 8 & all_results_cis_slice_genecount_onpatients_bytrial_cast$minPatients <= 9),]

# now you can filter by selected genes
# all_results_cis_slice_sharedgenesmin2pt <- all_results_cis_slice[which(all_results_cis_slice$GeneName %in% as.character(all_results_cis_slice_genecount_onpatients_bytrial_cast_slice_min2$GeneName)),]
all_results_cis_slice_sharedgenesmin3pt <- all_results_cis_slice[which(all_results_cis_slice$GeneName %in% as.character(all_results_cis_slice_genecount_onpatients_bytrial_cast_slice$GeneName)),]
# all_results_cis_slice_sharedgenesmin4pt <- all_results_cis_slice[which(all_results_cis_slice$GeneName %in% as.character(all_results_cis_slice_genecount_onpatients_bytrial_cast_slice_min4$GeneName)),]
# all_results_cis_slice_sharedgenesmin5pt <- all_results_cis_slice[which(all_results_cis_slice$GeneName %in% as.character(all_results_cis_slice_genecount_onpatients_bytrial_cast_slice_min5$GeneName)),]
# all_results_cis_slice_sharedgenesmin8pt <- all_results_cis_slice[which(all_results_cis_slice$GeneName %in% as.character(all_results_cis_slice_genecount_onpatients_bytrial_cast_slice_min8$GeneName)),]

# # no min N patients
# kruskal.test(neg_zscore_minus_log2_integration_freq_withtolerance ~ ClinicalTrial, data = all_results_cis_slice)
# pairwise.wilcox.test(all_results_cis_slice_sharedgenesmin3pt$neg_zscore_minus_log2_integration_freq_withtolerance, all_results_cis_slice_sharedgenesmin3pt$ClinicalTrial, p.adjust.method = "BH")
# # min 2
# kruskal.test(neg_zscore_minus_log2_integration_freq_withtolerance ~ ClinicalTrial, data = all_results_cis_slice_sharedgenesmin2pt)
# pairwise.wilcox.test(all_results_cis_slice_sharedgenesmin2pt$neg_zscore_minus_log2_integration_freq_withtolerance, all_results_cis_slice_sharedgenesmin2pt$ClinicalTrial, p.adjust.method = "BH")
# pairwise.wilcox.test(all_results_cis_slice_sharedgenesmin2pt$neg_zscore_minus_log2_integration_freq_withtolerance, all_results_cis_slice_sharedgenesmin2pt$ClinicalTrial, p.adjust.method = "fdr")
# min 3
kruskal.test(neg_zscore_minus_log2_integration_freq_withtolerance ~ ClinicalTrial, data = all_results_cis_slice_sharedgenesmin3pt)
pairwise.wilcox.test(all_results_cis_slice_sharedgenesmin3pt$neg_zscore_minus_log2_integration_freq_withtolerance, all_results_cis_slice_sharedgenesmin3pt$ClinicalTrial, p.adjust.method = "BH")
pairwise.wilcox.test(all_results_cis_slice_sharedgenesmin3pt$neg_zscore_minus_log2_integration_freq_withtolerance, all_results_cis_slice_sharedgenesmin3pt$ClinicalTrial, p.adjust.method = "fdr")
# # min 4
# kruskal.test(neg_zscore_minus_log2_integration_freq_withtolerance ~ ClinicalTrial, data = all_results_cis_slice_sharedgenesmin4pt)
# pairwise.wilcox.test(all_results_cis_slice_sharedgenesmin3pt$neg_zscore_minus_log2_integration_freq_withtolerance, all_results_cis_slice_sharedgenesmin3pt$ClinicalTrial, p.adjust.method = "BH")
# pairwise.wilcox.test(all_results_cis_slice_sharedgenesmin3pt$neg_zscore_minus_log2_integration_freq_withtolerance, all_results_cis_slice_sharedgenesmin3pt$ClinicalTrial, p.adjust.method = "fdr")
# # min 5
# kruskal.test(neg_zscore_minus_log2_integration_freq_withtolerance ~ ClinicalTrial, data = all_results_cis_slice_sharedgenesmin5pt)
# pairwise.wilcox.test(all_results_cis_slice_sharedgenesmin5pt$neg_zscore_minus_log2_integration_freq_withtolerance, all_results_cis_slice_sharedgenesmin3pt$ClinicalTrial, p.adjust.method = "BH")
# pairwise.wilcox.test(all_results_cis_slice_sharedgenesmin3pt$neg_zscore_minus_log2_integration_freq_withtolerance, all_results_cis_slice_sharedgenesmin3pt$ClinicalTrial, p.adjust.method = "fdr")
# # min 8
# # kruskal.test(neg_zscore_minus_log2_integration_freq_withtolerance ~ ClinicalTrial, data = all_results_cis_slice_sharedgenesmin8pt)

# ANOVA
# Compute the analysis of variance
res.aov <- aov(neg_zscore_minus_log2_integration_freq_withtolerance ~ ClinicalTrial, data = all_results_cis_slice_sharedgenesmin3pt)
# Summary of the analysis
summary(res.aov)

# PROVA CON GENI PER TRIAL 


###############################################################
## Gene Ontology
###############################################################
go_base_folder <- "analyses/go/"
### test GO
all_results_cis_avgByPatient_GOlist <- sqldf("select GeneName, max(avg_minus_log_p_fdr) as p from all_results_cis_avgByPatient where positive_outlier_and_significant like 1 group by GeneName order by p desc") # all_results_cis_avgByPatient[which(all_results_cis_avgByPatient$positive_outlier_and_significant), c("avg_minus_log_p_fdr", "GeneName")]
all_results_cis_avgByPatient_GOlist_ext <- sqldf("select GeneName, max(avg_minus_log_p_fdr) as p from all_results_cis_avgByPatient where 1 group by GeneName order by p desc limit 500") # all_results_cis_avgByPatient[which(all_results_cis_avgByPatient$positive_outlier_and_significant), c("avg_minus_log_p_fdr", "GeneName")]
all_results_cis_avgByPatient_GOlist$rank <- seq(1, nrow(all_results_cis_avgByPatient_GOlist))
gene <- as.character(all_results_cis_avgByPatient_GOlist$GeneName)

geneList <- all_results_cis_avgByPatient_GOlist[,2]
names(geneList) <- as.character(all_results_cis_avgByPatient_GOlist[,1])

library(clusterProfiler)
library(org.Hs.eg.db)
# gene.df <- bitr(gene, 
#                 fromType = "ENTREZID",
#                 toType = c("ENSEMBL", "SYMBOL"),
#                 OrgDb = org.Hs.eg.db)
gene.df <- bitr(gene, 
                fromType="SYMBOL", 
                # toType = c("ENTREZID", "ENSEMBL", "SYMBOL"),
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
ego_bp <- enrichGO(gene          = gene.df$ENTREZID,
                # universe      = names(geneList),
                # universe      = gene.df$SYMBOL,
                # keyType       = 'ENSEMBL',
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
ego_cc <- enrichGO(gene          = gene.df$ENTREZID,
                # universe      = names(geneList),
                # universe      = gene.df$SYMBOL,
                # keyType       = 'ENSEMBL',
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
ego_mf <- enrichGO(gene          = gene.df$ENTREZID,
                   # universe      = names(geneList),
                   # universe      = gene.df$SYMBOL,
                   # keyType       = 'ENSEMBL',
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

# ggo <- groupGO(gene     = gene.df$ENTREZID,
#                OrgDb    = org.Hs.eg.db,
#                ont      = "CC",
#                level    = 3,
#                readable = TRUE)
# ego2 <- enrichGO(gene          = gene.df$SYMBOL,
#                 # universe      = names(geneList),
#                 # keyType       = 'ENSEMBL',
#                 keyType       = 'SYMBOL',
#                 OrgDb         = org.Hs.eg.db,
#                 ont           = "CC",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.01,
#                 qvalueCutoff  = 0.05,
#                 readable      = TRUE)
# 
# ego2 <- enrichGO(gene          = gene.df$SYMBOL,
#                 # universe      = names(geneList),
#                 # keyType       = 'ENSEMBL',
#                 keyType       = 'SYMBOL',
#                 OrgDb         = org.Hs.eg.db,
#                 ont           = "CC",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.01,
#                 qvalueCutoff  = 0.05,
#                 readable      = TRUE)

library(msigdbr)
all_gene_sets <- msigdbr(species = "Homo sapiens")


library(DOSE)
library(enrichplot)

barplot(ego_bp, showCategory=50) 
barplot(ego_cc, showCategory=50)
barplot(ego_mf, showCategory=50) 


# # ego_gsedo <- gseDO(geneList = geneList)
# ego_gsedo <- gseDO(geneList = gene.df$ENTREZID)
dotplot(ego_bp, showCategory=30) + ggtitle("dotplot for ORA")
# dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")

## convert gene ID to Symbol
edox_bp <- setReadable(ego_bp, 'org.Hs.eg.db', 'ENTREZID')
p3_bp <- cnetplot(edox_bp, node_label="all") 
edox_cc <- setReadable(ego_cc, 'org.Hs.eg.db', 'ENTREZID')
p3_cc <- cnetplot(edox_cc, node_label="all") 
edox_mf <- setReadable(ego_mf, 'org.Hs.eg.db', 'ENTREZID')
p3_mf <- cnetplot(edox_mf, node_label="all") 

edo_pairs <- pairwise_termsim(ego_bp)
p1 <- emapplot(edo_pairs)
p2 <- emapplot(edo_pairs, cex_category=1.5)
p3 <- emapplot(edo_pairs, layout="kk")
p3 <- emapplot(edo_pairs, layout="kk", cex_category=1.5)

# p1 <- cnetplot(edox, foldChange=geneList)
# ## categorySize can be scaled by 'pvalue' or 'geneNum'
# p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
# p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
# cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

library(GOSemSim)
hsGO_cc <- godata('org.Hs.eg.db', ont="CC")
hsGO_mf <- godata('org.Hs.eg.db', ont="MF")
hsGO_bp <- godata('org.Hs.eg.db', ont="BP")

goSim("GO:0019886", "GO:0045637", semData=hsGO_bp, measure="Jiang")

# go1 <- as.data.frame(ego_bp)$ID
# mgoSim(go1, go2, semData=hsGO, measure="Wang", combine=NULL)
# mgoSim(go1, go2, semData=hsGO, measure="Jiang", combine=NULL)




# mld_GOlist <- sqldf("select GeneName, max(avg_minus_log_p_fdr) as p from all_results_cis_avgByPatient where positive_outlier_and_significant like 1 and ClinicalTrial like 'MLD' group by GeneName order by p desc")
# was_GOlist <- sqldf("select GeneName, max(avg_minus_log_p_fdr) as p from all_results_cis_avgByPatient where positive_outlier_and_significant like 1 and ClinicalTrial like 'WAS' group by GeneName order by p desc")

# mld_GOlist <- sqldf("select GeneName, max(avg_minus_log_p_fdr) as p from all_results_cis_avgByPatient where ClinicalTrial like 'MLD' group by GeneName order by p desc limit 20") 
# was_GOlist <- sqldf("select GeneName, max(avg_minus_log_p_fdr) as p from all_results_cis_avgByPatient where ClinicalTrial like 'WAS' group by GeneName order by p desc limit 20") 

# WAS gamma retro
bed_short_col_names <- c("chr", "integration_locus", "is_end", "hg19_annot_elem_chr", "annotsource", "hg19_annot_elem_type", "hg19_annot_elem_start", "hg19_annot_elem_end", "f1", "hg19_annot_elem_strand", "f2", "hg19_annot_elem_details", "hg19_annot_elem_distance")
bed_wasde <- read.csv(file = "source/WAS_DE/bed/annotated.new.WAS.Germany.P1.bed", header=F, fill=T, check.names = FALSE, sep = '\t')
names(bed_wasde) <- bed_short_col_names

bed_wasde_annot_geneid <- as.data.frame( t(as.data.frame(lapply(strsplit(as.character(bed_wasde$hg19_annot_elem_details), ';', fixed = T), function(x) {c(x)[1]}))) )
bed_wasde_annot_geneid$V1 <- gsub("gene_id ", "", bed_wasde_annot_geneid$V1)
names(bed_wasde_annot_geneid) <- "GeneName"
bed_wasde <- cbind(bed_wasde, bed_wasde_annot_geneid)

# do CIS
bed_wasde_df <- bed_wasde[c("chr", "integration_locus", "hg19_annot_elem_strand", "GeneName", "hg19_annot_elem_strand", "hg19_annot_elem_end")]
names(bed_wasde_df) <- c(id_cols, "WASDE_point")
bed_wasde_df$chr <- gsub("^chr", "", bed_wasde_df$chr)
gammawas_cis_df <- CISGrubbs(df = bed_wasde_df, genomic_annotation_genebased_file = "/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/TIGET/Workspace/Andrea/RefGene/hg19.refGene.oracle.tsv", add_standard_padjust = T, annotation_cols = id_cols)
gammawas_cis_df$minus_log_p_fdr <- -log(gammawas_cis_df$tdist_fdr, base = 10)



mld_GOlist <- sqldf("select GeneName, max(avg_minus_log_p_fdr) as p from all_results_cis_avgByPatient where avg_neg_zscore_minus_log2_int_freq_tolerance >= 0 and ClinicalTrial like 'MLD' group by GeneName order by p desc")
was_GOlist <- sqldf("select GeneName, max(avg_minus_log_p_fdr) as p from all_results_cis_avgByPatient where avg_neg_zscore_minus_log2_int_freq_tolerance >= 0 and ClinicalTrial like 'WAS' group by GeneName order by p desc")
bthal_GOlist <- sqldf("select GeneName, max(avg_minus_log_p_fdr) as p from all_results_cis_avgByPatient where avg_neg_zscore_minus_log2_int_freq_tolerance >= 0 and ClinicalTrial like 'BTHAL' group by GeneName order by p desc")
gammawas_GOlist <- sqldf("select GeneName, max(minus_log_p_fdr) as p from gammawas_cis_df where neg_zscore_minus_log2_integration_freq_withtolerance >= 0 group by GeneName order by p desc")

ontology <- "BP"

mld_gene <- as.character(mld_GOlist$GeneName)
mld_geneList <- mld_GOlist[,2]

names(mld_geneList) <- as.character(mld_GOlist[,1])

mld_gene.df <- bitr(mld_gene, 
                fromType="SYMBOL", 
                # toType = c("ENTREZID", "ENSEMBL", "SYMBOL"),
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
mld_ego_bp <- enrichGO(gene          = mld_gene.df$ENTREZID,
                   # universe      = names(geneList),
                   # universe      = gene.df$SYMBOL,
                   # keyType       = 'ENSEMBL',
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
mld_ego_cc <- enrichGO(gene          = mld_gene.df$ENTREZID,
                       # universe      = names(geneList),
                       # universe      = gene.df$SYMBOL,
                       # keyType       = 'ENSEMBL',
                       OrgDb         = org.Hs.eg.db,
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
mld_ego_mf <- enrichGO(gene          = mld_gene.df$ENTREZID,
                       # universe      = names(geneList),
                       # universe      = gene.df$SYMBOL,
                       # keyType       = 'ENSEMBL',
                       OrgDb         = org.Hs.eg.db,
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

was_gene <- as.character(was_GOlist$GeneName)
was_geneList <- was_GOlist[,2]
names(was_geneList) <- as.character(was_GOlist[,1])

was_gene.df <- bitr(was_gene, 
                    fromType="SYMBOL", 
                    # toType = c("ENTREZID", "ENSEMBL", "SYMBOL"),
                    toType = c("ENTREZID", "SYMBOL"),
                    OrgDb = org.Hs.eg.db)
was_ego_bp <- enrichGO(gene          = was_gene.df$ENTREZID,
                       # universe      = names(geneList),
                       # universe      = gene.df$SYMBOL,
                       # keyType       = 'ENSEMBL',
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",
                       # pAdjustMethod = "fdr",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       # pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
was_ego_cc <- enrichGO(gene          = was_gene.df$ENTREZID,
                       # universe      = names(geneList),
                       # universe      = gene.df$SYMBOL,
                       # keyType       = 'ENSEMBL',
                       OrgDb         = org.Hs.eg.db,
                       ont           = "CC",
                       # pAdjustMethod = "fdr",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       # pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
was_ego_mf <- enrichGO(gene          = was_gene.df$ENTREZID,
                       # universe      = names(geneList),
                       # universe      = gene.df$SYMBOL,
                       # keyType       = 'ENSEMBL',
                       OrgDb         = org.Hs.eg.db,
                       ont           = "MF",
                       # pAdjustMethod = "fdr",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       # pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)


bthal_gene <- as.character(bthal_GOlist$GeneName)
bthal_geneList <- bthal_GOlist[,2]
names(bthal_geneList) <- as.character(bthal_GOlist[,1])

bthal_gene.df <- bitr(bthal_gene, 
                    fromType="SYMBOL", 
                    # toType = c("ENTREZID", "ENSEMBL", "SYMBOL"),
                    toType = c("ENTREZID", "SYMBOL"),
                    OrgDb = org.Hs.eg.db)
bthal_ego_bp <- enrichGO(gene          = bthal_gene.df$ENTREZID,
                       # universe      = names(geneList),
                       # universe      = gene.df$SYMBOL,
                       # keyType       = 'ENSEMBL',
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
bthal_ego_cc <- enrichGO(gene          = bthal_gene.df$ENTREZID,
                       # universe      = names(geneList),
                       # universe      = gene.df$SYMBOL,
                       # keyType       = 'ENSEMBL',
                       OrgDb         = org.Hs.eg.db,
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
bthal_ego_mf <- enrichGO(gene          = bthal_gene.df$ENTREZID,
                       # universe      = names(geneList),
                       # universe      = gene.df$SYMBOL,
                       # keyType       = 'ENSEMBL',
                       OrgDb         = org.Hs.eg.db,
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

gammawas_gene <- as.character(gammawas_GOlist$GeneName)
gammawas_geneList <- gammawas_GOlist[,2]
names(gammawas_geneList) <- as.character(gammawas_GOlist[,1])

gammawas_gene.df <- bitr(gammawas_gene, 
                    fromType="SYMBOL", 
                    # toType = c("ENTREZID", "ENSEMBL", "SYMBOL"),
                    toType = c("ENTREZID", "SYMBOL"),
                    OrgDb = org.Hs.eg.db)
gammawas_ego_bp <- enrichGO(gene          = gammawas_gene.df$ENTREZID,
                       # universe      = names(geneList),
                       # universe      = gene.df$SYMBOL,
                       # keyType       = 'ENSEMBL',
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
gammawas_ego_cc <- enrichGO(gene          = gammawas_gene.df$ENTREZID,
                       # universe      = names(geneList),
                       # universe      = gene.df$SYMBOL,
                       # keyType       = 'ENSEMBL',
                       OrgDb         = org.Hs.eg.db,
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
gammawas_ego_mf <- enrichGO(gene          = gammawas_gene.df$ENTREZID,
                       # universe      = names(geneList),
                       # universe      = gene.df$SYMBOL,
                       # keyType       = 'ENSEMBL',
                       OrgDb         = org.Hs.eg.db,
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

library(GOSemSim)
mld_go_bp <- as.data.frame(mld_ego_bp)$ID
was_go_bp <- as.data.frame(was_ego_bp)$ID
bthal_go_bp <- as.data.frame(bthal_ego_bp)$ID
gammawas_ego_bp <- as.data.frame(gammawas_ego_bp)$ID
mw_mgoSim_compare_wang <- mgoSim(mld_go_bp, was_go_bp, semData=hsGO_bp, measure="Wang", combine=NULL)
mw_mgoSim_compare_jiang <- mgoSim(mld_go_bp, was_go_bp, semData=hsGO_bp, measure="Jiang", combine=NULL)

mb_mgoSim_compare_wang <- mgoSim(mld_go_bp, bthal_go_bp, semData=hsGO_bp, measure="Wang", combine=NULL)
mb_mgoSim_compare_jiang <- mgoSim(mld_go_bp, bthal_go_bp, semData=hsGO_bp, measure="Jiang", combine=NULL)

mgoSim_compare_BMA <- mgoSim(mld_go_bp, was_go_bp, semData=hsGO_bp, measure="Wang", combine="BMA")
# mgoSim_compare_BMA <- mw_mgoSim_compare_BMA
mw_mgoSim_compare_BMA_bp <- mgoSim(mld_go_bp, was_go_bp, semData=hsGO_bp, measure="Jiang", combine="BMA")
mb_mgoSim_compare_BMA_bp <- mgoSim(mld_go_bp, bthal_go_bp, semData=hsGO_bp, measure="Jiang", combine="BMA")
wb_mgoSim_compare_BMA_bp <- mgoSim(was_go_bp, bthal_go_bp, semData=hsGO_bp, measure="Jiang", combine="BMA")
gammaw_mgoSim_compare_BMA_bp <- mgoSim(gammawas_ego_bp, was_go_bp, semData=hsGO_bp, measure="Jiang", combine="BMA")
gammaw_mgoSim_compare_bp <- mgoSim(gammawas_ego_bp, was_go_bp, semData=hsGO_bp, measure="Jiang", combine=NULL)

# CC GO
mld_go_cc <- as.data.frame(mld_ego_cc)$ID
was_go_cc <- as.data.frame(was_ego_cc)$ID
bthal_go_cc <- as.data.frame(bthal_ego_cc)$ID
gammawas_go_cc <- as.data.frame(gammawas_ego_cc)$ID
mw_mgoSim_compare_BMA_cc <- mgoSim(mld_go_cc, was_go_cc, semData=hsGO_cc, measure="Jiang", combine="BMA")
mb_mgoSim_compare_BMA_cc <- mgoSim(mld_go_cc, bthal_go_cc, semData=hsGO_cc, measure="Jiang", combine="BMA")
wb_mgoSim_compare_BMA_cc <- mgoSim(was_go_cc, bthal_go_cc, semData=hsGO_cc, measure="Jiang", combine="BMA")
gammaw_mgoSim_compare_BMA_cc <- mgoSim(gammawas_go_cc, was_go_cc, semData=hsGO_cc, measure="Jiang", combine="BMA")
# MF GO
mld_go_mf <- as.data.frame(mld_ego_mf)$ID
was_go_mf <- as.data.frame(was_ego_mf)$ID
bthal_go_mf <- as.data.frame(bthal_ego_mf)$ID
gammawas_go_mf <- as.data.frame(gammawas_ego_mf)$ID
mw_mgoSim_compare_BMA_mf <- mgoSim(mld_go_mf, was_go_mf, semData=hsGO_mf, measure="Jiang", combine="BMA")
mb_mgoSim_compare_BMA_mf <- mgoSim(mld_go_mf, bthal_go_mf, semData=hsGO_mf, measure="Jiang", combine="BMA")
wb_mgoSim_compare_BMA_mf <- mgoSim(was_go_mf, bthal_go_mf, semData=hsGO_mf, measure="Jiang", combine="BMA")
gammaw_mgoSim_compare_BMA_mf <- mgoSim(was_go_mf, gammawas_go_mf, semData=hsGO_mf, measure="Jiang", combine="BMA")


# build up a full df
bp_allresults <- c(mw_mgoSim_compare_BMA_bp, mb_mgoSim_compare_BMA_bp, wb_mgoSim_compare_BMA_bp)
cc_allresults <- c(mw_mgoSim_compare_BMA_cc, mb_mgoSim_compare_BMA_cc, wb_mgoSim_compare_BMA_cc)
mf_allresults <- c(mw_mgoSim_compare_BMA_mf, mb_mgoSim_compare_BMA_mf, wb_mgoSim_compare_BMA_mf)

all_results_go <- data.frame(
  "BP" = bp_allresults,
  "CC" = cc_allresults,
  "MF" = mf_allresults
)
all_results_go <- data.frame(t(all_results_go))
names(all_results_go) <- c("MLD-WAS", "MLD-BTHAL", "BTHAL-WAS")
write.xlsx(x = all_results_go, file = paste0("analyses/go/", analysis_folder_date, "/", analysis_folder_date, ".GO.AllCategories.AmongTrials.xlsx"), sheetName = "Jiang BMA")
# 
# #### ternary
# myData <- all_results_go
# 
# TernaryPlot(
#   atip = "",
#   btip = "",
#   ctip = "",
#   alab = "Quartz",
#   blab = "Feldspar",
#   clab = "Rock.fragments",
#   lab.offset = 0.16,
#   lab.col = "#000000",
#   point = "up",
#   clockwise = TRUE,
#   xlim = NULL,
#   ylim = NULL,
#   lab.cex = 1,
#   lab.font = 1,
#   tip.cex = 1,
#   tip.font = 1,
#   tip.col = "#000000",
#   isometric = TRUE,
#   atip.rotate = NULL,
#   btip.rotate = NULL,
#   ctip.rotate = NULL,
#   atip.pos = NULL,
#   btip.pos = NULL,
#   ctip.pos = NULL,
#   padding = 0.08,
#   col = "#FFFFFF",
#   grid.lines = 10,
#   grid.col = "#A9A9A9",
#   grid.lty = "solid",
#   grid.lwd = 1,
#   grid.minor.lines = 4,
#   grid.minor.col = "#D3D3D3",
#   grid.minor.lty = "solid",
#   grid.minor.lwd = 1,
#   axis.lty = "solid",
#   axis.labels = TRUE,
#   axis.cex = 0.8,
#   axis.font = 1,
#   axis.rotate = TRUE,
#   axis.tick = TRUE,
#   axis.lwd = 1,
#   ticks.lwd = 1,
#   ticks.length = 0.025,
#   axis.col = "#000000",
#   ticks.col = "#A9A9A9"
# )
# 
# TernaryPoints(all_results_go[, 1:3],
#               type = "p",
#               cex = 1.8,
#               pch = 16,
#               lwd = 1,
#               lty = "solid",
#               col = "#222222"
# )
# 
# TernaryPlot(point = all_results_go, atip = 'A', btip = 'B', ctip = 'C',
#             alab = 'Aness', blab = 'Bness', clab = 'Cness')
# TernaryText(list(A = c(10, 1, 1), B = c(1, 10, 1), C = c(1, 1, 10)),
#             col = cbPalette8[4], font = 2)

library(ade4)
all_results_go <- read.csv(file = "/Users/calabria.andrea/Dropbox (HSR Global)/Project Clinical Trial Comparison/analyses/go/202104/202104.GO.AllCategories.AmongTrials.tsv", header=T, fill=T, sep='\t', check.names = F, row.names = 1)
# triangle.plot(all_results_go, label = row.names(all_results_go), clab = 1, labeltriangle = T)
wtriangleplot <- triangle.plot(all_results_go, label = row.names(all_results_go), 
              labeltriangle = T, cpoint = 3, 
              max3 = c(0.34, 0.34, 0.34), min3 = c(0.33, 0.33, 0.33), 
              clabel = 2, draw.line = T, 
              sub = "GO similarity", csub = 2, possub = "topright", box = F)
# points(wtriangleplot, col = "blue", pch = 20, cex = 5)
png(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.GOSimilarity.Triangularplot.figure1.png", sep = ""), height=8, width=8, units = "in", res = 300)
triangle.plot(all_results_go, label = row.names(all_results_go), 
              labeltriangle = T, cpoint = 3, 
              max3 = c(0.34, 0.34, 0.34), min3 = c(0.33, 0.33, 0.33), 
              clabel = 2, draw.line = T, 
              sub = "GO similarity", csub = 2, possub = "topright", box = F)
dev.off()
pdf(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.GOSimilarity.Triangularplot.figure1.pdf", sep = ""), height=8, width=8)
triangle.plot(all_results_go, label = row.names(all_results_go), 
              labeltriangle = T, cpoint = 3, 
              max3 = c(0.34, 0.34, 0.34), min3 = c(0.33, 0.33, 0.33), 
              clabel = 2, draw.line = T, 
              sub = "GO similarity", csub = 2, possub = "topright", box = F)
dev.off()




mgoSim_compare_matrix <- as.matrix(mgoSim_compare_wang)
mgoSim_compare_matrix <- as.matrix(mgoSim_compare_jiang)
# mgoSim_compare_matrix[mgoSim_compare_matrix==1] <- NA
pmap_go <- pheatmap(mgoSim_compare_matrix, scale = "none", 
                    main = paste0("GO similarity scores", ontology, " (overall score: ", mgoSim_compare_BMA, ")")
                    # cutree_rows = 2,
                    # cutree_cols = 2
                    )

gamma_pmap_go <- pheatmap(as.matrix(gammaw_mgoSim_compare_bp), scale = "none", 
                    main = paste0("GO similarity scores", ontology, " (overall score: ", gammaw_mgoSim_compare_BMA_bp, ")")
                    # cutree_rows = 2,
                    # cutree_cols = 2
                    )

# png(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, "SemSim_Jiang.figure1.png", sep = ""), height=8, width=8, units = "in", res = 300)
pheatmap(mgoSim_compare_matrix, scale = "none", 
                    main = paste0("GO similarity scores ", ontology, " (overall score: ", mgoSim_compare_BMA, ")"), 
                    filename = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, "SemSim_Jiang.figure2.png", sep = ""), 
                    show_rownames = F, show_colnames = F
                    # cutree_rows = 4,
                    # cutree_cols = 4
)
# dev.off()

pheatmap(mgoSim_compare_matrix, scale = "none", 
         main = paste0("GO similarity scores ", ontology, " (overall score: ", mgoSim_compare_BMA, ")"), 
         filename = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, "SemSim_Jiang.figure2.png", sep = ""), 
         show_rownames = F, show_colnames = F
         # cutree_rows = 4,
         # cutree_cols = 4
)
pheatmap(mgoSim_compare_matrix, scale = "none", 
         main = paste0("GO similarity scores ", ontology, " (overall score: ", mgoSim_compare_BMA, ")"), 
         filename = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, "SemSim_Jiang.figure2.pdf", sep = ""), 
         show_rownames = F, show_colnames = F
         # cutree_rows = 4,
         # cutree_cols = 4
)

pheatmap(as.matrix(gammaw_mgoSim_compare_bp), scale = "none", 
         main = paste0("GO similarity scores ", ontology, " (overall score: ", gammaw_mgoSim_compare_BMA_bp, ")"), 
         filename = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".WAS-Gamma_WAS-LV.SemSim_Jiang.figure2.pdf", sep = ""),
         show_rownames = F, show_colnames = F
         # cutree_rows = 4,
         # cutree_cols = 4
)
pheatmap(as.matrix(gammaw_mgoSim_compare_bp), scale = "none", 
         main = paste0("GO similarity scores ", ontology, " (overall score: ", gammaw_mgoSim_compare_BMA_bp, ")"), 
         filename = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".WAS-Gamma_WAS-LV.SemSim_Jiang.figure2.png", sep = ""),
         show_rownames = F, show_colnames = F
         # cutree_rows = 4,
         # cutree_cols = 4
)


venn(list("mld" = mld_go_bp, "was" = was_go_bp))

hsGO2 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="BP", computeIC=FALSE) 

# barplot(mld_ego_bp, showCategory=20) 
png(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".MLD.barplot.figure1.png", sep = ""), height=12, width=9, units = "in", res = 300)
barplot(mld_ego_bp, showCategory=50, title = paste("MLD -", ontology) ) 
dev.off()
pdf(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".MLD.barplot.figure1.pdf", sep = ""), height=12, width=9)
barplot(mld_ego_bp, showCategory=50, title = paste("MLD -", ontology)) 
dev.off()
# dotplot(ego_bp, showCategory=30) + ggtitle("dotplot for ORA")
# dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")

png(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".WAS.barplot.figure1.png", sep = ""), height=12, width=11.5, units = "in", res = 300)
barplot(was_ego_bp, showCategory=50, title = paste("WAS -", ontology) ) 
dev.off()
pdf(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".WAS.barplot.figure1.pdf", sep = ""), height=12, width=12)
barplot(was_ego_bp, showCategory=50, title = paste("WAS -", ontology) ) 
dev.off()

# barplot(bthal_ego_bp, showCategory=50) 
png(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".BTHAL.barplot.figure1.png", sep = ""), height=12, width=9, units = "in", res = 300)
barplot(bthal_ego_bp, showCategory=50, title = paste("BTHAL -", ontology) ) 
dev.off()
pdf(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".BTHAL.barplot.figure1.pdf", sep = ""), height=12, width=9)
barplot(bthal_ego_bp, showCategory=50, title = paste("BTHAL -", ontology)) 
dev.off()

### ----- now use only 15 classes ------
png(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".MLD.barplot.figure2.png", sep = ""), height=6, width=7, units = "in", res = 300)
barplot(mld_ego_bp, showCategory=15, title = paste("MLD -", ontology) ) 
dev.off()
pdf(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".MLD.barplot.figure2.pdf", sep = ""), height=6, width=7)
barplot(mld_ego_bp, showCategory=15, title = paste("MLD -", ontology)) 
dev.off()

png(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".WAS.barplot.figure2.png", sep = ""), height=6, width=7, units = "in", res = 300)
barplot(was_ego_bp, showCategory=15, title = paste("WAS -", ontology) ) 
dev.off()
pdf(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".WAS.barplot.figure2.pdf", sep = ""), height=6, width=7)
barplot(was_ego_bp, showCategory=15, title = paste("WAS -", ontology) ) 
dev.off()

png(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".BTHAL.barplot.figure2.png", sep = ""), height=6, width=7, units = "in", res = 300)
barplot(bthal_ego_bp, showCategory=15, title = paste("BTHAL -", ontology) ) 
dev.off()
pdf(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".BTHAL.barplot.figure2.pdf", sep = ""), height=6, width=7)
barplot(bthal_ego_bp, showCategory=15, title = paste("BTHAL -", ontology)) 
dev.off()


## convert gene ID to Symbol
mld_edox_bp <- setReadable(mld_ego_bp, 'org.Hs.eg.db', 'ENTREZID')
mld_p3_bp <- cnetplot(mld_edox_bp, node_label="all") 
png(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".MLD.cnetplot.figure1.png", sep = ""), height=8, width=8, units = "in", res = 300)
cnetplot(mld_edox_bp, node_label="all") 
dev.off()

was_edox_bp <- setReadable(was_ego_bp, 'org.Hs.eg.db', 'ENTREZID')
was_p3_bp <- cnetplot(was_edox_bp, node_label="all") 
png(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".WAS.cnetplot.figure1.png", sep = ""), height=8, width=8, units = "in", res = 300)
cnetplot(was_edox_bp, node_label="all") 
dev.off()


# edox_cc <- setReadable(ego_cc, 'org.Hs.eg.db', 'ENTREZID')
# p3_cc <- cnetplot(edox_cc, node_label="all") 
# edox_mf <- setReadable(ego_mf, 'org.Hs.eg.db', 'ENTREZID')
# p3_mf <- cnetplot(edox_mf, node_label="all") 

mld_edo_pairs <- pairwise_termsim(mld_ego_bp)
was_edo_pairs <- pairwise_termsim(was_ego_bp)
p1 <- emapplot(edo_pairs)
p2 <- emapplot(edo_pairs, cex_category=1.5)
p3 <- emapplot(edo_pairs, layout="kk")
p3 <- emapplot(edo_pairs, layout="kk", cex_category=1.5)

mld_pairs <- emapplot(mld_edo_pairs, layout="kk", cex_category=1.5)

png(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".MLD.emapplot.figure1.png", sep = ""), height=8, width=8, units = "in", res = 300)
emapplot(mld_edo_pairs, layout="kk", cex_category=1.5)
dev.off()
png(file = paste(go_base_folder, analysis_folder_date, "/", analysis_folder_date, ".GO.", ontology, ".WAS.emapplot.figure1.png", sep = ""), height=8, width=8, units = "in", res = 300)
emapplot(was_edo_pairs, layout="kk", cex_category=1.5)
dev.off()


# similarity among common genes
union_genes <- union(x = mld_gene, y = was_gene)
genesim_matrix <- mgeneSim(union_genes, semData=hsGO2, measure="Wang", combine="BMA", verbose=FALSE)
genesim_matrix <- mgeneSim(union_genes, semData=hsGO2, measure="Wang")




# # similarity between gen lists
# clusterSim(cluster1 = mld_gene, cluster2 = was_gene, semData=hsGO_bp, measure="Wang", combine="BMA")
# clusterSim(cluster1 = mld_gene, cluster2 = was_gene, semData=hsGO_bp, measure="Wang", combine=NULL)

library(pheatmap)
genesim_matrix_nodiagonal <- as.matrix(genesim_matrix)
genesim_matrix_nodiagonal[genesim_matrix_nodiagonal==1] <- NA
pmap <- pheatmap(genesim_matrix_nodiagonal, scale = "none", 
         cutree_rows = 4,
         cutree_cols = 4)


###############################################################
## Waves of Long lasting clones
###############################################################
markerlist <- c("CD13", "CD14", "CD15", "CD19", "CD3", "CD34", "GLY", "CD36", "GLYA")
wave_base_folder <- paste0("analyses/waves/")
dir.create(file.path(getwd(), wave_base_folder), showWarnings = FALSE)
dir.create(file.path(wave_base_folder, analysis_folder_date), showWarnings = FALSE)
min_is_for_sharing <- 50 # numver of min IS at a specific time point for considering the sample in the sharing

mld_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial MLD/"
bthal_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial bThalassemia/"
was_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial WAS/"
names_to_replace <- "AggregationGroup|SuperGroup"

# import all data
mld_wave <- read.csv(file = paste0(mld_wd, "analyses/15.waves/", analysis_folder_date, "/20220121_MLD_source_of_iss_SuperGroup_purityfilt.tsv.gz"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
mld_wave$ClinicalTrial <- "MLD"
names(mld_wave) <- gsub(names_to_replace, "CellMarker", colnames(mld_wave))
# colname_to_convertto_cellmarker <- "SuperGroup"
# if (colname_to_convertto_cellmarker %in% colnames(mld_wave) & !("CellMarker" %in% colnames(mld_wave))) {
#   mld_wave$CellMarker <- mld_wave[,colname_to_convertto_cellmarker]
# }

was_wave <- read.csv(file = paste0(was_wd, "analyses/15.waves/", analysis_folder_date, "/20220121_WAS_source_of_iss_SuperGroup_purityfilt_LAM-only.tsv.gz"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
was_wave$ClinicalTrial <- "WAS"
names(was_wave) <- gsub(names_to_replace, "CellMarker", colnames(was_wave))
# colname_to_convertto_cellmarker <- "SuperGroup"
# if (colname_to_convertto_cellmarker %in% colnames(was_wave) & !("CellMarker" %in% colnames(was_wave))) {
#   was_wave$CellMarker <- was_wave[,colname_to_convertto_cellmarker]
# }

bthal_wave <- read.csv(file = paste0(bthal_wd, "analyses/15.waves/", analysis_folder_date, "/20220121_BTHAL_source_of_iss_AggregationGroup_purityfilt.tsv.gz"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
bthal_wave$ClinicalTrial <- "BTHAL"
names(bthal_wave) <- gsub(names_to_replace, "CellMarker", colnames(bthal_wave))
# colname_to_convertto_cellmarker <- "AggregationGroup"
# if (colname_to_convertto_cellmarker %in% colnames(bthal_wave) & !("CellMarker" %in% colnames(bthal_wave))) {
#   bthal_wave$CellMarker <- bthal_wave[,colname_to_convertto_cellmarker]
# }

# merge
common_labels <- Reduce(intersect, list(colnames(mld_wave), colnames(was_wave), colnames(bthal_wave)))
all_results_wave <- Reduce(rbind, list(mld_wave[common_labels], was_wave[common_labels], bthal_wave[common_labels]))

all_results_wave$Study <- factor(all_results_wave$ClinicalTrial, levels = study_list)
all_results_wave[is.na(all_results_wave)] <- NA

# add info of lineages
library(openxlsx)
blood_lineages <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "blood_lineages")
# blood_lineages_colors <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "colors_noery")
blood_lineages_colors <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "colors")
blood_lineages <- merge(x = blood_lineages, y = blood_lineages_colors, all.x = T, by = c("CellType"))
rownames(blood_lineages) <- blood_lineages$CellMarker
all_results_wave$CellMarker <- all_results_wave$g1_CellMarker
all_results_wave <- merge(x = all_results_wave, y = blood_lineages, by = c("CellMarker"), all.x = T)
all_results_wave$CellLineage <- factor(all_results_wave$CellType, levels = celltype_list)
# rownames(all_results_wave) <- as.character(all_results_wave$gdf_names)

scale_color_manual_colors_sortedbyname <- blood_lineages_colors[order(blood_lineages_colors$CellType), "colorcode"]

# slice data
all_results_wave$SourceTime <- str_pad(all_results_wave$g1_TimepointMonths, 2, pad = "0")
# all_results_wave$g1 <- apply(all_results_wave[c("g1_SubjectID", "g1_CellMarker", "g1_Tissue", "g1_TimepointMonths")], 1, function(x) {paste(x[1], x[2], x[3], str_pad(round(as.numeric(x[4])), 2, pad = "0"), sep = '_')})
all_results_wave_slice <- all_results_wave[which(all_results_wave$shared > 0 &
                                                 (all_results_wave$g1 != all_results_wave$g2) &
                                                 (all_results_wave$count_g1 >= min_is_for_sharing & all_results_wave$count_g2 >= min_is_for_sharing) 
                                                 ),]

all_results_wave_slice_max <- all_results_wave_slice %>% group_by(g1_SubjectID, Study) %>% summarise(max_g2 = max(on_g2), g2_TimepointMonths = max(g2_TimepointMonths))

write.xlsx(x = all_results_wave_slice, file = paste(wave_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WavesIS.xlsx", sep = ""), sheetName = "All patients")

# do plots
plot_wave_barplot_perc <- ggplot(data = all_results_wave_slice[which(all_results_wave_slice$g1_SubjectID == "MLD14" 
                                                                     # & all_results_wave_slice$g1_Tissue == "PB"
                                                                     )
                                                               ,],  
                                 aes(x = g2_TimepointMonths, y = on_g2, colour = SourceTime, fill = SourceTime, alluvium = SourceTime)) + 
  # geom_bar() +
  # geom_bar(stat="identity", color="white", width = 2, position = position_stack(reverse = TRUE)) +
  geom_bar(stat="identity", color="white", width = 2, position = position_stack(reverse = F)) +
  # geom_alluvium(alpha = .75, decreasing = FALSE) +
  # geom_stratum() +
  # geom_text(aes(label=scales::percent(reads_percent, accuracy = 0.1) ), vjust=-.5, color="black") +
  # geom_text(aes(label=shared ), color="black", vjust=-.5, angle = 0, position = position_stack(vjust = 0.5)) +
  geom_text(aes(label=shared ), color="black", position = position_stack(vjust = 0.5)) +
  # geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray")
  theme_bw() +
  facet_wrap(g2_Tissue ~ g1_CellMarker, ncol = 3, scales = "free_y") +
  # scale_y_continuous(labels=scales::percent) +
  scale_x_continuous(breaks = seq(0, max(all_results_wave_slice$g2_TimepointMonths, na.rm = T), 6) ) +
  theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=14, angle = 0), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=22)) +
  labs(title = paste0("Analysis of waves of IS"), 
       x = "Months after GT", y = "Perc. IS", colour = "Source time", fill = "Source time", 
       subtitle = paste0("First observed clone (IS) is colored with the corresponding time point.\nNumbers within each bar represent the absolute number of shared IS of that time point.") ) 

pdf(file = paste(wave_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WavesIS.MLD14.pdf", sep = ""), height=11, width=14)
plot(plot_wave_barplot_perc)
dev.off()
png(file = paste(wave_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WavesIS.MLD14.png", sep = ""), height=11, width=14, units = "in", res = 300)
plot(plot_wave_barplot_perc)
dev.off()


plot_wave_splineplot_perc <- ggplot(data = all_results_wave_slice[which(all_results_wave_slice$g1_SubjectID == "WAS1001" & all_results_wave_slice$CellType == "Myeloid"),],  
                                    aes(x = g2_TimepointMonths, y = on_g2, colour = SourceTime, fill = SourceTime)) + 
  geom_point(alpha = 0.2, size = 2) +
  # geom_jitter(alpha = 0.8, size = 4) +
  geom_line(size=3, alpha = .7) +
  # geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  # geom_smooth(method = "loess", stat = "smooth", alpha = 0.6, level = 0.75) +
  # geom_text(aes(label=shared ), color="black", position = position_stack(vjust = 0.5)) +
  theme_bw() +
  facet_wrap(g2_Tissue ~ g1_CellMarker, ncol = 3, scales = "free_y") +
  # geom_label_repel(data = subset(all_results_wave_slice_max, max_g2>=10),
  #                  aes(label = g1_SubjectID),
  #                  color = 'white',
  #                  box.padding = unit(0.35, "lines"),
  #                  point.padding = unit(0.4, "lines"),
  #                  segment.color = 'black') +
  # scale_y_continuous(labels=scales::percent) +
  scale_x_continuous(breaks = seq(0, max(all_results_wave_slice$g2_TimepointMonths, na.rm = T), 6) ) +
  theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=14, angle = 0), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=22)) +
  labs(title = paste0("Analysis of waves of IS"), 
       x = "Months after GT", y = "Perc. IS", colour = "Source time", fill = "Source time", 
       subtitle = paste0("First observed IS will be colored with the corresponding time point.") ) 

# pdf(file = paste(cum_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".CumulativeNIS.ByPatient.pdf", sep = ""), height=5.5, width=13)
# plot(plot_cumulative_by_patient)
# dev.off()
# png(file = paste(cum_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".CumulativeNIS.ByPatient.png", sep = ""), height=5.5, width=13, units = "in", res = 300)
# plot(plot_cumulative_by_patient)
# dev.off()


# analysis, steps
# 1. get only late clones and counts
all_results_wave_slice_late <- all_results_wave_slice[which(all_results_wave_slice$g2_TimepointMonths >= 36),]
# 2. add phases to time points: 1-4:early, 9-18:mid, 24-30:steady
# cut(all_results_wave_slice_late$g1_TimepointMonths, breaks=c(0, 6.1, 9, Inf), labels=c("lo", "mid", "hi"))
all_results_wave_slice_late$SourceTime_Class <- case_when((all_results_wave_slice_late$g1_TimepointMonths > 0) & (all_results_wave_slice_late$g1_TimepointMonths < 5) ~ "early",
          (all_results_wave_slice_late$g1_TimepointMonths >= 9) & (all_results_wave_slice_late$g1_TimepointMonths <= 18) ~ "mid",
          (all_results_wave_slice_late$g1_TimepointMonths >= 24) & (all_results_wave_slice_late$g1_TimepointMonths <= 31) ~ "steady")
# 3. compute stats (avg, sd) by time class. qui fai media tra tutte le stacked bars >36 per il gruppo scelto (occhio che nel qui scegli per marker, dopo farai stats per classe di tempo; l'alternativa è mettere le medie direttamente sulla classe di tempo che livella le misurazioni ripetute)
all_results_wave_slice_late <- all_results_wave_slice_late %>% group_by(g1_SubjectID, g1_CellMarker, g1_Tissue, g1_TimepointMonths) %>% mutate(inLate_avg_on_g2 = mean(on_g2, na.rm = T), inLate_sd_on_g2 = sd(on_g2, na.rm = T))
write.xlsx(x = all_results_wave_slice_late, file = paste(wave_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WavesIS.statsFocus_LongLastingClones.xlsx", sep = ""), sheetName = "All patients")

# 4. do stats (first descriptive: boxplots)

plot_wave_slice_late_boxplot_test <- ggplot(data = all_results_wave_slice_late[which(all_results_wave_slice_late$g1_SubjectID == "MLD11" 
                                                                                # & all_results_wave_slice_late$CellType == "Myeloid"
                                                                                & !is.na(all_results_wave_slice_late$SourceTime_Class) ),],  
                                       aes(x = SourceTime_Class, y = on_g2, colour = SourceTime, fill = SourceTime), na.rm = T, se = TRUE) +
  # scale_color_manual(values = trials_colors) +
  # scale_fill_manual(values = trials_colors) +
  # geom_point(alpha = 0.2, size = 2) +
  geom_boxplot(alpha = .7, outlier.size = 0) +
  geom_jitter(alpha = 0.8, size = 4) +
  # geom_line(size=3, alpha = .7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  # geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  # scale_x_continuous(breaks = seq(0, max(all_results_h_slice$TimePoint, na.rm = T), 6) ) +
  facet_wrap(g2_Tissue ~ g1_CellMarker, ncol = 3, scales = "free_y") +
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
  labs(title = paste0("Source of last clones, backtracking in time"), 
       subtitle = "Filters: selected in timepoints >=36. Months split as follows: 1-4 early, 9-18 mid, 24-30 steady.\nGiven last time points, get avg of source time by frame and plot values.",
       x = "Source time points", 
       y = "Sharing perc. (average)", 
       colour = "Study", fill = "Study")

plot_wave_slice_late_boxplot_original <- ggplot(data = all_results_wave_slice_late[which(all_results_wave_slice_late$g1_SubjectID == "MLD14" 
                                                                                # & all_results_wave_slice_late$CellType == "Myeloid"
                                                                                & !is.na(all_results_wave_slice_late$SourceTime_Class) ),],  
                                       aes(x = SourceTime_Class, y = inLate_avg_on_g2, colour = SourceTime, fill = SourceTime, group = SourceTime_Class ), na.rm = T, se = TRUE) +
  # scale_color_manual(values = trials_colors) +
  # scale_fill_manual(values = trials_colors) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  geom_point(alpha = 0.6, size = 4) +
  # geom_jitter(alpha = 0.8, size = 4) +
  # geom_line(size=3, alpha = .7) +
  facet_wrap(g2_Tissue ~ g1_CellMarker, ncol = 3, scales = "free_y") +
  theme_bw() +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Source of last clones, backtracking in time"), 
       subtitle = "Filters: selected in timepoints >=36. Months split as follows: 1-4 early, 9-18 mid, 24-30 steady.\nGiven last time points, get avg of source time by frame and plot values.",
       x = "Source time points", 
       y = "Sharing perc. (average)", 
       colour = "Source time", fill = "Source time")

plot_wave_slice_late_boxplot <- ggplot(data = all_results_wave_slice_late[which(all_results_wave_slice_late$g1_SubjectID == "MLD14" 
                                                                                # & all_results_wave_slice_late$CellType == "Myeloid"
                                                                                & !is.na(all_results_wave_slice_late$SourceTime_Class) ),],  
                                       aes(x = SourceTime_Class, y = inLate_avg_on_g2, group = SourceTime_Class ), na.rm = T, se = TRUE) +
  # scale_color_manual(values = trials_colors) +
  # scale_fill_manual(values = trials_colors) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  geom_point(alpha = 0.6, size = 4) +
  # geom_jitter(alpha = 0.8, size = 4) +
  # geom_line(size=3, alpha = .7) +
  facet_wrap(g2_Tissue ~ g1_CellMarker, ncol = 3, scales = "free_y") +
  theme_bw() +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Source of last clones, backtracking in time"), 
       subtitle = "Filters: selected in timepoints >=36. Months split as follows: 1-4 early, 9-18 mid, 24-30 steady.\nGiven last time points, get avg of source time by frame and plot values.",
       x = "Source time points", 
       y = "Sharing perc. (average)", 
       colour = "Source time", fill = "Source time")

pdf(file = paste(wave_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WavesIS.MLD14.lateClones_sourceAvg.v2.pdf", sep = ""), height=11, width=14)
plot(plot_wave_slice_late_boxplot)
dev.off()
png(file = paste(wave_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WavesIS.MLD14.lateClones_sourceAvg.v2.png", sep = ""), height=11, width=14, units = "in", res = 300)
plot(plot_wave_slice_late_boxplot)
dev.off()

plot_wave_slice_late_boxplot_bystudy <- ggplot(data = all_results_wave_slice_late[which(!is.na(all_results_wave_slice_late$SourceTime_Class) 
                                                                                        & all_results_wave_slice_late$Study == "MLD"
                                                                                        ),],  
                                       aes(x = SourceTime_Class, y = inLate_avg_on_g2, colour = SourceTime, fill = SourceTime, group = SourceTime_Class ), na.rm = T, se = TRUE) +
  # scale_color_manual(values = trials_colors) +
  # scale_fill_manual(values = trials_colors) +
  # geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  geom_violin(colour = "darkblue", fill = "lightblue", alpha = 0.6) +
  # geom_line(size=3, alpha = .7) +
  facet_wrap(g2_Tissue ~ g1_CellMarker, ncol = 3, scales = "free_y") +
  theme_bw() +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Source of last clones, backtracking in time"), 
       subtitle = "Filters: selected in timepoints >=36. Months split as follows: 1-4 early, 9-18 mid, 24-30 steady.\nGiven last time points, get avg of source time by frame and plot values.",
       x = "Source time points", 
       y = "Sharing perc. (average)", 
       colour = "Source time", fill = "Source time")

# dataset for stats and plot
all_results_wave_slice_late_nodup <- unique(all_results_wave_slice_late[which(!is.na(all_results_wave_slice_late$SourceTime_Class) ), c("ClinicalTrial", "Study", "g1_SubjectID", "CellType", "Keywords", "HematoLineage", "SuperGroup", "AggregationGroup", "LineageByPurity", "HematoLineageMutations", "colorcode", "CellLineage", "SourceTime", "g2_Tissue", "g1_CellMarker", "inLate_avg_on_g2", "SourceTime_Class") ])
write.xlsx(x = all_results_wave_slice_late_nodup, file = paste(wave_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WavesIS.statsFocus_LongLastingClones.unique.xlsx", sep = ""), sheetName = "All patients")


my_comparisons <- list(
  c("early", "mid"),
  c("early", "steady"),
  c("mid", "steady")
)
plot_wave_slice_late_boxplot_bystudy_stats <- ggplot(data = all_results_wave_slice_late[which(!is.na(all_results_wave_slice_late$SourceTime_Class) 
                                                                                        & all_results_wave_slice_late$Study == "MLD"
                                                                                        #& all_results_wave_slice_late$g1_SubjectID == "MLD11"
                                                                                        ),],  
                                                    # aes(x = SourceTime_Class, y = inLate_avg_on_g2, colour = SourceTime, fill = SourceTime, group = SourceTime_Class ), na.rm = T, se = TRUE) +
                                                    aes(x = SourceTime_Class, y = inLate_avg_on_g2, group = SourceTime_Class ), na.rm = T, se = TRUE) +
  # scale_color_manual(values = trials_colors) +
  # scale_fill_manual(values = trials_colors) +
  # geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.4, size = 2) +
  geom_violin(colour = "darkblue", fill = "lightblue", alpha = 0.6) +
  # geom_line(size=3, alpha = .7) +
  facet_wrap(g2_Tissue ~ g1_CellMarker, ncol = 3) +
  theme_bw() +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(30, 35, 40), p.adjust.method = "fdr", label = "p.signif") +
  stat_compare_means(comparisons = my_comparisons, label.y = c(70, 60, 50)*1.1, label = "p.signif") +
  scale_y_continuous(limits = c(0,100)) +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Source of last clones, backtracking in time"), 
       subtitle = "Filters: selected in timepoints >=36. Months split as follows: 1-4 early, 9-18 mid, 24-30 steady.\nGiven last time points, get avg of source time by frame and plot values.",
       x = "Source time points", 
       y = "Sharing perc. (average)", 
       colour = "Source time", fill = "Source time")

pdf(file = paste(wave_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WavesIS.MLD.lateClones_sourceAvg.Stats.v2.pdf", sep = ""), height=14, width=14)
plot(plot_wave_slice_late_boxplot_bystudy_stats)
dev.off()
png(file = paste(wave_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WavesIS.MLD.lateClones_sourceAvg.Stats.v2.png", sep = ""), height=14, width=14, units = "in", res = 300)
plot(plot_wave_slice_late_boxplot_bystudy_stats)
dev.off()


plot_wave_slice_late_boxplot_bystudy_celltype_stats <- ggplot(data = all_results_wave_slice_late_nodup,  
                                                              # aes(x = SourceTime_Class, y = inLate_avg_on_g2, colour = SourceTime, fill = SourceTime, group = SourceTime_Class ), na.rm = T, se = TRUE) +
                                                              aes(x = SourceTime_Class, y = inLate_avg_on_g2, group = SourceTime_Class ), na.rm = T, se = TRUE) +
  # scale_color_manual(values = trials_colors) +
  # scale_fill_manual(values = trials_colors) +
  # geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.4, size = 2, color = "lightblue") +
  geom_violin(colour = "darkblue", fill = "steelblue2", alpha = 0.9) +
  # geom_line(size=3, alpha = .7) +
  # facet_wrap(. ~ CellType, ncol = 6) +
  # facet_grid(Study ~ CellLineage, scales = "free_y") +
  facet_grid(Study ~ CellLineage) +
  theme_bw() +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(30, 35, 40), p.adjust.method = "fdr", label = "p.signif") +
  stat_compare_means(comparisons = my_comparisons, label.y = c(50, 40, 30), label = "p.signif") +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Source of last clones, backtracking in time"), 
       subtitle = "Filters: selected in timepoints >=36. Months split as follows: 1-4 early, 9-18 mid, 24-30 steady.\nGiven last time points, get avg of source time by frame and plot values.",
       x = "Source time points", 
       y = "Sharing perc. (average)", 
       colour = "Source time", fill = "Source time")

pdf(file = paste(wave_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WavesIS.", paste(study_list, collapse = "_") ,".lateClones_sourceAvg.Stats.byLineage.v2.pdf", sep = ""), height=8, width=14)
plot(plot_wave_slice_late_boxplot_bystudy_celltype_stats)
dev.off()
png(file = paste(wave_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WavesIS.", paste(study_list, collapse = "_"), ".lateClones_sourceAvg.Stats.byLineage.v2.png", sep = ""), height=8, width=12, units = "in", res = 300)
plot(plot_wave_slice_late_boxplot_bystudy_celltype_stats)
dev.off()


###############################################################
## Sharing IS with WHOLE populations
###############################################################
markerlist <- c("CD13", "CD14", "CD15", "CD19", "CD3", "CD34", "GLY", "CD36", "GLYA")
wholesharing_base_folder <- paste0("analyses/sharing/")
dir.create(file.path(getwd(), wholesharing_base_folder), showWarnings = FALSE)
dir.create(file.path(wholesharing_base_folder, analysis_folder_date), showWarnings = FALSE)
min_is_for_sharing <- 50 # numver of min IS at a specific time point for considering the sample in the sharing

mld_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial MLD/"
bthal_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial bThalassemia/"
was_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial WAS/"
names_to_replace <- "AggregationGroup|SuperGroup"

# import all data
mld_whsharing <- read.csv(file = paste0(mld_wd, "analyses/12.sharing/", analysis_folder_date, "/20220121_MLD_WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_SuperGroup.tsv.gz"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
mld_whsharing$ClinicalTrial <- "MLD"
names(mld_whsharing) <- gsub(names_to_replace, "CellMarker", colnames(mld_whsharing))
# colname_to_convertto_cellmarker <- "SuperGroup"
# if (colname_to_convertto_cellmarker %in% colnames(mld_whsharing) & !("CellMarker" %in% colnames(mld_whsharing))) {
#   mld_whsharing$CellMarker <- mld_whsharing[,colname_to_convertto_cellmarker]
# }

was_whsharing <- read.csv(file = paste0(was_wd, "analyses/13.sharing/", analysis_folder_date, "/20220121_WAS_WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_SuperGroup_LAM-only.tsv.gz"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
was_whsharing$ClinicalTrial <- "WAS"
names(was_whsharing) <- gsub(names_to_replace, "CellMarker", colnames(was_whsharing))
# colname_to_convertto_cellmarker <- "SuperGroup"
# if (colname_to_convertto_cellmarker %in% colnames(was_whsharing) & !("CellMarker" %in% colnames(was_whsharing))) {
#   was_whsharing$CellMarker <- was_whsharing[,colname_to_convertto_cellmarker]
# }

bthal_whsharing <- read.csv(file = paste0(bthal_wd, "analyses/13.sharing/", analysis_folder_date, "/20220121_BTHAL_WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_AggregationGroup.tsv.gz"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
bthal_whsharing$ClinicalTrial <- "BTHAL"
names(bthal_whsharing) <- gsub(names_to_replace, "CellMarker", colnames(bthal_whsharing))
# colname_to_convertto_cellmarker <- "AggregationGroup"
# if (colname_to_convertto_cellmarker %in% colnames(bthal_whsharing) & !("CellMarker" %in% colnames(bthal_whsharing))) {
#   bthal_whsharing$CellMarker <- bthal_whsharing[,colname_to_convertto_cellmarker]
# }

# merge
common_labels <- Reduce(intersect, list(colnames(mld_whsharing), colnames(was_whsharing), colnames(bthal_whsharing)))
all_results_whsharing <- Reduce(rbind, list(mld_whsharing[common_labels], was_whsharing[common_labels], bthal_whsharing[common_labels]))

all_results_whsharing$Study <- factor(all_results_whsharing$ClinicalTrial, levels = study_list)
all_results_whsharing[is.na(all_results_whsharing)] <- NA

colnames_aftersplit <- c("SubjectID", "CellMarker", "Tissue")
all_results_whsharing_details_g1 <- as.data.frame( t(as.data.frame(lapply(strsplit(as.character(all_results_whsharing$g1), '_', fixed = T), function(x) {c(x)}))) )
all_results_whsharing_details_g2 <- as.data.frame( t(as.data.frame(lapply(strsplit(as.character(all_results_whsharing$g2), '_', fixed = T), function(x) {c(x)}))) )
names(all_results_whsharing_details_g1) <- paste0("g1_", c("SubjectID", "Tissue", "CellMarker"))
names(all_results_whsharing_details_g2) <- paste0("g2_", colnames_aftersplit)
all_results_whsharing <- Reduce(cbind, list(all_results_whsharing, all_results_whsharing_details_g1, all_results_whsharing_details_g2))
all_results_whsharing$CellMarker <- all_results_whsharing$g2_CellMarker

# add info of lineages
library(openxlsx)
blood_lineages <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "blood_lineages")
# blood_lineages_colors <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "colors_noery")
blood_lineages_colors <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "colors")
blood_lineages <- merge(x = blood_lineages, y = blood_lineages_colors, all.x = T, by = c("CellType"))
rownames(blood_lineages) <- blood_lineages$CellMarker
all_results_whsharing <- merge(x = all_results_whsharing, y = blood_lineages, by = c("CellMarker"), all.x = T)
all_results_whsharing$CellLineage <- factor(all_results_whsharing$CellType, levels = celltype_list)
# rownames(all_results_whsharing) <- as.character(all_results_whsharing$gdf_names)

scale_color_manual_colors_sortedbyname <- blood_lineages_colors[order(blood_lineages_colors$CellType), "colorcode"]

# slice data
all_results_whsharing_slice <- all_results_whsharing[which(all_results_whsharing$shared > 0 &
                                                           (all_results_whsharing$g1 != all_results_whsharing$g2) &
                                                           (all_results_whsharing$count_g1 >= min_is_for_sharing & all_results_whsharing$count_g2 >= min_is_for_sharing) 
                                                    ),]
all_results_whsharing_slice$SubjectID <- all_results_whsharing_slice$g1_SubjectID
# transform percentage in log odds
all_results_whsharing_slice$logit_perc_on_g1 <- log( ((all_results_whsharing_slice$on_g1)/100) / (1 - (all_results_whsharing_slice$on_g1)/100 ))
all_results_whsharing_slice$logit_perc_on_g2 <- log( ((all_results_whsharing_slice$on_g2)/100) / (1 - (all_results_whsharing_slice$on_g2)/100 ))
# do scaling by patient, tissue
all_results_whsharing_slice <- all_results_whsharing_slice %>% group_by(SubjectID, g1_Tissue, g2_Tissue) %>% mutate(on_g1_zscaled=scale(on_g1), on_g2_zscaled=scale(on_g2))

# all_results_whsharing_slice_max <- all_results_whsharing_slice %>% group_by(g1_SubjectID, Study) %>% summarise(max_g2 = max(on_g2), g2_TimepointMonths = max(g2_TimepointMonths))

write.xlsx(x = all_results_whsharing_slice, file = paste(wholesharing_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WholeSharing.CellMarkers.PB-BM.xlsx", sep = ""), sheetName = "All patients")

# remove some samples here
all_results_whsharing_slice_test <- all_results_whsharing_slice[which( (all_results_whsharing_slice$g2_Tissue == "PB") | 
                                                                         (all_results_whsharing_slice$g2_Tissue == "BM" & all_results_whsharing_slice$CellType == "Erythroid") ) ,]
all_results_whsharing_slice_test <- all_results_whsharing_slice_test[which( (all_results_whsharing_slice_test$g2_Tissue == all_results_whsharing_slice_test$g1_Tissue) ) ,]

# do plots
my_comparisons <- list(
  c("MLD", "WAS"),
  c("MLD", "BTHAL"),
  c("WAS", "BTHAL")
)
plot_whsharing_slice_boxplot_bymarker_distincttissues_stats <- ggplot(data = all_results_whsharing_slice_test,
                                                     aes(x = Study, y = on_g1, colour = CellMarker, fill = CellMarker, group = Study ), na.rm = T, se = TRUE) +
  # scale_color_manual(values = trials_colors) +
  # scale_fill_manual(values = trials_colors) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  # geom_violin(colour = "darkblue", fill = "lightblue", alpha = 0.6) +
  # geom_line(size=3, alpha = .7) +
  facet_wrap(g1_Tissue ~ CellLineage, ncol = 4) +
  theme_bw() +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(30, 35, 40), p.adjust.method = "fdr", label = "p.signif") +
  stat_compare_means(comparisons = my_comparisons, label.y = c(0.55, 0.5, 0.45)*100, label = "p.signif") +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Percentage of shared clones in whole populations."), 
       subtitle = "Filters: purity, PB and BM separated, markers separated.",
       x = "Studies", 
       y = "Sharing percentage", 
       colour = "CellMarker", fill = "CellMarker")

pdf(file = paste(wholesharing_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WholeSharing.CellMarkers.PB-BM.Stats.pdf", sep = ""), height=5, width=12)
plot(plot_whsharing_slice_boxplot_bymarker_distincttissues_stats)
dev.off()
png(file = paste(wholesharing_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WholeSharing.CellMarkers.PB-BM.Stats.png", sep = ""), height=5, width=12, units = "in", res = 300)
plot(plot_whsharing_slice_boxplot_bymarker_distincttissues_stats)
dev.off()

perc_spacer_stats <- -3
perc_starting_max_stats <- 25
plot_whsharing_slice_boxplot_bymarker_distincttissues_stats_scaled <- ggplot(data = all_results_whsharing_slice,
                                                                      aes(x = Study, y = on_g1_zscaled, colour = CellMarker, fill = CellMarker, group = Study, shape = g2_Tissue ), na.rm = T, se = TRUE) +
  # scale_color_manual(values = trials_colors) +
  # scale_fill_manual(values = trials_colors) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  # geom_violin(colour = "darkblue", fill = "lightblue", alpha = 0.6) +
  # geom_line(size=3, alpha = .7) +
  # facet_wrap(g1_Tissue ~ CellLineage, ncol = 4) +
  scale_y_continuous(limits = c(-2, 3)) + 
  facet_grid(g1_Tissue ~ CellLineage) +
  theme_bw() +
  # stat_compare_means(label.y = 6) +
  stat_compare_means(comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Percentage of shared clones in whole populations."), 
       subtitle = "Filters: purity, PB and BM separated, markers separated. Percentages rescaled by patient and tissue (Z-score).",
       x = "Studies", 
       y = "Sharing percentage (Z-Score rescaled)", 
       colour = "Cell Marker", fill = "Cell Marker", shape = "Lineage Tissue")

pdf(file = paste(wholesharing_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WholeSharing.CellMarkers.PB-BM.Stats.zscore.full.pdf", sep = ""), height=8, width=14)
plot(plot_whsharing_slice_boxplot_bymarker_distincttissues_stats_scaled)
dev.off()
png(file = paste(wholesharing_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WholeSharing.CellMarkers.PB-BM.Stats.zscore.full.png", sep = ""), height=8, width=14, units = "in", res = 300)
plot(plot_whsharing_slice_boxplot_bymarker_distincttissues_stats_scaled)
dev.off()
# plot_whsharing_slice_boxplot_bymarker_distincttissues_stats_test <- ggplot(data = all_results_whsharing_slice_test,
#                                                                       aes(x = Study, y = on_g1, colour = CellMarker, fill = CellMarker, group = Study ), na.rm = T, se = TRUE) +
#   # scale_color_manual(values = trials_colors) +
#   # scale_fill_manual(values = trials_colors) +
#   geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
#   # geom_point(alpha = 0.6, size = 4) +
#   geom_jitter(alpha = 0.6, size = 4) +
#   # geom_violin(colour = "darkblue", fill = "lightblue", alpha = 0.6) +
#   # geom_line(size=3, alpha = .7) +
#   facet_wrap(g1_Tissue ~ CellLineage, ncol = 4) +
#   theme_bw() +
#   # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
#   # stat_compare_means(label.y = 6) +
#   # stat_compare_means(comparisons = my_comparisons, label.y = c(30, 35, 40), p.adjust.method = "fdr", label = "p.signif") +
#   stat_compare_means(comparisons = my_comparisons, label.y = c(0.55, 0.5, 0.45)*100, label = "p.signif") +
#   theme(strip.text = element_text(face="bold", size=16)) +
#   # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
#   theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
#   # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
#   labs(title = paste0("Percentage of shared clones in whole populations."), 
#        subtitle = "Filters: purity, PB and BM separated, markers separated.",
#        x = "Studies", 
#        y = "Sharing percentage", 
#        colour = "CellMarker", fill = "CellMarker")

my_comparisons <- list(
  c("Erythroid", "Myeloid"),
  c("Erythroid", "B"),
  c("Erythroid", "T"),
  c("Myeloid", "B"),
  c("Myeloid", "T"),
  c("B", "T")
)
perc_spacer_stats <- -5
plot_whsharing_slice_boxplot_bymarker_distincttissues_stats_lineage <- ggplot(data = all_results_whsharing_slice_test,
                                                                      aes(x = CellLineage, y = on_g1, colour = CellMarker, fill = CellMarker, group = CellLineage ), na.rm = T, se = TRUE) +
  # scale_color_manual(values = trials_colors) +
  # scale_fill_manual(values = trials_colors) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  # geom_violin(colour = "darkblue", fill = "lightblue", alpha = 0.6) +
  # geom_line(size=3, alpha = .7) +
  facet_wrap(. ~ Study, ncol = 4) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 85), breaks = seq(from = 0, to = 100, by = 20)) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(30, 35, 40), p.adjust.method = "fdr", label = "p.signif") +
  # stat_compare_means(comparisons = my_comparisons, label.y = seq(from = 75, to = (75+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*1, label = "p.signif") +
  stat_compare_means(comparisons = my_comparisons, label.y = seq(from = 75, to = (75+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*1, p.adjust.method = "fdr", label = "p.signif") +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Percentage of shared clones in whole populations."), 
       subtitle = "Filters: purity, PB and BM separated, markers separated.",
       x = "Lineage", 
       y = "Sharing percentage", 
       colour = "Cell Marker", fill = "Cell Marker")

pdf(file = paste(wholesharing_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WholeSharing.CellMarkers.PB-BM.Stats.byStudy.pdf", sep = ""), height=5, width=12)
plot(plot_whsharing_slice_boxplot_bymarker_distincttissues_stats_lineage)
dev.off()
png(file = paste(wholesharing_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WholeSharing.CellMarkers.PB-BM.Stats.byStudy.png", sep = ""), height=5, width=12, units = "in", res = 300)
plot(plot_whsharing_slice_boxplot_bymarker_distincttissues_stats_lineage)
dev.off()


# ----------------- now try with all tissues together ---------------------
# import all data
# mld_whsharing <- read.csv(file = paste0(mld_wd, "analyses/12.sharing/", analysis_folder_date, "/20220121_MLD_nested-CD34|WholeMNC-MyErBT-tpCollapsed-PB-BM_SuperGroup.tsv.gz"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
mld_whsharing <- read.csv(file = paste0(mld_wd, "analyses/12.sharing/", analysis_folder_date, "/20211202_MLD_WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_CellTypeOnly.tsv.gz"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
mld_whsharing$ClinicalTrial <- "MLD"
names(mld_whsharing) <- gsub(names_to_replace, "CellMarker", colnames(mld_whsharing))
# colname_to_convertto_cellmarker <- "SuperGroup"
# if (colname_to_convertto_cellmarker %in% colnames(mld_whsharing) & !("CellMarker" %in% colnames(mld_whsharing))) {
#   mld_whsharing$CellMarker <- mld_whsharing[,colname_to_convertto_cellmarker]
# }

# was_whsharing <- read.csv(file = paste0(was_wd, "analyses/13.sharing/", analysis_folder_date, "/20220121_WAS_nested-CD34|WholeMNC-MyErBT-tpCollapsed-PB-BM_SuperGroup.tsv.gz"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
# was_whsharing <- read.csv(file = paste0(was_wd, "analyses/13.sharing/", analysis_folder_date, "/20220121_WAS_WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_SuperGroup_LAM-only.tsv.gz"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
was_whsharing <- read.csv(file = paste0(was_wd, "analyses/13.sharing/", analysis_folder_date, "/20211210_WAS_WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_CellTypeOnly.tsv.gz"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
was_whsharing$ClinicalTrial <- "WAS"
names(was_whsharing) <- gsub(names_to_replace, "CellMarker", colnames(was_whsharing))
# colname_to_convertto_cellmarker <- "SuperGroup"
# if (colname_to_convertto_cellmarker %in% colnames(was_whsharing) & !("CellMarker" %in% colnames(was_whsharing))) {
#   was_whsharing$CellMarker <- was_whsharing[,colname_to_convertto_cellmarker]
# }

# bthal_whsharing <- read.csv(file = paste0(bthal_wd, "analyses/13.sharing/", analysis_folder_date, "/20220121_BTHAL_nested-CD34|WholeMNC-MyErBT-tpCollapsed-PB-BM_AggregationGroup.tsv.gz"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
bthal_whsharing <- read.csv(file = paste0(bthal_wd, "analyses/13.sharing/", analysis_folder_date, "/20220121_BTHAL_WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_CellTypeOnly.tsv.gz"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
bthal_whsharing$ClinicalTrial <- "BTHAL"
names(bthal_whsharing) <- gsub(names_to_replace, "CellMarker", colnames(bthal_whsharing))
# colname_to_convertto_cellmarker <- "AggregationGroup"
# if (colname_to_convertto_cellmarker %in% colnames(bthal_whsharing) & !("CellMarker" %in% colnames(bthal_whsharing))) {
#   bthal_whsharing$CellMarker <- bthal_whsharing[,colname_to_convertto_cellmarker]
# }

# merge
common_labels <- Reduce(intersect, list(colnames(mld_whsharing), colnames(was_whsharing), colnames(bthal_whsharing)))
all_results_whsharing <- Reduce(rbind, list(mld_whsharing[common_labels], was_whsharing[common_labels], bthal_whsharing[common_labels]))

all_results_whsharing$Study <- factor(all_results_whsharing$ClinicalTrial, levels = study_list)
all_results_whsharing[is.na(all_results_whsharing)] <- NA

colnames_aftersplit <- c("SubjectID", "CellMarker", "Tissue")
all_results_whsharing_details_g1 <- as.data.frame( t(as.data.frame(lapply(strsplit(as.character(all_results_whsharing$g1), '_', fixed = T), function(x) {c(x)}))) )
all_results_whsharing_details_g2 <- as.data.frame( t(as.data.frame(lapply(strsplit(as.character(all_results_whsharing$g2), '_', fixed = T), function(x) {c(x)}))) )
# names(all_results_whsharing_details_g1) <- paste0("g1_", c("SubjectID", "CD34", "Whole", "CD34_Tissue", "Whole_Tissue"))
names(all_results_whsharing_details_g1) <- paste0("g1_", c("SubjectID", "Tissue", "CellMarker"))
names(all_results_whsharing_details_g2) <- paste0("g2_", colnames_aftersplit)
all_results_whsharing <- Reduce(cbind, list(all_results_whsharing, all_results_whsharing_details_g1, all_results_whsharing_details_g2))
# all_results_whsharing$CellMarker <- all_results_whsharing$g2_CellMarker

# add info of lineages
blood_lineages <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "blood_lineages")
blood_lineages_colors <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "colors")
blood_lineages <- merge(x = blood_lineages, y = blood_lineages_colors, all.x = T, by = c("CellType"))
rownames(blood_lineages) <- blood_lineages$CellMarker
all_results_whsharing <- merge(x = all_results_whsharing, y = blood_lineages, by = c("CellMarker"), all.x = T)
# all_results_whsharing$CellLineage <- factor(all_results_whsharing$CellType, levels = celltype_list)
all_results_whsharing$CellLineage <- factor(all_results_whsharing$g2_CellMarker, levels = celltype_list)
# rownames(all_results_whsharing) <- as.character(all_results_whsharing$gdf_names)

scale_color_manual_colors_sortedbyname <- blood_lineages_colors[order(blood_lineages_colors$CellType), "colorcode"]

# slice data
all_results_whsharing_slice <- all_results_whsharing[which(all_results_whsharing$shared > 0 &
                                                             (all_results_whsharing$g1 != all_results_whsharing$g2) &
                                                             (all_results_whsharing$count_g1 >= min_is_for_sharing & all_results_whsharing$count_g2 >= min_is_for_sharing) 
),]
# transform percentage in log odds
all_results_whsharing_slice$logit_perc_on_g1 <- log( ((all_results_whsharing_slice$on_g1)/100) / (1 - (all_results_whsharing_slice$on_g1)/100 ))
all_results_whsharing_slice$logit_perc_on_g2 <- log( ((all_results_whsharing_slice$on_g2)/100) / (1 - (all_results_whsharing_slice$on_g2)/100 ))

# all_results_whsharing_slice_max <- all_results_whsharing_slice %>% group_by(g1_SubjectID, Study) %>% summarise(max_g2 = max(on_g2), g2_TimepointMonths = max(g2_TimepointMonths))

write.xlsx(x = all_results_whsharing_slice, file = paste(wholesharing_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WholeSharing_NestedCD34.CellMarkers.tissuesCollapsed.xlsx", sep = ""), sheetName = "All patients")

# remove some samples here
all_results_whsharing_slice_test <- all_results_whsharing_slice[which( (all_results_whsharing_slice$g2_Tissue == "PB") | 
                                                                         (all_results_whsharing_slice$g2_Tissue == "BM" & all_results_whsharing_slice$CellLineage == "Erythroid") ) ,]

# do plots
my_comparisons <- list(
  c("MLD", "WAS"),
  c("MLD", "BTHAL"),
  c("WAS", "BTHAL")
)
plot_whsharing_slice_boxplot_bymarker_collapsedtissues_stats <- ggplot(data = all_results_whsharing_slice_test[which(!(all_results_whsharing_slice_test$g1_SubjectID %in% c("WAS1009"))),],
                                                                      aes(x = Study, y = on_g1, colour = CellLineage, fill = CellLineage, group = Study ), na.rm = T, se = TRUE) +
  # scale_color_manual(values = trials_colors) +
  # scale_fill_manual(values = trials_colors) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  # geom_violin(colour = "darkblue", fill = "lightblue", alpha = 0.6) +
  # geom_line(size=3, alpha = .7) +
  # facet_wrap(g1_Whole_Tissue ~ CellLineage, ncol = 4) +
  facet_grid(g1_Tissue ~ CellLineage) +
  # facet_wrap(g1_Tissue ~ CellLineage, ncol = 4) +
  theme_bw() +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(30, 35, 40), p.adjust.method = "fdr", label = "p.signif") +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(0.55, 0.5, 0.45)*100, label = "p.signif", p.adjust.method = "fdr") +
  stat_compare_means(comparisons = my_comparisons, label.y = c(0.55, 0.5, 0.45)*100, label = "p.signif") +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Percentage of shared clones in whole populations."), 
       subtitle = "Filters: purity, PB and BM separated, markers separated.",
       x = "Studies", 
       y = "Sharing percentage", 
       colour = "CellLineage", fill = "CellLineage")

pdf(file = paste(wholesharing_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WholeSharing.CellMarkers.PB-BM.Stats.pdf", sep = ""), height=5, width=12)
plot(plot_whsharing_slice_boxplot_bymarker_distincttissues_stats)
dev.off()
png(file = paste(wholesharing_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".WholeSharing.CellMarkers.PB-BM.Stats.png", sep = ""), height=5, width=12, units = "in", res = 300)
plot(plot_whsharing_slice_boxplot_bymarker_distincttissues_stats)
dev.off()



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

# write results in a file
write.xlsx(x = all_results_h_slice_scaled, file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.xlsx", sep = ""))
write.xlsx(x = all_results_h_slice_ngdna_scaled, file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.ngdna_notNA.xlsx", sep = ""))

plot_hindex <-
  ggplot(data = all_results_h_slice, aes(x = TimePoint, y = Hindex, fill = CellType, color = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(alpha = 0.2, size = 2) +
  # geom_jitter(alpha = 0.8, size = 4) +
  # geom_line(size=3, alpha = .7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  scale_x_continuous(breaks = seq(0, max(all_results_h_slice_ngdna$TimePoint, na.rm = T), 6) ) +
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
  labs(title = paste0("Diversity index (Shannon)"), 
            x = "Months after gene therapy", 
            y = "Shannon diversity index (H')", 
            colour = "CellMarker", fill = "CellMarker")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.pdf", sep = ""), height=8, width=12)
plot(plot_hindex)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.png", sep = ""), height=8, width=12, units = "in", res = 300)
plot(plot_hindex)
dev.off()

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


plot_jindex <-
  # ggplot(data = all_results_h_slice_ngdna_scaled[which(all_results_h_slice_ngdna_scaled$PCRMethod == "SLiM"),], aes(x = TimePoint, y = Jindex, fill = CellType, color = CellType), na.rm = T, se = TRUE) +
  ggplot(data = all_results_h_slice_ngdna_scaled, aes(x = TimePoint, y = Jindex, fill = CellType, color = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(alpha = 0.2, size = 2) +
  # geom_jitter(alpha = 0.8, size = 4) +
  # geom_line(size=3, alpha = .7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  scale_x_continuous(breaks = seq(0, max(all_results_h_slice$TimePoint, na.rm = T), 6) ) +
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
  labs(title = paste0("Equitability index (Pielou)"), 
       x = "Months after gene therapy", 
       y = "Pielou index (J)", 
       colour = "CellMarker", fill = "CellMarker")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Pielou.ByMarker.pdf", sep = ""), height=8, width=12)
plot(plot_jindex)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Pielou.ByMarker.png", sep = ""), height=8, width=12, units = "in", res = 300)
plot(plot_jindex)
dev.off()


plot_hindex_larger <-
  ggplot(data = all_results_h_slice, aes(x = TimePoint, y = Hindex, fill = CellType, color = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(alpha = 0.2, size = 2) +
  # geom_jitter(alpha = 0.8, size = 4) +
  # geom_line(size=3, alpha = .7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  scale_x_continuous(breaks = seq(0, max(all_results_h_slice$TimePoint, na.rm = T), 12) ) +
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
  labs(title = paste0("Diversity index (Shannon)"), 
       x = "Months after gene therapy", 
       y = "Shannon diversity index (H')", 
       colour = "CellMarker", fill = "CellMarker")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.larger.pdf", sep = ""), height=7, width=10)
plot(plot_hindex_larger)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.larger.png", sep = ""), height=7, width=10, units = "in", res = 300)
plot(plot_hindex_larger)
dev.off()

plot_hindex_zscaled <-
  ggplot(data = all_results_h_slice_scaled, aes(x = TimePoint, y = h_zscaled, fill = CellType, color = CellType), na.rm = T, se = TRUE) +
  # ggplot(data = all_results_h_slice_ngdna_scaled, aes(x = TimePoint, y = h_zscaled, fill = CellType, color = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(alpha = 0.2, size = 2) +
  # geom_jitter(alpha = 0.8, size = 4) +
  # geom_line(size=3, alpha = .7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  scale_x_continuous(breaks = seq(0, max(all_results_h_slice$TimePoint, na.rm = T), 6) ) +
  facet_grid(Tissue ~ Study, scales = "free_x", space = "free_x") +
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
  labs(title = paste0("Clonal population diversity over time"), subtitle = "Z-score by patient, tissue, timepoint.",
       x = "Months after gene therapy", 
       y = "Shannon diversity index (scaled H' Z-score)", 
       colour = "CellMarker", fill = "CellMarker")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Divresity.Rescaled.ByMarker.Zscore.pdf", sep = ""), height=8, width=12)
plot(plot_hindex_zscaled)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Divresity.Rescaled.ByMarker.Zscore.png", sep = ""), height=8, width=12, units = "in", res = 300)
plot(plot_hindex_zscaled)
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


plot_hindex_norm <-
  ggplot(data = all_results_h_slice_scaled, aes(x = TimePoint, y = Hindex_normalized, fill = CellType, color = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(alpha = 0.2, size = 2) +
  # geom_line(size=3, alpha = .7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  scale_x_continuous(breaks = seq(0, max(all_results_h_slice$TimePoint, na.rm = T), 6) ) +
  facet_grid( Tissue ~ Study, scales = "free_x", space = "free") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Diversity index (Shannon)"), subtitle = "Normalized LAM-PCR samples.",
       x = "Months after gene therapy", 
       y = "Shannon diversity index (H')", 
       colour = "CellMarker", fill = "CellMarker")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.Normalized.pdf", sep = ""), height=8, width=12)
plot(plot_hindex_norm)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.Normalized.png", sep = ""), height=8, width=12, units = "in", res = 300)
plot(plot_hindex_norm)
dev.off()

plot_hindex_norm_zscaled <-
  ggplot(data = all_results_h_slice_scaled, aes(x = TimePoint, y = h_normalized_zscaled, fill = CellType, color = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(alpha = 0.2, size = 2) +
  # geom_line(size=3, alpha = .7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  scale_x_continuous(breaks = seq(0, max(all_results_h_slice$TimePoint, na.rm = T), 6) ) +
  facet_grid( Tissue ~ Study, scales = "free_x", space = "free") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Diversity index (Shannon)"), subtitle = "Normalized LAM-PCR samples. Z-score by patient, tissue, and timepoint.",
       x = "Months after gene therapy", 
       y = "Shannon diversity index (H')", 
       colour = "CellMarker", fill = "CellMarker")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.Normalized.Zscore.pdf", sep = ""), height=8, width=12)
plot(plot_hindex_norm_zscaled)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.Normalized.Zscore.png", sep = ""), height=8, width=12, units = "in", res = 300)
plot(plot_hindex_norm_zscaled)
dev.off()


plot_hindex_norm_homo <-
  ggplot(data = all_results_h_slice_scaled_homo, aes(x = TimePoint, y = Hindex_normalized, fill = CellType, color = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(alpha = 0.2, size = 2) +
  # geom_line(size=3, alpha = .7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  scale_x_continuous(breaks = seq(0, max(all_results_h_slice$TimePoint, na.rm = T), 6) ) +
  facet_grid( Tissue ~ Study, scales = "free_x", space = "free") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Diversity index (Shannon)"), subtitle = "Normalized LAM-PCR samples.",
       x = "Months after gene therapy", 
       y = "Shannon diversity index (H')", 
       colour = "CellMarker", fill = "CellMarker")


### zscore stats
all_results_h_slice_scaled_homo_stable <- all_results_h_slice_scaled_homo[which(all_results_h_slice_scaled_homo$TimepointMonths >= 24 & all_results_h_slice_scaled_homo$TimepointMonths <= 62),]
all_results_h_slice_scaled_homo_stable <- all_results_h_slice_scaled_homo_stable %>% group_by(SubjectID, Tissue, TimepointMonths) %>% mutate(n_obs_stt_overtime = n())
all_results_h_slice_scaled_homo_stable_min3 <- all_results_h_slice_scaled_homo_stable %>% group_by(SubjectID, Tissue, TimepointMonths) %>% mutate(n_obs_stt_overtime = n()) %>% filter(n_obs_stt_overtime >= 3)


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






# try statistical comparisons
# tutorial on repeated measures in R and ANOVA: https://www.r-bloggers.com/2021/04/repeated-measures-of-anova-in-r-complete-tutorial/

# plot grid data raw
plot_hindex_by_marker <-
  ggplot(data = all_results_h_slice, aes(x = TimePoint, y = Hindex, fill = CellType, color = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(alpha = 0.2, size = 2) +
  # geom_jitter(alpha = 0.8, size = 4) +
  # geom_line(size=3, alpha = .7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  scale_x_continuous(breaks = seq(0, max(all_results_h_slice$TimePoint, na.rm = T), 6) ) +
  facet_grid( CellLineage ~ Study, scales = "free_x", space = "free") +
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
  labs(title = paste0("Diversity index (Shannon)"), 
       x = "Months after gene therapy", 
       y = "Shannon diversity index (H')", 
       colour = "CellMarker", fill = "CellMarker")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker_grid.pdf", sep = ""), height=12, width=12)
plot(plot_hindex_by_marker)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker_grid.png", sep = ""), height=12, width=12, units = "in", res = 300)
plot(plot_hindex_by_marker)
dev.off()

plot_hindex_by_marker_study <-
  ggplot(data = all_results_h_slice, aes(x = TimePoint, y = Hindex, fill = Study, color = Study), na.rm = T, se = TRUE) +
  scale_color_manual(values = trials_colors) +
  scale_fill_manual(values = trials_colors) +
  # geom_point(alpha = 0.2, size = 2) +
  # geom_jitter(alpha = 0.8, size = 4) +
  # geom_line(size=3, alpha = .7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  scale_x_continuous(breaks = seq(0, max(all_results_h_slice$TimePoint, na.rm = T), 6) ) +
  facet_grid( . ~ CellLineage, scales = "free_x", space = "free") +
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
  labs(title = paste0("Diversity index (Shannon)"), 
       x = "Months after gene therapy", 
       y = "Shannon diversity index (H')", 
       colour = "Study", fill = "Study")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker_grid.pdf", sep = ""), height=7, width=22)
plot(plot_hindex_by_marker_study)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker_grid.png", sep = ""), height=7, width=22, units = "in", res = 300)
plot(plot_hindex_by_marker_study)
dev.off()


### ------- group and repeated measures ------
# slice data at stability
all_results_h_slice_scaled_stable <- all_results_h_slice_scaled[which(all_results_h_slice_scaled$TimePoint >= 24 & all_results_h_slice_scaled$TimePoint <= 60),]
all_results_h_slice_scaled_stable <- all_results_h_slice_scaled_stable %>% group_by(SubjectID, CellType) %>% mutate(n_celltype_overtime = n())
all_results_h_slice_scaled_stable_min3 <- all_results_h_slice_scaled_stable %>% group_by(SubjectID, CellType) %>% mutate(n_celltype_overtime = n()) %>% filter(n_celltype_overtime >= 3)

plot_hindex_by_marker_study_boxplot_overtime <-
  ggplot(data = all_results_h_slice_scaled_stable_min3, aes(x = Study, y = Hindex, fill = Study, color = Study), na.rm = T, se = TRUE) +
  scale_color_manual(values = trials_colors) +
  scale_fill_manual(values = trials_colors) +
  # geom_point(alpha = 0.2, size = 2) +
  geom_boxplot(alpha = .7, outlier.size = 0) +
  geom_jitter(alpha = 0.8, size = 4) +
  # geom_line(size=3, alpha = .7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  # geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  # scale_x_continuous(breaks = seq(0, max(all_results_h_slice$TimePoint, na.rm = T), 6) ) +
  facet_grid( . ~ CellType, scales = "free_x", space = "free") +
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
  labs(title = paste0("Diversity index (Shannon)"), subtitle = "Filters: timepoints >=24 & <60; min observed tp = 3.",
       x = "Study", 
       y = "Shannon diversity index (H')", 
       colour = "Study", fill = "Study")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker_grid.stable_t.boxplot_allpointsovertime.pdf", sep = ""), height=7, width=14)
plot(plot_hindex_by_marker_study_boxplot_overtime)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker_grid.stable_t.boxplot_allpointsovertime.png", sep = ""), height=7, width=14, units = "in", res = 300)
plot(plot_hindex_by_marker_study_boxplot_overtime)
dev.off()

# summary stats
patient_h_summary <- all_results_h_slice_scaled_stable_min3 %>%
  filter(n_celltype_overtime >= 3) %>%
  group_by(SubjectID, CellType) %>%
  get_summary_stats(Hindex, type = "mean_sd") 

# outliers
patient_h_summary_outliers <- all_results_h_slice_scaled_stable_min3 %>%
  filter(n_celltype_overtime >= 3) %>%
  group_by(SubjectID, CellType) %>%
  identify_outliers(Hindex)
# filter data
all_results_h_slice_scaled_stable_min3_nooutliers <- all_results_h_slice_scaled_stable_min3 %>%
  anti_join(patient_h_summary_outliers, by = c("Study", "SubjectID", "CellType", "CellMarker", "Tissue", "TimepointMonths"))

# normality
patient_h_summary_normality <- all_results_h_slice_scaled_stable_min3_nooutliers %>%
  group_by(SubjectID, CellType) %>%
  shapiro_test(Hindex)
# remove data that are not normnally distributed (p>0.05):: Tested data was normally distributed at each time point, as assessed by Shapiro-Wilk’s test (p > 0.05).
all_results_h_slice_scaled_stable_min3_nooutliers_ext <- merge(x = all_results_h_slice_scaled_stable_min3_nooutliers, y = patient_h_summary_normality, by = c("SubjectID", "CellType"), all.x = T)
all_results_h_slice_scaled_stable_min3_nooutliers_normal <- all_results_h_slice_scaled_stable_min3_nooutliers_ext %>%
  filter(p > 0.5)

plot_hindex_by_marker_study_boxplot_overtime_postfiltering <-
  ggplot(data = all_results_h_slice_scaled_stable_min3_nooutliers_normal, aes(x = Study, y = Hindex, fill = Study, color = Study), na.rm = T, se = TRUE) +
  scale_color_manual(values = trials_colors) +
  scale_fill_manual(values = trials_colors) +
  # geom_point(alpha = 0.2, size = 2) +
  geom_boxplot(alpha = .7, outlier.size = 0) +
  geom_jitter(alpha = 0.8, size = 4) +
  # geom_line(size=3, alpha = .7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  # geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  # scale_x_continuous(breaks = seq(0, max(all_results_h_slice$TimePoint, na.rm = T), 6) ) +
  facet_grid( . ~ CellType, scales = "free_x", space = "free") +
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
  labs(title = paste0("Diversity index (Shannon)"), subtitle = "Filters: timepoints >=24 & <60; min observed tp = 3. Then removed outliers and not normal samples.",
       x = "Study", 
       y = "Shannon diversity index (H')", 
       colour = "Study", fill = "Study")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker_grid.stable_t.boxplot_allpointsovertime.postfiltering.pdf", sep = ""), height=7, width=14)
plot(plot_hindex_by_marker_study_boxplot_overtime_postfiltering)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker_grid.stable_t.boxplot_allpointsovertime.postfiltering.png", sep = ""), height=7, width=14, units = "in", res = 300)
plot(plot_hindex_by_marker_study_boxplot_overtime_postfiltering)
dev.off()

# summary stats
patient_h_summary_postfiltering <- all_results_h_slice_scaled_stable_min3_nooutliers_normal %>%
  group_by(Study, SubjectID, CellType) %>%
  get_summary_stats(Hindex, type = "mean_sd") 

write.xlsx(x = patient_h_summary_postfiltering, file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker.stable_t.postfiltering.summarystats.toPRISM.xlsx", sep = ""))
# se non trovi di meglio mettili in PRISM

plot_hindex_by_marker_study_boxplot_overtime_postfiltering_mean <-
  ggplot(data = patient_h_summary_postfiltering, aes(x = SubjectID, y = mean, fill = Study, color = Study), na.rm = T, se = TRUE) +
  scale_color_manual(values = trials_colors) +
  scale_fill_manual(values = trials_colors) +
  # geom_point(alpha = 0.2, size = 2) +
  geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd)) +
  # geom_boxplot(alpha = .7, outlier.size = 0) +
  # geom_jitter(alpha = 0.8, size = 4) +
  geom_line(size=3, alpha = .7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  # geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  # scale_x_continuous(breaks = seq(0, max(all_results_h_slice$TimePoint, na.rm = T), 6) ) +
  facet_grid( . ~ CellType, scales = "free_x", space = "free") +
  theme_bw() +
  # scale_y_log10() + 
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  # facet_grid(ProjectID ~ ., scales = "free_y", space = "free") +
  # facet_grid(ProjectID ~ ., scales = "free_y") +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16, angle=55, hjust=1, vjust=1), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Diversity index (Shannon)"), subtitle = "Filters: timepoints >=24 & <60; min observed tp = 3. Then removed outliers and not normal samples.",
       x = "Study", 
       y = "Shannon diversity index (H')", 
       colour = "Study", fill = "Study")

pdf(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker_grid.stable_t.scatterplot_mean_sd.postfiltering.pdf", sep = ""), height=7, width=18)
plot(plot_hindex_by_marker_study_boxplot_overtime_postfiltering_mean)
dev.off()
png(file = paste(h_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".Diversity.ByMarker_grid.stable_t.scatterplot_mean_sd.postfiltering.png", sep = ""), height=7, width=18, units = "in", res = 300)
plot(plot_hindex_by_marker_study_boxplot_overtime_postfiltering_mean)
dev.off()



##### in progress ........
# res<-anova_test(data=all_results_h_slice_scaled_stable_min3_nooutliers_normal, dv=Hindex, wid=c(SubjectID, CellType), within=TimepointMonths) 
# get_anova_table(res) 

# from https://www.sheffield.ac.uk/polopoly_fs/1.885219!/file/105_RepeatedANOVA.pdf
library(ez)
repeat1 <- ezANOVA(data=all_results_h_slice_scaled_stable_min3_nooutliers_normal,
                   dv=.(Hindex),
                   wid=.(SubjectID, CellType),
                   within=.(TimepointMonths),
                   type=3)


###############################################################
## Cumulative
###############################################################
markerlist <- c("CD13", "CD14", "CD15", "CD19", "CD3", "CD34", "GLY", "CD36", "GLYA")
cum_base_folder <- paste0("analyses/cumulative/")
dir.create(file.path(getwd(), cum_base_folder), showWarnings = FALSE)
dir.create(file.path(cum_base_folder, analysis_folder_date), showWarnings = FALSE)

mld_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial MLD/"
bthal_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial bThalassemia/"
was_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial WAS/"

## ----- per patient ------
mld_patcum <- read.csv(file = paste0(mld_wd, "analyses/03.cumulative/", analysis_folder_date, "/20220110_MLD_SubjectID_TimepointMonths_cumulated_counts.tsv"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
mld_patcum$ClinicalTrial <- "MLD"
mld_patcum <- mld_patcum[which(mld_patcum$SubjectID %in% grep("MLD", levels(factor(mld_patcum$SubjectID)), value = T)),]

was_patcum <- read.csv(file = paste0(was_wd, "analyses/03.cumulative/", analysis_folder_date, "/20221001_WAS_SubjectID_TimepointMonths_cumulated_counts.tsv"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
was_patcum$ClinicalTrial <- "WAS"
was_patcum <- was_patcum[which(was_patcum$SubjectID %in% grep("WAS", levels(factor(was_patcum$SubjectID)), value = T)),]

bthal_patcum <- read.csv(file = paste0(bthal_wd, "analyses/12.cumulative/", analysis_folder_date, "/20220110_BTHAL_SubjectID_TimepointMonths_cumulated_counts.tsv"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
bthal_patcum$ClinicalTrial <- "BTHAL"
bthal_patcum <- bthal_patcum[which(bthal_patcum$SubjectID %in% grep("BTHAL", levels(factor(bthal_patcum$SubjectID)), value = T)),]

# merge
common_labels <- Reduce(intersect, list(colnames(was_patcum), colnames(mld_patcum), colnames(bthal_patcum)))
all_results_cumulative <- Reduce(rbind, list(mld_patcum[common_labels], was_patcum[common_labels], bthal_patcum[common_labels]))

all_results_cumulative$Study <- factor(all_results_cumulative$ClinicalTrial, levels = study_list)
all_results_cumulative[is.na(all_results_cumulative)] <- NA
# all_results_cumulative$Hindex <- all_results_cumulative$fragmentEstimate_sum_shannon
# all_results_cumulative$gdf_names <- rownames(all_results_cumulative)
all_results_cumulative$TimePoint <- all_results_cumulative$TimepointMonths
all_results_cumulative_max <- all_results_cumulative %>% group_by(SubjectID, Study) %>% summarise(count = max(count), TimePoint = max(TimePoint))

write.xlsx(x = all_results_cumulative, file = paste(cum_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".CumulativeNIS.ByPatient.xlsx", sep = ""), sheetName = "per patient")

# do plots
plot_cumulative_by_patient <-
  ggplot(data = all_results_cumulative, aes(x = TimePoint, y = count, fill = Study, color = Study), na.rm = T, se = TRUE) +
  scale_color_manual(values = trials_colors) +
  scale_fill_manual(values = trials_colors) +
  geom_point(alpha = 0.2, size = 2) +
  # geom_jitter(alpha = 0.8, size = 4) +
  geom_line(aes(group = SubjectID), size=1.5, alpha = .2, colour = "gray") +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  # geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  geom_smooth(method = "loess", formula = y ~ log(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  scale_x_continuous(breaks = seq(0, max(all_results_cumulative$TimePoint, na.rm = T), 6) ) +
  facet_grid( . ~ Study, scales = "free_x", space = "free") +
  theme_bw() +
  geom_label_repel(data = subset(all_results_cumulative_max, count>100000),
                   aes(label = SubjectID),
                   color = 'white',
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.4, "lines"),
                   segment.color = 'black') +
  scale_y_continuous(labels = scales::comma) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  # facet_grid(ProjectID ~ ., scales = "free_y", space = "free") +
  # facet_grid(ProjectID ~ ., scales = "free_y") +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Cumulative number of IS"), 
       x = "Months after gene therapy", 
       y = "Cumulative n. IS", 
       colour = "Study", fill = "Study")

pdf(file = paste(cum_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".CumulativeNIS.ByPatient.pdf", sep = ""), height=5.5, width=13)
plot(plot_cumulative_by_patient)
dev.off()
png(file = paste(cum_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".CumulativeNIS.ByPatient.png", sep = ""), height=5.5, width=13, units = "in", res = 300)
plot(plot_cumulative_by_patient)
dev.off()

# do plots
plot_cumulative_by_patient_noLabs <-
  ggplot(data = all_results_cumulative, aes(x = TimePoint, y = count, fill = Study, color = Study), na.rm = T, se = TRUE) +
  scale_color_manual(values = trials_colors) +
  scale_fill_manual(values = trials_colors) +
  geom_point(alpha = 0.2, size = 2) +
  # geom_jitter(alpha = 0.8, size = 4) +
  geom_line(aes(group = SubjectID), size=1.5, alpha = .2, colour = "gray") +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  # geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  geom_smooth(method = "loess", formula = y ~ log(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  scale_x_continuous(breaks = seq(0, max(all_results_cumulative$TimePoint, na.rm = T), 12) ) +
  facet_grid( . ~ Study, scales = "free_x", space = "free") +
  theme_bw() +
  # geom_label_repel(data = subset(all_results_cumulative_max, count>100000),
  #                  aes(label = SubjectID),
  #                  color = 'white',
  #                  box.padding = unit(0.35, "lines"),
  #                  point.padding = unit(0.4, "lines"),
  #                  segment.color = 'black') +
  scale_y_continuous(labels = scales::comma) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  # facet_grid(ProjectID ~ ., scales = "free_y", space = "free") +
  # facet_grid(ProjectID ~ ., scales = "free_y") +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Cumulative number of IS"), 
       x = "Months after gene therapy", 
       y = "Cumulative n. IS", 
       colour = "Study", fill = "Study")

pdf(file = paste(cum_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".CumulativeNIS.ByPatient.larger_noLabs.pdf", sep = ""), height=4.5, width=10)
plot(plot_cumulative_by_patient_noLabs)
dev.off()
png(file = paste(cum_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".CumulativeNIS.ByPatient.larger_noLabs.png", sep = ""), height=4.5, width=10, units = "in", res = 300)
plot(plot_cumulative_by_patient_noLabs)
dev.off()


## ----- per cell marker ------
mld_patcum <- read.csv(file = paste0(mld_wd, "analyses/03.cumulative/", analysis_folder_date, "/20220110_MLD_SubjectID_CellMarker_Tissue_TimepointMonths_cumulated_counts.tsv"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
mld_patcum$ClinicalTrial <- "MLD"
mld_patcum <- mld_patcum[which(mld_patcum$SubjectID %in% grep("MLD", levels(factor(mld_patcum$SubjectID)), value = T)),]
colname_to_convertto_cellmarker <- "SuperGroup"
if (colname_to_convertto_cellmarker %in% colnames(mld_patcum) & !("CellMarker" %in% colnames(mld_patcum))) {
  mld_patcum$CellMarker <- mld_patcum[,colname_to_convertto_cellmarker]
}

was_patcum <- read.csv(file = paste0(was_wd, "analyses/03.cumulative/", analysis_folder_date, "/20221001_WAS_SubjectID_SuperGroup_Tissue_TimepointMonths_cumulated_counts.tsv"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
was_patcum$ClinicalTrial <- "WAS"
was_patcum <- was_patcum[which(was_patcum$SubjectID %in% grep("WAS", levels(factor(was_patcum$SubjectID)), value = T)),]
colname_to_convertto_cellmarker <- "SuperGroup"
if (colname_to_convertto_cellmarker %in% colnames(was_patcum) & !("CellMarker" %in% colnames(was_patcum))) {
  was_patcum$CellMarker <- was_patcum[,colname_to_convertto_cellmarker]
}

bthal_patcum <- read.csv(file = paste0(bthal_wd, "analyses/12.cumulative/", analysis_folder_date, "/20220110_BTHAL_SubjectID_AggregationGroup_Tissue_TimepointMonths_cumulated_counts.tsv"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
bthal_patcum$ClinicalTrial <- "BTHAL"
bthal_patcum <- bthal_patcum[which(bthal_patcum$SubjectID %in% grep("BTHAL", levels(factor(bthal_patcum$SubjectID)), value = T)),]
colname_to_convertto_cellmarker <- "AggregationGroup"
if (colname_to_convertto_cellmarker %in% colnames(bthal_patcum) & !("CellMarker" %in% colnames(bthal_patcum))) {
  bthal_patcum$CellMarker <- bthal_patcum[,colname_to_convertto_cellmarker]
}


# merge
common_labels <- Reduce(intersect, list(colnames(was_patcum), colnames(mld_patcum), colnames(bthal_patcum)))
all_results_cumulative <- Reduce(rbind, list(mld_patcum[common_labels], was_patcum[common_labels], bthal_patcum[common_labels]))

all_results_cumulative$Study <- factor(all_results_cumulative$ClinicalTrial, levels = study_list)
all_results_cumulative[is.na(all_results_cumulative)] <- NA
all_results_cumulative$TimePoint <- all_results_cumulative$TimepointMonths
all_results_cumulative_max <- all_results_cumulative %>% group_by(CellMarker, SubjectID, Study) %>% summarise(count = max(count), TimePoint = max(TimePoint))

# merge blood lineage data
blood_lineages <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "blood_lineages")
# blood_lineages_colors <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "colors_noery")
blood_lineages_colors <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "colors")
blood_lineages <- merge(x = blood_lineages, y = blood_lineages_colors, all.x = T, by = c("CellType"))
rownames(blood_lineages) <- blood_lineages$CellMarker
all_results_cumulative <- merge(x = all_results_cumulative, y = blood_lineages, by = c("CellMarker"), all.x = T)

# write all
write.xlsx(x = all_results_cumulative, file = paste(cum_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".CumulativeNIS.ByPatient_CellMarker.xlsx", sep = ""), sheetName = "per cellmarker")

scale_color_manual_colors_sortedbyname <- blood_lineages_colors[order(blood_lineages_colors$CellType), "colorcode"]
all_results_cumulative_slice <- all_results_cumulative[which(all_results_cumulative$CellMarker %in% markerlist),]

# do plots
plot_cumulative_by_patient_cellmarker <-
  ggplot(data = all_results_cumulative_slice, aes(x = TimePoint, y = count, fill = CellType, color = CellType), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(alpha = 0.2, size = 2) +
  # geom_jitter(alpha = 0.8, size = 4) +
  geom_line(aes(group = interaction(SubjectID, CellMarker, Tissue)), size=1.5, alpha = .2, colour = "gray") +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  # geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  geom_smooth(method = "loess", formula = y ~ log(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  scale_x_continuous(breaks = seq(0, max(all_results_cumulative$TimePoint, na.rm = T), 6) ) +
  facet_grid( . ~ Study, scales = "free_x", space = "free") +
  theme_bw() +
  # geom_label_repel(data = subset(all_results_cumulative_max, count>100000),
  #                  aes(label = SubjectID),
  #                  color = 'white',
  #                  box.padding = unit(0.35, "lines"),
  #                  point.padding = unit(0.4, "lines"),
  #                  segment.color = 'black') +
  scale_y_continuous(labels = scales::comma) +
  scale_y_break(c(50000, 50000), scales = 0.1, ticklabels = c(60000, 90000, 120000)) + #scale_y_break(c(50000, 120000), scales=0.2) +
  # scale_y_cut(breaks=c(50000, 80000), which=c(1, 3), scales=c(0.5, 0.5, 4)) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  # facet_grid(ProjectID ~ ., scales = "free_y", space = "free") +
  # facet_grid(ProjectID ~ ., scales = "free_y") +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Cumulative number of IS"), 
       x = "Months after gene therapy", 
       y = "Cumulative n. IS", 
       colour = "CellType", fill = "CellType")

pdf(file = paste(cum_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".CumulativeNIS.ByPatient_CellType.ybreak.pdf", sep = ""), height=6, width=13)
# plot(plot_cumulative_by_patient_cellmarker)
plot_cumulative_by_patient_cellmarker
dev.off()
png(file = paste(cum_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".CumulativeNIS.ByPatient_CellType.ybreak.png", sep = ""), height=6, width=13, units = "in", res = 300)
plot_cumulative_by_patient_cellmarker
dev.off()


###############################################################
## CD34 output
###############################################################
followup_limit_crosstrial <- 60
patients_summary <- read.xlsx(xlsxFile = "source/Clinical/Patient Summary.xlsx", sheet = "summary")
allpatients_shared34_stats <- read.csv(file = paste("analyses/cd34output/202112/MLD_WAS_BTHAL.SharingIS_CD34BM_SuperGroup.allpatients_shared34_stats.tsv", sep = ""), header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NA", ''))
allpatients_shared34_stats <- merge(x = allpatients_shared34_stats, y = patients_summary, by = c("SubjectID"), all.x = T)
# allpatients_shared34_stats$AgeGroup 
allpatients_shared34_stats <- allpatients_shared34_stats %>%
  mutate(AgeGroup = case_when(
    TreatmentAge.years <= 2 ~"0-2",
    TreatmentAge.years <= 15 ~"2-15",
    TreatmentAge.years > 15 ~"30+"
  ))
age_list <- c("0-2", "2-15", "30+")
allpatients_shared34_stats$AgeGroup_Sorted <- factor(allpatients_shared34_stats$AgeGroup, levels = age_list)

library(openxlsx)
blood_lineages <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "blood_lineages")
blood_lineages_colors <- read.xlsx(xlsxFile = blood_lineages_file, sheet = "colors")
blood_lineages <- merge(x = blood_lineages, y = blood_lineages_colors, all.x = T, by = c("CellType"))
rownames(blood_lineages) <- blood_lineages$CellMarker
scale_color_manual_colors_sortedbyname <- blood_lineages_colors[order(blood_lineages_colors$CellType), "colorcode"]
allpatients_shared34_stats <- merge(x = allpatients_shared34_stats, y = blood_lineages, by = c("CellMarker"), all.x = T)

### do plots and stats ###
sqldf("select SubjectID, min(TimePoint) as Min_CD34_InVivo_TP from allpatients_shared34_stats where CellMarker = 'CD34' and TimePoint > 0 group by SubjectID order by Min_CD34_InVivo_TP desc ") # with thi squery, you know the patients to exclude from the analysis
patients_to_exclude <- c("MLD22", "MLDCUP01", "MLDCUP02")

# patients_without_erythroid <- c("MLD08", "MLD11", "MLD14", "MLD15", "MLD16", "MLD17", "MLD20", "MLDC02", "MLDCUP04", "MLDHE01", "MLDCUP05", "MLDCUP03", "WAS1009")
patients_without_erythroid <- c("MLDCUP04", "MLDHE01", "MLDCUP05", "MLDCUP03", "MLDC02", "WAS1009")
patients_uncertainty <- c("MLD12")

allpatients_shared34_stats_noCD34 <- allpatients_shared34_stats[which(!(allpatients_shared34_stats$CellMarker %in% c("CD34", "Plasma")) & 
                                                                        !(is.na(allpatients_shared34_stats$colorcode)) ),]
allpatients_shared34_stats_noCD34$PercNIS_onOverallSharedCD34BM <- (allpatients_shared34_stats_noCD34$SharedCD34IS_NumIS * 100) / allpatients_shared34_stats_noCD34$Progenitor_nIS
allpatients_shared34_stats_noCD34$Study <- factor(allpatients_shared34_stats_noCD34$StudyID, levels = study_list)

allpatients_shared34_stats_noCD34 <- allpatients_shared34_stats_noCD34 %>% group_by(SubjectID, Tissue, TimepointMonths) %>% mutate(PercNIS_onOverallSharedCD34BM_zscaled=scale(PercNIS_onOverallSharedCD34BM))
allpatients_shared34_stats_noCD34[is.na(allpatients_shared34_stats_noCD34)] <- NA

# 
# allpatients_shared34_stats_noCD34 <- allpatients_shared34_stats[which(!(allpatients_shared34_stats$CellMarker %in% c("CD34", "Plasma")) ),]
# allpatients_shared34_stats_noCD34 <- allpatients_shared34_stats[which(!(allpatients_shared34_stats$CellMarker %in% c("CD34", "Plasma", "CD36")) & 
#                                                                         !(allpatients_shared34_stats$StudyID == "WAS" & allpatients_shared34_stats$PCRMethod %in% c("SLiM|LAM-PCR", "SLiM") ) & 
#                                                                         !(allpatients_shared34_stats_noCD34$SubjectID %in% c(patients_to_exclude, patients_without_erythroid, patients_uncertainty)) 
#                                                                       ),]

# # allpatients_shared34_stats_supergroup <- read.csv(file = paste("/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/Clinical Trial MLD/analyses/20.compare_trials/202004/202004.MLD_WAS_BTHAL.SharingIS_CD34BM.tsv", sep = ""), header=TRUE, fill=T, sep='\t', check.names = F, row.names = 1)
# allpatients_shared34_stats_supergroup_noCD34 <- allpatients_shared34_stats_noCD34[which(!(allpatients_shared34_stats_noCD34$CellMarker %in% c("CD34")) &
#                                                                                           !(allpatients_shared34_stats_noCD34$SubjectID %in% c(patients_to_exclude, patients_without_erythroid, patients_uncertainty)) ),]
# allpatients_shared34_stats_supergroup_noCD34 <- allpatients_shared34_stats_noCD34
# allpatients_shared34_stats_supergroup_noCD34$Study <- factor(allpatients_shared34_stats_supergroup_noCD34$StudyID, levels = study_list)
# 
# allpatients_shared34_stats_supergroup_pediatriconly_noCD34 <- allpatients_shared34_stats_noCD34[which(!(allpatients_shared34_stats_noCD34$CellMarker %in% c("CD34")) & 
#                                                                                                         !(allpatients_shared34_stats_noCD34$SubjectID %in% c(patients_to_exclude, patients_without_erythroid, patients_uncertainty)) ),]
# allpatients_shared34_stats_supergroup_pediatriconly_noCD34$PercNIS_onOverallSharedCD34BM <- (allpatients_shared34_stats_supergroup_pediatriconly_noCD34$SharedCD34IS_NumIS * 100) / allpatients_shared34_stats_supergroup_pediatriconly_noCD34$Progenitor_nIS
# allpatients_shared34_stats_supergroup_pediatriconly_noCD34 <- allpatients_shared34_stats_supergroup_pediatriconly_noCD34[which(allpatients_shared34_stats_supergroup_pediatriconly_noCD34$PercNIS_onOverallSharedCD34BM > 0),]

allpatients_shared34_stats_supergroup_noCD34 <- allpatients_shared34_stats_noCD34

## ----- sandbox -----
        # from: https://m-clark.github.io/mixed-models-with-R/random_intercepts.html 
        # try linear mixed model effects
        library(lme4)
        gpa_lm = lm(PercNIS_onOverallSharedCD34BM ~ FollowUp + CellType + Study, data = allpatients_shared34_stats_supergroup_noCD34)
        summary(gpa_lm)
        # gpa_mixed = lmer(PercNIS_onOverallSharedCD34BM ~ FollowUp + Study + (1 | CellType), data = allpatients_shared34_stats_supergroup_noCD34)
        gpa_mixed = lmer(PercNIS_onOverallSharedCD34BM ~ FollowUp + (1|Study) + (1 | CellType), data = allpatients_shared34_stats_supergroup_noCD34)
        summary(gpa_mixed)
        confint(gpa_mixed)
        library(merTools)
        predictInterval(gpa_mixed)   # for various model predictions, possibly with new data
        REsim(gpa_mixed)             # mean, median and sd of the random effect estimates
        plotREsim(REsim(gpa_mixed))  # plot the interval estimates
        predict_with_re = predict(gpa_mixed)
        mixed.lmer <- gpa_mixed
        qqnorm(resid(mixed.lmer))
        library(sjPlot)
## ---- end sandbox ----

# do plot of trends over time

scale_color_manual_colors_sortedbyname <- colors_lineages[1,]
# p_lines_mbte <- ggplot(data = allpatients_shared34_stats_supergroup_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = CellType), na.rm = T, se = TRUE) +
p_lines_mbte <- ggplot(data = allpatients_shared34_stats_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = CellType), na.rm = T, se = TRUE) +
  # p <- ggplot(data = allpatients_shared34_stats_supergroup_pediatriconly_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = HematoLineage), na.rm = T, se = TRUE) +
  # geom_point(size=4, alpha = .4) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(aes(shape = CellMarker), size=4, alpha = .6) +
  # geom_line(size=3, alpha = .3) +
  # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
  # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
  # stat_smooth(level = 2, size = 2) +
  # stat_smooth(aes(fill = CellType), level = 0.7, size = 2) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  # stat_smooth(aes(fill = CellType), method = lm, formula = y ~ splines::bs(x, 3), size = 2, level = 0.4, se = T) +
  # facet_wrap( ~ SubjectID, ncol = 4) +
  facet_grid(. ~ Study, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_noCD34$FollowUp, na.rm = T), 12) ) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 60) + # Add global p-value
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34.")

p_lines_mbte_identicalFU <- ggplot(data = allpatients_shared34_stats_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = CellType), na.rm = T, se = TRUE) +
  # p <- ggplot(data = allpatients_shared34_stats_supergroup_pediatriconly_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = HematoLineage), na.rm = T, se = TRUE) +
  # geom_point(size=4, alpha = .4) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(aes(shape = CellMarker), size=4, alpha = .6) +
  # geom_line(size=3, alpha = .3) +
  # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
  # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
  # stat_smooth(level = 2, size = 2) +
  # stat_smooth(aes(fill = CellType), level = 0.7, size = 2) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  # stat_smooth(aes(fill = CellType), method = lm, formula = y ~ splines::bs(x, 3), size = 2, level = 0.4, se = T) +
  # facet_wrap( ~ SubjectID, ncol = 4) +
  facet_grid(. ~ Study, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_noCD34$FollowUp, na.rm = T), 12), limits = c(0, followup_limit_crosstrial) ) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 60) + # Add global p-value
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34.")


p_lines_mbte_zscore <- ggplot(data = allpatients_shared34_stats_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM_zscaled, color = CellType), na.rm = T, se = TRUE) +
  # p <- ggplot(data = allpatients_shared34_stats_supergroup_pediatriconly_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = HematoLineage), na.rm = T, se = TRUE) +
  # geom_point(size=4, alpha = .4) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(aes(shape = CellMarker), size=4, alpha = .6) +
  # geom_line(size=3, alpha = .3) +
  # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
  # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
  # stat_smooth(level = 2, size = 2) +
  # stat_smooth(aes(fill = CellType), level = 0.7, size = 2) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  # stat_smooth(aes(fill = CellType), method = lm, formula = y ~ splines::bs(x, 3), size = 2, level = 0.4, se = T) +
  # facet_wrap( ~ SubjectID, ncol = 4) +
  facet_grid(. ~ Study, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_noCD34$FollowUp, na.rm = T), 12) ) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 60) + # Add global p-value
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), 
       x = "FollowUp months after GT", y = "Perc. of shared CD34 IS (Z-score)", 
       subtitle = "All clinical trials included, all patients included. Percentages are scaled (Z-score) by Patient, Tissue and Follow-up.")

p_lines_mbte_zscore_agegr <- ggplot(data = allpatients_shared34_stats_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM_zscaled, color = CellType), na.rm = T, se = TRUE) +
  # p <- ggplot(data = allpatients_shared34_stats_supergroup_pediatriconly_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = HematoLineage), na.rm = T, se = TRUE) +
  # geom_point(size=4, alpha = .4) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(aes(shape = CellMarker), size=4, alpha = .6) +
  # geom_line(size=3, alpha = .3) +
  # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
  # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
  # stat_smooth(level = 2, size = 2) +
  # stat_smooth(aes(fill = CellType), level = 0.7, size = 2) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  # stat_smooth(aes(fill = CellType), method = lm, formula = y ~ splines::bs(x, 3), size = 2, level = 0.4, se = T) +
  # facet_wrap( ~ SubjectID, ncol = 4) +
  facet_grid(AgeGroup_Sorted ~ Study, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_noCD34$FollowUp, na.rm = T), 12) ) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 60) + # Add global p-value
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), 
       x = "FollowUp months after GT", y = "Perc. of shared CD34 IS (Z-score)", 
       subtitle = "All clinical trials included, all patients included. Percentages are scaled (Z-score) by Patient, Tissue and Follow-up.")

p_lines_mbte_bypatient <- ggplot(data = allpatients_shared34_stats_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = CellType), na.rm = T, se = TRUE) +
  # p <- ggplot(data = allpatients_shared34_stats_supergroup_pediatriconly_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = HematoLineage), na.rm = T, se = TRUE) +
  # geom_point(size=4, alpha = .4) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(aes(shape = CellMarker), size=4, alpha = .6) +
  # geom_line(size=3, alpha = .3) +
  # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
  # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
  # stat_smooth(level = 2, size = 2) +
  stat_smooth(aes(fill = CellType), level = 0.7, size = 2) +
  # stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  # stat_smooth(aes(fill = CellType), method = lm, formula = y ~ splines::bs(x, 3), size = 2, level = 0.4, se = T) +
  facet_wrap( ~ SubjectID, ncol = 8, scales = "free_x") +
  # facet_grid(. ~ Study, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_supergroup_noCD34$FollowUp, na.rm = T), 12) ) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 60) + # Add global p-value
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34.")

p_lines_mbte_NIS <- ggplot(data = allpatients_shared34_stats_supergroup_noCD34, aes(x = FollowUp, y = nIS, color = CellType), na.rm = T, se = TRUE) +
  # p <- ggplot(data = allpatients_shared34_stats_supergroup_pediatriconly_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = HematoLineage), na.rm = T, se = TRUE) +
  # geom_point(size=4, alpha = .4) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(aes(shape = CellMarker), size=4, alpha = .6) +
  # geom_line(size=3, alpha = .3) +
  # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
  # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
  # stat_smooth(level = 2, size = 2) +
  # stat_smooth(aes(fill = CellType), level = 0.7, size = 2) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  # stat_smooth(aes(fill = CellType), method = lm, formula = y ~ splines::bs(x, 3), size = 2, level = 0.4, se = T) +
  # facet_wrap( ~ SubjectID, ncol = 4) +
  facet_grid(. ~ Study, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_supergroup_noCD34$FollowUp, na.rm = T), 12)) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 60) + # Add global p-value
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Number of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34.")


p_lines_mbte_agegr <- ggplot(data = allpatients_shared34_stats_supergroup_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = CellType), na.rm = T, se = TRUE) +
  # p <- ggplot(data = allpatients_shared34_stats_supergroup_pediatriconly_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = HematoLineage), na.rm = T, se = TRUE) +
  geom_point(size=4, alpha = .4) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(aes(shape = CellMarker), size=4, alpha = .6) +
  # geom_line(size=3, alpha = .3) +
  # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
  # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
  # stat_smooth(level = 2, size = 2) +
  # stat_smooth(aes(fill = CellType), level = 0.7, size = 2) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  # facet_wrap( ~ SubjectID, ncol = 4) +
  facet_grid(AgeGroup_Sorted ~ Study, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_supergroup_noCD34$FollowUp, na.rm = T), 12) ) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 60) + # Add global p-value
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34.")

p_lines_mbte_agegr_nopoints <- ggplot(data = allpatients_shared34_stats_supergroup_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = CellType), na.rm = T, se = TRUE) +
  # p <- ggplot(data = allpatients_shared34_stats_supergroup_pediatriconly_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = HematoLineage), na.rm = T, se = TRUE) +
  # geom_point(size=4, alpha = .4) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(aes(shape = CellMarker), size=4, alpha = .6) +
  # geom_line(size=3, alpha = .3) +
  # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
  # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
  # stat_smooth(level = 2, size = 2) +
  # stat_smooth(aes(fill = CellType), level = 0.7, size = 2) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  # facet_wrap( ~ SubjectID, ncol = 4) +
  facet_grid(AgeGroup_Sorted ~ Study, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_supergroup_noCD34$FollowUp, na.rm = T), 12) ) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 60) + # Add global p-value
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34.")

p_lines_mbte_agegr_nopoints_identicalFU <- ggplot(data = allpatients_shared34_stats_supergroup_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = CellType), na.rm = T, se = TRUE) +
  # p <- ggplot(data = allpatients_shared34_stats_supergroup_pediatriconly_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = HematoLineage), na.rm = T, se = TRUE) +
  # geom_point(size=4, alpha = .4) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(aes(shape = CellMarker), size=4, alpha = .6) +
  # geom_line(size=3, alpha = .3) +
  # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
  # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
  # stat_smooth(level = 2, size = 2) +
  # stat_smooth(aes(fill = CellType), level = 0.7, size = 2) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  # facet_wrap( ~ SubjectID, ncol = 4) +
  facet_grid(AgeGroup_Sorted ~ Study, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_supergroup_noCD34$FollowUp, na.rm = T), 12), limits = c(0, followup_limit_crosstrial) ) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 60) + # Add global p-value
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34.")

p_lines_mbte_agegr_nopoints_identicalFU_freey <- ggplot(data = allpatients_shared34_stats_supergroup_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = CellType), na.rm = T, se = TRUE) +
  # p <- ggplot(data = allpatients_shared34_stats_supergroup_pediatriconly_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = HematoLineage), na.rm = T, se = TRUE) +
  # geom_point(size=4, alpha = .4) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(aes(shape = CellMarker), size=4, alpha = .6) +
  # geom_line(size=3, alpha = .3) +
  # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
  # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
  # stat_smooth(level = 2, size = 2) +
  # stat_smooth(aes(fill = CellType), level = 0.7, size = 2) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  # facet_wrap( ~ SubjectID, ncol = 4) +
  facet_grid(AgeGroup_Sorted ~ Study, scales = "free", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_supergroup_noCD34$FollowUp, na.rm = T), 12), limits = c(0, followup_limit_crosstrial) ) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 60) + # Add global p-value
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34.")

p_lines_mbte_perpatient <- ggplot(data = allpatients_shared34_stats_supergroup_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = CellType), na.rm = T, se = TRUE) +
  # p <- ggplot(data = allpatients_shared34_stats_supergroup_pediatriconly_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = HematoLineage), na.rm = T, se = TRUE) +
  geom_point(size=4, alpha = .4) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(aes(shape = CellMarker), size=4, alpha = .6) +
  # geom_line(size=3, alpha = .3) +
  # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
  # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
  # stat_smooth(level = 2, size = 2) +
  # stat_smooth(aes(fill = CellType), level = 0.7, size = 2) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  facet_wrap( ~ SubjectID, ncol = 6) +
  # facet_grid(. ~ Study) +
  # scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_supergroup_noCD34$FollowUp, na.rm = T), 6) ) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 60) + # Add global p-value
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34.")


plot_footnote_details <- "All patients and all cell markers"
plot_filename <- paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.", sep = "")
plot_title <- paste0("CD34 BM outout over time")

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v2.wrap_trials.pdf", sep = ""), height=4.5, width=10)
print(p_lines_mbte)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v2.wrap_trials.png", sep = ""), height=4.5, width=10, units = "in", res = 300)
print(p_lines_mbte)
dev.off()

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v2.wrap_trials.60m.pdf", sep = ""), height=4.5, width=10)
print(p_lines_mbte_identicalFU)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v2.wrap_trials.60m.png", sep = ""), height=4.5, width=10, units = "in", res = 300)
print(p_lines_mbte_identicalFU)
dev.off()

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v2.wrap_trials.zscore.pdf", sep = ""), height=4.5, width=10)
print(p_lines_mbte_zscore)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v2.wrap_trials.zscore.png", sep = ""), height=4.5, width=10, units = "in", res = 300)
print(p_lines_mbte_zscore)
dev.off()

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v2.wrap_trials_age.pdf", sep = ""), 
    height=7, width=10)
print(p_lines_mbte_agegr)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v2.wrap_trials_age.png", sep = ""), 
    height=7, width=10, units = "in", res = 300)
print(p_lines_mbte_agegr)
dev.off()

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v2.zscoe.wrap_trials_age.pdf", sep = ""), 
    height=7, width=10)
print(p_lines_mbte_zscore_agegr)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v2.zscore.wrap_trials_age.png", sep = ""), 
    height=7, width=10, units = "in", res = 300)
print(p_lines_mbte_zscore_agegr)
dev.off()

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v2.wrap_trials_age.nopoints.pdf", sep = ""), 
    height=7, width=10)
print(p_lines_mbte_agegr_nopoints)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v2.wrap_trials_age.nopoints.png", sep = ""), 
    height=7, width=10, units = "in", res = 300)
print(p_lines_mbte_agegr_nopoints)
dev.off()

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v2.wrap_trials_age.nopoints.60m.pdf", sep = ""), 
    height=8, width=11)
print(p_lines_mbte_agegr_nopoints_identicalFU)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v2.wrap_trials_age.nopoints.60m.png", sep = ""), 
    height=8, width=11, units = "in", res = 300)
print(p_lines_mbte_agegr_nopoints_identicalFU)
dev.off()

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v2.wrap_trials_age.nopoints.60m.fy.pdf", sep = ""), 
    height=7, width=10)
print(p_lines_mbte_agegr_nopoints_identicalFU_freey)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v2.wrap_trials_age.nopoints.60m.fy.png", sep = ""), 
    height=7, width=10, units = "in", res = 300)
print(p_lines_mbte_agegr_nopoints_identicalFU_freey)
dev.off()


# png(file = plot_filename, height=20, width=12, units = "in", res = 300)
# grid.arrange(p_lines_mbte,
#              p_lines_mbte_agegr,
#              ncol = 1, 
#              top = textGrob(plot_title, gp = gpar(fontsize = 22)),
#              bottom = textGrob(plot_footnote_details, gp = gpar(fontface = 3, fontsize = 18), hjust=1, x = 1)
# )
# dev.off()


p_lines_mbte_trial <- ggplot(data = allpatients_shared34_stats_supergroup_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = Study), na.rm = T, se = TRUE) +
  # p <- ggplot(data = allpatients_shared34_stats_supergroup_pediatriconly_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = HematoLineage), na.rm = T, se = TRUE) +
  # geom_point(size=4, alpha = .4) +
  scale_color_manual(values = trials_colors) +
  scale_fill_manual(values = trials_colors) +
  # geom_point(aes(shape = CellMarker), size=4, alpha = .6) +
  # geom_line(size=3, alpha = .3) +
  # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
  # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
  # stat_smooth(level = 2, size = 2) +
  # stat_smooth(aes(fill = CellType), level = 0.7, size = 2) +
  stat_smooth(formula = y ~ log(x), level = 0.4, size = 2) +
  # stat_smooth(method = lm, formula = y ~ splines::bs(x, 3), size = 2, level = 0.4, se = T) +
  # facet_wrap( ~ SubjectID, ncol = 4) +
  facet_grid(. ~ CellType, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_supergroup_noCD34$FollowUp, na.rm = T), 6) ) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 60) + # Add global p-value
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16, angle = 45), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34.")


pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.byLineage.pdf", sep = ""), height=5, width=17)
print(p_lines_mbte_trial)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.byLineage.png", sep = ""), height=5, width=17, units = "in", res = 300)
print(p_lines_mbte_trial)
dev.off()


# p <- ggplot(data = allpatients_shared34_stats_supergroup_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = CellMarker), na.rm = T, se = TRUE)
# print(
#   p + 
#     geom_point(size=4, alpha = .7) +
#     geom_line(size=3, alpha = .7) +
#     # facet_wrap( ~ StudyID, ncol = 4) +
#     facet_grid(. ~ StudyID, scales = "free", space = "free") +
#     scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_supergroup_noCD34$FollowUp, na.rm = T), 6) ) +
#     # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
#     # stat_compare_means(label.y = 60) + # Add global p-value
#     theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
#     theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
#     labs(list(title = paste("Profile of CD34 BM output towards Myeloid/Lymphoid/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34."))
# )
# print(
#   p + 
#     # geom_point(size=4, alpha = .4) +
#     geom_point(aes(shape = CellMarker), size=4, alpha = .3) +
#     # geom_line(size=3, alpha = .3) +
#     # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
#     # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
#     stat_smooth(level = 2, size = 2) +
#     # facet_wrap( ~ SubjectID, ncol = 4) +
#     facet_grid(. ~ StudyID, scales = "free", space = "free") +
#     scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_supergroup_noCD34$FollowUp, na.rm = T), 6) ) +
#     # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
#     # stat_compare_means(label.y = 60) + # Add global p-value
#     theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
#     theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
#     labs(list(title = paste("Profile of CD34 BM output towards Myeloid/Lymphoid/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34."))
# ) 




### ------- group and repeated measures ------
# slice data at stability
allpatients_shared34_stats_supergroup_noCD34_stable <- allpatients_shared34_stats_supergroup_noCD34[which(allpatients_shared34_stats_supergroup_noCD34$TimepointMonths >= 24 & allpatients_shared34_stats_supergroup_noCD34$TimepointMonths <= 60),]
allpatients_shared34_stats_supergroup_noCD34_stable <- allpatients_shared34_stats_supergroup_noCD34_stable %>% group_by(SubjectID, Tissue, TimepointMonths) %>% mutate(n_obs_stt_overtime = n())
allpatients_shared34_stats_supergroup_noCD34_stable_min3 <- allpatients_shared34_stats_supergroup_noCD34_stable %>% group_by(SubjectID, Tissue, TimepointMonths) %>% mutate(n_obs_stt_overtime = n()) %>% filter(n_obs_stt_overtime >= 3)

p_lines_mbte_zscore_stable_min3 <-
  ggplot(data = allpatients_shared34_stats_supergroup_noCD34_stable_min3, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM_zscaled, color = CellType), na.rm = T, se = TRUE) +
  # p <- ggplot(data = allpatients_shared34_stats_supergroup_pediatriconly_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = HematoLineage), na.rm = T, se = TRUE) +
  # geom_point(size=2, alpha = .4) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(aes(shape = CellMarker), size=4, alpha = .6) +
  # geom_line(aes(SubjectID), size=1, alpha = .3) +
  # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
  # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
  # stat_smooth(level = 2, size = 2) +
  stat_smooth(aes(fill = CellType), level = 0.7, size = 2) +
  # stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  # stat_smooth(aes(fill = CellType), method = lm, formula = y ~ splines::bs(x, 3), size = 2, level = 0.4, se = T) +
  # facet_wrap( ~ SubjectID, ncol = 4) +
  facet_grid(. ~ Study, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_noCD34$FollowUp, na.rm = T), 12) ) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 60) + # Add global p-value
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), 
       x = "FollowUp months after GT", y = "Perc. of shared CD34 IS (Z-score)", 
       subtitle = "All clinical trials included, all patients included. Percentages are scaled (Z-score) by Patient, Tissue and Follow-up.")

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg.pdf", sep = ""), height=5, width=12)
plot(p_lines_mbte_zscore_stable_min3)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg.png", sep = ""), height=5, width=12, units = "in", res = 300)
plot(p_lines_mbte_zscore_stable_min3)
dev.off()

# summary stats
patient_shared34_summary <- allpatients_shared34_stats_supergroup_noCD34_stable_min3 %>%
  filter(n_obs_stt_overtime >= 3) %>%
  group_by(Study, SubjectID, Tissue, CellMarker, CellType, TreatmentAge.years, AgeGroup, AgeGroup_Sorted) %>%
  # get_summary_stats(PercNIS_onOverallSharedCD34BM_zscaled, type = "mean_sd")
  get_summary_stats(PercNIS_onOverallSharedCD34BM_zscaled, type = "full")

patient_shared34_summary_filtered <- patient_shared34_summary %>% filter(n >= 2)
patient_shared34_summary_filtered$CellLineage <- factor(patient_shared34_summary_filtered$CellType, levels = celltype_list)
patient_shared34_summary_filtered$Study <- factor(patient_shared34_summary_filtered$Study, levels = study_list)
patient_shared34_summary_filtered$AgeGroup_Sorted <- factor(patient_shared34_summary_filtered$AgeGroup, levels = age_list)

write.xlsx(x = patient_shared34_summary_filtered, file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg.patient_shared34_summary_filtered.full.xlsx", sep = ""))

# do plots
my_comparisons <- list(
  c("MLD", "WAS"),
  c("MLD", "BTHAL"),
  c("WAS", "BTHAL")
)
perc_spacer_stats <- -3
perc_starting_max_stats <- 18
p_lines_mbte_zscore_stable_min3_stats <- ggplot(data = patient_shared34_summary_filtered,
                                                aes(x = Study, y = mean, colour = CellMarker, fill = CellMarker, group = Study, shape = Tissue ), na.rm = T, se = TRUE) +
  # scale_color_manual(values = trials_colors) +
  # scale_fill_manual(values = trials_colors) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  # geom_violin(colour = "darkblue", fill = "lightblue", alpha = 0.6) +
  # geom_line(size=3, alpha = .7) +
  facet_wrap(. ~ CellLineage, ncol = 4) +
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
  labs(title = paste0("Percentage of shared clones in CD34 cells"), 
       subtitle = "Filters: purity, PB and BM separated, markers separated. Months under analysis: 24-60",
       x = "Studies", 
       y = "Sharing percentage", 
       colour = "CellMarker", fill = "CellMarker")

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg.byStudy.pdf", sep = ""), height=4.5, width=12)
plot(p_lines_mbte_zscore_stable_min3_stats)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg.byStudy.png", sep = ""), height=4.5, width=12, units = "in", res = 300)
plot(p_lines_mbte_zscore_stable_min3_stats)
dev.off()

p_lines_mbte_zscore_stable_min3_stat_median <- ggplot(data = patient_shared34_summary_filtered,
                                                aes(x = Study, y = median, colour = CellMarker, fill = CellMarker, group = Study, shape = Tissue ), na.rm = T, se = TRUE) +
  # scale_color_manual(values = trials_colors) +
  # scale_fill_manual(values = trials_colors) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  # geom_violin(colour = "darkblue", fill = "lightblue", alpha = 0.6) +
  # geom_line(size=3, alpha = .7) +
  facet_wrap(. ~ CellLineage, ncol = 4) +
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
  labs(title = paste0("Percentage of shared clones in CD34 cells"), 
       subtitle = "Filters: purity, PB and BM separated, markers separated. Months under analysis: 24-60",
       x = "Studies", 
       y = "Sharing percentage", 
       colour = "CellMarker", fill = "CellMarker")

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg_median.byStudy.pdf", sep = ""), height=4.5, width=12)
plot(p_lines_mbte_zscore_stable_min3_stat_median)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg_median.byStudy.png", sep = ""), height=4.5, width=12, units = "in", res = 300)
plot(p_lines_mbte_zscore_stable_min3_stat_median)
dev.off()

p_lines_mbte_zscore_stable_min3_stat_median_newcolors <- ggplot(data = patient_shared34_summary_filtered,
                                                      aes(x = Study, y = median, colour = Study, fill = Study, group = Study, shape = Tissue ), na.rm = T, se = TRUE) +
  scale_color_manual(values = trials_colors) +
  scale_fill_manual(values = trials_colors) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  # geom_violin(colour = "darkblue", fill = "lightblue", alpha = 0.6) +
  # geom_line(size=3, alpha = .7) +
  facet_wrap(. ~ CellLineage, ncol = 4) +
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
  labs(title = paste0("Percentage of shared clones in CD34 cells"), 
       subtitle = "Filters: purity, PB and BM separated, markers separated. Months under analysis: 24-60",
       x = "Studies", 
       y = "Sharing percentage", 
       colour = "Study", fill = "Study")

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg_median.byStudy.colors.pdf", sep = ""), height=6, width=12)
plot(p_lines_mbte_zscore_stable_min3_stat_median_newcolors)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg_median.byStudy.colors.png", sep = ""), height=6, width=12, units = "in", res = 300)
plot(p_lines_mbte_zscore_stable_min3_stat_median_newcolors)
dev.off()

p_lines_mbte_zscore_stable_min3_stats_2cols <- ggplot(data = patient_shared34_summary_filtered,
                                                # aes(x = Study, y = mean, colour = CellMarker, fill = CellMarker, group = Study, shape = Tissue ), na.rm = T, se = TRUE) +
                                                aes(x = Study, y = mean, group = Study, shape = Tissue ), na.rm = T, se = TRUE) +
  # scale_color_manual(values = trials_colors) +
  # scale_fill_manual(values = trials_colors) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  # geom_violin(colour = "darkblue", fill = "lightblue", alpha = 0.6) +
  # geom_line(size=3, alpha = .7) +
  facet_wrap(. ~ CellLineage, ncol = 2) +
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
  labs(title = paste0("Percentage of shared clones in CD34 cells"), 
       subtitle = "Filters: purity, PB and BM separated, markers separated. Months under analysis: 24-60",
       x = "Studies", 
       y = "Sharing percentage", 
       colour = "CellMarker", fill = "CellMarker")

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg.byStudy.2cols.v2.pdf", sep = ""), height=7, width=6)
plot(p_lines_mbte_zscore_stable_min3_stats_2cols)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg.byStudy.2cols.v2.png", sep = ""), height=7, width=6, units = "in", res = 300)
plot(p_lines_mbte_zscore_stable_min3_stats_2cols)
dev.off()

## age NOT working
# p_lines_mbte_zscore_stable_min3_stats_age <- ggplot(data = patient_shared34_summary_filtered,
#                                                 aes(x = Study, y = mean, colour = CellMarker, fill = CellMarker, group = Study ), na.rm = T, se = TRUE) +
#   # scale_color_manual(values = trials_colors) +
#   # scale_fill_manual(values = trials_colors) +
#   geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
#   # geom_point(alpha = 0.6, size = 4) +
#   geom_jitter(alpha = 0.6, size = 4) +
#   # geom_violin(colour = "darkblue", fill = "lightblue", alpha = 0.6) +
#   # geom_line(size=3, alpha = .7) +
#   facet_wrap(AgeGroup_Sorted ~ CellLineage, ncol = 4) +
#   scale_y_continuous(limits = c(-2, 2.3)) + 
#   theme_bw() +
#   # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
#   # stat_compare_means(label.y = 6) +
#   # stat_compare_means(comparisons = my_comparisons, label.y = c(30, 35, 40), p.adjust.method = "fdr", label = "p.signif") +
#   stat_compare_means(comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
#   theme(strip.text = element_text(face="bold", size=16)) +
#   # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
#   theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
#   # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
#   labs(title = paste0("Percentage of shared clones in CD34 cells"), 
#        subtitle = "Filters: purity, PB and BM separated, markers separated. Months under analysis: 24-60",
#        x = "Studies", 
#        y = "Sharing percentage", 
#        colour = "CellMarker", fill = "CellMarker")


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
p_lines_mbte_zscore_stable_min3_stats_lineage <- ggplot(data = patient_shared34_summary_filtered,
                                                            aes(x = CellLineage, y = mean, colour = CellMarker, fill = CellMarker, group = CellLineage, shape = Tissue), na.rm = T, se = TRUE) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  facet_wrap(. ~ Study, nrow = 1) +
  theme_bw() +
  scale_y_continuous(limits = c(-2, 3)) +
  stat_compare_means(comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16, angle = 0), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Percentage of shared clones in CD34 cells"), 
       subtitle = "Filters: purity, PB and BM separated, markers separated. Months under analysis: 24-60",
       x = "Lineage", 
       y = "Sharing percentage", 
       colour = "Cell Marker", fill = "Cell Marker")

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg.byLineage.pdf", sep = ""), height=5, width=14)
plot(p_lines_mbte_zscore_stable_min3_stats_lineage)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg.byLineage.png", sep = ""), height=5, width=14, units = "in", res = 300)
plot(p_lines_mbte_zscore_stable_min3_stats_lineage)
dev.off()

perc_spacer_stats <- -3
perc_starting_max_stats <- 26
p_lines_mbte_zscore_stable_min3_stats_lineage_age <- ggplot(data = patient_shared34_summary_filtered,
                                                        aes(x = CellLineage, y = mean, colour = CellMarker, fill = CellMarker, group = CellLineage, shape = Tissue), na.rm = T, se = TRUE) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  facet_wrap(AgeGroup_Sorted ~ Study, nrow = 2) +
  theme_bw() +
  scale_y_continuous(limits = c(-2, 3)) +
  stat_compare_means(comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16, angle = 20), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Percentage of shared clones in CD34 cells"), 
       subtitle = "Filters: purity, PB and BM separated, markers separated. Months under analysis: 24-60",
       x = "Lineage", 
       y = "Sharing percentage", 
       colour = "Cell Marker", fill = "Cell Marker")

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg.byLineage.age.v2.pdf", sep = ""), height=8, width=9)
plot(p_lines_mbte_zscore_stable_min3_stats_lineage_age)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg.byLineage.age.v2.png", sep = ""), height=8, width=9, units = "in", res = 300)
plot(p_lines_mbte_zscore_stable_min3_stats_lineage_age)
dev.off()

p_lines_mbte_zscore_stable_min3_stats_lineage_age_median <- ggplot(data = patient_shared34_summary_filtered,
                                                            # aes(x = CellLineage, y = median, colour = CellMarker, fill = CellMarker, group = CellLineage, shape = Tissue), na.rm = T, se = TRUE) +
                                                            aes(x = CellLineage, y = median, group = CellLineage, shape = Tissue), na.rm = T, se = TRUE) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  facet_wrap(AgeGroup_Sorted ~ Study, nrow = 2) +
  theme_bw() +
  scale_y_continuous(limits = c(-2, 3)) +
  stat_compare_means(comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16, angle = 20), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Percentage of shared clones in CD34 cells"), 
       subtitle = "Filters: purity, PB and BM separated, markers separated. Months under analysis: 24-60",
       x = "Lineage", 
       y = "Sharing percentage", 
       colour = "Cell Marker", fill = "Cell Marker")

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg_median.byLineage.age.v2.pdf", sep = ""), height=8, width=8)
plot(p_lines_mbte_zscore_stable_min3_stats_lineage_age_median)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg_median.byLineage.age.v2.png", sep = ""), height=8, width=8, units = "in", res = 300)
plot(p_lines_mbte_zscore_stable_min3_stats_lineage_age_median)
dev.off()



## find among lineages
my_comparisons <- list(
  c("0-2", "2-15")
)

perc_spacer_stats <- -3
perc_starting_max_stats <- 24
p_lines_mbte_zscore_stable_min3_stats_age_median_p1 <- ggplot(data = patient_shared34_summary_filtered,
                                                        aes(x = AgeGroup_Sorted, y = median, group = AgeGroup, shape = Tissue), na.rm = T, se = TRUE) +
                                                        # aes(x = AgeGroup_Sorted, y = median, colour = CellMarker, fill = CellMarker, group = AgeGroup, shape = Tissue), na.rm = T, se = TRUE) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  # facet_wrap(CellLineage ~ Study, nrow = 1) +
  facet_grid(CellLineage ~ Study, scales = "free_x") +
  theme_bw() +
  scale_y_continuous(limits = c(-2, 3)) +
  stat_compare_means(comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(2.5, 2), p.adjust.method = "fdr") +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16, angle = 0), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Percentage of shared clones in CD34 cells"), 
       subtitle = "Filters: purity, PB and BM separated, markers separated. Months under analysis: 24-60",
       x = "Lineage", 
       y = "Sharing percentage", 
       colour = "Cell Marker", fill = "Cell Marker")

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg_median.byLineage.p1.v2.pdf", sep = ""), height=7, width=7)
plot(p_lines_mbte_zscore_stable_min3_stats_age_median_p1)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg_median.byLineage.p1.v2.png", sep = ""), height=7, width=7, units = "in", res = 300)
plot(p_lines_mbte_zscore_stable_min3_stats_age_median_p1)
dev.off()

my_comparisons <- list(
  c("2-15", "30+")
)
p_lines_mbte_zscore_stable_min3_stats_age_median_p2 <- ggplot(data = patient_shared34_summary_filtered[which(patient_shared34_summary_filtered$Study == "BTHAL"),],
                                                              # aes(x = AgeGroup_Sorted, y = median, colour = CellMarker, fill = CellMarker, group = AgeGroup, shape = Tissue), na.rm = T, se = TRUE) +
                                                              aes(x = AgeGroup_Sorted, y = median, group = AgeGroup, shape = Tissue), na.rm = T, se = TRUE) +
  geom_boxplot(colour = "lightblue", fill = "lightblue", alpha = .3, outlier.size = 0) +
  # geom_point(alpha = 0.6, size = 4) +
  geom_jitter(alpha = 0.6, size = 4) +
  # facet_wrap(CellLineage ~ Study, nrow = 1) +
  facet_grid(CellLineage ~ Study, scales = "free_x") +
  theme_bw() +
  scale_y_continuous(limits = c(-2, 3)) +
  stat_compare_means(comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(2.5, 2), p.adjust.method = "fdr") +
  theme(strip.text = element_text(face="bold", size=16)) +
  # theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16, angle = 0), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Percentage of shared clones in CD34 cells"), 
       subtitle = "Filters: purity, PB and BM separated, markers separated. Months under analysis: 24-60",
       x = "Lineage", 
       y = "Sharing percentage", 
       colour = "Cell Marker", fill = "Cell Marker")

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg_median.byLineage.p2.v2.pdf", sep = ""), height=7, width=3)
plot(p_lines_mbte_zscore_stable_min3_stats_age_median_p2)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.v1.wrap_trials.zscore.zoom_stats_avg_median.byLineage.p2.v2.png", sep = ""), height=7, width=3, units = "in", res = 300)
plot(p_lines_mbte_zscore_stable_min3_stats_age_median_p2)
dev.off()



# --- test: do median from percentages -----
allpatients_shared34_stats_supergroup_noCD34_tmp <- allpatients_shared34_stats_supergroup_noCD34
allpatients_shared34_stats_supergroup_noCD34_tmp_stat <- allpatients_shared34_stats_supergroup_noCD34_tmp %>% group_by(StudyID, FollowUp, CellType) %>% summarise(Avg_PercNIS_onOverallSharedCD34BM = mean(PercNIS_onOverallSharedCD34BM, na.rm = TRUE), Median_PercNIS_onOverallSharedCD34BM = median(PercNIS_onOverallSharedCD34BM, na.rm = T))
write.xlsx(x = allpatients_shared34_stats_supergroup_noCD34_tmp_stat, file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.stats_avg_median_pertrial.xlsx", sep = ""), sheetName = "Avg Median per study")
allpatients_shared34_stats_supergroup_noCD34_tmp_stat$Study <- factor(allpatients_shared34_stats_supergroup_noCD34_tmp_stat$StudyID, levels = study_list)

scale_color_manual_colors_sortedbyname <- colors_lineages[1,]
p_lines_mbte_avg <- ggplot(data = allpatients_shared34_stats_supergroup_noCD34_tmp_stat, aes(x = FollowUp, y = Avg_PercNIS_onOverallSharedCD34BM, color = CellType), na.rm = T, se = TRUE) +
  # p <- ggplot(data = allpatients_shared34_stats_supergroup_pediatriconly_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = HematoLineage), na.rm = T, se = TRUE) +
  # geom_point(size=4, alpha = .4) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(aes(shape = CellMarker), size=4, alpha = .6) +
  # geom_line(size=3, alpha = .3) +
  # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
  # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
  # stat_smooth(level = 2, size = 2) +
  # stat_smooth(aes(fill = CellType), level = 0.7, size = 2) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  # stat_smooth(aes(fill = CellType), method = lm, formula = y ~ splines::bs(x, 3), size = 2, level = 0.4, se = T) +
  # facet_wrap( ~ SubjectID, ncol = 4) +
  facet_grid(. ~ Study, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_supergroup_noCD34_tmp_stat$FollowUp, na.rm = T), 6) ) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 60) + # Add global p-value
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34.")

allpatients_shared34_stats_supergroup_noCD34_tmp_stat_scaled <- allpatients_shared34_stats_supergroup_noCD34_tmp_stat %>% group_by(StudyID, FollowUp) %>% mutate(sum_Avg_PercNIS_onOverallSharedCD34BM = sum(Avg_PercNIS_onOverallSharedCD34BM, na.rm = TRUE), newPerc_Avg_PercNIS_onOverallSharedCD34BM = (Avg_PercNIS_onOverallSharedCD34BM/sum_Avg_PercNIS_onOverallSharedCD34BM))

# p_lines_mbte_avg_scaled <- ggplot(data = allpatients_shared34_stats_supergroup_noCD34_tmp_stat_scaled[which(allpatients_shared34_stats_supergroup_noCD34_tmp_stat_scaled$FollowUp < 90 & !(allpatients_shared34_stats_supergroup_noCD34_tmp_stat_scaled$StudyID == "BTHAL" & allpatients_shared34_stats_supergroup_noCD34_tmp_stat_scaled$FollowUp == 18)),], aes(x = FollowUp, y = newPerc_Avg_PercNIS_onOverallSharedCD34BM, color = CellType, fill = CellType), na.rm = T, se = TRUE) +
p_lines_mbte_avg_scaled <- ggplot(data = allpatients_shared34_stats_supergroup_noCD34_tmp_stat_scaled, aes(x = FollowUp, y = newPerc_Avg_PercNIS_onOverallSharedCD34BM, color = CellType, fill = CellType), na.rm = T, se = TRUE) +
  # p <- ggplot(data = allpatients_shared34_stats_supergroup_pediatriconly_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = HematoLineage), na.rm = T, se = TRUE) +
  # geom_point(size=4, alpha = .4) +
  # geom_bar(stat="identity") +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(aes(shape = CellMarker), size=4, alpha = .6) +
  # geom_line(size=3, alpha = .3) +
  # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
  # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
  # stat_smooth(level = 2, size = 2) +
  # stat_smooth(aes(fill = CellType), level = 0.4, size = 2) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  # stat_smooth(aes(fill = CellType), method = lm, formula = y ~ splines::bs(x, 3), size = 2, level = 0.4, se = T) +
  # facet_wrap( ~ SubjectID, ncol = 4) +
  facet_grid(. ~ Study, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent) + 
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_supergroup_noCD34_tmp_stat_scaled$FollowUp, na.rm = T), 12) ) +
  # scale_x_continuous(limits = c(0, 84) ) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 60) + # Add global p-value
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34. Percentages re-scaled to 100%.")

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.avg_scaled.nopoints.v2.pdf", sep = ""), 
    height=4.5, width=10)
print(p_lines_mbte_avg_scaled)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.avg_scaled.nopoints.v2.png", sep = ""), 
    height=4.5, width=10, units = "in", res = 300)
print(p_lines_mbte_avg_scaled)
dev.off()

p_lines_mbte_avg_scaled_identicalFU <- ggplot(data = allpatients_shared34_stats_supergroup_noCD34_tmp_stat_scaled, aes(x = FollowUp, y = newPerc_Avg_PercNIS_onOverallSharedCD34BM, color = CellType, fill = CellType), na.rm = T, se = TRUE) +
  # p <- ggplot(data = allpatients_shared34_stats_supergroup_pediatriconly_noCD34, aes(x = FollowUp, y = PercNIS_onOverallSharedCD34BM, color = HematoLineage), na.rm = T, se = TRUE) +
  # geom_point(size=4, alpha = .4) +
  # geom_bar(stat="identity") +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # geom_point(aes(shape = CellMarker), size=4, alpha = .6) +
  # geom_line(size=3, alpha = .3) +
  # geom_smooth(aes(color = HematoLineage, fill = HematoLineage)) +
  # stat_smooth(aes(shape = CellMarker), level = 2, size = 2) +
  # stat_smooth(level = 2, size = 2) +
  # stat_smooth(aes(fill = CellType), level = 0.4, size = 2) +
  stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2) +
  # stat_smooth(aes(fill = CellType), method = lm, formula = y ~ splines::bs(x, 3), size = 2, level = 0.4, se = T) +
  # facet_wrap( ~ SubjectID, ncol = 4) +
  facet_grid(. ~ Study, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent) + 
  scale_x_continuous(breaks = seq(0, max(allpatients_shared34_stats_supergroup_noCD34_tmp_stat_scaled$FollowUp, na.rm = T), 12), limits = c(0, followup_limit_crosstrial) ) +
  # scale_x_continuous(limits = c(0, 84) ) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(40, 50, 30)) + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 60) + # Add global p-value
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 16, colour = "darkred", angle = 270, face="bold")) +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=18), plot.title = element_text(size=20)) +
  labs(title = paste("Profile of CD34 BM output towards Myeloid/B/T/Erythroid lineages"), x = "FollowUp months after GT", y = "Perc. of IS shared with CD34 BM", 
       subtitle = "All clinical trials included, all patients included. WAS clinical trial includes early progenitors as supergroup of CD34. Percentages re-scaled to 100%.")

pdf(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.avg_scaled.nopoints.v2.60m.pdf", sep = ""), 
    height=4.5, width=10)
print(p_lines_mbte_avg_scaled_identicalFU)
dev.off()
png(file = paste0("analyses/cd34output/", analysis_folder_date, "/", paste(study_list, collapse = "_"), ".CD34output.allpatients.avg_scaled.nopoints.v2.60m.png", sep = ""), 
    height=4.5, width=10, units = "in", res = 300)
print(p_lines_mbte_avg_scaled_identicalFU)
dev.off()




###############################################################
## HSPC estimate
###############################################################
markerlist <- c("CD13", "CD14", "CD15", "CD19", "CD3", "CD34", "GLY", "CD36", "GLYA")
hspce_base_folder <- paste0("analyses/hspc_size/")
dir.create(file.path(getwd(), hspce_base_folder), showWarnings = FALSE)
dir.create(file.path(hspce_base_folder, analysis_folder_date), showWarnings = FALSE)

mld_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial MLD/"
bthal_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial bThalassemia/"
was_wd <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial WAS/"

## ----- per patient ------
mld_hspc <- read.csv(file = paste0(mld_wd, "analyses/11.stem_population/", analysis_folder_date, "/20220113_MLD_HSC_population_size_estimate.tsv"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
mld_hspc$ClinicalTrial <- "MLD"
mld_hspc <- mld_hspc[which(mld_hspc$SubjectID %in% grep("MLD", levels(factor(mld_hspc$SubjectID)), value = T)),]

was_hspc <- read.csv(file = paste0(was_wd, "analyses/11.stem_population/", analysis_folder_date, "/20220113_WAS_HSC_population_size_estimate.tsv"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
was_hspc$ClinicalTrial <- "WAS"
was_hspc <- was_hspc[which(was_hspc$SubjectID %in% grep("WAS", levels(factor(was_hspc$SubjectID)), value = T)),]

bthal_hspc <- read.csv(file = paste0(bthal_wd, "analyses/09.stem_population/", analysis_folder_date, "/20220113_BTHAL_HSC_population_size_estimate.tsv"), header=T, fill=T, check.names = FALSE, sep = '\t', na.strings = c("NA", ''))
bthal_hspc$ClinicalTrial <- "BTHAL"
bthal_hspc <- bthal_hspc[which(bthal_hspc$SubjectID %in% grep("BTHAL", levels(factor(bthal_hspc$SubjectID)), value = T)),]

# merge
common_labels <- Reduce(intersect, list(colnames(was_hspc), colnames(mld_hspc), colnames(bthal_hspc)))
all_results_hspc <- Reduce(rbind, list(mld_hspc[common_labels], was_hspc[common_labels], bthal_hspc[common_labels]))

all_results_hspc$Study <- factor(all_results_hspc$ClinicalTrial, levels = study_list)
all_results_hspc[is.na(all_results_hspc)] <- NA
# all_results_hspc_max <- all_results_hspc %>% group_by(SubjectID, Study) %>% summarise(count = max(abundance), TimePoint = min(TimePoint_to))

# compute avg VCN in vivo in myeloid cells using the df all_results_h (see diversity section)
all_results_h_avgVCNmyeloid12m <- all_results_h %>% filter(TimepointMonths >= 1, CellType == "Myeloid", Tissue == "PB") %>% group_by(SubjectID, Study, CellType, Tissue) %>% summarise(avg_invivo_VCN = mean(VCN_avg, na.rm = T), sd_invivo_VCN = sd(VCN_avg, na.rm = T), min_timepoint_invivo = min(TimepointMonths, na.rm = T), avg_invivo_Hindex = mean(Hindex, na.rm = T))

label_metadata_small <- label_metadata[c("ProjectID", "SubjectID", "VCN", "TimePoint", "Tissue", "CellMarker")]
all_results_h_avgVCNmyeloid12m <- sqldf("select  ProjectID as Study, SubjectID, AVG(VCN) as avg_invivo_VCN
      from label_metadata_small
      -- where TimePoint > 360 and CellMarker in ('CD13', 'CD14', 'CD15') and Tissue like 'PB'
      where Tissue like 'PB'
      group by ProjectID, SubjectID ")
      
# combine data of treatment by patient
all_results_hspc_extended <- merge(x = all_results_hspc, y = patients_summary, by = c("SubjectID"), all.x = T)
all_results_hspc_extended <- merge(x = all_results_hspc_extended, y = all_results_h_avgVCNmyeloid12m[c("SubjectID",  "avg_invivo_VCN")], by = c("SubjectID"), all.x = T)

# do parameter corrections by VCN
# Bushman 1: correction of CD34 by VCN od CD34
all_results_hspc_extended$CD34.pro.Kg_corrected <- ifelse(all_results_hspc_extended$VCN < 1, 
                                                          all_results_hspc_extended$CD34.pro.Kg*all_results_hspc_extended$VCN, 
                                                          all_results_hspc_extended$CD34.pro.Kg)
all_results_hspc_extended$abundance_corrected <- ifelse(all_results_hspc_extended$VCN > 1,
                                                        all_results_hspc_extended$abundance / all_results_hspc_extended$VCN,
                                                        all_results_hspc_extended$abundance)
all_results_hspc_extended$PopSize_corrected_VCNinvivo <- ifelse(all_results_hspc_extended$avg_invivo_VCN > 1,
                                                                all_results_hspc_extended$PopSize / all_results_hspc_extended$avg_invivo_VCN,
                                                                all_results_hspc_extended$PopSize)
all_results_hspc_extended$PopSize_corrected_VCNCD34 <- ifelse(all_results_hspc_extended$VCN > 1,
                                                              all_results_hspc_extended$PopSize / all_results_hspc_extended$VCN,
                                                              all_results_hspc_extended$PopSize)

# write all results
write.xlsx(x = all_results_hspc_extended, file = paste(hspce_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".HSPCsize.AllResults.ByPatient.xlsx", sep = ""), sheetName = "per patient")

# slice and plot
allstudies_hspc_slice <- all_results_hspc_extended[which((all_results_hspc_extended$Timepoints == "Consecutive") 
                                                           & (all_results_hspc_extended$Model == c("Mth Chao (LB)") )
                                                           & !is.na(all_results_hspc_extended$abundance)), ]

# few stats to report in the maunscript: avg per patient within 12 m and over
allstudies_hspc_slice_summary_lt12 <- allstudies_hspc_slice %>% filter(TimePoint_to < 12) %>% group_by(SubjectID, Study) %>% summarise(avg_PopSize_lt12m = mean(PopSize_corrected_VCNinvivo, na.rm = T), sd_PopSize_lt12m = sd(PopSize_corrected_VCNinvivo, na.rm = T))
allstudies_hspc_slice_summary_gt24 <- allstudies_hspc_slice %>% filter(TimePoint_to >= 24) %>% group_by(SubjectID, Study) %>% summarise(avg_PopSize_gt12m = mean(PopSize_corrected_VCNinvivo, na.rm = T), sd_PopSize_gt12m = sd(PopSize_corrected_VCNinvivo, na.rm = T))
write.xlsx(x = allstudies_hspc_slice_summary_lt12, file = paste(hspce_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".HSPCsize.AllResults.ByPatient.summarystats_lt12m.PopSize_corrected_VCNinvivo.xlsx", sep = ""), sheetName = "< 12m")
write.xlsx(x = allstudies_hspc_slice_summary_gt24, file = paste(hspce_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".HSPCsize.AllResults.ByPatient.summarystats_gt12m.PopSize_corrected_VCNinvivo.xlsx", sep = ""), sheetName = ">= 24m", )

# the same without correction VCN in vivo
allstudies_hspc_slice_summary_lt12 <- allstudies_hspc_slice %>% filter(TimePoint_to < 12) %>% group_by(SubjectID, Study) %>% summarise(avg_PopSize_lt12m = mean(PopSize_corrected_VCNCD34, na.rm = T), sd_PopSize_lt12m = sd(PopSize_corrected_VCNCD34, na.rm = T))
allstudies_hspc_slice_summary_gt24 <- allstudies_hspc_slice %>% filter(TimePoint_to >= 24) %>% group_by(SubjectID, Study) %>% summarise(avg_PopSize_gt12m = mean(PopSize_corrected_VCNCD34, na.rm = T), sd_PopSize_gt12m = sd(PopSize_corrected_VCNCD34, na.rm = T))
write.xlsx(x = allstudies_hspc_slice_summary_lt12, file = paste(hspce_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".HSPCsize.AllResults.ByPatient.summarystats_lt12m.PopSize_corrected_VCNCD34.xlsx", sep = ""), sheetName = "< 12m")
write.xlsx(x = allstudies_hspc_slice_summary_gt24, file = paste(hspce_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".HSPCsize.AllResults.ByPatient.summarystats_gt12m.PopSize_corrected_VCNCD34.xlsx", sep = ""), sheetName = ">= 24m", )


# plot_data_slice <- allstudies_hspc_slice
# plot_hspc_overtime <- ggplot(data = plot_data_slice, aes(x = FU_end, y = PopSize, fill = Study, color = Study), na.rm = T, se = TRUE) +
#   # geom_point(aes(fill = SubjectID, color = SubjectID), size = 3, alpha = 0.7) +
#   geom_pointrange(size = 3, alpha = 0.7) +
#   geom_line(aes(group = SubjectID, color = "gray"), size = 1, alpha = 0.2) +
#   # geom_errorbar(aes(ymin=(abundance-stderr), ymax=(abundance+stderr)), width=0.1) +
#   # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
#   geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
#   # scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
#   # scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
#   # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)) ) +
#   theme_bw() +
#   scale_y_continuous(labels = scales::comma, limits = c(0,100000)) +
#   scale_x_continuous(breaks = seq(0, max(plot_data_slice$FU_end, na.rm = T), 6), limits = c(0, 36)) +
#   facet_wrap( ~ ClinicalStudy, ncol = 3) +
#   theme(strip.text = element_text(face="bold", size=16)) +
#   theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
#   theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
#   labs(title = paste0("Estimate of population size of active and engrafted HSPC over time"), x = "Months after gene therapy", y = "HSPC estimate", colour = "Patient ID", fill = "Patient ID"
#        )

plot_hspc_overtime_by_patient <-
  # ggplot(data = allstudies_hspc_slice, aes(x = TimePoint_from, y = PopSize, fill = Study, color = Study), na.rm = T, se = TRUE) +
  ggplot(data = allstudies_hspc_slice, aes(x = TimePoint_from, y = PopSize, fill = Study, color = Study), na.rm = T, se = TRUE) +
  scale_color_manual(values = trials_colors) +
  scale_fill_manual(values = trials_colors) +
  # geom_point(alpha = 0.2, size = 2) +
  # geom_jitter(alpha = 0.8, size = 4) +
  geom_line(aes(group = SubjectID), size=1.5, alpha = .2, colour = "gray") +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  # geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  geom_smooth(method = "loess", formula = y ~ -log(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  scale_x_continuous(breaks = seq(0, max(allstudies_hspc_slice$TimePoint_to, na.rm = T), 6) ) +
  facet_grid( . ~ Study, scales = "free_x", space = "free") +
  theme_bw() +
  # geom_label_repel(data = subset(all_results_cumulative_max, count>100000),
  #                  aes(label = SubjectID),
  #                  color = 'white',
  #                  box.padding = unit(0.35, "lines"),
  #                  point.padding = unit(0.4, "lines"),
  #                  segment.color = 'black') +
  # scale_y_continuous(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  # scale_y_break(c(500000, 800000), scales = 0.1, ticklabels = c(60000, 90000, 120000)) + #scale_y_break(c(50000, 120000), scales=0.2) +
  # scale_y_cut(breaks=c(50000, 80000), which=c(1, 3), scales=c(0.5, 0.5, 4)) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  # facet_grid(ProjectID ~ ., scales = "free_y", space = "free") +
  # facet_grid(ProjectID ~ ., scales = "free_y") +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Estimated number of active and engrafted HSPC"), 
       subtitle = paste0("Filters: only PB Myeloid cells, with an overall number of reads >=3."),
       x = "Months after gene therapy", 
       y = "Estimated HSPC", 
       colour = "Study", fill = "Study")

pdf(file = paste(hspce_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".HSPCsize.AllResults.ByPatient.pdf", sep = ""), height=6, width=13)
# plot(plot_cumulative_by_patient_cellmarker)
plot_hspc_overtime_by_patient
dev.off()
png(file = paste(hspce_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".HSPCsize.AllResults.ByPatient.png", sep = ""), height=6, width=13, units = "in", res = 300)
plot_hspc_overtime_by_patient
dev.off()

plot_hspc_overtime_by_patient_larger <-
  # ggplot(data = allstudies_hspc_slice, aes(x = TimePoint_from, y = PopSize, fill = Study, color = Study), na.rm = T, se = TRUE) +
  ggplot(data = allstudies_hspc_slice, aes(x = TimePoint_from, y = PopSize_corrected_VCNinvivo, fill = Study, color = Study), na.rm = T, se = TRUE) +
  scale_color_manual(values = trials_colors) +
  scale_fill_manual(values = trials_colors) +
  # geom_point(alpha = 0.2, size = 2) +
  # geom_jitter(alpha = 0.8, size = 4) +
  geom_line(aes(group = SubjectID), size=1.5, alpha = .2, colour = "gray") +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.7) +
  # geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  geom_smooth(method = "loess", formula = y ~ -log(x), stat = "smooth", position = "identity", alpha = 0.6, level = 0.75) +
  scale_x_continuous(breaks = seq(0, max(allstudies_hspc_slice$TimePoint_to, na.rm = T), 12) ) +
  facet_grid( . ~ Study, scales = "free_x", space = "free") +
  theme_bw() +
  # geom_label_repel(data = subset(all_results_cumulative_max, count>100000),
  #                  aes(label = SubjectID),
  #                  color = 'white',
  #                  box.padding = unit(0.35, "lines"),
  #                  point.padding = unit(0.4, "lines"),
  #                  segment.color = 'black') +
  # scale_y_continuous(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  # scale_y_break(c(500000, 800000), scales = 0.1, ticklabels = c(60000, 90000, 120000)) + #scale_y_break(c(50000, 120000), scales=0.2) +
  # scale_y_cut(breaks=c(50000, 80000), which=c(1, 3), scales=c(0.5, 0.5, 4)) +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(7)) +
  # stat_compare_means(label.y = 6) +
  # facet_grid(ProjectID ~ ., scales = "free_y", space = "free") +
  # facet_grid(ProjectID ~ ., scales = "free_y") +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Estimated number of active and engrafted HSPC"), 
       subtitle = paste0("Filters: only PB Myeloid cells, with an overall number of reads >=3."),
       x = "Months after gene therapy", 
       y = "Estimated HSPC", 
       colour = "Study", fill = "Study")

pdf(file = paste(hspce_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".HSPCsize.AllResults.ByPatient.correctedVCNinvivo.larger.pdf", sep = ""), height=5, width=8)
# plot(plot_cumulative_by_patient_cellmarker)
plot_hspc_overtime_by_patient_larger
dev.off()
png(file = paste(hspce_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".HSPCsize.AllResults.ByPatient.correctedVCNinvivo.larger.png", sep = ""), height=5, width=8, units = "in", res = 300)
plot_hspc_overtime_by_patient_larger
dev.off()



# slice and plot
allstudies_hspc_slice_global <- all_results_hspc_extended[which((all_results_hspc_extended$Timepoints == "Stable") 
                                                         & (all_results_hspc_extended$ModelSetUp == c("mthchaobc") )
                                                         & !is.na(all_results_hspc_extended$abundance)), ]
# write all results
write.xlsx(x = allstudies_hspc_slice_global, file = paste(hspce_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".HSPCsize.Stable.ByPatient.xlsx", sep = ""), sheetName = "per patient")

model_popsize <- lm(PopSize_corrected_VCNCD34 ~ CD34.pro.Kg_corrected, data = allstudies_hspc_slice_global)
summary(model_popsize)

plot_hspc_global_by_patient <-
  ggplot(data = allstudies_hspc_slice_global, aes(x = CD34.pro.Kg_corrected, y = PopSize_corrected_VCNCD34), na.rm = T, se = TRUE) +
  scale_color_manual(values = trials_colors) +
  scale_fill_manual(values = trials_colors) +
  geom_pointrange(aes(fill = Study, color = Study, ymin=PopSize_corrected_VCNCD34-stderr, ymax=PopSize_corrected_VCNCD34+stderr)) +
  # geom_pointrange(aes(ymin=abundance_corrected-stderr, ymax=abundance_corrected+stderr)) +
  # facet_grid( . ~ Study, scales = "free_x", space = "free") +
  theme_bw() +
  scale_y_log10(labels = scales::comma) +
  geom_smooth(method = "lm", stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  # geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  # theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
  labs(title = paste0("Size of active and engrafted HSPC"), 
       subtitle = paste0("Filters: only PB Myeloid cells, with an overall number of reads and cells >=3.\nInfused dose of CD34 pro kg corrected by VCN of CD34 cells."),
       x = "Dose infused", 
       y = "Estimated HSPC", 
       colour = "Study", fill = "Study")

pdf(file = paste(hspce_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".HSPCsize.AllResults.Global.pdf", sep = ""), height=5, width=5)
# plot(plot_cumulative_by_patient_cellmarker)
plot_hspc_global_by_patient
dev.off()
png(file = paste(hspce_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".HSPCsize.AllResults.Global.png", sep = ""), height=6, width=6, units = "in", res = 300)
plot_hspc_global_by_patient
dev.off()






# fai plot per scale free
plot_hspc_overtime_allfu_normVCN <- ggplot(data = plot_data_slice_full, aes(x = FU_start, y = PopSize_ScaledVNC), na.rm = T, se = TRUE) +
  geom_point(aes(fill = SubjectID, color = SubjectID), size = 3, alpha = 0.7) +
  geom_line(aes(fill = SubjectID, color = SubjectID), size = 1, alpha = 0.2) + 
  # geom_errorbar(aes(ymin=(abundance-stderr), ymax=(abundance+stderr)), width=0.1) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  # scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)) ) +
  theme_bw() +
  # scale_y_continuous(labels = scales::comma, limits = c(0,100000)) +
  scale_x_continuous(breaks = seq(0, max(plot_data_slice$FU_start, na.rm = T), 6) ) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Estimate of population size of active and engrafted HSPC over time"), 
       subtitle = "Normalized on in vivo VCN (computed as average on stable tp and myeloid cells).",
       x = "Months after gene therapy", y = "HSPC estimate", 
       colour = "Patient ID", fill = "Patient ID"
  )

pdf(file = paste0("analyses/hspc_size/202004/", paste(study_list, collapse = "_"), ".HSPC_size_overtime.allFU.normVCN.pdf", sep = ""), height=7.5, width=12)
print(plot_hspc_overtime_allfu_normVCN)
dev.off()
png(file = paste0("analyses/hspc_size/202004/", paste(study_list, collapse = "_"), ".HSPC_size_overtime.allFU.normVCN.png", sep = ""), height=7.5, width=12, units = "in", res = 300)
print(plot_hspc_overtime_allfu_normVCN)
dev.off()



### try correlations

library(umap)
allstudies_hspc_slice_global_edited <- read.xlsx(xlsxFile = paste(hspce_base_folder, "/", analysis_folder_date, "/", analysis_folder_date, ".HSPCsize.Stable.ByPatient.addedVCNWAS.xlsx", sep = ""))
# slice_cols_pca <- c("H_inVivo", "HSPC_estimate", "avg_VCN_MyeloPB", "n..IS", "weight", "Vol", "CD34.10e6", "TransductionEfficiency", "VCN", "TreatmentAge.years")
slice_cols_pca <- c("avg_invivo_Hindex", "PopSize_corrected_VCNCD34", "TransductionEfficiency", "avg_invivo_VCN", "TreatmentAge.years", "CD34.pro.Kg_corrected", "SubjectID")
# CD34.pro.Kg_corrected, y = PopSize_corrected_VCNCD34
slice_data_topca <- allstudies_hspc_slice_global_edited[slice_cols_pca]
slice_data_topca_filtered <- slice_data_topca[which(!is.na(slice_data_topca$PopSize_corrected_VCNCD34) & !is.na(slice_data_topca$TransductionEfficiency) & !is.na(slice_data_topca$avg_invivo_Hindex)),]
rownames(slice_data_topca_filtered) <- slice_data_topca_filtered$SubjectID
slice_data_topca_filtered_umap <- umap(slice_data_topca_filtered[setdiff(colnames(slice_data_topca_filtered), "SubjectID")])

umap_pts <- as.data.frame(slice_data_topca_filtered_umap$layout)
umap_pts$SubjectID <- rownames(umap_pts)
# umap_pts_extended <- merge(x = umap_pts, y = allstudies_hspc_stable_extended, by = c("SubjectID"), all.x = T)
umap_pts_extended <- merge(x = umap_pts, y = allstudies_hspc_slice_global_edited, by = c("SubjectID"), all.x = T)

umap_pts_extended_t1 <-
  ggplot(data = umap_pts_extended, aes(x = V1, y = V2, color = Disease, fill = Disease), na.rm = T, se = TRUE) + 
  geom_point(alpha = .9) +
  # geom_errorbar(aes(ymin=(Avg_VCN_from12m_unsorted-Stdev_VCN_from12m_unsorted), ymax=(Avg_VCN_from12m_unsorted+Stdev_VCN_from12m_unsorted)), width=0.1) +
  # geom_line(size = 2, alpha = .7) +
  # facet_grid( PatientGroup ~ Transfusion ) +
  # scale_x_continuous(breaks = seq(0, max(stats$TimePoint), 6), limits = c(0, (max(stats$TimePoint)) ) ) +
  # scale_y_log10(labels = scales::comma) +
  theme_bw() +
  geom_label_repel(
    data = umap_pts_extended,
    # data = subset(all_results_cis, !is.na(Onco1_TS2) & minus_log2_integration_freq_withtolerance > 2.8),
    aes(label = SubjectID),
    color = 'white',
    size = 3.5,
    segment.color = 'black') +
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  theme(strip.text.x = element_text(size = 14, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 14, colour = "darkred", angle = 270, face="bold")) +
  labs(title = paste("UMAP PCA data integration"), x = "UMAP 1", y = "UMAP 2", 
            color = "Disease", fill = "Disease",
            subtitle = "In vivo data, from 12 months. Integrating molecular and treatment data.")

umap_pts_extended_t2 <-
  ggplot(data = umap_pts_extended, aes(x = V1, y = V2, color = GT.type, fill = GT.type), na.rm = T, se = TRUE) + 
  geom_point(alpha = .9) +
  # geom_errorbar(aes(ymin=(Avg_VCN_from12m_unsorted-Stdev_VCN_from12m_unsorted), ymax=(Avg_VCN_from12m_unsorted+Stdev_VCN_from12m_unsorted)), width=0.1) +
  # geom_line(size = 2, alpha = .7) +
  # facet_grid( PatientGroup ~ Transfusion ) +
  # scale_x_continuous(breaks = seq(0, max(stats$TimePoint), 6), limits = c(0, (max(stats$TimePoint)) ) ) +
  # scale_y_log10(labels = scales::comma) +
  theme_bw() +
  geom_label_repel(
    data = umap_pts_extended,
    # data = subset(all_results_cis, !is.na(Onco1_TS2) & minus_log2_integration_freq_withtolerance > 2.8),
    aes(label = SubjectID),
    color = 'white',
    size = 3.5,
    segment.color = 'black') +
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  theme(strip.text.x = element_text(size = 14, colour = "darkblue", angle = 0, face="bold"), strip.text.y = element_text(size = 14, colour = "darkred", angle = 270, face="bold")) +
  labs(title = paste("UMAP PCA data integration"), x = "UMAP 1", y = "UMAP 2", 
       color = "GT type", fill = "GT type",
       subtitle = "In vivo data, from 12 months. Integrating molecular and treatment data.")

# write.xlsx(x = umap_pts_extended, file = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".UMAP.t1.xlsx", sep = "") )
pdf(file = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".UMAP.t2.pdf", sep = ""), height=6, width=8)
print(umap_pts_extended_t1)
dev.off()
png(file = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".UMAP.t2.png", sep = ""), height=6, width=8, units = "in", res = 300)
print(umap_pts_extended_t1)
dev.off()
png(file = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".UMAP.t2.v2.png", sep = ""), height=7, width=9, units = "in", res = 300)
print(umap_pts_extended_t2)
dev.off()


my_data <- slice_data_topca_filtered[setdiff(colnames(slice_data_topca_filtered), "SubjectID")]
write.xlsx(x = slice_data_topca_filtered, file = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".my_data.full.xlsx", sep = ""))
# then editing patient names by hand and saving file ".colors". read it again
my_data <- read.xlsx(xlsxFile = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".my_data.colors.xlsx", sep = ""))
rownames(my_data) <- my_data$SubjectID
my_data <- my_data[c(-1,-2)]
names(my_data) <- c("H index", "HSPC size", "Tr.Efficiency", "VCN in vivo", "Age", "CD34+ 10e6/Kg")
library("PerformanceAnalytics")
chart.Correlation(my_data, histogram=TRUE, pch=3, method = c("pearson"))
# library(corrplot)
# res2 <- rcorr(as.matrix(my_data))
# corrplot(my_data, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)

# -----------------------
# try to extract the numbers below , first test on example in https://stackoverflow.com/questions/55589048/is-there-an-r-function-to-export-the-correlations-showed-in-a-correlation-matrix 
data_list <- split(s1, s1$media)
p_value <- lapply(data_list, function(x) corrplot::cor.mtest(x[, 4:12])[["p"]])
correlation <- lapply(data_list, function(x) cor(x[, 4:12], method = "spearman"))
# now my data
my_data <- read.xlsx(xlsxFile = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".my_data.colors.paperEdits.xlsx", sep = ""))
my_data$CTC <- "CTC"
data_list <- split(my_data, my_data$CTC)
# p_value <- lapply(data_list, function(x) corrplot::cor.mtest(x[, 3:8], conf.level = 0.95, adjust="bonferroni", alpha=.05)[["p"]])
p_value <- lapply(data_list, function(x) corrplot::cor.mtest(x[, 3:8], adjust="none", conf.level = 0.95, alpha=.05)[["p"]])
p_value <- lapply(data_list, function(x) corrplot::cor.mtest(x[, 3:8], adjust="none", alpha=.05)[["p"]])
full_correlations <- lapply(data_list, function(x) corrplot::cor.mtest(x[, 3:8], adjust="none", alpha=.05) )
correlation <- lapply(data_list, function(x) cor(x[, 3:8], method = "pearson"))
p_value_corrected_bonferroni <- as.data.frame(p_value) * 15 # since 15 are the combinations of tests
names(p_value_corrected_bonferroni) <- c("H index", "HSPC size", "Tr.Efficiency", "VCN in vivo", "Age", "CD34 proKg")
p_value_corrected_bonferroni$Variable <- colnames(p_value_corrected_bonferroni)
write.xlsx(x = p_value_corrected_bonferroni, file = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".my_data.colors.paperEdits.Pvals.Bonferroni.xlsx", sep = "") )
# -----------------------


pdf(file = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".Corrlations.t2.pdf", sep = ""), height=8, width=8)
correlations_t1 <- chart.Correlation(my_data, histogram=TRUE, pch=19, method = c("pearson"))
dev.off()
png(file = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".Corrlations.t2.png", sep = ""), height=8, width=8, units = "in", res = 300)
correlations_t1 <- chart.Correlation(my_data, histogram=TRUE, pch=19, method = c("pearson"))
dev.off()

# library("Hmisc")
# res2 <- rcorr(as.matrix(my_data))
# plot(res2)

library(factoextra)
res.pca <- prcomp(my_data, scale = TRUE)
# res.pca <- prcomp(my_data, scale = F)

pca_residuals <- fviz_eig(res.pca)
pca_cos_dist <- fviz_pca_ind(res.pca,
             title = "PCA - Patients",
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
pca_contrib_dist <- fviz_pca_var(res.pca,
             title = "PCA - Variables",
             col.var = "contrib", # Color by contributions to the PC
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             gradient.cols = c("orange", "darkred"),
             repel = TRUE     # Avoid text overlapping
)
pca_biplot_eig <- fviz_pca_biplot(res.pca, 
                title = "PCA - Biplot",
                repel = TRUE,
                max.overlaps = 5,
                col.var = "darkred", # Variables color
                col.ind = "gray40"  # Individuals color
)

pdf(file = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".PCA.t2.newnames.pdf", sep = ""), height=8, width=8)
pca_biplot_eig
pca_cos_dist
pca_contrib_dist
dev.off()

slice_data_topca_filtered_extended <- merge(x = slice_data_topca_filtered, y = allstudies_hspc_slice_global_edited[c("SubjectID", "GT.type")], by = c("SubjectID"), all.x = T)
write.xlsx(x = slice_data_topca_filtered_extended, file = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".slice_data_topca_filtered_extended.xlsx", sep = ""))
# then editing patient names by hand and saving file back and read it again
# slice_data_topca_filtered_extended <- read.xlsx(xlsxFile = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".slice_data_topca_filtered_extended.xlsx", sep = ""))
slice_data_topca_filtered_extended <- read.xlsx(xlsxFile = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".slice_data_topca_filtered_extended.newcodes.xlsx", sep = ""))
rownames(my_data) <- my_data$SubjectID
# groups <- as.factor(decathlon2$Competition[1:23])
groups <- as.factor(slice_data_topca_filtered_extended$GT.type)
pca_bp_gr1 <- fviz_pca_ind(res.pca,
             title = "PCA - Biplot with clustering",
             col.ind = groups, # color by groups
             palette = c("orange",  "forestgreen", "lightblue"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)
pca_bp_gr2 <- fviz_pca_biplot(res.pca, 
                title = "PCA - Biplot with clustering",
                max.overlaps = 5,
                # col.var = "darkred", # Variables color
                col.ind = groups, # color by groups
                palette = c("orange",  "forestgreen", "mediumpurple"),
                addEllipses = TRUE, # Concentration ellipses
                ellipse.type = "confidence",
                legend.title = "GT source",
                repel = TRUE
)

pdf(file = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".PCA.t2.groups.newnames.pdf", sep = ""), height=7, width=8)
pca_bp_gr2
pca_bp_gr1
dev.off()

pdf(file = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".PCA.t2.groups.newcodes.pdf", sep = ""), height=7, width=8)
pca_bp_gr2
pca_bp_gr1
dev.off()

png(file = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".Corrlations.t2.png", sep = ""), height=8, width=8, units = "in", res = 300)
pca_bp_gr1
dev.off()

# try to add the disease data in the plot
slice_data_topca_filtered_extended2 <- read.xlsx(xlsxFile = paste("analyses/pca/", analysis_folder_date, "/", analysis_folder_date, ".", ProjectID, ".slice_data_topca_filtered_extended.newcodes2.xlsx", sep = ""))
groups <- as.factor(slice_data_topca_filtered_extended2$GT.type)
diseaseconditions <- as.factor(slice_data_topca_filtered_extended2$Disease)
pca_bp_gr3 <- fviz_pca_biplot(res.pca, 
                              title = "PCA - Biplot with clustering",
                              max.overlaps = 5,
                              # geom = "point",
                              # col.var = "darkred", # Variables color
                              col.ind = groups, # color by groups
                              palette = c("orange",  "forestgreen", "mediumpurple"),
                              addEllipses = TRUE, # Concentration ellipses
                              ellipse.type = "confidence",
                              legend.title = "GT source",
                              repel = TRUE
) + geom_point(aes(shape = diseaseconditions)) + labs(shape = "Disease")


###############################################################
## Multi-to-Uni PLASTICITY COMMITTMENT 
###############################################################
# sharing data for each lineage
# add all trials data
global_allstudies_profile_skewing_fulldata_notNA_stats <- NULL
global_patient_skewing_summary_filtered <- NULL

#####------ BM ---------------------
faceting_order <- c("MLD", "WAS", "BTHAL")
lineage_faceting_order <- c("Multilineage", "Erythroid", "Myeloid", "B", "T")
trials_colors <- c("darkblue", "forestgreen", "firebrick")
prefix_trials_outname <- paste(faceting_order, collapse = "_")
results_folder_name <- "multi_uni_sharing/"
this_run_foldername <- "MyBTEry34_allMarkers_noLy_BM_SC3/" # "MyBTEry34_strictMarkers_fromwhole_noLy_AllMarkers/" # "MyBTEry34/" "MyBTEry34_ISgt6/" "MyBTEry34_BMonly/" "MyBTEry34_PBonly/"
# this_run_foldername_in_was <- "MyBTEry34_strictMarkers/" # "MyBTEry34/" "MyBTEry34_ISgt6/" "MyBTEry34_BMonly/" "MyBTEry34_PBonly/"
threshold_abundance <- 1 # < threshold_abundance will be removed
out_xlsxfile_infix <- "gdf_hlfu_filtered"
out_file_infix <- ".hlfu_34MyBTEry"
prefix_dest_folder <- "AllPts_"
suffix_plotfile <- "" # ".noTransfDep"
followup_limit_crosstrial <- 60 # plot up to this month after GT
excel_bagcases_foldername <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial MLD/analyses/13.hsc_profiles/"
excel_bagcases_filename <- "MLD.hlfu_34MyBTEry_noLy_no34exclusive.flag_matrix.bagcases.xlsx" # "MLD.hlfu_34MyBTEry.flag_matrix.bagcases.xlsx"  # "MLD.hlfu_34PBMyBT.flag_matrix.bagcases.xlsx" # MLD.hlfu_34MyBTEry.flag_matrix.bagcases.xlsx MLD.hlfu_34MyBTEry_noLy.flag_matrix.bagcases.xlsx MLD.hlfu_34PBMyBT_noLy.flag_matrix.bagcases.xlsx
# patients_to_exclude <- c("BTHAL001", "BTHAL003", "BTHAL004", "BTHAL006", "BTHAL010")
patients_to_exclude <- c()
name_list_trials <- prefix_trials_outname
# patients_to_exclude <- c()

color_schema_df <- read.xlsx(xlsxFile = paste(excel_bagcases_foldername, excel_bagcases_filename, sep = ""), sheet = "color_schema")
color_schema_df_sortedbyname <- color_schema_df[order(color_schema_df$ClassName),]
color_schema_df_sortedbyplotlabel <- color_schema_df[order(color_schema_df$Label),]

# add all trials data
mld_allpatients_flag_matrix_relabeled_correctedDC_freqs_melt_merge <- read.csv(file = paste0("/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial MLD/analyses/13.hsc_profiles/202112/", this_run_foldername, "MLD.AllPatients.hlfu_excl34MyBTEry.flag_matrix_relabeled_correctedDC_freqs_melt_merge.tsv"), header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", ""))
was_allpatients_flag_matrix_relabeled_correctedDC_freqs_melt_merge <- read.csv(file = paste0("/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial WAS/analyses/14.hsc_profiles/202112/", this_run_foldername, "WAS.AllPatients.hlfu_34MyBTEry.flag_matrix_relabeled_correctedDC_freqs_melt_merge.tsv"), header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", ""))
# bthal_allpatients_flag_matrix_relabeled_correctedDC_freqs_melt_merge <- read.csv(file = paste0("/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial bThalassemia/analyses/14.hsc_profiles/202112/", this_run_foldername, "BTHAL.AllPatients.hlfu_34MyBTEry.flag_matrix_relabeled_correctedDC_freqs_melt_merge.tsv"), header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", ""))
bthal_allpatients_flag_matrix_relabeled_correctedDC_freqs_melt_merge <- read.csv(file = paste0("/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial bThalassemia/analyses/14.hsc_profiles/202112/", this_run_foldername, "BTHAL.AllPatients.hlfu_34MyBTEry.flag_matrix_relabeled_correctedDC_freqs_melt_merge.tsv"), header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", ""))
# bthal_allpatients_flag_matrix_relabeled_correctedDC_freqs_melt_merge <- read.csv(file = paste0("/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial bThalassemia/analyses/14.hsc_profiles/202112/MyBTEry34_allMarkers_noLy_BM_SC3/BTHAL.AllPatients.hlfu_34MyBTEry.flag_matrix_relabeled_correctedDC_freqs_melt_merge.tsv"), header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", ""))

allstudies_profile_skewing_fulldata <- Reduce(rbind, list(mld_allpatients_flag_matrix_relabeled_correctedDC_freqs_melt_merge, was_allpatients_flag_matrix_relabeled_correctedDC_freqs_melt_merge, bthal_allpatients_flag_matrix_relabeled_correctedDC_freqs_melt_merge))
allstudies_profile_skewing_fulldata$ClassName <- allstudies_profile_skewing_fulldata$Class
allstudies_profile_skewing_fulldata <- merge(x = allstudies_profile_skewing_fulldata, y = color_schema_df[c("ClassName", "colorcode", "ColorSchema", "PlotLabel")], by = c("ClassName"))
allstudies_profile_skewing_fulldata$ClinicalStudy <- factor(allstudies_profile_skewing_fulldata$ProjectID, 
                                                            levels = faceting_order)
allstudies_profile_skewing_fulldata$PlotLabel <- factor(allstudies_profile_skewing_fulldata$PlotLabel, 
                                                        levels = lineage_faceting_order)

allstudies_profile_skewing_fulldata_notNA <- allstudies_profile_skewing_fulldata[which(!is.na(allstudies_profile_skewing_fulldata$SonR)),]

# adjustments of values
allstudies_profile_skewing_fulldata_notNA$CSonCR <- ifelse((allstudies_profile_skewing_fulldata_notNA$CR > 0 & allstudies_profile_skewing_fulldata_notNA$CU > 0), allstudies_profile_skewing_fulldata_notNA$CS / allstudies_profile_skewing_fulldata_notNA$CR, NA)
allstudies_profile_skewing_fulldata_notNA$CSonCU <- ifelse((allstudies_profile_skewing_fulldata_notNA$CU > 0), allstudies_profile_skewing_fulldata_notNA$CS *100 / allstudies_profile_skewing_fulldata_notNA$CU, NA)
allstudies_profile_skewing_fulldata_notNA$CRonCU <- ifelse((allstudies_profile_skewing_fulldata_notNA$CU > 0), allstudies_profile_skewing_fulldata_notNA$CR *100 / allstudies_profile_skewing_fulldata_notNA$CU, NA)
allstudies_profile_skewing_fulldata_notNA$CSonU <- ifelse((allstudies_profile_skewing_fulldata_notNA$U > 0), allstudies_profile_skewing_fulldata_notNA$CS *100 / allstudies_profile_skewing_fulldata_notNA$U, NA)
allstudies_profile_skewing_fulldata_notNA$CRonU <- ifelse((allstudies_profile_skewing_fulldata_notNA$U > 0), allstudies_profile_skewing_fulldata_notNA$CR *100 / allstudies_profile_skewing_fulldata_notNA$U, NA)
allstudies_profile_skewing_fulldata_notNA$CSonS <- ifelse((allstudies_profile_skewing_fulldata_notNA$S > 0), allstudies_profile_skewing_fulldata_notNA$CS *100 / allstudies_profile_skewing_fulldata_notNA$S, NA)
allstudies_profile_skewing_fulldata_notNA$CRonR <- ifelse((allstudies_profile_skewing_fulldata_notNA$R > 0), allstudies_profile_skewing_fulldata_notNA$CR *100 / allstudies_profile_skewing_fulldata_notNA$R, NA)

# patients_to_exclude <- c()
allstudies_profile_skewing_fulldata_notNA <- allstudies_profile_skewing_fulldata_notNA[which(!(allstudies_profile_skewing_fulldata_notNA$SubjectID %in% patients_to_exclude) ),]

allstudies_profile_skewing_fulldata_notNA_stats <- merge(x = allstudies_profile_skewing_fulldata_notNA, y = patients_summary, by = c("SubjectID"), all.x = T)
# allpatients_shared34_stats$AgeGroup <- 
allstudies_profile_skewing_fulldata_notNA_stats <- allstudies_profile_skewing_fulldata_notNA_stats %>%
  mutate(AgeGroup = case_when(
    TreatmentAge.years <= 2 ~"0-2",
    TreatmentAge.years <= 15 ~"2-15",
    TreatmentAge.years > 15 ~"30+"
  ))
age_list <- c("0-2", "2-15", "30+")
allstudies_profile_skewing_fulldata_notNA_stats$AgeGroup_Sorted <- factor(allstudies_profile_skewing_fulldata_notNA_stats$AgeGroup, levels = age_list)

# add the Tissue only here
allstudies_profile_skewing_fulldata_notNA_stats$Tissue <- "BM"
allstudies_profile_skewing_fulldata_notNA_stats <- allstudies_profile_skewing_fulldata_notNA_stats %>% group_by(ProjectID, SubjectID, Tissue, FollowUp) %>% mutate(CRonR_zscaled = scale(CRonR), CSonU_zscaled = scale(CSonU), CRonU_zscaled = scale(CRonU), CSonS_zscaled = scale(CSonS))
allstudies_profile_skewing_fulldata_notNA_stats[is.na(allstudies_profile_skewing_fulldata_notNA_stats)] <- NA
# copy to the global df
if (length(global_allstudies_profile_skewing_fulldata_notNA_stats) > 0) {
  global_allstudies_profile_skewing_fulldata_notNA_stats <- rbind(global_allstudies_profile_skewing_fulldata_notNA_stats, allstudies_profile_skewing_fulldata_notNA_stats)
} else {
  global_allstudies_profile_skewing_fulldata_notNA_stats <- allstudies_profile_skewing_fulldata_notNA_stats
}
# slice data at stability
allstudies_profile_skewing_fulldata_notNA_stats_stable <- allstudies_profile_skewing_fulldata_notNA_stats[which(allstudies_profile_skewing_fulldata_notNA_stats$FollowUp >= 24 & allstudies_profile_skewing_fulldata_notNA_stats$FollowUp <= 60),]
allstudies_profile_skewing_fulldata_notNA_stats_stable <- allstudies_profile_skewing_fulldata_notNA_stats %>% group_by(SubjectID, Tissue, Class) %>% mutate(n_obs_stt_overtime = n())
allstudies_profile_skewing_fulldata_notNA_stats_stable_min3 <- allstudies_profile_skewing_fulldata_notNA_stats %>% group_by(SubjectID, Tissue, Class) %>% mutate(n_obs_stt_overtime = n()) %>% filter(n_obs_stt_overtime >= 3)


# summary stats
patient_skewing_summary <- allstudies_profile_skewing_fulldata_notNA_stats_stable_min3 %>%
  filter(n_obs_stt_overtime >= 3) %>%
  group_by(ClinicalStudy, SubjectID, Tissue, Class, TreatmentAge.years, AgeGroup, AgeGroup_Sorted, colorcode, ColorSchema, PlotLabel) %>%
  get_summary_stats(CRonR_zscaled, type = "full") 

patient_skewing_summary_filtered <- patient_skewing_summary %>% filter(n >= 2)
# patient_skewing_summary_filtered$CellLineage <- factor(patient_skewing_summary_filtered$CellType, levels = celltype_list)
patient_skewing_summary_filtered$Study <- factor(patient_skewing_summary_filtered$ClinicalStudy, levels = study_list)
patient_skewing_summary_filtered$AgeGroup_Sorted <- factor(patient_skewing_summary_filtered$AgeGroup, levels = age_list)

write.xlsx(x = patient_skewing_summary_filtered, file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.StableFU.CRonR.xlsx"))
global_patient_skewing_summary_filtered <- patient_skewing_summary_filtered


allstudies_profile_skewing_fulldata_byt <- sqldf("select ProjectID, ClinicalStudy, SubjectID, Timepoint, FollowUp, SumOverT_OverallNIS, SumOverT_SingletoneNIS, SumOverT_SharingNIS, U, S, R, SumOverT_Singletone_OnTotal, SumOverT_Sharing_OnTotal, SumOverT_SingletonOnSharing, SonU, RonU, SonR from allstudies_profile_skewing_fulldata_notNA where R > 0 group by ProjectID, SubjectID, Timepoint, FollowUp")

scale_color_manual_colors_sortedbyname <- as.character(color_schema_df_sortedbyname[which(color_schema_df_sortedbyname$ClassName %in% allstudies_profile_skewing_fulldata_notNA$Class), "colorcode"])
plot_CRonR_global <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CRonR, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 12) ) +
  facet_wrap( ~ ClinicalStudy, ncol = 3) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance) ))

plot_CRonR_global_allFU <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CRonR, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 12) ) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ))

pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.v2.pdf"), height=7, width=12)
print(plot_CRonR_global_allFU)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.v2.png"), height=7, width=12, units = "in", res = 300)
print(plot_CRonR_global_allFU)
dev.off()


plot_CRonR_global_allFU_nopoints <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CRonR, fill = Class, color = Class), na.rm = T, se = TRUE) +
  # geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  # geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  geom_smooth(method = "loess", formula = y ~ log(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 12) ) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ))

pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.nopoints.v2.pdf"), height=5, width=10)
print(plot_CRonR_global_allFU_nopoints)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.nopoints.v2.png"), height=5, width=10, units = "in", res = 300)
print(plot_CRonR_global_allFU_nopoints)
dev.off()

plot_CRonR_zscaled_global_allFU_nopoints <- ggplot(allstudies_profile_skewing_fulldata_notNA_stats, aes(x = FollowUp, y = CRonR_zscaled, fill = Class, color = Class), na.rm = T, se = TRUE) +
  # geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  # geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  geom_smooth(method = "loess", formula = y ~ log(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 12) ) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ),
       x = "Months after gene therapy", y = "% of IS (Z-score)", 
       colour = "Class")

pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.zscaled.nopoints.v2.pdf"), height=5, width=10)
print(plot_CRonR_zscaled_global_allFU_nopoints)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.zscaled.nopoints.v2.png"), height=5, width=10, units = "in", res = 300)
print(plot_CRonR_zscaled_global_allFU_nopoints)
dev.off()

plot_CRonR_global_allFU_byage <- ggplot(allstudies_profile_skewing_fulldata_notNA_stats, aes(x = FollowUp, y = CRonR, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 12) ) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(AgeGroup_Sorted ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ))

pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.byAge.pdf"), height=10, width=12)
print(plot_CRonR_global_allFU_byage)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.byAge.png"), height=10, width=12, units = "in", res = 300)
print(plot_CRonR_global_allFU_byage)
dev.off()


plot_CSonCU_global <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CSonCU, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 12) ) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Singletons IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CS on CU (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance) ))

plot_CRonCU_global <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CRonCU, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 12) ) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CR on CU (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance) ))

plot_CSonCR_global <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CSonCR, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.2) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  # scale_y_log10(
  #   breaks = scales::trans_breaks("log10", function(x) 10^x),
  #   labels = scales::trans_format("log10", scales::math_format(10^.x))
  # ) +
  # annotation_logticks(sides = "lr") +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 12) ) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Fold change Singletons vs Recaptured"), x = "Months after gene therapy", y = "FC on N.IS", colour = "Class", subtitle = paste0("CS on CR (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance) ))

plot_CSonS_global <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CSonS, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 12) ) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Singletons IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CS on S (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance) ))

plot_CRonR_global <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CRonR, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 12) ) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance) ))

plot_CSonU_global <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CSonU, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 12) ) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Singletons IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CS on U (C:Class, S:Singletons, R:Recaptured, U:All). Equivalent to rescaled CS/S: a*(CS/S), where a=S/U.", " Abundance threshold = ", (threshold_abundance) ))

plot_CSonU_global_nopoints <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CSonU, fill = Class, color = Class), na.rm = T, se = TRUE) +
  # geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 12) ) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Singletons IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CS on U (C:Class, S:Singletons, R:Recaptured, U:All). Equivalent to rescaled CS/S: a*(CS/S), where a=S/U.", " Abundance threshold = ", (threshold_abundance) ))

plot_CRonU_global <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CRonU, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 12) ) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CR on U (C:Class, S:Singletons, R:Recaptured, U:All). Equivalent to rescaled CR/R: a*(CR/R), where a=R/U.", " Abundance threshold = ", (threshold_abundance) ))

blank <- grid.rect(gp=gpar(col="white"))

pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".flag_matrix_relabeled_correctedDC_freqs_melt_merge.Classes.SingletonsVSMulti.facetStudies", suffix_plotfile, ".pdf"), height=7, width=16)
print(plot_CSonU_global) 
print(plot_CRonU_global)
print(plot_CSonS_global) 
print(plot_CRonR_global) 
dev.off()

plot_footnote_details <- "Corrected DC"
plot_title <- paste0("Short vs Long -term IS, by classes using Flags, relabeled and corrected DC")
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".flag_matrix_relabeled_correctedDC_freqs_melt_merge.Classes.SingletonsVSMulti.facetStudies", suffix_plotfile, ".png"), height=20, width=12, units = "in", res = 300)
grid.arrange(plot_CSonU_global,
             plot_CRonU_global,
             plot_CSonS_global,
             plot_CRonR_global,
             ncol = 1, 
             top = textGrob(plot_title, gp = gpar(fontsize = 22)),
             bottom = textGrob(plot_footnote_details, gp = gpar(fontface = 3, fontsize = 18), hjust=1, x = 1)
)
dev.off()


######

allstudies_profile_skewing_fulldata_notNA_avgont <- allstudies_profile_skewing_fulldata_notNA_stats %>% group_by(ProjectID, Class, FollowUp) %>% summarise(Avg_CRonR = mean(CRonR, na.rm = TRUE))

allstudies_profile_skewing_fulldata_notNA_avgont_cast <- dcast(data = allstudies_profile_skewing_fulldata_notNA_avgont, fun.aggregate = mean, formula = FollowUp ~ Class + ProjectID, value.var = "Avg_CRonR")
allstudies_profile_skewing_fulldata_notNA_avgont_cast[is.na(allstudies_profile_skewing_fulldata_notNA_avgont_cast)] <- NA
write.xlsx(x = allstudies_profile_skewing_fulldata_notNA_avgont_cast, file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, "allstudies_profile_skewing_fulldata_notNA_avgont_cast.xlsx"))

scale_color_manual_colors_sortedbyname <- as.character(color_schema_df_sortedbyname[which(color_schema_df_sortedbyname$ClassName %in% allstudies_profile_skewing_fulldata_notNA$Class), "colorcode"])

my_comparisons <- list(
  # c("Multilineage", "UniErythroid"),
  # c("Multilineage", "UniLyB"),
  # c("Multilineage", "UniLyT"),
  # c("Multilineage", "UniMyeloid"),
  c("UniMyeloid", "UniErythroid"),
  c("UniMyeloid", "UniLyB"),
  c("UniMyeloid", "UniLyT"),
  c("UniLyT", "UniErythroid"),
  c("UniLyB", "UniErythroid")
  #UniLyB UniLymphoid UniLyT UniMyeloid
)
# p <- ggplot(data = stats, aes(x = factor(TestGroup), y = Perc_onSharedCD34BM), na.rm = T, se = TRUE)
data_slice <- allstudies_profile_skewing_fulldata_notNA_avgont[which(!(allstudies_profile_skewing_fulldata_notNA_avgont$Class %in% c("UniLymphoid", "Multilineage") )),]
plot_stats_multiuni_anova1 <- ggplot(data = data_slice, aes(x = Class, y = Avg_CRonR), color = ProjectID, na.rm = T, se = TRUE) +
  # geom_boxplot(alpha = .7, outlier.size = 0) +
  # geom_violin() +
  geom_jitter(aes(colour = Class), height = 0, width = 0.1, alpha = 0.7, size = 3) +
  facet_grid(. ~ ProjectID) +
  # geom_point(size=4, alpha = .7, position = "jitter") +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(1, 1, 2, 2))+
  stat_compare_means(comparisons = my_comparisons, label.y = seq(from = 60, to = 95, by = 5), label = "p.signif") +
  theme_bw() +
  stat_compare_means(label.y = 100) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16, angle = 30), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Multilineage potential over time"), x = "Lineages", y = "Sharing (%)", colour = "Lineage", fill = "Lineage"
  )

pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, "allstudies_profile_skewing_fulldata_notNA_avgont_cast.testv0.pdf"), height=6, width=10)
print(plot_stats_multiuni_anova1)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, "allstudies_profile_skewing_fulldata_notNA_avgont_cast.testv0.png"), height=6, width=10, units = "in", res = 300)
print(plot_stats_multiuni_anova1)
dev.off()

my_comparisons <- list(
  c("MLD", "WAS"),
  c("MLD", "BTHAL"),
  c("BTHAL", "WAS")
)
# p <- ggplot(data = stats, aes(x = factor(TestGroup), y = Perc_onSharedCD34BM), na.rm = T, se = TRUE)
# data_slice <- allstudies_profile_skewing_fulldata_notNA_avgont[which(!(allstudies_profile_skewing_fulldata_notNA_avgont$Class %in% c("UniLymphoid", "Multilineage") )),]
data_slice <- allstudies_profile_skewing_fulldata_notNA_avgont[which(!(allstudies_profile_skewing_fulldata_notNA_avgont$Class %in% c("UniLymphoid") )),]
data_slice$Study <- factor(data_slice$ProjectID, levels = study_list)

n_fun <- function(x){
  return(data.frame(y = c(100),
                    label = paste0("N=", length(x)) ))
}

plot_stats_multiuni_anova2 <- ggplot(data = data_slice, aes(x = Study, y = Avg_CRonR), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  geom_boxplot(alpha = .5, fill = "white", outlier.size = 0) +
  # geom_violin() +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size = 5) +
  geom_jitter(aes(colour = Class), height = 0, width = 0.1, alpha = 0.7, size = 3) +
  facet_grid(. ~ Class) +
  # geom_point(size=4, alpha = .7, position = "jitter") +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(1, 1, 2, 2))+
  # stat_compare_means(comparisons = my_comparisons, label.y = seq(from = 60, to = 95, by = 5), label = "p.signif") +
  stat_compare_means(comparisons = my_comparisons, label.y = seq(from = 35, to = 95, by = 5)) +
  theme_bw() +
  # stat_compare_means(label.y = 100) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16, angle = 45), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("BM - Lineage commitment over time"), x = "Disease", y = "Sharing (%)", colour = "Lineage", fill = "Lineage"
  )

# ggplot(data = data_slice, aes(x = FollowUp, y = Avg_CRonR, group = Class, color = Class, fill = Class), na.rm = T, se = TRUE) +
#   # geom_boxplot(alpha = .7, outlier.size = 0) +
#   # geom_violin() +
#   # geom_jitter(aes(colour = Class), height = 0, width = 0.1, alpha = 0.7, size = 3) +
#   facet_grid(. ~ ProjectID) +
#   geom_point(size=4, alpha = .7, position = "jitter") +
#   geom_line() +
#   # stat_compare_means(comparisons = my_comparisons, label.y = c(1, 1, 2, 2))+
#   # stat_compare_means(comparisons = my_comparisons, label.y = seq(from = 60, to = 95, by = 5), label = "p.signif") +
#   theme_bw() +
#   # stat_compare_means(label.y = 100) +
#   theme(strip.text = element_text(face="bold", size=16)) +
#   theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
#   theme(axis.text.x = element_text(size=16, angle = 30), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
#   labs(title = paste0("Multilineage potential over time"), x = "Lineages", y = "Sharing (%)", colour = "Lineage", fill = "Lineage"
#   )

pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, "allstudies_profile_skewing_fulldata_notNA_avgont_cast.testv2.pdf"), height=5.5, width=11)
print(plot_stats_multiuni_anova2)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, "allstudies_profile_skewing_fulldata_notNA_avgont_cast.testv2.png"), height=5.5, width=11, units = "in", res = 300)
print(plot_stats_multiuni_anova2)
dev.off()



# compare zscaled values
my_comparisons <- list(
  c("BTHAL", "WAS"),
  c("MLD", "BTHAL"),
  c("MLD", "WAS")
)

# n_fun <- function(x){
#   return(data.frame(y = c(max(patient_skewing_summary_filtered$mean)*1.1),
#                     label = paste0("N=", length(x)) ))
# }
n_fun <- function(x){
  return(data.frame(y = c(min(patient_skewing_summary_filtered$mean)*1.1),
                    label = paste0("N=", length(x)) ))
}
perc_starting_max_stats <- 20
perc_spacer_stats <- 2
plot_stats_multiuni_zscaled_anova2 <- ggplot(data = global_patient_skewing_summary_filtered, aes(x = Study, y = mean), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  geom_boxplot(alpha = .5, fill = "white", outlier.size = 0) +
  # geom_violin() +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size = 5) +
  geom_jitter(aes(colour = Class), height = 0, width = 0.1, alpha = 0.7, size = 3) +
  facet_grid(Tissue ~ Class) +
  # geom_point(size=4, alpha = .7, position = "jitter") +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(1, 1, 2, 2))+
  stat_compare_means(comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
  # stat_compare_means(comparisons = my_comparisons, label.y = seq(from = 2.5, to = 1.5, by = -0.5), label = "p.signif", p.adjust.method = "fdr") +
  theme_bw() +
  # stat_compare_means(label.y = 100) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16, angle = 0), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Lineage commitment over time"), 
       x = "Disease", y = "Percentage scaled Z-score", colour = "Lineage", fill = "Lineage"
  )

pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, "patient_skewing_summary_filtered.CRonR.zscore.pdf"), height=5, width=10)
print(plot_stats_multiuni_zscaled_anova2)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, "patient_skewing_summary_filtered.CRonR.zscore.png"), height=5, width=10, units = "in", res = 300)
print(plot_stats_multiuni_zscaled_anova2)
dev.off()




#####------ PB ---------------------
faceting_order <- c("MLD", "WAS", "BTHAL")
lineage_faceting_order <- c("Multilineage", "Myeloid", "B", "T")
trials_colors <- c("darkblue", "forestgreen", "firebrick")
prefix_trials_outname <- paste(faceting_order, collapse = "_")
results_folder_name <- "multi_uni_sharing/"
this_run_foldername <- "MyBTEry34_allMarkers_noLy_PB_SC3/" # "MyBTEry34_strictMarkers_fromwhole_noLy_AllMarkers/" # "MyBTEry34/" "MyBTEry34_ISgt6/" "MyBTEry34_BMonly/" "MyBTEry34_PBonly/"
# this_run_foldername_in_was <- "MyBTEry34_strictMarkers/" # "MyBTEry34/" "MyBTEry34_ISgt6/" "MyBTEry34_BMonly/" "MyBTEry34_PBonly/"
threshold_abundance <- 1 # < threshold_abundance will be removed
out_xlsxfile_infix <- "gdf_hlfu_filtered"
out_file_infix <- ".hlfu_34MyBTEry"
prefix_dest_folder <- "AllPts_"
suffix_plotfile <- "" # ".noTransfDep"
followup_limit_crosstrial <- 36 # plot up to this month after GT
excel_bagcases_foldername <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial MLD/analyses/13.hsc_profiles/"
excel_bagcases_filename <- "MLD.hlfu_34PBMyBT_noLy_no34exclusive.flag_matrix.bagcases.xlsx" # "MLD.hlfu_34MyBTEry.flag_matrix.bagcases.xlsx"  # "MLD.hlfu_34PBMyBT.flag_matrix.bagcases.xlsx" # MLD.hlfu_34MyBTEry.flag_matrix.bagcases.xlsx MLD.hlfu_34MyBTEry_noLy.flag_matrix.bagcases.xlsx MLD.hlfu_34PBMyBT_noLy.flag_matrix.bagcases.xlsx
# patients_to_exclude <- c("BTHAL001", "BTHAL003", "BTHAL004", "BTHAL006", "BTHAL010")
patients_to_exclude <- c("MLD13")
name_list_trials <- prefix_trials_outname
# patients_to_exclude <- c()

color_schema_df <- read.xlsx(xlsxFile = paste(excel_bagcases_foldername, excel_bagcases_filename, sep = ""), sheet = "color_schema")
color_schema_df_sortedbyname <- color_schema_df[order(color_schema_df$ClassName),]
color_schema_df_sortedbyplotlabel <- color_schema_df[order(color_schema_df$Label),]

# add all trials data
mld_allpatients_flag_matrix_relabeled_correctedDC_freqs_melt_merge <- read.csv(file = paste0("/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial MLD/analyses/13.hsc_profiles/202112/", this_run_foldername, "MLD.AllPatients.hlfu_excl34MyBTEry.flag_matrix_relabeled_correctedDC_freqs_melt_merge.tsv"), header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", ""))
was_allpatients_flag_matrix_relabeled_correctedDC_freqs_melt_merge <- read.csv(file = paste0("/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial WAS/analyses/14.hsc_profiles/202112/", this_run_foldername, "WAS.AllPatients.hlfu_34MyBTEry.flag_matrix_relabeled_correctedDC_freqs_melt_merge.tsv"), header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", ""))
# bthal_allpatients_flag_matrix_relabeled_correctedDC_freqs_melt_merge <- read.csv(file = paste0("/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial bThalassemia/analyses/14.hsc_profiles/202103/", this_run_foldername, "/BTHAL.AllPatients.hlfu_34MyBTEry.flag_matrix_relabeled_correctedDC_freqs_melt_merge.tsv"), header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", ""))
bthal_allpatients_flag_matrix_relabeled_correctedDC_freqs_melt_merge <- read.csv(file = paste0("/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial bThalassemia/analyses/14.hsc_profiles/202112/", this_run_foldername, "BTHAL.AllPatients.hlfu_34MyBTEry.flag_matrix_relabeled_correctedDC_freqs_melt_merge.tsv"), header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", ""))


allstudies_profile_skewing_fulldata <- Reduce(rbind, list(mld_allpatients_flag_matrix_relabeled_correctedDC_freqs_melt_merge, was_allpatients_flag_matrix_relabeled_correctedDC_freqs_melt_merge, bthal_allpatients_flag_matrix_relabeled_correctedDC_freqs_melt_merge))
allstudies_profile_skewing_fulldata$ClassName <- allstudies_profile_skewing_fulldata$Class
allstudies_profile_skewing_fulldata <- merge(x = allstudies_profile_skewing_fulldata, y = color_schema_df[c("ClassName", "colorcode", "ColorSchema", "PlotLabel")], by = c("ClassName"))
allstudies_profile_skewing_fulldata$ClinicalStudy <- factor(allstudies_profile_skewing_fulldata$ProjectID, 
                                                            levels = faceting_order)
allstudies_profile_skewing_fulldata$PlotLabel <- factor(allstudies_profile_skewing_fulldata$PlotLabel, 
                                                            levels = lineage_faceting_order)

allstudies_profile_skewing_fulldata_notNA <- allstudies_profile_skewing_fulldata[which(!is.na(allstudies_profile_skewing_fulldata$SonR)),]


# adjustments of values
allstudies_profile_skewing_fulldata_notNA$CSonCR <- ifelse((allstudies_profile_skewing_fulldata_notNA$CR > 0 & allstudies_profile_skewing_fulldata_notNA$CU > 0), allstudies_profile_skewing_fulldata_notNA$CS / allstudies_profile_skewing_fulldata_notNA$CR, NA)
allstudies_profile_skewing_fulldata_notNA$CSonCU <- ifelse((allstudies_profile_skewing_fulldata_notNA$CU > 0), allstudies_profile_skewing_fulldata_notNA$CS *100 / allstudies_profile_skewing_fulldata_notNA$CU, NA)
allstudies_profile_skewing_fulldata_notNA$CRonCU <- ifelse((allstudies_profile_skewing_fulldata_notNA$CU > 0), allstudies_profile_skewing_fulldata_notNA$CR *100 / allstudies_profile_skewing_fulldata_notNA$CU, NA)
allstudies_profile_skewing_fulldata_notNA$CSonU <- ifelse((allstudies_profile_skewing_fulldata_notNA$U > 0), allstudies_profile_skewing_fulldata_notNA$CS *100 / allstudies_profile_skewing_fulldata_notNA$U, NA)
allstudies_profile_skewing_fulldata_notNA$CRonU <- ifelse((allstudies_profile_skewing_fulldata_notNA$U > 0), allstudies_profile_skewing_fulldata_notNA$CR *100 / allstudies_profile_skewing_fulldata_notNA$U, NA)
allstudies_profile_skewing_fulldata_notNA$CSonS <- ifelse((allstudies_profile_skewing_fulldata_notNA$S > 0), allstudies_profile_skewing_fulldata_notNA$CS *100 / allstudies_profile_skewing_fulldata_notNA$S, NA)
allstudies_profile_skewing_fulldata_notNA$CRonR <- ifelse((allstudies_profile_skewing_fulldata_notNA$R > 0), allstudies_profile_skewing_fulldata_notNA$CR *100 / allstudies_profile_skewing_fulldata_notNA$R, NA)

# patients_to_exclude <- c()
allstudies_profile_skewing_fulldata_notNA <- allstudies_profile_skewing_fulldata_notNA[which(!(allstudies_profile_skewing_fulldata_notNA$SubjectID %in% patients_to_exclude) ),]

allstudies_profile_skewing_fulldata_notNA_stats <- merge(x = allstudies_profile_skewing_fulldata_notNA, y = patients_summary, by = c("SubjectID"), all.x = T)
# allpatients_shared34_stats$AgeGroup <- 
allstudies_profile_skewing_fulldata_notNA_stats <- allstudies_profile_skewing_fulldata_notNA_stats %>%
  mutate(AgeGroup = case_when(
    TreatmentAge.years <= 2 ~"0-2",
    TreatmentAge.years <= 15 ~"2-15",
    TreatmentAge.years > 15 ~"30+"
  ))
age_list <- c("0-2", "2-15", "30+")
allstudies_profile_skewing_fulldata_notNA_stats$AgeGroup_Sorted <- factor(allstudies_profile_skewing_fulldata_notNA_stats$AgeGroup, levels = age_list)


# add the Tissue only here
allstudies_profile_skewing_fulldata_notNA_stats$Tissue <- "PB"
allstudies_profile_skewing_fulldata_notNA_stats <- allstudies_profile_skewing_fulldata_notNA_stats %>% group_by(ProjectID, SubjectID, Tissue, FollowUp) %>% mutate(CRonR_zscaled = scale(CRonR), CSonU_zscaled = scale(CSonU), CRonU_zscaled = scale(CRonU), CSonS_zscaled = scale(CSonS))
allstudies_profile_skewing_fulldata_notNA_stats[is.na(allstudies_profile_skewing_fulldata_notNA_stats)] <- NA
# copy to the global df
if (length(global_allstudies_profile_skewing_fulldata_notNA_stats) > 0) {
  global_allstudies_profile_skewing_fulldata_notNA_stats <- rbind(global_allstudies_profile_skewing_fulldata_notNA_stats, allstudies_profile_skewing_fulldata_notNA_stats)
} else {
  global_allstudies_profile_skewing_fulldata_notNA_stats <- allstudies_profile_skewing_fulldata_notNA_stats
}
# slice data at stability
allstudies_profile_skewing_fulldata_notNA_stats_stable <- allstudies_profile_skewing_fulldata_notNA_stats[which(allstudies_profile_skewing_fulldata_notNA_stats$FollowUp >= 24 & allstudies_profile_skewing_fulldata_notNA_stats$FollowUp <= 60),]
allstudies_profile_skewing_fulldata_notNA_stats_stable <- allstudies_profile_skewing_fulldata_notNA_stats %>% group_by(SubjectID, Tissue, Class) %>% mutate(n_obs_stt_overtime = n())
allstudies_profile_skewing_fulldata_notNA_stats_stable_min3 <- allstudies_profile_skewing_fulldata_notNA_stats %>% group_by(SubjectID, Tissue, Class) %>% mutate(n_obs_stt_overtime = n()) %>% filter(n_obs_stt_overtime >= 3)


# summary stats
patient_skewing_summary <- allstudies_profile_skewing_fulldata_notNA_stats_stable_min3 %>%
  filter(n_obs_stt_overtime >= 3) %>%
  group_by(ClinicalStudy, SubjectID, Tissue, Class, TreatmentAge.years, AgeGroup, AgeGroup_Sorted, colorcode, ColorSchema, PlotLabel) %>%
  get_summary_stats(CRonR_zscaled, type = "full") 

patient_skewing_summary_filtered <- patient_skewing_summary %>% filter(n >= 2)
# patient_skewing_summary_filtered$CellLineage <- factor(patient_skewing_summary_filtered$CellType, levels = celltype_list)
patient_skewing_summary_filtered$Study <- factor(patient_skewing_summary_filtered$ClinicalStudy, levels = study_list)
patient_skewing_summary_filtered$AgeGroup_Sorted <- factor(patient_skewing_summary_filtered$AgeGroup, levels = age_list)
write.xlsx(x = patient_skewing_summary_filtered, file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.StableFU.CRonR.xlsx"))

global_patient_skewing_summary_filtered <- rbind(global_patient_skewing_summary_filtered, patient_skewing_summary_filtered)


allstudies_profile_skewing_fulldata_byt <- sqldf("select ProjectID, ClinicalStudy, SubjectID, Timepoint, FollowUp, SumOverT_OverallNIS, SumOverT_SingletoneNIS, SumOverT_SharingNIS, U, S, R, SumOverT_Singletone_OnTotal, SumOverT_Sharing_OnTotal, SumOverT_SingletonOnSharing, SonU, RonU, SonR from allstudies_profile_skewing_fulldata_notNA where R > 0 group by ProjectID, SubjectID, Timepoint, FollowUp")

write.xlsx(x = allstudies_profile_skewing_fulldata_notNA, file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, ".allstudies_profile_skewing_fulldata_notNA.xlsx"))
write.xlsx(x = allstudies_profile_skewing_fulldata_byt, file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, ".allstudies_profile_skewing_fulldata_byt.xlsx"))

scale_color_manual_colors_sortedbyname <- as.character(color_schema_df_sortedbyname[which(color_schema_df_sortedbyname$ClassName %in% levels(factor(allstudies_profile_skewing_fulldata_notNA$Class))), "colorcode"])
scale_color_manual_colors_sortedbyplotlabel <- as.character(color_schema_df_sortedbyplotlabel[which(color_schema_df_sortedbyplotlabel$PlotLabel %in% levels(factor(allstudies_profile_skewing_fulldata_notNA$PlotLabel))), "colorcode"])
# plot_CRonR_global_allFU <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CRonR, fill = Class, color = Class), na.rm = T, se = TRUE) +
plot_CRonR_global_allFU <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CRonR, fill = PlotLabel, color = PlotLabel), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyplotlabel) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyplotlabel) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 6) ) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), x = "Months after gene therapy", y = "% of IS", fill = "Lineage direction", colour = "Lineage direction", 
       subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ))

pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.pdf"), height=7, width=12)
print(plot_CRonR_global_allFU)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.png"), height=7, width=12, units = "in", res = 300)
print(plot_CRonR_global_allFU)
dev.off()


plot_CRonR_global_allFU_nopoints <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CRonR, fill = Class, color = Class), na.rm = T, se = TRUE) +
  # geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 6) ) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ))

pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.nopoints.pdf"), height=5, width=10)
print(plot_CRonR_global_allFU_nopoints)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.nopoints.png"), height=5, width=10, units = "in", res = 300)
print(plot_CRonR_global_allFU_nopoints)
dev.off()

plot_CRonR_zscaled_global_allFU_nopoints <- ggplot(allstudies_profile_skewing_fulldata_notNA_stats, aes(x = FollowUp, y = CRonR_zscaled, fill = Class, color = Class), na.rm = T, se = TRUE) +
  # geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 6) ) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ),
       x = "Months after gene therapy", y = "% of IS (Z-score)", 
       colour = "Class")
pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.zscaled.nopoints.pdf"), height=5, width=10)
print(plot_CRonR_zscaled_global_allFU_nopoints)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.zscaled.nopoints.png"), height=5, width=10, units = "in", res = 300)
print(plot_CRonR_zscaled_global_allFU_nopoints)
dev.off()

plot_CRonR_global_allFU_byage <- ggplot(allstudies_profile_skewing_fulldata_notNA_stats, aes(x = FollowUp, y = CRonR, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 6) ) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(AgeGroup_Sorted ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ))

pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.byAge.pdf"), height=10, width=12)
print(plot_CRonR_global_allFU_byage)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.byAge.png"), height=10, width=12, units = "in", res = 300)
print(plot_CRonR_global_allFU_byage)
dev.off()


plot_CRonR_global_allFU_byage <- ggplot(allstudies_profile_skewing_fulldata_notNA_stats, aes(x = FollowUp, y = CRonR, fill = PlotLabel, color = PlotLabel), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyplotlabel) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyplotlabel) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 6) ) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(AgeGroup_Sorted ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), x = "Months after gene therapy", y = "% of IS", fill = "Lineage direction", colour = "Lineage direction", 
       subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ))

pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.byAge.pdf"), height=10, width=12)
print(plot_CRonR_global_allFU_byage)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.byAge.png"), height=10, width=12, units = "in", res = 300)
print(plot_CRonR_global_allFU_byage)
dev.off()


plot_CSonCU_global <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CSonCU, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 6) ) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Singletons IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CS on CU (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance) ))

plot_CRonCU_global <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CRonCU, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 6) ) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CR on CU (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance) ))

plot_CSonCR_global <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CSonCR, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.2) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  # scale_y_log10(
  #   breaks = scales::trans_breaks("log10", function(x) 10^x),
  #   labels = scales::trans_format("log10", scales::math_format(10^.x))
  # ) +
  # annotation_logticks(sides = "lr") +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 6) ) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Fold change Singletons vs Recaptured"), x = "Months after gene therapy", y = "FC on N.IS", colour = "Class", subtitle = paste0("CS on CR (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance) ))

plot_CSonS_global <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CSonS, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 6) ) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Singletons IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CS on S (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance) ))

plot_CRonR_global <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CRonR, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 6) ) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance) ))

plot_CSonU_global <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CSonU, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 6) ) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Singletons IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CS on U (C:Class, S:Singletons, R:Recaptured, U:All). Equivalent to rescaled CS/S: a*(CS/S), where a=S/U.", " Abundance threshold = ", (threshold_abundance) ))

plot_CRonU_global <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CRonU, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_point(size = 2, alpha = 0.7) +
  # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
  geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 6) ) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CR on U (C:Class, S:Singletons, R:Recaptured, U:All). Equivalent to rescaled CR/R: a*(CR/R), where a=R/U.", " Abundance threshold = ", (threshold_abundance) ))

blank <- grid.rect(gp=gpar(col="white"))

pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".flag_matrix_relabeled_correctedDC_freqs_melt_merge.Classes.SingletonsVSMulti.facetStudies", suffix_plotfile, ".pdf"), height=7, width=16)
print(plot_CSonU_global) 
print(plot_CRonU_global)
print(plot_CSonS_global) 
print(plot_CRonR_global) 
dev.off()

plot_footnote_details <- "Corrected DC"
plot_title <- paste0("Short vs Long -term IS, by classes using Flags, relabeled and corrected DC")
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, paste(study_list, collapse = "_"), ".flag_matrix_relabeled_correctedDC_freqs_melt_merge.Classes.SingletonsVSMulti.facetStudies", suffix_plotfile, ".png"), height=20, width=12, units = "in", res = 300)
grid.arrange(plot_CSonU_global,
             plot_CRonU_global,
             plot_CSonS_global,
             plot_CRonR_global,
             ncol = 1, 
             top = textGrob(plot_title, gp = gpar(fontsize = 22)),
             bottom = textGrob(plot_footnote_details, gp = gpar(fontface = 3, fontsize = 18), hjust=1, x = 1)
)
dev.off()



# scale_color_manual_colors_sortedbyname <- as.character(color_schema_df_sortedbyname[which(color_schema_df_sortedbyname$ClassName %in% allstudies_profile_skewing_fulldata_notNA$Class), "colorcode"])
# plot_CRonR_global <- ggplot(allstudies_profile_skewing_fulldata_notNA, aes(x = FollowUp, y = CRonR, fill = Class, color = Class), na.rm = T, se = TRUE) +
#   geom_point(size = 2, alpha = 0.7) +
#   # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
#   geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
#   # scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
#   # scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
#   # scale_y_continuous(limits = c(0, 100)) +
#   scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 6) ) +
#   facet_wrap( ~ ClinicalStudy, ncol = 3) +
#   theme(strip.text = element_text(face="bold", size=16)) +
#   theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
#   theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
#   labs(title = paste0("Recaptured IS in percentage"), x = "Months after gene therapy", y = "% of IS", colour = "Class", subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance) ))


#### stats! ---- mean among all patients over time -> track the profile
# allstudies_profile_skewing_fulldata_notNA_avgont <- allstudies_profile_skewing_fulldata_notNA %>% group_by(ProjectID, Class, FollowUp) %>% mutate(Avg_CRonR = mean(CRonR, na.rm = TRUE))
allstudies_profile_skewing_fulldata_notNA_avgont <- allstudies_profile_skewing_fulldata_notNA %>% group_by(ProjectID, Class, PlotLabel, FollowUp) %>% summarise(Avg_CRonR = mean(CRonR, na.rm = TRUE), n_patients = n())

allstudies_profile_skewing_fulldata_notNA_avgont_cast <- dcast(data = allstudies_profile_skewing_fulldata_notNA_avgont, fun.aggregate = mean, formula = FollowUp ~ Class + PlotLabel + ProjectID, value.var = "Avg_CRonR")
allstudies_profile_skewing_fulldata_notNA_avgont_cast[is.na(allstudies_profile_skewing_fulldata_notNA_avgont_cast)] <- NA
write.xlsx(x = allstudies_profile_skewing_fulldata_notNA_avgont_cast, file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, ".allstudies_profile_skewing_fulldata_notNA_avgont_cast.xlsx"))

scale_color_manual_colors_sortedbyname <- as.character(color_schema_df_sortedbyname[which(color_schema_df_sortedbyname$ClassName %in% levels(factor(allstudies_profile_skewing_fulldata_notNA$Class))), "colorcode"])
scale_color_manual_colors_sortedbyplotlabel <- as.character(color_schema_df_sortedbyplotlabel[which(color_schema_df_sortedbyplotlabel$PlotLabel %in% levels(factor(allstudies_profile_skewing_fulldata_notNA$PlotLabel))), "colorcode"])

my_comparisons <- list(
  # c("Multilineage", "UniErythroid"),
  # c("Multilineage", "UniLyB"),
  # c("Multilineage", "UniLyT"),
  # c("Multilineage", "UniMyeloid"),
  c("UniMyeloid", "UniLyB"),
  c("UniMyeloid", "UniLyT"),
  c("UniLyB", "UniLyT")
  #UniLyB UniLymphoid UniLyT UniMyeloid
)
# p <- ggplot(data = stats, aes(x = factor(TestGroup), y = Perc_onSharedCD34BM), na.rm = T, se = TRUE)
data_slice <- allstudies_profile_skewing_fulldata_notNA_avgont[which(!(allstudies_profile_skewing_fulldata_notNA_avgont$Class %in% c("UniLymphoid") ) 
                                                                     & allstudies_profile_skewing_fulldata_notNA_avgont$n_patients > 2),]
plot_stats_multiuni_anova1 <- ggplot(data = data_slice, aes(x = Class, y = Avg_CRonR), color = ProjectID, na.rm = T, se = TRUE) +
  # geom_boxplot(alpha = .7, outlier.size = 0) +
  # geom_violin() +
  geom_jitter(aes(colour = Class), height = 0, width = 0.1, alpha = 0.7, size = 3) +
  facet_grid(. ~ ProjectID) +
  # geom_point(size=4, alpha = .7, position = "jitter") +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(1, 1, 2, 2))+
  stat_compare_means(comparisons = my_comparisons, label.y = seq(from = 60, to = 95, by = 5), label = "p.signif") +
  theme_bw() +
  stat_compare_means(label.y = 100) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16, angle = 30), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Multilineage potential over time"), x = "Lineages", y = "Sharing (%)", colour = "Lineage", fill = "Lineage"
  )
  
  pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, "allstudies_profile_skewing_fulldata_notNA_avgont_cast.testv0.pdf"), height=6, width=8)
  print(plot_stats_multiuni_anova1)
  dev.off()
  png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, "allstudies_profile_skewing_fulldata_notNA_avgont_cast.testv0.png"), height=6, width=8, units = "in", res = 300)
  print(plot_stats_multiuni_anova1)
  dev.off()

  
my_comparisons <- list(
  c("MLD", "BTHAL"),
  c("MLD", "WAS"),
  c("BTHAL", "WAS")
)
# p <- ggplot(data = stats, aes(x = factor(TestGroup), y = Perc_onSharedCD34BM), na.rm = T, se = TRUE)
# data_slice <- allstudies_profile_skewing_fulldata_notNA_avgont[which(!(allstudies_profile_skewing_fulldata_notNA_avgont$Class %in% c("UniLymphoid", "Multilineage") )),]
data_slice <- allstudies_profile_skewing_fulldata_notNA_avgont[which(!(allstudies_profile_skewing_fulldata_notNA_avgont$Class %in% c("UniLymphoid") ) 
                                                                     & allstudies_profile_skewing_fulldata_notNA_avgont$n_patients > 2),]
data_slice$Study <- factor(data_slice$ProjectID, levels = study_list)

n_fun <- function(x){
  return(data.frame(y = c(100),
                    label = paste0("N=", length(x)) ))
}

plot_stats_multiuni_anova2 <- ggplot(data = data_slice, aes(x = Study, y = Avg_CRonR), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyplotlabel) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyplotlabel) +
  geom_boxplot(alpha = .5, fill = "white", outlier.size = 0) +
  # geom_violin() +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size = 5) +
  geom_jitter(aes(colour = PlotLabel), height = 0, width = 0.1, alpha = 0.7, size = 3) +
  facet_grid(. ~ PlotLabel) +
  # geom_point(size=4, alpha = .7, position = "jitter") +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(1, 1, 2, 2))+
  # stat_compare_means(comparisons = my_comparisons, label.y = seq(from = 60, to = 95, by = 5), label = "p.signif") +
  stat_compare_means(comparisons = my_comparisons, p.adjust.method = "fdr", label.y = seq(from = 35, to = 95, by = 5)) +
  theme_bw() +
  # stat_compare_means(label.y = 100) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16, angle = 0), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("PB - Lineage commitment over time"), x = "Disease", y = "Sharing (%)", colour = "Lineage", fill = "Lineage",
       subtitle = paste0("Test on means, p-value corrected (FDR). Data filtering: points with N patients < 2.")
  )


pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, "allstudies_profile_skewing_fulldata_notNA_avgont_cast.testv2.pdf"), height=5.5, width=9.5)
print(plot_stats_multiuni_anova2)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, "allstudies_profile_skewing_fulldata_notNA_avgont_cast.testv2.png"), height=5.5, width=9.5, units = "in", res = 300)
print(plot_stats_multiuni_anova2)
dev.off()


# now try with the age stratification
allstudies_profile_skewing_fulldata_notNA_avgont_age <- allstudies_profile_skewing_fulldata_notNA_stats %>% group_by(ProjectID, Class, PlotLabel, FollowUp, AgeGroup) %>% summarise(Avg_CRonR = mean(CRonR, na.rm = TRUE), n_patients = n())

allstudies_profile_skewing_fulldata_notNA_avgont_age_cast <- dcast(data = allstudies_profile_skewing_fulldata_notNA_avgont_age, fun.aggregate = mean, formula = FollowUp ~ Class + PlotLabel + AgeGroup + ProjectID, value.var = "Avg_CRonR")
allstudies_profile_skewing_fulldata_notNA_avgont_age_cast[is.na(allstudies_profile_skewing_fulldata_notNA_avgont_age_cast)] <- NA
write.xlsx(x = allstudies_profile_skewing_fulldata_notNA_avgont_age_cast, file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, ".allstudies_profile_skewing_fulldata_notNA_avgont_age_cast.xlsx"))

scale_color_manual_colors_sortedbyname <- as.character(color_schema_df_sortedbyname[which(color_schema_df_sortedbyname$ClassName %in% levels(factor(allstudies_profile_skewing_fulldata_notNA$Class))), "colorcode"])


my_comparisons <- list(
  # c("MLD", "BTHAL"),
  c("MLD", "WAS")
  # c("BTHAL", "WAS")
)
# p <- ggplot(data = stats, aes(x = factor(TestGroup), y = Perc_onSharedCD34BM), na.rm = T, se = TRUE)
# data_slice <- allstudies_profile_skewing_fulldata_notNA_avgont_age[which(!(allstudies_profile_skewing_fulldata_notNA_avgont_age$Class %in% c("UniLymphoid") )),]
data_slice <- allstudies_profile_skewing_fulldata_notNA_avgont_age[which(!(allstudies_profile_skewing_fulldata_notNA_avgont_age$Class %in% c("UniLymphoid") ) 
                                                                     & allstudies_profile_skewing_fulldata_notNA_avgont_age$n_patients > 2),]
data_slice$Study <- factor(data_slice$ProjectID, levels = study_list)

n_fun <- function(x){
  return(data.frame(y = c(100),
                    label = paste0("N=", length(x)) ))
}

plot_stats_multiuni_anova3 <- ggplot(data = data_slice[which(data_slice$ProjectID %in% c("MLD", "WAS")),], aes(x = Study, y = Avg_CRonR), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyplotlabel) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyplotlabel) +
  geom_boxplot(alpha = .5, fill = "white", outlier.size = 0) +
  # geom_violin() +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size = 5) +
  geom_jitter(aes(colour = PlotLabel), height = 0, width = 0.1, alpha = 0.7, size = 3) +
  facet_grid(AgeGroup ~ PlotLabel) +
  # geom_point(size=4, alpha = .7, position = "jitter") +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(1, 1, 2, 2))+
  # stat_compare_means(comparisons = my_comparisons, label.y = seq(from = 60, to = 95, by = 5), label = "p.signif") +
  stat_compare_means(comparisons = my_comparisons, p.adjust.method = "fdr", label.y = seq(from = 35, to = 95, by = 5)) +
  theme_bw() +
  # stat_compare_means(label.y = 100) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16, angle = 0), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("PB - Lineage commitment over time"), x = "Disease", y = "Sharing (%)", colour = "Lineage", fill = "Lineage",
       subtitle = paste0("Test on means, p-value corrected (FDR). Data filtering: points with N patients < 2.")
  )


pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, "allstudies_profile_skewing_fulldata_notNA_avgont_age.testv2.pdf"), height=8, width=9.5)
print(plot_stats_multiuni_anova3)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, "allstudies_profile_skewing_fulldata_notNA_avgont_age.testv2.png"), height=8, width=9.5, units = "in", res = 300)
print(plot_stats_multiuni_anova3)
dev.off()


my_comparisons <- list(
  c("0-2", "3-21")
  # c("3-21", ">21"),
  # c("0-2", ">21")
)
plot_stats_multiuni_anova4 <- ggplot(data = data_slice[which(data_slice$ProjectID %in% c("MLD", "WAS")),], aes(x = AgeGroup, y = Avg_CRonR), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyplotlabel) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyplotlabel) +
  geom_boxplot(alpha = .5, fill = "white", outlier.size = 0) +
  # geom_violin() +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size = 5) +
  geom_jitter(aes(colour = PlotLabel), height = 0, width = 0.1, alpha = 0.7, size = 3) +
  facet_grid(Study ~ PlotLabel) +
  # geom_point(size=4, alpha = .7, position = "jitter") +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(1, 1, 2, 2))+
  # stat_compare_means(comparisons = my_comparisons, label.y = seq(from = 60, to = 95, by = 5), label = "p.signif") +
  stat_compare_means(comparisons = my_comparisons, p.adjust.method = "fdr", label.y = seq(from = 35, to = 95, by = 5)) +
  theme_bw() +
  # stat_compare_means(label.y = 100) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16, angle = 0), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("PB - Lineage commitment over time"), x = "Disease", y = "Sharing (%)", colour = "Lineage", fill = "Lineage",
       subtitle = paste0("Test on means, p-value corrected (FDR). Data filtering: points with N patients < 2.")
  )


pdf(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, "allstudies_profile_skewing_fulldata_notNA_avgont_age_intratrial.testv2.pdf"), height=8, width=9.5)
print(plot_stats_multiuni_anova4)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", this_run_foldername, prefix_trials_outname, "allstudies_profile_skewing_fulldata_notNA_avgont_age_intratrial.testv2.png"), height=8, width=9.5, units = "in", res = 300)
print(plot_stats_multiuni_anova4)
dev.off()






## ========= using global analysis and data --- re-apply stats just to use plot labels =========
# summary stats
excel_bagcases_foldername <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial MLD/analyses/13.hsc_profiles/"
excel_bagcases_filename <- "MLD.hlfu_34MyBTEry_noLy_no34exclusive.flag_matrix.bagcases.xlsx" # "MLD.hlfu_34MyBTEry.flag_matrix.bagcases.xlsx"  # "MLD.hlfu_34PBMyBT.flag_matrix.bagcases.xlsx" # MLD
color_schema_df <- read.xlsx(xlsxFile = paste(excel_bagcases_foldername, excel_bagcases_filename, sep = ""), sheet = "color_schema")
color_schema_df_sortedbyname <- color_schema_df[order(color_schema_df$ClassName),]
color_schema_df_sortedbyplotlabel <- color_schema_df[order(color_schema_df$PaperOrder),]
scale_color_manual_colors_sortedbyname <- as.character(color_schema_df_sortedbyplotlabel[which(color_schema_df_sortedbyplotlabel$ClassName %in% levels(factor(global_patient_skewing_summary_filtered$Class))), "colorcode"])

lineage_faceting_order <- c("Multilineage", "Erythroid", "Myeloid", "B", "T")
global_patient_skewing_summary_filtered$PlotLabel <- factor(global_patient_skewing_summary_filtered$PlotLabel, levels = lineage_faceting_order)
write.xlsx(x = global_patient_skewing_summary_filtered, file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.StableFU.CRonR.xlsx"))

global_patient_skewing_summary_filtered_bk <- global_patient_skewing_summary_filtered

# compare zscaled values
my_comparisons <- list(
  c("BTHAL", "WAS"),
  c("MLD", "BTHAL"),
  c("MLD", "WAS")
)

# n_fun <- function(x){
#   return(data.frame(y = c(max(global_patient_skewing_summary_filtered$mean)*1.1),
#                     label = paste0("N=", length(x)) ))
# }
n_fun <- function(x){
  return(data.frame(y = c(min(global_patient_skewing_summary_filtered$mean)*1.1),
                    label = paste0("N=", length(x)) ))
}
perc_starting_max_stats <- 20
perc_spacer_stats <- 2.2
plot_stats_multiuni_zscaled_anova2 <- ggplot(data = global_patient_skewing_summary_filtered, aes(x = Study, y = median), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  geom_boxplot(alpha = .5, fill = "white", outlier.size = 0) +
  # geom_violin() +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size = 5) +
  geom_jitter(aes(colour = PlotLabel), height = 0, width = 0.1, alpha = 0.7, size = 3) +
  facet_grid(Tissue ~ PlotLabel) +
  scale_y_continuous(limits = c(min(global_patient_skewing_summary_filtered$mean)*1.1, max(global_patient_skewing_summary_filtered$mean)*1.5)) +
  # geom_point(size=4, alpha = .7, position = "jitter") +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(1, 1, 2, 2))+
  stat_compare_means(comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
  # stat_compare_means(comparisons = my_comparisons, label.y = seq(from = 2.5, to = 1.5, by = -0.5), label = "p.signif", p.adjust.method = "fdr") +
  theme_bw() +
  # stat_compare_means(label.y = 100) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16, angle = 0), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Lineage commitment over time"), 
       x = "Disease", y = "Percentage scaled Z-score", colour = "Lineage", fill = "Lineage"
  )

pdf(file = paste0("analyses/multi_uni_sharing/202112/", prefix_trials_outname, ".patient_skewing_summary_filtered.CRonR.zscore.pdf"), height=8, width=11)
print(plot_stats_multiuni_zscaled_anova2)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", prefix_trials_outname, ".patient_skewing_summary_filtered.CRonR.zscore.png"), height=8, width=11, units = "in", res = 300)
print(plot_stats_multiuni_zscaled_anova2)
dev.off()



## find among lineages, global and then by tissue
perc_spacer_stats <- -3
perc_starting_max_stats <- 24

n_fun <- function(x){
  return(data.frame(y = c(min(global_patient_skewing_summary_filtered$mean)*1.1),
                    label = paste0("N=", length(x)) ))
}

my_comparisons <- list(
  c("0-2", "2-15"),
  c("2-15", "30+")
)
plot_stats_multiuni_zscaled_anova3_p1 <- ggplot(data = global_patient_skewing_summary_filtered, aes(x = AgeGroup_Sorted, y = median), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  geom_boxplot(alpha = .5, fill = "white", outlier.size = 0) +
  # geom_violin() +
  # stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size = 5) +
  geom_jitter(aes(colour = PlotLabel), height = 0, width = 0.1, alpha = 0.7, size = 3) +
  facet_grid(Study ~ PlotLabel) +
  scale_y_continuous(limits = c(min(global_patient_skewing_summary_filtered$mean)*1.1, max(global_patient_skewing_summary_filtered$mean)*1.5)) +
  # geom_point(size=4, alpha = .7, position = "jitter") +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(1, 1, 2, 2))+
  stat_compare_means(aes(group = AgeGroup), comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
  # stat_compare_means(comparisons = my_comparisons, label.y = seq(from = 2.5, to = 1.5, by = -0.5), label = "p.signif", p.adjust.method = "fdr") +
  theme_bw() +
  # stat_compare_means(label.y = 100) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16, angle = 0), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Lineage commitment over time"), 
       x = "Age", y = "Percentage scaled Z-score", colour = "Lineage", fill = "Lineage"
  )

my_comparisons <- list(
  c("0-2", "2-15")
)
plot_stats_multiuni_zscaled_anova3_p2 <- ggplot(data = global_patient_skewing_summary_filtered, aes(x = AgeGroup_Sorted, y = median), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  geom_boxplot(alpha = .5, fill = "white", outlier.size = 0) +
  # geom_violin() +
  # stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size = 5) +
  geom_jitter(aes(colour = PlotLabel), height = 0, width = 0.1, alpha = 0.7, size = 3) +
  facet_grid(Study ~ PlotLabel) +
  scale_y_continuous(limits = c(min(global_patient_skewing_summary_filtered$mean)*1.1, max(global_patient_skewing_summary_filtered$mean)*1.5)) +
  # geom_point(size=4, alpha = .7, position = "jitter") +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(1, 1, 2, 2))+
  stat_compare_means(aes(group = AgeGroup), comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
  # stat_compare_means(comparisons = my_comparisons, label.y = seq(from = 2.5, to = 1.5, by = -0.5), label = "p.signif", p.adjust.method = "fdr") +
  theme_bw() +
  # stat_compare_means(label.y = 100) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16, angle = 0), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Lineage commitment over time"), 
       x = "Age", y = "Percentage scaled Z-score", colour = "Lineage", fill = "Lineage"
  )

pdf(file = paste0("analyses/multi_uni_sharing/202112/", prefix_trials_outname, ".patient_skewing_summary_filtered.CRonR.zscore.age.median.p1.pdf"), height=8, width=11)
print(plot_stats_multiuni_zscaled_anova3_p1)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", prefix_trials_outname, ".patient_skewing_summary_filtered.CRonR.zscore.age.median.p1.png"), height=8, width=11, units = "in", res = 300)
print(plot_stats_multiuni_zscaled_anova3_p1)
dev.off()
pdf(file = paste0("analyses/multi_uni_sharing/202112/", prefix_trials_outname, ".patient_skewing_summary_filtered.CRonR.zscore.age.median.p2.pdf"), height=8, width=11)
print(plot_stats_multiuni_zscaled_anova3_p2)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", prefix_trials_outname, ".patient_skewing_summary_filtered.CRonR.zscore.age.median.p2.png"), height=8, width=11, units = "in", res = 300)
print(plot_stats_multiuni_zscaled_anova3_p2)
dev.off()

## now by tissue
my_comparisons <- list(
  c("0-2", "2-15"),
  c("2-15", "30+")
)
plot_stats_multiuni_zscaled_anova4_p1 <- ggplot(data = global_patient_skewing_summary_filtered, aes(x = AgeGroup_Sorted, y = median), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  geom_boxplot(alpha = .5, fill = "white", outlier.size = 0) +
  # geom_violin() +
  # stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size = 5) +
  geom_jitter(aes(colour = PlotLabel), height = 0, width = 0.1, alpha = 0.7, size = 3) +
  facet_grid(PlotLabel + Tissue ~ Study) +
  scale_y_continuous(limits = c(min(global_patient_skewing_summary_filtered$mean)*1.1, max(global_patient_skewing_summary_filtered$mean)*1.5)) +
  # geom_point(size=4, alpha = .7, position = "jitter") +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(1, 1, 2, 2))+
  stat_compare_means(aes(group = AgeGroup), comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
  # stat_compare_means(comparisons = my_comparisons, label.y = seq(from = 2.5, to = 1.5, by = -0.5), label = "p.signif", p.adjust.method = "fdr") +
  theme_bw() +
  # stat_compare_means(label.y = 100) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16, angle = 0), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Lineage commitment over time"), 
       x = "Age", y = "Percentage scaled Z-score", colour = "Lineage", fill = "Lineage"
  )

my_comparisons <- list(
  c("0-2", "2-15")
)
plot_stats_multiuni_zscaled_anova4_p2 <- ggplot(data = global_patient_skewing_summary_filtered, aes(x = AgeGroup_Sorted, y = median), na.rm = T, se = TRUE) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  geom_boxplot(alpha = .5, fill = "white", outlier.size = 0) +
  # geom_violin() +
  # stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, size = 5) +
  geom_jitter(aes(colour = PlotLabel), height = 0, width = 0.1, alpha = 0.7, size = 3) +
  facet_grid(PlotLabel + Tissue ~ Study) +
  scale_y_continuous(limits = c(min(global_patient_skewing_summary_filtered$mean)*1.1, max(global_patient_skewing_summary_filtered$mean)*1.5)) +
  # geom_point(size=4, alpha = .7, position = "jitter") +
  # stat_compare_means(comparisons = my_comparisons, label.y = c(1, 1, 2, 2))+
  stat_compare_means(aes(group = AgeGroup), comparisons = my_comparisons, label.y = seq(from = perc_starting_max_stats, to = (perc_starting_max_stats+(length(my_comparisons)*perc_spacer_stats)), by = perc_spacer_stats)*10^-1, p.adjust.method = "fdr", label = "p.signif") +
  # stat_compare_means(comparisons = my_comparisons, label.y = seq(from = 2.5, to = 1.5, by = -0.5), label = "p.signif", p.adjust.method = "fdr") +
  theme_bw() +
  # stat_compare_means(label.y = 100) +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
  theme(axis.text.x = element_text(size=16, angle = 0), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Lineage commitment over time"), 
       x = "Age", y = "Percentage scaled Z-score", colour = "Lineage", fill = "Lineage"
  )

pdf(file = paste0("analyses/multi_uni_sharing/202112/", prefix_trials_outname, ".patient_skewing_summary_filtered.CRonR.zscore.age-tissue.median.p1.v2.pdf"), height=16, width=8)
print(plot_stats_multiuni_zscaled_anova4_p1)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", prefix_trials_outname, ".patient_skewing_summary_filtered.CRonR.zscore.age-tissue.median.p1.v2.png"), height=16, width=8, units = "in", res = 300)
print(plot_stats_multiuni_zscaled_anova4_p1)
dev.off()
pdf(file = paste0("analyses/multi_uni_sharing/202112/", prefix_trials_outname, ".patient_skewing_summary_filtered.CRonR.zscore.age-tissue.median.p2.v2.pdf"), height=16, width=8)
print(plot_stats_multiuni_zscaled_anova4_p2)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", prefix_trials_outname, ".patient_skewing_summary_filtered.CRonR.zscore.age-tissue.median.p2.v2.png"), height=16, width=8, units = "in", res = 300)
print(plot_stats_multiuni_zscaled_anova4_p2)
dev.off()






# ------ global results ----------
faceting_order <- c("MLD", "WAS", "BTHAL")
lineage_faceting_order <- c("Multilineage", "Erythroid", "Myeloid", "B", "T")
trials_colors <- c("darkblue", "forestgreen", "firebrick")
prefix_trials_outname <- paste(faceting_order, collapse = "_")
results_folder_name <- "multi_uni_sharing/"
out_file_infix <- ".hlfu_34MyBTEry"
prefix_dest_folder <- "AllPts_"
suffix_plotfile <- "" # ".noTransfDep"
followup_limit_crosstrial <- 60 # plot up to this month after GT
excel_bagcases_foldername <- "/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial MLD/analyses/13.hsc_profiles/"
excel_bagcases_filename <- "MLD.hlfu_34MyBTEry_noLy_no34exclusive.flag_matrix.bagcases.xlsx" # "MLD.hlfu_34MyBTEry.flag_matrix.bagcases.xlsx"  # "MLD.hlfu_34PBMyBT.flag_matrix.bagcases.xlsx" # MLD.hlfu_34MyBTEry
color_schema_df <- read.xlsx(xlsxFile = paste(excel_bagcases_foldername, excel_bagcases_filename, sep = ""), sheet = "color_schema")
color_schema_df_sortedbyname <- color_schema_df[order(color_schema_df$ClassName),]
color_schema_df_sortedbyplotlabel <- color_schema_df[order(color_schema_df$Label),]
scale_color_manual_colors_sortedbyname <- as.character(color_schema_df_sortedbyname[which(color_schema_df_sortedbyname$ClassName %in% levels(factor(global_allstudies_profile_skewing_fulldata_notNA_stats$Class))), "colorcode"])
scale_color_manual_colors_sortedbyplotlabel <- as.character(color_schema_df_sortedbyplotlabel[which(color_schema_df_sortedbyplotlabel$PlotLabel %in% levels(factor(global_allstudies_profile_skewing_fulldata_notNA_stats$PlotLabel))), "colorcode"])
write.xlsx(x = global_allstudies_profile_skewing_fulldata_notNA_stats, file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.global_allstudies_profile_skewing_fulldata_notNA_stats.xlsx"))


plot_CRonR_allFU_nopoints_global_identicalFU <- ggplot(global_allstudies_profile_skewing_fulldata_notNA_stats, aes(x = FollowUp, y = CRonR, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_smooth(method = "loess", formula = y ~ log(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(global_allstudies_profile_skewing_fulldata_notNA_stats$FollowUp, na.rm = T), 12), limits = c(0, followup_limit_crosstrial)) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(Tissue ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), 
       subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ),
       x = "Months after gene therapy", y = "% of IS", colour = "Class", fill = "Class")

pdf(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonR.nopoints.maxFU", followup_limit_crosstrial, ".v2.pdf"), height=8, width=10)
print(plot_CRonR_allFU_nopoints_global_identicalFU)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonR.nopoints.maxFU", followup_limit_crosstrial, ".v2.png"), height=8, width=10, units = "in", res = 300)
print(plot_CRonR_allFU_nopoints_global_identicalFU)
dev.off()

plot_CRonR_allFU_nopoints_global <- ggplot(global_allstudies_profile_skewing_fulldata_notNA_stats, aes(x = FollowUp, y = CRonR, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_smooth(method = "loess", formula = y ~ log(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(global_allstudies_profile_skewing_fulldata_notNA_stats$FollowUp, na.rm = T), 12) ) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(Tissue ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), 
       subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ),
       x = "Months after gene therapy", y = "% of IS", colour = "Class", fill = "Class")

pdf(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonR.nopoints.v2.pdf"), height=8, width=10)
print(plot_CRonR_allFU_nopoints_global)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonR.nopoints.v2.png"), height=8, width=10, units = "in", res = 300)
print(plot_CRonR_allFU_nopoints_global)
dev.off()

plot_CRonR_allFU_nopoints_global_byage_PB <- ggplot(global_allstudies_profile_skewing_fulldata_notNA_stats[which(global_allstudies_profile_skewing_fulldata_notNA_stats$Tissue == "PB"),], aes(x = FollowUp, y = CRonR, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_smooth(method = "loess", formula = y ~ log(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname[-2]) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname[-2]) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(global_allstudies_profile_skewing_fulldata_notNA_stats$FollowUp, na.rm = T), 12) ) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(AgeGroup_Sorted ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), 
       subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ),
       x = "Months after gene therapy", y = "% of IS", colour = "Class", fill = "Class")

pdf(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonR.nopoints.PB_byage.v2.pdf"), height=9, width=10)
print(plot_CRonR_allFU_nopoints_global_byage_PB)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonR.nopoints.PB_byage.v2.png"), height=9, width=10, units = "in", res = 300)
print(plot_CRonR_allFU_nopoints_global_byage_PB)
dev.off()

plot_CRonR_allFU_nopoints_global_byage_PB_identicalFU <- ggplot(global_allstudies_profile_skewing_fulldata_notNA_stats[which(global_allstudies_profile_skewing_fulldata_notNA_stats$Tissue == "PB"),], aes(x = FollowUp, y = CRonR, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_smooth(method = "loess", formula = y ~ log(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname[-2]) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname[-2]) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(global_allstudies_profile_skewing_fulldata_notNA_stats$FollowUp, na.rm = T), 12), limits = c(0, followup_limit_crosstrial)) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(AgeGroup_Sorted ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage, PB"), 
       subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ),
       x = "Months after gene therapy", y = "% of IS", colour = "Class", fill = "Class")

pdf(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonR.nopoints.PB_byage.maxFU", followup_limit_crosstrial, ".v2.pdf"), height=10, width=10)
print(plot_CRonR_allFU_nopoints_global_byage_PB_identicalFU)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonR.nopoints.PB_byage.maxFU", followup_limit_crosstrial, ".v2.png"), height=10, width=10, units = "in", res = 300)
print(plot_CRonR_allFU_nopoints_global_byage_PB_identicalFU)
dev.off()

plot_CRonR_allFU_nopoints_global_byage_BM_identicalFU <- ggplot(global_allstudies_profile_skewing_fulldata_notNA_stats[which(global_allstudies_profile_skewing_fulldata_notNA_stats$Tissue == "BM"),], aes(x = FollowUp, y = CRonR, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_smooth(method = "loess", formula = y ~ log(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(global_allstudies_profile_skewing_fulldata_notNA_stats$FollowUp, na.rm = T), 12), limits = c(0, followup_limit_crosstrial)) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(AgeGroup_Sorted ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage, BM"), 
       subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ),
       x = "Months after gene therapy", y = "% of IS", colour = "Class", fill = "Class")

pdf(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonR.nopoints.BM_byage.maxFU", followup_limit_crosstrial, ".v2.pdf"), height=10, width=10)
print(plot_CRonR_allFU_nopoints_global_byage_BM_identicalFU)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonR.nopoints.BM_byage.maxFU", followup_limit_crosstrial, ".v2.png"), height=10, width=10, units = "in", res = 300)
print(plot_CRonR_allFU_nopoints_global_byage_BM_identicalFU)
dev.off()

plot_CRonR_zscaled_global_allFU_nopoints_global <- ggplot(global_allstudies_profile_skewing_fulldata_notNA_stats, aes(x = FollowUp, y = CRonR_zscaled, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_smooth(method = "loess", formula = y ~ log(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 6) ) +
  facet_grid(Tissue ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), 
       subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ),
       x = "Months after gene therapy", y = "% of IS (Z-score)", 
       colour = "Class")

pdf(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonR_zscaled.nopoints.pdf"), height=8, width=12)
print(plot_CRonR_zscaled_global_allFU_nopoints_global)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonR_zscaled.nopoints.png"), height=8, width=12, units = "in", res = 300)
print(plot_CRonR_zscaled_global_allFU_nopoints_global)
dev.off()


plot_CRonU_allFU_nopoints_global <- ggplot(global_allstudies_profile_skewing_fulldata_notNA_stats, aes(x = FollowUp, y = CRonU, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_smooth(method = "loess", formula = y ~ log(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(global_allstudies_profile_skewing_fulldata_notNA_stats$FollowUp, na.rm = T), 12) ) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(Tissue ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), 
       subtitle = paste0("CR on U (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ),
       x = "Months after gene therapy", y = "% of IS", colour = "Class", fill = "Class")

pdf(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonU.nopoints.v2.pdf"), height=8, width=10)
print(plot_CRonU_allFU_nopoints_global)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonU.nopoints.v2.png"), height=8, width=10, units = "in", res = 300)
print(plot_CRonU_allFU_nopoints_global)
dev.off()

plot_CRonU_zscaled_global_allFU_nopoints_global <- ggplot(global_allstudies_profile_skewing_fulldata_notNA_stats, aes(x = FollowUp, y = CRonU_zscaled, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_smooth(method = "loess", formula = y ~ log(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 12) ) +
  facet_grid(Tissue ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), 
       subtitle = paste0("CR on U (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ),
       x = "Months after gene therapy", y = "% of IS (Z-score)", 
       colour = "Class")

pdf(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonU_zscaled.nopoints.v2.pdf"), height=8, width=10)
print(plot_CRonU_zscaled_global_allFU_nopoints_global)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonU_zscaled.nopoints.v2.png"), height=8, width=10, units = "in", res = 300)
print(plot_CRonU_zscaled_global_allFU_nopoints_global)
dev.off()



plot_CSonU_allFU_nopoints_global <- ggplot(global_allstudies_profile_skewing_fulldata_notNA_stats, aes(x = FollowUp, y = CSonU, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_smooth(method = "loess", formula = y ~ log(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(global_allstudies_profile_skewing_fulldata_notNA_stats$FollowUp, na.rm = T), 6), limits = c(0, 24)) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(. ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Singletons IS in percentage"), 
       subtitle = paste0("CS on U (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ),
       x = "Months after gene therapy", y = "% of IS", colour = "Class", fill = "Class")

pdf(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CSonU.nopoints.0-24m.pdf"), height=6, width=7)
print(plot_CSonU_allFU_nopoints_global)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CSonU.nopoints.0-24m.png"), height=6, width=7, units = "in", res = 300)
print(plot_CSonU_allFU_nopoints_global)
dev.off()

plot_CSonU_allFU_nopoints_global_byage <- ggplot(global_allstudies_profile_skewing_fulldata_notNA_stats, aes(x = FollowUp, y = CSonU, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_smooth(method = "loess", formula = y ~ log(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, max(global_allstudies_profile_skewing_fulldata_notNA_stats$FollowUp, na.rm = T), 6), limits = c(0, 18)) +
  # facet_wrap( ~ ClinicalStudy, ncol = 3) +
  facet_grid(AgeGroup_Sorted ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Singletons IS in percentage"), 
       subtitle = paste0("CS on U (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ),
       x = "Months after gene therapy", y = "% of IS", colour = "Class", fill = "Class")

pdf(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CSonU.nopoints.0-12m.byage.v2.pdf"), height=10, width=5)
print(plot_CSonU_allFU_nopoints_global_byage)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CSonU.nopoints.0-12m.byage.v2.png"), height=10, width=5, units = "in", res = 300)
print(plot_CSonU_allFU_nopoints_global_byage)
dev.off()

plot_CSonU_zscaled_global_allFU_nopoints_global <- ggplot(global_allstudies_profile_skewing_fulldata_notNA_stats, aes(x = FollowUp, y = CSonU_zscaled, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_smooth(method = "loess", formula = y ~ log(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_x_continuous(breaks = seq(0, max(allstudies_profile_skewing_fulldata_notNA$FollowUp, na.rm = T), 6) ) +
  facet_grid(Tissue ~ ClinicalStudy, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Singletons IS in percentage"), 
       subtitle = paste0("CS on U (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ),
       x = "Months after gene therapy", y = "% of IS (Z-score)", 
       colour = "Class")

pdf(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonU_zscaled.nopoints.pdf"), height=8, width=12)
print(plot_CRonU_zscaled_global_allFU_nopoints_global)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonU_zscaled.nopoints.png"), height=8, width=12, units = "in", res = 300)
print(plot_CRonU_zscaled_global_allFU_nopoints_global)
dev.off()




###############################################################
## ANNOTATIONS
###############################################################
#### feature annotations
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)

# prendi un test BE
bed_colnames <- c("chr","start","end","value1", "score", "strand")
bed_was1 <- read.csv(file = "source/WAS/bed/WAS1001.SC-FE.sampleFiltered.noDup.bed", header=F, fill=T, check.names = FALSE, sep = '\t')
names(bed_was1) <- bed_colnames
bed_mld1 <- read.csv(file = "source/MLD/bed/IS_matrix_classic_strand_specific_method_mld_202004_MLD01.no0.annotated.bed", header=F, fill=T, check.names = FALSE, sep = '\t')
names(bed_mld1) <- bed_colnames
bed_bthal1 <- read.csv(file = "source/BTHAL/bed/BTHAL001.collDate.filt.gdf.byPMTF.SS.invivo.bed", header=F, fill=T, check.names = FALSE, sep = '\t')
names(bed_bthal1) <- bed_colnames

#### WAS ----- WAS1001
peak <- readPeakFile("source/WAS/bed/WAS1001.SC-FE.sampleFiltered.noDup.bed")
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)

plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 20)

plotAvgProf2(peak = peak, 
              conf = 0.95, 
              TxDb = txdb, weightCol = "V5")

tagMatrix_binning <- getTagMatrix(peak = peak, TxDb = txdb, 
                                  upstream = 500, downstream = 500, 
                                  type = "start_site", by = "gene", 
                                  weightCol = "V5")

# upset plot
peakAnno <- annotatePeak("source/WAS/bed/WAS1001.SC-FE.sampleFiltered.noDup.bed", tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
# vennpie(peakAnno)
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie=TRUE)

pdf(file = paste0("analyses/circos/", paste(study_list, collapse = "_"), ".UpSetPlot.pdf"), height=5, width=8)
upsetplot(peakAnno)
dev.off()

#### MLD ----- MLD01
# upset plot
peakAnno <- annotatePeak("source/MLD/bed/IS_matrix_classic_strand_specific_method_mld_202004_MLD01.no0.annotated.bed", tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

pdf(file = paste0("analyses/circos/", paste(study_list, collapse = "_"), ".UpSetPlot.MLD01.pdf"), height=5, width=8)
upsetplot(peakAnno)
dev.off()

#### BTHAL ----- BTHAL001
# upset plot
peakAnno <- annotatePeak("source/BTHAL/bed/BTHAL001.collDate.filt.gdf.byPMTF.SS.invivo.bed", tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

pdf(file = paste0("analyses/circos/", paste(study_list, collapse = "_"), ".UpSetPlot.BTHAL001.pdf"), height=5, width=8)
upsetplot(peakAnno)
dev.off()

