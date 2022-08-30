###--- MLD - analysis script ---###
## Latest update: 02/02/2022
library(ISAnalytics)
library(dplyr)
library(ggplot2)

date <- "20220509"
date_month <- "202112"
project_name <- "MLD"

# FOLDER STRUCTURE
project_folder <- fs::path(fs::path_wd(), date_month)
matrix_folder <- fs::path(project_folder, "matrices")
matrix_cumulative_folder <- fs::path(matrix_folder, "cumulative")
matrix_agg_folder <- fs::path(matrix_folder, "aggregated")
matrix_purity_folder <- fs::path(matrix_folder, "purity_filtered")
report_folder <- fs::path(project_folder, "report")
cis_folder <- fs::path(project_folder, "CIS")
alluv_fold <- fs::path(project_folder, "alluvials")
sharing_fold <- fs::path(project_folder, "sharing")
plots_fold <- fs::path(project_folder, "plots")

for (fold in c(project_folder, matrix_folder, matrix_agg_folder,
                 matrix_cumulative_folder, report_folder, cis_folder, 
                 alluv_fold, sharing_fold, plots_fold, matrix_purity_folder)) {
  fs::dir_create(fold) # does nothing if folder already exists
}

# UTILITIES
source("../per_patient_matrix_save.R") # saves matrices per patient
source("../prototype_stats_summary.R")

# PRELIMINARY STEPS ----
## AF loading
af_path <- fs::path("/mnt/OracleResults/AssociationFiles", 
                    "asso.complete.220427.tsv") # latest version

### * Warning: this detailed process is required for this project in particular
### * due to project-specific issues. Workflow requires separate import for
### * SLiM and LAM
af_slim_only <- import_association_file(af_path,
                                        root = "/mnt/OracleResults", 
                                        filter_for = list(ProjectID = "MLD",
                                                          PCRMethod = "SLiM"), 
                                        report_path = report_folder)

pools_to_exclude <- c("ETPoolMLDLowVCN-2", "SLiMMLDVCN", "SLiMMLDVCN-1",
                      "MLDHE01-BTHAL001-2", 
                      "MLDHE01-BTHAL001-3", "20190128-MLD-AAPOOL35-1",
                      "20190128-MLD-AAPOOL34-1", 
                      "20211228-MLDslim68-AApool110-1", 
                      "20211228-MLDslim68-AApool111-1", 
                      "20211124-MLDB67-AAPOOL107-1",
                      "20211124-MLDA67-AAPOOL106-1")

af_slim_only <- af_slim_only %>%
  filter(!concatenatePoolIDSeqRun %in% pools_to_exclude)

slim_data <- import_parallel_Vispa2Matrices(
  af_slim_only, quantification_type = c("seqCount", "fragmentEstimate"),
  workers = 4, report_path = report_folder, mode = "AUTO"
)

slim_iss_stats <- import_Vispa2_stats(af_slim_only, 
                                      file_prefixes = c(default_iss_file_prefixes(),
                                                        "stats\\.mld"),
                                      report_path = report_folder
                                      )

af_all <- import_association_file(af_path, 
                                  root = NULL, 
                                  filter_for = list(ProjectID = "MLD"), 
                                  report_path = NULL)

lam_all_p <- "/mnt/OracleResults/Matrixes/Human/MLD-VISPA2_allLAM_seqCount_matrix_oldFormat.no0.annotated.tsv.gz"
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072*3)
lam_all <- import_single_Vispa2Matrix(lam_all_p)
lam_all <- lam_all %>%
  filter(CompleteAmplificationID != "all")
lam_with_meta <- lam_all %>%
  dplyr::left_join(af_all, by = "CompleteAmplificationID")

## These samples are from WAS (WAS-MLD mixed matrix) - we're excluding these
LAM_cAmps <- c("MLD01_CD13_BM_01", "MLD01_CD13_BM_03", "MLD01_CD13_BM_12", "MLD01_CD13_BM_18",
               "MLD01_CD14_PB_01", "MLD01_CD14_PB_03", "MLD01_CD14_PB_06", "MLD01_CD14_PB_09",
               "MLD01_CD14_PB_12", "MLD01_CD14_PB_18", "MLD01_CD15_BM_01", "MLD01_CD15_BM_03",
               "MLD01_CD15_BM_12", "MLD01_CD15_BM_18", "MLD01_CD15_PB_01", "MLD01_CD15_PB_03",
               "MLD01_CD15_PB_06", "MLD01_CD15_PB_09", "MLD01_CD15_PB_12", "MLD01_CD15_PB_18",
               "MLD01_CD19_BM_01", "MLD01_CD19_BM_03", "MLD01_CD19_BM_12", "MLD01_CD19_BM_18",
               "MLD01_CD19_PB_01", "MLD01_CD19_PB_03", "MLD01_CD19_PB_06", "MLD01_CD19_PB_09",
               "MLD01_CD19_PB_12", "MLD01_CD19_PB_18", "MLD01_CD3_BM_03", "MLD01_CD3_BM_12",
               "MLD01_CD3_BM_18", "MLD01_CD3_PB_03", "MLD01_CD3_PB_06", "MLD01_CD3_PB_09",
               "MLD01_CD3_PB_12", "MLD01_CD3_PB_18", "MLD01_CD34_BM_00", "MLD01_CD34_BM_01",
               "MLD01_CD34_BM_03", "MLD01_CD34_BM_12", "MLD01_CD34_BM_18", "MLD01_CD56_BM_01",
               "MLD01_CD56_BM_03", "MLD01_CD56_BM_12", "MLD01_CD56_BM_18", "MLD01_CD56_PB_01",
               "MLD01_CD56_PB_03", "MLD01_CD56_PB_06", "MLD01_CD56_PB_09", "MLD01_CD56_PB_12",
               "MLD01_CD56_PB_18", "MLD01_CD61_BM_01", "MLD01_CD61_BM_03", "MLD01_CD61_BM_12",
               "MLD01_CD61_BM_18", "MLD01_CFC_BM_00", "MLD01_CFC_BM_03", "MLD01_CFC_BM_12",
               "MLD01_CFC_BM_18", "MLD01_CFC_PB_06", "MLD01_CFC_PB_12", "MLD01_CFC_PB_18",
               "MLD01_GLY_BM_01", "MLD01_GLY_BM_03", "MLD01_GLY_BM_12", "MLD01_GLY_BM_18",
               "MLD01_MNC_BM_01", "MLD01_MNC_BM_03", "MLD01_MNC_BM_12", "MLD01_MNC_BM_18",
               "MLD01_MNC_PB_01", "MLD01_MNC_PB_03", "MLD01_MNC_PB_06", "MLD01_MNC_PB_09",
               "MLD01_MNC_PB_12", "MLD01_MNC_PB_18", "MLD01_WHOLE-BM_BM_01", "MLD01_WHOLE-BM_BM_03",
               "MLD01_WHOLE-BM_BM_12", "MLD01_WHOLE-BM_BM_18", "MLD01_WHOLE-PB_PB_01",
               "MLD01_WHOLE-PB_PB_03", "MLD01_WHOLE-PB_PB_06", "MLD01_WHOLE-PB_PB_09",
               "MLD01_WHOLE-PB_PB_12", "MLD01_WHOLE-PB_PB_18", "MLD02_CD13_BM_01", "MLD02_CD13_BM_03",
               "MLD02_CD13_BM_06", "MLD02_CD13_BM_12", "MLD02_CD13_BM_18", "MLD02_CD14_PB_01",
               "MLD02_CD14_PB_03", "MLD02_CD14_PB_06", "MLD02_CD14_PB_09", "MLD02_CD14_PB_12",
               "MLD02_CD14_PB_18", "MLD02_CD15_BM_01", "MLD02_CD15_BM_03", "MLD02_CD15_BM_06",
               "MLD02_CD15_BM_12", "MLD02_CD15_BM_18", "MLD02_CD15_PB_01", "MLD02_CD15_PB_03",
               "MLD02_CD15_PB_06", "MLD02_CD15_PB_09", "MLD02_CD15_PB_12", "MLD02_CD15_PB_18",
               "MLD02_CD19_BM_01", "MLD02_CD19_BM_03", "MLD02_CD19_BM_06", "MLD02_CD19_BM_12",
               "MLD02_CD19_BM_18", "MLD02_CD19_PB_01", "MLD02_CD19_PB_03", "MLD02_CD19_PB_06",
               "MLD02_CD19_PB_09", "MLD02_CD19_PB_12", "MLD02_CD19_PB_18", "MLD02_CD3_BM_01",
               "MLD02_CD3_BM_03", "MLD02_CD3_BM_06", "MLD02_CD3_BM_12", "MLD02_CD3_PB_01",
               "MLD02_CD3_PB_03", "MLD02_CD3_PB_06", "MLD02_CD3_PB_09", "MLD02_CD3_PB_12",
               "MLD02_CD3_PB_18", "MLD02_CD34_BM_00", "MLD02_CD34_BM_01", "MLD02_CD34_BM_03",
               "MLD02_CD34_BM_06", "MLD02_CD34_BM_12", "MLD02_CD34_BM_18", "MLD02_CD56_BM_01",
               "MLD02_CD56_BM_03", "MLD02_CD56_BM_06", "MLD02_CD56_BM_18", "MLD02_CD56_PB_01",
               "MLD02_CD56_PB_03", "MLD02_CD56_PB_06", "MLD02_CD56_PB_09", "MLD02_CD56_PB_12",
               "MLD02_CD61_BM_01", "MLD02_CD61_BM_03", "MLD02_CD61_BM_06", "MLD02_CD61_BM_12",
               "MLD02_CD61_BM_18", "MLD02_CFC_BM_00", "MLD02_CFC_BM_01", "MLD02_CFC_BM_03",
               "MLD02_CFC_BM_06", "MLD02_CFC_BM_12", "MLD02_CFC_BM_18", "MLD02_CFC_PB_03",
               "MLD02_CFC_PB_06", "MLD02_CFC_PB_09", "MLD02_CFC_PB_12", "MLD02_CFC_PB_18",
               "MLD02_GLY_BM_01", "MLD02_GLY_BM_03", "MLD02_GLY_BM_06", "MLD02_GLY_BM_12",
               "MLD02_GLY_BM_18", "MLD02_MNC_BM_01", "MLD02_MNC_BM_03", "MLD02_MNC_BM_06",
               "MLD02_MNC_BM_12", "MLD02_MNC_BM_18", "MLD02_MNC_PB_01", "MLD02_MNC_PB_03",
               "MLD02_MNC_PB_06", "MLD02_MNC_PB_09", "MLD02_MNC_PB_12",
               "MLD02_MNC_PB_18", "MLD02_WHOLE-BM_BM_01", "MLD02_WHOLE-BM_BM_03",
               "MLD02_WHOLE-BM_BM_06", "MLD02_WHOLE-BM_BM_12", "MLD02_WHOLE-BM_BM_18",
               "MLD02_WHOLE-PB_PB_01", "MLD02_WHOLE-PB_PB_03", "MLD02_WHOLE-PB_PB_06",
               "MLD02_WHOLE-PB_PB_09",
               "MLD02_WHOLE-PB_PB_12", "MLD02_WHOLE-PB_PB_18", "MLD03_CD13_BM_01",
               "MLD03_CD13_BM_03", "MLD03_CD13_BM_06", "MLD03_CD13_BM_12", "MLD03_CD14_PB_01",
               "MLD03_CD14_PB_03", "MLD03_CD14_PB_06", "MLD03_CD14_PB_09", "MLD03_CD14_PB_12",
               "MLD03_CD15_BM_01", "MLD03_CD15_BM_03", "MLD03_CD15_BM_06", "MLD03_CD15_BM_12",
               "MLD03_CD15_PB_01", "MLD03_CD15_PB_03", "MLD03_CD15_PB_06", "MLD03_CD15_PB_09",
               "MLD03_CD19_BM_01", "MLD03_CD19_BM_03", "MLD03_CD19_BM_06", "MLD03_CD19_BM_12",
               "MLD03_CD19_PB_01", "MLD03_CD19_PB_03", "MLD03_CD19_PB_06", "MLD03_CD19_PB_09",
               "MLD03_CD19_PB_12", "MLD03_CD3_BM_01", "MLD03_CD3_BM_03", "MLD03_CD3_BM_06",
               "MLD03_CD3_BM_12", "MLD03_CD3_PB_01", "MLD03_CD3_PB_03", "MLD03_CD3_PB_06",
               "MLD03_CD3_PB_09", "MLD03_CD3_PB_12", "MLD03_CD34_BM_00", "MLD03_CD34_BM_01",
               "MLD03_CD34_BM_03", "MLD03_CD34_BM_06", "MLD03_CD34_BM_12", "MLD03_CD56_BM_01",
               "MLD03_CD56_BM_03", "MLD03_CD56_BM_06", "MLD03_CD56_BM_12", "MLD03_CD56_PB_01",
               "MLD03_CD56_PB_03", "MLD03_CD56_PB_06", "MLD03_CD56_PB_09", "MLD03_CD56_PB_12",
               "MLD03_CD61_BM_01", "MLD03_CD61_BM_03", "MLD03_CD61_BM_06", "MLD03_CD61_BM_12",
               "MLD03_CFC_BM_00", "MLD03_CFC_BM_03", "MLD03_CFC_BM_12", "MLD03_CFC_PB_12",
               "MLD03_GLY_BM_01", "MLD03_GLY_BM_03", "MLD03_GLY_BM_06", "MLD03_GLY_BM_12",
               "MLD03_MNC_BM_01", "MLD03_MNC_BM_03", "MLD03_MNC_BM_06", "MLD03_MNC_BM_12",
               "MLD03_MNC_PB_01", "MLD03_MNC_PB_03", "MLD03_MNC_PB_06", "MLD03_MNC_PB_09",
               "MLD03_MNC_PB_12", "MLD03_WHOLE-BM_BM_01", "MLD03_WHOLE-BM_BM_03", "MLD03_WHOLE-BM_BM_06",
               "MLD03_WHOLE-BM_BM_12", "MLD03_WHOLE-PB_PB_01", "MLD03_WHOLE-PB_PB_03",
               "MLD03_WHOLE-PB_PB_06", "MLD03_WHOLE-PB_PB_09", "MLD03_WHOLE-PB_PB_12")
lam_with_meta <- lam_with_meta %>%
  filter(!is.na(ProjectID))
### * Check all cAMPs of LAM samples are included
all(LAM_cAmps %in% lam_with_meta$CompleteAmplificationID)
lam_all <- lam_with_meta %>%
  select(all_of(c(mandatory_IS_vars(), annotation_IS_vars())), 
         CompleteAmplificationID, Value)

af_lam <- lam_with_meta %>%
  select(!all_of(c(mandatory_IS_vars(), annotation_IS_vars())), 
         -Value) %>%
  distinct()

#stats_lam_path <- "/storage/d4/workspace/MLD_ISAnalytics/stats.vispa2.MLD.202004.tsv"
stats_lam_path <- "http://172.25.39.2/projects/MLD_ISAnalytics/stats.vispa2.MLD.202004.tsv"

stats_lam <- data.table::fread(stats_lam_path, sep = "\t")
stats_lam <- stats_lam %>% mutate(TAG = stringr::str_replace(TAG, "\\.", ""))

af_lam <- af_lam %>%
  left_join(stats_lam, by = c("concatenatePoolIDSeqRun" = "POOL", "TagSequence" = "TAG"))

final_af <- af_lam %>%
  bind_rows(slim_iss_stats) %>%
  distinct()
total_matrix <- lam_all %>% rename(seqCount = "Value") %>% bind_rows(slim_data)

### ** Free some space
rm(af_lam, stats_lam, LAM_cAmps, stats_lam, lam_all, lam_with_meta, af_slim_only,
   pools_to_exclude, slim_data, slim_iss_stats, af_all, lam_all_p, fold)
gc()

agg_key <- c("SubjectID", "CellMarker", "Tissue", "TimepointMonths")
agg_key_days <- c("SubjectID", "CellMarker", "Tissue", "TimePoint")
patients <- unique(final_af$SubjectID[!outliers_removed$SubjectID %in% 
                                            c("UTR", "FB", "UT", "CEM37")])

step_1_summary <- compute_intermediate_stats(total_matrix, final_af, patients,
                                             agg_key_days, step_suffix = "input")
step_1_summary %>%
  readr::write_tsv(file = fs::path(project_folder, paste(
    date, project_name, "summary_post_import.tsv", sep = "_"
  )))

## RECALIBRATION
recalibr <- compute_near_integrations(total_matrix, file_path = report_folder, 
                                      max_workers = 10)

step_2_summary <- compute_intermediate_stats(recalibr, final_af, patients,
                                             agg_key_days, step_suffix = "post_recalibr")
step_2_summary %>%
  readr::write_tsv(file = fs::path(project_folder, paste(
    date, project_name, "summary_post_recalibr.tsv", sep = "_"
  )))

## COLLISIONS
coll <- remove_collisions(recalibr, association_file = final_af, 
                          report_path = NULL, max_workers = 15)

step_3_summary <- compute_intermediate_stats(coll, final_af, patients,
                                             agg_key_days, 
                                             step_suffix = "post_coll")
step_3_summary %>%
  readr::write_tsv(file = fs::path(project_folder, paste(
    date, project_name, "summary_post_coll.tsv", sep = "_"
  )))

## OUTLIER REMOVAL
outliers_removed <- outlier_filter(final_af, report_path = report_folder)

final_matrix <- coll %>%
  semi_join(outliers_removed, by = "CompleteAmplificationID") %>%
  distinct()

step_4_summary <- compute_intermediate_stats(final_matrix, outliers_removed, 
                                             patients,
                                             agg_key_days, 
                                             step_suffix = "post_outliers")
step_4_summary %>%
  readr::write_tsv(file = fs::path(project_folder, paste(
    date, project_name, "summary_post_outliers.tsv", sep = "_"
  )))

## DUPLICATE REMOVAL - Removing possible duplicates between LAM and SLiM
grouping_keys <- c(
  "SubjectID",
  "CellMarker",
  "Tissue",
  "TimePoint",
  "PCRMethod"
)

joint <- final_matrix %>%
  left_join(outliers_removed, by = "CompleteAmplificationID") %>%
  select(all_of(c("CompleteAmplificationID", grouping_keys)))

to_keep <- joint %>% group_by(across(all_of(grouping_keys[grouping_keys != "PCRMethod"]))) %>% 
  distinct() %>%
  group_modify(~ {
    if (length(unique(.x$PCRMethod)) == 1) {
      return(.x)
    } else {
      return(.x %>% filter(PCRMethod == "SLiM"))
    }
  }) %>%
  ungroup()

final_matrix <- final_matrix %>% 
  filter(CompleteAmplificationID %in% to_keep$CompleteAmplificationID)

### LAM samples always lack fragment estimate: replacing missing values with
### sequence count value (mixed quantification)
final_matrix <- final_matrix %>%
  mutate(fragmentEstimate = if_else(is.na(fragmentEstimate),
                                    seqCount, fragmentEstimate))

step_5_summary <- compute_intermediate_stats(final_matrix, outliers_removed, 
                                             patients,
                                             agg_key_days, 
                                             step_suffix = "post_rem_dupl")
step_5_summary %>%
  readr::write_tsv(file = fs::path(project_folder, paste(
    date, project_name, "summary_post_duplicate_removal.tsv", sep = "_"
  )))

overall_summary <- step_1_summary %>%
  dplyr::full_join(step_2_summary, by = "SubjectID") %>%
  dplyr::full_join(step_3_summary, by = "SubjectID") %>%
  dplyr::full_join(step_4_summary, by = "SubjectID") %>%
  dplyr::full_join(step_5_summary, by = "SubjectID")

overall_summary %>%
  readr::write_tsv(file = fs::path(project_folder, paste(
    date, project_name, "summary_overall.tsv", sep = "_"
  )))

### **Free some space
rm(joint, to_keep, final_af)
gc()

# AGGREGATION
concat_agg_key <- paste0(agg_key, collapse = "_")
agg <- aggregate_values_by_key(final_matrix, 
                               association_file = outliers_removed, 
                               value_cols = c("seqCount", "fragmentEstimate"),
                               key = agg_key)
agg_meta <- aggregate_metadata(outliers_removed,
                               grouping_keys = agg_key)

agg_sc_filt <- agg %>%
  filter(seqCount_sum >= 3)

agg_key_subj <- c("SubjectID", "TimepointMonths")
concat_agg_key_subj <- paste0(agg_key_subj, collapse = "_")
agg_subj_tp <- aggregate_values_by_key(final_matrix,
                                       association_file = outliers_removed, 
                                       value_cols = c("seqCount", "fragmentEstimate"),
                                       key = agg_key_subj)

agg_key_sg <- c("SubjectID", "SuperGroup", "Tissue", "TimepointMonths")
concat_agg_key_sg <- paste0(agg_key_sg, collapse = "_")
agg_sg <- aggregate_values_by_key(final_matrix,
                                   association_file = outliers_removed %>%
                                    left_join(blood_lineages_default(), by = "CellMarker"), 
                                   value_cols = c("seqCount", "fragmentEstimate"),
                                   key = agg_key_sg)

# AGG MATRIX SAVE - PER PATIENT
### * NOTE: framentEstimate matrices will contain mixed quantifications!
matrix_file_prefix <- paste(date, project_name, paste0(agg_key, collapse = "_"),
                            sep = "_")
per_patient_sparse(agg,
                   quantific = c(seqCount = "seqCount_sum",
                                 fragmentEstimate = "fragmentEstimate_sum"),
                   prefix = matrix_file_prefix,
                   dir_name = matrix_agg_folder, 
                   row_totals = TRUE)

per_patient_sparse(agg_sc_filt,
                   quantific = c(seqCount = "seqCount_sum",
                                 fragmentEstimate = "fragmentEstimate_sum"),
                   prefix = paste0(matrix_file_prefix, "_SC-filtered"),
                   dir_name = matrix_agg_folder, 
                   row_totals = TRUE)

per_patient_sparse(agg_sg,
                   quantific = c(seqCount = "seqCount_sum",
                                 fragmentEstimate = "fragmentEstimate_sum"),
                   prefix = paste(date, project_name, 
                                  concat_agg_key_sg,
                                  sep = "_"),
                   dir_name = matrix_agg_folder, 
                   row_totals = TRUE)

# CUMULATED MATRICES
cumulated <- cumulative_is(agg,
                           key = agg_key,
                           timepoint_col = "TimepointMonths",
                           keep_og_is = FALSE,
                           expand = TRUE)
cumulated_per_patient <- cumulated %>%
  group_by(SubjectID) %>%
  group_split()
for (df in cumulated_per_patient) {
  subj <- df$SubjectID[1]
  print(paste("Saving:", subj, "..."))
  file_name <- paste0(subj, ".",
                        paste(date, project_name, 
                               paste0(agg_key, collapse = "_"),
                               sep = "_"),
                        "_cumulated_is.tsv.gz")
  df %>%
    readr::write_tsv(file = fs::path(matrix_cumulative_folder, file_name), 
                     na = "")
}

rm(cumulated_per_patient, cumulated, subj, file_name, df)
gc()

cumulated_subj <- cumulative_is(agg_subj_tp,
                           key = agg_key_subj,
                           timepoint_col = "TimepointMonths",
                           keep_og_is = FALSE,
                           expand = TRUE)
cumulated_per_patient_subj <- cumulated_subj %>%
  group_by(SubjectID) %>%
  group_split()
for (df in cumulated_per_patient_subj) {
  subj <- df$SubjectID[1]
  print(paste("Saving:", subj, "..."))
  file_name <- paste0(subj, ".",
                      paste(date, project_name, 
                            concat_agg_key_subj,
                            sep = "_"),
                      "_cumulated_is.tsv.gz")
  df %>%
    readr::write_tsv(file = fs::path(matrix_cumulative_folder, file_name), 
                     na = "")
}

cumulated_counts <- cumulative_count_union(agg, timepoint_column = "TimepointMonths",
                                           key = agg_key, zero = "00")

cumulated_counts %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         concat_agg_key,
                                         "cumulated_counts.tsv",
                                         sep = "_")))

cumulated_counts_subj <- cumulative_count_union(agg_subj_tp, 
                                                timepoint_column = "TimepointMonths",
                                                key = agg_key_subj, zero = "00")

cumulated_counts_subj %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         concat_agg_key_subj,
                                         "cumulated_counts.tsv",
                                         sep = "_")))

cumulated_from_12 <- cumulative_is(agg %>%
                                     filter(as.numeric(TimepointMonths) >= 12),
                           key = agg_key,
                           timepoint_col = "TimepointMonths",
                           keep_og_is = FALSE,
                           expand = TRUE)

cumulated_per_patient_from_12 <- cumulated_from_12 %>%
  group_by(SubjectID) %>%
  group_split()
for (df in cumulated_per_patient_from_12) {
  subj <- df$SubjectID[1]
  print(paste("Saving:", subj, "..."))
  file_name <- paste0(subj, ".",
                      paste(date, project_name, 
                            paste0(agg_key, collapse = "_"),
                            sep = "_"),
                      "_cumulated_is_FROM-TP-12.tsv.gz")
  df %>%
    readr::write_tsv(file = fs::path(matrix_cumulative_folder, file_name), 
                     na = "")
}

cumulated_counts_from_12 <- cumulative_count_union(agg %>%
                                                     filter(as.numeric(TimepointMonths) >= 12), 
                                                   timepoint_column = "TimepointMonths",
                                                  key = agg_key, zero = "00")

cumulated_counts_from_12 %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         concat_agg_key,
                                         "cumulated_counts_FROM-TP-12.tsv",
                                         sep = "_")))

cumulated_counts_subj_from_12 <- cumulative_count_union(agg_subj_tp %>%
                                                          filter(as.numeric(TimepointMonths) >= 12), 
                                                timepoint_column = "TimepointMonths",
                                                key = agg_key_subj, zero = "00")

cumulated_counts_subj_from_12 %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         concat_agg_key_subj,
                                         "cumulated_counts_FROM-TP-12.tsv",
                                         sep = "_")))

cumulated_subj_from_12 <- cumulative_is(agg_subj_tp %>%
                                          filter(as.numeric(TimepointMonths) >= 12),
                                key = agg_key_subj,
                                timepoint_col = "TimepointMonths",
                                keep_og_is = FALSE,
                                expand = TRUE)
cumulated_per_patient_subj_from_12  <- cumulated_subj_from_12 %>%
  group_by(SubjectID) %>%
  group_split()
for (df in cumulated_per_patient_subj_from_12) {
  subj <- df$SubjectID[1]
  print(paste("Saving:", subj, "..."))
  file_name <- paste0(subj, ".",
                      paste(date, project_name, 
                            concat_agg_key_subj,
                            sep = "_"),
                      "_cumulated_is_FROM-TP-12.tsv.gz")
  df %>%
    readr::write_tsv(file = fs::path(matrix_cumulative_folder, file_name), 
                     na = "")
}

# DESCRIPTIVE STATS
## Single PCR
stats_file_prefix <- paste(date, project_name, "descriptive-stats", sep = "_")
single_pcr_stats <- sample_statistics(final_matrix, 
                                      metadata = outliers_removed,
                                      value_columns = c("seqCount", "fragmentEstimate"))
single_pcr_stats$metadata %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste0(stats_file_prefix, "_single-pcr.tsv.gz")),
                   na = "")

## Aggregated
### ** Note: temporary procedure to add PCRMethod, will be changed in ISAnalytics in
### ** future update
# concat_pcr_methods <- outliers_removed %>%
#   select(all_of(c(agg_key, "PCRMethod", "NGSTechnology"))) %>%
#   group_by(across(all_of(agg_key))) %>%
#   summarise(PCRMethod = paste0(unique(PCRMethod), collapse = "|"),
#             NGSTechnology = paste0(unique(NGSTechnology), collapse = "|"), 
#             .groups = "drop")

agg_stats <- sample_statistics(agg,
                               agg_meta,
                               sample_key = agg_key,
                               value_columns = c("seqCount_sum",
                                                 "fragmentEstimate_sum"))
# agg_stats <- agg_stats$metadata %>%
#   left_join(concat_pcr_methods, by = agg_key)
agg_stats$metadata %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste0(stats_file_prefix, 
                                          "_",
                                          concat_agg_key,
                                          ".tsv.gz")),
                   na = "")

agg_stats_sc_filt <- sample_statistics(agg_sc_filt,
                                       agg_meta,
                                       sample_key = agg_key,
                                       value_columns = c("seqCount_sum",
                                                         "fragmentEstimate_sum"))
agg_stats_sc_filt$metadata %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste0(stats_file_prefix, 
                                          "_",
                                          concat_agg_key,
                                          "SC-filtered.tsv.gz")),
                   na = "")

# ABUNDANCE
abund <- compute_abundance(agg, key = agg_key)
readr::write_tsv(abund, 
                 file = fs::path(project_folder, 
                                 paste0(paste(date, project_name,
                                        "abundance", 
                                        concat_agg_key,
                                        sep = "_"), ".tsv.gz")
                                 ),
                 na = "")

# CIS (in vivo)
cis <- CIS_grubbs(agg %>% 
                    filter(TimepointMonths != "00"), 
                  by = "SubjectID")
purrr::walk2(cis[patients], names(cis[patients]), ~ {
  filename <- paste0(.y, ".", paste(date, project_name, 
                                    "CISGrubbs_inVivo.tsv.gz",
                                    sep = "_"))
  readr::write_tsv(.x, file = fs::path(cis_folder, filename), na = "")
})

volcanos <- purrr::map2(cis[patients], names(cis[patients]), ~ {
  plot_title <- paste0("MLD - ", .y, "\n")
  CIS_volcano_plot(.x, title_prefix = plot_title)
})

purrr::walk2(volcanos, names(volcanos), ~ {
  prefix <- .y
  filename <- paste0(.y, ".", paste(date, project_name, 
                                    "volcano-plot",
                                    sep = "_"))
  ggplot2::ggsave(plot = .x, filename = paste0(filename, ".pdf"), 
                  path = cis_folder, 
                  width = 8, height = 6)
  ggplot2::ggsave(plot = .x, filename = paste0(filename, ".png"), 
                  path = cis_folder, 
                  width = 8, height = 6)
})

# ALLUVIALS
alluv <- integration_alluvial_plot(x = abund %>%
                                        filter(SubjectID %in% patients), 
                                   top_abundant_tbl = TRUE, 
                                   plot_x = "TimepointMonths")
purrr::walk2(alluv, names(alluv), ~ {
  prettified_id <- stringr::str_replace_all(.y, "_", ", ")
  p <- .x$plot +
    theme_bw() +
    theme(legend.position = "none") +
    ggplot2::labs(title = paste0("MLD - IS Abundance over time - ",
                                 prettified_id),
                  subtitle = paste0("Colored flows indicate IS having ",
                                    "a relative abundance >= 1%\n",
                                    "in at least 1 time point"),
                  x = "Months after GT",
                  y = "Abundace (%)")
  file_prefix <- paste0(.y, ".",
                        paste(date, project_name, sep = "_"))
  ggplot2::ggsave(plot = p, path = alluv_fold, 
                  filename = paste0(file_prefix, "_alluvial-plot.pdf"), 
                  width = 8, height = 5)
  ggplot2::ggsave(plot = p, path = alluv_fold, 
                  filename = paste0(file_prefix, "_alluvial-plot.png"), 
                  width = 8, height = 5)
  pdf(file = fs::path(alluv_fold, 
                      paste0(file_prefix, "_top-abund-tbl.pdf")), 
      width = 13, height = 7)
  gridExtra::grid.arrange(.x$tables)
  dev.off()
})

# TOP ABUNDANT IS/TOP TARGETED GENES
top20_abund_agg <- top_integrations(abund, key = agg_key)
top20_abund_agg %>%
  readr::write_tsv(file = fs::path(project_folder, 
                                   paste0(paste(date, project_name,
                                         "top20-abundant-iss",
                                         concat_agg_key, sep = "_"),
                                         ".tsv.gz"
                                   )),
                   na = "")

n_is_by_gene_agg <- agg %>%
  group_by(across(all_of(c(annotation_IS_vars(), agg_key)))) %>%
  summarise(n_IS = n_distinct(chr, integration_locus, strand), 
            .groups = "drop")

top_20_targeted_agg <- n_is_by_gene_agg  %>% 
  group_by(across(all_of(agg_key))) %>%
  arrange(desc(n_IS)) %>% 
  slice_head(n = 20) %>%
  ungroup()

top_20_targeted_agg %>%
  readr::write_tsv(file = fs::path(project_folder, 
                                   paste0(paste(date, project_name,
                                                "top20-targeted-genes",
                                                concat_agg_key, sep = "_"),
                                          ".tsv.gz"
                                   )),
                   na = "")


# SHARING
## CD34 output: for each patient, apply the purity filter and compute sharing
## between CD34+ and lineages Herithroid, Myeloid and Lymphoid.
## Purity filter is applied by timepoint only on groups of interest
agg_by_patient <- agg %>%
  filter(SubjectID %in% patients) %>%
  group_by(SubjectID) %>%
  group_split()

groups_of_interest <- c("CD34", "CD13", "CD14", "CD15",
                        "CD19", "CD3", "CD4", "CD8", 
                        "CD36", "GLY", "GLYA")

subj_names <- purrr::map_chr(agg_by_patient, ~ .x$SubjectID[1])

purity_filtered <- purrr::map(agg_by_patient, 
                              ~ purity_filter(.x, 
                                              selected_groups = groups_of_interest,
                                              aggregation_key = agg_key, 
                                              timepoint_column = "TimepointMonths"))
purity_filtered <- purity_filtered %>% purrr::set_names(subj_names)
purrr::walk2(purity_filtered, names(purity_filtered), ~ {
  .x %>%
    readr::write_tsv(file = fs::path(
      matrix_purity_folder,
      paste0(
        .y, ".",
        paste(date, project_name, "purity-filtered", sep = "_"),
        ".tsv.gz"
      )
    ), na = "")
})

cd34_output <- purrr::map2_df(purity_filtered, names(purity_filtered),
  ~ {
    temp <- .x %>% mutate(SubjectID = .y) %>%
      dplyr::left_join(blood_lineages_default(),
                       by = "CellMarker")
    cd34_only <- temp %>%
      filter(CellMarker == "CD34")
    e_l_m <- temp %>%
      filter(HematoLineage %in% c("Erythroid",
                                  "Myeloid",
                                  "Lymphoid"))
    if (nrow(cd34_only) == 0 || nrow(e_l_m) == 0) {
      return(NULL)
    }
    shared <- is_sharing(cd34_only, e_l_m,
                         group_keys = list(g1 = c("SubjectID", "SuperGroup", "Tissue"),
                                           g2 = c("SubjectID","SuperGroup", "Tissue", "TimepointMonths")),
                         keep_genomic_coord = TRUE)
    shared
  })

cd34_output %>%
  select(-is_coord) %>%
  readr::write_tsv(
    file = fs::path(sharing_fold,
                    paste(date, project_name,
                          "CD34-Output_SuperGroup.tsv.gz", sep = "_")),
    na = ""
  )

to_plot_cd34_output <- cd34_output %>%
  tidyr::separate(col = "g2", into = c("SubjectID", "CellMarker",
                                       "Tissue", "TimepointMonths"),
                  convert = TRUE) %>%
  # filter(!SubjectID %in% c("MLDC02", "MLDCO2", "MLDCUP01", "MLDCUP02",
  #                          "MLDCUP03", "MLDCUP04", "MLDCUP05", "MLDHE01",
  #                          "MLDHE02","MLDHE03")) %>%
  left_join(blood_lineages_default(), by = "CellMarker")


cd34_output_plot <- ggplot(to_plot_cd34_output,
       aes(x = TimepointMonths,
           y = on_g1, fill = CellType, group = CellType,
           color = CellType, shape = CellMarker)) +
  # stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2)
  stat_smooth(size = 2, level = 0.4, se = T) +
  labs(title = "MLD - CD34 Output all patients",
       y = "% IS Shared with CD34", x = "Months after GT") +
  theme_bw()
ggsave(plot = cd34_output_plot,
       path = sharing_fold,
       filename = paste(date, project_name, "CD34-Output-plot.pdf", sep = "_"),
       width = 9, height = 5)
ggsave(plot = cd34_output_plot,
       path = sharing_fold,
       filename = paste(date, project_name, "CD34-Output-plot.png", sep = "_"),
       width = 9, height = 5)


cd34_output_plasma <- purrr::map2_df(purity_filtered, names(purity_filtered),
  ~ {
    temp <- .x %>% mutate(SubjectID = .y) %>%
      dplyr::left_join(blood_lineages_default(),
                       by = "CellMarker")
    cd34_only <- temp %>%
      filter(CellMarker == "CD34")
    plasma <- temp %>%
      filter(CellMarker == "Plasma")
    if (nrow(cd34_only) == 0 || nrow(plasma) == 0) {
      return(NULL)
    }
    shared <- is_sharing(cd34_only, plasma,
                         group_keys = list(g1 = c("SubjectID", "SuperGroup", "Tissue"),
                                           g2 = c("SubjectID","SuperGroup", "Tissue", "TimepointMonths")),
                         keep_genomic_coord = TRUE)
    shared
  })

cd34_output_plasma %>%
  select(-is_coord) %>%
  readr::write_tsv(
    file = fs::path(sharing_fold,
                    paste(date, project_name,
                          "CD34-Output_PLASMA_SuperGroup.tsv.gz", sep = "_")),
    na = ""
  )

plasma_output <- purrr::map2_df(purity_filtered, names(purity_filtered),
  ~ {
    temp <- .x %>% mutate(SubjectID = .y) %>%
      dplyr::left_join(blood_lineages_default(),
                       by = "CellMarker")
    plasma <- temp %>%
      filter(CellMarker == "Plasma")
    e_l_m_cd34 <- temp %>%
      filter(HematoLineage %in% c("Erythroid",
                                  "Myeloid",
                                  "Lymphoid",
                                  "CD34"))
    if (nrow(plasma) == 0 || nrow(e_l_m_cd34) == 0) {
      return(NULL)
    }
    shared <- is_sharing(plasma, e_l_m_cd34,
                         group_keys = list(g1 = c("SubjectID", "SuperGroup", "Tissue"),
                                           g2 = c("SubjectID","SuperGroup", "Tissue", "TimepointMonths")),
                         keep_genomic_coord = TRUE)
    shared
  })

plasma_output %>%
  select(-is_coord) %>%
  readr::write_tsv(
    file = fs::path(sharing_fold,
                    paste(date, project_name,
                          "PLASMA-Output_SuperGroup.tsv.gz", sep = "_")),
    na = ""
  )

## Sharing CD34 vs. Myeloid, Erythro, B, T
### ** Only in vivo tp (no 0); all tp collapsed together;
### ** only PB & BM (together and separated)
cd34_vs_lineages_sepTissue <- purrr::map2_df(purity_filtered,
                                             names(purity_filtered),
  ~ {
    print(paste0("Processing ", .y, "..."))
    invivo <- .x %>%
      mutate(SubjectID = .y) %>%
      filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
      left_join(blood_lineages_default(),  by = "CellMarker")
     cd34_only <- invivo %>%
       filter(CellMarker == "CD34")
     other_lineages <- invivo %>%
       filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
     if (nrow(cd34_only) == 0 || nrow(other_lineages) == 0) {
       return(NULL)
     }
     is_sharing(cd34_only, other_lineages,
                group_key = c("SubjectID", "CellMarker", "Tissue"))
     })

cd34_vs_lineages_sepTissue %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-PB-BM_withCellMarker.tsv.gz", sep = "_"))
                   )

cd34_vs_lineages_sepTissue_ct <- purrr::map2_df(purity_filtered, 
     names(purity_filtered),
     ~ {
       print(paste0("Processing ", .y, "..."))
       invivo <- .x %>% 
         mutate(SubjectID = .y) %>%
         filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
         left_join(blood_lineages_default(),  by = "CellMarker")
       cd34_only <- invivo %>%
         filter(CellMarker == "CD34")
       other_lineages <- invivo %>%
         filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
       if (nrow(cd34_only) == 0 || nrow(other_lineages) == 0) {
         return(NULL)
       }
       is_sharing(cd34_only, other_lineages,
                  group_key = c("SubjectID", "CellType", "Tissue"), 
                  table_for_venn = TRUE)
     })
cd34_vs_lineages_sepTissue_ct %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-PB-BM_CellTypeOnly.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_sepTissue_sg <- purrr::map2_df(purity_filtered, 
      names(purity_filtered),
      ~ {
        print(paste0("Processing ", .y, "..."))
        invivo <- .x %>% 
          mutate(SubjectID = .y) %>%
          filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
          left_join(blood_lineages_default(),  by = "CellMarker")
        cd34_only <- invivo %>%
          filter(CellMarker == "CD34")
        other_lineages <- invivo %>%
          filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
        if (nrow(cd34_only) == 0 || nrow(other_lineages) == 0) {
          return(NULL)
        }
        is_sharing(cd34_only, other_lineages,
                   group_key = c("SubjectID", "SuperGroup", "Tissue"), 
                   table_for_venn = TRUE)
      })

cd34_vs_lineages_sepTissue_sg %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-PB-BM_SuperGroup.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_sepTissue_ge12 <- purrr::map2_df(purity_filtered, 
  names(purity_filtered),
  ~ {
    print(paste0("Processing ", .y, "..."))
    invivo <- .x %>% 
      mutate(SubjectID = .y) %>%
      filter(as.numeric(TimepointMonths) >= 12, Tissue %in% c("PB", "BM")) %>%
      left_join(blood_lineages_default(),  by = "CellMarker")
    cd34_only <- invivo %>%
      filter(CellMarker == "CD34")
    other_lineages <- invivo %>%
      filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
    if (nrow(cd34_only) == 0 || nrow(other_lineages) == 0) {
      return(NULL)
    }
    is_sharing(cd34_only, other_lineages,
               group_key = c("SubjectID", "CellMarker", "Tissue"), 
               table_for_venn = TRUE)
  })

cd34_vs_lineages_sepTissue_ct_ge12 <- purrr::map2_df(purity_filtered, 
   names(purity_filtered),
   ~ {
     print(paste0("Processing ", .y, "..."))
     invivo <- .x %>% 
       mutate(SubjectID = .y) %>%
       filter(as.numeric(TimepointMonths) >= 12, Tissue %in% c("PB", "BM")) %>%
       left_join(blood_lineages_default(),  by = "CellMarker")
     cd34_only <- invivo %>%
       filter(CellMarker == "CD34")
     other_lineages <- invivo %>%
       filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
     if (nrow(cd34_only) == 0 || nrow(other_lineages) == 0) {
       return(NULL)
     }
     is_sharing(cd34_only, other_lineages,
                group_key = c("SubjectID", "CellType", "Tissue"), 
                table_for_venn = TRUE)
   })

cd34_vs_lineages_sepTissue_sg_ge12 <- purrr::map2_df(purity_filtered, 
   names(purity_filtered),
   ~ {
     print(paste0("Processing ", .y, "..."))
     invivo <- .x %>% 
       mutate(SubjectID = .y) %>%
       filter(as.numeric(TimepointMonths) >= 12, Tissue %in% c("PB", "BM")) %>%
       left_join(blood_lineages_default(),  by = "CellMarker")
     cd34_only <- invivo %>%
       filter(CellMarker == "CD34")
     other_lineages <- invivo %>%
       filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
     if (nrow(cd34_only) == 0 || nrow(other_lineages) == 0) {
       return(NULL)
     }
     is_sharing(cd34_only, other_lineages,
                group_key = c("SubjectID", "SuperGroup", "Tissue"), 
                table_for_venn = TRUE)
   })

cd34_vs_lineages_sepTissue_ge12 %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-PB-BM_withCellMarker_FROM-TP-12.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_sepTissue_ct_ge12 %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-PB-BM_CellTypeOnly_FROM-TP-12.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_sepTissue_sg_ge12 %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-PB-BM_SuperGroup_FROM-TP-12.tsv.gz", sep = "_"))
  )


cd34_vs_lineages_collTissue <- purrr::map2_df(purity_filtered,
   names(purity_filtered),
   ~ {
     print(paste0("Processing ", .y, "..."))
     invivo <- .x %>%
       mutate(SubjectID = .y) %>%
       filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
       left_join(blood_lineages_default(),  by = "CellMarker")
     cd34_only <- invivo %>%
       filter(CellMarker == "CD34")
     other_lineages <- invivo %>%
       filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
     if (nrow(cd34_only) == 0 || nrow(other_lineages) == 0) {
       return(NULL)
     }
     is_sharing(cd34_only, other_lineages,
                group_keys = list(g1 = c("SubjectID", "CellMarker", "Tissue"),
                                  g2 = c("SubjectID", "CellMarker")))
   })

cd34_vs_lineages_collTissue %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-tissueCollapsed_withCellMarker.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_collTissue_ct <- purrr::map2_df(purity_filtered, 
  names(purity_filtered),
  ~ {
    print(paste0("Processing ", .y, "..."))
    invivo <- .x %>% 
      mutate(SubjectID = .y) %>%
      filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
      left_join(blood_lineages_default(),  by = "CellMarker")
    cd34_only <- invivo %>%
      filter(CellMarker == "CD34")
    other_lineages <- invivo %>%
      filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
    if (nrow(cd34_only) == 0 || nrow(other_lineages) == 0) {
      return(NULL)
    }
    is_sharing(cd34_only, other_lineages,
               group_keys = list(g1 = c("SubjectID", "CellType", "Tissue"),
                                 g2 = c("SubjectID", "CellType")), 
               table_for_venn = TRUE)
  })

cd34_vs_lineages_collTissue_ct %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-tissueCollapsed_CellTypeOnly.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_collTissue_sg <- purrr::map2_df(purity_filtered, 
     names(purity_filtered),
     ~ {
       print(paste0("Processing ", .y, "..."))
       invivo <- .x %>% 
         mutate(SubjectID = .y) %>%
         filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
         left_join(blood_lineages_default(),  by = "CellMarker")
       cd34_only <- invivo %>%
         filter(CellMarker == "CD34")
       other_lineages <- invivo %>%
         filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
       if (nrow(cd34_only) == 0 || nrow(other_lineages) == 0) {
         return(NULL)
       }
       is_sharing(cd34_only, other_lineages,
                  group_keys = list(g1 = c("SubjectID", "SuperGroup", "Tissue"),
                                    g2 = c("SubjectID", "SuperGroup")), 
                  table_for_venn = TRUE)
     })

cd34_vs_lineages_collTissue_sg %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-tissueCollapsed_SuperGroup.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_collTissue_ge12 <- purrr::map2_df(purity_filtered,
   names(purity_filtered),
   ~ {
     print(paste0("Processing ", .y, "..."))
     invivo <- .x %>%
       mutate(SubjectID = .y) %>%
       filter(as.numeric(TimepointMonths) >= 12, Tissue %in% c("PB", "BM")) %>%
       left_join(blood_lineages_default(),  by = "CellMarker")
     cd34_only <- invivo %>%
       filter(CellMarker == "CD34")
     other_lineages <- invivo %>%
       filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
     if (nrow(cd34_only) == 0 || nrow(other_lineages) == 0) {
       return(NULL)
     }
     is_sharing(cd34_only, other_lineages,
                group_keys = list(g1 = c("SubjectID", "CellMarker", "Tissue"),
                                  g2 = c("SubjectID", "CellMarker")))
   })

cd34_vs_lineages_collTissue_ct_ge12 <- purrr::map2_df(purity_filtered, 
  names(purity_filtered),
  ~ {
    print(paste0("Processing ", .y, "..."))
    invivo <- .x %>% 
      mutate(SubjectID = .y) %>%
      filter(as.numeric(TimepointMonths) >= 12, Tissue %in% c("PB", "BM")) %>%
      left_join(blood_lineages_default(),  by = "CellMarker")
    cd34_only <- invivo %>%
      filter(CellMarker == "CD34")
    other_lineages <- invivo %>%
      filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
    if (nrow(cd34_only) == 0 || nrow(other_lineages) == 0) {
      return(NULL)
    }
    is_sharing(cd34_only, other_lineages,
               group_keys = list(g1 = c("SubjectID", "CellType", "Tissue"),
                                 g2 = c("SubjectID", "CellType")), 
               table_for_venn = TRUE)
  })

cd34_vs_lineages_collTissue_sg_ge12 <- purrr::map2_df(purity_filtered, 
  names(purity_filtered),
  ~ {
    print(paste0("Processing ", .y, "..."))
    invivo <- .x %>% 
      mutate(SubjectID = .y) %>%
      filter(as.numeric(TimepointMonths) >= 12, Tissue %in% c("PB", "BM")) %>%
      left_join(blood_lineages_default(),  by = "CellMarker")
    cd34_only <- invivo %>%
      filter(CellMarker == "CD34")
    other_lineages <- invivo %>%
      filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
    if (nrow(cd34_only) == 0 || nrow(other_lineages) == 0) {
      return(NULL)
    }
    is_sharing(cd34_only, other_lineages,
               group_keys = list(g1 = c("SubjectID", "SuperGroup", "Tissue"),
                                 g2 = c("SubjectID", "SuperGroup")), 
               table_for_venn = TRUE)
  })

cd34_vs_lineages_collTissue_ge12 %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-tissueCollapsed_withCellMarker_FROM-TP-12.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_collTissue_ct_ge12 %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-tissueCollapsed_CellTypeOnly_FROM-TP-12.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_collTissue_sg_ge12 %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-tissueCollapsed_SuperGroup_FROM-TP-12.tsv.gz", sep = "_"))
  )

## Sharing WHOLE/MNC vs. Myeloid, Erythro, B, T
### ** Only in vivo tp (no 0); all tp collapsed together;
### ** only PB & BM (together and separated)
sep_whole_other <- function(df, sub) {
  invivo <- df %>% 
    mutate(SubjectID = sub) %>%
    filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
    left_join(blood_lineages_default(),  by = "CellMarker")
  whole_mnc_only <- invivo %>%
    filter(CellMarker %in% c("Whole", "MNC"))
  other_lineages <- invivo %>%
    filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
  return(list(whole = whole_mnc_only, other = other_lineages))
}

insert_id_whole <- function(ids) {
  splits <- stringr::str_split(ids, "_")
  purrr::map_chr(splits, ~ {
    .x[3] <- .x[2]
    .x[2] <- "Whole"
    paste0(.x, collapse = "_")
  })
}

whole_vs_lineages_sepTissue <- purrr::map2_df(purity_filtered,
  names(purity_filtered),
  ~ {
   print(paste0("Processing ", .y, "..."))
   sep <- sep_whole_other(.x, .y)
   if (nrow(sep$whole) == 0 || nrow(sep$other) == 0) {
     return(NULL)
   }
   is_sharing(sep$whole, sep$other,
              group_keys = list(g1 = c("SubjectID", "Tissue"),
                                g2 = c("SubjectID", "CellMarker", "Tissue"))
              ) %>%
     mutate(g1 = insert_id_whole(g1))
  })

whole_vs_lineages_sepTissue %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_withCellMarker.tsv.gz", sep = "_"))
  )

whole_vs_lineages_sepTissue_ct <- purrr::map2_df(purity_filtered, 
    names(purity_filtered),
    ~ {
      print(paste0("Processing ", .y, "..."))
      sep <- sep_whole_other(.x, .y)
      if (nrow(sep$whole) == 0 || nrow(sep$other) == 0) {
        return(NULL)
      }
      is_sharing(sep$whole, sep$other,
                 group_keys = list(g1 = c("SubjectID", "Tissue"),
                                   g2 = c("SubjectID", "CellType", "Tissue")),
                 table_for_venn = TRUE
      ) %>%
        mutate(g1 = insert_id_whole(g1))
    })

whole_vs_lineages_sepTissue_ct %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_CellTypeOnly.tsv.gz", sep = "_"))
  )

whole_vs_lineages_sepTissue_sg <- purrr::map2_df(purity_filtered, 
     names(purity_filtered),
     ~ {
       print(paste0("Processing ", .y, "..."))
       sep <- sep_whole_other(.x, .y)
       if (nrow(sep$whole) == 0 || nrow(sep$other) == 0) {
         return(NULL)
       }
       is_sharing(sep$whole, sep$other,
                  group_keys = list(g1 = c("SubjectID", "Tissue"),
                                    g2 = c("SubjectID", "SuperGroup", "Tissue")),
                  table_for_venn = TRUE
       ) %>%
         mutate(g1 = insert_id_whole(g1))
     })

whole_vs_lineages_sepTissue_sg %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_SuperGroup.tsv.gz", sep = "_"))
  )

sep_whole_other_ge12 <- function(df, sub) {
  invivo <- df %>% 
    mutate(SubjectID = sub) %>%
    filter(as.numeric(TimepointMonths) >= 12, Tissue %in% c("PB", "BM")) %>%
    left_join(blood_lineages_default(),  by = "CellMarker")
  whole_mnc_only <- invivo %>%
    filter(CellMarker %in% c("Whole", "MNC"))
  other_lineages <- invivo %>%
    filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
  return(list(whole = whole_mnc_only, other = other_lineages))
}

whole_vs_lineages_sepTissue_ge12 <- purrr::map2_df(purity_filtered,
 names(purity_filtered),
 ~ {
   print(paste0("Processing ", .y, "..."))
   sep <- sep_whole_other_ge12(.x, .y)
   if (nrow(sep$whole) == 0 || nrow(sep$other) == 0) {
     return(NULL)
   }
   is_sharing(sep$whole, sep$other,
              group_keys = list(g1 = c("SubjectID", "Tissue"),
                                g2 = c("SubjectID", "CellMarker", "Tissue"))
   ) %>%
     mutate(g1 = insert_id_whole(g1))
 })

whole_vs_lineages_sepTissue_ct_ge12 <- purrr::map2_df(purity_filtered, 
  names(purity_filtered),
  ~ {
    print(paste0("Processing ", .y, "..."))
    sep <- sep_whole_other_ge12(.x, .y)
    if (nrow(sep$whole) == 0 || nrow(sep$other) == 0) {
      return(NULL)
    }
    is_sharing(sep$whole, sep$other,
               group_keys = list(g1 = c("SubjectID", "Tissue"),
                                 g2 = c("SubjectID", "CellType", "Tissue")),
               table_for_venn = TRUE
    ) %>%
      mutate(g1 = insert_id_whole(g1))
  })

whole_vs_lineages_sepTissue_sg_ge12 <- purrr::map2_df(purity_filtered, 
  names(purity_filtered),
  ~ {
    print(paste0("Processing ", .y, "..."))
    sep <- sep_whole_other_ge12(.x, .y)
    if (nrow(sep$whole) == 0 || nrow(sep$other) == 0) {
      return(NULL)
    }
    is_sharing(sep$whole, sep$other,
               group_keys = list(g1 = c("SubjectID", "Tissue"),
                                 g2 = c("SubjectID", "SuperGroup", "Tissue")),
               table_for_venn = TRUE
    ) %>%
      mutate(g1 = insert_id_whole(g1))
  })

whole_vs_lineages_sepTissue_ge12 %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_withCellMarker_FROM-TP-12.tsv.gz", sep = "_"))
  )

whole_vs_lineages_sepTissue_ct_ge12 %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_CellTypeOnly_FROM-TP-12.tsv.gz", sep = "_"))
  )

whole_vs_lineages_sepTissue_sg_ge12 %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_SuperGroup_FROM-TP-12.tsv.gz", sep = "_"))
  )

whole_vs_lineages_collTissue <- purrr::map2_df(purity_filtered,
    names(purity_filtered),
    ~ {
      print(paste0("Processing ", .y, "..."))
      sep <- sep_whole_other(.x, .y)
      if (nrow(sep$whole) == 0 || nrow(sep$other) == 0) {
        return(NULL)
      }
      is_sharing(sep$whole, sep$other,
                 group_keys = list(g1 = c("SubjectID"),
                                   g2 = c("SubjectID", "CellMarker")),
                 table_for_venn = TRUE
      ) %>%
        mutate(g1 = paste0(g1, "_Whole"))
    })

whole_vs_lineages_collTissue %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-tissueCollapsed_withCellMarker.tsv.gz", sep = "_"))
  )

whole_vs_lineages_collTissue_ct <- purrr::map2_df(purity_filtered, 
   names(purity_filtered),
   ~ {
     print(paste0("Processing ", .y, "..."))
     sep <- sep_whole_other(.x, .y)
     if (nrow(sep$whole) == 0 || nrow(sep$other) == 0) {
       return(NULL)
     }
     is_sharing(sep$whole, sep$other,
                group_keys = list(g1 = c("SubjectID"),
                                  g2 = c("SubjectID", "CellType")),
                table_for_venn = TRUE
     ) %>%
       mutate(g1 = paste0(g1, "_Whole"))
   })

whole_vs_lineages_collTissue_ct %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-tissueCollapsed_CellTypeOnly.tsv.gz", sep = "_"))
  )

whole_vs_lineages_collTissue_sg <- purrr::map2_df(purity_filtered, 
  names(purity_filtered),
  ~ {
    print(paste0("Processing ", .y, "..."))
    sep <- sep_whole_other(.x, .y)
    if (nrow(sep$whole) == 0 || nrow(sep$other) == 0) {
      return(NULL)
    }
    is_sharing(sep$whole, sep$other,
               group_keys = list(g1 = c("SubjectID"),
                                 g2 = c("SubjectID", "SuperGroup")),
               table_for_venn = TRUE
    ) %>%
      mutate(g1 = paste0(g1, "_Whole"))
  })

whole_vs_lineages_collTissue_sg %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-tissueCollapsed_SuperGroup.tsv.gz", sep = "_"))
  )

whole_vs_lineages_collTissue_ge12 <- purrr::map2_df(purity_filtered,
  names(purity_filtered),
  ~ {
    print(paste0("Processing ", .y, "..."))
    sep <- sep_whole_other_ge12(.x, .y)
    if (nrow(sep$whole) == 0 || nrow(sep$other) == 0) {
      return(NULL)
    }
    is_sharing(sep$whole, sep$other,
               group_keys = list(g1 = c("SubjectID"),
                                 g2 = c("SubjectID", "CellMarker")),
               table_for_venn = TRUE
    ) %>%
      mutate(g1 = paste0(g1, "_Whole"))
  })

whole_vs_lineages_collTissue_ct_ge12 <- purrr::map2_df(purity_filtered, 
   names(purity_filtered),
   ~ {
     print(paste0("Processing ", .y, "..."))
     sep <- sep_whole_other_ge12(.x, .y)
     if (nrow(sep$whole) == 0 || nrow(sep$other) == 0) {
       return(NULL)
     }
     is_sharing(sep$whole, sep$other,
                group_keys = list(g1 = c("SubjectID"),
                                  g2 = c("SubjectID", "CellType")),
                table_for_venn = TRUE
     ) %>%
       mutate(g1 = paste0(g1, "_Whole"))
   })

whole_vs_lineages_collTissue_sg_ge12 <- purrr::map2_df(purity_filtered, 
   names(purity_filtered),
   ~ {
     print(paste0("Processing ", .y, "..."))
     sep <- sep_whole_other_ge12(.x, .y)
     if (nrow(sep$whole) == 0 || nrow(sep$other) == 0) {
       return(NULL)
     }
     is_sharing(sep$whole, sep$other,
                group_keys = list(g1 = c("SubjectID"),
                                  g2 = c("SubjectID", "SuperGroup")),
                table_for_venn = TRUE
     ) %>%
       mutate(g1 = paste0(g1, "_Whole"))
   })

whole_vs_lineages_collTissue_ge12 %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-tissueCollapsed_withCellMarker_FROM-TP-12.tsv.gz", sep = "_"))
  )

whole_vs_lineages_collTissue_ct_ge12 %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-tissueCollapsed_CellTypeOnly_FROM-TP-12.tsv.gz", sep = "_"))
  )

whole_vs_lineages_collTissue_sg_ge12 %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-tissueCollapsed_SuperGroup_FROM-TP-12.tsv.gz", sep = "_"))
  )


## NESTED SHARING (cd34 vs WholeMNC) vs. Myeloid, Erythro, B, T
### ** Only in vivo tp (no 0); all tp collapsed together;
### ** only PB & BM (together and separated)
cd34_vs_wholemnc_sepTissue <- purrr::map2(purity_filtered, 
  names(purity_filtered),
  ~ {
    print(paste0("Processing ", .y, "..."))
    invivo <- .x %>% 
      mutate(SubjectID = .y) %>%
      filter(TimepointMonths != "00", Tissue %in% c("PB", "BM"))
    whole_mnc_only <- invivo %>%
      filter(CellMarker %in% c("Whole", "MNC"))
    cd34_only <- invivo %>%
      filter(CellMarker == "CD34")
    if (nrow(cd34_only) == 0 || nrow(whole_mnc_only) == 0) {
      return(NULL)
    }
    is_sharing(cd34_only, whole_mnc_only,
               group_keys = list(g1 = c("SubjectID", "CellMarker", "Tissue"),
                                g2 = c("SubjectID", "Tissue")), 
               keep_genomic_coord = TRUE) %>%
      mutate(g2 = insert_id_whole(g2))
  })

# purrr::reduce(cd34_vs_wholemnc_sepTissue, bind_rows) %>% 
#   select(-is_coord) %>%
#   readr::write_tsv(
#   file = fs::path(sharing_fold,
#                   paste(date, project_name,
#                         "CD34-WholeMNC-inVivo-tpCollapsed-PB-BM.tsv.gz", sep = "_"))
# )
  
regen <- function(df) {
  if (is.null(df)) {
    return(NULL)
  }
  df <- df %>%
    tidyr::separate(col = g1,
                    into = c("SubjectID_g1", "CellMarker_g1", "Tissue_g1")) %>%
    tidyr::separate(col = g2,
                    into = c("SubjectID_g2", "CellMarker_g2", "Tissue_g2")) %>%
    tidyr::unite(SubjectID_g1, CellMarker_g1, Tissue_g1, CellMarker_g2,
                 Tissue_g2, col = "id") %>%
    select(id, is_coord) %>%
    tidyr::unnest(is_coord) %>%
    distinct()
}

shared_nested_sepTissue <- purrr::map2_df(cd34_vs_wholemnc_sepTissue,
                                          names(cd34_vs_wholemnc_sepTissue),
  ~ {
    print(paste0("Processing ", .y, "..."))
    regenerated <- regen(.x)
    if (is.null(regenerated)) {
      return(NULL)
    }
    other_lineages <- purity_filtered[[.y]] %>%
      mutate(SubjectID = .y) %>%
      filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
      left_join(blood_lineages_default(),  by = "CellMarker") %>%
      filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
    if (nrow(other_lineages) == 0) {
      return(NULL)
    }
    is_sharing(regenerated, other_lineages,
               group_keys = list(g1 = c("id"),
                                 g2 = c("SubjectID", "CellMarker", "Tissue")),
               table_for_venn = TRUE)
  })

shared_nested_sepTissue %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-PB-BM_withCellMarker.tsv.gz", sep = "_"))
  )

shared_nested_sepTissue_ct <- purrr::map2_df(cd34_vs_wholemnc_sepTissue,
                                          names(cd34_vs_wholemnc_sepTissue),
    ~ {
      print(paste0("Processing ", .y, "..."))
      regenerated <- regen(.x)
      if (is.null(regenerated)) {
        return(NULL)
      }
      other_lineages <- purity_filtered[[.y]] %>%
        mutate(SubjectID = .y) %>%
        filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
        left_join(blood_lineages_default(),  by = "CellMarker") %>%
        filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
      if (nrow(other_lineages) == 0) {
        return(NULL)
      }
      is_sharing(regenerated, other_lineages,
                 group_keys = list(g1 = c("id"),
                                   g2 = c("SubjectID", "CellType", "Tissue")),
                 table_for_venn = TRUE)
    })

shared_nested_sepTissue_ct %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-PB-BM_CellTypeOnly.tsv.gz", sep = "_"))
  )

shared_nested_sepTissue_sg <- purrr::map2_df(cd34_vs_wholemnc_sepTissue,
     names(cd34_vs_wholemnc_sepTissue),
     ~ {
       print(paste0("Processing ", .y, "..."))
       regenerated <- regen(.x)
       if (is.null(regenerated)) {
         return(NULL)
       }
       other_lineages <- purity_filtered[[.y]] %>%
         mutate(SubjectID = .y) %>%
         filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
         left_join(blood_lineages_default(),  by = "CellMarker") %>%
         filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
       if (nrow(other_lineages) == 0) {
         return(NULL)
       }
       is_sharing(regenerated, other_lineages,
                  group_keys = list(g1 = c("id"),
                                    g2 = c("SubjectID", "SuperGroup", "Tissue")),
                  table_for_venn = TRUE)
     })

shared_nested_sepTissue_sg %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-PB-BM_SuperGroup.tsv.gz", sep = "_"))
  )

cd34_vs_wholemnc_sepTissue_ge12 <- purrr::map2(purity_filtered,
   names(purity_filtered),
   ~ {
     print(paste0("Processing ", .y, "..."))
     invivo <- .x %>%
       mutate(SubjectID = .y) %>%
       filter(as.numeric(TimepointMonths) >= 12, Tissue %in% c("PB", "BM"))
     whole_mnc_only <- invivo %>%
       filter(CellMarker %in% c("Whole", "MNC"))
     cd34_only <- invivo %>%
       filter(CellMarker == "CD34")
     if (nrow(cd34_only) == 0 || nrow(whole_mnc_only) == 0) {
       return(NULL)
     }
     is_sharing(cd34_only, whole_mnc_only,
                group_keys = list(g1 = c("SubjectID", "CellMarker", "Tissue"),
                                  g2 = c("SubjectID", "Tissue")),
                keep_genomic_coord = TRUE) %>%
       mutate(g2 = insert_id_whole(g2))
   })

shared_nested_sepTissue_ge12 <- purrr::map2_df(cd34_vs_wholemnc_sepTissue_ge12,
 names(cd34_vs_wholemnc_sepTissue_ge12),
 ~ {
   print(paste0("Processing ", .y, "..."))
   regenerated <- regen(.x)
   if (is.null(regenerated)) {
     return(NULL)
   }
   other_lineages <- purity_filtered[[.y]] %>%
     mutate(SubjectID = .y) %>%
     filter(as.numeric(TimepointMonths) >= 12, Tissue %in% c("PB", "BM")) %>%
     left_join(blood_lineages_default(),  by = "CellMarker") %>%
     filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
   if (nrow(other_lineages) == 0) {
     return(NULL)
   }
   is_sharing(regenerated, other_lineages,
              group_keys = list(g1 = c("id"),
                                g2 = c("SubjectID", "CellMarker", "Tissue")),
              table_for_venn = TRUE)
 })

shared_nested_sepTissue_ct_ge12 <- purrr::map2_df(cd34_vs_wholemnc_sepTissue_ge12,
  names(cd34_vs_wholemnc_sepTissue_ge12),
  ~ {
    print(paste0("Processing ", .y, "..."))
    regenerated <- regen(.x)
    if (is.null(regenerated)) {
      return(NULL)
    }
    other_lineages <- purity_filtered[[.y]] %>%
      mutate(SubjectID = .y) %>%
      filter(as.numeric(TimepointMonths) >= 12, Tissue %in% c("PB", "BM")) %>%
      left_join(blood_lineages_default(),  by = "CellMarker") %>%
      filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
    if (nrow(other_lineages) == 0) {
      return(NULL)
    }
    is_sharing(regenerated, other_lineages,
               group_keys = list(g1 = c("id"),
                                 g2 = c("SubjectID", "CellType", "Tissue")),
               table_for_venn = TRUE)
  })

shared_nested_sepTissue_sg_ge12 <- purrr::map2_df(cd34_vs_wholemnc_sepTissue_ge12,
  names(cd34_vs_wholemnc_sepTissue_ge12),
  ~ {
    print(paste0("Processing ", .y, "..."))
    regenerated <- regen(.x)
    if (is.null(regenerated)) {
      return(NULL)
    }
    other_lineages <- purity_filtered[[.y]] %>%
      mutate(SubjectID = .y) %>%
      filter(as.numeric(TimepointMonths) >= 12, Tissue %in% c("PB", "BM")) %>%
      left_join(blood_lineages_default(),  by = "CellMarker") %>%
      filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
    if (nrow(other_lineages) == 0) {
      return(NULL)
    }
    is_sharing(regenerated, other_lineages,
               group_keys = list(g1 = c("id"),
                                 g2 = c("SubjectID", "SuperGroup", "Tissue")),
               table_for_venn = TRUE)
  })

shared_nested_sepTissue_ge12 %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-PB-BM_withCellMarker_FROM-TP-12.tsv.gz", sep = "_"))
  )

shared_nested_sepTissue_ct_ge12 %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-PB-BM_CellTypeOnly_FROM-TP-12.tsv.gz", sep = "_"))
  )

shared_nested_sepTissue_sg_ge12 %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-PB-BM_SuperGroup_FROM-TP-12.tsv.gz", sep = "_"))
  )


cd34_vs_wholemnc_collTissue <- purrr::map2(purity_filtered, 
                                          names(purity_filtered),
    ~ {
      print(paste0("Processing ", .y, "..."))
      invivo <- .x %>% 
        mutate(SubjectID = .y) %>%
        filter(TimepointMonths != "00", Tissue %in% c("PB", "BM"))
      whole_mnc_only <- invivo %>%
        filter(CellMarker %in% c("Whole", "MNC"))
      cd34_only <- invivo %>%
        filter(CellMarker == "CD34")
      if (nrow(cd34_only) == 0 || nrow(whole_mnc_only) == 0) {
        return(NULL)
      }
      is_sharing(cd34_only, whole_mnc_only,
                 group_keys = list(g1 = c("SubjectID", "CellMarker"),
                                   g2 = c("SubjectID")), 
                 keep_genomic_coord = TRUE) %>%
        mutate(g2 = paste0(g2, "_Whole"))
    })

# purrr::reduce(cd34_vs_wholemnc_collTissue, bind_rows) %>% 
#   select(-is_coord) %>%
#   readr::write_tsv(
#     file = fs::path(sharing_fold,
#                     paste(date, project_name,
#                           "CD34-WholeMNC-inVivo-tpCollapsed-tissueCollapsed.tsv.gz", sep = "_"))
#   )

regen2 <- function(df) {
  if (is.null(df)) {
    return(NULL)
  }
  df <- df %>%
    tidyr::separate(col = g1,
                    into = c("SubjectID_g1", "CellMarker_g1")) %>%
    tidyr::separate(col = g2,
                    into = c("SubjectID_g2", "CellMarker_g2")) %>%
    tidyr::unite(SubjectID_g1, CellMarker_g1, CellMarker_g2, col = "id") %>%
    select(id, is_coord) %>%
    tidyr::unnest(is_coord) %>%
    distinct()
}

shared_nested_collTissue <- purrr::map2_df(cd34_vs_wholemnc_collTissue,
                                          names(cd34_vs_wholemnc_collTissue),
  ~ {
    print(paste0("Processing ", .y, "..."))
    regenerated <- regen2(.x)
    if (is.null(regenerated)) {
      return(NULL)
    }
    other_lineages <- purity_filtered[[.y]] %>%
      mutate(SubjectID = .y) %>%
      filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
      left_join(blood_lineages_default(),  by = "CellMarker") %>%
      filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
    if (nrow(other_lineages) == 0) {
      return(NULL)
    }
    is_sharing(regenerated, other_lineages,
               group_keys = list(g1 = c("id"),
                                 g2 = c("SubjectID", "CellMarker")),
               table_for_venn = TRUE)
  })

shared_nested_collTissue %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-tissueCollapsed_withCellMarker.tsv.gz", sep = "_"))
  )

shared_nested_collTissue_ct <- purrr::map2_df(cd34_vs_wholemnc_collTissue,
                                             names(cd34_vs_wholemnc_collTissue),
   ~ {
     print(paste0("Processing ", .y, "..."))
     regenerated <- regen2(.x)
     if (is.null(regenerated)) {
       return(NULL)
     }
     other_lineages <- purity_filtered[[.y]] %>%
       mutate(SubjectID = .y) %>%
       filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
       left_join(blood_lineages_default(),  by = "CellMarker") %>%
       filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
     if (nrow(other_lineages) == 0) {
       return(NULL)
     }
     is_sharing(regenerated, other_lineages,
                group_keys = list(g1 = c("id"),
                                  g2 = c("SubjectID", "CellType")),
                table_for_venn = TRUE)
   })

shared_nested_collTissue_ct %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-tissueCollapsed_CellTypeOnly.tsv.gz", sep = "_"))
  )

shared_nested_collTissue_sg <- purrr::map2_df(cd34_vs_wholemnc_collTissue,
      names(cd34_vs_wholemnc_collTissue),
      ~ {
        print(paste0("Processing ", .y, "..."))
        regenerated <- regen2(.x)
        if (is.null(regenerated)) {
          return(NULL)
        }
        other_lineages <- purity_filtered[[.y]] %>%
          mutate(SubjectID = .y) %>%
          filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
          left_join(blood_lineages_default(),  by = "CellMarker") %>%
          filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
        if (nrow(other_lineages) == 0) {
          return(NULL)
        }
        is_sharing(regenerated, other_lineages,
                   group_keys = list(g1 = c("id"),
                                     g2 = c("SubjectID", "SuperGroup")),
                   table_for_venn = TRUE)
      })

shared_nested_collTissue_sg %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-tissueCollapsed_SuperGroup.tsv.gz", sep = "_"))
  )

cd34_vs_wholemnc_collTissue_ge12 <- purrr::map2(purity_filtered,
    names(purity_filtered),
    ~ {
      print(paste0("Processing ", .y, "..."))
      invivo <- .x %>%
        mutate(SubjectID = .y) %>%
        filter(as.numeric(TimepointMonths) >= 12, Tissue %in% c("PB", "BM"))
      whole_mnc_only <- invivo %>%
        filter(CellMarker %in% c("Whole", "MNC"))
      cd34_only <- invivo %>%
        filter(CellMarker == "CD34")
      if (nrow(cd34_only) == 0 || nrow(whole_mnc_only) == 0) {
        return(NULL)
      }
      is_sharing(cd34_only, whole_mnc_only,
                 group_keys = list(g1 = c("SubjectID", "CellMarker"),
                                   g2 = c("SubjectID")),
                 keep_genomic_coord = TRUE) %>%
        mutate(g2 = paste0(g2, "_Whole"))
    })

shared_nested_collTissue_ge12 <- purrr::map2_df(cd34_vs_wholemnc_collTissue_ge12,
  names(cd34_vs_wholemnc_collTissue_ge12),
  ~ {
    print(paste0("Processing ", .y, "..."))
    regenerated <- regen2(.x)
    if (is.null(regenerated)) {
      return(NULL)
    }
    other_lineages <- purity_filtered[[.y]] %>%
      mutate(SubjectID = .y) %>%
      filter(as.numeric(TimepointMonths) >= 12, Tissue %in% c("PB", "BM")) %>%
      left_join(blood_lineages_default(),  by = "CellMarker") %>%
      filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
    if (nrow(other_lineages) == 0) {
      return(NULL)
    }
    is_sharing(regenerated, other_lineages,
               group_keys = list(g1 = c("id"),
                                 g2 = c("SubjectID", "CellMarker")),
               table_for_venn = TRUE)
  })

shared_nested_collTissue_ct_ge12 <- purrr::map2_df(cd34_vs_wholemnc_collTissue_ge12,
  names(cd34_vs_wholemnc_collTissue_ge12),
  ~ {
    print(paste0("Processing ", .y, "..."))
    regenerated <- regen2(.x)
    if (is.null(regenerated)) {
      return(NULL)
    }
    other_lineages <- purity_filtered[[.y]] %>%
      mutate(SubjectID = .y) %>%
      filter(as.numeric(TimepointMonths) >= 12, Tissue %in% c("PB", "BM")) %>%
      left_join(blood_lineages_default(),  by = "CellMarker") %>%
      filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
    if (nrow(other_lineages) == 0) {
      return(NULL)
    }
    is_sharing(regenerated, other_lineages,
               group_keys = list(g1 = c("id"),
                                 g2 = c("SubjectID", "CellType")),
               table_for_venn = TRUE)
  })

shared_nested_collTissue_sg_ge12 <- purrr::map2_df(cd34_vs_wholemnc_collTissue_ge12,
    names(cd34_vs_wholemnc_collTissue_ge12),
    ~ {
      print(paste0("Processing ", .y, "..."))
      regenerated <- regen2(.x)
      if (is.null(regenerated)) {
        return(NULL)
      }
      other_lineages <- purity_filtered[[.y]] %>%
        mutate(SubjectID = .y) %>%
        filter(as.numeric(TimepointMonths) >= 12, Tissue %in% c("PB", "BM")) %>%
        left_join(blood_lineages_default(),  by = "CellMarker") %>%
        filter(CellType %in% c("Myeloid", "Erythroid", "B", "T"))
      if (nrow(other_lineages) == 0) {
        return(NULL)
      }
      is_sharing(regenerated, other_lineages,
                 group_keys = list(g1 = c("id"),
                                   g2 = c("SubjectID", "SuperGroup")),
                 table_for_venn = TRUE)
    })

shared_nested_collTissue_ge12 %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-tissueCollapsed_withCellMarker_FROM-TP-12.tsv.gz", sep = "_"))
  )

shared_nested_collTissue_ct_ge12 %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-tissueCollapsed_CellTypeOnly_FROM-TP-12.tsv.gz", sep = "_"))
  )

shared_nested_collTissue_sg_ge12 %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-tissueCollapsed_SuperGroup_FROM-TP-12.tsv.gz", sep = "_"))
  )

## Sharing Plasma vs. Myeloid, Erythro, B, T
### ** Only in vivo tp (no 0); all tp collapsed together;
### ** only PB & BM (separated)

plasma_vs_lineages_sepTissue <- purrr::map2_df(purity_filtered,
                                             names(purity_filtered),
  ~ {
    print(paste0("Processing ", .y, "..."))
    invivo <- .x %>%
      mutate(SubjectID = .y) %>%
      filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
      left_join(blood_lineages_default(),  by = "CellMarker")
     plasma_only <- invivo %>%
       filter(CellMarker == "Plasma")
     other_lineages <- invivo %>%
       filter(CellType %in% c("Myeloid", "Erythroid", "B", "T", "CD34"))
     if (nrow(plasma_only) == 0 || nrow(other_lineages) == 0) {
       return(NULL)
     }
     is_sharing(plasma_only, other_lineages,
                group_key = c("SubjectID", "CellMarker", "Tissue"))
     })

plasma_vs_lineages_sepTissue %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "Plasma-MyErBT-inVivo-tpCollapsed-PB-BM_withCellMarker.tsv.gz", sep = "_"))
                   )

cd34_vs_all_lineages_sepTissue_ct <- purrr::map2_df(purity_filtered, 
     names(purity_filtered),
     ~ {
       print(paste0("Processing ", .y, "..."))
       invivo <- .x %>% 
         mutate(SubjectID = .y) %>%
         filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
         left_join(blood_lineages_default(),  by = "CellMarker")
       cd34_only <- invivo %>%
         filter(CellMarker == "CD34")
       sharing_df <- NULL
       for (tissue in c("PB", "BM")) {
         myelo <- invivo %>%
          filter(CellType == "Myeloid", Tissue == tissue)
         ery <- invivo %>%
           filter(CellType == "Erythroid", Tissue == tissue)
         b <- invivo %>%
           filter(CellType == "B", Tissue == tissue)
         t <- invivo %>%
           filter(CellType == "T", Tissue == tissue)
         lineages_list <- list(cd34_only, myelo, ery, b, t)
         groups_list <- list()
         for (df in lineages_list) {
           if (nrow(df) > 0) {
             groups_list <- append(groups_list, list(df))
           }
         }
         if (length(groups_list) == 0) {
         return(NULL)
         }
         sh <- is_sharing(!!!groups_list,
                    group_key = c("SubjectID", "CellType", "Tissue"), 
                    table_for_venn = TRUE)
         if (is.null(sharing_df)) {
           sharing_df <- sh
         } else {
           sharing_df <- bind_rows(sharing_df, sh)
         }
       }
       sharing_df
     })

prettify <- function(venn_obj) {
  # Extract rownames
  row_names <- rownames(venn_obj$ellipses)
  split <- stringr::str_split(row_names, "_")
  patient <- split[[1]][1]
  tissue <- if (split[[1]][2] != "CD34") {
    stringr::str_replace(split[[1]][3], "\\(g[1-9]\\)", "") 
  } else {
    stringr::str_replace(split[[2]][3], "\\(g[1-9]\\)", "") 
  }
  markers <- purrr::map_chr(split, ~ .x[2])
  rownames(venn_obj$ellipses) <- markers
  venn_obj$Patient <- patient
  venn_obj$Tissue <- tissue
  venn_obj
}

venns <- sharing_venn(cd34_vs_all_lineages_sepTissue_ct, euler = FALSE)
eulers <- sharing_venn(cd34_vs_all_lineages_sepTissue_ct, euler = TRUE)

venns <- purrr::map(venns, prettify)
eulers <- purrr::map(eulers, prettify)
colorscale <- tibble::tribble(
  ~ color, ~ marker,
  "green4",	"CD34",
  "orange",	"Myeloid",
  "red3",	"Erythroid",
  "dodgerblue3",	"B",
  "deepskyblue",	"T"
)

for (venn in venns) {
  colors <- tibble::as_tibble(list(marker = rownames(venn$ellipses))) %>%
    left_join(colorscale) %>%
    pull(color)
  pdf(file = fs::path(sharing_fold, 
                      paste0(venn$Patient,
                             "_", venn$Tissue, ".",
                             paste(date, project_name, 
                                   "sharing_venn.pdf", 
                                    sep = "_"))), 
      width = 8, height = 7)
  Sys.sleep(0)
  print(plot(venn, 
             quantities = TRUE,
             legend = TRUE, 
             labels = FALSE, 
             fills = colors,
             main = paste(venn$Patient, venn$Tissue)))
  Sys.sleep(0)
  dev.off()
}

for (eul in eulers) {
  colors <- tibble::as_tibble(list(marker = rownames(eul$ellipses))) %>%
    left_join(colorscale) %>%
    pull(color)
  pdf(file = fs::path(sharing_fold, 
                      paste0(eul$Patient,
                             "_", eul$Tissue, ".",
                             paste(date, project_name, 
                                   "sharing_euler.pdf", 
                                    sep = "_"))), 
      width = 8, height = 7)
  Sys.sleep(0)
  print(plot(eul, quantities = TRUE,
     legend = TRUE, 
     labels = FALSE, 
     fills = colors,
     main = paste(eul$Patient, eul$Tissue)))
  Sys.sleep(0)
  dev.off()
}

## HSC POPULATION ESTIMATE
stable_tps <- c(9, 12)

HSC_estimate <- HSC_population_size_estimate(
  x = agg, 
  metadata = agg_meta,
  aggregation_key = agg_key,
  timepoint_column = "TimepointMonths", 
  stable_timepoints = stable_tps
)

HSC_estimate %>%
  readr::write_tsv(file = fs::path(project_folder, 
                                   paste(date, project_name, "HSC_population_size_estimate.tsv", 
                                   sep = "_"), na = ""))

hsc_plot <- HSC_population_plot(HSC_estimate, project_name = "MLD")
hsc_plot <- hsc_plot +
  theme_bw()
ggsave(plot = hsc_plot, path = plots_fold,
       paste(date, project_name, 
             "HSC_population_size_estimate_plot.pdf", sep = "_"),
       width = 8, height = 5)
ggsave(plot = hsc_plot, path = plots_fold,
       paste(date, project_name, 
             "HSC_population_size_estimate_plot.png", sep = "_"),
       width = 8, height = 5)

# SOURCE OF ISS
## For groups CD34 - Myelo - Ery - T - B, at single marker level, intramarker
selected_markers <- blood_lineages_default() %>%
  filter(CellType %in% c("CD34", "Myeloid", "Erythroid", "B", "T")) %>%
  distinct(CellMarker) %>%
  pull(CellMarker)

source_of_iss <- list()

for (marker in selected_markers) {
  for (patient_df in agg_by_patient) {
    print(paste0("Processing marker ", marker, " for patient ", 
                 patient_df$SubjectID[1]))
    patient_df <- patient_df %>%
      filter(CellMarker == marker)
    if (nrow(patient_df) > 0) {
      iss_s <- iss_source(reference = patient_df, selection = patient_df, 
                          ref_group_key = agg_key, selection_group_key = agg_key, 
                          timepoint_column = "TimepointMonths")
      patient_source <- source_of_iss[[patient_df$SubjectID[1]]]
      if (is.null(patient_source)) {
        source_of_iss[[patient_df$SubjectID[1]]] <- iss_s[[1]]
      } else {
        source_of_iss[[patient_df$SubjectID[1]]] <- bind_rows(
          source_of_iss[[patient_df$SubjectID[1]]], iss_s[[1]])
      }
    }
  }
}

source_all_table <- purrr::reduce(source_of_iss, bind_rows)
source_all_table <- source_all_table %>%
  filter(g1_Tissue == g2_Tissue)

source_all_table %>%
  select(-is_coord) %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         "source_of_iss.tsv.gz",
                                         sep = "_")))

selected_supergroups <- blood_lineages_default() %>%
  filter(CellType %in% c("CD34", "Myeloid", "Erythroid", "B", "T")) %>%
  distinct(SuperGroup) %>%
  pull(SuperGroup)

agg_sg_per_patient <- agg_sg %>%
  filter(as.numeric(TimepointMonths) > 0, SubjectID %in% patients) %>%
  group_by(SubjectID) %>%
  group_split()

source_of_iss_sg <- list()
for (marker in selected_supergroups) {
  for (patient_df in agg_sg_per_patient) {
    print(paste0("Processing marker ", marker, " for patient ", 
                 patient_df$SubjectID[1]))
    patient_df <- patient_df %>%
      filter(SuperGroup == marker)
    if (nrow(patient_df) > 0) {
      iss_s <- iss_source(reference = patient_df, selection = patient_df, 
                          ref_group_key = agg_key_sg, selection_group_key = agg_key_sg, 
                          timepoint_column = "TimepointMonths")
      patient_source <- source_of_iss_sg[[patient_df$SubjectID[1]]]
      if (is.null(patient_source)) {
        source_of_iss_sg[[patient_df$SubjectID[1]]] <- iss_s[[1]]
      } else {
        source_of_iss_sg[[patient_df$SubjectID[1]]] <- bind_rows(
          source_of_iss_sg[[patient_df$SubjectID[1]]], iss_s[[1]])
      }
    }
  }
}

source_all_table_sg <- purrr::reduce(source_of_iss_sg, bind_rows)
source_all_table_sg <- source_all_table_sg %>%
  filter(g1_Tissue == g2_Tissue)

source_all_table_sg %>%
  select(-is_coord) %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         "source_of_iss_SuperGroup.tsv.gz",
                                         sep = "_")))

source_of_iss_sg_purityfilt <- list()
purity_mod <- purrr::map2(purity_filtered, names(purity_filtered), ~ {
  .x %>%
    mutate(SubjectID = .y) %>%
    left_join(blood_lineages_default(), by = "CellMarker")
})

for (marker in selected_supergroups) {
  for (patient_df in purity_mod) {
    print(paste0("Processing marker ", marker, " for patient ", 
                 patient_df$SubjectID[1]))
    patient_df <- patient_df %>%
      filter(SuperGroup == marker, as.numeric(TimepointMonths) > 0)
    if (nrow(patient_df) > 0) {
      iss_s <- iss_source(reference = patient_df, selection = patient_df, 
                          ref_group_key = agg_key_sg, 
                          selection_group_key = agg_key_sg, 
                          timepoint_column = "TimepointMonths")
      patient_source <- source_of_iss_sg_purityfilt[[patient_df$SubjectID[1]]]
      if (is.null(patient_source)) {
        source_of_iss_sg_purityfilt[[patient_df$SubjectID[1]]] <- iss_s[[1]]
      } else {
        source_of_iss_sg_purityfilt[[patient_df$SubjectID[1]]] <- bind_rows(
          source_of_iss_sg_purityfilt[[patient_df$SubjectID[1]]], iss_s[[1]])
      }
    }
  }
}

source_all_table_sg_purityfilt <- purrr::reduce(source_of_iss_sg_purityfilt, bind_rows)
source_all_table_sg_purityfilt <- source_all_table_sg_purityfilt %>%
  filter(g1_Tissue == g2_Tissue)

source_all_table_sg_purityfilt %>%
  select(-is_coord) %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         "source_of_iss_SuperGroup_purityfilt.tsv.gz",
                                         sep = "_")))

#------------------------------------------------------------------------------
# WHOLE + MNC Workflow
## ** In these steps MNC cell marker is replaced both in matrix and af by 
## ** Whole and sensitive steps are re-run.
## ** NOTE: all sharing steps are NOT affected by this since a re-aggregation
## ** step is performed prior calculations. Also not affected: CIS

# SUBSTITUTION
no_mnc_af <- outliers_removed %>%
  mutate(CellMarker = if_else(CellMarker == "MNC", "Whole", CellMarker))
no_mnc_agg <- aggregate_values_by_key(final_matrix,
                                      association_file = no_mnc_af, 
                                      value_cols = c("seqCount", "fragmentEstimate"),
                                      key = agg_key)
no_mnc_agg_meta <- aggregate_metadata(no_mnc_af, grouping_keys = agg_key)

# STATS
no_mnc_stats <- sample_statistics(no_mnc_agg, no_mnc_agg_meta,
                                  sample_key = agg_key,
                                  value_columns = c("seqCount_sum",
                                                    "fragmentEstimate_sum"))

# Filtering: we apply a filter for IS 
# with sc < 20 (only on markers that are NOT 'Plasma')
to_remove <- no_mnc_stats$metadata %>%
  filter(nIS < 20, CellMarker != "Plasma")

no_mnc_agg_filtered <- no_mnc_agg %>%
  anti_join(to_remove, by = agg_key)

no_mnc_stats <- sample_statistics(no_mnc_agg_filtered, no_mnc_agg_meta,
                                  sample_key = agg_key,
                                  value_columns = c("seqCount_sum",
                                                    "fragmentEstimate_sum"))

no_mnc_stats$metadata %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste0(stats_file_prefix, 
                                          "_",
                                          concat_agg_key,
                                          "_NO-MNC.tsv.gz")),
                   na = "")

# ABUNDANCE
no_mnc_abund <- compute_abundance(no_mnc_agg_filtered, key = agg_key)
readr::write_tsv(no_mnc_abund, 
                 file = fs::path(project_folder, 
                                 paste0(paste(date, project_name,
                                              "abundance", 
                                              concat_agg_key,
                                              sep = "_"), "_NO-MNC.tsv.gz")
                 ),
                 na = "")

# ALLUVIAL - plotting Whole only for each patient
alluv_nomnc <- integration_alluvial_plot(x = no_mnc_abund %>%
                                     filter(SubjectID %in% patients,
                                            CellMarker == "Whole"), 
                                   top_abundant_tbl = TRUE, 
                                   alluvia_plot_y_threshold = 1,
                                   plot_x = "TimepointMonths")
purrr::walk2(alluv_nomnc, names(alluv_nomnc), ~ {
  prettified_id <- stringr::str_replace_all(.y, "_", ", ")
  p <- .x$plot +
    theme_bw() +
    theme(legend.position = "none") +
    ggplot2::labs(title = paste0("MLD - IS Abundance over time - ",
                                 prettified_id),
                  subtitle = paste0("Colored flows indicate IS having ",
                                    "a relative abundance >= 1%\n",
                                    "in at least 1 time point"),
                  x = "Months after GT",
                  y = "Abundace (%)")
  file_prefix <- paste0(.y, ".",
                        paste(date, project_name, sep = "_"))
  ggplot2::ggsave(plot = p, path = alluv_fold, 
                  filename = paste0(file_prefix, "_alluvial-plot_NO-MNC.pdf"), 
                  width = 8, height = 5)
  ggplot2::ggsave(plot = p, path = alluv_fold, 
                  filename = paste0(file_prefix, "_alluvial-plot_NO-MNC.png"), 
                  width = 8, height = 5)
  pdf(file = fs::path(alluv_fold, 
                      paste0(file_prefix, "_top-abund-tbl_NO-MNC.pdf")), 
      width = 13, height = 7)
  gridExtra::grid.arrange(.x$tables)
  dev.off()
})

# DIVERSITY PLOT - WHOLE/MNC
h_index_data <- no_mnc_stats$metadata %>%
  mutate(TimepointMonths = as.numeric(TimepointMonths)) %>%
  filter(CellMarker %in% c("Whole"),
         TimepointMonths > 0,
         SubjectID %in% patients) %>%
  select(SubjectID, Tissue, TimepointMonths, 
         seqCount_sum_shannon, fragmentEstimate_sum_shannon)

h_index_plot <- ggplot(data = h_index_data, aes(x = TimepointMonths,
                                                y = seqCount_sum_shannon,
                                                group = SubjectID, 
                                                color = SubjectID)) +
  geom_point(size = 2.5) +
  geom_line(size = 1.5, alpha = .7) +
  facet_wrap(~Tissue) +
  #geom_smooth(aes(group=Tissue)) +
  theme_bw() +
  labs(y = "H index", x = "Time point (months after GT)",
       title = "BTHAL Clonal population diversity over time - Whole/MNC") +
  scale_x_continuous(breaks = seq(0, max(h_index_data$TimepointMonths, 
                                         na.rm = T), 6) ) +
  scale_y_continuous(limits = c(0, 10)) +
  theme(axis.text.x = element_text(size=16), 
        axis.text.y = element_text(size=16), 
        axis.title = element_text(size=16), 
        plot.title = element_text(size=20, face = "bold"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.key.width = unit(0.2, "cm")) +
  guides(color = guide_legend(nrow = 2))

ggsave(plot = h_index_plot, filename = paste0(paste(date, project_name, 
                                                    "hindex-plot-WHOLE-MNC-withSmooth",
                                                    sep = "_"), ".pdf"),
       width = 13, height = 8, path = plots_fold)
ggsave(plot = h_index_plot, filename = paste0(paste(date, project_name, 
                                                    "hindex-plot-WHOLE-MNC-withSmooth",
                                                    sep = "_"), ".png"),
       width = 13, height = 8, path = plots_fold)
ggsave(plot = h_index_plot, filename = paste0(paste(date, project_name, 
                                                    "hindex-plot-WHOLE-MNC-NOSmooth",
                                                    sep = "_"), ".pdf"),
       width = 13, height = 8, path = plots_fold)
ggsave(plot = h_index_plot, filename = paste0(paste(date, project_name, 
                                                    "hindex-plot-WHOLE-MNC-NOSmooth",
                                                    sep = "_"), ".png"),
       width = 13, height = 8, path = plots_fold)

# CUMULATED MATRICES
cumulated_no_mnc <- cumulative_is(no_mnc_agg_filtered,
                                  key = agg_key,
                                  timepoint_col = "TimepointMonths",
                                  keep_og_is = FALSE,
                                  expand = TRUE)
cumulated_per_patient_no_mnc <- cumulated_no_mnc %>%
  group_by(SubjectID) %>%
  group_split()
for (df in cumulated_per_patient_no_mnc) {
  subj <- df$SubjectID[1]
  print(paste("Saving:", subj, "..."))
  file_name <- paste0(subj, ".",
                      paste(date, project_name, 
                            paste0(agg_key, collapse = "_"),
                            sep = "_"),
                      "_cumulated_is_NO-MNC.tsv.gz")
  df %>%
    readr::write_tsv(file = fs::path(matrix_cumulative_folder, file_name), 
                     na = "")
}

rm(cumulated_per_patient_no_mnc, cumulated_no_mnc, subj, file_name, df)
gc()



#-----------------------------------------------------------------------------
# Analyses with timeframes instead of timepoints
#-----------------------------------------------------------------------------
af_timeframes <- outliers_removed %>%
  mutate(TimeFrame = if_else(
    condition = as.numeric(TimepointMonths) %in% c(1:6),
    true = "early",
    false = if_else(
      condition = as.numeric(TimepointMonths) %in% c(12:24),
      true = "mid",
      false = if_else(
        condition = as.numeric(TimepointMonths) >= 36,
        true = "steady",
        false = NA_character_
      )
    )
  ))

timeframes_key <- c("SubjectID", "CellMarker", "Tissue", "TimeFrame")

agg_timeframes <- aggregate_values_by_key(
  final_matrix, 
  association_file = af_timeframes, 
  value_cols = "seqCount",
  key = timeframes_key
)

concat_tf_key <- paste0(timeframes_key, collapse = "_")
matrix_file_prefix_tf <- paste(date, project_name, concat_tf_key,
                            sep = "_")
per_patient_sparse(agg_timeframes,
                   quantific = c(seqCount = "seqCount_sum",
                                 fragmentEstimate = "fragmentEstimate_sum"),
                   prefix = matrix_file_prefix_tf,
                   dir_name = matrix_agg_folder, 
                   row_totals = TRUE)

per_patient_sparse(agg_timeframes %>% filter(seqCount_sum >= 3),
                   quantific = c(seqCount = "seqCount_sum",
                                 fragmentEstimate = "fragmentEstimate_sum"),
                   prefix = paste0(matrix_file_prefix_tf, "_SC-filtered"),
                   dir_name = matrix_agg_folder, 
                   row_totals = TRUE)


cd34_output_tf <- purrr::map2_df(purity_filtered, names(purity_filtered),
                              ~ {
                                temp <- .x %>% 
                                  mutate(SubjectID = .y) %>%
                                  left_join(af_timeframes %>%
                                              select("SubjectID", "CellMarker", 
                                                     "Tissue", "TimepointMonths",
                                                     "TimeFrame")) %>%
                                  dplyr::left_join(blood_lineages_default(), 
                                                   by = "CellMarker")
                                cd34_only <- temp %>%
                                  filter(CellMarker == "CD34")
                                e_l_m <- temp %>%
                                  filter(HematoLineage %in% c("Erythroid",
                                                              "Myeloid",
                                                              "Lymphoid"))
                                if (nrow(cd34_only) == 0 || nrow(e_l_m) == 0) {
                                  return(NULL)
                                }
                                shared <- is_sharing(cd34_only, e_l_m,
                                                     group_keys = list(g1 = c("SubjectID", "CellMarker", "Tissue"),
                                                                       g2 = c("SubjectID","CellMarker", "Tissue", "TimeFrame")),
                                                     keep_genomic_coord = TRUE)
                                shared
                              })

cd34_output_tf %>%
  select(-is_coord) %>%
  readr::write_tsv(
    file = fs::path(sharing_fold,
                    paste(date, project_name,
                          "CD34-Output_TimeFrame.tsv.gz", sep = "_")),
    na = ""
  )

to_plot_cd34_output_tf <- cd34_output_tf %>%
  tidyr::separate(col = "g2", into = c("SubjectID", "CellMarker", 
                                       "Tissue", "TimeFrame"),
                  convert = TRUE) %>%
  # filter(!SubjectID %in% c("MLDC02", "MLDCO2", "MLDCUP01", "MLDCUP02", 
  #                          "MLDCUP03", "MLDCUP04", "MLDCUP05", "MLDHE01",  
  #                          "MLDHE02","MLDHE03")) %>%
  left_join(blood_lineages_default(), by = "CellMarker")


cd34_output_plot_tf <- ggplot(to_plot_cd34_output_tf %>%
                                filter(!is.na(TimeFrame)), 
                           aes(x = TimeFrame, 
                               y = on_g1, fill = CellType, group = CellType,
                               color = CellType, shape = CellMarker)) +
  # stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2)
  stat_smooth(size = 2, level = 0.4, se = T) +
  labs(title = "MLD - CD34 Output all patients",
       y = "% IS Shared with CD34", x = "Months after GT") +
  theme_bw()
ggsave(plot = cd34_output_plot_tf,
       path = sharing_fold,
       filename = paste(date, project_name, "CD34-Output-plot_TimeFrame.pdf", sep = "_"),
       width = 9, height = 5)
ggsave(plot = cd34_output_plot_tf,
       path = sharing_fold,
       filename = paste(date, project_name, "CD34-Output-plot_TimeFrame.png", sep = "_"),
       width = 9, height = 5)
