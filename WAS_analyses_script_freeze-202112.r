###--- WAS - analysis script ---###
## Latest update: 16/02/2022
library(ISAnalytics)
library(dplyr)
library(ggplot2)

date <- "20220509"
date_month <- "202112"
project_name <- "WAS"

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
source(".../per_patient_matrix_save.R") # saves matrices per patient

# PRELIMINARY STEPS ----
## AF loading
af_path <- "/home/giulia/Asso.all.projects.final.formatting.csv" # latest version

## AF IMPORT
af <- import_association_file(af_path, root = "",
                              separator = ',',
                              filter_for = list(ProjectID = "WAS"),
                              report_path = report_folder,
                              import_iss = TRUE,
                              file_prefixes = c(default_iss_file_prefixes(),
                                                "stats\\.was."),
                              dates_format = 'mdy')
patients <- unique(af$SubjectID[!af$SubjectID %in% 
                                c("UTR", "FB", "UT", "CEM37",
                                  "CEM41.W27", "HD")])
agg_key <- c("SubjectID", "CellMarker", "Tissue", "TimepointMonths")
agg_key_days <- c("SubjectID", "CellMarker", "Tissue", "TimePoint")

## MATRICES IMPORT
matrices <- import_parallel_Vispa2Matrices(af, 
                                           quantification_type = c("seqCount", 
                                                                   "fragmentEstimate"), 
                                           workers = 5, mode = "AUTO", 
                                           report_path = report_folder)

step_1_summary <- compute_intermediate_stats(matrices, af, patients,
                                             agg_key_days, step_suffix = "input")
step_1_summary %>%
  readr::write_tsv(file = fs::path(project_folder, paste(
    date, project_name, "summary_post_import.tsv", sep = "_"
  )))

## RECALIBRATION
recalibr <- compute_near_integrations(matrices, file_path = report_folder)

step_2_summary <- compute_intermediate_stats(recalibr, af, patients,
                                             agg_key_days, step_suffix = "post_recalibr")
step_2_summary %>%
  readr::write_tsv(file = fs::path(project_folder, paste(
    date, project_name, "summary_post_recalibr.tsv", sep = "_"
  )))

## COLLISIONS
coll <- remove_collisions(recalibr, 
                          association_file = af, 
                          max_workers = 15,
                          report_path = report_folder)

step_3_summary <- compute_intermediate_stats(coll, af, patients,
                                             agg_key_days, 
                                             step_suffix = "post_coll")
step_3_summary %>%
  readr::write_tsv(file = fs::path(project_folder, paste(
    date, project_name, "summary_post_coll.tsv", sep = "_"
  )))

## OUTLIER REMOVAL
outliers_removed <- outlier_filter(af, report_path = report_folder)

final_matrix <- coll %>% 
  filter(CompleteAmplificationID %in% outliers_removed$CompleteAmplificationID)

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

### LAM samples always lack fragment estimate: replacing missing values with
### sequence count value (mixed quantification)
final_matrix <- final_matrix %>%
  mutate(fragmentEstimate = if_else(is.na(fragmentEstimate),
                                    seqCount, fragmentEstimate))

# AGGREGATION
concat_agg_key <- paste0(agg_key, collapse = "_")
agg <- aggregate_values_by_key(final_matrix, 
                               association_file = outliers_removed, 
                               value_cols = c("seqCount", "fragmentEstimate"),
                               key = agg_key)
agg_meta <- aggregate_metadata(outliers_removed,
                               grouping_keys = agg_key)


agg_key_supergroup <- c("SubjectID", "SuperGroup", "Tissue", "TimepointMonths")
concat_agg_key_supergroup <- paste0(agg_key_supergroup, collapse = "_")
af_blood_l <- outliers_removed %>%
  select(-Keywords) %>%
  left_join(blood_lineages_default(), by = "CellMarker")

agg_supergroup <- aggregate_values_by_key(final_matrix, 
                                          association_file = af_blood_l, 
                                          value_cols = c("seqCount", "fragmentEstimate"),
                                          key = agg_key_supergroup)

agg_meta_super <- aggregate_metadata(af_blood_l,
                               grouping_keys = agg_key_supergroup)

agg_key_subj <- c("SubjectID", "TimepointMonths")
concat_agg_key_subj <- paste0(agg_key_subj, collapse = "_")

agg_subj_tp <- aggregate_values_by_key(final_matrix, 
                                          association_file = af_blood_l, 
                                          value_cols = c("seqCount", "fragmentEstimate"),
                                          key = agg_key_subj)

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

matrix_file_prefix_super <- paste(date, project_name, 
                                  paste0(agg_key_supergroup, collapse = "_"),
                            sep = "_")
per_patient_sparse(agg_supergroup,
                   quantific = c(seqCount = "seqCount_sum",
                                 fragmentEstimate = "fragmentEstimate_sum"),
                   prefix = matrix_file_prefix_super,
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


cumulated_sg <- cumulative_is(agg_supergroup,
                             key = agg_key_supergroup,
                             timepoint_col = "TimepointMonths",
                             keep_og_is = FALSE,
                             expand = TRUE)

cumulated_per_patient_sg <- cumulated_sg %>%
  group_by(SubjectID) %>%
  group_split()
for (df in cumulated_per_patient_sg) {
  subj <- df$SubjectID[1]
  print(paste("Saving:", subj, "..."))
  file_name <- paste0(subj, ".",
                      paste(date, project_name, 
                            concat_agg_key_supergroup,
                            sep = "_"),
                      "_cumulated_is.tsv.gz")
  df %>%
    readr::write_tsv(file = fs::path(matrix_cumulative_folder, file_name), 
                     na = "")
}

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

cumulate_count_sg <- cumulative_count_union(agg_supergroup, 
                                         key = agg_key_supergroup,
                                         timepoint_column = "TimepointMonths", 
                                         zero = "00")

cumulate_count_sg %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         concat_agg_key_supergroup,
                                         "cumulated_counts.tsv", 
                                         sep = "_")))

cumulate_count_sub <- cumulative_count_union(agg_subj_tp, 
                                            key = agg_key_subj,
                                            timepoint_column = "TimepointMonths", 
                                            zero = "00")
cumulate_count_sub %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         concat_agg_key_subj,
                                         "cumulated_counts.tsv", 
                                         sep = "_")))

cumulated_sg_from_12 <- cumulative_is(agg_supergroup %>%
                                        filter(as.numeric(TimepointMonths) >= 12),
                              key = agg_key_supergroup,
                              timepoint_col = "TimepointMonths",
                              keep_og_is = FALSE,
                              expand = TRUE)

cumulated_per_patient_sg_from_12 <- cumulated_sg_from_12 %>%
  group_by(SubjectID) %>%
  group_split()
for (df in cumulated_per_patient_sg_from_12) {
  subj <- df$SubjectID[1]
  print(paste("Saving:", subj, "..."))
  file_name <- paste0(subj, ".",
                      paste(date, project_name, 
                            concat_agg_key_supergroup,
                            sep = "_"),
                      "_cumulated_is_FROM-TP-12.tsv.gz")
  df %>%
    readr::write_tsv(file = fs::path(matrix_cumulative_folder, file_name), 
                     na = "")
}

cumulate_count_sg_from_12 <- cumulative_count_union(agg_supergroup %>%
                                                      filter(as.numeric(TimepointMonths) >= 12), 
                                            key = agg_key_supergroup,
                                            timepoint_column = "TimepointMonths", 
                                            zero = "00")

cumulate_count_sg_from_12 %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         concat_agg_key_supergroup,
                                         "cumulated_counts_FROM-TP-12.tsv", 
                                         sep = "_")))

cumulated_subj_from_12 <- cumulative_is(agg_subj_tp %>%
                                          filter(as.numeric(TimepointMonths) >= 12),
                                key = agg_key_subj,
                                timepoint_col = "TimepointMonths",
                                keep_og_is = FALSE,
                                expand = TRUE)

cumulated_per_patient_subj_from_12 <- cumulated_subj_from_12 %>%
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

cumulate_count_sub_from_12 <- cumulative_count_union(agg_subj_tp %>%
                                                       filter(as.numeric(TimepointMonths) >= 12), 
                                             key = agg_key_subj,
                                             timepoint_column = "TimepointMonths", 
                                             zero = "00")
cumulate_count_sub_from_12 %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         concat_agg_key_subj,
                                         "cumulated_counts_FROM-TP-12.tsv", 
                                         sep = "_")))

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
concat_pcr_methods <- outliers_removed %>%
  select(all_of(c(agg_key, "PCRMethod", "NGSTechnology"))) %>%
  group_by(across(all_of(agg_key))) %>%
  summarise(PCRMethod = paste0(unique(PCRMethod), collapse = "|"),
            NGSTechnology = paste0(unique(NGSTechnology), collapse = "|"), 
            .groups = "drop")

agg_stats <- sample_statistics(agg,
                               agg_meta,
                               sample_key = agg_key,
                               value_columns = c("seqCount_sum",
                                                 "fragmentEstimate_sum"))
agg_stats <- agg_stats$metadata %>%
  left_join(concat_pcr_methods, by = agg_key)
agg_stats %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste0(stats_file_prefix, 
                                          "_",
                                          concat_agg_key,
                                          ".tsv.gz")),
                   na = "")

agg_stats_super <- sample_statistics(agg_supergroup,
                               agg_meta_super,
                               sample_key = agg_key_supergroup,
                               value_columns = c("seqCount_sum",
                                                 "fragmentEstimate_sum"))

agg_stats_super$metadata %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste0(stats_file_prefix, 
                                          "_",
                                          paste0(agg_key_supergroup, collapse = "_"),
                                          ".tsv.gz")),
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
  plot_title <- paste0("WAS - ", .y, "\n")
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
    ggplot2::labs(title = paste0("WAS - IS Abundance over time - ",
                                 prettified_id),
                  subtitle = paste0("Colored flows indicate IS having ",
                                    "a relative abundance >= 1%\n",
                                    "in at least 1 time point"),
                  x = "Months after GT",
                  y = "Abundance (%)")
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

# cd34_output <- purrr::map2_df(purity_filtered, names(purity_filtered),
#                               ~ {
#                                 temp <- .x %>% mutate(SubjectID = .y) %>%
#                                   dplyr::left_join(blood_lineages_default(), 
#                                                    by = "CellMarker")
#                                 cd34_only <- temp %>%
#                                   filter(CellMarker == "CD34")
#                                 e_l_m <- temp %>%
#                                   filter(HematoLineage %in% c("Erythroid",
#                                                               "Myeloid",
#                                                               "Lymphoid"))
#                                 if (nrow(cd34_only) == 0 || nrow(e_l_m) == 0) {
#                                   return(NULL)
#                                 }
#                                 shared <- is_sharing(cd34_only, e_l_m,
#                                                      group_keys = list(g1 = c("SubjectID", "CellMarker", "Tissue"),
#                                                                        g2 = c("SubjectID","CellMarker", "Tissue", "TimepointMonths")),
#                                                      keep_genomic_coord = TRUE)
#                                 shared
#                               })
# 
# cd34_output %>%
#   select(-is_coord) %>%
#   readr::write_tsv(
#     file = fs::path(sharing_fold,
#                     paste(date, project_name,
#                           "CD34-Output.tsv.gz", sep = "_")),
#     na = ""
#   )
# 
# to_plot_cd34_output <- cd34_output %>%
#   tidyr::separate(col = "g2", into = c("SubjectID", "CellMarker", 
#                                        "Tissue", "TimepointMonths"),
#                   convert = TRUE) %>%
#   # filter(!SubjectID %in% c("MLDC02", "MLDCO2", "MLDCUP01", "MLDCUP02", 
#   #                          "MLDCUP03", "MLDCUP04", "MLDCUP05", "MLDHE01",  
#   #                          "MLDHE02","MLDHE03")) %>%
#   left_join(blood_lineages_default(), by = "CellMarker")
# 
# 
# cd34_output_plot <- ggplot(to_plot_cd34_output, 
#                            aes(x = TimepointMonths, 
#                                y = on_g1, fill = CellType, group = CellType,
#                                color = CellType, shape = CellMarker)) +
#   # stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2)
#   stat_smooth(size = 2, level = 0.4, se = T) +
#   labs(title = "WAS - CD34 Output all patients",
#        y = "% IS Shared with CD34", x = "Months after GT") +
#   theme_bw()
# ggsave(plot = cd34_output_plot,
#        path = sharing_fold,
#        filename = paste(date, project_name, "CD34-Output-plot.pdf", sep = "_"),
#        width = 9, height = 5)
# ggsave(plot = cd34_output_plot,
#        path = sharing_fold,
#        filename = paste(date, project_name, "CD34-Output-plot.png", sep = "_"),
#        width = 9, height = 5)

cd34_output_super <- purrr::map2_df(purity_filtered, names(purity_filtered),
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

cd34_output_super %>%
  select(-is_coord) %>%
  readr::write_tsv(
    file = fs::path(sharing_fold,
                    paste(date, project_name,
                          "CD34-Output_by-SuperGroup.tsv.gz", sep = "_")),
    na = ""
  )

to_plot_cd34_output_super <- cd34_output_super %>%
  tidyr::separate(col = "g2", into = c("SubjectID", "SuperGroup", 
                                       "Tissue", "TimepointMonths"),
                  convert = TRUE) %>%
  # filter(!SubjectID %in% c("MLDC02", "MLDCO2", "MLDCUP01", "MLDCUP02", 
  #                          "MLDCUP03", "MLDCUP04", "MLDCUP05", "MLDHE01",  
  #                          "MLDHE02","MLDHE03")) %>%
  left_join(blood_lineages_default(), by = c("SuperGroup" = "CellMarker"))


cd34_output_plot_super <- ggplot(to_plot_cd34_output_super, 
                                     aes(x = TimepointMonths, 
                                         y = on_g1, fill = CellType, group = CellType,
                                         color = CellType)) +
  # stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2)
  stat_smooth(size = 2, level = 0.4, se = T) +
  labs(title = "WAS - CD34 Output all patients",
       y = "% IS Shared with CD34", x = "Months after GT") +
  theme_bw()
ggsave(plot = cd34_output_plot_super,
       path = sharing_fold,
       filename = paste(date, project_name, "CD34-Output-plot_by-SuperGroup.pdf", sep = "_"),
       width = 9, height = 5)
ggsave(plot = cd34_output_plot_super,
       path = sharing_fold,
       filename = paste(date, project_name, "CD34-Output-plot_by-SuperGroup.png", sep = "_"),
       width = 9, height = 5)

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
#     file = fs::path(sharing_fold,
#                     paste(date, project_name,
#                           "CD34-WholeMNC-inVivo-tpColapsed-PB-BM.tsv.gz", sep = "_"))
#   )

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


cd34_vs_all_lineages_sepTissue_ct <- purrr::map2_df(purity_filtered, 
     names(purity_filtered),
     ~ {
       print(paste0("Processing ", .y, "..."))
       invivo <- .x %>% 
         mutate(SubjectID = .y) %>%
         filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
         left_join(blood_lineages_default(),  by = "CellMarker")
       cd34_only <- invivo %>%
         filter(CellMarker == "CD34", Tissue == "BM")
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
         if (length(groups_list) < 2) {
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

hsc_plot <- HSC_population_plot(HSC_estimate, project_name = "WAS")
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

agg_sg_per_patient <- agg_supergroup %>%
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
                          ref_group_key = agg_key_supergroup, 
                          selection_group_key = agg_key_supergroup, 
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
                          ref_group_key = agg_key_supergroup, 
                          selection_group_key = agg_key_supergroup, 
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
    ggplot2::labs(title = paste0("WAS - IS Abundance over time - ",
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
         SubjectID %in% patients,
         Tissue %in% c("PB", "BM")) %>%
  select(SubjectID, Tissue, TimepointMonths, 
         seqCount_sum_shannon, fragmentEstimate_sum_shannon)

h_index_plot <- ggplot(data = h_index_data, aes(x = TimepointMonths,
                                                y = seqCount_sum_shannon,
                                                group = SubjectID, 
                                                color = SubjectID)) +
  geom_point(size = 2.5) +
  geom_line(size = 1.5, alpha = .7) +
  facet_wrap(~Tissue) +
  # geom_smooth(aes(group=Tissue)) +
  theme_bw() +
  labs(y = "H index", x = "Time point (months after GT)",
       title = "WAS Clonal population diversity over time - Whole/MNC") +
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

#------------------------------------------------------------------------------
# WAS - only LAM workflow (requested)
#------------------------------------------------------------------------------
lam_matrix_fold <- fs::path(matrix_folder, "LAM_only")
lam_matrix_agg <- fs::path(lam_matrix_fold, "aggregated")
lam_cumulative_folder <- fs::path(lam_matrix_fold, "cumulative")
lam_purity_folder <- fs::path(lam_matrix_fold, "purity_filtered")
fs::dir_create(lam_matrix_fold)
fs::dir_create(lam_matrix_agg)
fs::dir_create(lam_cumulative_folder)
fs::dir_create(lam_purity_folder)

## AF IMPORT
af_lam_only <- import_association_file(af_path, root = "",
                                separator = ',',
                                filter_for = list(ProjectID = "WAS",
                                                  PCRMethod = "LAM-PCR"),
                                report_path = report_folder,
                                import_iss = TRUE,
                                file_prefixes = c(default_iss_file_prefixes(),
                                                  "stats\\.was."),
                                dates_format = 'mdy')
## MATRICES IMPORT
matrices_lam <- import_parallel_Vispa2Matrices(af_lam_only,
                                               quantification_type = "seqCount",
                                               workers = 5, mode = "AUTO", 
                                               report_path = report_folder)


## RECALIBRATION
recalibr <- compute_near_integrations(matrices_lam, 
                                      value_columns = "seqCount",
                                      file_path = report_folder)

## COLLISIONS
coll_lam <- remove_collisions(recalibr, 
                          association_file = af_lam_only, 
                          quant_cols = c(seqCount = "seqCount"),
                          max_workers = 10,
                          report_path = report_folder)

## OUTLIER REMOVAL
outliers_removed_lam <- outlier_filter(af_lam_only, 
                                       report_path = report_folder)

final_matrix_lam <- coll_lam %>% 
  filter(CompleteAmplificationID %in% outliers_removed_lam$CompleteAmplificationID)

# AGGREGATION
agg_lam <- aggregate_values_by_key(final_matrix_lam, 
                               association_file = outliers_removed_lam, 
                               value_cols = c("seqCount"),
                               key = agg_key)
#### Filter for sc < 3
agg_lam <- agg_lam %>%
  filter(seqCount_sum >= 3)
agg_meta_lam <- aggregate_metadata(outliers_removed_lam,
                               grouping_keys = agg_key)


af_blood_l_lam <- outliers_removed_lam %>%
  select(-Keywords) %>%
  left_join(blood_lineages_default(), by = "CellMarker")

agg_supergroup_lam <- aggregate_values_by_key(final_matrix_lam, 
                                          association_file = af_blood_l_lam, 
                                          value_cols = c("seqCount"),
                                          key = agg_key_supergroup)

agg_supergroup_lam_sc_filt <- agg_supergroup_lam %>%
  filter(seqCount_sum >= 3)
agg_meta_super_lam <- aggregate_metadata(af_blood_l_lam,
                                     grouping_keys = agg_key_supergroup)

agg_subj_tp_lam <- aggregate_values_by_key(final_matrix_lam, 
                                           association_file = outliers_removed_lam,
                                           value_cols = c("seqCount"),
                                           key = agg_key_subj)

# AGG MATRIX SAVE - PER PATIENT
### * NOTE: framentEstimate matrices will contain mixed quantifications!
matrix_file_prefix_lam_nosc <- paste(date, project_name, 
                                concat_agg_key, 
                                "LAM-only",
                                sep = "_")
per_patient_sparse(agg_lam,
                   quantific = c(seqCount = "seqCount_sum"),
                   prefix = matrix_file_prefix_lam_nosc,
                   dir_name = lam_matrix_agg, 
                   row_totals = TRUE)

matrix_file_prefix_lam <- paste(date, project_name, 
                                concat_agg_key, 
                                "LAM-only_SC-filtered",
                                sep = "_")
per_patient_sparse(agg_lam,
                   quantific = c(seqCount = "seqCount_sum"),
                   prefix = matrix_file_prefix_lam,
                   dir_name = lam_matrix_agg, 
                   row_totals = TRUE)

matrix_file_prefix_super_lam <- paste(date, project_name, 
                                      concat_agg_key_supergroup,
                                      "LAM-only",
                                      sep = "_")
per_patient_sparse(agg_supergroup_lam,
                   quantific = c(seqCount = "seqCount_sum"),
                   prefix = matrix_file_prefix_super_lam,
                   dir_name = lam_matrix_agg, 
                   row_totals = TRUE)

matrix_file_prefix_super_lam <- paste(date, project_name, 
                                  concat_agg_key_supergroup,
                                  "LAM-only_SC-filtered",
                                  sep = "_")
per_patient_sparse(agg_supergroup_lam_sc_filt,
                   quantific = c(seqCount = "seqCount_sum"),
                   prefix = matrix_file_prefix_super_lam,
                   dir_name = lam_matrix_agg, 
                   row_totals = TRUE)

# CUMULATED MATRICES
cumulated_lam <- cumulative_is(agg_lam,
                           key = agg_key,
                           timepoint_col = "TimepointMonths",
                           keep_og_is = FALSE,
                           expand = TRUE)
cumulated_per_patient_lam <- cumulated_lam %>%
  group_by(SubjectID) %>%
  group_split()
for (df in cumulated_per_patient_lam) {
  subj <- df$SubjectID[1]
  print(paste("Saving:", subj, "..."))
  file_name <- paste0(subj, ".",
                      paste(date, project_name, 
                            paste0(agg_key, collapse = "_"),
                            sep = "_"),
                      "_cumulated_is_LAM-only.tsv.gz")
  df %>%
    readr::write_tsv(file = fs::path(lam_cumulative_folder, file_name), 
                     na = "")
}

rm(cumulated_per_patient_lam, cumulated_lam, subj, file_name, df)
gc()

cumulated_lam_sg <- cumulative_is(agg_supergroup_lam,
                               key = agg_key_supergroup,
                               timepoint_col = "TimepointMonths",
                               keep_og_is = FALSE,
                               expand = TRUE)

cumulated_per_patient_lam_sg <- cumulated_lam_sg %>%
  group_by(SubjectID) %>%
  group_split()
for (df in cumulated_per_patient_lam_sg) {
  subj <- df$SubjectID[1]
  print(paste("Saving:", subj, "..."))
  file_name <- paste0(subj, ".",
                      paste(date, project_name, 
                            concat_agg_key_supergroup,
                            sep = "_"),
                      "_cumulated_is_LAM-only.tsv.gz")
  df %>%
    readr::write_tsv(file = fs::path(lam_cumulative_folder, file_name), 
                     na = "")
}

cumulated_lam_sub <- cumulative_is(agg_subj_tp_lam,
                                  key = agg_key_subj,
                                  timepoint_col = "TimepointMonths",
                                  keep_og_is = FALSE,
                                  expand = TRUE)

cumulated_per_patient_lam_sub <- cumulated_lam_sub %>%
  group_by(SubjectID) %>%
  group_split()
for (df in cumulated_per_patient_lam_sub) {
  subj <- df$SubjectID[1]
  print(paste("Saving:", subj, "..."))
  file_name <- paste0(subj, ".",
                      paste(date, project_name, 
                            concat_agg_key_subj,
                            sep = "_"),
                      "_cumulated_is_LAM-only.tsv.gz")
  df %>%
    readr::write_tsv(file = fs::path(lam_cumulative_folder, file_name), 
                     na = "")
}

cumulate_count_sg_lam <- cumulative_count_union(agg_supergroup_lam,
                                                timepoint_column = "TimepointMonths",
                                                key = agg_key_supergroup,
                                                zero = "00")
cumulate_count_sg_lam %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         concat_agg_key_supergroup,
                                         "cumulated_counts_LAM-only.tsv", 
                                         sep = "_")))

cumulate_count_sub_lam <- cumulative_count_union(agg_subj_tp_lam,
                                                timepoint_column = "TimepointMonths",
                                                key = agg_key_subj,
                                                zero = "00")
cumulate_count_sub_lam %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         concat_agg_key_subj,
                                         "cumulated_counts_LAM-only.tsv", 
                                         sep = "_")))

cumulated_lam_sg_from_12 <- cumulative_is(agg_supergroup_lam %>%
                                            filter(as.numeric(TimepointMonths) >= 12),
                                  key = agg_key_supergroup,
                                  timepoint_col = "TimepointMonths",
                                  keep_og_is = FALSE,
                                  expand = TRUE)

cumulated_per_patient_lam_sg_from_12 <- cumulated_lam_sg_from_12 %>%
  group_by(SubjectID) %>%
  group_split()
for (df in cumulated_per_patient_lam_sg_from_12) {
  subj <- df$SubjectID[1]
  print(paste("Saving:", subj, "..."))
  file_name <- paste0(subj, ".",
                      paste(date, project_name, 
                            concat_agg_key_supergroup,
                            sep = "_"),
                      "_cumulated_is_LAM-only_FROM-TP-12.tsv.gz")
  df %>%
    readr::write_tsv(file = fs::path(lam_cumulative_folder, file_name), 
                     na = "")
}

cumulate_count_sg_lam_from_12 <- cumulative_count_union(agg_supergroup_lam %>%
                                                  filter(as.numeric(TimepointMonths) >= 12),
                                                timepoint_column = "TimepointMonths",
                                                key = agg_key_supergroup,
                                                zero = "00")
cumulate_count_sg_lam_from_12 %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         concat_agg_key_supergroup,
                                         "cumulated_counts_LAM-only_FROM-TP-12.tsv", 
                                         sep = "_")))

cumulated_lam_sub_from_12 <- cumulative_is(agg_subj_tp_lam %>%
                                             filter(as.numeric(TimepointMonths) >= 12),
                                   key = agg_key_subj,
                                   timepoint_col = "TimepointMonths",
                                   keep_og_is = FALSE,
                                   expand = TRUE)

cumulated_per_patient_lam_sub_from_12  <- cumulated_lam_sub_from_12  %>%
  group_by(SubjectID) %>%
  group_split()
for (df in cumulated_per_patient_lam_sub_from_12) {
  subj <- df$SubjectID[1]
  print(paste("Saving:", subj, "..."))
  file_name <- paste0(subj, ".",
                      paste(date, project_name, 
                            concat_agg_key_subj,
                            sep = "_"),
                      "_cumulated_is_LAM-only_FROM-TP-12.tsv.gz")
  df %>%
    readr::write_tsv(file = fs::path(lam_cumulative_folder, file_name), 
                     na = "")
}

cumulate_count_sub_lam_from_12 <- cumulative_count_union(agg_subj_tp_lam  %>%
                                                           filter(as.numeric(TimepointMonths) >= 12),
                                                 timepoint_column = "TimepointMonths",
                                                 key = agg_key_subj,
                                                 zero = "00")
cumulate_count_sub_lam_from_12 %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         concat_agg_key_subj,
                                         "cumulated_counts_LAM-only_FROM-TP-12.tsv", 
                                         sep = "_")))

# DESCRIPTIVE STATS
## Single PCR
stats_file_prefix_lam <- paste(date, project_name, "descriptive-stats", sep = "_")
single_pcr_stats_lam <- sample_statistics(final_matrix_lam, 
                                      metadata = outliers_removed_lam,
                                      value_columns = c("seqCount"))
single_pcr_stats_lam$metadata %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste0(stats_file_prefix, "_single-pcr_LAM-only.tsv.gz")),
                   na = "")

## Aggregated
agg_stats_lam <- sample_statistics(agg_lam,
                               agg_meta_lam,
                               sample_key = agg_key,
                               value_columns = c("seqCount_sum"))
agg_stats_lam$metadata %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste0(stats_file_prefix, 
                                          "_",
                                          concat_agg_key,
                                          "_LAM-only.tsv.gz")),
                   na = "")

agg_stats_super_lam <- sample_statistics(agg_supergroup_lam,
                                     agg_meta_super_lam,
                                     sample_key = agg_key_supergroup,
                                     value_columns = c("seqCount_sum"))

agg_stats_super_lam$metadata %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste0(stats_file_prefix, 
                                          "_",
                                          paste0(agg_key_supergroup, collapse = "_"),
                                          "_LAM-only.tsv.gz")),
                   na = "")

agg_stats_super_lam_filt <- sample_statistics(agg_supergroup_lam_sc_filt,
                                              agg_meta_super_lam, 
                                              sample_key = agg_key_supergroup,
                                              value_columns = "seqCount_sum")

agg_stats_super_lam_filt$metadata %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste0(stats_file_prefix, 
                                          "_",
                                          paste0(agg_key_supergroup, collapse = "_"),
                                          "_LAM-only_SC-filtered.tsv.gz")),
                   na = "")

# ABUNDANCE
abund_lam <- compute_abundance(agg_lam, key = agg_key, columns = "seqCount_sum")
readr::write_tsv(abund_lam, 
                 file = fs::path(project_folder, 
                                 paste0(paste(date, project_name,
                                              "abundance", 
                                              concat_agg_key,
                                              sep = "_"), "_LAM-only.tsv.gz")
                 ),
                 na = "")

# CIS (in vivo)
cis_lam <- CIS_grubbs(agg_lam %>% 
                    filter(TimepointMonths != "00"), 
                  by = "SubjectID")
purrr::walk2(cis_lam[patients], names(cis_lam[patients]), ~ {
  filename <- paste0(.y, ".", paste(date, project_name, 
                                    "CISGrubbs_inVivo_LAM-only.tsv.gz",
                                    sep = "_"))
  readr::write_tsv(.x, file = fs::path(cis_folder, filename), na = "")
})

volcanos_lam <- purrr::map2(cis_lam[patients], names(cis_lam[patients]), ~ {
  if (!is.null(.x)) {
    plot_title <- paste0("WAS - ", .y, "\n")
    CIS_volcano_plot(.x, title_prefix = plot_title)
  }
})

purrr::walk2(volcanos_lam, names(volcanos_lam), ~ {
  if (!is.null(.x)) {
    prefix <- .y
    filename <- paste0(.y, ".", paste(date, project_name, 
                                      "volcano-plot_LAM-only",
                                      sep = "_"))
    ggplot2::ggsave(plot = .x, filename = paste0(filename, ".pdf"), 
                    path = cis_folder, 
                    width = 8, height = 6)
    ggplot2::ggsave(plot = .x, filename = paste0(filename, ".png"), 
                    path = cis_folder, 
                    width = 8, height = 6)
  }
})

# ALLUVIALS
alluv_lam <- integration_alluvial_plot(x = abund_lam %>%
                                     filter(SubjectID %in% patients), 
                                     top_abundant_tbl = TRUE, 
                                     plot_y = "seqCount_sum_PercAbundance",
                                     plot_x = "TimepointMonths")
purrr::walk2(alluv_lam, names(alluv_lam), ~ {
  prettified_id <- stringr::str_replace_all(.y, "_", ", ")
  p <- .x$plot +
    theme_bw() +
    theme(legend.position = "none") +
    ggplot2::labs(title = paste0("WAS - IS Abundance over time - ",
                                 prettified_id),
                  subtitle = paste0("Colored flows indicate IS having ",
                                    "a relative abundance >= 1%\n",
                                    "in at least 1 time point"),
                  x = "Months after GT",
                  y = "Abundance (%)")
  file_prefix <- paste0(.y, ".",
                        paste(date, project_name, sep = "_"))
  ggplot2::ggsave(plot = p, path = alluv_fold, 
                  filename = paste0(file_prefix, "_alluvial-plot_LAM-only.pdf"), 
                  width = 8, height = 5)
  ggplot2::ggsave(plot = p, path = alluv_fold, 
                  filename = paste0(file_prefix, "_alluvial-plot_LAM-only.png"), 
                  width = 8, height = 5)
  pdf(file = fs::path(alluv_fold, 
                      paste0(file_prefix, "_top-abund-tbl_LAM-only.pdf")), 
      width = 13, height = 7)
  gridExtra::grid.arrange(.x$tables)
  dev.off()
})

# TOP ABUNDANT IS/TOP TARGETED GENES
top20_abund_agg_lam <- top_integrations(abund_lam, key = agg_key, 
                                        columns = "seqCount_sum_RelAbundance")
top20_abund_agg_lam %>%
  readr::write_tsv(file = fs::path(project_folder, 
                                   paste0(paste(date, project_name,
                                                "top20-abundant-iss",
                                                concat_agg_key, sep = "_"),
                                          "_LAM-only.tsv.gz"
                                   )),
                   na = "")

n_is_by_gene_agg_lam <- agg_lam %>%
  group_by(across(all_of(c(annotation_IS_vars(), agg_key)))) %>%
  summarise(n_IS = n_distinct(chr, integration_locus, strand), 
            .groups = "drop")

top_20_targeted_agg_lam <- n_is_by_gene_agg_lam  %>% 
  group_by(across(all_of(agg_key))) %>%
  arrange(desc(n_IS)) %>% 
  slice_head(n = 20) %>%
  ungroup()

top_20_targeted_agg_lam %>%
  readr::write_tsv(file = fs::path(project_folder, 
                                   paste0(paste(date, project_name,
                                                "top20-targeted-genes",
                                                concat_agg_key, sep = "_"),
                                          "_LAM-only.tsv.gz"
                                   )),
                   na = "")

# SHARING
## CD34 output: for each patient, apply the purity filter and compute sharing
## between CD34+ and lineages Herithroid, Myeloid and Lymphoid.
## Purity filter is applied by timepoint only on groups of interest
agg_by_patient_lam <- agg_lam %>%
  filter(SubjectID %in% patients) %>%
  group_by(SubjectID) %>%
  group_split()

subj_names_lam <- purrr::map_chr(agg_by_patient_lam, ~ .x$SubjectID[1])

purity_filtered_lam <- purrr::map(agg_by_patient_lam, 
                              ~ purity_filter(.x, 
                                              selected_groups = groups_of_interest,
                                              aggregation_key = agg_key, 
                                              timepoint_column = "TimepointMonths"))
purity_filtered_lam <- purity_filtered_lam %>% purrr::set_names(subj_names_lam)
purrr::walk2(purity_filtered_lam, names(purity_filtered_lam), ~ {
  .x %>%
    readr::write_tsv(file = fs::path(
      lam_purity_folder,
      paste0(
        .y, ".",
        paste(date, project_name, "purity-filtered", sep = "_"),
        "_LAM-only.tsv.gz"
      )
    ), na = "")
})

cd34_output_lam <- purrr::map2_df(purity_filtered_lam, names(purity_filtered_lam),
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
                                                     group_keys = list(g1 = c("SubjectID", "CellMarker", "Tissue"),
                                                                       g2 = c("SubjectID","CellMarker", "Tissue", "TimepointMonths")),
                                                     keep_genomic_coord = TRUE)
                                shared
                              })

cd34_output_lam %>%
  select(-is_coord) %>%
  readr::write_tsv(
    file = fs::path(sharing_fold,
                    paste(date, project_name,
                          "CD34-Output_LAM-only.tsv.gz", sep = "_")),
    na = ""
  )

to_plot_cd34_output_lam <- cd34_output_lam %>%
  tidyr::separate(col = "g2", into = c("SubjectID", "CellMarker", 
                                       "Tissue", "TimepointMonths"),
                  convert = TRUE) %>%
  # filter(!SubjectID %in% c("MLDC02", "MLDCO2", "MLDCUP01", "MLDCUP02", 
  #                          "MLDCUP03", "MLDCUP04", "MLDCUP05", "MLDHE01",  
  #                          "MLDHE02","MLDHE03")) %>%
  left_join(blood_lineages_default(), by = "CellMarker")


cd34_output_plot_lam <- ggplot(to_plot_cd34_output_lam, 
                           aes(x = TimepointMonths, 
                               y = on_g1, fill = CellType, group = CellType,
                               color = CellType, shape = CellMarker)) +
  # stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2)
  stat_smooth(size = 2, level = 0.4, se = T) +
  labs(title = "WAS - CD34 Output all patients",
       y = "% IS Shared with CD34", x = "Months after GT") +
  theme_bw()
ggsave(plot = cd34_output_plot_lam,
       path = sharing_fold,
       filename = paste(date, project_name, "CD34-Output-plot_LAM-only.pdf", sep = "_"),
       width = 9, height = 5)
ggsave(plot = cd34_output_plot_lam,
       path = sharing_fold,
       filename = paste(date, project_name, "CD34-Output-plot_LAM-only.png", sep = "_"),
       width = 9, height = 5)


cd34_output_lam_super <- purrr::map2_df(purity_filtered_lam, names(purity_filtered_lam),
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

cd34_output_lam_super %>%
  select(-is_coord) %>%
  readr::write_tsv(
    file = fs::path(sharing_fold,
                    paste(date, project_name,
                          "CD34-Output_by-SuperGroup_LAM-only.tsv.gz", sep = "_")),
    na = ""
  )

to_plot_cd34_output_lam_super <- cd34_output_lam_super %>%
  tidyr::separate(col = "g2", into = c("SubjectID", "SuperGroup", 
                                       "Tissue", "TimepointMonths"),
                  convert = TRUE) %>%
  # filter(!SubjectID %in% c("MLDC02", "MLDCO2", "MLDCUP01", "MLDCUP02", 
  #                          "MLDCUP03", "MLDCUP04", "MLDCUP05", "MLDHE01",  
  #                          "MLDHE02","MLDHE03")) %>%
  left_join(blood_lineages_default(), by = c("SuperGroup" = "CellMarker"))


cd34_output_plot_lam_super <- ggplot(to_plot_cd34_output_lam_super, 
                               aes(x = TimepointMonths, 
                                   y = on_g1, fill = CellType, group = CellType,
                                   color = CellType)) +
  # stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2)
  stat_smooth(size = 2, level = 0.4, se = T) +
  labs(title = "WAS - CD34 Output all patients",
       y = "% IS Shared with CD34", x = "Months after GT") +
  theme_bw()
ggsave(plot = cd34_output_plot_lam_super,
       path = sharing_fold,
       filename = paste(date, project_name, "CD34-Output-plot_by-SuperGroup_LAM-only.pdf", sep = "_"),
       width = 9, height = 5)
ggsave(plot = cd34_output_plot_lam_super,
       path = sharing_fold,
       filename = paste(date, project_name, "CD34-Output-plot_by-SuperGroup_LAM-only.png", sep = "_"),
       width = 9, height = 5)


## Sharing CD34 vs. Myeloid, Erythro, B, T
### ** Only in vivo tp (no 0); all tp collapsed together;
### ** only PB & BM (together and separated)
cd34_vs_lineages_sepTissue_lam <- purrr::map2_df(purity_filtered_lam,
     names(purity_filtered_lam),
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

cd34_vs_lineages_sepTissue_lam %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-PB-BM_withCellMarker_LAM-only.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_sepTissue_ct_lam <- purrr::map2_df(purity_filtered_lam, 
        names(purity_filtered_lam),
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

cd34_vs_lineages_sepTissue_ct_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-PB-BM_CellTypeOnly_LAM-only.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_sepTissue_sg_lam <- purrr::map2_df(purity_filtered_lam, 
            names(purity_filtered_lam),
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

cd34_vs_lineages_sepTissue_sg_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-PB-BM_SuperGroup_LAM-only.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_sepTissue_ge12_lam <- purrr::map2_df(purity_filtered_lam, 
    names(purity_filtered_lam),
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

cd34_vs_lineages_sepTissue_ct_ge12_lam <- purrr::map2_df(purity_filtered_lam, 
     names(purity_filtered_lam),
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

cd34_vs_lineages_sepTissue_sg_ge12_lam <- purrr::map2_df(purity_filtered_lam, 
   names(purity_filtered_lam),
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

cd34_vs_lineages_sepTissue_ge12_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-PB-BM_withCellMarker_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_sepTissue_ct_ge12_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-PB-BM_CellTypeOnly_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_sepTissue_sg_ge12_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-PB-BM_SuperGroup_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )


cd34_vs_lineages_collTissue_lam <- purrr::map2_df(purity_filtered_lam,
                                              names(purity_filtered_lam),
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

cd34_vs_lineages_collTissue_lam %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-tissueCollapsed_withCellMarker_LAM-only.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_collTissue_ct_lam <- purrr::map2_df(purity_filtered_lam, 
                                                 names(purity_filtered_lam),
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

cd34_vs_lineages_collTissue_ct_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-tissueCollapsed_CellTypeOnly_LAM-only.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_collTissue_sg_lam <- purrr::map2_df(purity_filtered_lam, 
     names(purity_filtered_lam),
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

cd34_vs_lineages_collTissue_sg_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-tissueCollapsed_SuperGroup_LAM-only.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_collTissue_ge12_lam <- purrr::map2_df(purity_filtered_lam,
  names(purity_filtered_lam),
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

cd34_vs_lineages_collTissue_ct_ge12_lam <- purrr::map2_df(purity_filtered_lam, 
  names(purity_filtered_lam),
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

cd34_vs_lineages_collTissue_sg_ge12_lam <- purrr::map2_df(purity_filtered_lam, 
  names(purity_filtered_lam),
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

cd34_vs_lineages_collTissue_ge12_lam %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-tissueCollapsed_withCellMarker_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_collTissue_ct_ge12_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-tissueCollapsed_CellTypeOnly_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )

cd34_vs_lineages_collTissue_sg_ge12_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "CD34-MyErBT-inVivo-tpCollapsed-tissueCollapsed_SuperGroup_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )


## Sharing WHOLE/MNC vs. Myeloid, Erythro, B, T
### ** Only in vivo tp (no 0); all tp collapsed together;
### ** only PB & BM (together and separated)
whole_vs_lineages_sepTissue_lam <- purrr::map2_df(purity_filtered_lam,
  names(purity_filtered_lam),
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

whole_vs_lineages_sepTissue_lam %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_withCellMarker_LAM-only.tsv.gz", sep = "_"))
  )

whole_vs_lineages_sepTissue_ct_lam <- purrr::map2_df(purity_filtered_lam, 
     names(purity_filtered_lam),
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

whole_vs_lineages_sepTissue_ct_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_CellTypeOnly_LAM-only.tsv.gz", sep = "_"))
  )

whole_vs_lineages_sepTissue_sg_lam <- purrr::map2_df(purity_filtered_lam, 
     names(purity_filtered_lam),
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

whole_vs_lineages_sepTissue_sg_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_SuperGroup_LAM-only.tsv.gz", sep = "_"))
  )

whole_vs_lineages_collTissue_lam <- purrr::map2_df(purity_filtered_lam,
     names(purity_filtered_lam),
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

whole_vs_lineages_collTissue_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-tissueCollapsed_withCellMarker_LAM-only.tsv.gz", sep = "_"))
  )

whole_vs_lineages_collTissue_ct_lam <- purrr::map2_df(purity_filtered_lam, 
      names(purity_filtered_lam),
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

whole_vs_lineages_collTissue_ct_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-tissueCollapsed_CellTypeOnly_LAM-only.tsv.gz", sep = "_"))
  )

whole_vs_lineages_collTissue_sg_lam <- purrr::map2_df(purity_filtered_lam, 
      names(purity_filtered_lam),
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

whole_vs_lineages_collTissue_sg_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-tissueCollapsed_SuperGroup_LAM-only.tsv.gz", sep = "_"))
  )

whole_vs_lineages_sepTissue_ge12_lam <- purrr::map2_df(purity_filtered_lam,
 names(purity_filtered_lam),
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

whole_vs_lineages_sepTissue_ct_ge12_lam <- purrr::map2_df(purity_filtered_lam, 
  names(purity_filtered_lam),
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

whole_vs_lineages_sepTissue_sg_ge12_lam <- purrr::map2_df(purity_filtered_lam, 
  names(purity_filtered_lam),
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

whole_vs_lineages_sepTissue_ge12_lam %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_withCellMarker_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )

whole_vs_lineages_sepTissue_ct_ge12_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_CellTypeOnly_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )

whole_vs_lineages_sepTissue_sg_ge12_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-PB-BM_SuperGroup_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )

whole_vs_lineages_collTissue_ge12_lam <- purrr::map2_df(purity_filtered_lam,
  names(purity_filtered_lam),
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

whole_vs_lineages_collTissue_ct_ge12_lam <- purrr::map2_df(purity_filtered_lam, 
   names(purity_filtered_lam),
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

whole_vs_lineages_collTissue_sg_ge12_lam <- purrr::map2_df(purity_filtered_lam, 
   names(purity_filtered_lam),
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

whole_vs_lineages_collTissue_ge12_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-tissueCollapsed_withCellMarker_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )

whole_vs_lineages_collTissue_ct_ge12_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-tissueCollapsed_CellTypeOnly_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )

whole_vs_lineages_collTissue_sg_ge12_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "WholeMNC-MyErBT-inVivo-tpCollapsed-tissueCollapsed_SuperGroup_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )


## NESTED SHARING (cd34 vs WholeMNC) vs. Myeloid, Erythro, B, T
### ** Only in vivo tp (no 0); all tp collapsed together;
### ** only PB & BM (together and separated)
cd34_vs_wholemnc_sepTissue_lam <- purrr::map2(purity_filtered_lam,
    names(purity_filtered_lam),
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

# purrr::reduce(cd34_vs_wholemnc_sepTissue_lam, bind_rows) %>%
#   select(-is_coord) %>%
#   readr::write_tsv(
#     file = fs::path(sharing_fold,
#                     paste(date, project_name,
#                           "CD34-WholeMNC-inVivo-tpColapsed-PB-BM_LAM-only.tsv.gz", sep = "_"))
#   )

shared_nested_sepTissue_lam <- purrr::map2_df(cd34_vs_wholemnc_sepTissue_lam,
    names(cd34_vs_wholemnc_sepTissue_lam),
    ~ {
      print(paste0("Processing ", .y, "..."))
      regenerated <- regen(.x)
      if (is.null(regenerated)) {
        return(NULL)
      }
      other_lineages <- purity_filtered_lam[[.y]] %>%
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

shared_nested_sepTissue_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-PB-BM_withCellMarker_LAM-only.tsv.gz", sep = "_"))
  )

shared_nested_sepTissue_ct_lam <- purrr::map2_df(cd34_vs_wholemnc_sepTissue_lam,
 names(cd34_vs_wholemnc_sepTissue_lam),
 ~ {
   print(paste0("Processing ", .y, "..."))
   regenerated <- regen(.x)
   if (is.null(regenerated)) {
     return(NULL)
   }
   other_lineages <- purity_filtered_lam[[.y]] %>%
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

shared_nested_sepTissue_ct_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-PB-BM_CellTypeOnly_LAM-only.tsv.gz", sep = "_"))
  )

shared_nested_sepTissue_sg_lam <- purrr::map2_df(cd34_vs_wholemnc_sepTissue_lam,
     names(cd34_vs_wholemnc_sepTissue_lam),
     ~ {
       print(paste0("Processing ", .y, "..."))
       regenerated <- regen(.x)
       if (is.null(regenerated)) {
         return(NULL)
       }
       other_lineages <- purity_filtered_lam[[.y]] %>%
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

shared_nested_sepTissue_sg_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-PB-BM_SuperGroup_LAM-only.tsv.gz", sep = "_"))
  )

cd34_vs_wholemnc_sepTissue_ge12_lam <- purrr::map2(purity_filtered_lam,
   names(purity_filtered_lam),
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

shared_nested_sepTissue_ge12_lam <- purrr::map2_df(cd34_vs_wholemnc_sepTissue_ge12_lam,
 names(cd34_vs_wholemnc_sepTissue_ge12_lam),
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

shared_nested_sepTissue_ct_ge12_lam <- purrr::map2_df(cd34_vs_wholemnc_sepTissue_ge12_lam,
  names(cd34_vs_wholemnc_sepTissue_ge12_lam),
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

shared_nested_sepTissue_sg_ge12_lam <- purrr::map2_df(cd34_vs_wholemnc_sepTissue_ge12_lam,
  names(cd34_vs_wholemnc_sepTissue_ge12_lam),
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

shared_nested_sepTissue_ge12_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-PB-BM_withCellMarker_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )

shared_nested_sepTissue_ct_ge12_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-PB-BM_CellTypeOnly_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )

shared_nested_sepTissue_sg_ge12_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-PB-BM_SuperGroup_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )

cd34_vs_wholemnc_collTissue_lam <- purrr::map2(purity_filtered_lam, 
     names(purity_filtered_lam),
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

# purrr::reduce(cd34_vs_wholemnc_collTissue_lam, bind_rows) %>% 
#   select(-is_coord) %>%
#   readr::write_tsv(
#     file = fs::path(sharing_fold,
#                     paste(date, project_name,
#                           "CD34-WholeMNC-inVivo-tpCollapsed-tissueCollapsed_LAM-only.tsv.gz", sep = "_"))
#   )

shared_nested_collTissue_lam <- purrr::map2_df(cd34_vs_wholemnc_collTissue_lam,
     names(cd34_vs_wholemnc_collTissue_lam),
     ~ {
       print(paste0("Processing ", .y, "..."))
       regenerated <- regen2(.x)
       if (is.null(regenerated)) {
         return(NULL)
       }
       other_lineages <- purity_filtered_lam[[.y]] %>%
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

shared_nested_collTissue_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-tissueCollapsed_withCellMarker_LAM-only.tsv.gz", sep = "_"))
  )

shared_nested_collTissue_ct_lam <- purrr::map2_df(cd34_vs_wholemnc_collTissue_lam,
                                             names(cd34_vs_wholemnc_collTissue_lam),
                                             ~ {
                                               print(paste0("Processing ", .y, "..."))
                                               regenerated <- regen2(.x)
                                               if (is.null(regenerated)) {
                                                 return(NULL)
                                               }
                                               other_lineages <- purity_filtered_lam[[.y]] %>%
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

shared_nested_collTissue_ct_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-tissueCollapsed_CellTypeOnly_LAM-only.tsv.gz", sep = "_"))
  )

shared_nested_collTissue_sg_lam <- purrr::map2_df(cd34_vs_wholemnc_collTissue_lam,
                                                 names(cd34_vs_wholemnc_collTissue_lam),
                                                 ~ {
                                                   print(paste0("Processing ", .y, "..."))
                                                   regenerated <- regen2(.x)
                                                   if (is.null(regenerated)) {
                                                     return(NULL)
                                                   }
                                                   other_lineages <- purity_filtered_lam[[.y]] %>%
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

shared_nested_collTissue_sg_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-tissueCollapsed_SuperGroup_LAM-only.tsv.gz", sep = "_"))
  )

cd34_vs_wholemnc_collTissue_ge12_lam <- purrr::map2(purity_filtered_lam,
    names(purity_filtered_lam),
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

shared_nested_collTissue_ge12_lam <- purrr::map2_df(cd34_vs_wholemnc_collTissue_ge12_lam,
  names(cd34_vs_wholemnc_collTissue_ge12_lam),
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

shared_nested_collTissue_ct_ge12_lam <- purrr::map2_df(cd34_vs_wholemnc_collTissue_ge12_lam,
  names(cd34_vs_wholemnc_collTissue_ge12_lam),
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

shared_nested_collTissue_sg_ge12_lam <- purrr::map2_df(cd34_vs_wholemnc_collTissue_ge12_lam,
    names(cd34_vs_wholemnc_collTissue_ge12_lam),
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

shared_nested_collTissue_ge12_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-tissueCollapsed_withCellMarker_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )

shared_nested_collTissue_ct_ge12_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-tissueCollapsed_CellTypeOnly_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )

shared_nested_collTissue_sg_ge12_lam %>%
  select(-truth_tbl_venn) %>%
  readr::write_tsv(na = "",
                   file = fs::path(sharing_fold,
                                   paste(date, project_name,
                                         "nested-CD34|WholeMNC-MyErBT-tpCollapsed-tissueCollapsed_SuperGroup_LAM-only_FROM-TP-12.tsv.gz", sep = "_"))
  )

cd34_vs_all_lineages_sepTissue_ct_lam <- purrr::map2_df(purity_filtered_lam, 
     names(purity_filtered_lam),
     ~ {
       print(paste0("Processing ", .y, "..."))
       invivo <- .x %>% 
         mutate(SubjectID = .y) %>%
         filter(TimepointMonths != "00", Tissue %in% c("PB", "BM")) %>%
         left_join(blood_lineages_default(),  by = "CellMarker")
       cd34_only <- invivo %>%
         filter(CellMarker == "CD34", Tissue == "BM")
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
         if (length(groups_list) < 2) {
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


venns_lam <- sharing_venn(cd34_vs_all_lineages_sepTissue_ct_lam, euler = FALSE)
eulers_lam <- sharing_venn(cd34_vs_all_lineages_sepTissue_ct_lam, euler = TRUE)

venns_lam <- purrr::map(venns_lam, prettify)
eulers_lam <- purrr::map(eulers_lam, prettify)

for (venn in venns_lam) {
  colors <- tibble::as_tibble(list(marker = rownames(venn$ellipses))) %>%
    left_join(colorscale) %>%
    pull(color)
  pdf(file = fs::path(sharing_fold, 
                      paste0(venn$Patient,
                             "_", venn$Tissue, ".",
                             paste(date, project_name, 
                                   "sharing_venn_LAM-only.pdf", 
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

for (eul in eulers_lam) {
  colors <- tibble::as_tibble(list(marker = rownames(eul$ellipses))) %>%
    left_join(colorscale) %>%
    pull(color)
  pdf(file = fs::path(sharing_fold, 
                      paste0(eul$Patient,
                             "_", eul$Tissue, ".",
                             paste(date, project_name, 
                                   "sharing_euler_LAM-only.pdf", 
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

## HSC population size estimate
HSC_estimate_lam <- HSC_population_size_estimate(
  x = agg_lam, 
  metadata = agg_meta_lam,
  aggregation_key = agg_key,
  timepoint_column = "TimepointMonths", 
  stable_timepoints = stable_tps
)

HSC_estimate_lam %>%
  readr::write_tsv(file = fs::path(project_folder, 
                                   paste(date, project_name, 
                                         "HSC_population_size_estimate_LAM-only.tsv", 
                                         sep = "_"), na = ""))

hsc_plot_lam <- HSC_population_plot(HSC_estimate_lam, project_name = "WAS")
hsc_plot_lam <- hsc_plot_lam +
  theme_bw()
ggsave(plot = hsc_plot_lam, path = plots_fold,
       paste(date, project_name, 
             "HSC_population_size_estimate_plot_LAM-only.pdf", sep = "_"),
       width = 8, height = 5)
ggsave(plot = hsc_plot_lam, path = plots_fold,
       paste(date, project_name, 
             "HSC_population_size_estimate_plot_LAM-only.png", sep = "_"),
       width = 8, height = 5)

# SOURCE OF ISS
agg_patient_split_lam <- agg_lam %>%
  filter(SubjectID %in% patients) %>%
  group_by(SubjectID) %>%
  group_split()

source_of_iss_lam <- list()

for (marker in selected_markers) {
  for (patient_df in agg_patient_split_lam) {
    print(paste0("Processing marker ", marker, " for patient ", 
                 patient_df$SubjectID[1]))
    patient_df <- patient_df %>%
      filter(CellMarker == marker)
    if (nrow(patient_df) > 0) {
      iss_s <- iss_source(reference = patient_df, selection = patient_df, 
                          ref_group_key = agg_key, selection_group_key = agg_key, 
                          timepoint_column = "TimepointMonths")
      patient_source <- source_of_iss_lam[[patient_df$SubjectID[1]]]
      if (is.null(patient_source)) {
        source_of_iss_lam[[patient_df$SubjectID[1]]] <- iss_s[[1]]
      } else {
        source_of_iss_lam[[patient_df$SubjectID[1]]] <- bind_rows(
          source_of_iss_lam[[patient_df$SubjectID[1]]], iss_s[[1]])
      }
    }
  }
}

source_all_table_lam <- purrr::reduce(source_of_iss_lam, bind_rows)
source_all_table_lam <- source_all_table_lam %>%
  filter(g1_Tissue == g2_Tissue)

source_all_table_lam %>%
  select(-is_coord) %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         "source_of_iss_LAM-only.tsv.gz",
                                         sep = "_")))

agg_sg_per_patient_lam <- agg_supergroup_lam %>%
  filter(as.numeric(TimepointMonths) > 0, SubjectID %in% patients) %>%
  group_by(SubjectID) %>%
  group_split()

source_of_iss_sg_lam <- list()
for (marker in selected_supergroups) {
  for (patient_df in agg_sg_per_patient_lam) {
    print(paste0("Processing marker ", marker, " for patient ", 
                 patient_df$SubjectID[1]))
    patient_df <- patient_df %>%
      filter(SuperGroup == marker)
    if (nrow(patient_df) > 0) {
      iss_s <- iss_source(reference = patient_df, selection = patient_df, 
                          ref_group_key = agg_key_supergroup, 
                          selection_group_key = agg_key_supergroup, 
                          timepoint_column = "TimepointMonths")
      patient_source <- source_of_iss_sg_lam[[patient_df$SubjectID[1]]]
      if (is.null(patient_source)) {
        source_of_iss_sg_lam[[patient_df$SubjectID[1]]] <- iss_s[[1]]
      } else {
        source_of_iss_sg_lam[[patient_df$SubjectID[1]]] <- bind_rows(
          source_of_iss_sg_lam[[patient_df$SubjectID[1]]], iss_s[[1]])
      }
    }
  }
}

source_all_table_sg_lam <- purrr::reduce(source_of_iss_sg_lam, bind_rows)
source_all_table_sg_lam <- source_all_table_sg_lam %>%
  filter(g1_Tissue == g2_Tissue)

source_all_table_sg_lam %>%
  select(-is_coord) %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         "source_of_iss_SuperGroup_LAM-only.tsv.gz",
                                         sep = "_")))

source_of_iss_sg_purityfilt_lam <- list()
purity_mod_lam <- purrr::map2(purity_filtered_lam, names(purity_filtered_lam), ~ {
  .x %>%
    mutate(SubjectID = .y) %>%
    left_join(blood_lineages_default(), by = "CellMarker")
})

for (marker in selected_supergroups) {
  for (patient_df in purity_mod_lam) {
    print(paste0("Processing marker ", marker, " for patient ", 
                 patient_df$SubjectID[1]))
    patient_df <- patient_df %>%
      filter(SuperGroup == marker, as.numeric(TimepointMonths) > 0)
    if (nrow(patient_df) > 0) {
      iss_s <- iss_source(reference = patient_df, selection = patient_df, 
                          ref_group_key = agg_key_supergroup, 
                          selection_group_key = agg_key_supergroup, 
                          timepoint_column = "TimepointMonths")
      patient_source <- source_of_iss_sg_purityfilt_lam[[patient_df$SubjectID[1]]]
      if (is.null(patient_source)) {
        source_of_iss_sg_purityfilt_lam[[patient_df$SubjectID[1]]] <- iss_s[[1]]
      } else {
        source_of_iss_sg_purityfilt_lam[[patient_df$SubjectID[1]]] <- bind_rows(
          source_of_iss_sg_purityfilt_lam[[patient_df$SubjectID[1]]], iss_s[[1]])
      }
    }
  }
}

source_all_table_sg_purityfilt_lam <- purrr::reduce(source_of_iss_sg_purityfilt_lam, bind_rows)
source_all_table_sg_purityfilt_lam <- source_all_table_sg_purityfilt_lam %>%
  filter(g1_Tissue == g2_Tissue)

source_all_table_sg_purityfilt_lam %>%
  select(-is_coord) %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste(date, project_name,
                                         "source_of_iss_SuperGroup_purityfilt_LAM-only.tsv.gz",
                                         sep = "_")))


#------------------------------------------------------------------------------
# WHOLE + MNC Workflow
## ** In these steps MNC cell marker is replaced both in matrix and af by 
## ** Whole and sensitive steps are re-run.
## ** NOTE: all sharing steps are NOT affected by this since a re-aggregation
## ** step is performed prior calculations. Also not affected: CIS

# SUBSTITUTION
no_mnc_af_lam <- outliers_removed_lam %>%
  mutate(CellMarker = if_else(CellMarker == "MNC", "Whole", CellMarker))
no_mnc_agg_lam <- aggregate_values_by_key(final_matrix_lam,
                                      association_file = no_mnc_af_lam, 
                                      value_cols = c("seqCount"),
                                      key = agg_key)
no_mnc_agg_meta_lam <- aggregate_metadata(no_mnc_af_lam, grouping_keys = agg_key)

# STATS
no_mnc_stats_lam <- sample_statistics(no_mnc_agg_lam, no_mnc_agg_meta_lam,
                                  sample_key = agg_key,
                                  value_columns = c("seqCount_sum"))

# Filtering: we apply a filter for IS 
# with sc < 20 (only on markers that are NOT 'Plasma')
to_remove_lam <- no_mnc_stats_lam$metadata %>%
  filter(nIS < 20, CellMarker != "Plasma")

no_mnc_agg_filtered_lam <- no_mnc_agg_lam %>%
  anti_join(to_remove_lam, by = agg_key)

no_mnc_stats_lam <- sample_statistics(no_mnc_agg_filtered_lam, no_mnc_agg_meta_lam,
                                  sample_key = agg_key,
                                  value_columns = c("seqCount_sum"))

no_mnc_stats_lam$metadata %>%
  readr::write_tsv(file = fs::path(project_folder,
                                   paste0(stats_file_prefix, 
                                          "_",
                                          concat_agg_key,
                                          "_NO-MNC_LAM-only.tsv.gz")),
                   na = "")

# ABUNDANCE
no_mnc_abund_lam <- compute_abundance(no_mnc_agg_filtered_lam, 
                                      key = agg_key, 
                                      columns = "seqCount_sum")
readr::write_tsv(no_mnc_abund_lam, 
                 file = fs::path(project_folder, 
                                 paste0(paste(date, project_name,
                                              "abundance", 
                                              concat_agg_key,
                                              sep = "_"), "_NO-MNC_LAM-only.tsv.gz")
                 ),
                 na = "")

# ALLUVIAL - plotting Whole only for each patient
alluv_no_mnc_lam <- integration_alluvial_plot(x = no_mnc_abund_lam %>%
                                     filter(SubjectID %in% patients,
                                            CellMarker == "Whole"), 
                                   top_abundant_tbl = TRUE, 
                                   alluvia_plot_y_threshold = 1,
                                   plot_x = "TimepointMonths", 
                                   plot_y = "seqCount_sum_PercAbundance")
purrr::walk2(alluv_no_mnc_lam, names(alluv_no_mnc_lam), ~ {
  prettified_id <- stringr::str_replace_all(.y, "_", ", ")
  p <- .x$plot +
    theme_bw() +
    theme(legend.position = "none") +
    ggplot2::labs(title = paste0("WAS - IS Abundance over time - ",
                                 prettified_id),
                  subtitle = paste0("Colored flows indicate IS having ",
                                    "a relative abundance >= 1%\n",
                                    "in at least 1 time point"),
                  x = "Months after GT",
                  y = "Abundace (%)")
  file_prefix <- paste0(.y, ".",
                        paste(date, project_name, sep = "_"))
  ggplot2::ggsave(plot = p, path = alluv_fold, 
                  filename = paste0(file_prefix, "_alluvial-plot_NO-MNC_LAM-only.pdf"), 
                  width = 8, height = 5)
  ggplot2::ggsave(plot = p, path = alluv_fold, 
                  filename = paste0(file_prefix, "_alluvial-plot_NO-MNC_LAM-only.png"), 
                  width = 8, height = 5)
  pdf(file = fs::path(alluv_fold, 
                      paste0(file_prefix, "_top-abund-tbl_NO-MNC_LAM-only.pdf")), 
      width = 13, height = 7)
  gridExtra::grid.arrange(.x$tables)
  dev.off()
})

# DIVERSITY PLOT - WHOLE/MNC
h_index_data_lam <- no_mnc_stats_lam$metadata %>%
  mutate(TimepointMonths = as.numeric(TimepointMonths)) %>%
  filter(CellMarker %in% c("Whole"),
         TimepointMonths > 0,
         SubjectID %in% patients,
         Tissue %in% c("PB", "BM")) %>%
  select(SubjectID, Tissue, TimepointMonths, 
         seqCount_sum_shannon)

h_index_plot_lam <- ggplot(data = h_index_data_lam, aes(x = TimepointMonths,
                                                y = seqCount_sum_shannon,
                                                group = SubjectID, 
                                                color = SubjectID)) +
  geom_point(size = 2.5) +
  geom_line(size = 1.5, alpha = .7) +
  facet_wrap(~Tissue) +
  # geom_smooth(aes(group=Tissue)) +
  theme_bw() +
  labs(y = "H index", x = "Time point (months after GT)",
       title = "WAS Clonal population diversity over time - Whole/MNC") +
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

ggsave(plot = h_index_plot_lam, filename = paste0(paste(date, project_name, 
                                                    "hindex-plot-WHOLE-MNC-withSmooth_LAM-only",
                                                    sep = "_"), ".pdf"),
       width = 13, height = 8, path = plots_fold)
ggsave(plot = h_index_plot_lam, filename = paste0(paste(date, project_name, 
                                                    "hindex-plot-WHOLE-MNC-withSmooth_LAM-only",
                                                    sep = "_"), ".png"),
       width = 13, height = 8, path = plots_fold)
ggsave(plot = h_index_plot_lam, filename = paste0(paste(date, project_name, 
                                                    "hindex-plot-WHOLE-MNC-NOSmooth_LAM-only",
                                                    sep = "_"), ".pdf"),
       width = 13, height = 8, path = plots_fold)
ggsave(plot = h_index_plot_lam, filename = paste0(paste(date, project_name, 
                                                    "hindex-plot-WHOLE-MNC-NOSmooth_LAM-only",
                                                    sep = "_"), ".png"),
       width = 13, height = 8, path = plots_fold)

# CUMULATED MATRICES
cumulated_no_mnc_lam <- cumulative_is(no_mnc_agg_filtered_lam,
                                  key = agg_key,
                                  timepoint_col = "TimepointMonths",
                                  keep_og_is = FALSE,
                                  expand = TRUE)
cumulated_per_patient_no_mnc_lam <- cumulated_no_mnc_lam %>%
  group_by(SubjectID) %>%
  group_split()
for (df in cumulated_per_patient_no_mnc_lam) {
  subj <- df$SubjectID[1]
  print(paste("Saving:", subj, "..."))
  file_name <- paste0(subj, ".",
                      paste(date, project_name, 
                            paste0(agg_key, collapse = "_"),
                            sep = "_"),
                      "_cumulated_is_NO-MNC_LAM-only.tsv.gz")
  df %>%
    readr::write_tsv(file = fs::path(lam_cumulative_folder, file_name), 
                     na = "")
}

rm(cumulated_per_patient_no_mnc_lam, cumulated_no_mnc_lam, subj, file_name, df)
gc()

#-----------------------------------------------------------------------------
# Analyses with timeframes instead of timepoints
#-----------------------------------------------------------------------------
af_timeframes <- af_blood_l_lam %>%
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

timeframes_key <- c("SubjectID", "SuperGroup", "Tissue", "TimeFrame")

agg_timeframes <- aggregate_values_by_key(
  final_matrix_lam, 
  association_file = af_timeframes, 
  value_cols = "seqCount",
  key = timeframes_key
)

concat_tf_key <- paste0(timeframes_key, collapse = "_")

matrix_file_prefix_tf <- paste(date, project_name, 
                               concat_tf_key,
                                "LAM-only",
                                sep = "_")

per_patient_sparse(agg_timeframes,
                   quantific = c(seqCount = "seqCount_sum"),
                   prefix = matrix_file_prefix_tf,
                   dir_name = lam_matrix_agg, 
                   row_totals = TRUE)

matrix_file_prefix_tf_filt <- paste(date, project_name, 
                                    concat_tf_key,
                                    "LAM-only_SC-filtered",
                                    sep = "_")
per_patient_sparse(agg_timeframes %>% filter(seqCount_sum >= 3),
                   quantific = c(seqCount = "seqCount_sum"),
                   prefix = matrix_file_prefix_tf_filt,
                   dir_name = lam_matrix_agg, 
                   row_totals = TRUE)

cd34_output_timeframe <- purrr::map2_df(purity_filtered_lam, 
                                        names(purity_filtered_lam),
  ~ {
    temp <- .x %>% mutate(SubjectID = .y) %>%
      dplyr::left_join(af_timeframes %>%
                         select(SubjectID, CellMarker, Tissue, TimepointMonths,
                                SuperGroup, TimeFrame, HematoLineage), 
                       by = c("SubjectID", "CellMarker", 
                               "Tissue", "TimepointMonths"))
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
                                           g2 = c("SubjectID","SuperGroup", "Tissue", "TimeFrame")),
                         keep_genomic_coord = TRUE)
    shared
  })

cd34_output_timeframe %>%
  select(-is_coord) %>%
  readr::write_tsv(
    file = fs::path(sharing_fold,
                    paste(date, project_name,
                          "CD34-Output_by-SuperGroup-TimeFrame_LAM-only.tsv.gz", sep = "_")),
    na = ""
  )

to_plot_cd34_output_tf <- cd34_output_timeframe %>%
  tidyr::separate(col = "g2", into = c("SubjectID", "SuperGroup", 
                                       "Tissue", "TimeFrame"),
                  convert = TRUE) %>%
  left_join(blood_lineages_default(), by = c("SuperGroup" = "CellMarker"))


cd34_output_plot_tf <- ggplot(to_plot_cd34_output_tf %>%
                                filter(!is.na(TimeFrame)), 
                               aes(x = TimeFrame, 
                                   y = on_g1, fill = CellType, group = CellType,
                                   color = CellType)) +
  # stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2)
  stat_smooth(size = 2, level = 0.4, se = T) +
  labs(title = "WAS - CD34 Output all patients",
       y = "% IS Shared with CD34", x = "Months after GT") +
  theme_bw()
ggsave(plot = cd34_output_plot_tf,
       path = sharing_fold,
       filename = paste(date, project_name, "CD34-Output-plot_by-SuperGroup-TimeFrame_LAM-only.pdf", sep = "_"),
       width = 9, height = 5)
ggsave(plot = cd34_output_plot_tf,
       path = sharing_fold,
       filename = paste(date, project_name, "CD34-Output-plot_by-SuperGroup-TimeFrame_LAM-only.png", sep = "_"),
       width = 9, height = 5)

