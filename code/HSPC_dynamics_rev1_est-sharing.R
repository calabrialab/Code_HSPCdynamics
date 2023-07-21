################################################################################
# HSPC Dynamics additional analyses - Estimates and sharing
################################################################################
# 2023-07-21
library(ISAnalytics)
library(dplyr)
library(ggplot2)

progressr::handlers(global = TRUE)
progressr::handlers("cli")

proj_folder <- fs::path_wd()
est_folder <- fs::path(proj_folder, "population_estimates")
sharing_folder <- fs::path(proj_folder, "sharing")

mld_agg_af_path <- fs::path("/mnt/CerberoWorkspaceISA/MLD/archive/202112",
                            "20211202_MLD_descriptive-stats_SubjectID_CellMarker_Tissue_TimepointMonthsSC-filtered.tsv.gz")
mld_agg_af <- readr::read_tsv(mld_agg_af_path, guess_max = Inf)

bthal_agg_af_path <- fs::path("/mnt/CerberoWorkspaceISA/BTHAL/archive/202112",
                              "20211206_BTHAL_descriptive-stats_SubjectID_CellMarker_Tissue_TimepointMonths_SC-filtered.tsv.gz")
bthal_agg_af <- readr::read_tsv(bthal_agg_af_path, guess_max = Inf)

was_agg_af_path <- fs::path("/mnt/CerberoWorkspaceISA/WAS/archive/202112",
                            "20211210_WAS_descriptive-stats_SubjectID_CellMarker_Tissue_TimepointMonths.tsv.gz")
was_agg_af <- readr::read_tsv(was_agg_af_path, guess_max = Inf)

included_markers <- c(
  "CD34",
  "CD3",
  "CD4",
  "CD8",
  "CD19",
  "CD36",
  "CD38",
  "GLY",
  "CD13",
  "CD14",
  "CD15",
  "WHOLE",
  "MNC"
)

# Population estimates ---------------------------------------------------------
datasets_for_popul <- list(
  MLD = mld_data |>
    filter(as.numeric(TimepointMonths) > 0,
           CellMarker %in% included_markers),
  BTHAL = bthal_data |>
    filter(as.numeric(TimepointMonths) > 0,
           CellMarker %in% included_markers),
  WAS = was_data |>
    filter(as.numeric(TimepointMonths) > 0,
           CellMarker %in% included_markers)
)

afs <- list(
  MLD = mld_agg_af,
  BTHAL = bthal_agg_af,
  WAS = was_agg_af
)

agg_key <- c("SubjectID", "CellMarker", "Tissue", "TimepointMonths")

HSC_estimate_ge24 <- purrr::map2(datasets_for_popul, afs,
                                ~ HSC_population_size_estimate(
                                  x = .x |>
                                    filter(as.numeric(TimepointMonths) >= 24),
                                  metadata = .y, 
                                  stable_timepoints = seq(24, 100), 
                                  aggregation_key = agg_key,
                                  timepoint_column = "TimepointMonths", 
                                  seqCount_column = "seqCount",
                                  fragmentEstimate_column = NULL
                                ))
HSC_estimate_less24 <- purrr::map2(datasets_for_popul, afs,
                                   ~ HSC_population_size_estimate(
                                     x = .x |>
                                       filter(as.numeric(TimepointMonths) < 24),
                                     metadata = .y, 
                                     stable_timepoints = seq(9, 12), 
                                     aggregation_key = agg_key,
                                     timepoint_column = "TimepointMonths", 
                                     seqCount_column = "seqCount",
                                     fragmentEstimate_column = NULL
                                   ))

purrr::iwalk(HSC_estimate_ge24, ~ {
  file_name <- paste0(.y, ".", "20230721_HSPC-estimates-TPge24.tsv")
  .x$est |>
    readr::write_tsv(
      file = fs::path(est_folder, file_name),
      na = ""
    )
})

purrr::iwalk(HSC_estimate_less24, ~ {
  file_name <- paste0(.y, ".", "20230721_HSPC-estimates-TPless24.tsv")
  .x$est |>
    readr::write_tsv(
      file = fs::path(est_folder, file_name),
      na = ""
    )
})

# Sharing ----------------------------------------------------------------------
sharing_same_p <- function(df, patient) {
  df_less24 <- df |>
    filter(as.numeric(TimepointMonths) < 24,
           SubjectID == patient)
  df_ge24 <- df |>
    filter(as.numeric(TimepointMonths) >= 24,
           SubjectID == patient)
  shared <- is_sharing(
    df_less24, df_ge24, group_key = "SubjectID", keep_genomic_coord = TRUE,
    table_for_venn = TRUE
  )
  return(shared)
}

sharing_results <- purrr::map(datasets_for_popul, ~ {
  df <- .x
  patients <- unique(.x$SubjectID)
  purrr::map(patients, ~ sharing_same_p(df, .x)) |>
    purrr::list_rbind()
})

purrr::iwalk(sharing_results, ~ {
  file_name <- paste0(.y, ".", "20230721_sharing-TP24.tsv")
  .x |>
    select(-is_coord, -truth_tbl_venn) |>
    readr::write_tsv(
      file = fs::path(sharing_folder, file_name),
      na = ""
    )
})

save.image(file = "./HSPC_dynamics_rev1_est-sharing.RData")
