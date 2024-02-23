################################################################################
# HSPC Dynamics additional analyses - XSCID dataset
################################################################################
# 2023-08-29
library(ISAnalytics)
library(dplyr)
library(openxlsx)
library(ggplot2)

source("../per_patient_matrix_save.R")

progressr::handlers(global = TRUE)
progressr::handlers("cli")

proj_folder <- fs::path_wd()
xscid_fold <- fs::path(proj_folder, "xscid_dataset")
fs::dir_create(xscid_fold)

data_path <- fs::path("/giuliahome/HSPC_dynamics_data",
                      "41467_2022_31344_MOESM6_ESM.xlsx")

xscid_dataset <- purrr::map(3:9, ~ {
  read.xlsx(data_path, sheet = .x, startRow = 2) |>
    tibble::as_tibble() |>
    mutate(
      patient = stringr::str_remove(getSheetNames(data_path)[.x], 
                                    "Combined_allsites")
    )
}) |> purrr::list_rbind()

xscid_dataset_melt <- xscid_dataset |>
  select(-Gene_Start, -Gene_end,
         -distance_to_gene_or_exon_intron, -cumulative.frequency,
         -UCSC_link, -in_out_gene, -RefSeq_ID) |>
  tidyr::pivot_longer(cols = !c(chr, strand, bp, patient, Gene, Gene_Strand), 
                      names_to = "sample", values_to = "seqCount", 
                      values_drop_na = TRUE)
xscid_dataset_melt <- xscid_dataset_melt |>
  rename(integration_locus = "bp", GeneName = "Gene",
         GeneStrand = "Gene_Strand") |>
  tidyr::separate(col = sample, sep = "_", 
                  into = c("CellMarker", "TimepointMonths"))
xscid_dataset_melt <- xscid_dataset_melt |>
  mutate(chr = stringr::str_remove(chr, "chr"),
         integration_locus = as.integer(integration_locus),
         CellMarker = stringr::str_replace(CellMarker, "-", ""),
         TimepointMonths = stringr::str_remove(TimepointMonths, "m(R)*(S)*")
         )

# SHARING ----
## CD34 output: for each patient, apply the purity filter and compute sharing
## between CD34+ and lineages Erythroid, Myeloid and Lymphoid.
## Purity filter is applied by timepoint only on groups of interest

groups_of_interest <- c("CD34", "CD13", "CD14", "CD15",
                        "CD19", "CD3", "CD4", "CD8", 
                        "CD36", "GLY", "GLYA")

filtered_ds <- xscid_dataset_melt |>
  filter(seqCount >= 3) |>
  rename(SubjectID = "patient") |>
  group_by(chr, strand, integration_locus, GeneName, GeneStrand, 
           SubjectID, CellMarker, TimepointMonths) |>
  summarise(seqCount = sum(seqCount), .groups = "drop")

per_patient_sparse(filtered_ds, dir_name = fs::path(xscid_fold, 
                                                    "matrices",
                                                    "CellMarker_TimepointMonths"), 
                   prefix = paste("20230830", "HSPC-dynamics",
                                  "xscid-dataset",
                                  "CellMarker_TimepointMonths",
                                  "sc-filt3",
                                  sep = "_"))

filtered_ds_supergroup <- filtered_ds |>
  left_join(blood_lineages_default() |> select(CellMarker, SuperGroup),
            by = "CellMarker") |>
  group_by(chr, strand, integration_locus, GeneName, GeneStrand, 
           SubjectID, SuperGroup, TimepointMonths) |>
  summarise(seqCount = sum(seqCount), .groups = "drop")

per_patient_sparse(filtered_ds_supergroup, dir_name = fs::path(xscid_fold, 
                                                    "matrices",
                                                    "SuperGroup_TimepointMonths"), 
                   prefix = paste("20230830", "HSPC-dynamics",
                                  "xscid-dataset",
                                  "SuperGroup_TimepointMonths",
                                  "sc-filt3",
                                  sep = "_"))

agg_by_patient <- filtered_ds |>
  group_by(SubjectID) |>
  group_split()

subj_names <- purrr::map_chr(agg_by_patient, ~ .x$SubjectID[1])

purity_filtered <- purrr::map(agg_by_patient, 
                              ~ purity_filter(.x, 
                                              selected_groups = groups_of_interest,
                                              group_key = "CellMarker",
                                              aggregation_key = c("SubjectID",
                                                                  "CellMarker", "TimepointMonths"), 
                                              timepoint_column = "TimepointMonths", 
                                              value_column = "seqCount"))
purity_filtered <- purity_filtered |> purrr::set_names(subj_names)

purity_filtered_all <- purrr::imap(purity_filtered, ~ .x |>
                                     mutate(SubjectID = .y)) |>
  purrr::list_rbind()

# filter_stats <- purrr::map2(purity_filtered, names(purity_filtered), ~{
#   perc_1 <- (.x |>
#     filter(Value >= 1) |>
#     nrow()) / nrow(.x)
#   perc_1_2 <- (.x |>
#                  filter(Value >= 2) |>
#                  nrow()) / nrow(.x)
#   perc_1_2_3 <- (.x |>
#                    filter(Value >= 3) |>
#                    nrow()) / nrow(.x)
#   perc_1_CD34 <- (.x |>
#                     filter(Value >= 1, CellMarker == "CD34") |>
#                     nrow()) / nrow(.x |>
#                                      filter(CellMarker == "CD34"))
#   perc_1_2_CD34 <- (.x |>
#                     filter(Value >= 2, CellMarker == "CD34") |>
#                     nrow()) / nrow(.x |>
#                                      filter(CellMarker == "CD34"))
#   perc_1_2_3_CD34 <- (.x |>
#                       filter(Value >= 3, CellMarker == "CD34") |>
#                       nrow()) / nrow(.x |>
#                                        filter(CellMarker == "CD34"))
#   subj <- .y
#   tibble::tribble(
#     ~ Patient, ~ Perc_sc_1, ~ Perc_sc_ge_2, ~ Perc_sc_ge_3, ~ Perc_CD34_sc_1,
#     ~ Perc_CD34_sc_ge_2, ~ Perc_CD34_sc_ge_3,
#     subj, perc_1*100, perc_1_2*100, perc_1_2_3*100, perc_1_CD34*100,
#     perc_1_2_CD34*100, perc_1_2_3_CD34*100
#   )
# }) |> purrr::list_rbind()
# 
# filter_stats_toplot <- filter_stats |>
#   tidyr::pivot_longer(cols = !c(Patient), names_to = "type", 
#                       values_to = "value")
# 
# plot_filter_stats <- ggplot(filter_stats_toplot, aes(x = Patient, y = value, fill = type)) +
#   geom_col(position = "dodge") +
#   theme_bw() +
#   labs(y = "%") +
#   scale_fill_discrete(
#     name = "", 
#     labels = c("% CD34 with SC >= 1",
#                "% CD34 with SC >= 2",
#                "% CD34 with SC >= 3",
#                "% Overall with SC >= 1",
#                "% Overall with SC >= 2",
#                "% Overall with SC >= 3")
#   )


purrr::walk2(purity_filtered, names(purity_filtered), ~ {
  .x |>
    readr::write_tsv(file = fs::path(
      fs::path(xscid_fold, 
               "matrices",
               "Purity-sc-filtered"),
      paste0(
        .y, ".",
        paste("20230830", "HSPC-dynamics", "xscid-dataset",
              "CellMarker_TimepointMonths",
              "purity-filtered", sep = "_"),
        ".tsv.gz"
      )
    ), na = "")
})

cd34_output_super <- purrr::map2_df(purity_filtered, names(purity_filtered),
  ~ {
    temp <- .x |> mutate(SubjectID = .y) %>%
      dplyr::left_join(blood_lineages_default(), 
                       by = "CellMarker")
    cd34_only <- temp |>
      filter(CellMarker == "CD34")
    e_l_m <- temp |>
      filter(HematoLineage %in% c("Erythroid",
                                  "Myeloid",
                                  "Lymphoid"))
    if (nrow(cd34_only) == 0 || nrow(e_l_m) == 0) {
      return(NULL)
    }
    shared <- is_sharing(cd34_only, e_l_m,
                         group_keys = list(g1 = c("SubjectID", "SuperGroup"),
                                           g2 = c("SubjectID","SuperGroup", "TimepointMonths")),
                         keep_genomic_coord = TRUE)
    shared
  })

cd34_output_super <- cd34_output_super |>
  mutate(
    g1 = stringr::str_remove(g1, "-1$"),
    g2 = stringr::str_remove(g2, "-2$")
  )

cd34_output_super |>
  select(-is_coord) |>
  readr::write_tsv(
    file = fs::path(xscid_fold,
                    "cd34_output",
                    paste("20230830", "HSPC-dynamics", "xscid-dataset",
                          "CD34-Output_by-SuperGroup.tsv.gz", sep = "_")),
    na = ""
  )

# to_plot_cd34_output_super <- cd34_output_super %>%
#   tidyr::separate(col = "g2", into = c("SubjectID", "SuperGroup", "TimepointMonths", "n"),
#                   convert = TRUE) %>%
#   select(-n) |>
#   left_join(blood_lineages_default(), by = c("SuperGroup" = "CellMarker"))
# 
# 
# cd34_output_plot_super <- ggplot(to_plot_cd34_output_super, 
#                                  aes(x = TimepointMonths, 
#                                      y = on_g1, fill = CellType, group = CellType,
#                                      color = CellType)) +
#   # stat_smooth(aes(fill = CellType), formula = y ~ log(x), level = 0.4, size = 2)
#   stat_smooth(linewidth = 2, level = 0.4, se = T) +
#   labs(title = "XSCID - CD34 Output all patients",
#        y = "% IS Shared with CD34", x = "Months after GT") +
#   theme_bw()

joined_data <- filtered_ds |>
  full_join(purity_filtered_all) |>
  left_join(blood_lineages_default(), by = "CellMarker") |>
  filter(!is.na(Value), !is.na(SuperGroup)) |>
  group_by(chr, integration_locus, strand, GeneName, GeneStrand,
           SubjectID, SuperGroup, TimepointMonths) |>
  summarise(
    seqCount = sum(seqCount),
    .groups = "drop"
  )

cd34_output_iss_details <- joined_data  |>
  semi_join(
    cd34_output_super |>
      select(g1, is_coord) |>
      tidyr::unnest(is_coord) |>
      distinct() |>
      tidyr::separate(col = g1, 
                      into = c("SubjectID", "SuperGroup"))
  ) |>
  bind_rows(
    joined_data |>
      semi_join(
        cd34_output_super |>
          select(g2, is_coord) |>
          tidyr::unnest(is_coord) |>
          distinct() |>
          tidyr::separate(col = g2, 
                          into = c("SubjectID", "SuperGroup",
                                   "TimepointMonths"))
      )
  ) |> 
  distinct()

cd34_output_iss_details |>
  readr::write_tsv(
    file = fs::path(save_folder, "XSCID", 
                    "20231006_XSCID_CD34-output-IS-selection.tsv.gz")
  )

## STATS ----
agg_stats <- sample_statistics(
  x = filtered_ds,
  metadata = filtered_ds |> distinct(SubjectID, CellMarker, TimepointMonths),
  value_columns = "seqCount",
  sample_key = c("SubjectID", "CellMarker", "TimepointMonths")
)

agg_stats$metadata |>
  readr::write_tsv(
    file = fs::path(xscid_fold,
                    "20231019_XSCID_descriptive-stats_SubjectID-CellMarker-TimepointMonths-SC3.tsv"),
    na = ""
  )

agg_stats_su <- sample_statistics(
  x = filtered_ds_supergroup,
  metadata = filtered_ds_supergroup |> distinct(SubjectID, SuperGroup, TimepointMonths),
  value_columns = "seqCount",
  sample_key = c("SubjectID", "SuperGroup", "TimepointMonths")
)

agg_stats_su$metadata |>
  readr::write_tsv(
    file = fs::path(xscid_fold,
                    "20231019_XSCID_descriptive-stats_SubjectID-SuperGroup-TimepointMonths-SC3.tsv"),
    na = ""
  )

## ABUNDANCE ----
xscid_abundance_classic <- compute_abundance(filtered_ds,
                                             key = c("SubjectID", "CellMarker",
                                                     "TimepointMonths"),
                                             columns = "seqCount")

xscid_abundance_classic |>
  readr::write_tsv(
    file = fs::path(xscid_fold, 
                    paste("XSCID.20231030",
                          "abundance-IS-selection",
                          "SubjectID_CellMarker_TimepointMonths.tsv.gz",
                          sep = "_"))
  )

save.image(file = "./HSPC_dynamics_rev1_xscid.RData")




