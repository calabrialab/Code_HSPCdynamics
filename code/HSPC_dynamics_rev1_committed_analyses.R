################################################################################
# HSPC Dynamics additional analyses - Committed analyses
################################################################################
# 2023-09-05
library(ISAnalytics)
library(dplyr)
library(ggplot2)
library(stringr)

progressr::handlers(global = TRUE)
progressr::handlers("cli")

proj_folder <- fs::path_wd()
comm_fold <- fs::path(proj_folder, "committed_analyses")
for (fold in c(comm_fold)) {
  fs::dir_create(fold)
}
sharing_folder <- fs::path(proj_folder, "sharing")

load("/projects/HSPC_dynamics_paper/HSPC_dynamics_rev1_abundance.RData")

rm(list = ls()[!ls() %in% c("mld_flag_data", "bthal_flag_data", 
                            "was_flag_data")])

flag_data <- list(
  mld = mld_flag_data,
  bthal = bthal_flag_data,
  was = was_flag_data
)

is_per_late_label <- purrr::map(flag_data, ~ .x |> 
                                  filter(Late_StateLabel != "Multi") |> 
                                  group_by(Late_StateLabel) |> 
                                  summarise(n = n_distinct(chr, 
                                                           integration_locus, 
                                                           strand)))

SAMPLE_SIZE <- 80
REPS <- 1000
set.seed(456893)

# Functions --------------------------------------------------------------------
extract_data_test <- function(dataset) {
  is_distinct_per_gene <- dataset |>
    group_by(GeneName) |>
    summarise(nIS = n_distinct(chr, integration_locus, strand))
  great_tot <- sum(is_distinct_per_gene$nIS)
  is_distinct_per_gene <- is_distinct_per_gene |>
    mutate(tot_minus_target = great_tot - .data$nIS,
           totIS = great_tot)
  is_distinct_per_gene
}

fisher_per_gene <- function(df1, df2, threshold = 0.05) {
  p <- progressr::progressor(steps = length(union(df1$GeneName, df2$GeneName)))
  tmp <- purrr::map(union(df1$GeneName, df2$GeneName), ~ {
    gene <- .x
    gene_is_1 <- df1 |>
      filter(GeneName == gene) |>
      pull(nIS)
    if (purrr::is_empty(gene_is_1)) {
      return(NULL) # don't count genes that are not in common
    }
    gene_is_2 <- df2 |>
      filter(GeneName == gene) |>
      pull(nIS)
    if (purrr::is_empty(gene_is_2)) {
      return(NULL) # don't count genes that are not in common
    }
    other_is_1 <- df1 |>
      filter(GeneName == gene) |>
      pull(tot_minus_target)
    if (purrr::is_empty(other_is_1)) {
      other_is_1 <- df1$totIS[1]
    }
    other_is_2 <- df2 |>
      filter(GeneName == gene) |>
      pull(tot_minus_target)
    if (purrr::is_empty(other_is_2)) {
      other_is_2 <- df2$totIS[1]
    }
    contingency <- data.frame(
      "g1" = c(gene_is_1, other_is_1),
      "g2" = c(gene_is_2, other_is_2),
      stringsAsFactors = FALSE
    )
    colnames(contingency) <-
      c("in gene", "not in gene")
    test_res <-
      fisher.test(contingency, conf.level = 1 - threshold)
    p()
    return(
      tibble::tibble(
        GeneName = gene,
        nIS_1 = gene_is_1,
        n_IS_2 = gene_is_2,
        other_IS_1 = other_is_1,
        other_IS_2 = other_is_2,
        p_value = test_res$p.value
      )
    )
  }) |> 
    purrr::list_rbind() 
  
  if (!purrr::is_empty(tmp)) {
    tmp <- tmp |>
      mutate(
        p_value_fdr = p.adjust(p_value, method = "fdr", 
                               n = length(p_value)),
        p_value_bonferroni = p.adjust(p_value, method = "bonferroni", 
                                      n = length(p_value)),
        p_value_benjamini = p.adjust(p_value, method = "BY", 
                                     n = length(p_value))
      )
    return(tmp)
  } else {
    return(NULL)
  }
}

# Intra-trial function ---------------------------------------------------------
# To compute fisher tests ON THE SAME TRIAL, between different lineages 
# (no info on patient and tissue is preserved)

intra_trial_comp <- function(trial = c("MLD", "BTHAL", "WAS")) {
  trial <- str_to_lower(rlang::arg_match(trial)) 
  
  comparison_table <- gtools::combinations(
    n = length(is_per_late_label[[trial]]$Late_StateLabel), 
    r = 2, v = is_per_late_label[[trial]]$Late_StateLabel, set = TRUE)
  colnames(comparison_table) <- c("g1", "g2")
  comparison_table <- tibble::as_tibble(comparison_table)
  
  ## NO SAMPLING ----
  print("Performing simple Fisher tests...")
  fisher_tests_simple <- purrr::pmap(comparison_table, function(g1, g2, ...) {
    df_g1 <- flag_data[[trial]] |> 
      filter(Late_StateLabel == g1) |>
      extract_data_test()
    df_g2 <- flag_data[[trial]] |> 
             filter(Late_StateLabel == g2) |>
             extract_data_test()
    fisher_per_gene(df_g1, df_g2)
  }) |> purrr::set_names(
    comparison_table |> tidyr::unite(col = "id", g1, g2) |> pull("id")
  )
  
  ## WITH SAMPLING ----
  print("Performing Fisher tests with sampling...")
  sampling_tests <- function(g1, g2, ...) {
    results_df <- NULL
    for (i in 1:REPS) {
      print(paste("Rep", i, "of", REPS))
      df_g1 <- flag_data[[trial]] |> 
        filter(Late_StateLabel == g1) |>
        slice_sample(n = SAMPLE_SIZE, replace = FALSE) |>
        extract_data_test()
      df_g2 <- flag_data[[trial]] |> 
        filter(Late_StateLabel == g2) |>
        slice_sample(n = SAMPLE_SIZE, replace = FALSE) |>
        extract_data_test()
      res <- fisher_per_gene(df_g1, df_g2)
      if (!is.null(res)) {
        if (is.null(results_df)) {
          results_df <- res |>
            mutate(rep_n = i)
        } else {
          results_df <- results_df |>
            bind_rows(
              res |>
                mutate(rep_n = i)
            )
        }
      }
    }
    return(results_df)
  }
  fisher_sampling <- purrr::pmap(comparison_table, sampling_tests) |> 
    purrr::set_names(
    comparison_table |> tidyr::unite(col = "id", g1, g2) |> pull("id")
  )
  
  return(list(simple = fisher_tests_simple, sampling = fisher_sampling))
}

mld_intra_trial_res <- intra_trial_comp("MLD")
mld_intra_trial_simple <- purrr::map2(mld_intra_trial_res$simple,
                                      names(mld_intra_trial_res$simple),
                                        ~ {
                                          name_sep <- str_split(.y, "_")
                                          .x |>
                                            mutate(
                                              Group1 = name_sep[[1]][1],
                                              Group2 = name_sep[[1]][2]
                                            )
                                        }) |> purrr::list_rbind()

mld_intra_trial_simple |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("MLD.", 
                           "20230905_INTRA-trial-comparison-lineages-committed.tsv")),
    na = ""
  )

mld_intra_trial_sampling <- purrr::map2(mld_intra_trial_res$sampling,
                                        names(mld_intra_trial_res$sampling),
                                       ~ {
                                         groups <- str_split(
                                           .y, "_"
                                         )
                                         .x |> mutate(
                                           Group1 = groups[[1]][1],
                                           Group2 = groups[[1]][1]
                                         )
                                       }) |> 
  purrr::list_rbind()

mld_intra_p_values_stats <- mld_intra_trial_sampling |>
  group_by(GeneName, Group1, Group2) |>
  summarise(
    p_value_min = min(p_value, na.rm = TRUE),
    p_value_max = max(p_value, na.rm = TRUE),
    p_value_mean = mean(p_value, na.rm = TRUE),
    p_value_sd = sd(p_value, na.rm = TRUE),
    p_value_count = n(),
    .groups = "drop"
  )

mld_intra_p_values_stats |>
  arrange(p_value_mean) |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("MLD.", 
                           "20230906_INTRA-trial-comparison-lineages-committed_SAMPLING-PVALUE-STATS.tsv")),
    na = ""
  )

mld_intra_trial_sampling |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("MLD.", 
                           "20230906_INTRA-trial-comparison-lineages-committed_SAMPLING-REP-DETAILS.tsv")),
    na = ""
  )

bthal_intra_trial_res <- intra_trial_comp("BTHAL")
bthal_intra_trial_simple <- purrr::map2(bthal_intra_trial_res$simple,
                                      names(bthal_intra_trial_res$simple),
                                      ~ {
                                        name_sep <- str_split(.y, "_")
                                        .x |>
                                          mutate(
                                            Group1 = name_sep[[1]][1],
                                            Group2 = name_sep[[1]][2]
                                          )
                                      }) |> purrr::list_rbind()
bthal_intra_trial_simple |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("BTHAL.", 
                           "20230905_INTRA-trial-comparison-lineages-committed.tsv")),
    na = ""
  )

bthal_intra_trial_sampling <- purrr::map2(bthal_intra_trial_res$sampling,
                                        names(bthal_intra_trial_res$sampling),
                                        ~ {
                                          groups <- str_split(
                                            .y, "_"
                                          )
                                          .x |> mutate(
                                            Group1 = groups[[1]][1],
                                            Group2 = groups[[1]][1]
                                          )
                                        }) |> 
  purrr::list_rbind()

bthal_intra_p_values_stats <- bthal_intra_trial_sampling |>
  group_by(GeneName, Group1, Group2) |>
  summarise(
    p_value_min = min(p_value, na.rm = TRUE),
    p_value_max = max(p_value, na.rm = TRUE),
    p_value_mean = mean(p_value, na.rm = TRUE),
    p_value_sd = sd(p_value, na.rm = TRUE),
    p_value_count = n(),
    .groups = "drop"
  )

bthal_intra_p_values_stats |>
  arrange(p_value_mean) |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("BTHAL.", 
                           "20230906_INTRA-trial-comparison-lineages-committed_SAMPLING-PVALUE-STATS.tsv")),
    na = ""
  )

bthal_intra_trial_sampling |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("BTHAL.", 
                           "20230906_INTRA-trial-comparison-lineages-committed_SAMPLING-REP-DETAILS.tsv")),
    na = ""
  )


was_intra_trial_res <- intra_trial_comp("WAS")
was_intra_trial_simple <- purrr::map2(was_intra_trial_res$simple,
                                        names(was_intra_trial_res$simple),
                                        ~ {
                                          name_sep <- str_split(.y, "_")
                                          .x |>
                                            mutate(
                                              Group1 = name_sep[[1]][1],
                                              Group2 = name_sep[[1]][2]
                                            )
                                        }) |> purrr::list_rbind()
was_intra_trial_simple |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("WAS.", 
                           "20230905_INTRA-trial-comparison-lineages-committed.tsv")),
    na = ""
  )

was_intra_trial_sampling <- purrr::map2(was_intra_trial_res$sampling,
                                          names(was_intra_trial_res$sampling),
                                          ~ {
                                            groups <- str_split(
                                              .y, "_"
                                            )
                                            .x |> mutate(
                                              Group1 = groups[[1]][1],
                                              Group2 = groups[[1]][1]
                                            )
                                          }) |> 
  purrr::list_rbind()

was_intra_p_values_stats <- was_intra_trial_sampling |>
  group_by(GeneName, Group1, Group2) |>
  summarise(
    p_value_min = min(p_value, na.rm = TRUE),
    p_value_max = max(p_value, na.rm = TRUE),
    p_value_mean = mean(p_value, na.rm = TRUE),
    p_value_sd = sd(p_value, na.rm = TRUE),
    p_value_count = n(),
    .groups = "drop"
  )

was_intra_p_values_stats |>
  arrange(p_value_mean) |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("WAS.", 
                           "20230906_INTRA-trial-comparison-lineages-committed_SAMPLING-PVALUE-STATS.tsv")),
    na = ""
  )

was_intra_trial_sampling |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("WAS.", 
                           "20230906_INTRA-trial-comparison-lineages-committed_SAMPLING-REP-DETAILS.tsv")),
    na = ""
  )

# Inter-trial function ---------------------------------------------------------
# To compute fisher tests BETWEEN DIFFERENT TRIALS, between different lineages 
# (no info on patient and tissue is preserved)

inter_trial_comp <- function(trial1 = c("MLD", "BTHAL", "WAS"),
                             trial2 = c("MLD", "BTHAL", "WAS")) {
  trial1 <- str_to_lower(rlang::arg_match(trial1)) 
  trial2 <- str_to_lower(rlang::arg_match(trial2))
  
  comparison_table <- tidyr::expand_grid(
    g1 = paste0(str_to_upper(trial1), "-", is_per_late_label[[trial1]]$Late_StateLabel),
    g2 = paste0(str_to_upper(trial2), "-", is_per_late_label[[trial2]]$Late_StateLabel)
  )
  
  ## NO SAMPLING ----
  print("Performing simple Fisher tests...")
  fisher_tests_simple <- purrr::pmap(comparison_table, function(g1, g2, ...) {
    df_g1 <- flag_data[[trial1]] |> 
      filter(Late_StateLabel == str_remove(g1, paste0(str_to_upper(trial1), "-"))) |>
      extract_data_test()
    df_g2 <- flag_data[[trial2]] |> 
      filter(Late_StateLabel == str_remove(g2, paste0(str_to_upper(trial2), "-"))) |>
      extract_data_test()
    fisher_per_gene(df_g1, df_g2)
  }) |> purrr::set_names(
    comparison_table |> tidyr::unite(col = "id", g1, g2) |> pull("id")
  )
  
  ## WITH SAMPLING ----
  print("Performing Fisher tests with sampling...")
  sampling_tests <- function(g1, g2, ...) {
    results_df <- NULL
    for (i in 1:REPS) {
      print(paste("Rep", i, "of", REPS))
      df_g1 <- flag_data[[trial1]] |> 
        filter(Late_StateLabel == str_remove(g1, paste0(str_to_upper(trial1), "-"))) |>
        slice_sample(n = SAMPLE_SIZE, replace = FALSE) |>
        extract_data_test()
      df_g2 <- flag_data[[trial2]] |> 
        filter(Late_StateLabel == str_remove(g2, paste0(str_to_upper(trial2), "-"))) |>
        slice_sample(n = SAMPLE_SIZE, replace = FALSE) |>
        extract_data_test()
      res <- fisher_per_gene(df_g1, df_g2)
      if (!is.null(res)) {
        if (is.null(results_df)) {
          results_df <- res |>
            mutate(rep_n = i)
        } else {
          results_df <- results_df |>
            bind_rows(
              res |>
                mutate(rep_n = i)
            )
        }
      }
    }
    return(results_df)
  }
  fisher_sampling <- purrr::pmap(comparison_table, sampling_tests) |> 
    purrr::set_names(
      comparison_table |> tidyr::unite(col = "id", g1, g2) |> pull("id")
    )
  
  return(list(simple = fisher_tests_simple, sampling = fisher_sampling))
}

mld_bthal_inter_trial_res <- inter_trial_comp(trial1 = "MLD", trial2 = "BTHAL")
mld_bthal_inter_trial_simple <- purrr::map2(mld_bthal_inter_trial_res$simple,
                                        names(mld_bthal_inter_trial_res$simple),
                                        ~ {
                                          name_sep <- str_split(.y, "_")
                                          .x |>
                                            mutate(
                                              Group1 = name_sep[[1]][1],
                                              Group2 = name_sep[[1]][2]
                                            )
                                        }) |> purrr::list_rbind()
mld_bthal_inter_trial_simple |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("MLD-BTHAL.", 
                           "20230906_INTER-trial-comparison-lineages-committed.tsv")),
    na = ""
  )

mld_bthal_inter_trial_sampling <- purrr::map2(mld_bthal_inter_trial_res$sampling,
                                        names(mld_bthal_inter_trial_res$sampling),
                                        ~ {
                                          groups <- str_split(
                                            .y, "_"
                                          )
                                          if (!is.null(.x)) {
                                            return(.x |> mutate(
                                              Group1 = groups[[1]][1],
                                              Group2 = groups[[1]][1]
                                            ))
                                          }
                                        }) |> 
  purrr::list_rbind()

mld_bthal_inter_p_values_stats <- mld_bthal_inter_trial_sampling |>
  group_by(GeneName, Group1, Group2) |>
  summarise(
    p_value_min = min(p_value, na.rm = TRUE),
    p_value_max = max(p_value, na.rm = TRUE),
    p_value_mean = mean(p_value, na.rm = TRUE),
    p_value_sd = sd(p_value, na.rm = TRUE),
    p_value_count = n(),
    .groups = "drop"
  )

mld_bthal_inter_p_values_stats |>
  arrange(p_value_mean) |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("MLD-BTHAL.", 
                           "20230906_INTER-trial-comparison-lineages-committed_SAMPLING-PVALUE-STATS.tsv")),
    na = ""
  )

mld_bthal_inter_trial_sampling |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("MLD-BTHAL.", 
                           "20230906_INTER-trial-comparison-lineages-committed_SAMPLING-REP-DETAILS.tsv")),
    na = ""
  )

mld_was_inter_trial_res <- inter_trial_comp(trial1 = "MLD", trial2 = "WAS")
mld_was_inter_trial_simple <- purrr::map2(mld_was_inter_trial_res$simple,
                                            names(mld_was_inter_trial_res$simple),
                                            ~ {
                                              name_sep <- str_split(.y, "_")
                                              .x |>
                                                mutate(
                                                  Group1 = name_sep[[1]][1],
                                                  Group2 = name_sep[[1]][2]
                                                )
                                            }) |> purrr::list_rbind()
mld_was_inter_trial_simple |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("MLD-WAS.", 
                           "20230906_INTER-trial-comparison-lineages-committed.tsv")),
    na = ""
  )

mld_was_inter_trial_sampling <- purrr::map2(mld_was_inter_trial_res$sampling,
                                              names(mld_was_inter_trial_res$sampling),
                                              ~ {
                                                groups <- str_split(
                                                  .y, "_"
                                                )
                                                if (!is.null(.x)) {
                                                  return(.x |> mutate(
                                                    Group1 = groups[[1]][1],
                                                    Group2 = groups[[1]][1]
                                                  ))
                                                }
                                              }) |> 
  purrr::list_rbind()

mld_was_inter_p_values_stats <- mld_was_inter_trial_sampling |>
  group_by(GeneName, Group1, Group2) |>
  summarise(
    p_value_min = min(p_value, na.rm = TRUE),
    p_value_max = max(p_value, na.rm = TRUE),
    p_value_mean = mean(p_value, na.rm = TRUE),
    p_value_sd = sd(p_value, na.rm = TRUE),
    p_value_count = n(),
    .groups = "drop"
  )

mld_was_inter_p_values_stats |>
  arrange(p_value_mean) |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("MLD-WAS.", 
                           "20230906_INTER-trial-comparison-lineages-committed_SAMPLING-PVALUE-STATS.tsv")),
    na = ""
  )

mld_was_inter_trial_sampling |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("MLD-WAS.", 
                           "20230906_INTER-trial-comparison-lineages-committed_SAMPLING-REP-DETAILS.tsv")),
    na = ""
  )


bthal_was_inter_trial_res <- inter_trial_comp(trial1 = "BTHAL", trial2 = "WAS")
bthal_was_inter_trial_simple <- purrr::map2(bthal_was_inter_trial_res$simple,
                                          names(bthal_was_inter_trial_res$simple),
                                          ~ {
                                            name_sep <- str_split(.y, "_")
                                            .x |>
                                              mutate(
                                                Group1 = name_sep[[1]][1],
                                                Group2 = name_sep[[1]][2]
                                              )
                                          }) |> purrr::list_rbind()
bthal_was_inter_trial_simple |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("BTHAL-WAS.", 
                           "20230906_INTER-trial-comparison-lineages-committed.tsv")),
    na = ""
  )

bthal_was_inter_trial_sampling <- purrr::map2(bthal_was_inter_trial_res$sampling,
                                            names(bthal_was_inter_trial_res$sampling),
                                            ~ {
                                              groups <- str_split(
                                                .y, "_"
                                              )
                                              if (!is.null(.x)) {
                                                return(.x |> mutate(
                                                  Group1 = groups[[1]][1],
                                                  Group2 = groups[[1]][1]
                                                ))
                                              }
                                            }) |> 
  purrr::list_rbind()

bthal_was_inter_p_values_stats <- bthal_was_inter_trial_sampling |>
  group_by(GeneName, Group1, Group2) |>
  summarise(
    p_value_min = min(p_value, na.rm = TRUE),
    p_value_max = max(p_value, na.rm = TRUE),
    p_value_mean = mean(p_value, na.rm = TRUE),
    p_value_sd = sd(p_value, na.rm = TRUE),
    p_value_count = n(),
    .groups = "drop"
  )

bthal_was_inter_p_values_stats |>
  arrange(p_value_mean) |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("BTHAL-WAS.", 
                           "20230906_INTER-trial-comparison-lineages-committed_SAMPLING-PVALUE-STATS.tsv")),
    na = ""
  )

bthal_was_inter_trial_sampling |>
  readr::write_tsv(
    file = fs::path(comm_fold, 
                    "comparisons",
                    paste0("BTHAL-WAS.", 
                           "20230906_INTER-trial-comparison-lineages-committed_SAMPLING-REP-DETAILS.tsv")),
    na = ""
  )

# Sharing early/late with no filtering -----------------------------------------
## MLD ----
mld_data_folder <- fs::path("/mnt/CerberoWorkspaceISA/MLD/archive/202112")
mld_matrix_regex <- paste0(
  mld_data_folder,
  "/matrices",
  "/aggregated",
  "/MLD.+\\..+SubjectID_CellMarker_Tissue_TimepointMonths",
  "_fragmentEstimate"
)
mld_matrices_paths <- fs::dir_ls(fs::path(mld_data_folder, "matrices",
                                          "aggregated"), regexp = mld_matrix_regex)
### Import IS matrices ----
mld_data <- purrr::map(mld_matrices_paths, ~ {
  matrix <- readr::read_tsv(
    .x,
    col_types = readr::cols(
      .default = "d",
      chr = "c",
      strand = "c",
      GeneName = "c",
      GeneStrand = "c",
      integration_locus = "i"
    )
  )
  matrix |> tidyr::pivot_longer(
    cols = !c(mandatory_IS_vars(), annotation_IS_vars()),
    names_to = c("SubjectID", "CellMarker", "Tissue", "TimepointMonths"),
    names_sep = "_",
    values_to = "fragmentEstimate",
    values_drop_na = TRUE
  )
}) |> purrr::list_rbind()

bthal_data_folder <-
  fs::path("/mnt/CerberoWorkspaceISA/BTHAL/archive/202112")
bthal_matrix_regex <- paste0(
  "BTHAL[0-9]{3}.+_SubjectID_CellMarker_Tissue_",
  "TimepointMonths_fragmentEstimate"
)
bthal_matrices_paths <- fs::dir_ls(
  fs::path(bthal_data_folder, "matrices", "aggregated"), 
  regexp = bthal_matrix_regex)

### Import IS matrices ----
bthal_data <- purrr::map(bthal_matrices_paths, ~ {
  matrix <- readr::read_tsv(
    .x,
    col_types = readr::cols(
      .default = "d",
      chr = "c",
      strand = "c",
      GeneName = "c",
      GeneStrand = "c",
      integration_locus = "i"
    )
  )
  matrix |> 
    tidyr::pivot_longer(
      cols = !c(mandatory_IS_vars(), annotation_IS_vars()),
      names_to = c("SubjectID", "CellMarker", "Tissue", "TimepointMonths"),
      names_sep = "_",
      values_to = "fragmentEstimate",
      values_drop_na = TRUE
    )
}) |> purrr::list_rbind()


was_data_folder <-
  fs::path("/mnt/CerberoWorkspaceISA/WAS/archive/202112")
was_matrix_regex <- paste0(
  was_data_folder,
  "/matrices",
  "/aggregated",
  "/WAS.+\\..+SubjectID_CellMarker_Tissue_TimepointMonths_fragmentEstimate"
)
was_matrices_paths <-
  fs::dir_ls(fs::path(was_data_folder, "matrices",
                      "aggregated"),
             regexp = was_matrix_regex)

### Import IS matrices ----
was_data <- purrr::map(was_matrices_paths, ~ {
  matrix <- readr::read_tsv(
    .x,
    col_types = readr::cols(
      .default = "d",
      chr = "c",
      strand = "c",
      GeneName = "c",
      GeneStrand = "c",
      integration_locus = "i"
    )
  )
  matrix |> 
    tidyr::pivot_longer(
      cols = !c(mandatory_IS_vars(), annotation_IS_vars()),
      names_to = c("SubjectID", "CellMarker", "Tissue", "TimepointMonths"),
      names_sep = "_",
      values_to = "fragmentEstimate",
      values_drop_na = TRUE
    )
}) |> purrr::list_rbind()

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

datasets <- list(
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

sharing_results <- purrr::map(datasets, ~ {
  df <- .x
  patients <- unique(.x$SubjectID)
  purrr::map(patients, ~ sharing_same_p(df, .x)) |>
    purrr::list_rbind()
})

purrr::iwalk(sharing_results, ~ {
  file_name <- paste0(.y, ".", "202300906_sharing-TP24_NO-FILTERS.tsv")
  .x |>
    select(-is_coord, -truth_tbl_venn) |>
    readr::write_tsv(
      file = fs::path(sharing_folder, file_name),
      na = ""
    )
})

save.image(file = "./HSPC_dynamics_rev1_committed_analyses.RData")
