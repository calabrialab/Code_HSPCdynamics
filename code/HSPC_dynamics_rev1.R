################################################################################
# HSPC Dynamics additional analyses
################################################################################
# 2023-06-06

library(ISAnalytics)
library(dplyr)

progressr::handlers(global = TRUE)
progressr::handlers("cli")

proj_folder <- fs::path_wd()
track_comp_folder <- fs::path(proj_folder, "track_comparison")

for (fold in c(track_comp_folder)) {
  fs::dir_create(fold)
  
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

disease_genes <- c("ARSA", "WAS", "HBB", "HBE1")

# Functions --------------------------------------------------------------------
extract_data_test <- function(dataset) {
  is_distinct_per_gene <- dataset |>
    group_by(GeneName) |>
    summarise(nIS = n_distinct(chr, integration_locus, strand))
  great_tot <- sum(is_distinct_per_gene$nIS)
  is_distinct_per_gene <- is_distinct_per_gene |>
    mutate(tot_minus_target = great_tot - .data$nIS,
           totIS = great_tot)
}

fisher_per_gene <- function(df1, df2, threshold = 0.05) {
  p <- progressr::progressor(steps = length(union(df1$GeneName, df2$GeneName)))
  purrr::map(union(df1$GeneName, df2$GeneName), ~ {
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
    purrr::list_rbind() |>
    mutate(
      p_value_fdr = p.adjust(p_value, method = "fdr", 
                             n = length(p_value)),
      p_value_bonferroni = p.adjust(p_value, method = "bonferroni", 
                                    n = length(p_value)),
      p_value_benjamini = p.adjust(p_value, method = "BY", 
                                   n = length(p_value))
    )
}


# BTHAL ------------------------------------------------------------------------
bthal_data_folder <-
  fs::path("/mnt/CerberoWorkspaceISA/BTHAL/archive/202112")
# rstudioapi::navigateToFile(file = fs::path(
#   bthal_data_folder, "BTHAL_analyses_script_freeze-202112.R"))
bthal_matrix_regex <- paste0(
    "BTHAL[0-9]{3}.+_SubjectID_CellMarker_Tissue_",
    "TimepointMonths_SC-filtered_seqCount"
  )
bthal_matrices_paths <- fs::dir_ls(fs::path(bthal_data_folder, "matrices",
                      "aggregated"), regexp = bthal_matrix_regex)

## Import matrices ----
bthal_data <- purrr::map(bthal_matrices_paths, ~ {
  matrix <- readr::read_tsv(
    .x,
    col_types = readr::cols(
      .default = "i",
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
    values_to = "seqCount",
    values_drop_na = TRUE
  )
}) |> purrr::list_rbind()

## Only in-vivo data ----
bthal_data_invivo <- bthal_data |>
  filter(as.numeric(TimepointMonths) > 0,
         CellMarker %in% included_markers,
         !GeneName %in% disease_genes)


## Extract data for Fisher test ----
bthal_fisher_data <- extract_data_test(bthal_data_invivo)

p_values_stats <- function(ds, gene) {
  gene_ds <- ds |>
    filter(GeneName == gene)
  result <- tibble::tibble(
    GeneName = gene,
    p_value_fdr_min = min(gene_ds$Fisher_p_value_fdr, na.rm = TRUE),
    p_value_fdr_max = max(gene_ds$Fisher_p_value_fdr, na.rm = TRUE),
    p_value_fdr_mean = mean(gene_ds$Fisher_p_value_fdr, na.rm = TRUE),
    p_value_bonferroni_min = min(gene_ds$Fisher_p_value_bonferroni, na.rm = TRUE),
    p_value_bonferroni_max = max(gene_ds$Fisher_p_value_bonferroni, na.rm = TRUE),
    p_value_bonferroni_mean = mean(gene_ds$Fisher_p_value_bonferroni, na.rm = TRUE),
    p_value_benjamini_min = min(gene_ds$Fisher_p_value_benjamini, na.rm = TRUE),
    p_value_benjamini_max = max(gene_ds$Fisher_p_value_benjamini, na.rm = TRUE),
    p_value_benjamini_mean = mean(gene_ds$Fisher_p_value_benjamini, na.rm = TRUE)
  )
  return(result)
}


# WAS --------------------------------------------------------------------------
was_data_folder <-
  fs::path("/mnt/CerberoWorkspaceISA/WAS/archive/202112")
# rstudioapi::navigateToFile(file = fs::path(
#   was_data_folder, "WAS_analyses_script_freeze-202112.R"))
was_matrix_regex <- paste0(
  was_data_folder,
  "/matrices",
  "/aggregated",
  "/WAS.+\\..+SubjectID_CellMarker_Tissue_TimepointMonths_seqCount"
)
was_matrices_paths <-
  fs::dir_ls(fs::path(was_data_folder, "matrices",
                      "aggregated"),
             regexp = was_matrix_regex)
## Import matrices ----
was_data <- purrr::map(was_matrices_paths, ~ {
  matrix <- readr::read_tsv(
    .x,
    col_types = readr::cols(
      .default = "i",
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
    values_to = "seqCount",
    values_drop_na = TRUE
  )
}) |> purrr::list_rbind()

## Only in-vivo data ----
was_data_invivo <- was_data |>
  filter(
    as.numeric(TimepointMonths) > 0,
    CellMarker %in% included_markers,
    Tissue %in% c("BM", "PB"),
    !GeneName %in% disease_genes
  )

## Extract data for Fisher test ----
was_fisher_data <- extract_data_test(was_data_invivo)


# MLD --------------------------------------------------------------------------
mld_data_folder <-
  fs::path("/mnt/CerberoWorkspaceISA/MLD/archive/202112")
# rstudioapi::navigateToFile(file = fs::path(
#   mld_data_folder, "MLD_analyses_script_freeze-202112.R"))
mld_matrix_regex <- paste0(
  mld_data_folder,
  "/matrices",
  "/aggregated",
  "/MLD.+\\..+SubjectID_CellMarker_Tissue_TimepointMonths_SC-filtered_seqCount"
)
mld_matrices_paths <-
  fs::dir_ls(fs::path(mld_data_folder, "matrices",
                      "aggregated"),
             regexp = mld_matrix_regex)
## Import matrices ----
mld_data <- purrr::map(mld_matrices_paths, ~ {
  matrix <- readr::read_tsv(
    .x,
    col_types = readr::cols(
      .default = "i",
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
    values_to = "seqCount",
    values_drop_na = TRUE
  )
}) |> purrr::list_rbind()

## Only in-vivo data ----
mld_data_invivo <- mld_data |>
  filter(as.numeric(TimepointMonths) > 0,
         CellMarker %in% included_markers,
         !GeneName %in% disease_genes)

## Extract data for Fisher test ----
mld_fisher_data <- extract_data_test(mld_data_invivo)

# Per-gene Fisher test (all patients) ------------------------------------------
## BTHAL vs. WAS ----
fisher_bthal_was <- fisher_per_gene(df1 = bthal_fisher_data,
                                    df2 = was_fisher_data)
fisher_bthal_was |>
  arrange(p_value_fdr, p_value_bonferroni) |>
  readr::write_tsv(
    fs::path(track_comp_folder, 
             "20230606_Fisher-gene-targeting_BTHAL-vs-WAS.tsv")
  )

## BTHAL vs. MLD ----
fisher_bthal_mld <- fisher_per_gene(df1 = bthal_fisher_data,
                                    df2 = mld_fisher_data)
fisher_bthal_mld |>
  arrange(p_value_fdr, p_value_bonferroni) |>
  readr::write_tsv(
    fs::path(track_comp_folder, 
             "20230606_Fisher-gene-targeting_BTHAL-vs-MLD.tsv")
  )

## WAS vs. MLD ----
fisher_was_mld <- fisher_per_gene(df1 = was_fisher_data,
                                  df2 = mld_fisher_data)

fisher_was_mld |>
  arrange(p_value_fdr, p_value_bonferroni) |>
  readr::write_tsv(
    fs::path(track_comp_folder, 
             "20230606_Fisher-gene-targeting_WAS-vs-MLD.tsv")
  )

# Sampling method --------------------------------------------------------------
sampling_dir <- fs::path(proj_folder, "sampling_method")
fs::dir_create(sampling_dir)

set.seed(436)
sample_size <- 100000
n_trials <- 10
sample_and_test <- function(n_trial, data_list, threshold = 0.05) {
  print(sprintf("Processing trial %i...", n_trial))
  # Get samples from datasets
  bthal_sample <- data_list$bthal |>
    slice_sample(n = sample_size, replace = FALSE)
  was_sample <- data_list$bthal |>
    slice_sample(n = sample_size, replace = FALSE)
  mld_sample <- data_list$mld |>
    slice_sample(n = sample_size, replace = FALSE)
  # Extract data
  bthal_sample_data <- extract_data_test(bthal_sample)
  was_sample_data <- extract_data_test(was_sample)
  mld_sample_data <- extract_data_test(mld_sample)
  # Perform tests
  ## BTHAL vs. WAS
  print("--- Fisher bthal vs. was")
  sample_fisher_bthal_was <- fisher_per_gene(
    df1 = bthal_sample_data,
    df2 = was_sample_data,
    threshold = threshold
  )
  ## BTHAL vs. MLD
  print("--- Fisher bthal vs. mld")
  sample_fisher_bthal_mld <- fisher_per_gene(
    df1 = bthal_sample_data,
    df2 = mld_sample_data,
    threshold = threshold
  )
  ## WAS vs. MLD
  print("--- Fisher was vs. mld")
  sample_fisher_was_mld <- fisher_per_gene(
    df1 = was_sample_data,
    df2 = mld_sample_data,
    threshold = threshold
  )
  return(list(
    bthal_was = sample_fisher_bthal_was,
    bthal_mld = sample_fisher_bthal_mld,
    was_mld = sample_fisher_was_mld
  ))
}

data_list <- list(
  bthal = bthal_data_invivo,
  was = was_data_invivo,
  mld = mld_data_invivo
)

sampling_results <- purrr::map(seq_len(n_trials), 
                               ~ sample_and_test(.x, data_list = data_list))

purrr::walk(seq_along(sampling_results), ~ {
  print(sprintf("Saving trial %i...", .x))
  file_name <- sprintf("0%i.20230606_sampling_trial-fisher-results.xlsx",
                       .x)
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(
    wb = wb,
    sheetName = "BTHAL vs. WAS"
  )
  openxlsx::addWorksheet(
    wb = wb,
    sheetName = "BTHAL vs. MLD"
  )
  openxlsx::addWorksheet(
    wb = wb,
    sheetName = "WAS vs. MLD"
  )
  openxlsx::writeData(
    wb = wb,
    sheet = "BTHAL vs. WAS",
    x = sampling_results[[.x]]$bthal_was |>
      arrange(p_value_fdr, p_value_bonferroni)
  )
  openxlsx::writeData(
    wb = wb,
    sheet = "BTHAL vs. MLD",
    x = sampling_results[[.x]]$bthal_mld |>
      arrange(p_value_fdr, p_value_bonferroni)
  )
  openxlsx::writeData(
    wb = wb,
    sheet = "WAS vs. MLD",
    x = sampling_results[[.x]]$was_mld |>
      arrange(p_value_fdr, p_value_bonferroni)
  )
  openxlsx::saveWorkbook(
    wb = wb,
    file = fs::path(sampling_dir, file_name),
    overwrite = TRUE
  )
})

## Saving single big matrix
reshape_frames <- function(n_trial, bthal_was, bthal_mld, was_mld) {
  bthal_was_df <- bthal_was |>
    rename(nIS_bthal = "nIS_1",
           nIS_was = "n_IS_2",
           other_IS_bthal = "other_IS_1",
           other_IS_was = "other_IS_2") |>
    tidyr::pivot_longer(
      cols = starts_with("p_value"),
      names_to = "p_value_type",
      values_to = "p_value"
    ) |>
    mutate(
      comparison = "BTHAL_WAS"
    ) |>
    mutate(
      nIS_mld = NA_integer_,
      .after = nIS_was
    ) |>
    mutate(
      other_IS_mld = NA_integer_,
      .after = other_IS_was
    )
  
  bthal_mld_df <- bthal_mld |>
    rename(nIS_bthal = "nIS_1",
           nIS_mld = "n_IS_2",
           other_IS_bthal = "other_IS_1",
           other_IS_mld = "other_IS_2") |>
    tidyr::pivot_longer(
      cols = starts_with("p_value"),
      names_to = "p_value_type",
      values_to = "p_value"
    ) |>
    mutate(
      comparison = "BTHAL_MLD"
    ) |>
    mutate(
      nIS_was = NA_integer_,
      .after = nIS_bthal
    ) |>
    mutate(
      other_IS_was = NA_integer_,
      .after = other_IS_bthal
    )
  
  was_mld_df <- was_mld |>
    rename(nIS_was = "nIS_1",
           nIS_mld = "n_IS_2",
           other_IS_was = "other_IS_1",
           other_IS_mld = "other_IS_2") |>
    tidyr::pivot_longer(
      cols = starts_with("p_value"),
      names_to = "p_value_type",
      values_to = "p_value"
    ) |>
    mutate(
      comparison = "WAS_MLD"
    ) |>
    mutate(
      nIS_bthal = NA_integer_,
      .before = nIS_was
    ) |>
    mutate(
      other_IS_was = NA_integer_,
      .before = other_IS_mld
    )
  
  bthal_was_df |>
    full_join(bthal_mld_df) |>
    full_join(was_mld_df) |>
    mutate(
      n_trial = n_trial
    )
}

sampling_all_matrix <- purrr::map(seq_along(sampling_results), ~ {
  reshape_frames(
    n_trial = .x,
    bthal_was = sampling_results[[.x]]$bthal_was,
    bthal_mld = sampling_results[[.x]]$bthal_mld,
    was_mld = sampling_results[[.x]]$was_mld
  )
}) |> purrr::list_rbind()

sampling_all_matrix |>
  readr::write_tsv(
    fs::path(sampling_dir, 
             "20230606_sampling_ALL-TRIALS.tsv.gz"),
    na = ""
  )

sampling_all_matrix_benj <- sampling_all_matrix |>
  filter(p_value_type == "p_value_benjamini")

p_values_stats <- sampling_all_matrix_benj |>
  group_by(GeneName, comparison) |>
  summarise(
    p_value_min = min(p_value, na.rm = TRUE),
    p_value_max = max(p_value, na.rm = TRUE),
    p_value_mean = mean(p_value, na.rm = TRUE),
    p_value_sd = sd(p_value, na.rm = TRUE),
    p_value_count = n(),
    .groups = "drop"
  )

p_values_stats |>
  readr::write_tsv(
    fs::path(sampling_dir, 
             "20230705_sampling_P-VALUES-STATS.tsv.gz"),
    na = ""
  )

# Genome binning ---------------------------------------------------------------
library(QDNAseq.hg19)

bins_100 <- getBinAnnotations(binSize = 100, genome = "hg19")
hg19_100kbp_bins <- Biobase::pData(bins_100) |>
  tibble::as_tibble() |>
  dplyr::rename(seqnames = "chromosome") |>
  dplyr::mutate(seqnames = paste0("chr", seqnames)) |>
  plyranges::as_granges()

bthal_distinct_is <- bthal_data_invivo |>
  distinct(chr, integration_locus, strand) |>
  mutate(chr = paste0("chr", chr)) |>
  dplyr::rename(seqnames = "chr", start = "integration_locus") |>
  mutate(end = start) |>
  plyranges::as_granges()

was_distinct_is <- was_data_invivo |>
  distinct(chr, integration_locus, strand) |>
  mutate(chr = paste0("chr", chr)) |>
  dplyr::rename(seqnames = "chr", start = "integration_locus") |>
  mutate(end = start) |>
  plyranges::as_granges()

mld_distinct_is <- mld_data_invivo |>
  distinct(chr, integration_locus, strand) |>
  mutate(chr = paste0("chr", chr)) |>
  dplyr::rename(seqnames = "chr", start = "integration_locus") |>
  mutate(end = start) |>
  plyranges::as_granges()

overlaps <- hg19_100kbp_bins |>
  mutate(
    bthal_olap = plyranges::count_overlaps(hg19_100kbp_bins, 
                                           bthal_distinct_is,
                                           maxgap = -1L,
                                           minoverlap = 0L),
    was_olap = plyranges::count_overlaps(hg19_100kbp_bins, 
                                         was_distinct_is,
                                         maxgap = -1L,
                                         minoverlap = 0L),
    mld_olap = plyranges::count_overlaps(hg19_100kbp_bins, 
                                         mld_distinct_is,
                                         maxgap = -1L,
                                         minoverlap = 0L)
  )

overlaps_tbl <- tibble::as_tibble(overlaps) |>
  select(-bases, -gc, -mappability, -blacklist, -residual, -use) |>
  rowwise() |>
  mutate(
    totIS_bin = sum(bthal_olap, was_olap, mld_olap)
  ) |>
  ungroup()

## Apply Fisher tests rowwise ----
.compute_test <- function(g1, g2, tot, threshold = 0.05) {
  contingency_tbl <- data.frame(
    "g1" = c(g1, tot - g1),
    "g2" = c(g2, tot - g2),
    stringsAsFactors = FALSE
  )
  fisher.test(contingency_tbl, conf.level = 1 - threshold)$p.value
}

overlaps_tbl_tests <- overlaps_tbl |>
  rowwise() |>
  mutate(
    bthal_was.pvalue = .compute_test(bthal_olap, was_olap, totIS_bin),
    bthal_mld.pvalue = .compute_test(bthal_olap, mld_olap, totIS_bin),
    was_mld.pvalue = .compute_test(was_olap, mld_olap, totIS_bin)
  )
overlaps_tbl_tests <- overlaps_tbl_tests |>
  ungroup() |>
  mutate(
    across(
      ends_with("pvalue"),
      .fns = list(
        fdr = ~ p.adjust(.x, method = "fdr", n = length(.x)),
        bonferroni = ~ p.adjust(.x, method = "bonferroni", 
                                n = length(.x)),
        benjamini = ~ p.adjust(.x, method = "BY",
                                n = length(.x))
      ),
      .names = "{.col}_{.fn}"
    )
  ) |>
  arrange(across(ends_with("pvalue_fdr")))

overlaps_tbl_tests |>
  readr::write_tsv(
    fs::path(track_comp_folder, 
             "20230606_binning_100kbp_fisher.tsv")
  )







