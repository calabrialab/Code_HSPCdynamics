################################################################################
# HSPC Dynamics additional analyses - Abundance curves
################################################################################
# 2023-06-07
library(ISAnalytics)
library(dplyr)
library(ggplot2)
library(ggridges)

progressr::handlers(global = TRUE)
progressr::handlers("cli")

proj_folder <- fs::path_wd()
ab_curves_fold <- fs::path(proj_folder, "abundance_curves")
for (fold in c(ab_curves_fold)) {
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

# Functions --------------------------------------------------------------------
get_plots <- function(data, project, x_scale = c("def", "log"),
                      stat_curves = c("density_ridges", "density")) {
  x_scale <- rlang::arg_match(x_scale)
  stat_curves <- rlang::arg_match(stat_curves)
  patients <- data |> distinct(SubjectID) |> pull(SubjectID)
  purrr::map(patients, ~ {
    if (stat_curves == "density_ridges") {
      plt <- ggplot(data |>
                      filter(SubjectID == .x),
                    aes(x = AB_mean, y = FU, fill = FU)) +
        geom_density_ridges(alpha = 0.8, size = 0.3) +
        facet_grid(Tissue ~ StateChangeLabel) +
        theme_ridges() +
        theme(legend.position = "none")
      if (x_scale == "log") {
        plt <- plt +
          labs(y = "", x = "Mean % abundance (log10)", 
               title = paste0(project, " - Patient: ", .x)) +
          scale_x_log10(labels = scales::scientific)
      } else {
        plt <- plt +
          labs(y = "", x = "Mean % abundance", 
               title = paste0(project, " - Patient: ", .x))
      }
      return(plt)
    } else {
      plt <- ggplot(data |>
                      filter(SubjectID == .x),
                    aes(x = AB_mean, y = FU, fill = FU,
                        height = after_stat(density))) +
        geom_density_ridges(alpha = 0.8, size = 0.3,
                            stat = "density") +
        facet_grid(Tissue ~ StateChangeLabel) +
        theme_ridges() +
        theme(legend.position = "none")
      if (x_scale == "log") {
        plt <- plt +
          labs(y = "", x = "Mean % abundance (log10)", 
               title = paste0(project, " - Patient: ", .x)) +
          scale_x_log10(labels = scales::scientific)
      } else {
        plt <- plt +
          labs(y = "", x = "Mean % abundance", 
               title = paste0(project, " - Patient: ", .x))
      }
      return(plt)
    }
  }) |> purrr::set_names(patients)
}

save_plots <- function(plot_list, folder, width, height, 
                       type = c("def", "log"), 
                       plot_type = c("abund", "deltas")) {
  type <- rlang::arg_match(type)
  plot_type <- rlang::arg_match(plot_type)
  if (plot_type == "abund") {
    purrr::walk2(plot_list, names(plot_list), ~ {
      file_name <- if (type == "def") {
        paste0(.y, ".abundance_dist_plot")
      } else {
        paste0(.y, ".abundance_dist_plot_log-10-x")
      }
      try(ggsave(
        plot = .x,
        filename = paste0(file_name, ".pdf"),
        path = folder,
        width = width,
        height = height
      ))
      try(ggsave(
        plot = .x,
        filename = paste0(file_name, ".png"),
        path = folder,
        width = width,
        height = height
      ))
    })
    return()
  }
  purrr::walk2(plot_list, names(plot_list), ~ {
    file_names <- if (type == "def") {
      c(paste0(.y, ".delta_dist_plot"), paste0(.y, ".fold-change_dist_plot"))
    } else {
      c(paste0(.y, ".delta_dist_plot_log-10-x"),
        paste0(.y, ".fold-change_dist_plot_log-10-x"))
    }
    try(ggsave(
      plot = .x$delta,
      filename = paste0(file_names[1], ".pdf"),
      path = folder,
      width = width,
      height = height
    ))
    try(ggsave(
      plot = .x$delta,
      filename = paste0(file_names[1], ".png"),
      path = folder,
      width = width,
      height = height
    ))
    try(ggsave(
      plot = .x$fc,
      filename = paste0(file_names[2], ".pdf"),
      path = folder,
      width = width,
      height = height
    ))
    try(ggsave(
      plot = .x$fc,
      filename = paste0(file_names[2], ".png"),
      path = folder,
      width = width,
      height = height
    ))
  })
}

plot_delta_fc <- function(data, project, x_scale = c("def", "log")) {
  x_scale <- rlang::arg_match(x_scale)
  patients <- data |> distinct(SubjectID) |> pull(SubjectID)
  purrr::map(patients, ~ {
    plt_delta <- ggplot(data |>
                          filter(SubjectID == .x),
                        aes(x = delta)
    ) +
      geom_density() +
      theme_bw() +
      labs(x = "Delta % mean abundance (late - early)",
           y = "",
           title = paste0(project, " - Patient: ", .x)) +
      facet_grid(Tissue ~ StateChangeLabel)
    if (x_scale == "log") {
      plt_delta <- plt_delta +
        labs(x = "Log10 Delta % mean abundance (late - early)") +
        scale_x_log10(labels = scales::scientific)
    }
    plt_fc <- ggplot(data |>
                       filter(SubjectID == .x),
                     aes(x = fold_change)
    ) +
      geom_density() +
      theme_bw() +
      labs(x = "Fold change % mean abundance (late/early)",
           y = "",
           title = paste0(project, " - Patient: ", .x)) +
      facet_grid(Tissue ~ StateChangeLabel)
    if (x_scale == "log") {
      plt_fc <- plt_fc +
        labs(x = "Log10 fold change % mean abundance (late/early)") +
        scale_x_log10(labels = scales::scientific)
    }
    return(list(delta = plt_delta, fc = plt_fc))
  }) |> purrr::set_names(patients)
}



# Abundance distributions ------------------------------------------------------
## MLD ----
mld_data_folder <- fs::path("/mnt/CerberoWorkspaceISA/MLD/archive/202112")
mld_matrix_regex <- paste0(
  mld_data_folder,
  "/matrices",
  "/aggregated",
  "/MLD.+\\..+SubjectID_CellMarker_Tissue_TimepointMonths_SC-filtered",
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

mld_data_mod <- mld_data |>
  filter(as.numeric(TimepointMonths) > 0,
         CellMarker %in% included_markers) |>
  mutate(
    FU = if_else(as.numeric(TimepointMonths) < 24,
                 "Early", "Late")
  )

mld_data_reagg <- mld_data_mod |>
  left_join(blood_lineages_default(), by = "CellMarker") |>
  group_by(across(
    all_of(c(mandatory_IS_vars(), annotation_IS_vars(),
             "SubjectID", "Tissue", "CellType", "FU")),
  )) |>
  summarise(
    FE = sum(fragmentEstimate, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(
    CellType %in% c("Myeloid", "B", "T", "CD34", "Erythroid")
  )

### Flag matrices ----
mld_flag_folder <- fs::path("/home/giulia/HSPC_dynamics_data",
                            "MLD_flag_matrices")

mld_flag_data_BM <- purrr::map(fs::dir_ls(fs::path(mld_flag_folder, "BM")), ~{
  try({
    read.delim(file = .x, row.names = 1, sep = "\t") |>
      tibble::as_tibble()
  })
})
mld_flag_data_BM <- mld_flag_data_BM[purrr::map_lgl(mld_flag_data_BM, 
                                                    ~is.data.frame(.x))]
mld_flag_data_BM <- mld_flag_data_BM[-10] # MLD11 excluded - missing info

mld_flag_data_BM_mod <- mld_flag_data_BM |>
  purrr::list_rbind() |>
  mutate(
    Tissue = "BM"
  ) |>
  tidyr::separate(
    col = CloneID,
    into = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"),
    convert = TRUE, sep = "_"
  ) |>
  mutate(
    chr = stringr::str_remove(chr, "chr")
  )

mld_flag_data_PB <- purrr::map(fs::dir_ls(fs::path(mld_flag_folder, "PB")), ~{
  try({
    read.delim(file = .x, row.names = 1, sep = "\t") |>
      tibble::as_tibble()
  })
})
mld_flag_data_PB <- mld_flag_data_PB[purrr::map_lgl(mld_flag_data_PB, 
                                                    ~is.data.frame(.x))]

mld_flag_data_PB_mod <- mld_flag_data_PB |>
  purrr::list_rbind() |>
  mutate(
    Tissue = "PB"
  ) |>
  tidyr::separate(
    col = CloneID,
    into = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"),
    convert = TRUE, sep = "_"
  ) |>
  mutate(
    chr = stringr::str_remove(chr, "chr")
  )
mld_flag_data_PB_mod[45145, ]$GeneName <- "XLOC_008559"
mld_flag_data_PB_mod[45145, ]$GeneStrand <- "+"

mld_flag_data <- mld_flag_data_BM_mod |>
  bind_rows(mld_flag_data_PB_mod)

### Abundance ----
mld_abund <- compute_abundance(
  mld_data_reagg, columns = "FE",
  key = c("SubjectID", "Tissue", "CellType", "FU")
)
mld_abund_classic <- compute_abundance(
  mld_data_mod |>
    select(-FU),
  columns = "fragmentEstimate",
  key = c("SubjectID", "Tissue", "CellMarker", "TimepointMonths")
)
mld_abund_classic <- mld_abund_classic |>
  semi_join(mld_flag_data, by = mandatory_IS_vars())

mld_abund_classic |>
  readr::write_tsv(
    file = fs::path(ab_curves_fold, 
                    paste("MLD.20230713",
                          "abundance-IS-selection",
                          "SubjectID_CellMarker_Tissue_TimepointMonths.tsv.gz",
                          sep = "_"))
  )

### Statistics per IS ----
mld_ab_stats <- mld_abund |>
  group_by(across(all_of(c(mandatory_IS_vars(), annotation_IS_vars(),
                           "SubjectID", "Tissue", "FU")))) |>
  summarise(
    AB_min = min(FE_PercAbundance),
    AB_max = max(FE_PercAbundance),
    AB_mean = mean(FE_PercAbundance),
    AB_median = median(FE_PercAbundance),
    AB_sd = sd(FE_PercAbundance),
    AB_se = sd(FE_PercAbundance) / sqrt(length(FE_PercAbundance)),
    .groups = "drop"
  )

### Merge flag info and plot ----
mld_ab_flags <- mld_flag_data |>
  left_join(mld_ab_stats)
mld_ab_flags |>
  readr::write_tsv(
    file = fs::path(ab_curves_fold, 
                    "MLD.20230607_abundance_stats_ALL-patients-tissues.tsv"),
    na = ""
  )

mld_plots_patients <- get_plots(mld_ab_flags, "MLD")
save_plots(mld_plots_patients, ab_curves_fold, width = 20, height = 5)

mld_plots_patients_log <- get_plots(mld_ab_flags, "MLD", x_scale = "log")
save_plots(mld_plots_patients_log, ab_curves_fold, width = 20, height = 5,
           type = "log")

mld_plot_overview <- ggplot(mld_ab_flags, 
                            aes(x = AB_mean, y = FU, fill = FU)) +
  geom_density_ridges(alpha = 0.8, size = 0.3) +
  facet_grid(Tissue ~ StateChangeLabel) +
  labs(y = "", x = "Mean % abundance", 
       title = "MLD - all patients") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_ridges() +
  theme(legend.position = "none")

ggsave(
  plot = mld_plot_overview,
  filename = paste0("MLD.abundance_curves-all_patients", ".pdf"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)
ggsave(
  plot = mld_plot_overview,
  filename = paste0("MLD.abundance_curves-all_patients", ".png"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)

mld_plot_overview_log10 <- ggplot(mld_ab_flags, 
                            aes(x = AB_mean, y = FU, fill = FU)) +
  geom_density_ridges(alpha = 0.8, size = 0.3) +
  facet_grid(Tissue ~ StateChangeLabel) +
  labs(y = "", x = "Mean % abundance (log10)", 
       title = "MLD - all patients") +
  scale_x_log10() +
  theme_ridges() +
  theme(legend.position = "none")
ggsave(
  plot = mld_plot_overview_log10,
  filename = paste0("MLD.abundance_curves-all_patients_log10-x", ".pdf"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)
ggsave(
  plot = mld_plot_overview_log10,
  filename = paste0("MLD.abundance_curves-all_patients_log10-x", ".png"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)

### Delta and fold change ----
mld_mean_deltaFC <- mld_ab_flags |>
  select(-AB_min, -AB_max, -AB_median, -AB_sd, -AB_se) |>
  tidyr::pivot_wider(names_from = FU, values_from = AB_mean) |>
  mutate(
    delta = Late - Early,
    fold_change = Late / Early
  )
mld_mean_deltaFC |>
  readr::write_tsv(
    file = fs::path(ab_curves_fold, 
                    "MLD.20230609_mean-abund-change_ALL-pt-tissues.tsv"),
    na = ""
  )
mld_plots_delta_fc <- plot_delta_fc(mld_mean_deltaFC, "MLD", x_scale = "log")
save_plots(mld_plots_delta_fc, ab_curves_fold, width = 20, height = 5,
           type = "log", plot_type = "deltas")

mld_plot_overview_delta_log10 <- ggplot(mld_mean_deltaFC,
                    aes(x = delta)
) +
  geom_density() +
  theme_bw() +
  labs(x = "Log10 Delta % mean abundance (late - early)",
       y = "",
       title = "MLD - all patients") +
  scale_x_log10() + 
  facet_grid(Tissue ~ StateChangeLabel)

mld_plot_overview_fc_log10 <- ggplot(mld_mean_deltaFC,
                                        aes(x = fold_change)
) +
  geom_density() +
  theme_bw() +
  labs(x = "Log10 fold change % mean abundance (late/early)",
       y = "",
       title = "MLD - all patients") +
  scale_x_log10() + 
  facet_grid(Tissue ~ StateChangeLabel)
ggsave(
  plot = mld_plot_overview_delta_log10,
  filename = paste0("MLD.delta-curves_-all_patients_log10-x", ".pdf"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)
ggsave(
  plot = mld_plot_overview_delta_log10,
  filename = paste0("MLD.delta-curves_-all_patients_log10-x", ".png"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)
ggsave(
  plot = mld_plot_overview_fc_log10,
  filename = paste0("MLD.FC-curves_-all_patients_log10-x", ".pdf"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)
ggsave(
  plot = mld_plot_overview_fc_log10,
  filename = paste0("MLD.FC-curves_-all_patients_log10-x", ".png"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)

## BTHAL ----
bthal_data_folder <-
  fs::path("/mnt/CerberoWorkspaceISA/BTHAL/archive/202112")
bthal_matrix_regex <- paste0(
  "BTHAL[0-9]{3}.+_SubjectID_CellMarker_Tissue_",
  "TimepointMonths_SC-filtered_fragmentEstimate"
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

bthal_data_mod <- bthal_data |>
  filter(as.numeric(TimepointMonths) > 0,
         CellMarker %in% included_markers) |>
  mutate(
    FU = if_else(as.numeric(TimepointMonths) < 24,
                 "Early", "Late")
  )

bthal_data_reagg <- bthal_data_mod |>
  left_join(blood_lineages_default(), by = "CellMarker") |>
  group_by(across(
    all_of(c(mandatory_IS_vars(), annotation_IS_vars(),
             "SubjectID", "Tissue", "CellType", "FU")),
  )) |>
  summarise(
    FE = sum(fragmentEstimate, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(
    CellType %in% c("Myeloid", "B", "T", "CD34", "Erythroid")
  )

### Flag matrices ----
bthal_flag_folder <- fs::path("/home/giulia/HSPC_dynamics_data",
                            "BTHAL_flag_matrices")

bthal_flag_data_BM <- purrr::map(fs::dir_ls(
  fs::path(bthal_flag_folder, "BM")), ~{
  try({
    read.delim(file = .x, row.names = 1, sep = "\t") |>
      tibble::as_tibble()
  })
})

bthal_flag_data_BM_mod <- bthal_flag_data_BM |>
  purrr::list_rbind() |>
  mutate(
    Tissue = "BM"
  ) |>
  tidyr::separate(
    col = CloneID,
    into = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"),
    convert = TRUE, sep = "_"
  ) |>
  mutate(
    chr = stringr::str_remove(chr, "chr")
  )

bthal_flag_data_PB <- purrr::map(fs::dir_ls(
  fs::path(bthal_flag_folder, "PB")), ~{
  try({
    read.delim(file = .x, row.names = 1, sep = "\t") |>
      tibble::as_tibble()
  })
})

bthal_flag_data_PB_mod <- bthal_flag_data_PB |>
  purrr::list_rbind() |>
  mutate(
    Tissue = "PB"
  ) |>
  tidyr::separate(
    col = CloneID,
    into = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"),
    convert = TRUE, sep = "_"
  ) |>
  mutate(
    chr = stringr::str_remove(chr, "chr")
  )

bthal_flag_data <- bthal_flag_data_BM_mod |>
  bind_rows(bthal_flag_data_PB_mod)

### Abundance ----
bthal_abund <- compute_abundance(
  bthal_data_reagg, columns = "FE",
  key = c("SubjectID", "Tissue", "CellType", "FU")
)

bthal_abund_classic <- compute_abundance(
  bthal_data_mod |>
    select(-FU),
  columns = "fragmentEstimate",
  key = c("SubjectID", "Tissue", "CellMarker", "TimepointMonths")
)
bthal_abund_classic <- bthal_abund_classic |>
  semi_join(bthal_flag_data, by = mandatory_IS_vars())

bthal_abund_classic |>
  readr::write_tsv(
    file = fs::path(ab_curves_fold, 
                    paste("BTHAL.20230713",
                          "abundance-IS-selection",
                          "SubjectID_CellMarker_Tissue_TimepointMonths.tsv.gz",
                          sep = "_"))
  )

### Statistics per IS ----
bthal_ab_stats <- bthal_abund |>
  group_by(across(all_of(c(mandatory_IS_vars(), annotation_IS_vars(),
                           "SubjectID", "Tissue", "FU")))) |>
  summarise(
    AB_min = min(FE_PercAbundance),
    AB_max = max(FE_PercAbundance),
    AB_mean = mean(FE_PercAbundance),
    AB_median = median(FE_PercAbundance),
    AB_sd = sd(FE_PercAbundance),
    AB_se = sd(FE_PercAbundance) / sqrt(length(FE_PercAbundance)),
    .groups = "drop"
  )

### Merge flag info and plot ----
bthal_ab_flags <- bthal_flag_data |>
  left_join(bthal_ab_stats)
bthal_ab_flags |>
  readr::write_tsv(
    file = fs::path(ab_curves_fold, 
                    "BTHAL.20230607_abundance_stats_ALL-patients-tissues.tsv"),
    na = ""
  )

bthal_plots_patients <- get_plots(bthal_ab_flags, "BTHAL")
save_plots(bthal_plots_patients, ab_curves_fold, width = 20, height = 5)

bthal_plots_patients_log <- get_plots(bthal_ab_flags, "BTHAL", x_scale = "log")
save_plots(bthal_plots_patients_log, ab_curves_fold, width = 25, height = 5,
           type = "log")

bthal_plot_overview <- ggplot(bthal_ab_flags, 
                            aes(x = AB_mean, y = FU, fill = FU)) +
  geom_density_ridges(alpha = 0.8, size = 0.3) +
  facet_grid(Tissue ~ StateChangeLabel) +
  labs(y = "", x = "Mean % abundance", 
       title = "BTHAL - all patients") +
  scale_x_continuous(limits = c(0,1)) +
  theme_ridges() +
  theme(legend.position = "none")

ggsave(
  plot = bthal_plot_overview,
  filename = paste0("BTHAL.abundance_curves-all_patients", ".pdf"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)
ggsave(
  plot = bthal_plot_overview,
  filename = paste0("BTHAL.abundance_curves-all_patients", ".png"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)

bthal_plot_overview_log10 <- ggplot(bthal_ab_flags, 
                              aes(x = AB_mean, y = FU, fill = FU)) +
  geom_density_ridges(alpha = 0.8, size = 0.3) +
  facet_grid(Tissue ~ StateChangeLabel) +
  labs(y = "", x = "Mean % abundance (log10)", 
       title = "BTHAL - all patients") +
  scale_x_log10() +
  theme_ridges() +
  theme(legend.position = "none")

ggsave(
  plot = bthal_plot_overview_log10,
  filename = paste0("BTHAL.abundance_curves-all_patients_log10-x", ".pdf"),
  path = ab_curves_fold,
  width = 25,
  height = 5
)
ggsave(
  plot = bthal_plot_overview_log10,
  filename = paste0("BTHAL.abundance_curves-all_patients_log10-x", ".png"),
  path = ab_curves_fold,
  width = 25,
  height = 5
)

### Delta and fold change ----
bthal_mean_deltaFC <- bthal_ab_flags |>
  select(-AB_min, -AB_max, -AB_median, -AB_sd, -AB_se) |>
  tidyr::pivot_wider(names_from = FU, values_from = AB_mean) |>
  mutate(
    delta = Late - Early,
    fold_change = Late / Early
  )
bthal_mean_deltaFC |>
  readr::write_tsv(
    file = fs::path(ab_curves_fold, 
                    "BTHAL.20230609_mean-abund-change_ALL-pt-tissues.tsv"),
    na = ""
  )
bthal_plots_delta_fc <- plot_delta_fc(bthal_mean_deltaFC, "BTHAL", x_scale = "log")
save_plots(bthal_plots_delta_fc, ab_curves_fold, width = 20, height = 5,
           type = "log", plot_type = "deltas")

bthal_plot_overview_delta_log10 <- ggplot(bthal_mean_deltaFC,
                                        aes(x = delta)
) +
  geom_density() +
  theme_bw() +
  labs(x = "Log10 Delta % mean abundance (late - early)",
       y = "",
       title = "BTHAL - all patients") +
  scale_x_log10() + 
  facet_grid(Tissue ~ StateChangeLabel)

bthal_plot_overview_fc_log10 <- ggplot(bthal_mean_deltaFC,
                                     aes(x = fold_change)
) +
  geom_density() +
  theme_bw() +
  labs(x = "Log10 fold change % mean abundance (late/early)",
       y = "",
       title = "BTHAL - all patients") +
  scale_x_log10() + 
  facet_grid(Tissue ~ StateChangeLabel)
ggsave(
  plot = bthal_plot_overview_delta_log10,
  filename = paste0("BTHAL.delta-curves_-all_patients_log10-x", ".pdf"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)
ggsave(
  plot = bthal_plot_overview_delta_log10,
  filename = paste0("BTHAL.delta-curves_-all_patients_log10-x", ".png"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)
ggsave(
  plot = bthal_plot_overview_fc_log10,
  filename = paste0("BTHAL.FC-curves_-all_patients_log10-x", ".pdf"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)
ggsave(
  plot = bthal_plot_overview_fc_log10,
  filename = paste0("BTHAL.FC-curves_-all_patients_log10-x", ".png"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)


## WAS ----
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

was_data_mod <- was_data |>
  filter(as.numeric(TimepointMonths) > 0,
         CellMarker %in% included_markers) |>
  mutate(
    FU = if_else(as.numeric(TimepointMonths) < 24,
                 "Early", "Late")
  )

was_data_reagg <- was_data_mod |>
  left_join(blood_lineages_default(), by = "CellMarker") |>
  group_by(across(
    all_of(c(mandatory_IS_vars(), annotation_IS_vars(),
             "SubjectID", "Tissue", "CellType", "FU")),
  )) |>
  summarise(
    FE = sum(fragmentEstimate, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(
    CellType %in% c("Myeloid", "B", "T", "CD34", "Erythroid")
  )

### Flag matrices ----
was_flag_folder <- fs::path("/home/giulia/HSPC_dynamics_data",
                              "WAS_flag_matrices")

was_flag_data_BM <- purrr::map(fs::dir_ls(
  fs::path(was_flag_folder, "BM")), ~{
    try({
      read.delim(file = .x, row.names = 1, sep = "\t") |>
        tibble::as_tibble()
    })
  })


was_flag_data_BM_mod <- was_flag_data_BM |>
  purrr::list_rbind() |>
  mutate(
    Tissue = "BM"
  ) |>
  tidyr::separate(
    col = CloneID,
    into = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"),
    convert = TRUE, sep = "_"
  ) |>
  mutate(
    chr = stringr::str_remove(chr, "chr")
  )


was_flag_data_PB <- purrr::map(fs::dir_ls(
  fs::path(was_flag_folder, "PB")), ~{
    try({
      read.delim(file = .x, row.names = 1, sep = "\t") |>
        tibble::as_tibble()
    })
  })

was_flag_data_PB_mod <- was_flag_data_PB |>
  purrr::list_rbind() |>
  mutate(
    Tissue = "PB"
  ) |>
  tidyr::separate(
    col = CloneID,
    into = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"),
    convert = TRUE, sep = "_"
  ) |>
  mutate(
    chr = stringr::str_remove(chr, "chr")
  )
was_flag_data_PB_mod[1735, ]$GeneName <- "C4B_2"
was_flag_data_PB_mod[1735, ]$GeneStrand <- "+"

was_flag_data <- was_flag_data_BM_mod |>
  bind_rows(was_flag_data_PB_mod)

### Abundance ----
was_abund <- compute_abundance(
  was_data_reagg, columns = "FE",
  key = c("SubjectID", "Tissue", "CellType", "FU")
)

was_abund_classic <- compute_abundance(
  was_data_mod |>
    select(-FU),
  columns = "fragmentEstimate",
  key = c("SubjectID", "Tissue", "CellMarker", "TimepointMonths")
)
was_abund_classic <- was_abund_classic |>
  semi_join(was_flag_data, by = mandatory_IS_vars())

was_abund_classic |>
  readr::write_tsv(
    file = fs::path(ab_curves_fold, 
                    paste("WAS.20230713",
                          "abundance-IS-selection",
                          "SubjectID_CellMarker_Tissue_TimepointMonths.tsv.gz",
                          sep = "_"))
  )

### Statistics per IS ----
was_ab_stats <- was_abund |>
  group_by(across(all_of(c(mandatory_IS_vars(), annotation_IS_vars(),
                           "SubjectID", "Tissue", "FU")))) |>
  summarise(
    AB_min = min(FE_PercAbundance),
    AB_max = max(FE_PercAbundance),
    AB_mean = mean(FE_PercAbundance),
    AB_median = median(FE_PercAbundance),
    AB_sd = sd(FE_PercAbundance),
    AB_se = sd(FE_PercAbundance) / sqrt(length(FE_PercAbundance)),
    .groups = "drop"
  )

### Merge flag info and plot ----
was_ab_flags <- was_flag_data |>
  left_join(was_ab_stats)
was_ab_flags |>
  readr::write_tsv(
    file = fs::path(ab_curves_fold, 
                    "WAS.20230607_abundance_stats_ALL-patients-tissues.tsv"),
    na = ""
  )

was_ab_flags <- was_ab_flags |>
  filter(!is.na(FU))

was_ab_flags |>
  readr::write_tsv(
    file = fs::path(ab_curves_fold, 
                    "WAS.20230609_abundance_stats_ALL-patients-tissues_NO-NA.tsv"),
    na = ""
  )

was_plots_patients <- get_plots(was_ab_flags, "WAS")
save_plots(was_plots_patients, ab_curves_fold, width = 20, height = 5)

was_plots_patients_log <- get_plots(was_ab_flags, "WAS", x_scale = "log",
                                    stat_curves = "density")
save_plots(was_plots_patients_log, ab_curves_fold, width = 20, height = 5,
           type = "log")

was_plot_overview <- ggplot(was_ab_flags, 
                              aes(x = AB_mean, y = FU, fill = FU)) +
  geom_density_ridges(alpha = 0.8, size = 0.3) +
  facet_grid(Tissue ~ StateChangeLabel) +
  labs(y = "", x = "Mean % abundance", 
       title = "WAS - all patients") +
  scale_x_continuous(limits = c(0,1)) +
  theme_ridges() +
  theme(legend.position = "none")

ggsave(
  plot = was_plot_overview,
  filename = paste0("WAS.abundance_curves-all_patients", ".pdf"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)
ggsave(
  plot = was_plot_overview,
  filename = paste0("WAS.abundance_curves-all_patients", ".png"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)

was_plot_overview_log10 <- ggplot(was_ab_flags, 
                            aes(x = AB_mean, y = FU, fill = FU,
                                height = after_stat(density))) +
  geom_density_ridges(alpha = 0.8, size = 0.3, stat = "density") +
  facet_grid(Tissue ~ StateChangeLabel) +
  labs(y = "", x = "Mean % abundance (log10)", 
       title = "WAS - all patients") +
  scale_x_log10() +
  theme_ridges() +
  theme(legend.position = "none")

ggsave(
  plot = was_plot_overview_log10,
  filename = paste0("WAS.abundance_curves-all_patients_log10-x", ".pdf"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)
ggsave(
  plot = was_plot_overview_log10,
  filename = paste0("WAS.abundance_curves-all_patients_log10-x", ".png"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)

### Delta and fold change ----
was_mean_deltaFC <- was_ab_flags |>
  select(-AB_min, -AB_max, -AB_median, -AB_sd, -AB_se) |>
  tidyr::pivot_wider(names_from = FU, values_from = AB_mean) |>
  mutate(
    delta = Late - Early,
    fold_change = Late / Early
  )
was_mean_deltaFC |>
  readr::write_tsv(
    file = fs::path(ab_curves_fold, 
                    "WAS.20230609_mean-abund-change_ALL-pt-tissues.tsv"),
    na = ""
  )
was_plots_delta_fc <- plot_delta_fc(was_mean_deltaFC, "WAS", x_scale = "log")
save_plots(was_plots_delta_fc, ab_curves_fold, width = 20, height = 5,
           type = "log", plot_type = "deltas")

was_plot_overview_delta_log10 <- ggplot(was_mean_deltaFC,
                                        aes(x = delta)
) +
  geom_density() +
  theme_bw() +
  labs(x = "Log10 Delta % mean abundance (late - early)",
       y = "",
       title = "WAS - all patients") +
  scale_x_log10() + 
  facet_grid(Tissue ~ StateChangeLabel)

was_plot_overview_fc_log10 <- ggplot(was_mean_deltaFC,
                                     aes(x = fold_change)
) +
  geom_density() +
  theme_bw() +
  labs(x = "Log10 fold change % mean abundance (late/early)",
       y = "",
       title = "WAS - all patients") +
  scale_x_log10() + 
  facet_grid(Tissue ~ StateChangeLabel)
ggsave(
  plot = was_plot_overview_delta_log10,
  filename = paste0("WAS.delta-curves_-all_patients_log10-x", ".pdf"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)
ggsave(
  plot = was_plot_overview_delta_log10,
  filename = paste0("WAS.delta-curves_-all_patients_log10-x", ".png"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)
ggsave(
  plot = was_plot_overview_fc_log10,
  filename = paste0("WAS.FC-curves_-all_patients_log10-x", ".pdf"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)
ggsave(
  plot = was_plot_overview_fc_log10,
  filename = paste0("WAS.FC-curves_-all_patients_log10-x", ".png"),
  path = ab_curves_fold,
  width = 20,
  height = 5
)

# Abundances comparison --------------------------------------------------------
# 2023-07-10
library(ggpubr)
perc_spacer_stats <- 0.2
perc_starting_max_stats <- 0.75
### Classes with UNI will contain all uni states without distinction
classes <- c("MULTI-MULTI", "MULTI-UNI", "UNI-UNI")
comparisons <- gtools::combinations(length(classes), r = 2, v = classes) |> 
  t() |> as.data.frame() |> as.list()

## MLD -----------------
mld_ab_classes <- mld_ab_flags |>
  mutate(
    class = if_else(
      Early_StateLabel == "Multi" & Late_StateLabel == "Multi",
      classes[1],
      if_else(
        Early_StateLabel == "Multi" & 
          stringr::str_starts(Late_StateLabel, "Uni"),
        classes[2],
        classes[3]
      )
    )
  )

mld_ab_classes_means <- mld_ab_classes |>
  group_by(across(all_of(c("SubjectID",
                           "class", "Tissue", "FU")))) |>
  summarise(
    mean_AB = mean(AB_mean, na.rm = TRUE),
    .groups = "drop"
  )

mld_ab_classes_plot <- ggplot(mld_ab_classes_means, 
                              aes(x = class, y = mean_AB,
                                 shape = class)) +
  geom_jitter(alpha = 0.6, size = 4) +
  scale_y_continuous(limits = c(0, 1.5),
                     labels = scales::label_percent(scale = 1)) +
  theme_bw() +
  facet_grid(Tissue ~ FU) +
  labs(x = "", y = "Mean abundance per patient",
       title = "MLD abundances comparison") +
  theme(legend.position = "none", text = element_text(size = 16),
        strip.text = element_text(face = "bold", size = 16)) +
  geom_pwc(method = "t_test", p.adjust.method = "fdr",
           label = "p.adj", y.position = perc_starting_max_stats,
           step.increase = 0.3)

mld_ab_classes_plot <- ggadjust_pvalue(mld_ab_classes_plot,
                                       p.adjust.method = "fdr",
                                       label = "p.adj.signif")

ggsave(
  plot = mld_ab_classes_plot,
  path = ab_curves_fold,
  filename = "MLD.abundances_means_comparison.png",
  width = 9,
  height = 7
)
ggsave(
  plot = mld_ab_classes_plot,
  path = ab_curves_fold,
  filename = "MLD.abundances_means_comparison.pdf",
  width = 9,
  height = 7
)

mld_ab_classes_plot2 <- ggplot(mld_ab_classes_means, 
                               aes(x = FU, y = mean_AB,
                                 fill = FU,
                                 color = FU)) +
  geom_jitter(alpha = 0.6, size = 4) +
  scale_y_continuous(limits = c(0, 1.5), 
                     labels = scales::label_percent(scale = 1)) +
  theme_bw() +
  facet_grid(Tissue ~ class) +
  labs(x = "", y = "Mean abundance per patient", 
       title = "MLD abundances comparison") +
  theme(legend.position = "none", text = element_text(size = 16),
        strip.text = element_text(face = "bold", size = 16)) +
  geom_pwc(method = "t_test", p.adjust.method = "fdr",
           label = "p.adj", y.position = perc_starting_max_stats + 0.5,
           step.increase = 0.3)

mld_ab_classes_plot2 <- ggadjust_pvalue(mld_ab_classes_plot2,
                                        p.adjust.method = "fdr",
                                        label = "p.adj.signif")

ggsave(
  plot = mld_ab_classes_plot2,
  path = ab_curves_fold,
  filename = "MLD.abundances_means_comparison_FU-x.png",
  width = 10,
  height = 7
)
ggsave(
  plot = mld_ab_classes_plot2,
  path = ab_curves_fold,
  filename = "MLD.abundances_means_comparison_FU-x.pdf",
  width = 10,
  height = 7
)

## BTHAL ---------------
bthal_ab_classes <- bthal_ab_flags |>
  mutate(
    class = if_else(
      Early_StateLabel == "Multi" & Late_StateLabel == "Multi",
      classes[1],
      if_else(
        Early_StateLabel == "Multi" & 
          stringr::str_starts(Late_StateLabel, "Uni"),
        classes[2],
        classes[3]
      )
    )
  )

bthal_ab_classes_means <- bthal_ab_classes |>
  group_by(across(all_of(c("SubjectID",
                           "class", "Tissue", "FU")))) |>
  summarise(
    mean_AB = mean(AB_mean, na.rm = TRUE),
    .groups = "drop"
  )

bthal_ab_classes_plot <- ggplot(bthal_ab_classes_means, 
                                aes(x = class, y = mean_AB,
                                    shape = class)) +
  geom_jitter(alpha = 0.6, size = 4) +
  scale_y_continuous(limits = c(0, 1.5), 
                     labels = scales::label_percent(scale = 1)) +
  theme_bw() +
  facet_grid(Tissue ~ FU) +
  labs(x = "", y = "Mean abundance per patient", 
       title = "BTHAL abundances comparison") +
  theme(legend.position = "none", text = element_text(size = 16),
        strip.text = element_text(face = "bold", size = 16)) +
  geom_pwc(method = "t_test", p.adjust.method = "fdr",
           label = "p.adj", y.position = perc_starting_max_stats,
           step.increase = 1)

bthal_ab_classes_plot <- ggadjust_pvalue(bthal_ab_classes_plot,
                                       p.adjust.method = "fdr",
                                       label = "p.adj.signif")

ggsave(
  plot = bthal_ab_classes_plot,
  path = ab_curves_fold,
  filename = "BTHAL.abundances_means_comparison.png",
  width = 9,
  height = 7
)
ggsave(
  plot = bthal_ab_classes_plot,
  path = ab_curves_fold,
  filename = "BTHAL.abundances_means_comparison.pdf",
  width = 9,
  height = 7
)

bthal_ab_classes_plot2 <- ggplot(bthal_ab_classes_means, 
                                 aes(x = FU, y = mean_AB,
                                     fill = FU,
                                     color = FU)) +
  geom_jitter(alpha = 0.6, size = 4) +
  scale_y_continuous(limits = c(0, 1.5), 
                     labels = scales::label_percent(scale = 1)) +
  theme_bw() +
  facet_grid(Tissue ~ class) +
  labs(x = "", y = "Mean abundance per patient", 
       title = "BTHAL abundances comparison") +
  theme(legend.position = "none", text = element_text(size = 16),
        strip.text = element_text(face = "bold", size = 16)) +
  
  geom_pwc(method = "t_test", p.adjust.method = "fdr",
           label = "p.adj", y.position = perc_starting_max_stats,
           step.increase = 1)
bthal_ab_classes_plot2 <- ggadjust_pvalue(bthal_ab_classes_plot2,
                                         p.adjust.method = "fdr",
                                         label = "p.adj.signif")

ggsave(
  plot = bthal_ab_classes_plot2,
  path = ab_curves_fold,
  filename = "BTHAL.abundances_means_comparison_FU-x.png",
  width = 10,
  height = 7
)
ggsave(
  plot = bthal_ab_classes_plot2,
  path = ab_curves_fold,
  filename = "BTHAL.abundances_means_comparison_FU-x.pdf",
  width = 10,
  height = 7
)

## WAS -----------------
was_ab_classes <- was_ab_flags |>
  mutate(
    class = if_else(
      Early_StateLabel == "Multi" & Late_StateLabel == "Multi",
      classes[1],
      if_else(
        Early_StateLabel == "Multi" & 
          stringr::str_starts(Late_StateLabel, "Uni"),
        classes[2],
        classes[3]
      )
    )
  )

was_ab_classes_means <- was_ab_classes |>
  group_by(across(all_of(c("SubjectID",
                           "class", "Tissue", "FU")))) |>
  summarise(
    mean_AB = mean(AB_mean, na.rm = TRUE),
    .groups = "drop"
  ) 


was_ab_classes_plot <- ggplot(was_ab_classes_means, 
                              aes(x = class, y = mean_AB,
                                  shape = class)) +
  geom_jitter(alpha = 0.6, size = 4) +
  scale_y_continuous(limits = c(0, 1.5),
                     labels = scales::label_percent(scale = 1)) +
  theme_bw() +
  facet_grid(Tissue ~ FU) +
  labs(x = "", y = "Mean abundance per patient",
       title = "WAS abundances comparison") +
  theme(legend.position = "none", text = element_text(size = 16),
        strip.text = element_text(face = "bold", size = 16)) +
  geom_pwc(method = "t_test", p.adjust.method = "fdr",
           label = "p.adj", y.position = perc_starting_max_stats,
           step.increase = 1)

was_ab_classes_plot <- ggadjust_pvalue(was_ab_classes_plot, 
                        p.adjust.method = "fdr", 
                        label = "p.adj.signif")


ggsave(
  plot = was_ab_classes_plot,
  path = ab_curves_fold,
  filename = "WAS.abundances_means_comparison.png",
  width = 9,
  height = 7
)
ggsave(
  plot = was_ab_classes_plot,
  path = ab_curves_fold,
  filename = "WAS.abundances_means_comparison.pdf",
  width = 9,
  height = 7
)

was_ab_classes_plot2 <- ggplot(was_ab_classes_means, aes(x = FU, y = mean_AB,
                                                             fill = FU,
                                                             color = FU)) +
  geom_jitter(alpha = 0.6, size = 4) +
  scale_y_continuous(limits = c(0, 1.5), 
                     labels = scales::label_percent(scale = 1)) +
  theme_bw() +
  facet_grid(Tissue ~ class) +
  labs(x = "", y = "Mean abundance per patient", 
       title = "WAS abundances comparison") +
  theme(legend.position = "none", text = element_text(size = 16),
        strip.text = element_text(face = "bold", size = 16)) +
  geom_pwc(method = "t_test", p.adjust.method = "fdr",
           label = "p.adj", y.position = perc_starting_max_stats,
           step.increase = 1)

was_ab_classes_plot2 <- ggadjust_pvalue(was_ab_classes_plot2, 
                                       p.adjust.method = "fdr", 
                                       label = "p.adj.signif")

ggsave(
  plot = was_ab_classes_plot2,
  path = ab_curves_fold,
  filename = "WAS.abundances_means_comparison_FU-x.png",
  width = 10,
  height = 7
)
ggsave(
  plot = was_ab_classes_plot2,
  path = ab_curves_fold,
  filename = "WAS.abundances_means_comparison_FU-x.pdf",
  width = 10,
  height = 7
)


## ALL TRIALS -----------------
ab_classes_means <- mld_ab_classes_means |>
  mutate(
    trial = "MLD"
  ) |>
  bind_rows(bthal_ab_classes_means |>
              mutate(
                trial = "BTHAL"
              )) |>
  bind_rows(was_ab_classes_means |>
              mutate(
                trial = "WAS"
              ))

ab_classes_plot <- ggplot(ab_classes_means, 
                          aes(x = class, y = mean_AB, 
                              shape = class)) +
  geom_jitter(aes(color = trial, fill = trial), alpha = 0.6, size = 4) +
  scale_y_continuous(limits = c(0, 1.5), 
                     labels = scales::label_percent(scale = 1)) +
  scale_color_manual(values = c(MLD = "darkblue", WAS = "forestgreen",
                                BTHAL = "firebrick")) +
  scale_fill_manual(values = c(MLD = "darkblue", WAS = "forestgreen",
                                BTHAL = "firebrick")) +
  theme_bw() +
  facet_grid(Tissue ~ FU) +
  labs(x = "", y = "Mean abundance per patient", 
       title = "All trials abundances comparison",
       fill = "Clinical trial", color = "Clinical trial") +
  theme(legend.position = "bottom", text = element_text(size = 16),
        strip.text = element_text(face = "bold", size = 16)) +
  guides(shape = "none") +
  geom_pwc(method = "t_test", p.adjust.method = "fdr",
           label = "p.adj", y.position = perc_starting_max_stats,
           step.increase = 0.3, group.by = "x.var")

ab_classes_plot <- ggadjust_pvalue(ab_classes_plot, 
                                        p.adjust.method = "fdr", 
                                        label = "p.adj.signif")

ggsave(
  plot = ab_classes_plot,
  path = ab_curves_fold,
  filename = "ALL.abundances_means_comparison.png",
  width = 9,
  height = 7
)
ggsave(
  plot = ab_classes_plot,
  path = ab_curves_fold,
  filename = "ALL.abundances_means_comparison.pdf",
  width = 9,
  height = 7
)

ab_classes_plot2 <- ggplot(ab_classes_means, aes(x = FU, y = mean_AB)) +
  geom_jitter(aes(fill = trial,
                  color = trial), alpha = 0.6, size = 4) +
  scale_y_continuous(limits = c(0, 1.5), 
                     labels = scales::label_percent(scale = 1)) +
  scale_color_manual(values = c(MLD = "darkblue", WAS = "forestgreen",
                                BTHAL = "firebrick")) +
  scale_fill_manual(values = c(MLD = "darkblue", WAS = "forestgreen",
                               BTHAL = "firebrick")) +
  theme_bw() +
  facet_grid(Tissue ~ class) +
  labs(x = "", y = "Mean abundance per patient", 
       title = "All trials abundances comparison",
       fill = "Clinical trial", color = "Clinical trial") +
  theme(legend.position = "bottom", text = element_text(size = 16),
        strip.text = element_text(face = "bold", size = 16)) +
  geom_pwc(method = "t_test", p.adjust.method = "fdr",
           label = "p.adj", y.position = perc_starting_max_stats,
           step.increase = 0.3, group.by = "x.var")

ab_classes_plot2 <- ggadjust_pvalue(ab_classes_plot2, 
                                   p.adjust.method = "fdr", 
                                   label = "p.adj.signif")

ggsave(
  plot = ab_classes_plot2,
  path = ab_curves_fold,
  filename = "ALL.abundances_means_comparison_FU-x.png",
  width = 10,
  height = 7
)
ggsave(
  plot = ab_classes_plot2,
  path = ab_curves_fold,
  filename = "ALL.abundances_means_comparison_FU-x.pdf",
  width = 10,
  height = 7
)

save.image(file = "./HSPC_dynamics_rev1_abundance.RData")
