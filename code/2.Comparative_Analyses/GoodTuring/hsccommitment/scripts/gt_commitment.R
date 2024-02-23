main_path = "~/Dropbox/shared/HSC dynamics/"

setwd(main_path)
source(paste0(main_path, "r_functions/gt_helper_fns.R"))
library(tidyverse)

# patients_without_erythroid = c("MLDCUP04", "MLDHE01", "MLDCUP05", "MLDCUP03", "MLDC02", "WAS1009")
# patients_uncertainty = c("MLD12")

## Load datasets ####

old_stats = readxl::read_xlsx(paste0(main_path, "hsccommitment/data/MLD_WAS_BTHAL.AllPatients.hlfu_34MyBTEry.global_allstudies_profile_skewing_fulldata_notNA_stats.xlsx")) %>%
  dplyr::rename(TimepointMonths=Timepoint) %>% 
  dplyr::bind_rows(
    read.csv(paste0(main_path, "hsccommitment/data/XSCID.AllPatients.hlfu_excl34MyBTEry.flag_matrix_relabeled_correctedDC_freqs_melt_merge.tsv"), sep="\t") %>%
      dplyr::rename(TimepointMonths=Timepoint, ClinicalStudy=ProjectID) %>% 
      dplyr::mutate(Tissue="PB")
  )

## collect classnames by timepoint
classification = read.csv(paste0(main_path, "hsccommitment/datasets/flags_timepoint_wXSCID.csv")) %>% 
  tibble::as_tibble() %>% 
  
  ## XSCID only PB? no UniErythroid
  dplyr::mutate(Tissue=dplyr::case_when(Tissue=="" ~ "PB", .default=Tissue)) %>% 
  
  # add classnames
  dplyr::full_join(
    readxl::read_xlsx(paste0(main_path, "hsccommitment/data/MLD.hlfu_34MyBTEry_noLy_no34exclusive.flag_matrix.bagcases.xlsx")) %>% 
      dplyr::mutate(Tissue="BM") %>% 
      dplyr::bind_rows(
        readxl::read_xlsx(paste0(main_path, "hsccommitment/data/MLD.hlfu_34PBMyBT_noLy_no34exclusive.flag_matrix.bagcases.xlsx")) %>% 
          dplyr::mutate(Tissue="PB")
      ) %>% dplyr::rename(classcode=Label) %>% 
      dplyr::select(classcode, ClassName, Tissue, colorcode)
  )

## load abundances for all ISs
abundances = read.csv(paste0(main_path, "hsccommitment/datasets/abundances_all_ISs_wXSCID.csv")) %>% 
  dplyr::mutate(IS=paste0("chr", paste(chr, integration_locus, strand, GeneName, GeneStrand, sep="_")),
                Tissue=dplyr::case_when(is.na(Tissue)~"PB", .default=Tissue)) %>% 
  dplyr::select(-chr, -integration_locus, -strand, -GeneName, - GeneStrand) %>% 
  dplyr::select(IS, dplyr::everything()) %>% tibble::as_tibble() %>% 
  
  ## add celltypes -> some markers removed
  dplyr::left_join(
    readxl::read_xlsx(paste0(main_path, "paper docs/Calabria et al - HSC dynamics - Extended Data Tables 1-2.xlsx"), sheet=1) %>% 
      dplyr::select(CellMarker, Lineage) %>% unique() %>% dplyr::rename(CellType=Lineage),
    by="CellMarker"
  ) %>% 
  dplyr::mutate(abundance=dplyr::case_when(!is.na(seqCount) ~ seqCount,
                                           .default=fragmentEstimate))


input_perc = classification %>% dplyr::filter(ClassName!="Unobserved") %>% 
  ## add abundances -> removes some ISs from "classification"
  dplyr::inner_join(abundances, multiple="all", 
                   by=c("IS","SubjectID","ClinicalStudy",
                        "Tissue","TimepointMonths")) %>% 

  dplyr::select(IS, SubjectID, ClinicalStudy, Tissue, TimepointMonths, 
                CellType, CellMarker, abundance, ClassName, colorcode) %>% 

  ## sum abundance across markers
  dplyr::group_by(IS, SubjectID, ClinicalStudy, Tissue, TimepointMonths, CellType) %>% 
  dplyr::summarise(abundance=sum(abundance), 
                   ClassName=unique(ClassName), 
                   colorcode=unique(colorcode)) %>% 
  dplyr::ungroup() %>% unique()

## celltyper per class
celltype_list = input_perc %>% dplyr::select(ClassName) %>% unique() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(celltype_list=dplyr::case_when(
    ClassName=="Multilineage" ~ list(c("B","T","Myeloid","Erythroid","CD34")),
    .default = list(c("CD34",stringr::str_replace_all(ClassName, "Uni|_|Ly",""))))) %>%
  dplyr::ungroup()


res_perc = input_perc %>% 
  dplyr::inner_join(celltype_list) %>% 
  tidyr::pivot_wider(names_from="CellType",
                     values_from="abundance", 
                     values_fill=0) %>% 
  dplyr::group_by(SubjectID, ClinicalStudy, Tissue, TimepointMonths, 
                  ClassName, celltype_list, colorcode) %>%
  tidyr::nest() %>%
  dplyr::summarise(gt_richness=compute_richness(assemblages_df=data[[1]], gt=TRUE,
                                                subset_cols=unique(unlist(celltype_list[[1]]))),
                   emp_richness=compute_richness(assemblages_df=data[[1]], gt=FALSE,
                                                 subset_cols=unique(unlist(celltype_list[[1]])))) %>% 
  dplyr::ungroup()


## get denominator
tot_IS = input_perc %>% 
  # sum abundance over classes
  dplyr::group_by(IS, SubjectID, ClinicalStudy, Tissue, TimepointMonths) %>% 
  dplyr::summarise(abundance=sum(abundance)) %>% 
  dplyr::ungroup() %>% 
  
  dplyr::group_by(SubjectID, ClinicalStudy, Tissue, TimepointMonths) %>%
  tidyr::nest() %>% 
  dplyr::summarise(gt_ntot=compute_richness(assemblages_df=data[[1]],
                                            subset_cols="abundance", gt=T),
                   emp_ntot=compute_richness(assemblages_df=data[[1]],
                                             subset_cols="abundance", gt=F))


## Save dataframe ####
ratio_df = res_perc %>% 
  dplyr::full_join(tot_IS) %>% 
  dplyr::mutate(ratio_chao=gt_richness / gt_ntot*100,
                ratio_emp=emp_richness / emp_ntot*100) %>% 
  dplyr::select(-celltype_list)

write.csv(ratio_df, paste0(main_path, "hsccommitment/datasets/adj_ratios_commitment_wXSCID.csv"))

# ratio_df2 = read.csv(paste0(main_path, "hsccommitment/datasets/adj_ratios_commitment.csv"))

## Plot ####
ratio_df = read.csv(paste0(main_path, "hsccommitment/datasets/adj_ratios_commitment_wXSCID.csv"), row.names=1)
color_pal_perc = classification$colorcode %>% unique() %>%
  setNames(classification$ClassName %>% unique())

ratio_df %>% 
  tidyr::pivot_longer(cols=c("ratio_chao","ratio_emp"), 
                      names_to="method", values_to="ratio") %>% 

  ggplot() +
  
  geom_smooth(aes(x=TimepointMonths, y=ratio,
                  color=ClassName, fill=ClassName), 
              
              method="loess", 
              formula=y~log(x),
              stat="smooth",
              position="identity", alpha=0.4, level=0.75, se=T) +
  
  scale_color_manual(values=color_pal_perc) +
  scale_fill_manual(values=color_pal_perc) +
  
  ggh4x::facet_nested(method + Tissue ~ 
                        factor(ClinicalStudy, levels=c("MLD","WAS","BTHAL","XSCID")),
                      scales="free", space="free_x") +
  scale_x_continuous(breaks = seq(0, max(ratio_df$TimepointMonths, na.rm = T), 12)) +
  theme_bw() +
  theme(strip.text=element_text(face="bold", size=16)) +
  theme(legend.position="right") +
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16),
        axis.title=element_text(size=16), plot.title=element_text(size=20)) +
  labs(x="Months after gene therapy",
       y="% IS", 
       colour="ClassName", fill="ClassName")

ggsave(paste0(main_path, "hsccommitment/plots/gt_perc_ISs.pdf"), height=10, width=12)
ggsave(paste0(main_path, "hsccommitment/plots/gt_perc_ISs.png"), height=10, width=12, dpi=600)



color_pal_perc = old_stats$colorcode %>% unique() %>%
  setNames(old_stats$ClassName %>% unique())
old_stats %>% 
  ggplot() +
  
  geom_smooth(aes(x=TimepointMonths, y=CRonR, fill=Class, color=Class),
              method="loess",
              formula=y~log(x),
              stat="smooth",
              position="identity", alpha=0.4, level=0.75, se=T) +
  
  scale_color_manual(values=color_pal_perc) +
  scale_fill_manual(values=color_pal_perc) +
  
  ggh4x::facet_nested(Tissue ~ 
                        factor(ClinicalStudy, levels=c("MLD","WAS","BTHAL","XSCID")),
                      scales="free", space="free_x") +
  scale_x_continuous(breaks = seq(0, max(old_stats$TimepointMonths, na.rm = T), 12)) +
  theme_bw() +
  theme(strip.text=element_text(face="bold", size=16)) +
  theme(legend.position="right") +
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16),
        axis.title=element_text(size=16), plot.title=element_text(size=20)) +
  labs(x="Months after gene therapy",
       y="% IS", 
       colour="ClassName", fill="ClassName")

ggsave(paste0(main_path, "hsccommitment/plots/main.pdf"), height=10, width=12)
ggsave(paste0(main_path, "hsccommitment/plots/main.png"), height=10, width=12, dpi=600)




