main_path = "~/Dropbox/shared/HSC dynamics/"
save_path = "~/Dropbox/shared/HSC dynamics/cd34output/"
setwd(main_path)
source(paste0(main_path, "r_functions/gt_helper_fns.R"))

library(magrittr)
library(ggplot2)

# # Compute ratios ####
# patients_to_exclude = c("MLD22", "MLDCUP01", "MLDCUP02", "BTHAL001", "MLDC02")
patients_to_exclude = c()

old_stats = read.csv(paste0(main_path, "cd34output/datasets/cd34output_ISAnalytics.csv")) %>%
  dplyr::filter(grepl("BM$|CD34$", g1)) %>%
  tibble::as_tibble()


patients_infos = read.csv(paste0(main_path, "datasets/patients_age.csv"))


## load celltypes and colors
celltypes_markers = readxl::read_xlsx(paste0(main_path, "paper docs/blood_lineages_update.xlsx"), sheet=1) %>%
  dplyr::select(SuperGroup, AggregationGroup, CellType) %>% unique() %>%

  ## add colors
  dplyr::full_join(
    readxl::read_xlsx(paste0(main_path, "diversity/data/202112.Diversity.ByMarker.xlsx")) %>%
      dplyr::select(CellType, colorcode) %>% unique()
  ) %>% dplyr::filter(!is.na(colorcode))

celltypes_markers = celltypes_markers %>%
  dplyr::rename(AggrSuperGroup=SuperGroup) %>%
  dplyr::select(-AggregationGroup) %>%
  dplyr::bind_rows(celltypes_markers %>%
                     dplyr::rename(AggrSuperGroup=AggregationGroup) %>%
                     dplyr::select(-SuperGroup)) %>%
  unique()



## load abundances
## fragmentEstimate for XSCID is actually seqCount
abundances = read.csv(paste0(main_path, "cd34output/datasets/abundances_CD34.csv")) %>%
  dplyr::select(-chr, -integration_locus, -strand, -GeneName, - GeneStrand, 
                -SuperGroup, -AggregationGroup, -seqCount) %>%
  dplyr::select(ClinicalStudy, SubjectID, IS, dplyr::everything()) %>% tibble::as_tibble() %>%
  dplyr::inner_join(celltypes_markers)


input_sr = abundances %>%
  ## add counts of CD34
  dplyr::inner_join(
    old_stats %>%
      dplyr::select(ClinicalStudy, SubjectID, count_g1) %>%
      dplyr::rename(nIS_progenitor=count_g1) %>% unique()
  ) %>%
  ## remove patients to exclude
  dplyr::filter(!SubjectID %in% patients_to_exclude) %>%
  dplyr::select(IS, SubjectID, ClinicalStudy, Tissue, TimepointMonths,
                AggrSuperGroup, CellType, fragmentEstimate, nIS_progenitor, colorcode) %>%
  dplyr::rename(abundance=fragmentEstimate)


## only ISs from CD34
input_cd34 = input_sr %>% dplyr::filter(CellType == "CD34") %>%
  dplyr::group_by(IS, SubjectID, ClinicalStudy, AggrSuperGroup) %>%
  # sum abundance over timepoints
  dplyr::summarise(abundance=sum(abundance)) %>% dplyr::ungroup()
input_cd34 = input_cd34[!duplicated(input_cd34$IS), ]  # remove duplicated ISs
input_cd34_wide = input_cd34 %>%
  tidyr::pivot_wider(names_from="AggrSuperGroup", values_from="abundance", values_fill=0)
## ISs no CD34
input_nocd34 = input_sr %>% dplyr::filter(CellType != "CD34")


richness_cd34 = input_cd34 %>%
  tidyr::pivot_wider(names_from="AggrSuperGroup", values_from="abundance", values_fill=0) %>%
  dplyr::group_by(SubjectID, ClinicalStudy) %>%
  tidyr::nest() %>%
  dplyr::summarise(richness=compute_richness(data[[1]], subset_cols="CD34", gt=TRUE))



compute_gt_emp_ratios = function(abund, marker, tissue) {
  abund %>%
    dplyr::group_by(SubjectID, ClinicalStudy, TimepointMonths) %>%
    tidyr::nest() %>%
    dplyr::summarise(gt_shared=gt_shared_species(data[[1]], subset_cols=ct, ratio=F),
                     emp_shared=compute_ratio(assemblages_df=data[[1]],
                                              subset_cols=ct, ratio=F)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(AggrSuperGroup=marker, Tissue=tissue)
}


# Compute estimations
res_ab.withNAs = res_nIS = res_ab = tibble::tibble()
for (ct in unique(input_nocd34$AggrSuperGroup)) {
  for (tid in unique(input_nocd34$Tissue)) {
    m = input_nocd34 %>% dplyr::filter(AggrSuperGroup == ct, Tissue==tid) %>% nrow()
    if (m == 0) next

    subs_nocd34 = input_nocd34 %>%
      dplyr::select(-colorcode, -CellType) %>%
      dplyr::filter(AggrSuperGroup==ct, Tissue==tid) %>%
      dplyr::select(-Tissue) %>%
      tidyr::pivot_wider(names_from="AggrSuperGroup", values_from="abundance", values_fill=0)

    res_ab.withNAs = res_ab.withNAs %>% dplyr::bind_rows(
      subs_nocd34 %>%
        dplyr::left_join(input_cd34_wide) %>%
        tidyr::drop_na() %>%
        compute_gt_emp_ratios(marker=ct, tissue=tid)
    )
    
    res_ab = res_ab %>% dplyr::bind_rows(
      subs_nocd34 %>%
        dplyr::left_join(input_cd34_wide) %>%
        dplyr::mutate(CD34=replace(CD34, is.na(CD34), 0)) %>%
        compute_gt_emp_ratios(marker=ct, tissue=tid)
    )

    res_nIS = res_nIS %>% dplyr::bind_rows(
      subs_nocd34 %>%
        dplyr::rename(CD34=nIS_progenitor) %>%
        compute_gt_emp_ratios(marker=ct, tissue=tid)
    )
  }
}


ratio_ab.withNAs = res_ab.withNAs %>%
  dplyr::left_join(richness_cd34) %>%
  dplyr::mutate(gt_ratio=gt_shared / richness,
                emp_ratio=emp_shared / richness)


## ratio_ab.1 - IS Analytics as denominator
ratio_ab.1 = res_ab %>%
  dplyr::left_join(
    old_stats %>%
      dplyr::select(ClinicalStudy, SubjectID, count_g1) %>%
      dplyr::rename(richness=count_g1) %>% unique()) %>%
  dplyr::mutate(gt_ratio=gt_shared / richness,
                emp_ratio=emp_shared / richness)

## ratio_ab.2 - GT richness as denominator
ratio_ab.2 = res_ab %>%
  dplyr::left_join(richness_cd34) %>%
  dplyr::mutate(gt_ratio=gt_shared / richness,
                emp_ratio=emp_shared / richness)

ratio_nIS = res_nIS %>%
  dplyr::left_join(
      old_stats %>%
        dplyr::select(ClinicalStudy, SubjectID, g1, count_g1) %>%
        dplyr::rename(nIS_progenitor=count_g1) %>% unique()) %>%
    dplyr::mutate(gt_ratio=gt_shared / nIS_progenitor,
                  emp_ratio=emp_shared / nIS_progenitor)

ratio_all = ratio_ab.1 %>%
  dplyr::rename(nIS_progenitor=richness) %>%
  dplyr::mutate(type="GT/ISAnalytics") %>%
  dplyr::add_row(
    ratio_ab.2 %>%
      dplyr::rename(nIS_progenitor=richness) %>%
      dplyr::mutate(type="GT/GT")
  ) %>% 
  dplyr::add_row(
    ratio_nIS %>% dplyr::select(-g1) %>%
      dplyr::mutate(type="IsAnalytics/IsAnalytics")
    ) %>%

  dplyr::left_join(patients_infos %>%
                     dplyr::mutate(AgeGroup = dplyr::case_when(
                       age <= 2 ~"0-2", age <= 15 ~"2-15", age > 15 ~"30+"
  ))) %>%
  dplyr::left_join(celltypes_markers)


## Save dataframes
write.csv(ratio_all, paste0(save_path, "datasets/gtratio_all_cd34shared.csv"), row.names=FALSE)
write.csv(ratio_ab.2 %>% dplyr::select(-dplyr::starts_with("emp")) %>%
            dplyr::left_join(patients_infos %>%
                               dplyr::mutate(AgeGroup = dplyr::case_when(
                                 age <= 2 ~"0-2", age <= 15 ~"2-15", age > 15 ~"30+"
                               ))) %>%
            dplyr::left_join(celltypes_markers),
          paste0(save_path, "datasets/gtratio_gt:gt_cd34shared.csv"), row.names=FALSE)

# Read dataframes #####
ratio_all = read.csv(paste0(save_path, "datasets/gtratio_all_cd34shared.csv"))

# PLOTS ####
color_pal = ratio_all$colorcode %>% unique() %>%
  setNames(ratio_all$CellType %>% unique())


plot_main(ratio_all)
plot_main(ratio_ab)


# IS Analytics as denominator
fname = "gtratio_gt:gt_NAs"
plot_main(ratio_ab.withNAs)
ggsave(paste0(save_path, "plots/pdfs/", fname, "_cd34output.pdf"), height=6, width=12)
ggsave(paste0(save_path, "plots/", fname, "_cd34output.png"), height=6, width=12, dpi=600)

plot_main(ratio_ab.withNAs %>% dplyr::filter(ClinicalStudy!="XSCID"), by_age=TRUE)
ggsave(paste0(save_path, "plots/pdfs/", fname, "_byage_cd34output.pdf"), height=10, width=12)
ggsave(paste0(save_path, "plots/", fname, "_byage_cd34output.png"), height=10, width=12, dpi=600)


# IS Analytics as denominator
fname = "gtratio_gt:IS"
plot_main(ratio_ab.1)
ggsave(paste0(save_path, "plots/pdfs/", fname, "_cd34output.pdf"), height=6, width=12)
ggsave(paste0(save_path, "plots/", fname, "_cd34output.png"), height=6, width=12, dpi=600)

plot_main(ratio_ab.1 %>% dplyr::filter(ClinicalStudy!="XSCID"), by_age=TRUE)
ggsave(paste0(save_path, "plots/pdfs/", fname, "_byage_cd34output.pdf"), height=10, width=12)
ggsave(paste0(save_path, "plots/", fname, "_byage_cd34output.png"), height=10, width=12, dpi=600)


fname = "gtratio_gt:gt"
plot_main(ratio_ab.2)
ggsave(paste0(save_path, "plots/pdfs/", fname, "_cd34output.pdf"), height=6, width=12)
ggsave(paste0(save_path, "plots/", fname, "_cd34output.png"), height=6, width=12, dpi=600)

plot_main(ratio_ab.1 %>% dplyr::filter(ClinicalStudy!="XSCID"), by_age=TRUE)
ggsave(paste0(save_path, "plots/pdfs/", fname, "_byage_cd34output.pdf"), height=10, width=12)
ggsave(paste0(save_path, "plots/", fname, "_byage_cd34output.png"), height=10, width=12, dpi=600)


fname = "gtratio_IS:IS"
plot_main(ratio_nIS)
ggsave(paste0(save_path, "plots/pdfs/", fname, "_cd34output.pdf"), height=6, width=12)
ggsave(paste0(save_path, "plots/", fname, "_cd34output.png"), height=6, width=12, dpi=600)

plot_main(ratio_nIS %>% dplyr::filter(ClinicalStudy!="XSCID"), by_age=TRUE)
ggsave(paste0(save_path, "plots/pdfs/", fname, "_byage_cd34output.pdf"), height=10, width=12)
ggsave(paste0(save_path, "plots/", fname, "_byage_cd34output.png"), height=10, width=12, dpi=600)



## plots ISAnalytics ####

old_stats %>%
  dplyr::filter(!SubjectID %in% patients_to_exclude) %>%
  dplyr::left_join(patients_infos %>% dplyr::mutate(AgeGroup = dplyr::case_when(
    age <= 2 ~"0-2", age <= 15 ~"2-15", age > 15 ~"30+"
  ))) %>%
  dplyr::inner_join(celltypes_markers) %>%
  dplyr::filter(TimepointMonths <= 60, !AggrSuperGroup %in% c("CD34","Plasma")) %>%
  
  ggplot() +

  stat_smooth(aes(x=TimepointMonths, y=shared/count_g1 * 100,
                  color=CellType, fill=CellType),
              formula=y~log(x),
              se=T, level=.4) +

  scale_color_manual(values=color_pal) +
  scale_fill_manual(values=color_pal) +

  ggh4x::facet_nested(. ~ factor(ClinicalStudy, levels=c("MLD","WAS","BTHAL","XSCID")),
                      scales="free", space="free_x") +
  scale_x_continuous(breaks = seq(0, max(old_stats$TimepointMonths, na.rm = T), 12)) +
  theme_bw() +
  theme(strip.text=element_text(face="bold", size=16)) +
  theme(legend.direction="horizontal", legend.position="bottom", legend.box="horizontal") +
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16),
        axis.title=element_text(size=16), plot.title=element_text(size=20)) +
  labs(x="FollowUp months after gene therapy",
       y="Perc. of IS shared with CD34 BM",
       colour="CellType", fill="CellType")
ggsave(paste0(save_path, "plots/main_cd34output.pdf"), height=4, width=12)
ggsave(paste0(save_path, "plots/main_cd34output.png"), height=4, width=12, dpi=600)


old_stats %>%
  dplyr::filter(!SubjectID %in% patients_to_exclude,
                ClinicalStudy != "XSCID") %>%
  dplyr::left_join(patients_infos %>% dplyr::mutate(AgeGroup = dplyr::case_when(
    age <= 2 ~"0-2", age <= 15 ~"2-15", age > 15 ~"30+"
  ))) %>%
  dplyr::inner_join(celltypes_markers) %>%
  dplyr::filter(TimepointMonths <= 60, !AggrSuperGroup %in% c("CD34","Plasma")) %>%
  
  ggplot() +
  
  stat_smooth(aes(x=TimepointMonths, y=shared/count_g1 * 100,
                  color=CellType, fill=CellType),
              formula=y~log(x),
              se=T, level=.4) +
  
  scale_color_manual(values=color_pal) +
  scale_fill_manual(values=color_pal) +
  
  ggh4x::facet_nested(AgeGroup ~ factor(ClinicalStudy, levels=c("MLD","WAS","BTHAL","XSCID")),
                      scales="free", space="free_x") +
  scale_x_continuous(breaks = seq(0, max(old_stats$TimepointMonths, na.rm = T), 12)) +
  theme_bw() +
  theme(strip.text=element_text(face="bold", size=16)) +
  theme(legend.direction="horizontal", legend.position="bottom", legend.box="horizontal") +
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16),
        axis.title=element_text(size=16), plot.title=element_text(size=20)) +
  labs(x="FollowUp months after gene therapy",
       y="Perc. of IS shared with CD34 BM",
       colour="CellType", fill="CellType")
ggsave(paste0(save_path, "plots/main_byage_cd34output.pdf"), height=8, width=12)
ggsave(paste0(save_path, "plots/main_byage_cd34output.png"), height=8, width=12, dpi=600)



# Compare new-old ####
all = ratio_ab %>%
  dplyr::full_join(
    old_stats %>% tibble::rownames_to_column() %>%
      dplyr::filter(!AggrSuperGroup %in% c("CD34","Plasma"), !SubjectID %in% patients_to_exclude) %>%
      dplyr::select(SubjectID, ClinicalStudy, TimepointMonths, Tissue,
                    AggrSuperGroup, count_g1, shared)) %>%
  dplyr::left_join(patients_infos %>% dplyr::mutate(AgeGroup = dplyr::case_when(
    age <= 2 ~"0-2", age <= 15 ~"2-15", age > 15 ~"30+"
  ))) %>%
  dplyr::filter(!SubjectID %in% patients_to_exclude) %>%
  dplyr::inner_join(celltypes_markers)


all %>%
  dplyr::filter(ClinicalStudy == "BTHAL", AggrSuperGroup %in% c("CD36","GLY")) %>%
  dplyr::filter(TimepointMonths <= 60) %>%
  ggplot() +
  geom_point(aes(x=emp_shared/count_g1*100,
                 y=shared/count_g1*100,
                 color=AggrSuperGroup)) +
  ggh4x::facet_nested(TimepointMonths ~ Tissue + AgeGroup) +
  geom_abline() +
  theme_bw()


all %>%
  dplyr::filter(ClinicalStudy == "BTHAL", AggrSuperGroup %in% c("CD36","GLY"), AgeGroup == "30+") %>%
  ggplot() +
  geom_point(aes(x=gt_shared, y=shared, color=AggrSuperGroup)) +
  ggh4x::facet_nested(TimepointMonths ~ Tissue + AgeGroup) +
  # ggh4x::facet_nested(Tissue ~ TimepointMonths) +
  geom_abline() +
  theme_bw()


all %>%
  dplyr::filter(ClinicalStudy == "BTHAL", AggrSuperGroup %in% c("CD36","GLY"), AgeGroup == "30+") %>%
  ggplot() +
  geom_point(aes(x=count_g1, y=richness, color=AggrSuperGroup)) +
  ggh4x::facet_nested(TimepointMonths ~ Tissue + AgeGroup) +
  # ggh4x::facet_nested(TimepointMonths ~ Tissue) +
  geom_abline() +
  theme_bw()







### OLD 
# res_sr = lapply(cellmark, function(ct) {
#   lapply(tissues, function(tid) {
#     # filter tissue and marker
#     subs = input_nocd34 %>%
#       dplyr::select(-colorcode, -CellType) %>%
#       dplyr::filter(AggrSuperGroup == ct, Tissue==tid)
# 
#     if (nrow(subs) == 0) return()
#     
#     subs %>% dplyr::select(-Tissue, -nIS_progenitor) %>%
#       tidyr::pivot_wider(names_from="AggrSuperGroup", values_from="abundance", values_fill=0) %>%
# 
#       dplyr::left_join(input_cd34 %>%
#                          tidyr::pivot_wider(names_from="AggrSuperGroup",
#                                             values_from="abundance", values_fill=0)) %>%
#       # dplyr::rename(CD34=nIS_progenitor) %>%
# 
#       dplyr::group_by(SubjectID, ClinicalStudy, TimepointMonths) %>%
#       tidyr::nest() %>%
#       dplyr::summarise(gt_shared=gt_shared_species(data[[1]], subset_cols=ct, ratio=F),
#                        emp_shared=compute_ratio(assemblages_df=data[[1]],
#                                                 subset_cols=ct, ratio=F)) %>%
#       dplyr::ungroup() %>%
#       dplyr::mutate(AggrSuperGroup=ct, Tissue=tid)
#     }) %>% do.call(rbind, .)
#   }) %>% do.call(rbind, .) %>%
#   
#   dplyr::left_join(
#     old_stats %>%
#       dplyr::select(ClinicalStudy, SubjectID, g1, count_g1) %>%
#       dplyr::rename(nIS_progenitor=count_g1) %>% unique()) %>%
#   dplyr::mutate(gt_ratio=gt_shared / nIS_progenitor,
#                 emp_ratio=emp_shared / nIS_progenitor)


