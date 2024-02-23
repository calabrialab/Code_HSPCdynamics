main_path = "~/Dropbox/shared/HSC dynamics/"
setwd(main_path)

## Load the functions to generate the input df
source(paste0(main_path, "r_functions/regression_helper_fns.R"))
## Load the sklearn model functions
reticulate::source_python(paste0(main_path, "r_functions/sklearn_pipe.py"))

library(tidyverse)
require(ggplot2)

# Create the dataset and exploratory plots ####

patients_infos = read.csv(paste0(main_path, "datasets/patients_age.csv"))

full.data = read.csv(paste0(main_path, "hsccommitment/datasets/adj_ratios_commitment_wXSCID.csv"), row.names=1) %>%
  dplyr::left_join(patients_infos) %>%
  tibble::as_tibble() %>% 
  dplyr::rename(Timepoint=TimepointMonths)

color_pal = full.data$colorcode %>% unique() %>% setNames(unique(full.data$ClassName))

div.all = openxlsx::read.xlsx(paste0(main_path, "diversity/data/202112.Diversity.ByMarker.xlsx")) %>%
  tibble::as_tibble() %>% 
  dplyr::rename(ClinicalStudy=ClinicalTrial) %>% 
  dplyr::mutate(NGSTechnology=tolower(NGSTechnology),
                PCRMethod=tolower(PCRMethod)) %>% 
  dplyr::rename(Timepoint=TimePoint)

ng_scaled = div.all %>%
  dplyr::left_join(patients_infos) %>% 
  dplyr::filter(CellMarker %in% c("CD34", "CD13", "CD14", "CD15", "GLY",
                                  "CD36", "CD19", "CD3", "CD4", "CD8")) %>%
  dplyr::select(SubjectID, ClinicalStudy, Tissue, Timepoint, CellType,
                ng.DNA.corrected_sum, infused_cells, age, NGSTechnology, 
                Gender, PCRMethod) %>% 

  dplyr::group_by(SubjectID, ClinicalStudy, Tissue, Timepoint, CellType) %>%
  dplyr::mutate(mean_ng=mean(ng.DNA.corrected_sum)) %>% unique() %>% 
  dplyr::ungroup() %>% 
  
  dplyr::filter(CellType != "CD34") %>% 
  dplyr::mutate(ClassName=dplyr::case_when(
    CellType == "T" ~ "UniLyT",
    CellType == "B" ~ "UniLyB",
    CellType == "Myeloid" ~ "UniMyeloid",
    CellType == "Erythroid" ~ "UniErythroid"
  )) %>% dplyr::select(-CellType, -ng.DNA.corrected_sum) %>% 
  dplyr::group_by(SubjectID, Tissue, Timepoint, ClassName) %>% 
  dplyr::mutate(NGSTechnology=paste(unique(NGSTechnology),collapse="|"),
                PCRMethod=paste(unique(PCRMethod),collapse="|")) %>% 
  dplyr::ungroup() %>% unique()


ng_scaled2 = div.all %>%
  dplyr::left_join(patients_infos) %>% 
  dplyr::filter(CellMarker %in% c("CD34", "CD13", "CD14", "CD15", "GLY",
                                  "CD36", "CD19", "CD3", "CD4", "CD8")) %>%
  dplyr::select(SubjectID, ClinicalStudy, Tissue, Timepoint,
                ng.DNA.corrected_sum, infused_cells, age, NGSTechnology, 
                Gender, PCRMethod) %>% unique() %>% 

  dplyr::group_by(SubjectID, ClinicalStudy, Tissue, Timepoint) %>%
  dplyr::mutate(mean_ng=mean(ng.DNA.corrected_sum)) %>%
  dplyr::mutate(ClassName="Multilineage") %>%
  dplyr::ungroup() %>% 
  dplyr::select(-ng.DNA.corrected_sum) %>% 
  dplyr::group_by(SubjectID, Tissue, Timepoint, ClassName) %>% 
  dplyr::mutate(NGSTechnology=paste(unique(NGSTechnology),collapse="|"),
                PCRMethod=paste(unique(PCRMethod),collapse="|")) %>% 
  dplyr::ungroup() %>% unique()



nISs_VCN = readxl::read_xlsx(paste0(main_path, "hsccommitment/data/MLD_WAS_BTHAL.AllPatients.hlfu_34MyBTEry.global_allstudies_profile_skewing_fulldata_notNA_stats.xlsx")) %>%
  # dplyr::left_join(patients_infos) %>%
  dplyr::select(SubjectID, ClinicalStudy, Tissue, Timepoint,
                ClassName, CRonR, OverallNIS, VCN)

subs = full.data %>%
  dplyr::left_join(ng_scaled %>% dplyr::add_row(ng_scaled2)) %>%
  dplyr::left_join(nISs_VCN) %>% 
  dplyr::select(ClinicalStudy, SubjectID, ClassName, Timepoint, Tissue,
                dplyr::contains("ratio"), colorcode, age, infused_cells,
                dplyr::contains("_ng"), OverallNIS, VCN, NGSTechnology, 
                PCRMethod, Gender) %>% 
  dplyr::mutate(AgeGroup=dplyr::case_when(
    age <= 2 ~"0-2", age <= 15 ~"2-15", age > 15 ~"30+"))


# Fit the model ####

cols_fit = list(
  c("Timepoint", "mean_ng", "nIS", "VCN_avg", "infused_cells", 
    "NGSTechnology", "PCRMethod", "Gender"),
  c("Timepoint", "mean_ng", "nIS", "VCN_avg", "infused_cells", 
    "NGSTechnology", "PCRMethod", "Gender", "age")
)


## nested df
nested.fit = subs %>%
  # tidyr::drop_na() %>%
  dplyr::group_by(Tissue, ClinicalStudy, ClassName) %>%
  tidyr::nest() %>% 
  do_all_fits(y_col="ratio_chao", cols_fit=cols_fit) %>% 
  dplyr::mutate(y_col="ratio_chao")


## Unnest preds and save ####
unnested = nested.fit %>%
  dplyr::ungroup() %>%
  tidyr::pivot_longer(cols=tidyr::starts_with("fit"), names_to="fitname") %>%
  tidyr::unnest_wider(value) %>%
  tidyr::unnest(cols=c("fitname","SubjectID","Timepoint",
                       "mean_ng", "ratio_chao", "ratio_emp", "OverallNIS", "AgeGroup",
                       "VCN", "age", "infused_cells", "NGSTechnology", "PCRMethod", 
                       "Gender", "colorcode", "fit", "se")) %>%

  dplyr::rename(pred_chao=fit, pred_chao_se=se) %>%
  dplyr::mutate(regression_formula=dplyr::case_when(
    fitname == "fit0" ~ paste0(y_col, " ~ TimePoint"),
    fitname == "fit1" ~ paste0(y_col, " ~ ", paste(cols_fit[[1]], collapse=" + ")),
    fitname == "fit2" ~ paste0(y_col, " ~ ", paste(cols_fit[[2]], collapse=" + "))
  )) %>% 
  # dplyr::mutate(final_fit=(fitname=="fit2"))
  dplyr::mutate(final_noage=(ClinicalStudy=="XSCID"&fitname=="fit0")|
                  (ClinicalStudy!="XSCID"&fitname=="fit1"),
                final_age=(ClinicalStudy=="XSCID"&fitname=="fit0")|
                  (ClinicalStudy!="XSCID"&fitname=="fit2"))


## nested df
nested.fit.byage = subs %>%
  # tidyr::drop_na() %>%
  dplyr::group_by(Tissue, ClinicalStudy, ClassName, AgeGroup) %>%
  tidyr::nest() %>% 
  dplyr::filter(!is.na(AgeGroup)) %>% 
  do_all_fits(y_col="ratio_chao", cols_fit=cols_fit) %>% 
  dplyr::mutate(y_col="ratio_chao")


## Unnest preds and save ####
unnested.byage = nested.fit.byage %>%
  dplyr::ungroup() %>%
  tidyr::pivot_longer(cols=tidyr::starts_with("fit"), names_to="fitname") %>%
  tidyr::unnest_wider(value) %>%
  tidyr::unnest(cols=c("fitname","SubjectID","Timepoint",
                       "mean_ng", "ratio_chao", "ratio_emp", "OverallNIS", "AgeGroup",
                       "VCN", "age", "infused_cells", "NGSTechnology", "PCRMethod", 
                       "Gender", "colorcode", "fit", "se")) %>%
  
  dplyr::rename(pred_chao=fit, pred_chao_se=se) %>%
  dplyr::mutate(regression_formula=dplyr::case_when(
    fitname == "fit0" ~ paste0(y_col, " ~ TimePoint"),
    fitname == "fit1" ~ paste0(y_col, " ~ ", paste(cols_fit[[1]], collapse=" + ")),
    fitname == "fit2" ~ paste0(y_col, " ~ ", paste(cols_fit[[2]], collapse=" + "))
  )) %>% 
  # dplyr::mutate(final_fit=(fitname=="fit2"))
  dplyr::mutate(final_noage=(ClinicalStudy=="XSCID"&fitname=="fit0")|
                  (ClinicalStudy!="XSCID"&fitname=="fit1"),
                final_age=(ClinicalStudy=="XSCID"&fitname=="fit0")|
                  (ClinicalStudy!="XSCID"&fitname=="fit2"))


joined_fit = rbind(unnested %>% dplyr::mutate(fittype="noAge"), 
                   unnested.byage %>% dplyr::mutate(fittype="ageGroup"))

## Save RDS and csv
saveRDS(joined_fit, paste0(main_path, "hsccommitment/datasets/PRED_commitment_wXSCID.Rds"))

joined_fit %>%
  dplyr::select(-data) %>%
  write.csv(paste0(main_path, "hsccommitment/datasets/PRED_commitment_wXSCID.csv"),
            row.names=F)

joined_fit %>%
  dplyr::select(-data) %>%
  write.csv(paste0(main_path, "regression_datasets/PRED_commitment_wXSCID.csv"),
            row.names=F)





# Plot #####
joined_fit = readRDS(paste0(main_path, "hsccommitment/datasets/PRED_commitment_wXSCID.Rds"))

fittype = "noAge"

joined_fit %>%
  dplyr::filter(fittype == !!fittype) %>% 
  
  dplyr::group_by(Timepoint, Tissue, ClinicalStudy, ClassName, fitname) %>%

  ggplot() +
  geom_smooth(data=subs,
              aes(x=Timepoint, y=ratio_chao), fill="orange", color="orange",
              alpha=.3, level=.75, linetype="solid", linewidth=0.7) +

  geom_smooth(aes(x=Timepoint, y=pred_chao,
                  color=regression_formula,
                  fill=regression_formula),
              se=TRUE, level=.75) +

  ggh4x::facet_nested(Tissue + ClassName ~ ClinicalStudy, scales="free") +
  theme_bw() + theme(legend.position="right") +
  labs(color="Regression model", fill="Regression model")

ggsave(paste0(main_path, "hsccommitment/plots/PRED_gt_compare_fits_perc_ISs.pdf"), height=12, width=12)
ggsave(paste0(main_path, "hsccommitment/plots/PRED_gt_compare_fits_perc_ISs.png"), height=12, width=12, dpi=600)



regr_formula = joined_fit %>% dplyr::filter(final_age) %>% 
  dplyr::pull(regression_formula) %>% unique()

joined_fit %>%
  dplyr::filter(fittype == !!fittype, final_noage) %>% 
  dplyr::group_by(Timepoint, Tissue, ClinicalStudy, ClassName, fitname) %>%

  ggplot() +
  geom_smooth(aes(x=Timepoint, y=pred_chao, fill=ClassName, color=ClassName),
              se=TRUE, level=.75, alpha=0.4, formula=y~log(x), 
              stat="smooth", position="identity") +
  
  # geom_errorbar(aes(x=Timepoint, ymin=pred_chao-pred_chao_se, 
  #                   ymax=pred_chao+pred_chao_se, color=ClassName)) +
  
  scale_color_manual(values=color_pal) +
  scale_fill_manual(values=color_pal) +
  
  scale_x_continuous(breaks = seq(0, max(joined_fit$Timepoint, na.rm=T), 12) ) +
  facet_grid(Tissue ~ factor(ClinicalStudy, levels=c("MLD","WAS","BTHAL","XSCID")), 
             scales="free_x", space="free_x") +
  theme_bw() +
  theme(strip.text=element_text(face="bold", size=16)) +
  theme(legend.direction="horizontal", legend.position="bottom", 
        legend.box="horizontal") +
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16), 
        axis.title=element_text(size=16), plot.title=element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), 
       subtitle = paste0("Regression: ", regr_formula),
       x="Months after gene therapy", y="% of IS", 
       colour="Class", fill="Class")

ggsave(paste0(main_path, "hsccommitment/plots/PRED_commitment_ISs.pdf"), height=8, width=12)
ggsave(paste0(main_path, "hsccommitment/plots/PRED_commitment_ISs.png"), height=8, width=12, dpi=600)


joined_fit %>%
  dplyr::filter(fittype == !!fittype, Timepoint <= 60, final_noage,
                ClinicalStudy!="XSCID") %>% 
  dplyr::group_by(Timepoint, Tissue, ClinicalStudy, ClassName, fitname) %>%

  ggplot() +
  geom_smooth(aes(x=Timepoint, y=pred_chao, fill=ClassName, color=ClassName),
              se=TRUE, level=.75, alpha=0.4, formula=y~log(x), 
              stat="smooth", position="identity") +
  
  scale_color_manual(values=color_pal) +
  scale_fill_manual(values=color_pal) +
  
  scale_x_continuous(breaks = seq(0, max(joined_fit$Timepoint, na.rm=T), 12) ) +
  ggh4x::facet_nested(Tissue + AgeGroup ~ factor(ClinicalStudy, levels=c("MLD","WAS","BTHAL")), 
             scales="free_x", space="free_x") +
  theme_bw() +
  theme(strip.text=element_text(face="bold", size=16)) +
  theme(legend.direction="horizontal", legend.position="bottom", 
        legend.box="horizontal") +
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16), 
        axis.title=element_text(size=16), plot.title=element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), 
       subtitle = paste0("Regression: ", regr_formula),
       x="Months after gene therapy", y="% of IS", 
       colour="Class", fill="Class")

ggsave(paste0(main_path, "hsccommitment/plots/PRED_commitment_ISs_byage.pdf"), height=12, width=12)
ggsave(paste0(main_path, "hsccommitment/plots/PRED_commitment_ISs_byage.png"), height=12, width=12, dpi=600)




