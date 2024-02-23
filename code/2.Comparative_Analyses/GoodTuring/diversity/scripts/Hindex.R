main_path = "~/Dropbox/shared/HSC dynamics/"
setwd(main_path)
save_path = paste0(main_path, "diversity/")

## Load the functions to generate the input df
source(paste0(main_path, "r_functions/regression_helper_fns.R"))
## Load the sklearn model functions ####
reticulate::source_python(paste0(main_path, "r_functions/sklearn_pipe.py"))

require(tidyverse)
require(ggplot2)

# Import the dataset ####
div.all = openxlsx::read.xlsx(paste0(main_path, "diversity/data/202112.Diversity.ByMarker.xlsx")) %>%
  dplyr::filter(PCRMethod == "SLiM" | ClinicalTrial == "WAS")

# get the color palette
color_pal = div.all %>% dplyr::select(colorcode, CellType) %>% unique()
color_pal = color_pal$colorcode %>% setNames(color_pal$CellType)


#### H-Index ####
patients_infos = read.csv(paste0(main_path, "datasets/patients_age.csv"))

subs.tab = div.all %>%
  tibble::as_tibble() %>%
  dplyr::filter(!(Tissue=="PB" & CellType %in% c("CD34","Erythroid"))) %>%
  # dplyr::select(Tissue, Study, CellMarker, CellType, SubjectID, TimePoint, VCN_avg,
  #               ng.DNA.corrected_sum, NGSTechnology, Hindex, Hindex_normalized) %>%
  
  dplyr::left_join(patients_infos) %>% dplyr::filter(!is.na(infused_cells)) %>% 

  # to make the names uniform (both Novaseq and NovaSeq are present)
  dplyr::mutate(NGSTechnology=tolower(NGSTechnology),
                PCRMethod=tolower(PCRMethod)) %>%
  
  # fix NAs for VCN -> as avg of timepoint + tissue
  dplyr::group_by(TimePoint, Tissue) %>% 
  dplyr::mutate(VCN_avg=replace(VCN_avg, is.na(VCN_avg), mean(VCN_avg, na.rm=T))) %>% 
  dplyr::ungroup() %>% 
  # if still NAs -> avg for marker
  dplyr::group_by(CellMarker) %>% 
  dplyr::mutate(VCN_avg=replace(VCN_avg, is.na(VCN_avg), mean(VCN_avg, na.rm=T))) %>% 
  dplyr::ungroup() %>% 
  
  dplyr::group_by(SubjectID, ClinicalStudy, Tissue, TimePoint, CellType) %>%
  # dplyr::mutate(mean_ng=mean(ng.DNA.corrected_sum)) %>% unique() %>%
  dplyr::mutate(mean_ng=max(0.01, mean(ng.DNA.corrected_sum, na.rm=T))) %>%
  
  dplyr::group_by(SubjectID, Tissue, TimePoint) %>%
  dplyr::mutate(h_zscaled=as.numeric(scale(Hindex))) %>% 
  dplyr::ungroup() %>% 
  
  dplyr::rename(Timepoint=TimePoint) %>%
  
  dplyr::select(SubjectID, ClinicalStudy, Tissue, Timepoint, CellMarker, CellType,
                Hindex, mean_ng, infused_cells, age, nIS, 
                # Hindex_normalized, h_zscaled,
                VCN_avg, NGSTechnology, Gender, PCRMethod, colorcode)
  
  

# Create a dataframe with only the Tissue/Study/CellType and the rest of the dataset nested
# Column "data" contains a tibble with the data
nested.tab = subs.tab %>%
  dplyr::group_by(Tissue, ClinicalStudy, CellType) %>%
  tidyr::nest() %>%
  dplyr::ungroup()



# Fit - Hindex ####
cols_fit = list(
  c("Timepoint", "mean_ng", "nIS", "VCN_avg", "infused_cells", 
    "NGSTechnology", "PCRMethod", "Gender"),
  c("Timepoint", "mean_ng", "nIS", "VCN_avg", "infused_cells", 
    "NGSTechnology", "PCRMethod", "Gender", "age")
)

nested.fit = nested.tab %>%
  dplyr::group_by(Tissue, ClinicalStudy, CellType) %>%
  do_all_fits(y_col="Hindex", cols_fit=cols_fit) %>% 
  dplyr::mutate(y_col="Hindex")




## Unnest preds and save ####

unnested = nested.fit %>%
  dplyr::ungroup() %>%
  tidyr::pivot_longer(cols=tidyr::starts_with("fit"), names_to="fitname") %>%
  tidyr::unnest_wider(value) %>%
  tidyr::unnest(cols=c("fitname","CellMarker","SubjectID","Timepoint",
                       "CellType", "nIS", "mean_ng", "NGSTechnology",
                       "VCN_avg", "age", "infused_cells", "PCRMethod", "Gender",
                       # "Hindex_normalized", "h_zscaled", 
                       "Hindex", "fit", "se")) %>%
  dplyr::rename(pred_Hindex=fit, pred_Hindex_se=se) %>%
  dplyr::mutate(regression_formula=dplyr::case_when(
    fitname == "fit0" ~ paste0(y_col, " ~ Timepoint"),
    fitname == "fit1" ~ paste0(y_col, " ~ ", paste(cols_fit[[1]], collapse=" + ")),
    fitname == "fit2" ~ paste0(y_col, " ~ ", paste(cols_fit[[2]], collapse=" + "))
  )) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(SubjectID, Tissue, Timepoint) %>%
  dplyr::mutate(pred_h_zscaled=as.numeric(scale(pred_Hindex)),
                h_zscaled=as.numeric(scale(Hindex))) %>% 
  dplyr::ungroup()

## Save to RDS and csv
saveRDS(unnested %>% 
          dplyr::mutate(final_fit=fitname=="fit2"), 
        paste0(save_path, "datasets/PRED_202112.Diversity.ByMarker.Hindex.Rds"))

unnested %>%
  dplyr::select(-data) %>%
  dplyr::mutate(final_fit=fitname=="fit2") %>% 
  write.csv(paste0(save_path, "datasets/PRED_202112.Diversity.ByMarker.Hindex.csv"))

unnested %>%
  dplyr::select(-data) %>%
  dplyr::mutate(final_fit=fitname=="fit2") %>% 
  write.csv(paste0(main_path, "regression_datasets/PRED_202112.Diversity.ByMarker.Hindex.csv"))



## Comparison among fits ####

## Hindex
unnested %>%
  ggplot() +
  geom_smooth(aes(x=Timepoint, y=Hindex), fill="orange", color="orange",
              alpha=.3, level=.75, linetype="solid", linewidth=0.7) +
  
  geom_smooth(aes(x=Timepoint, y=pred_Hindex, fill=regression_formula, color=regression_formula),
              alpha=.3, level=.75, linetype="solid", linewidth=0.7) +
  
  ggh4x::facet_nested(Tissue + CellType ~ ClinicalStudy, scales="free_x") +
  theme_bw() + theme(legend.position="right") +
  labs(color="Regression model", fill="Regression model")


## Hindex scaled
unnested %>%

  ggplot() +
  geom_smooth(aes(x=Timepoint, y=h_zscaled), fill="orange", color="orange",
              alpha=.3, level=.75, linetype="solid", linewidth=0.7) +
  
  geom_smooth(aes(x=Timepoint, y=pred_h_zscaled,
                  fill=regression_formula, color=regression_formula),
              alpha=.3, level=.75, linetype="solid", linewidth=0.7) +
  
  ggh4x::facet_nested(Tissue + CellType ~ ClinicalStudy, scales="free_x") +
  theme_bw() + theme(legend.position="bottom") +
  labs(color="Regression model", fill="Regression model")







## Slice diversity data ####
library(ggpubr); library(rstatix)

study_list = c("MLD", "WAS", "BTHAL")
celltype_list = c("CD34", "Erythroid", "Myeloid", "B", "T")
faceting_order = c("MLD", "WAS", "BTHAL")
markerlist = c("CD13", "CD14", "CD15", "CD19", "CD3", "CD34", "GLY", "CD36", "GLYA")
trials_colors = c("darkblue", "forestgreen", "firebrick") %>% setNames(ClinicalStudy_list)

patient_h_scaled_summary_filtered = unnested %>% 
  dplyr::filter(CellMarker %in% markerlist, Tissue %in% c("PB", "BM"), 
                !is.na(pred_Hindex), Timepoint > 0) %>% 
  # dplyr::group_by(SubjectID, Tissue, Timepoint) %>% 
  # dplyr::mutate(h_normalized_zscaled=scale(Hindex_normalized)) %>% 
  # dplyr::filter(PCRMethod == "SLiM" | ClinicalStudy == "WAS") %>%   # already done
  dplyr::filter(Timepoint >= 24 & Timepoint <= 62) %>% 
  dplyr::group_by(SubjectID, Tissue, Timepoint, fitname) %>% 
  dplyr::mutate(n_obs_stt_overtime = dplyr::n()) %>% 
  dplyr::filter(n_obs_stt_overtime >= 3) %>% 
  
  dplyr::mutate(CellType=factor(CellType, levels=celltype_list),
                ClinicalStudy=factor(ClinicalStudy, levels=study_list)) %>% 
  
  dplyr::group_by(ClinicalStudy, SubjectID, Tissue, CellMarker, CellType, fitname) %>%
  ggpubr::get_summary_stats(pred_h_zscaled, type = "full") %>% 
  dplyr::filter(n >=2)


write.csv(patient_h_scaled_summary_filtered, 
          paste0(save_path, "datasets/pred.diversity.zscore.patient_h_scaled_summary_filtered.csv"))






# #### H-Index scaled
# subs.tab.scaled = subs.tab %>%
#   
#   dplyr::group_by(SubjectID, Tissue, Timepoint) %>%
#   dplyr::mutate(h_zscaled=as.numeric(scale(Hindex))) %>%
#   
#   dplyr::group_by(SubjectID, Tissue, Timepoint) %>%
#   dplyr::mutate(h_normalized_zscaled=as.numeric(scale(Hindex_normalized))) %>%
#   
#   dplyr::ungroup()
# 
# 
# nested.tab.scaled = subs.tab.scaled %>%
#   dplyr::group_by(Tissue, ClinicalStudy, CellType) %>%
#   tidyr::nest() %>%
#   dplyr::ungroup()

# Fit - Hindex Scaled
# y_cols = "h_zscaled"
# 
# nested.fit.h_zscaled = nested.tab.scaled %>%
#   tidyr::unnest(data) %>%
#   dplyr::filter(!is.na(h_zscaled) & !is.na(Timepoint) & !is.na(ng.DNA.corrected_sum) & !is.na(NGSTechnology)) %>%
# 
#   dplyr::group_by(Tissue, ClinicalStudy, CellType) %>%
#   tidyr::nest() %>%
#   dplyr::summarise(fit1 = list(fit_pipeline( get_x(data, to_py=T),
#                                              get_y(data, column=y_cols, to_py=T),
#                                              get_x(data, to_py=T),
#                                              subset=reticulate::tuple(c("Timepoint")))),
#                    fit2 = list(fit_pipeline( get_x(data, to_py=T),
#                                              get_y(data, column=y_cols, to_py=T),
#                                              get_x(data, to_py=T),
#                                              subset=c("Timepoint", "ng.DNA.corrected_sum"))),
#                    fit3 = list(fit_pipeline( get_x(data, to_py=T),
#                                              get_y(data, column=y_cols, to_py=T),
#                                              get_x(data, to_py=T),
#                                              subset=c("Timepoint", "ng.DNA.corrected_sum", "NGSTechnology"))),
#                    data = data)



## Unnest preds and save

# unnested.h_zscaled = nested.fit.h_zscaled %>%
#   dplyr::ungroup() %>%
#   tidyr::pivot_longer(cols=tidyr::starts_with("fit"), names_to="fitname") %>%
#   tidyr::unnest_wider(value) %>%
#   tidyr::unnest(cols=c("fitname","CellMarker","SubjectID","Timepoint",
#                        "ng.DNA.corrected_sum","NGSTechnology","Hindex",
#                        "Hindex_normalized","h_zscaled","h_normalized_zscaled",
#                        "fit", "se")) %>%
#   dplyr::rename(pred_h_zscaled=fit, pred_h_zscaled_se=se) %>%
#   dplyr::mutate(regression_formula=dplyr::case_when(
#     fitname == "fit1" ~ "h_zscaled ~ Timepoint",
#     fitname == "fit2" ~ "h_zscaled ~ Timepoint + ng.DNA.corrected_sum",
#     fitname == "fit3" ~ "h_zscaled ~ Timepoint + ng.DNA.corrected_sum + NGSTechnology"
#   ))
# 
# 
# ## Save RDS and csv
# saveRDS(unnested.h_zscaled, paste0(save_path, "predicted_h_zscaled.Rds"))
# 
# unnested.h_zscaled %>%
#   dplyr::select(-data) %>%
#   write.csv(paste0(save_path, "pred_202112.Diversity.ByMarker.h_zscaled.csv"))


