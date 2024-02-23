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


full.data = read.csv(paste0(main_path, "cd34output/datasets/gtratio_gt:gt_cd34shared.csv")) %>%
  dplyr::left_join(patients_infos) %>% 
  tibble::as_tibble() %>% 
  dplyr::rename(Timepoint=TimepointMonths)

# color_pal = full.data$colorcode %>% unique() %>% setNames(unique(full.data$ClassName))

div.all = openxlsx::read.xlsx(paste0(main_path, "diversity/data/202112.Diversity.ByMarker.xlsx")) %>%
  tibble::as_tibble() %>% 
  dplyr::rename(ClinicalStudy=ClinicalTrial)


ng_scaled = div.all %>%
  tibble::as_tibble() %>%

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

  dplyr::group_by(SubjectID, ClinicalStudy, Tissue, TimePoint, CellMarker) %>%
  dplyr::mutate(mean_ng=mean(ng.DNA.corrected_sum)) %>% unique() %>%

  dplyr::rename(Timepoint=TimePoint) %>%
  dplyr::ungroup() %>%

  dplyr::select(SubjectID, ClinicalStudy, Tissue, Timepoint, CellMarker,
                mean_ng, infused_cells, age, nIS, VCN_avg,
                NGSTechnology, Gender, PCRMethod, colorcode)
  

subs = full.data %>% dplyr::rename(CellMarker=AggrSuperGroup) %>% 
  dplyr::left_join(ng_scaled)



# Fit the models

cols_fit = list(
  c("Timepoint", "mean_ng", "nIS", "VCN_avg", "infused_cells", 
    "NGSTechnology", "PCRMethod", "Gender"),
  c("Timepoint", "mean_ng", "nIS", "VCN_avg", "infused_cells", 
    "NGSTechnology", "PCRMethod", "Gender", "age")
)


## Regression ####
nested.fit = subs %>%
  # tidyr::drop_na() %>%
  dplyr::group_by(Tissue, ClinicalStudy, CellMarker) %>%
  tidyr::nest() %>% 
  do_all_fits(y_col="gt_ratio", cols_fit=cols_fit) %>% 
  dplyr::mutate(y_col="gt_ratio")


## Unnest preds and save ####
unnested = nested.fit %>%
  tidyr::pivot_longer(cols=tidyr::starts_with("fit"), names_to="fitname") %>%
  tidyr::unnest_wider(value) %>%
  tidyr::unnest(cols=c("fitname","SubjectID","Timepoint",
                       "gt_shared", "richness", "AgeGroup", "CellType",
                       "gt_ratio", "nIS", "mean_ng", "NGSTechnology",
                       "VCN_avg", "age", "infused_cells", "PCRMethod", "Gender",
                       "colorcode", "fit", "se")) %>%
  
  dplyr::rename(pred=fit, pred_se=se) %>%
  dplyr::mutate(regression_formula=dplyr::case_when(
    fitname == "fit0" ~ paste0(y_col, " ~ TimePoint"),
    fitname == "fit1" ~ paste0(y_col, " ~ ", paste(cols_fit[[1]], collapse=" + ")),
    fitname == "fit2" ~ paste0(y_col, " ~ ", paste(cols_fit[[2]], collapse=" + "))
  )) %>% 
  dplyr::filter(!(ClinicalStudy=="XSCID" & fitname != "fit0")) %>% 
  
  dplyr::mutate(final_noage=(ClinicalStudy=="XSCID"&fitname=="fit0")|
                  (ClinicalStudy!="XSCID"&fitname=="fit1"),
                final_age=(ClinicalStudy=="XSCID"&fitname=="fit0")|
                  (ClinicalStudy!="XSCID"&fitname=="fit2")) %>% 
  
  dplyr::left_join(div.all %>% dplyr::select(CellMarker, CellType) %>% unique())


## Regression by age #####
nested.fit.byage = subs %>%
  dplyr::group_by(Tissue, ClinicalStudy, CellMarker, AgeGroup) %>%
  tidyr::nest() %>% 
  dplyr::filter(!is.na(AgeGroup)) %>% 
  do_all_fits(y_col="gt_ratio", cols_fit=cols_fit) %>% 
  dplyr::mutate(y_col="gt_ratio")


## Unnest preds and save ####
unnested.byage = nested.fit.byage %>%
  tidyr::pivot_longer(cols=tidyr::starts_with("fit"), names_to="fitname") %>%
  tidyr::unnest_wider(value) %>%
  tidyr::unnest(cols=c("fitname","SubjectID","Timepoint",
                       "gt_shared", "richness", "CellType",
                       "gt_ratio", "nIS", "mean_ng", "NGSTechnology",
                       "VCN_avg", "age", "infused_cells", "PCRMethod", 
                       "Gender", "colorcode", "fit", "se")) %>%
  
  dplyr::rename(pred=fit, pred_se=se) %>%
  dplyr::mutate(regression_formula=dplyr::case_when(
    fitname == "fit0" ~ paste0(y_col, " ~ TimePoint"),
    fitname == "fit1" ~ paste0(y_col, " ~ ", paste(cols_fit[[1]], collapse=" + ")),
    fitname == "fit2" ~ paste0(y_col, " ~ ", paste(cols_fit[[2]], collapse=" + "))
  )) %>% 
  dplyr::filter(!(ClinicalStudy=="XSCID" & fitname != "fit0")) %>% 
  
  dplyr::mutate(final_noage=(ClinicalStudy=="XSCID"&fitname=="fit0")|
                  (ClinicalStudy!="XSCID"&fitname=="fit1"),
                final_age=(ClinicalStudy=="XSCID"&fitname=="fit0")|
                  (ClinicalStudy!="XSCID"&fitname=="fit2")) %>% 
  
  dplyr::left_join(div.all %>% dplyr::select(CellMarker, CellType) %>% unique())



joined_fit = rbind(unnested %>% dplyr::mutate(fittype="noAge"), 
                   unnested.byage %>% dplyr::mutate(fittype="ageGroup"))


## Save RDS and csv
saveRDS(joined_fit, paste0(main_path, "cd34output/datasets/PRED_gtratio_cd34shared.Rds"))

joined_fit %>%
  dplyr::select(-data) %>%
  write.csv(paste0(main_path, "cd34output/datasets/PRED_gtratio_cd34shared.csv"),
            row.names=F)

joined_fit %>%
  dplyr::select(-data) %>%
  write.csv(paste0(main_path, "regression_datasets/PRED_gtratio_cd34shared.csv"),
            row.names=F)



# Plot #####
joined_fit = readRDS(paste0(main_path, "cd34output/datasets/PRED_gtratio_cd34shared.Rds"))

fittype = "noAge"

joined_fit %>%
  dplyr::filter(fittype == !!fittype) %>% 
  ggplot() +
  geom_smooth(data=subs %>% 
                dplyr::left_join(div.all %>% dplyr::select(CellMarker, CellType) %>% unique()),
              aes(x=Timepoint, y=gt_ratio), fill="orange", color="orange",
              alpha=.3, level=.75, linetype="solid", linewidth=0.7) +
  
  geom_smooth(aes(x=Timepoint, y=pred,
                  color=regression_formula,
                  fill=regression_formula),
              se=TRUE, level=.75) +
  
  ggh4x::facet_nested(Tissue + CellType ~ ClinicalStudy, scales="free") +
  theme_bw() + theme(legend.position="right") +
  labs(color="Regression model", fill="Regression model")

ggsave(paste0(main_path, "cd34output/plots/PRED_compare_fits_gtratio_cd34output.pdf"), height=12, width=12)
ggsave(paste0(main_path, "cd34output/plots/PRED_compare_fits_gtratio_cd34output.png"), height=12, width=12, dpi=600)



color_pal = joined_fit$colorcode %>% unique() %>%
  setNames(joined_fit$CellType %>% unique())


regr_formula_noage = joined_fit %>% 
  dplyr::filter(fittype==!!fittype, final_noage) %>% 
  dplyr::pull(regression_formula) %>% unique()
regr_formula_age = joined_fit %>% 
  dplyr::filter(fittype==!!fittype, final_age) %>% 
  dplyr::pull(regression_formula) %>% unique()

joined_fit %>%
  dplyr::filter(fittype==!!fittype, Timepoint <= 60) %>% 
  dplyr::left_join(div.all %>% dplyr::select(CellMarker, CellType) %>% unique()) %>% 
  dplyr::filter(final_age) %>%
  
  ggplot() +
  geom_smooth(aes(x=Timepoint, y=pred*100, fill=CellType, color=CellType),
              formula=y~log(x), method="loess", se=T, level=.4) +
  
  scale_color_manual(values=color_pal) +
  scale_fill_manual(values=color_pal) +
  
  scale_x_continuous(breaks = seq(0, max(joined_fit$Timepoint, na.rm=T), 12) ) +
  facet_grid(. ~ factor(ClinicalStudy, levels=c("MLD","WAS","BTHAL","XSCID")), 
             scales="free_x", space="free_x") +
  theme_bw() +
  theme(strip.text=element_text(face="bold", size=16)) +
  theme(legend.direction="horizontal", legend.position="bottom", 
        legend.box="horizontal") +
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16), 
        axis.title=element_text(size=16), plot.title=element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), 
       subtitle = paste0("Regression: ", regr_formula_age),
       x="Months after gene therapy", y="% of IS", 
       colour="Class", fill="Class")

ggsave(paste0(main_path, "cd34output/plots/pdfs/PRED_age_gtratio_cd34output.pdf"), height=5, width=12)
ggsave(paste0(main_path, "cd34output/plots/PRED_age_gtratio_cd34output.png"), height=5, width=12, dpi=600)




joined_fit %>%
  dplyr::filter(fittype==!!fittype, Timepoint <= 60) %>% 
  dplyr::left_join(div.all %>% dplyr::select(CellMarker, CellType) %>% unique()) %>% 
  dplyr::filter(final_age, ClinicalStudy!="XSCID") %>%
  
  ggplot() +
  geom_smooth(aes(x=Timepoint, y=pred*100, fill=CellType, color=CellType),
              formula=y~log(x), method="loess", se=T, level=.4) +
  
  scale_color_manual(values=color_pal) +
  scale_fill_manual(values=color_pal) +
  
  scale_x_continuous(breaks = seq(0, max(joined_fit$Timepoint, na.rm=T), 12) ) +
  facet_grid(AgeGroup ~ factor(ClinicalStudy, levels=c("MLD","WAS","BTHAL")), 
             scales="free_x", space="free_x") +
  theme_bw() +
  theme(strip.text=element_text(face="bold", size=16)) +
  theme(legend.direction="horizontal", legend.position="bottom", 
        legend.box="horizontal") +
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16), 
        axis.title=element_text(size=16), plot.title=element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), 
       subtitle = paste0("Regression: ", regr_formula_age),
       x="Months after gene therapy", y="% of IS", 
       colour="Class", fill="Class")

ggsave(paste0(main_path, "cd34output/plots/pdfs/PRED_age_gtratio_byage_cd34output.pdf"), height=8, width=12)
ggsave(paste0(main_path, "cd34output/plots/PRED_age_gtratio_byage_cd34output.png"), height=8, width=12, dpi=600)


