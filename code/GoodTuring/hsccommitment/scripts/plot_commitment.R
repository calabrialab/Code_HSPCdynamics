
# global_allstudies_profile_skewing_fulldata_notNA_stats <- read the excel file

global_allstudies_profile_skewing_fulldata_notNA_stats = readxl::read_xlsx(
  "MLD_WAS_BTHAL.AllPatients.hlfu_34MyBTEry.global_allstudies_profile_skewing_fulldata_notNA_stats.xlsx"
)
plot_CRonR_allFU_nopoints_global <- 
  ggplot(global_allstudies_profile_skewing_fulldata_notNA_stats, 
         aes(x = FollowUp, y = CRonR, fill = Class, color = Class), na.rm = T, se = TRUE) +
  geom_smooth(method = "loess", 
              formula = y ~ log(x),
              stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
  # scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
  # scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
  scale_x_continuous(breaks = seq(0, max(global_allstudies_profile_skewing_fulldata_notNA_stats$FollowUp, na.rm = T), 12) ) +
  facet_grid(Tissue ~ factor(ClinicalStudy, levels=c("MLD","WAS","BTHAL")), 
             scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=16)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
  labs(title = paste0("Recaptured IS in percentage"), 
       # subtitle = paste0("CR on R (C:Class, S:Singletons, R:Recaptured, U:All).", " Abundance threshold = ", (threshold_abundance), ". All patients included." ),
       x = "Months after gene therapy", y = "% of IS", colour = "Class", fill = "Class")

pdf(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonR.nopoints.v2.pdf"), height=8, width=10)
print(plot_CRonR_allFU_nopoints_global)
dev.off()
png(file = paste0("analyses/multi_uni_sharing/202112/", paste(study_list, collapse = "_"), ".AllPatients.hlfu_34MyBTEry.allFU.CRonR.nopoints.v2.png"), height=8, width=10, units = "in", res = 300)
print(plot_CRonR_allFU_nopoints_global)
dev.off()
