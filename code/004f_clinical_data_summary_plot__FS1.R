source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library("patchwork")

load(file = pff("patient_clinical_data.rsav"))
load(file = pff("expr_mat.rsav"))
patient_clinical_data$matching_RNAseq = patient_clinical_data$patient  %in% colnames(expr_mat)
patient_clinical_data$group = gsub("-CPC-GENE", "", patient_clinical_data$group)

# helper funtions to add median or bin counts
ggplot_count_func = function(x) 
		geom_text(stat = "count", aes(label = after_stat(count)), size = 3, y = 0, color = "purple")
ggplot_median_func_part1 = function(median_val) 
		geom_vline(xintercept = median_val, linetype = "dashed")
ggplot_median_func_part2 = function(median_val)
		geom_text(label = median_val, x = median_val, y = 0.5)
		
plot_theme2 = function() plot_theme() + theme(axis.title.x = element_blank())


median_age = median(patient_clinical_data$age_at_tumour_collection, na.rm = TRUE)
plt_age = ggplot(patient_clinical_data, aes(age_at_tumour_collection)) + 
		geom_histogram() + 
		ggplot_median_func_part1(median_age) + ggplot_median_func_part2(median_age) + 
		geom_vline(xintercept = 55) + 
		ggtitle (NULL, "Donor age") + 
		plot_theme2()

plt_country = ggplot(patient_clinical_data, aes(group)) +
		geom_bar(stat = "count") +
		ggtitle (NULL, "Country") + 
		ggplot_count_func() + 
		plot_theme2()
		
plt_RNAseq = ggplot(patient_clinical_data, aes(matching_RNAseq)) +
		geom_bar(stat = "count") +
		ggtitle (NULL, "Matching RNAseq") + 
		ggplot_count_func() + 
		plot_theme2()

# psa at tumour collection, cap at 100 for better visualisation
PSA_CAP_VALUE = 100
patient_clinical_data$psa_at_tumour_collection_cap100 = patient_clinical_data$psa_at_tumour_collection
patient_clinical_data$psa_at_tumour_collection_cap100 [patient_clinical_data$psa_at_tumour_collection_cap100 > PSA_CAP_VALUE] = PSA_CAP_VALUE
median_PSA = median(patient_clinical_data$psa_at_tumour_collection, na.rm = TRUE)

plt_psa = ggplot(patient_clinical_data, aes(psa_at_tumour_collection_cap100)) +
		geom_histogram() + 
		ggplot_median_func_part1(median_PSA) + ggplot_median_func_part2(median_PSA) + 
		ggtitle (NULL, paste("PSA at tumor collection, cap", PSA_CAP_VALUE)) + 
		plot_theme2()

# pathology T stage
plt_Tstage = ggplot(patient_clinical_data, aes(T_stage)) +
		geom_bar(stat = "count") + 
		ggtitle (NULL, "Pathology T-stage") + 
		ggplot_count_func() + 
		plot_theme2()

# gleason grade group
plt_ggg = ggplot(patient_clinical_data, aes(gleason_grade_group)) +
		geom_bar(stat = "count") +
		ggtitle (NULL, "Gleason Grade Group") + 
		ggplot_count_func() + 
		plot_theme2()

# ancestry, predicted
plt_ancestry = ggplot(patient_clinical_data, aes(predicted_ancestry)) + 
		geom_bar(stat = "count") +
		ggtitle (NULL, "Ancestry, predicted") + 
		ggplot_count_func() + 
		plot_theme2()

# 10 year relapse
plt_relapse_event10y = ggplot(patient_clinical_data, aes(relapse_event_10year)) + 
		geom_bar(stat = "count") +
		ggtitle (NULL, "Relapse status, 10 years") + 
		ggplot_count_func() + 
		plot_theme2()

# 10 year mets
plt_mets_event10y = ggplot(patient_clinical_data, aes(mets_event_10year)) + 
		geom_bar(stat = "count") +
		ggtitle (NULL, "Metastasis status, 10 years") + 
		ggplot_count_func() + 
		plot_theme2()
		
# germline risk genes
plt_germline_risk_genes = ggplot(patient_clinical_data, aes(germline_risk_gene)) + 
		geom_bar(stat = "count") +
		ggtitle (NULL, "Germline risk genes") + 
		ggplot_count_func() + 
		plot_theme2()
		
# percent genome altered
median_PGA = signif(median(patient_clinical_data$PGA, na.rm = TRUE), 2)
plt_pga = ggplot(patient_clinical_data, aes(PGA)) +
		geom_histogram() + 
		ggtitle (NULL, "Percent genome altered") + 
		ggplot_median_func_part1(median_PGA) + ggplot_median_func_part2(median_PGA) + 
		plot_theme2()
		
# SNV/indel count
median_GTMB_log1p = signif(median(patient_clinical_data$GTMB_log1p), 2)
plt_gtmb_log1p = ggplot(patient_clinical_data, aes(GTMB_log1p)) +
		geom_histogram() + 
		ggplot_median_func_part1(median_GTMB_log1p) + ggplot_median_func_part2(median_GTMB_log1p) + 
		ggtitle (NULL, "#SNVs/indels, log1p") + 
		plot_theme2()
		
# sequencing coverage, tumor and normal
load(pff("WGS_coverage_tumor_normal.rsav"))
TUM_DEP_CAP = 200
NORM_DEP_CAP = 100
WGS_coverage_tumor_normal$tumor_total_depth [WGS_coverage_tumor_normal$tumor_total_depth > TUM_DEP_CAP] = TUM_DEP_CAP
WGS_coverage_tumor_normal$normal_total_depth [WGS_coverage_tumor_normal$normal_total_depth > NORM_DEP_CAP] = NORM_DEP_CAP
median_depth_tum = signif(median(WGS_coverage_tumor_normal$tumor_total_depth), 2)
median_depth_norm = signif(median(WGS_coverage_tumor_normal$normal_total_depth), 2)

plt_coverage_tumor = ggplot(WGS_coverage_tumor_normal, aes(tumor_total_depth)) +
		geom_histogram() + 
		ggplot_median_func_part1(median_depth_tum) + ggplot_median_func_part2(median_depth_tum) + 
		ggtitle (NULL, paste("tumor total depth cap", TUM_DEP_CAP)) + 
		plot_theme2()

plt_coverage_normal = ggplot(WGS_coverage_tumor_normal, aes(normal_total_depth)) +
		geom_histogram() + 
		ggplot_median_func_part1(median_depth_norm) + ggplot_median_func_part2(median_depth_norm) + 
		ggtitle (NULL, paste("normal total depth cap", NORM_DEP_CAP)) + 
		plot_theme2()

# combine all plots
plt_combined = 
		(plt_country + plt_age + plt_psa + plt_Tstage + plt_ggg + plt_ancestry + plt_germline_risk_genes + 
		plt_mets_event10y + plt_relapse_event10y + plt_pga + plt_gtmb_log1p + 
		plt_RNAseq + plt_coverage_tumor + plt_coverage_normal) +  plot_layout(ncol = 3)

fname = pff("figures/Combined_clinical_summaries_TableS1.pdf")
ggsave(plt_combined, file = fname, width = 9, height = 12)
file_open_call2(fname)
