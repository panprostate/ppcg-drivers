source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
source(paste0("../bin/", date_tag, "/004xx_functions_for_clinical_analyses.R"))
library(patchwork)

load(pff("all_patients.rsav"))
MIN_DRIVER_FREQ = 0.01
MIN_MUT_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)
MIN_MUT_PATIENTS

load(file = pff("patient_clinical_data.rsav"))
load(file = pff("patient_sets_for_drivers.rsav"))
patient_sets_for_drivers = patient_sets_for_drivers[sapply(patient_sets_for_drivers, length) >= MIN_MUT_PATIENTS]

# compute clinical stats with age stage grade as covariates
fname = pff(paste0("figures/KMplots_drivers_ALL__", MFS_ANALYSIS_TAG, ".pdf"))
pdf(fname, width = 7, height = 7)
cat("todo:", length(patient_sets_for_drivers), "\n")
driver_clinical_stats = do.call(rbind, 
		lapply(1:length(patient_sets_for_drivers), 
				clinical_analysis, patient_sets_for_drivers, patient_clinical_data, do_plot = TRUE, MIN_MUT_PATIENTS, MFS_COVARIATES, 
				precomp_clin_stats = NULL))
dev.off()
file_open_call2(fname)

driver_clinical_stats$analysis = MFS_ANALYSIS_TAG
save(driver_clinical_stats, file = pff(c("driver_clinical_stats__", MFS_ANALYSIS_TAG, ".rsav")))



# compute clinical stats with no covariates, these are used to preselect variables for the multivariate model
NULL_COVARIATES = c()
NULL_ANALYSIS_TAG = "PPCG__No_covar"

fname = pff(paste0("figures/KMplots_drivers_ALL__", NULL_ANALYSIS_TAG, ".pdf"))
pdf(fname, width = 7, height = 7)

cat("todo:", length(patient_sets_for_drivers), "\n")
driver_clinical_stats = do.call(rbind, 
		lapply(1:length(patient_sets_for_drivers), 
				clinical_analysis, patient_sets_for_drivers, patient_clinical_data, do_plot = TRUE, MIN_MUT_PATIENTS, NULL_COVARIATES, 
				precomp_clin_stats = NULL))
dev.off()
file_open_call2(fname)

driver_clinical_stats$analysis = NULL_ANALYSIS_TAG
save(driver_clinical_stats, file = pff(c("driver_clinical_stats__", NULL_ANALYSIS_TAG, ".rsav")))


# merge different types of alterations per gene, repeat clinical analyses using covariates
load(file = pff("patient_sets_for_drivers.rsav"))
length(patient_sets_for_drivers)

# collapse across alteration types
RGA_stack = data.frame(as.matrix(stack(patient_sets_for_drivers)), stringsAsFactors = FALSE)
colnames(RGA_stack) = c("patient", "RGA")
RGA_stack$RGA = gsub("(.+)__(.+)", "\\1", RGA_stack$RGA)
RGA_list = by(RGA_stack$patient, RGA_stack$RGA, unique)

fname = pff(paste0("figures/KMplots_drivers_COLLAPSED__", MFS_ANALYSIS_TAG, ".pdf"))
pdf(fname, width = 7, height = 7)

cat("todo:", length(RGA_list), "\n")
driver_clinical_stats = do.call(rbind, 
		lapply(1:length(RGA_list), 
				clinical_analysis, RGA_list, patient_clinical_data, do_plot = TRUE, MIN_MUT_PATIENTS, MFS_COVARIATES, 
				precomp_clin_stats = NULL))
dev.off()
file_open_call2(fname)

driver_clinical_stats$analysis = MFS_ANALYSIS_TAG
save(driver_clinical_stats, file = pff(c("driver_clinical_stats_collapsed__", MFS_ANALYSIS_TAG, ".rsav")))


# merge different types of alterations per gene, repeat clinical analyses using covariates
# compute clinical stats with age stage grade as covariates, ERG-neg tumors only

load(file = pff("patient_sets_for_drivers.rsav"))
load(file = pff("patient_clinical_data.rsav"))

ETS_pos_tumors = patient_sets_for_drivers$ETS__SV


# collapse across alteration types
RGA_stack = data.frame(as.matrix(stack(patient_sets_for_drivers)), stringsAsFactors = FALSE)
colnames(RGA_stack) = c("patient", "RGA")
RGA_stack$RGA = gsub("(.+)__(.+)", "\\1", RGA_stack$RGA)
RGA_list = by(RGA_stack$patient, RGA_stack$RGA, unique)

RGA_list = lapply(RGA_list, function(x) setdiff(x, ETS_pos_tumors) )
RGA_list = RGA_list[sapply(RGA_list, length) > MIN_MUT_PATIENTS]

patient_clinical_data = patient_clinical_data[!patient_clinical_data$patient %in% ETS_pos_tumors, ]

fname = pff(paste0("figures/KMplots_drivers_ETSneg__", MFS_ANALYSIS_TAG, ".pdf"))
pdf(fname, width = 7, height = 7)

cat("todo:", length(RGA_list), "\n")
driver_clinical_stats_ETSneg = do.call(rbind, 
		lapply(1:length(RGA_list), 
				clinical_analysis, RGA_list, patient_clinical_data, do_plot = TRUE, MIN_MUT_PATIENTS, MFS_COVARIATES, 
				precomp_clin_stats = NULL))
dev.off()
file_open_call2(fname)

driver_clinical_stats_ETSneg$analysis = MFS_ANALYSIS_TAG
save(driver_clinical_stats_ETSneg, file = pff(c("driver_clinical_stats_ETSneg__", MFS_ANALYSIS_TAG, ".rsav")))




# supp figure for all MFS associations of RGAs
library(gtools)
library(patchwork)

load(pff("all_patients.rsav"))
MIN_DRIVER_FREQ = 0.01
MIN_MUT_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)
MIN_MUT_PATIENTS

FDR_CUTOFF = 0.05

load(file = pff(c("driver_clinical_stats__", MFS_ANALYSIS_TAG, ".rsav")))
driver_clinical_stats = driver_clinical_stats[driver_clinical_stats$clin_var %in% c("mets_time"), ]
driver_clinical_stats = driver_clinical_stats[driver_clinical_stats$n_mut > MIN_MUT_PATIENTS,]
driver_clinical_stats$id = gsub("_", " ", driver_clinical_stats$id)
driver_clinical_stats$fdr = p.adjust(driver_clinical_stats$pval, method = "fdr")
driver_clinical_stats$is_fdr_signf = driver_clinical_stats$fdr < FDR_CUTOFF

combined_p = by(driver_clinical_stats$pval, driver_clinical_stats$id, prod, na.rm = T)
driver_clinical_stats$id = factor(driver_clinical_stats$id, levels = names(sort(-log10(combined_p))))
driver_clinical_stats$fdr_score = -log10(driver_clinical_stats$pval) * sign(log2(driver_clinical_stats$effect))

driver_clinical_stats_MFS_selected = driver_clinical_stats
save(driver_clinical_stats_MFS_selected, file = pff(c("driver_clinical_stats_MFS_selected__", MFS_ANALYSIS_TAG, ".rsav")))


plt_title = paste0("MFS & BCR, FDR < ", FDR_CUTOFF, "; \nmin_mut=", MIN_MUT_PATIENTS, "\n", paste(MFS_COVARIATES, collapse= ","))

plt1 = ggplot(driver_clinical_stats, aes(clin_var, id, fill = fdr_score, label = stars.pval(pval), color = is_fdr_signf)) + 
		geom_tile(aes(color = is_fdr_signf)) + 
		geom_text() + 
		scale_fill_gradient2(high = "darkred", mid = "white", low = "darkblue") + 
		scale_color_manual(values = c("TRUE" = "black", "FALSE" = NA)) + 
		plot_theme() + 
		theme(legend.position = "bottom") + 
		ggtitle(NULL, plt_title)

driver_clinical_stats$HR_label = paste0(
		signif(driver_clinical_stats$effect, 2), " [", 
		signif(driver_clinical_stats$effect_lo, 2), "--", 
		signif(driver_clinical_stats$effect_hi, 2), "]") 

plt2 = ggplot(driver_clinical_stats, aes(id, effect, ymin = effect_lo, ymax = effect_hi, label = HR_label, color = clin_var)) + 
		geom_point(position = position_dodge(width = 0.9)) + 
		geom_errorbar(position = position_dodge(width = 0.9)) + 
		geom_text() + 
		coord_flip() +
		geom_hline(yintercept = 1) + 
		plot_theme() + 
		theme(legend.position = "bottom")


plt_combined  = (plt1 | plt2) + plot_layout(widths = c(1, 3))
fname = pff("figures/MFS_RGA_tile_plot.pdf")
ggsave(plt_combined, file = fname, height = 12, width = 12)
file_open_call2(fname)
