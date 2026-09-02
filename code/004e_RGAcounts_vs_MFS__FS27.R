source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
source(paste0("../bin/", date_tag, "/004xx_functions_for_clinical_analyses.R"))
library(gplots)
library(survival)
library(survminer)
library(patchwork)
library(ggsignif)
library(gtools)

load(pff("all_patients.rsav"))
MIN_DRIVER_FREQ = 0.01
MIN_MUT_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)

load(file = pff("patient_clinical_data.rsav"))
load(file = pff("patient_sets_for_drivers.rsav"))

patient_driver_stack = stack(patient_sets_for_drivers)
patient_driver_stack[, 2] = as.character(patient_driver_stack[, 2])
patient_driver_stack$mut_type = sapply(strsplit(patient_driver_stack[,2], split = "__"), '[[', 2)
patient_driver_stack$mut_type[patient_driver_stack$mut_type %in% c("gain", "loss", "hAMP")] = "CNA"
patient_driver_stack$mut_type[patient_driver_stack$mut_type %in% c("SNV_CDS", "SNV_NC")] = "SNV_indel"
colnames(patient_driver_stack) = c("patient", "rga", "mut_type")

# separate subcategory for all RGAs
patient_driver_stack_for_all = patient_driver_stack
patient_driver_stack_for_all$mut_type = "ALL_RGA"
patient_driver_stack = rbind(patient_driver_stack, patient_driver_stack_for_all)

# group patients by RGA counts, median-based
get_patient_sets_above_RGA_median = function(mut_type, patient_driver_stack, all_patients) {
	
	stack_here = patient_driver_stack[patient_driver_stack$mut_type == mut_type,]
	dfr = c(unlist(table(stack_here$patient)))[all_patients]
	names(dfr) = all_patients
	dfr[is.na(dfr)] = 0
	dfr = data.frame(n_rgas = dfr, patient = names(dfr), stringsAsFactors = FALSE)
	
	median_rga_count = median(dfr$n_rgas)
	dfr$above_median = dfr$n_rgas > median_rga_count
	dfr$median_rga_count = median_rga_count
	
	patients_above_median = dfr[dfr$above_median, "patient"]	
	res = list(patients_above_median)
	names(res) = paste0(mut_type, "_above_median_", median_rga_count)
	res
}

patient_sets_by_median_RGA_counts = lapply(unique(patient_driver_stack$mut_type), get_patient_sets_above_RGA_median, patient_driver_stack, all_patients)
patient_set_names = sapply(patient_sets_by_median_RGA_counts, function(x) names(x[1]))
patient_sets_by_median_RGA_counts = lapply(patient_sets_by_median_RGA_counts, '[[', 1)
names(patient_sets_by_median_RGA_counts) = patient_set_names


# compute and visualise median-RGA stats with age stage grade as covariates
fname = pff(paste0("figures/clinical_by_RGAcounts_binary_comparision__", MFS_ANALYSIS_TAG, ".pdf"))
pdf(fname, width = 7, height = 7)

cat("todo:", length(patient_sets_by_median_RGA_counts), "\n")
clinical_by_RGAcounts_stats = do.call(rbind, 
		lapply(1:length(patient_sets_by_median_RGA_counts), 
				clinical_analysis, patient_sets_by_median_RGA_counts, patient_clinical_data, do_plot = TRUE, MIN_MUT_PATIENTS, MFS_COVARIATES, 
				precomp_clin_stats = NULL))
dev.off()
file_open_call2(fname)

clinical_by_RGAcounts_stats$analysis = MFS_ANALYSIS_TAG
save(clinical_by_RGAcounts_stats, file = pff(c("clinical_by_RGAcounts_stats__", MFS_ANALYSIS_TAG, ".rsav")))
