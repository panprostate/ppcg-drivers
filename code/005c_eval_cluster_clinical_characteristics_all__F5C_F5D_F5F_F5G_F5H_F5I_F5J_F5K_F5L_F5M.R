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
input_data_list = list()
cluster_method = SELECT_CLUSTER_METHOD
this_k = SELECT_N_CLUSTERS

# add cluster identity to clinical dataset
load(pff("patient_cluster_membership.rsav"))
clin1 = merge(patient_clinical_data, patient_cluster_membership, by = "patient")

# add driver counts into the combination plot
load(pff("patient_sets_for_drivers.rsav"))
drivers_stack = data.frame(stack(patient_sets_for_drivers), stringsAsFactors = FALSE)
drivers_stack = data.frame(drivers_stack, stringsAsFactors = FALSE)
colnames(drivers_stack) = c("patient", "rga")
drivers_stack = stack(by(drivers_stack$rga, drivers_stack$patient, function(x) length(unique(x))))
colnames(drivers_stack) = c("n_rgas", "patient")
drivers_stack$patient = as.character(drivers_stack$patient)
clin1_1 = merge(clin1, drivers_stack, by = "patient")

plot_tag = paste0(cluster_method, "___", this_k)
group_colors = CLUSTER_COLORS[1:this_k]
names(group_colors) = paste0("cluster_", 1:this_k)

input_data_list[[plot_tag]][["group_colors"]] = group_colors
input_data_list[[plot_tag]][["clin1"]] = clin1_1

combo_plot_and_stats = get_multi_group_report(plot_tag, input_data_list, MFS_COVARIATES, MIN_MUT_PATIENTS)
combo_plot = combo_plot_and_stats[[1]]
cluster_KMplots_survival_stats = combo_plot_and_stats[[2]]
cluster_KMplots_additional_stats = combo_plot_and_stats[[3]]

fname = pff(c("figures/cluster_KMplots_TOP_OPTIONS_NEW__", plot_tag, "__", MFS_ANALYSIS_TAG, ".pdf"))
pdf(fname, width = 12, height = 12)
print(combo_plot)
dev.off()
file_open_call2(fname)

save(cluster_KMplots_survival_stats, 
		file = pff(c("cluster_KMplots_TOP_OPTIONS_NEW_survival_stats_", plot_tag, "__", MFS_ANALYSIS_TAG, ".rsav")))
save(cluster_KMplots_additional_stats, 
		file = pff(c("cluster_KMplots_TOP_OPTIONS_NEW_additional_stats_", plot_tag, "__", MFS_ANALYSIS_TAG, ".rsav")))