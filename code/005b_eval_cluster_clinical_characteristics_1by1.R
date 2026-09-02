source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
source(paste0("../bin/", date_tag, "/004xx_functions_for_clinical_analyses.R"))
library(gplots)
library(survival)
library(survminer)

load(pff("all_patients.rsav"))
MIN_DRIVER_FREQ = 0.01
MIN_MUT_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)
clust_to_try = SELECT_CLUSTER_METHOD
n_clusters = SELECT_N_CLUSTERS

load(file = pff("patient_clinical_data.rsav"))
load(pff("patient_cluster_membership.rsav"))
members_this_clust = split(patient_cluster_membership$patient, patient_cluster_membership$group_id)
combined_cluster_clinical_stats = combined_cluster_clinical_stats_merged = list()

# cluster size to print
cluster_sizes = sapply(members_this_clust, length)
cluster_sizes = data.frame(cluster_id = gsub("K_(.+)n_(.+)", "\\2", names(cluster_sizes)), 
		cluster_id2 = names(cluster_sizes), 
		n_patients = cluster_sizes, 
		stringsAsFactors = FALSE)
cluster_sizes$cluster_id = paste0("cluster_", cluster_sizes$cluster_id)

# plot stats for clusters one by one
fname = pff(paste0("figures/Clin_cluster_eval__TOP_OPTIONS__", clust_to_try, "__K_", n_clusters, ".pdf"))
pdf(fname, width = 7, height = 8)
driver_clinical_stats2 = do.call(rbind, 
		lapply(1:length(members_this_clust), 
				clinical_analysis, members_this_clust, patient_clinical_data, do_plot = TRUE, MIN_MUT_PATIENTS, MFS_COVARIATES, 
				precomp_clin_stats = NULL))
dev.off()
file_open_call2(fname)

n_clusters_tag = paste0("n_", n_clusters)
combined_cluster_clinical_stats[[clust_to_try]][[n_clusters_tag]] = driver_clinical_stats2
save(combined_cluster_clinical_stats, file = pff("combined_cluster_clinical_stats.rsav"))