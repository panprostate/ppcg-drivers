source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
load(pff("all_patients.rsav"))

FDR_CUTOFF = 0.05
FDR_CAP = 1e-16
LOG_OR_CAP = 4
clust_to_try = SELECT_CLUSTER_METHOD
n_clusters = SELECT_N_CLUSTERS

load(pff("patient_cluster_membership.rsav"))
RGA_clusters_list = split(patient_cluster_membership[,1], patient_cluster_membership[,2])
names(RGA_clusters_list) = gsub("^cluster", "RGA", names(RGA_clusters_list))


mutproc_clusters = read.delim(
		paste0("DATA_USED__", this_timestamp, "/JW_mutProcassignment_2025-11-05/PPCG_mutsig_IMF.tsv"), 
		stringsAsFactors = FALSE)
mutproc_clusters$patient = mutproc_clusters$PPCG_Donor_ID
mutproc_clusters$Cluster_ID = mutproc_clusters$IMF

mutproc_clusters_list = split(mutproc_clusters$patient, mutproc_clusters$Cluster_ID)
names(mutproc_clusters_list) = gsub("Cluster-", "MP_", names(mutproc_clusters_list))
mutproc_clusters_list = lapply(mutproc_clusters_list, unique)
mutproc_clusters_list$MP_noCIN = setdiff(all_patients, unique(unlist(mutproc_clusters_list)))

todo_list = as.matrix(expand.grid(names(mutproc_clusters_list), names(RGA_clusters_list)))


test_cluster_intersection = function(i, todo_list, mutproc_clusters_list, RGA_clusters_list, all_patients) { 
	MP_cluster = todo_list[i, 1]
	RGA_cluster = todo_list[i, 2]
	
	MP_patients = mutproc_clusters_list[[MP_cluster]]
	RGA_patients = RGA_clusters_list[[RGA_cluster]]
	
	common_patients = intersect(MP_patients, RGA_patients)
	n_common = length(common_patients)
	n_MP = length(MP_patients)
	n_RGA = length(RGA_patients)
	
	
	ftest = fisher.test(all_patients %in% MP_patients, all_patients %in% RGA_patients)
	pval = ftest$p.value
	or = ftest$estimate
	
	percent_MP = length(common_patients) / n_MP
	percent_RGA = length(common_patients) / n_RGA
	
	data.frame(i, RGA_cluster, MP_cluster, pval, or, percent_RGA, percent_MP, stringsAsFactors = FALSE)
}


cluster_ixn_stats = do.call(rbind, lapply(1:nrow(todo_list), 
		test_cluster_intersection, todo_list, mutproc_clusters_list, RGA_clusters_list, all_patients))
cluster_ixn_stats$fdr = p.adjust(cluster_ixn_stats$pval, method = "fdr")
cluster_ixn_stats_signf = cluster_ixn_stats[cluster_ixn_stats$fdr < FDR_CUTOFF,]

# cap values
cluster_ixn_stats_signf$fdr_cap = pmax(cluster_ixn_stats_signf$fdr, FDR_CAP)
cluster_ixn_stats_signf$log2_or_cap = log2(cluster_ixn_stats_signf$or)
cluster_ixn_stats_signf$log2_or_cap = pmax(cluster_ixn_stats_signf$log2_or_cap, -LOG_OR_CAP)
cluster_ixn_stats_signf$log2_or_cap = pmin(cluster_ixn_stats_signf$log2_or_cap, LOG_OR_CAP)
cluster_ixn_stats_signf$label_RGA_pct = paste0(100 * signif(cluster_ixn_stats_signf$percent_RGA, 2), "%")

plot_title = paste0("clusters: drivers vs mutProcesses;\nFDR<", FDR_CUTOFF, "; FDR_cap=", FDR_CAP, "; log2OR_cap=", LOG_OR_CAP)
plt = ggplot(cluster_ixn_stats_signf, aes(MP_cluster, RGA_cluster, size = -log10(fdr_cap), color = log2_or_cap, label = label_RGA_pct)) + 
		geom_point() + 
		geom_text(color = "darkgrey", size = 4, vjust = 2) + 
		geom_point(color = "black", shape = 21) + 
		scale_color_distiller(palette = "RdBu") + 
		plot_theme() + 
		ggtitle(NULL, plot_title)

fname = pff(c("figures/MutProc_clusters_comparison", SELECT_CLUSTER_METHOD, "_k", SELECT_N_CLUSTERS, ".pdf"))
ggsave(plt, file = fname)
file_open_call2(fname)
