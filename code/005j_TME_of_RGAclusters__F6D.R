source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(gtools)

load(pff("all_patients.rsav"))
load(pff("patient_cluster_membership.rsav"))
RGA_clusters_list = split(patient_cluster_membership[,1], patient_cluster_membership[,2])
names(RGA_clusters_list) = gsub("^cluster", "RGA", names(RGA_clusters_list))

RGA_clusters_dfr = data.frame(
		patient = patient_cluster_membership$patient, 
		cluster_id = gsub("cluster", "RGA", patient_cluster_membership$group_id), 
		stringsAsFactors = FALSE)

TME_scores = read.csv("DATA_USED__2025-09-24/TME_Lena_2025-11-13/cell_type_deconvolution_V4_2025-11-13.csv", 
		stringsAsFactors = FALSE)
TME_scores$patient = gsub("(PPCG....)(.)_(.+)", "\\1", TME_scores$PPCG_RNA_Assay_ID)

# use the same samples as in RNA-seq analyses
load(file = pff("expr_sample_meta.rsav"))
TME_scores = TME_scores[TME_scores$PPCG_RNA_Assay_ID %in% expr_sample_meta$PPCG_RNA_Assay_ID,]
save(TME_scores, file = pff("TME_scores.rsav"))

TME_scores_with_RGAcluster = merge(RGA_clusters_dfr, TME_scores, by = "patient")

cell_types = colnames(TME_scores_with_RGAcluster)[-1:-3]
cluster_ids = unique(TME_scores_with_RGAcluster[,2])
todo_list = as.matrix(expand.grid(cell_types, cluster_ids))

# i = 1; 	# cell_type = "t_cd8"; cluster_id = "RGA_5"
test_diff = function(i, todo_list, TME_scores_with_RGAcluster) {

	cell_type = todo_list[i, 1]
	cluster_id = todo_list[i, 2]
	
	this_TME_scores = data.frame (cluster_id = TME_scores_with_RGAcluster$cluster_id, cell_score = TME_scores_with_RGAcluster[, cell_type], 
			stringsAsFactors = FALSE)
	this_TME_scores$cell_score_ranknorm = rank(this_TME_scores$cell_score) / nrow(this_TME_scores)
	this_TME_scores$cell_score_ranknorm_z = (this_TME_scores$cell_score_ranknorm - 0.5) / sd(this_TME_scores$cell_score_ranknorm)
	
	this_vals = this_TME_scores[this_TME_scores$cluster_id == cluster_id, "cell_score"]
	other_vals = this_TME_scores[this_TME_scores$cluster_id != cluster_id, "cell_score"]
	mean_rank_z_score = mean(this_TME_scores[TME_scores_with_RGAcluster$cluster_id == cluster_id, "cell_score_ranknorm_z"])
	
	tested = wilcox.test(this_vals, other_vals)
	pval = tested$p.value
	
	data.frame(cell_type, cluster_id, pval, mean_rank_z_score, stringsAsFactors = FALSE)	
}


TME_results = do.call(rbind, lapply(1:nrow(todo_list), test_diff, todo_list, TME_scores_with_RGAcluster))
TME_results$fdr = p.adjust(TME_results$pval, method = "fdr")
TME_results$fdr_stars = stars.pval(TME_results$fdr)

top_cell_types = names(sort(by(-log10(TME_results$pval), TME_results$cell_type, sum)))
TME_results$cell_type = factor(TME_results$cell_type, levels = top_cell_types)

plt = ggplot(TME_results, aes(cluster_id, cell_type, fill = mean_rank_z_score, label = fdr_stars)) + 
		scale_fill_gradient2(high = "red", low = "blue", mid = "white") + 
		geom_tile() + 
		geom_text() + 
		ggtitle("RGA clusters vs TME features", "WX-test on raw, FDR") +
		plot_theme()

fname = pff(c("figures/TME_RGAcluster_heatmap_rankscore_Wx_on_raw_vals_", SELECT_CLUSTER_METHOD, "_k", SELECT_N_CLUSTERS, ".pdf"))
ggsave(plt, file = fname)
file_open_call2(fname)
