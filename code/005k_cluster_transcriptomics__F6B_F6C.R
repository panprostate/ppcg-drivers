source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library("ActivePathways")
library("umap")

FDR_CUTOFF = 0.05
clust_to_try = SELECT_CLUSTER_METHOD
n_clusters = SELECT_N_CLUSTERS

# keep only the patients with RNAseq data
load(pff("expr_mat.rsav"))
load(pff("patient_cluster_membership.rsav"))
patient_cluster_map = structure(names = patient_cluster_membership$patient, patient_cluster_membership$group_id)
patient_cluster_map = patient_cluster_map[names(patient_cluster_map) %in% colnames(expr_mat)]
# exclude a subset of samples that have no cluster (ie no RGAs)
expr_mat = expr_mat[,colnames(expr_mat) %in% names(patient_cluster_map)]

# eval if gene is upregulated in given cluster compared to rest
test_difexp = function(gene, expr_mat, samples_here, cluster_id) {
	
	if (runif(1)<0.001) cat("o")
	
	vals = expr_mat[gene,]
	cluster_vals =  vals[samples_here]
	out_vals = vals[setdiff(names(vals), samples_here)]	
	stat_test = wilcox.test(cluster_vals, out_vals, alter = "g")
	pval = stat_test$p.value
	stat = stat_test$statistic
	mean_clust = mean(cluster_vals)
	mean_out = mean(out_vals)
	fc = mean_clust - mean_out
	data.frame(gene, cluster_id, pval, fc, stat, mean_clust, mean_out)
}

# cluster_id = "cluster_2"
difexp_per_cluster = function(cluster_id, patient_cluster_map, expr_mat) {
	cat("\n", cluster_id, "\n")
	samples_here = names(which(patient_cluster_map == cluster_id))
	all_genes = rownames(expr_mat)
	expr_res = do.call(rbind, mclapply(all_genes, test_difexp, expr_mat, samples_here, cluster_id, mc.cores = 8))
	expr_res = expr_res[order(expr_res$pval),]
	expr_res
}

all_clusters = unique(patient_cluster_map)
difexpr_results = do.call(rbind, lapply(all_clusters, difexp_per_cluster, patient_cluster_map, expr_mat))
difexpr_results$fdr = p.adjust(difexpr_results$pval, method = "fdr")
save(difexpr_results, file = pff(c("difexpr_results__", clust_to_try, "__K", n_clusters, ".rsav")))



# run pathway analysis for each cluster separately
# get pvalues for each cluster as a 1d matrix
difexpr_pval_mat = dcast(gene ~ cluster_id, data = difexpr_results, value.var = "pval")
difexpr_pval_mat1 = as.matrix(difexpr_pval_mat[,-1])
rownames(difexpr_pval_mat1) = difexpr_pval_mat[,1]
save(difexpr_pval_mat1, file = pff(c("difexpr_pval_mat1__", clust_to_try, "__K", n_clusters, ".rsav")))

# custom bg from pvalue matrix and gmt
GMT_fname = pff("gobp_reac.gmt")
custom_bg = readLines(GMT_fname)
custom_bg = strsplit(custom_bg, split = "\t")
custom_bg = unique(unlist(sapply(custom_bg, '[', -1:-2)))
custom_bg = intersect(custom_bg, rownames(difexpr_pval_mat1))


# loop APW for each cluster, no multi test correction here; correct for full set of tests below.
cluster_set = colnames(difexpr_pval_mat1)
APW_results_cluster_expr_by_cluster = mclapply(cluster_set, function(this_cluster) {
		cat(this_cluster, ": \n\n")
		res = ActivePathways(difexpr_pval_mat1[, this_cluster, drop = FALSE], GMT_fname, 
				geneset_filter = c(25, 250),
				background = custom_bg,
				correction_method = "none",
				significant = 1)
		res$cluster = this_cluster
		res

}, mc.cores = 8)

APW_results_cluster_expr_by_cluster = do.call(rbind, APW_results_cluster_expr_by_cluster)
APW_results_cluster_expr_by_cluster$fwer = p.adjust(APW_results_cluster_expr_by_cluster$adjusted_p_val, method = "holm")
APW_res_per_cl = as.data.frame(APW_results_cluster_expr_by_cluster[APW_results_cluster_expr_by_cluster$fwer < FDR_CUTOFF,])
save(APW_res_per_cl, file = pff("APW_res_per_cl.rsav"))


# replicate the file structure for activePathways enrichment map
# pathway list with term_id	term_name	adjusted_p_val
pathways_table_for_enrichment_map = do.call(rbind, by(APW_res_per_cl, APW_res_per_cl$term_id, 
		function(x) x[which.min(x$adjusted_p_val), c("term_id", "term_name", "adjusted_p_val")]))

# GMT file with relevant term_names 
GMT_fname = pff("gobp_reac.gmt")
GMT_content = readLines(GMT_fname)
GMT_content = strsplit(GMT_content, split = "\t")

select_selected = which(sapply(GMT_content, '[[', 1) %in% pathways_table_for_enrichment_map$term_id)
GMT_content = GMT_content[select_selected]
GMT_content = sapply(GMT_content, paste, collapse = "\t")


# subgroups file with pathway/cluster mappings and instruct
cluster_onehot = do.call(rbind, by(APW_res_per_cl$cluster, APW_res_per_cl$term_id, function(x) 0 + (cluster_set %in% x) ))
colnames(cluster_onehot) = cluster_set
cluster_onehot = data.frame(term_id = rownames(cluster_onehot), cluster_onehot, stringsAsFactors = FALSE)
cluster_onehot$instruct = paste0(
		'piechart: attributelist="', 
		 paste(cluster_set, collapse = ","),
		'" colorlist="',
		paste( CLUSTER_COLORS[1:length(cluster_set)], collapse = ","),
		'" showlabels=FALSE'
)

# color legend pdf
dfr = data.frame(cl = cluster_set, fill = as.numeric(gsub("cluster_", "", cluster_set)), stringsAsFactors = FALSE)
cluster_colors = CLUSTER_COLORS[as.numeric(gsub("cluster_", "", cluster_set))]
names(cluster_colors) = cluster_set

plt = ggplot(dfr, aes(cl, factor(1), color = cl)) + 
		scale_color_manual (values = cluster_colors) +
		geom_point() + 
		plot_theme()



new_CS_file_tag = pff(c("figures/ActivePathways__by_cluster_expr__", clust_to_try, "__K", n_clusters, "__"))

fname = paste0(new_CS_file_tag, "pathways.txt")
write.table(pathways_table_for_enrichment_map, file = fname, quote = F, sep = "\t", row.names = FALSE)
file_open_call2(fname)

fname = paste0(new_CS_file_tag, "GMT.gmt")
writeLines(GMT_content, fname)
file_open_call2(fname)

fname = paste0(new_CS_file_tag, "subgroups.txt")
write.table(cluster_onehot, file = fname, quote = F, sep = "\t", row.names = FALSE)
file_open_call2(fname)

fname = paste0(new_CS_file_tag, "colors.pdf")
ggsave(plt, file = fname)
file_open_call2(fname)
