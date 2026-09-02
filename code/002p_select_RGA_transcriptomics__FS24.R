source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library("gtools")
library("ActivePathways")
library("ggrepel")


FDR_CUTOFF = 0.05


load(pff("cgc2024.rsav"))
load(pff("prostate_cancer_genes.rsav"))
cancer_genes = c(cgc2024, prostate_cancer_genes)



test_difexp_per_gene = function(gene, expr_mat, tumors_case, tumors_ctrl, RGA_to_test) {
	
	if (runif(1)<0.001) cat("o")
	
	vals_case = unlist(expr_mat[gene, tumors_case, drop = TRUE])
	vals_clrl = unlist(expr_mat[gene, tumors_ctrl, drop = TRUE])
	
	stat_test = wilcox.test(vals_case, vals_clrl)
	pval = stat_test$p.value
	stat = stat_test$statistic

	mean_case = mean(vals_case)
	mean_ctrl = mean(vals_clrl)
	fc = mean_case - mean_ctrl

	data.frame(gene, pval, fc, stat, mean_case, mean_ctrl, RGA_to_test, stringsAsFactors = FALSE)
}


test_difexp = function(RGA_to_test, patient_sets_for_drivers, expr_mat) {
	
	tumors_case = intersect(patient_sets_for_drivers[[RGA_to_test]], colnames(expr_mat))
	tumors_ctrl = setdiff(colnames(expr_mat), patient_sets_for_drivers[[RGA_to_test]])
	
	difexp_values = do.call(rbind, mclapply(rownames(expr_mat), test_difexp_per_gene, expr_mat, tumors_case, tumors_ctrl, RGA_to_test, mc.cores = 8))
	difexp_values$fdr = p.adjust(difexp_values$pval, method = "fdr")
	difexp_values
}


test_pathways = function(select_RGA_difexp_values, GMT_fname, APW_filename_tag) {
	
	difexp_matrix = as.matrix(select_RGA_difexp_values[, "pval", drop = FALSE])
	rownames(difexp_matrix) = select_RGA_difexp_values$gene
	
	# custom bg from pvalue matrix and gmt
	custom_bg = readLines(GMT_fname)
	custom_bg = strsplit(custom_bg, split = "\t")
	custom_bg = unique(unlist(sapply(custom_bg, '[', -1:-2)))
	custom_bg = intersect(custom_bg, rownames(difexp_matrix))

	APW_res = ActivePathways(
			difexp_matrix, 
			GMT_fname, 
			significant = FDR_CUTOFF, 
			cytoscape_file_tag = APW_filename_tag,
			geneset_filter = c(25, 250),
			background = custom_bg)
	APW_res
}


plot_volcano = function(select_RGA_difexp_values, RGA_to_test) {
	
	select_RGA_difexp_values$is_significant = select_RGA_difexp_values$fdr < FDR_CUTOFF
	select_RGA_difexp_values$is_signif_cancer_gene = select_RGA_difexp_values$is_significant & select_RGA_difexp_values$gene %in% cancer_genes
	
	select_RGA_difexp_values$color = "not significant"
	select_RGA_difexp_values$color[select_RGA_difexp_values$is_significant] = "significant"
	select_RGA_difexp_values$color[select_RGA_difexp_values$is_signif_cancer_gene] = "significant, cancer gene"
	point_colors = c("not significant" = "lightgrey", "significant" = "azure4", "significant, cancer gene" = "darkred")
	
	select_RGA_difexp_values$label = select_RGA_difexp_values$gene
	select_RGA_difexp_values[!select_RGA_difexp_values$gene %in% cancer_genes, "label"] = NA
	select_RGA_difexp_values[!select_RGA_difexp_values$is_significant, "label"] = NA
	
	plt_title = paste0(RGA_to_test, "; FDR < ", FDR_CUTOFF, 
			"; n = ", sum(select_RGA_difexp_values$fdr < 0.05),
			"; n_c = ", sum(!is.na(select_RGA_difexp_values$label)))

	plt = ggplot(select_RGA_difexp_values, aes(fc, -log10(fdr), color = color, label = label)) + 
			geom_point() + 
			scale_color_manual(values = point_colors) + 
			plot_theme() + 
			scale_x_continuous("fold-change, log2") +
			scale_y_continuous("FDR, log10") +
			geom_text_repel(size = 3, color = "purple", max.overlaps = Inf) + 
			ggtitle(NULL, plt_title) + 
			theme(legend.position = "bottom")
	plt
	
}




RGA_to_test = "COL5A1__gain"

load(pff("expr_mat.rsav"))
load(pff("patient_sets_for_drivers.rsav"))

# remove from expression matrix the genes located in the locus of the CNA
load(pff("CNA_modules.rsav"))
locus_target_genes = CNA_modules[CNA_modules$TGT_GENE == "COL5A1", "all_genes"]
locus_target_genes = strsplit(locus_target_genes, split = ",")[[1]]
expr_mat = expr_mat[!rownames(expr_mat) %in% locus_target_genes,]

# difexp genes
select_RGA_difexp_values = test_difexp(RGA_to_test, patient_sets_for_drivers, expr_mat)
save(select_RGA_difexp_values, file = pff(c("select_RGA_difexp_values__", RGA_to_test, ".rsav")))

# pathway enrichment
GMT_fname = pff("gobp_reac.gmt")
APW_filename_tag = pff(c("figures/select_RGA_pathway_analysis__", RGA_to_test, "__"))
select_RGA_enriched_pathways = test_pathways(select_RGA_difexp_values, GMT_fname, APW_filename_tag)
save(select_RGA_enriched_pathways, file = pff(c("select_RGA_enriched_pathways__", RGA_to_test, ".rsav")))

# volcano plot
plt = plot_volcano(select_RGA_difexp_values, RGA_to_test)
fname = pff(c("figures/select_RGA_difexp_volcano__", RGA_to_test, ".pdf"))
ggsave(plt, file = fname)
file_open_call2(fname)



RGA_to_test = "ZNF280D__SNV_NC"

load(pff("expr_mat.rsav"))
load(pff("patient_sets_for_drivers.rsav"))

# difexp genes
select_RGA_difexp_values = test_difexp(RGA_to_test, patient_sets_for_drivers, expr_mat)
save(select_RGA_difexp_values, file = pff(c("select_RGA_difexp_values__", RGA_to_test, ".rsav")))

# pathway enrichment
GMT_fname = pff("gobp_reac.gmt")
APW_filename_tag = pff(c("figures/select_RGA_pathway_analysis__", RGA_to_test, "__"))
select_RGA_enriched_pathways = test_pathways(select_RGA_difexp_values, GMT_fname, APW_filename_tag)
save(select_RGA_enriched_pathways, file = pff(c("select_RGA_enriched_pathways__", RGA_to_test, ".rsav")))

# volcano plot
plt = plot_volcano(select_RGA_difexp_values, RGA_to_test)
fname = pff(c("figures/select_RGA_difexp_volcano__", RGA_to_test, ".pdf"))
ggsave(plt, file = fname)
file_open_call2(fname)

