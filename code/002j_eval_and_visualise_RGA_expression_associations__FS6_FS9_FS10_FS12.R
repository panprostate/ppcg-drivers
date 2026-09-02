source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(ggpubr)
library(ggrepel)
library(gtools)
library(patchwork)

load(pff("expr_mat.rsav"))
load(pff("cgc2024.rsav"))
load(pff("prostate_cancer_genes.rsav"))
cancer_genes = c(cgc2024, prostate_cancer_genes)

FDR_CUTOFF = 0.05
N_MIN_SPLS = 3


test_gene_expression = function(id, patient_sets, expr_mat, mut_tag, N_MIN_SPLS, cancer_genes, do_plot = FALSE) {
	
	cat(id, " ")
	mut_patients = patient_sets[[id]]
	gene = gsub("(.+)__(.+)__(.+)", "\\1", id)
	mut_type = gsub("(.+)__(.+)__(.+)", "\\2", id)
	origin = gsub("(.+)__(.+)__(.+)", "\\3", id)
		
	# also need to have mrna
	mut_patients = intersect(colnames(expr_mat), mut_patients)
	notmut_patients = setdiff(colnames(expr_mat), mut_patients)

	if (length(mut_patients) < N_MIN_SPLS) {
		cat(gene, "has only", length(mut_patients), "samples; skipping\n")
		return(list())
	}
	if (!gene %in% rownames(expr_mat)) {
		cat(gene, "not shown in exprMat; skipping\n")
		return(list())		
	}

	dfr = rbind(
			data.frame(patient = mut_patients, expr = expr_mat[gene,  mut_patients], type = "mut", 
					stringsAsFactors = FALSE),
			data.frame(patient = notmut_patients, expr = expr_mat[gene,  notmut_patients], type = "wt", 
					stringsAsFactors = FALSE))
	plt = NULL
	if (do_plot) {
		ggtitl = paste0(id, "\nn=", length(mut_patients), "\n")
		if (gene %in% cancer_genes)
			ggtitl = paste0(ggtitl, "###")
		plt = ggplot(dfr, aes(type, expr, fill = type)) + 
				scale_y_continuous("expression, log2") +
				scale_fill_manual(values = c("mut" = "darkslategrey", "wt" = "lightgrey")) +
				geom_boxplot(width = 0.75, outlier.shape = NA) + 
				geom_jitter(width = 0.1, alpha = 0.5, color = "darkgrey") + 
				stat_compare_means(method = "wilcox.test", comparisons = list(c("mut", "wt"))) + 
				plot_theme() + 
				ggtitle(ggtitl) + 
				theme(plot.title = element_text(size = 8), plot.subtitle = element_text(size = 8), legend.position = "none")
	}
	
	mut_vals =  dfr[dfr$type == "mut", "expr"]
	notmut_vals =  dfr[dfr$type == "wt", "expr"]
	mut_mean =  mean(mut_vals)
	notmut_mean = mean(notmut_vals)
	log2FC = mut_mean - notmut_mean
	pval = wilcox.test(mut_vals, notmut_vals)$p.value
	patients_concat = paste(mut_patients, collapse = ",")
	
	dfr = data.frame(gene, pval, log2FC, mut_mean, notmut_mean, n = length(mut_vals), mut_type, mut_tag, id, origin, 
			donors = patients_concat, 
			stringsAsFactors = FALSE)
	
	return(list(dfr, plt))
}	


clean_res_list = function(res_list) {
	res_list = res_list[sapply(res_list, length) == 2]
	res_dfr = do.call(rbind, lapply(res_list, '[[', 1))
	res_dfr = res_dfr[order(res_dfr$pval),]
	res_dfr$fdr = p.adjust(res_dfr$pval, method = "fdr")
	res_dfr
}



#
# test RGA / expression associations
#
load(pff("patient_sets_for_drivers.rsav"))

# exclude protein-coding mutations
patient_sets_for_drivers = patient_sets_for_drivers[grep("__SNV_CDS", names(patient_sets_for_drivers), invert = TRUE)]

# add ETS factors from the consensus dataset
load(file = pff("dat_ETS.rsav"))
ETS_genes = c("ERG", "ETV1", "ETV4", "ETV5", "FLI1")
ETS_patients = sapply(ETS_genes, function(g)  dat_ETS[dat_ETS[, g] == TRUE, "PPCG_Donor_ID"])
names(ETS_patients) = paste0(names(ETS_patients), sep = "__SV")

patient_sets_for_drivers = c(patient_sets_for_drivers, ETS_patients)

# add 3rd token to gene set name: same gene as above, to label direct targets and distal targets
names(patient_sets_for_drivers) = 
		paste0(names(patient_sets_for_drivers), "__", gsub("(.+)__(.+)", "\\1", names(patient_sets_for_drivers)))


res_list = lapply(names(patient_sets_for_drivers), test_gene_expression, patient_sets_for_drivers, 
		expr_mat, "all_muts", N_MIN_SPLS, cancer_genes)
res_dfr = clean_res_list(res_list)

res_dfr$gene_type = "other"
res_dfr[res_dfr$gene %in% cgc2024, "gene_type"] = "CGC"
res_dfr[res_dfr$gene %in% prostate_cancer_genes, "gene_type"] = "PrCa"

drivers_with_expression_changes = res_dfr
save(drivers_with_expression_changes, file = pff("drivers_with_expression_changes.rsav"))




## barplot
res_dfr_here = res_dfr[res_dfr$fdr < FDR_CUTOFF, ]
res_dfr_here$gene = factor(res_dfr_here$gene, levels = unique(res_dfr_here[order(res_dfr_here$log2FC),"gene"]))
res_dfr_here$mut_type = factor(res_dfr_here$mut_type, levels = c("loss", "BI", "gain", "hAMP", "SV", "SNV_NC"))

gene_colors = c("other" = "cornflowerblue", "PrCa" = "darkred", "CGC" = "darkorange")
plot_subtitle = paste0("n=", nrow(res_dfr_here), "; FDR < ", FDR_CUTOFF, "; min",  N_MIN_SPLS, " spl")

ggplt = ggplot(res_dfr_here, aes(gene, log2FC, label = stars.pval(fdr), fill = gene_type)) + 
		geom_bar(stat = "identity") + 
		geom_text(size = 3) + 
		geom_text(aes(label = n), y = 0, size = 3) + 
		facet_grid(mut_type ~ 1, space = "free", scales = "free_y") + 
		scale_fill_manual(values = gene_colors) +
		coord_flip() + 
		plot_theme() +
		ggtitle("mRNA associations of RGAs", plot_subtitle) +
		theme(axis.text.x = element_text(size = 6), legend.position = "right")


fname = pff("figures/mRNA_vs_muts__LOG2FC_barplot.pdf")
ggsave(ggplt, file = fname, height = 8)
file_open_call2(fname)


# boxplots for top hits
genes_to_plot = res_dfr_here$id
res_list2 = lapply(genes_to_plot, test_gene_expression, patient_sets_for_drivers, expr_mat, "all_muts", 
		N_MIN_SPLS, cancer_genes, do_plot = TRUE)

fname = pff(paste0("figures/mRNA_boxplots_ALL__FDR", FDR_CUTOFF, ".pdf"))
pdf(fname, width = 3, height = 5)

plt_list_smaller = list()

for (i in 1:length(res_list2)) {
	cat(i, " ")
	# take FDR from previous data frame
	this_id = res_list2[[i]][[1]]$id
	fdr = signif(res_dfr[res_dfr$id == this_id, "fdr"], 2)
	n_patients = res_list2[[i]][[1]]$n
	
	plt = res_list2[[i]][[2]] + 
			ggtitle(this_id, paste0("n = ", n_patients, "; fdr = ", fdr )) + 
			theme(legend.position = "none") + 
			scale_x_discrete(NULL) +
			scale_y_continuous(NULL)
	print(plt)
	
	plt_list_smaller[[i]] = plt
}
dev.off()
file_open_call2(fname)



#
# expression of genes connected by chromatin loops
#
load(pff("patient_sets_for_loops.rsav"))
N_MIN_SPLS = 3
FDR_CUTOFF = 0.2

loops_res_list = lapply(names(patient_sets_for_loops), test_gene_expression, patient_sets_for_loops, 
		expr_mat, "all_muts", N_MIN_SPLS, cancer_genes)
loops_res_list = clean_res_list(loops_res_list)

drivers_looped_with_expression_changes = loops_res_list
save(drivers_looped_with_expression_changes, file = pff("drivers_looped_with_expression_changes.rsav"))

# boxplots for top hits
genes_to_plot = loops_res_list[loops_res_list$fdr < FDR_CUTOFF, "id"]
loops_res_list2 = lapply(genes_to_plot, test_gene_expression, patient_sets_for_loops, expr_mat, "loops_muts", 
		N_MIN_SPLS, cancer_genes, do_plot = T)

fname = pff(paste0("figures/mRNA_boxplots_LOOPS__FDR", FDR_CUTOFF, ".pdf"))
pdf(fname, width = 3, height = 5)


for (i in 1:length(loops_res_list2)) {
	cat(i, " ")
	# take FDR from previous data frame
	this_id = loops_res_list2[[i]][[1]]$id
	fdr = signif(loops_res_list[loops_res_list$id == this_id, "fdr"], 2)
	n_patients = loops_res_list2[[i]][[1]]$n
	
	plt = loops_res_list2[[i]][[2]] + 
			ggtitle(this_id, paste0("n = ", n_patients, "; fdr = ", fdr )) + 
			theme(legend.position = "none") + 
			scale_x_discrete(NULL) +
			scale_y_continuous(NULL)
	print(plt)
}
dev.off()
file_open_call2(fname)


# lncRNA expression
load(file = pff("expr_mat_with_LNCs.rsav"))
load(file = pff("prepared_lncRNA.rsav"))
lncRNA_symbols = unique(gsub("(.+)::(.+)::(.+)", "\\2", prepared_lncRNA$id))

FDR_CUTOFF = 1.0
N_MIN_SPLS = 3


load(pff("patient_sets_for_drivers.rsav"))
patient_sets_for_drivers = patient_sets_for_drivers[grep("__hAMP|__BI|__SNV_NC|__SV", names(patient_sets_for_drivers))]
names(patient_sets_for_drivers) = paste0(names(patient_sets_for_drivers), "__", gsub("(.+)__(.+)", "\\1", names(patient_sets_for_drivers)))

# keep only lncRNAs here
patient_sets_for_drivers = patient_sets_for_drivers[gsub("(.+)__(.+)__(.+)", "\\1", names(patient_sets_for_drivers)) %in% lncRNA_symbols]

res_list = lapply(names(patient_sets_for_drivers), test_gene_expression, patient_sets_for_drivers, 
		expr_mat, "all_muts", N_MIN_SPLS, cancer_genes)
res_dfr = clean_res_list(res_list)


drivers_lncrnas_with_expression_changes = res_dfr
save(drivers_lncrnas_with_expression_changes, file = pff("drivers_lncrnas_with_expression_changes.rsav"))


res_dfr$gene_type = "other"
res_dfr[res_dfr$gene %in% cgc2024, "gene_type"] = "CGC"
res_dfr[res_dfr$gene %in% prostate_cancer_genes, "gene_type"] = "PrCa"

# boxplots for top hits
genes_to_plot = res_dfr[res_dfr$fdr < FDR_CUTOFF, "id"]
res_list2 = lapply(genes_to_plot, test_gene_expression, patient_sets_for_drivers, expr_mat, "all_muts", 
		N_MIN_SPLS, cancer_genes, do_plot = T)


fname = pff(paste0("figures/mRNA_boxplots_lncRNA__FDR", FDR_CUTOFF, ".pdf"))
pdf(fname, width = 3, height = 5)

plt_list_smaller = list()

for (i in 1:length(res_list2)) {
	cat(i, " ")
	# take FDR from previous data frame
	this_id = res_list2[[i]][[1]]$id
	fdr = signif(res_dfr[res_dfr$id == this_id, "fdr"], 2)
	n_patients = res_list2[[i]][[1]]$n
	
	plt = res_list2[[i]][[2]] + 
			ggtitle(this_id, paste0("n = ", n_patients, "; fdr = ", fdr )) + 
			theme(legend.position = "none") + 
			scale_x_discrete(NULL) +
			scale_y_continuous(NULL)
	print(plt)
	
	plt_list_smaller[[i]] = plt
}
dev.off()
file_open_call2(fname)

