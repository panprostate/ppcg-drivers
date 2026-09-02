source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(patchwork)
library(gtools)
library(eulerr)


FDR_CUTOFF = 0.05
MIN_N_PATIENTS = 10

load(pff("cgc2024.rsav"))
load(pff("prostate_cancer_genes.rsav"))
cancer_genes = unique(c(cgc2024, prostate_cancer_genes))

# stat test: do RGAs cooccur or show mutual exclusivity
test_interaction = function(i, todo_list, combined_sample_lists, all_patients) {
	if (i %% 100 == 0) cat(i, " ")
	
	driver1 = todo_list[i, 2]
	driver2 = todo_list[i, 1]
	ixn_tag = paste(sort(c(driver1, driver2)), collapse = ",")
	tumors1 = unique(combined_sample_lists[[driver1]])
	tumors2 = unique(combined_sample_lists[[driver2]])
	common_tumors = intersect(tumors1, tumors2)
	common_tumors_string = paste(common_tumors, collapse = ";")
	
	n_common_tumors = length(common_tumors)
	n1 = length(tumors1)
	n2 = length(tumors2)
	prob1 = length(tumors1) / length(all_patients)
	prob2 = length(tumors2) / length(all_patients)
	exp_n_common = prob1 * prob2 * length(all_patients)
	
	ft = fisher.test (
			factor(all_patients %in% tumors1, levels = c("TRUE","FALSE")), 
			factor(all_patients %in% tumors2, levels = c("TRUE","FALSE")))
	pp = ft$p.value[[1]]
	or = ft$estimate[[1]]
	
	data.frame(driver1, driver2, pp, or, n_common_tumors, n1, n2, common_tumors_string, ixn_tag, exp_n_common,
			stringsAsFactors = FALSE)
}



load(pff("all_patients.rsav"))
load(pff("patient_sets_for_drivers.rsav"))
# exclude GISTIC CNAs
patient_sets_for_drivers = patient_sets_for_drivers[grep("__loss$$", names(patient_sets_for_drivers), invert = TRUE)]

# collapse patient sets for different mutation types
stack_driver_sets = as.matrix(stack(patient_sets_for_drivers))
stack_driver_sets = cbind(stack_driver_sets, gene = gsub("(.+)__(.+)", "\\1", stack_driver_sets[, "ind"]))

patient_sets_for_drivers_collapsed = split(stack_driver_sets[, "values"], stack_driver_sets[, "gene"])
patient_sets_for_drivers_collapsed = lapply(patient_sets_for_drivers_collapsed, unique)
patient_sets_for_drivers_collapsed = patient_sets_for_drivers_collapsed[sapply(patient_sets_for_drivers_collapsed, length) >= MIN_N_PATIENTS]
todo_list = t(combn(names(patient_sets_for_drivers_collapsed), 2))

CNA_SNV_ixns = do.call(rbind, mclapply(1:nrow(todo_list), test_interaction, 
		todo_list, patient_sets_for_drivers_collapsed, all_patients, mc.cores = 8))

# adjust for multiple testing
CNA_SNV_ixns$fdr = p.adjust(CNA_SNV_ixns$pp, method = "fdr")
CNA_SNV_ixns$edge_dir = c("-1"="low", "1"="high")[as.character(sign(log(CNA_SNV_ixns$or)))]
CNA_SNV_ixns = CNA_SNV_ixns[order(CNA_SNV_ixns$pp),]
CNA_SNV_ixns = CNA_SNV_ixns[CNA_SNV_ixns$fdr < FDR_CUTOFF,]

save(CNA_SNV_ixns, file = pff("CNA_SNV_ixns.rsav"))
save(patient_sets_for_drivers_collapsed, file = pff("patient_sets_for_drivers_collapsed.rsav"))

system(paste("mkdir", pff("Zoe_2025-01-13/")))
save(CNA_SNV_ixns, file = paste0(pff("Zoe_2025-01-13/"), "CNA_SNV_ixns.rsav"))
save(patient_sets_for_drivers_collapsed, file = paste0(pff("Zoe_2025-01-13/"), file = "patient_sets_for_drivers_collapsed.rsav"))


# cap FDR and OR values
CNA_SNV_ixns$or_capped = CNA_SNV_ixns$or
CNA_SNV_ixns$fdr_capped = CNA_SNV_ixns$fdr
CNA_SNV_ixns [CNA_SNV_ixns$or > 16, "or_capped"] = 16
CNA_SNV_ixns [CNA_SNV_ixns$or < 1/16, "or_capped"] = 1/16
CNA_SNV_ixns [CNA_SNV_ixns$fdr < 1e-16, "fdr_capped"] = 1e-16
CNA_SNV_ixns$score = -log10(CNA_SNV_ixns$fdr_capped) * sign(log2(CNA_SNV_ixns$or_capped))

# make the resulting matrix symmetrical for visualisation purposes
CNA_SNV_ixns_mirror = CNA_SNV_ixns[, c(2,1, 3:ncol(CNA_SNV_ixns))]
colnames(CNA_SNV_ixns_mirror) = colnames(CNA_SNV_ixns)
CNA_SNV_ixns_symm = rbind(CNA_SNV_ixns, CNA_SNV_ixns_mirror)

# cluster genes by that matrix and order genes for heatmap
this_dist = function(x) proxy::dist(x, by_rows = TRUE, method = "Jaccard")
matrix_for_hclust = dcast(CNA_SNV_ixns_symm, driver1~driver2, value.var = "score", fun.aggregate = sum)
rownames(matrix_for_hclust) = matrix_for_hclust[,1]
matrix_for_hclust = matrix_for_hclust[,-1]
hclustered = hclust(this_dist(as.matrix(matrix_for_hclust)))

# assign gene ordering
gene_ordering_in_heatmap = rownames(matrix_for_hclust)[hclustered$order]
CNA_SNV_ixns_symm$driver1 = factor(CNA_SNV_ixns_symm$driver1, levels = gene_ordering_in_heatmap)
CNA_SNV_ixns_symm$driver2 = factor(CNA_SNV_ixns_symm$driver2, levels = gene_ordering_in_heatmap)

CNA_SNV_ixns_symm$n_spl_ixnd = cut(CNA_SNV_ixns_symm$n_common_tumors, breaks = c(-1, 10, 50, 100, 200, 10000))

ixn_plt = ggplot(CNA_SNV_ixns_symm, aes(driver1, driver2, fill = score)) + 
		geom_tile(color = "black") + 
		geom_point(aes(shape = n_spl_ixnd), size = 1) +
		scale_fill_gradient2( mid = "white", high = "darkorange", low = "darkblue", limits = c(-16,16), breaks = c(-16,-8, 0, 8, 16)) + 
		scale_x_discrete(NULL) +
		scale_y_discrete(NULL) +
		plot_theme() + 
		geom_abline() +
		theme(axis.text.x = element_text(size = 4), axis.text.y = element_text(size = 3)) + 		
		ggtitle (paste0("Driver stat interactions, n = ", length(unique(CNA_SNV_ixns_symm$ixn_tag)), "; FDR < ", FDR_CUTOFF) )
		
dfr_label_cancer_genes = data.frame(driver1_s = gene_ordering_in_heatmap, symbol = gene_ordering_in_heatmap, cancer_genes = "none", 
		stringsAsFactors = FALSE)
dfr_label_cancer_genes$driver1_s = factor(dfr_label_cancer_genes$driver1_s, levels = gene_ordering_in_heatmap)

# label cancer genes here - if they include cgc or prca genes
dfr_label_cancer_genes$symbol = gsub("(.+)__(.+)", "\\1", dfr_label_cancer_genes$symbol)
is_this_cgc_gene = sapply(strsplit(dfr_label_cancer_genes$symbol, split = ","), function(x) any(cgc2024 %in% x))
is_this_prca_gene = sapply(strsplit(dfr_label_cancer_genes$symbol, split = ","), function(x) any(prostate_cancer_genes %in% x))
dfr_label_cancer_genes$cancer_genes [is_this_cgc_gene] = "CGC"
dfr_label_cancer_genes$cancer_genes [is_this_prca_gene] = "PrCa"
dfr_label_cancer_genes[dfr_label_cancer_genes$symbol == "ETS", "cancer_genes"] = "PrCa"

cancer_gene_symbol = c("PrCa" = 1, "CGC" = 4)

label_cancer_genes_plt = ggplot(dfr_label_cancer_genes, aes(factor(1), driver1_s, shape = cancer_genes)) + 
		geom_tile(fill = "white") + 
		geom_point() + 
		scale_fill_manual (values = c("CGC" = "grey", "PrCa" = "black", "none" = "whitesmoke"))  +
		scale_shape_manual(values = cancer_gene_symbol) + 
		plot_theme() + 
		theme(axis.text.y = element_text(size = 3))
		
plt_combined = ixn_plt + label_cancer_genes_plt + plot_layout(widths = unit(c(1, 0.1), 'null'))
		

fname = pff("figures/CNA_SNV_ixn_dotplot.pdf")
ggsave(plt_combined, file = fname, width = 8.95, height = 5)
file_open_call2(fname)


# oncoprints for RGA-RGA interactions
get_oncoprint_plot = function(driver1, CNA_SNV_ixns, driver_patient_sets) {
	
	cat(driver1, "\n")
	other_drivers1 = as.character(CNA_SNV_ixns[CNA_SNV_ixns$driver1 == driver1, "driver2"])
	other_drivers2 = as.character(CNA_SNV_ixns[CNA_SNV_ixns$driver2 == driver1, "driver1"])
	other_drivers = unique(c(other_drivers1, other_drivers2))
	
	ixn_direction = unique(CNA_SNV_ixns[CNA_SNV_ixns$driver1 == driver1, c("driver2", "or", "fdr")])
	ixn_direction$or = as.character(sign(log(ixn_direction$or)))
	ixn_direction = rbind(ixn_direction, c(driver1, 0, 1))
	ixn_direction$or = c("-1" = "exclusive", "1"="co-occur", "0" = "driver")[ ixn_direction$or ]
	ixn_direction = ixn_direction[order(ixn_direction$or, ixn_direction$fdr),]
	
	this_driver_stack = data.frame(as.matrix(stack(driver_patient_sets[driver1])), stringsAsFactors = FALSE)
	other_drivers_stack = data.frame(as.matrix(stack(driver_patient_sets[other_drivers])), stringsAsFactors = FALSE)
	this_driver_stack$type = "this"
	other_drivers_stack$type = "other"	
	driver_stack = rbind(this_driver_stack, other_drivers_stack)
	colnames(driver_stack) = c("tumor", "driver", "type")
	
	driver_stack = merge(driver_stack, ixn_direction, by.x = "driver", by.y = "driver2")
	driver_mat = dcast(driver_stack, driver~tumor, value.var = "tumor", fun.aggregate = length)
	rownames(driver_mat) = driver_mat[,1]
	# use the oncoprint ordering of rows and columns
	driver_mat1 = memoSort(driver_mat[,-1])
	
	driver1_tgts = unique(driver_stack[driver_stack$driver == driver1,"tumor"])
	other_tgts = setdiff(unique(driver_stack[driver_stack$driver != driver1,"tumor"]), driver1_tgts)

	patient_order = colnames(driver_mat1)
	driver_stack$tumor = factor(driver_stack$tumor, levels = patient_order)

	driver_order = c(driver1, 
		intersect(rownames(driver_mat1), ixn_direction[ixn_direction$or == "exclusive", "driver2"]),
		intersect(rownames(driver_mat1), ixn_direction[ixn_direction$or == "co-occur", "driver2"]))
	driver_order = rev(driver_order)
	driver_stack$driver = factor(driver_stack$driver, levels = driver_order)
	
	driver_stack$pval_stars = stars.pval(as.numeric(driver_stack$fdr))
	driver_stack$pval_stars[duplicated(driver_stack$driver)] = NA
	
	plt = ggplot(driver_stack, aes(tumor, driver, 1, fill = or, label = pval_stars)) + 
			geom_tile() + 
			scale_x_discrete("tumors", breaks = NULL, labels = NULL) +
			plot_theme() + 
			geom_text(x = 1) + 
			theme(legend.position = "bottom") + 
			scale_fill_manual(NULL, values = c("driver" = "darkgray", "co-occur" = "red3", "exclusive" = "royalblue")) + 
			ggtitle(NULL, paste0(driver1, "\n n_ixn=", length(driver_order)-1, "; n_spl=", length(patient_order)))
}	


all_drivers = unique(as.character(unlist(CNA_SNV_ixns[, c("driver1", "driver2")])))


fname = pff("CNA_SNV_ixn__oncoprints.pdf")
pdf(fname, width = 5, height = 3)
lapply(all_drivers, get_oncoprint_plot, CNA_SNV_ixns_symm, patient_sets_for_drivers_collapsed)
dev.off()
file_open_call2(fname)


load(file = pff("CNA_SNV_ixns.rsav"))
load(file = pff("patient_sets_for_drivers_collapsed.rsav"))

fname = pff("figures/RGA_RGA_interactions_Venn_diagrams.pdf")
pdf(fname)

# i = 1 
for (i in 1:nrow(CNA_SNV_ixns)) {
	cat(i, "\n")
	
	driver1 = CNA_SNV_ixns[i, "driver1"]
	driver2 = CNA_SNV_ixns[i, "driver2"]
	fdr = signif(CNA_SNV_ixns[i, "fdr"], 2)
	pval = signif(CNA_SNV_ixns[i, "pp"], 2)
	or = signif(CNA_SNV_ixns[i, "or"], 2)
	n_expt = signif(CNA_SNV_ixns[i, "exp_n_common"], 2)
	n_common = CNA_SNV_ixns[i, "n_common_tumors"]
	n_dr1 = CNA_SNV_ixns[i, "n1"]
	n_dr2 = CNA_SNV_ixns[i, "n2"]
	
	driver1_patients = patient_sets_for_drivers_collapsed[[driver1]]
	driver2_patients = patient_sets_for_drivers_collapsed[[driver2]]
	all_drivers = unique(c(driver2_patients, driver1_patients))

	# euler diagram for sites
	dfr_for_venn = data.frame(
			all_drivers, 
			driver1_patients = all_drivers %in% driver1_patients,
			driver2_patients = all_drivers %in% driver2_patients,
			stringsAsFactors = FALSE)
	colnames(dfr_for_venn) = c("all", driver1, driver2)
	
	venn_colors = c("#99CCFF", "#6666CC")
	
	plt_title = paste0("FDR=", fdr, "; OR=", or, "; p=", pval, "\n",
			"n_", driver1, "=", n_dr1, "; n_", driver2, "=", n_dr2, "; n_12=", n_common, "; n_exp=", n_expt )
	
	fit2 = euler(dfr_for_venn[,-1])
	print(plot(fit2, quantities = TRUE, main = plt_title, fills = list(fill = venn_colors)))
}
	
dev.off()
file_open_call2(fname)
