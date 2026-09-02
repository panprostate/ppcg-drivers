source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(scales)
library(patchwork)

load(pff("cgc2024.rsav"))
load(pff("prostate_cancer_genes.rsav"))
cancer_genes = unique(cgc2024, prostate_cancer_genes)


build_pathway_oncoprint = function(i, APW_results, drivers_mat, mut_type_colors) {
	
	cat(i, "/", nrow(APW_results), "\n")
	this_pw_res = APW_results[i,]
	cols_to_get_genes = c("overlap", grep("Genes_", colnames(APW_results), value = T))
	
	pw_genes = setdiff(unique(unlist(this_pw_res[, cols_to_get_genes])), NA)
	pw_id = this_pw_res[1, "term_id"][[1]][[1]]
	pw_name = this_pw_res[1, "term_name"][[1]][[1]]
	pw_fdr = signif(this_pw_res[1, "adjusted_p_val"][[1]][[1]], 2)
	
	this_drivers_mat = drivers_mat[drivers_mat$gene %in% pw_genes,]
	
	# use the oncoprint ordering of rows and columns
	# convert to numeric matrix for sorting
	this_drivers_mat1 = memoSort( this_drivers_mat[,-1] != "")
	patient_order = colnames(this_drivers_mat1)
	
	gene_order = apply(this_drivers_mat[,-1], 1, function(x) sum(x!=""))
	gene_order = this_drivers_mat[,1][order(gene_order)]
	
	this_stack = data.frame(as.matrix(melt(this_drivers_mat, id.var = "gene")), stringsAsFactors = FALSE)
	colnames(this_stack) = c("gene", "patient", "mut")
	this_stack$patient = factor(this_stack$patient, levels = patient_order)
	this_stack$gene = factor(this_stack$gene, levels = gene_order)
	this_stack = this_stack[this_stack$mut != "",]
	n_pat = length(unique(this_stack$patient))
	
	plt_oncoprint = ggplot(this_stack, aes(patient, gene, fill = mut)) + 
			geom_tile() + 
			scale_x_discrete("patients", breaks = NULL, labels = NULL) +
			scale_fill_manual(values = mut_type_colors) + 
			plot_theme() + 
			theme(legend.position = "none", plot.subtitle = element_text(size = 12)) + 
			ggtitle(NULL, paste0(pw_id, "\nFDR=", pw_fdr, "; n_spl=", n_pat, "\n", pw_name))
			
	plot_mutcounts = ggplot(this_stack, aes(gene, fill = mut)) + 
			geom_bar() + 
			coord_flip() + 
			scale_x_discrete(NULL, labels = NULL) + 
			scale_fill_manual(values = mut_type_colors) + 
			plot_theme() 
			
	plt_combined = plt_oncoprint | plot_mutcounts +	plot_layout(widths = unit(c(3, 1), 'null'))

	return(plt_combined)
}	


## gene-based visuals for evidence-collapsed pathway analysis
load(file = pff("APW_results_integrative.rsav"))
APW_results_integrative = APW_results_integrative[order(APW_results_integrative$adjusted_p_val),]

load(file = pff("driver_patient_sets_for_pathway_analysis.rsav"))
drivers_mat = dcast(driver_patient_sets_for_pathway_analysis, gene ~ patient, 
		value.var = "mut_type", fun.aggregate = function(x) paste(sort(x), collapse = ","))

# collapse some mutation types using a spreadsheet annotation
pw_mut_events_tab = table(unlist(drivers_mat[,-1]))
pw_mut_events_tab = data.frame(mut_type = names(pw_mut_events_tab), n_events = c(pw_mut_events_tab), ANNOT = names(pw_mut_events_tab), 
		stringsAsFactors = FALSE, row.names = NULL)

fname = pff("figures/Pathway_mutation_events___APW_results_integrative.csv")
write.csv(pw_mut_events_tab, file = fname)
file_open_call2(fname)

# annotate events in a csv file
# load it in again
pw_mut_events_tab = read.csv(pff("figures/Pathway_mutation_events___APW_results_integrative_EDIT.csv"), stringsAsFactors = FALSE)

# replace the new mutation types into the matrix
for (i in 1:nrow(pw_mut_events_tab)) {
	drivers_mat[,-1][drivers_mat[,-1] == pw_mut_events_tab[i, "mut_type"]] = pw_mut_events_tab[i, "ANNOT"]
}


mut_type_colors = c(
	"gain" = "hotpink",
	"hAMP" = "darkred",
	"loss" = "steelblue",
	"BI" = "darkblue",
	"SNV_CDS" = "gold",
	"SNV_NC" = "darkorange",
	"SV" = "darkgreen"
)
mut_type_colors_with_MULTIPLE = c(mut_type_colors, "MULTIPLE" = "grey24")
	

fname = pff("figures/pathways_with_drivers_integrative__oncoprints.pdf")
pdf(fname, height = 3.5, width = 6.5)

lapply(1:nrow(APW_results_integrative), build_pathway_oncoprint, APW_results_integrative, drivers_mat, mut_type_colors)

dev.off()
file_open_call2(fname)



# pathway-level oncoprint
get_samples_genes_per_pathway = function (i, APW_res, drivers_mat, cols_to_get_genes) {
	cat(i," ")
	
	this_APW_result = APW_res[i, ]
	genes_here = setdiff(unlist(this_APW_result[, cols_to_get_genes]), c("", NA))
	drivers_mat_here = drivers_mat[drivers_mat$gene %in% genes_here, ]
	
	# return samples per pathway as result 1
	drivers_mat_collapsed = apply(drivers_mat_here[,-1], 2, function(x) paste(setdiff(unique(x), ""), collapse = ","))
	drivers_mat_collapsed[grep(",", drivers_mat_collapsed)] = "MULTIPLE"
	result_samples = data.frame(cbind(this_APW_result, t(drivers_mat_collapsed)), stringsAsFactors = FALSE)
	
	# return gene counts per pathway as result 2
	gene_counts = sapply(drivers_mat_here$gene, function(x) sum(drivers_mat_here[drivers_mat_here$gene == x, -1] != ""))
	result_counts = data.frame(pathway = this_APW_result$term_name, gene = names(gene_counts), n_patients = gene_counts, stringsAsFactors = FALSE)
	
	list(result_samples, result_counts)
}


# take manually selected pathways, print these into an oncoprint plot
selected_APW_result = read.csv(pff("figures/APW_results_integrative_EDIT.csv"), stringsAsFactors = FALSE)
selected_APW_result = selected_APW_result[selected_APW_result$SELECT == "YES",]

load(pff("APW_results_integrative.rsav"))
APW_results_integrative = APW_results_integrative[APW_results_integrative$term_id %in% selected_APW_result$term_id, ]
cols_to_get_genes = c("overlap", grep("Genes_", colnames(APW_results_integrative), value = TRUE))


# gene counts and patient sets are now derived using the same function
APW_genes_patients__res_list = 
		lapply(1:nrow(APW_results_integrative), get_samples_genes_per_pathway, APW_results_integrative, drivers_mat, cols_to_get_genes)
selected_APW_result_with_samples = do.call(rbind, lapply(APW_genes_patients__res_list, '[[', 1))
selected_APW_result_gene_counts = do.call(rbind, lapply(APW_genes_patients__res_list, '[[', 2))

# order pathways by most to least frequent
pw_sample_matrix = selected_APW_result_with_samples[,grep("PPCG", colnames(selected_APW_result_with_samples))]
pw_sample_matrix = cbind(term_name = selected_APW_result_with_samples[,2], pw_sample_matrix)
pw_ordered = pw_sample_matrix[order(apply(pw_sample_matrix[,-1], 1, function(x) sum(x != ""))), "term_name"]
pw_sample_matrix_memosort = memoSort( pw_sample_matrix[,-1] != "")
patient_ordered = colnames(pw_sample_matrix_memosort)

this_stack = data.frame(as.matrix(melt(pw_sample_matrix, id.var = "term_name")), stringsAsFactors = FALSE)
colnames(this_stack) = c("pathway", "patient", "mut")
this_stack$patient = factor(this_stack$patient, levels = patient_ordered)
this_stack$pathway = factor(this_stack$pathway, levels = pw_ordered)
this_stack = this_stack[this_stack$mut != "",]
n_pat = length(unique(this_stack$patient))

# gene counts from most to least mutated
selected_APW_result_gene_counts$pathway = factor(selected_APW_result_gene_counts$pathway, levels = pw_ordered)
gene_ordered = rev(names(sort(by(selected_APW_result_gene_counts$n_patients, selected_APW_result_gene_counts$gene, min))))
selected_APW_result_gene_counts$gene = factor(selected_APW_result_gene_counts$gene, levels = gene_ordered)

plt_oncoprint = ggplot(this_stack, aes(patient, pathway, fill = mut)) + 
		geom_tile() + 
		scale_x_discrete("patients", breaks = NULL, labels = NULL) +
		scale_fill_manual(values = mut_type_colors_with_MULTIPLE) + 
		plot_theme() + 
		theme(legend.position = "bottom", plot.subtitle = element_text(size = 12)) + 
		ggtitle(NULL, paste0("n_spl = ", n_pat))
		
plt_mut_stack = ggplot(this_stack, aes(pathway, fill = mut)) + 
		geom_bar() + 
		scale_x_discrete(NULL, labels = NULL) + 
		scale_fill_manual(values = mut_type_colors_with_MULTIPLE) + 
		coord_flip() + 
		plot_theme() +
		theme(legend.position = "none", plot.subtitle = element_text(size = 12))
		
plt_gene_dots = ggplot(selected_APW_result_gene_counts, aes(gene, pathway, fill = n_patients)) + 
		scale_y_discrete(NULL, breaks = NULL, labels = NULL) +
		geom_tile() + 
		scale_fill_gradient( trans = "log10", low = "pink", high = "darkred") + 
		plot_theme() +
		theme(legend.position = "bottom", plot.subtitle = element_text(size = 12))
		

plt_combined = plt_oncoprint + plt_mut_stack + plt_gene_dots + plot_layout(widths = unit(c(3,1,4), 'null'))
fname = pff("key_pathways_oncoprint__APW_results_integrative.pdf")
ggsave(plt_combined, file = fname, width = 12, height = 5)
file_open_call2(fname)

samples_mutated_in_main_pathways = this_stack
save(samples_mutated_in_main_pathways, file = pff("samples_mutated_in_main_pathways__APW_results_integrative.rsav"))
