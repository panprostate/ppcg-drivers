source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")

FDR_CUTOFF = 0.05

system(paste("mkdir -p", pff("figures")))
load(pff("cgc2024.rsav"))
load(pff("prostate_cancer_genes.rsav"))
cancer_genes_combined = unique(c(cgc2024, prostate_cancer_genes))

# what tasks were required to do initially
load(file = pff(paste0("ADWGS_tmp/todo_list.0.rsav")))
all_todo_tags = apply(todo_list[,c("element", "muts_dataset")], 1, paste, collapse = "__")

# what tasks have been completed
fnames = list.files(pff("ADWGS_tmp"), pattern = "^result_iter_", full.names = TRUE)
results = unique(do.call(rbind, mclapply(fnames, function(fname) {load(fname); cat("."); result }, mc.cores = 16)))
all_done_tags = apply(results[,c("id", "ds")], 1, paste, collapse = "__")

# all tasks have been completed
all_analyzed = all(all_todo_tags %in% all_done_tags)
if (!all_analyzed) (stop("some results not analyzed yet\n\n"))

# none of the results are duplicated
no_dups = !any(duplicated(all_done_tags))
if (!no_dups) (stop("some duplicates in results\n\n"))

# add silent mutations for which no analysis was performed 
load(file = pff("ADWGS_tmp/todo_list_silent.rsav"))
ADWGS_silent = data.frame(
		id = todo_list_silent[,"element"], pp_element = 1, 
		element_muts_obs = NA, element_muts_exp = NA, 
		element_enriched = NA, pp_site = NA, 
		site_muts_obs = NA, site_muts_exp = NA, 
		site_enriched = NA, 
		ds = todo_list_silent[, "muts_dataset"],
		row_id = (1:nrow(todo_list_silent)) + nrow(results),
		stringsAsFactors = FALSE)

results = rbind(results, ADWGS_silent)
results[is.na(results$pp_element), "pp_element"] = 1
rownames(results) = NULL

results$element_type = gsub("(.+)::(.+)::(.+)::(.+)", "\\1", results$id)
results$symbol = gsub("(.+)::(.+)::(.+)::(.+)", "\\1::\\3", results$id)
results$symbol = gsub("CDS::|lncRNA::|CDSgene::", "", results$symbol)

results$logFC = log1p(results$element_muts_obs) - log1p(results$element_muts_exp)
results$is_cancer_gene = results$symbol %in% cgc2024
results$is_prca_gene = results$symbol %in% prostate_cancer_genes

# one multiple testing correction for all elements; CDS correction comes from dNdScv
results$fdr = p.adjust(results$pp_element, method = "fdr")
results = results[order(results$pp_element),]
save(results, file = pff("results_ADWGS.rsav"))

results_signf = results[results$fdr < FDR_CUTOFF,]

# add gene annotations to non-coding sites, all genes and cancer genes
annot_site_to_gene = function(results_signf, annots, nth_element, filter_cgc = NULL) {
	
	this_annots = annots[results_signf$id]
	this_annots = lapply(this_annots, function(x) {
		if (is.null(x)) return(x)
		sapply(strsplit(x, s = "::"), '[[', nth_element)
	})
	
	this_annots = lapply(this_annots, unique)
	
	if (!is.null(filter_cgc)) { 
		this_annots = lapply(this_annots, function(x) intersect(x, filter_cgc))
	}
	
	this_annots = sapply(this_annots, paste, collapse = ",")
	this_annots
}


load(pff("nc_to_gene_flank.rsav"))
load(pff("nc_to_gene_loops.rsav"))

annots_flank = nc_to_gene_flank[results_signf$id]
results_signf$annots_flank = annot_site_to_gene(results_signf, annots_flank, 3)
results_signf$annots_flank_CGC = annot_site_to_gene(results_signf, annots_flank, 3, cancer_genes_combined)

annots_loops = nc_to_gene_loops[results_signf$id]
results_signf$annots_loops = annot_site_to_gene(results_signf, annots_loops, 3)
results_signf$annots_loops_CGC = annot_site_to_gene(results_signf, annots_loops, 3, cancer_genes_combined)


# list of tumour samples concatenated for each element
load(file = pff("variants_to_elements.rsav"))

altered_patients = lapply(1:nrow(results_signf), function(i) 
		unique(variants_to_elements[[results_signf$ds[i]]][[results_signf$id[[i]]]]$patient))
		
results_signf$patient_ids_concat = sapply(altered_patients, paste, collapse = ",")
rownames(results_signf) = NULL

save(results_signf, file = pff("results_signf.rsav"))