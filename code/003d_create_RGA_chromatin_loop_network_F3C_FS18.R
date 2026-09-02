source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")

load(pff("results_signf_merged_annot.rsav"))
results_signf_merged_annot = results_signf_merged_annot[
		!is.na(results_signf_merged_annot$annots_loops) & 
		results_signf_merged_annot$annots_loops != "",]
# exclude protein-coding drivers, ETS SVs
results_signf_merged_annot = results_signf_merged_annot[results_signf_merged_annot$mut_type != "SNV_CDS",]
results_signf_merged_annot = results_signf_merged_annot[results_signf_merged_annot$annots_MAIN != "ETS",]

# chromatin loop annotations of non-coding sites and genes
get_loops_from_results = function(element, results_signf_merged_annot) { 
	
	res_here = results_signf_merged_annot[results_signf_merged_annot$annots_MAIN == element, ]
	loop_tgts = setdiff(res_here$annots_loops, c(NA, ""))
	loop_tgts = unique(strsplit(loop_tgts, split = ",")[[1]])
	mut_type = paste(sort(unique(res_here$mut_type)), collapse = ",")
	data.frame(driver = element, loop_tgts, mut_type, stringsAsFactors = FALSE)
	
}

main_genes = results_signf_merged_annot$annots_MAIN
loops_combined = do.call(rbind, lapply(main_genes, get_loops_from_results, results_signf_merged_annot))
loops_combined = unique(loops_combined)
n_patients_per_element = split(results_signf_merged_annot$patient_ids, results_signf_merged_annot$annots_MAIN)
n_patients_per_element = sapply(n_patients_per_element, function(x) length(unique(unlist(strsplit(x, s = ",")) )))

fname = pff("figures/FMRE_loops_network_edges.txt")
write.table(loops_combined, file = fname, quote = FALSE, sep = "\t", row.names = FALSE)
file_open_call2(fname)
FMRE_loops_network_edges = loops_combined
save(FMRE_loops_network_edges, file = pff("FMRE_loops_network_edges.rsav"))


all_elements = unique(unlist(loops_combined[,c(1, 2)]))
elements_loops_combined = data.frame(el = all_elements, n_patients = n_patients_per_element[all_elements], stringsAsFactors = FALSE)
elements_loops_combined$is_driver = elements_loops_combined$el %in% main_genes
fname = pff("figures/FMRE_loops_network_nodes.txt")
write.table(elements_loops_combined, file = fname, quote = FALSE, sep = "\t", row.names = FALSE)
file_open_call2(fname)
FMRE_loops_network_nodes = elements_loops_combined
save(FMRE_loops_network_nodes, file = pff("FMRE_loops_network_nodes.rsav"))


