source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")

# merge unique elements to make sure the same region does not get counted multiple times
merge_elements_to_islands = function(results_signf_this, gr_prepared_elements, tag) {

	# merge elements into element islands with overlap >= 1
	gr_elements = gr_prepared_elements[gr_prepared_elements$id %in% results_signf_this$id]
	element_coverage = coverage(gr_elements)
	
	chromosomes = as.character(unique(seqnames(gr_elements)))
	list_element_islands = lapply(chromosomes, function(chr) slice(element_coverage, lower = 1)[[chr]])
	names(list_element_islands) = chromosomes
	
	gr_element_islands = do.call(c, lapply(chromosomes, function (chr) {
		GRanges(chr, IRanges(start(list_element_islands[[chr]]), end = end(list_element_islands[[chr]])))
	}))
	
	
	el_chromosome = seqnames(gr_element_islands)
	el_starts = start(gr_element_islands)
	el_ends = end(gr_element_islands) + 1
	el_mid = ceiling( (el_starts + el_ends)/2)
	gr_element_islands$id = paste0(tag, "_", el_chromosome, "_", el_mid)

	ov_islands_vs_elements = findOverlaps(gr_element_islands, gr_elements)
	dfr_islands_vs_elements = cbind(
			island_id = gr_element_islands$id[queryHits(ov_islands_vs_elements)],
			element_id = gr_elements$id[subjectHits(ov_islands_vs_elements)])
	
	island_to_element = split(dfr_islands_vs_elements[,"element_id"], dfr_islands_vs_elements[,"island_id"])
	island_to_element = lapply(island_to_element, unique)
	list(gr_element_islands, island_to_element)
}


collapse_annots = function(ids, nc_to_gene, cancer_genes) {
	
	annots = setdiff(unique(unlist(nc_to_gene[ids])), "")
	annots = gsub("(.+)::(.+)::(.+)::(.+)", "\\3", annots)
	if (is.null(cancer_genes[[1]])) {
		cancer_genes = annots
	}
	annots = annots[annots %in% cancer_genes]
	annots = unique(annots)
	paste(annots, collapse = ",")
}


# take merged elements and assign patients, stats and original element IDs to them in a data frame similar to original results
annotate_merged_element = function(island_id, island_to_element, results_signf_noncoding, 
		nc_to_gene_flank, nc_to_gene_loops, cancer_genes_combined, gr_prepared_elements) {
	element_ids = island_to_element[[island_id]]
	results_signf_here = results_signf_noncoding[results_signf_noncoding$id %in% element_ids,]

	fdr = min(results_signf_here$fdr)
	pval = min(results_signf_here$pp_element)
	
	gr_el_combined = gr_prepared_elements[gr_prepared_elements$id %in%  island_to_element[[island_id]]]
	location_chr = as.character(seqnames(gr_el_combined))[[1]]
	location_start = min(start(gr_el_combined))
	location_end = max(end(gr_el_combined))
	location_width = location_end - location_start
	
	patient_ids = sort(unique(unlist(strsplit(results_signf_here$patient_ids_concat, s = ","))))
	element_types = sort(unique(results_signf_here$element_type))
	element_ids = sort(unique(results_signf_here$id))
	ds_types = sort(unique(results_signf_here$ds))
	n_patients = length(patient_ids)
	
	element_ids = paste(element_ids, collapse = ",")
	patient_ids = paste(patient_ids, collapse = ",")
	element_type = paste(element_types, collapse = ",")
	ds = paste(ds_types, collapse = ",")

	annots_flank = collapse_annots(results_signf_here$id, nc_to_gene_flank, NULL)
	annots_loops = collapse_annots(results_signf_here$id, nc_to_gene_loops, NULL)
	annots_flank_CGC = collapse_annots(results_signf_here$id, nc_to_gene_flank, cancer_genes_combined)
	annots_loops_CGC = collapse_annots(results_signf_here$id, nc_to_gene_loops, cancer_genes_combined)

	dfr = data.frame(id = island_id, fdr, pval, patient_ids, element_type, n_patients, ds, element_ids, 
			location_chr, location_start, location_end, location_width,
			annots_flank, annots_flank_CGC, annots_loops, annots_loops_CGC,
			stringsAsFactors = FALSE)
	rownames(dfr) = NULL
	dfr
}

# get protein-coding mutations from dndscv, reformat consistently with others
format_results_from_genes = function(genes_results, standard_colnames, gene_CDS_id_map, dnds_results_annotations, gr_prepared_elements) {
		
	id = gene_CDS_id_map[genes_results$gene_name]
	
	# get gene locus to match other types of elements
	gr_this = gr_prepared_elements[gr_prepared_elements$id %in% id]
	gr_this_split = split(gr_this, gr_this$id)
	this_chr = sapply(gr_this_split, function(x) as.character(seqnames(x))[[1]])
	this_start = sapply(gr_this_split, function(x) min(start(x)))
	this_end = sapply(gr_this_split, function(x) max(end(x)))
	this_width = sapply(gr_this_split, function(x) max(end(x)) - min(start(x)))
	
	fdr = genes_results$qglobal_cv
	pp_element = genes_results$pglobal_cv
	
	patients_here = lapply(split(dnds_results_annotations$sampleID, dnds_results_annotations$gene), unique)
	patient_ids_concat = sapply(patients_here, function(x) paste(sort(x), collapse = ","))
	n_patients = sapply(patients_here, length)
	element_ids = id
	names(patient_ids_concat) = names(n_patients) = element_ids[names(patients_here)]

	merged_genes_results = data.frame(
			id = id, 
			fdr = fdr,
			pp_element = pp_element,
			patient_ids_concat = patient_ids_concat[id],

			element_type = "CDS",
			n_patients = n_patients[id],
			ds = "SNV_indel",
			element_ids = id, 
			location_chr = this_chr[id], 
			location_start = this_start[id],
			location_end = this_end[id], 
			location_width = this_width[id],
			
			annots_flank = genes_results$gene_name,
			annots_flank_CGC = sapply(genes_results$gene_name, function(x) ifelse(x %in% cancer_genes_combined, x, NA)),
			annots_loops = NA,
			annots_loops_CGC = NA,
			
			stringsAsFactors = FALSE
		)	
		
	colnames(merged_genes_results) = standard_colnames
	merged_genes_results
}

load(file = pff("results_signf.rsav"))
load(file = pff("gr_prepared_elements.rsav"))
load(file = pff("nc_to_gene_flank.rsav"))
load(file = pff("nc_to_gene_loops.rsav"))
load(file = pff("prepared_variants__SNV_indel.rsav"))

print(load(pff("cgc2024.rsav")))
print(load(pff("prostate_cancer_genes.rsav")))
cancer_genes_combined = c(cgc2024, prostate_cancer_genes)

results_signf_ncSNV = results_signf[results_signf$ds == "SNV_indel",]
results_signf_SV = results_signf[results_signf$ds == "SV",]

# merge elements that overlap into consecutive regions
ncSNV_merged_islands_results = merge_elements_to_islands(results_signf_ncSNV, gr_prepared_elements, "ncSNV")
SV_merged_islands_results = merge_elements_to_islands(results_signf_SV, gr_prepared_elements, "SV")

# two objects returned from merging: list of merged elements and their mapping to initial elements
gr_element_islands_ncSNV = ncSNV_merged_islands_results[[1]]
island_to_element_ncSNV = ncSNV_merged_islands_results[[2]]
gr_element_islands_SV = SV_merged_islands_results[[1]]
island_to_element_SV = SV_merged_islands_results[[2]]

save(island_to_element_ncSNV, file = pff("island_to_element_ncSNV.rsav"))
save(island_to_element_SV, file = pff("island_to_element_SV.rsav"))

# annotate merged elements to keep the original column structure in results_signf
results_signf_merged_ncSNV = do.call(rbind, lapply(gr_element_islands_ncSNV$id, 
		annotate_merged_element, island_to_element_ncSNV, results_signf_ncSNV, nc_to_gene_flank, nc_to_gene_loops, 
		cancer_genes_combined, gr_prepared_elements))
results_signf_merged_SV = do.call(rbind, lapply(gr_element_islands_SV$id, 
		annotate_merged_element, island_to_element_SV, results_signf_SV, nc_to_gene_flank, nc_to_gene_loops, 
		cancer_genes_combined, gr_prepared_elements))

# make CDS SNV results formatted the same as the rest
# keep only patients that have protein-coding variants, removing silent variants
load(pff("dnds_results.rsav"))
load(pff("dnds_results_annotations.rsav"))

# map between symbols and longer ids
gene_CDS_id_map = unique(grep("CDS::", gr_prepared_elements$id, value = T))
gene_CDS_id_map = structure(gene_CDS_id_map, names = gsub("CDS::gencode::(.+)::(.+)", "\\1", gene_CDS_id_map))

results_signf_merged_cdsSNV = format_results_from_genes(
		dnds_results, colnames(results_signf_merged_ncSNV), gene_CDS_id_map, dnds_results_annotations, gr_prepared_elements)

results_signf_merged_cdsSNV$mut_type = "SNV_CDS"
results_signf_merged_ncSNV$mut_type = "SNV_NC"
results_signf_merged_SV$mut_type = "SV"

results_signf_merged = rbind(results_signf_merged_cdsSNV, results_signf_merged_ncSNV, results_signf_merged_SV)
results_signf_merged = results_signf_merged[order(results_signf_merged$pval),]

# gene annotations using either CGC, or if missing, all flanks, or if missing, the ID
results_signf_merged$annots_MAIN = results_signf_merged$annots_flank_CGC
results_signf_merged[results_signf_merged$annots_MAIN == "", "annots_MAIN"] = 
		results_signf_merged[results_signf_merged$annots_MAIN == "", "annots_flank"]
results_signf_merged[results_signf_merged$annots_MAIN == "", "annots_MAIN"] = 
		results_signf_merged[results_signf_merged$annots_MAIN == "", "id"]

# save file, organise main genes manually
save(results_signf_merged, file = pff("results_signf_merged.rsav"))
fname = pff("figures/results_signf_merged.csv")
write.csv(results_signf_merged, file = fname, row.names = FALSE)
file_open_call2(fname)

## manual annotation of candidate drivers based on main gene
## load the annotations from the CSV file that was prepared
results_signf_merged_annot = read.csv(file = pff("figures/results_signf_merged_EDIT.csv"), stringsAsFactors = FALSE)
results_signf_merged_annot = results_signf_merged_annot[order(results_signf_merged_annot$fdr),]

# remove ETS SVs from the SV results; handled separately based on the consensus dataset
load(file = pff("ETS_target_genes.rsav"))
load(file = pff("ETS_patients.rsav"))
ETS_target_genes = c(ETS_target_genes, "ETS")
index_ETS_results = results_signf_merged_annot$annots_MAIN %in% ETS_target_genes & results_signf_merged_annot$mut_type == "SV"
ets_based_drivers = results_signf_merged_annot[index_ETS_results,]
results_signf_merged_annot = results_signf_merged_annot[!index_ETS_results,]
save(ets_based_drivers, file = pff("ets_based_drivers.rsav"))

ETS_dfr = data.frame(
		id = "ETS",
		fdr = 1e-300,
		pval = 1e-300,
		patient_ids = paste(ETS_patients, collapse = ","),
		element_type = NA,
		n_patients = length(ETS_patients),
		ds = "SV",
		element_ids = NA, 
		location_chr = NA, 
		location_start = NA,
		location_end = NA, 
		location_width = NA,
		annots_flank = "ETS",
		annots_flank_CGC = "ETS",
		annots_loops = "ETS",
		annots_loops_CGC = "ETS",
		mut_type = "SV",
		annots_MAIN = "ETS",
		stringsAsFactors = FALSE
	)
results_signf_merged_annot = rbind(ETS_dfr, results_signf_merged_annot)	
save(results_signf_merged_annot, file = pff("results_signf_merged_annot.rsav"))
