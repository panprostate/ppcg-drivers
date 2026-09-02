source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(igraph)
library(patchwork)

load(pff("all_patients.rsav"))
MIN_DRIVER_FREQ = 0.01
MIN_N_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)
MIN_N_PATIENTS
#[1] 9.61

HOTSPOT_PADDING = 1

find_hotspot_variants = function(gr_prepared_variants, padding_bp, N_MIN_PATIENTS) {
	
	gr_padded = GRanges(seqnames(gr_prepared_variants), IRanges(start(gr_prepared_variants) - padding_bp, end(gr_prepared_variants) + padding_bp))
	ov = findOverlaps(gr_padded, gr_padded)
	ov_dfr = as.data.frame(ov)
	
	# model as connected components - this will merge overlapping sets of indices
	graph = graph_from_data_frame(ov_dfr, directed = FALSE)
	connected_components = components(graph)$membership
	groups_of_variants = split(names(connected_components), connected_components)
	groups_of_variants = lapply(groups_of_variants, as.numeric)
	groups_of_variants = groups_of_variants[sapply(groups_of_variants, length) > N_MIN_PATIENTS]
	groups_of_variants
}



# x = variant_index_sets_flank1[[2]]
get_group_variants = function(x, prepared_variants, dnds_results_annotations_all, gr_RGA_elements, gr_prepared_elements, MIN_N_PATIENTS) {
	cat(".")
	var_here = prepared_variants[x, ]
	var_here1 = merge(var_here, dnds_results_annotations_all[, c("chr", "pos", "sampleID", "gene", "impact", "aachange")], 
			by.x = c("chr", "pos1", "patient"), by.y = c("chr", "pos", "sampleID"), all.x = TRUE)

	gr_vars_here1 = GRanges(var_here1$chr, IRanges(var_here1$pos1, var_here1$pos2))
	ov_RGA_elements = findOverlaps(gr_vars_here1, gr_RGA_elements)
	
	var_here1$RGA_ids = NA
	if (length(ov_RGA_elements) > 0) {
		var_here1$RGA_ids = unique(gr_RGA_elements[subjectHits(ov_RGA_elements)]$RGA_id)	
	}
	
	cds_gene = unique(var_here1$gene)
	cds_impact = paste(names(rev(sort(table(var_here1$impact)))), collapse = ",")
	cds_change = paste(names(rev(sort(table(var_here1$aachange)))), collapse = ",")
	noncoding_rga = unique(var_here1$RGA_ids)
	n_unq_patients = length(unique(var_here1$patient))
	coord = paste0( gsub("chr", "", var_here1$chr[1]), ":", min(var_here1$pos1), "-", max(var_here1$pos2))
	annovar_annot = paste(unique(paste(var_here$Gene.refGene, var_here$Func.refGene, sep = " ")), collapse = ",")
	
	# label based on if hotspot was tested as part of ADWGS, or not tested. exclude CDSgene as these are only used for SV analyses
	ov_with_elements = findOverlaps(gr_prepared_elements, gr_vars_here1)
	overlapping_element_ids = unique(gr_prepared_elements[queryHits(ov_with_elements)]$id)
	overlapping_element_ids = grep("^CDSgene", overlapping_element_ids, invert = TRUE, value = TRUE)
	
	global_annot = "not RGA, not tested"
	global_label = annovar_annot
	
	if (length(overlapping_element_ids) > 0) {
		global_annot = "not RGA, tested & not signifcant"	
		global_label = paste(overlapping_element_ids, collapse = ";")
	}
	
	if (!is.na(noncoding_rga)) {
		global_annot = "non-coding RGA"
		global_label = noncoding_rga	
	}
	if (!is.na(cds_gene)) {
		global_annot = "CDS RGA"
		global_label = paste(cds_gene, cds_impact, cds_change, sep = ",")
	}
	
	if (n_unq_patients < MIN_N_PATIENTS) {
		return(NULL)
	}
	
	var_summary = data.frame(global_annot, n_unq_patients, coord, cds_gene, cds_impact, noncoding_rga, annovar_annot, global_label, 
			stringsAsFactors = FALSE)
	
	list(var_here1, var_summary)
}




load(pff("gr_prepared_variants__SNV_indel.rsav"))
load(pff("prepared_variants__SNV_indel.rsav"))

# annotate mutation hotspots with either non-coding RGAs or with SNV/indel findings from DnDs-cv
load(pff("dnds_results_annotations_all.rsav"))
dnds_results_annotations_all$chr = paste0("chr", dnds_results_annotations_all$chr)
gr_dnds = GRanges(dnds_results_annotations_all$chr, IRanges(dnds_results_annotations_all$pos, dnds_results_annotations_all$pos), 
		mcols = dnds_results_annotations_all[, c("gene", "aachange", "impact")])

load(pff("results_signf_merged_annot.rsav"))
annot_RGAs = results_signf_merged_annot[results_signf_merged_annot$mut_type == "SNV_NC",]
RGAs_to_elements = as.matrix(stack(by(annot_RGAs, annot_RGAs$annots_MAIN, function(x) unique(unlist(strsplit(x$element_ids, sp =","))))))
RGAs_to_elements_map = structure(names = RGAs_to_elements[,1], RGAs_to_elements[,2])

load(file = pff("gr_prepared_elements.rsav"))
gr_RGA_elements = gr_prepared_elements[gr_prepared_elements$id %in% RGAs_to_elements[,1]]
gr_RGA_elements$RGA_id = RGAs_to_elements_map[gr_RGA_elements$id]


# use variant overlaps with padding to find overlapping variant sets, merge these using graph connected components
# padding 1 - 1bp upstream and downstream
variant_index_sets_flank1 = find_hotspot_variants(gr_prepared_variants, HOTSPOT_PADDING, MIN_N_PATIENTS)
groups_of_variants_annots = lapply(variant_index_sets_flank1, get_group_variants, 
		prepared_variants, dnds_results_annotations_all, gr_RGA_elements, gr_prepared_elements, MIN_N_PATIENTS)

hotspot_summaries = do.call(rbind, lapply(groups_of_variants_annots, '[[', 2))
hotspot_summaries = hotspot_summaries[order(hotspot_summaries$n_unq_patients, decreasing = TRUE),]
hotspot_variants =  do.call(rbind, lapply(groups_of_variants_annots, '[[', 1))
hotspot_summaries$coord = factor(hotspot_summaries$coord, levels = rev(hotspot_summaries$coord))

colors_for_annotation_types = c(
		"CDS RGA" = "gold", 
		"non-coding RGA" = "darkorange", 
		"not RGA, not tested" = "lightgrey",
		"not RGA, tested & not signifcant" = "darkgrey")


titl = paste0("SNV/indel hotspots, padding = ", HOTSPOT_PADDING, 
		"\n n= ", nrow(hotspot_summaries), "; n_var = ", nrow(hotspot_variants), "; min_patients=", MIN_N_PATIENTS)

fname = pff(paste0("figures/hotspot_SNVs_padding", HOTSPOT_PADDING, ".pdf"))
plt_pad1 = ggplot(hotspot_summaries, aes(coord, n_unq_patients, fill = global_annot, label = global_label)) + 
		geom_bar(stat = "identity", color = "black") + 
		geom_text(size = 2, y = 1, angle = 0, hjust = 0) + 
		plot_theme() +
		coord_flip() + 
		scale_fill_manual (values = colors_for_annotation_types) + 
		theme(legend.position = "bottom") + 
		ggtitle(NULL, titl)


ggsave(plt_pad1, file = fname, width = 9, height = 7)
file_open_call2(fname)


save(hotspot_summaries, file = pff("hotspot_summaries.rsav"))
save(hotspot_variants, file = pff("hotspot_variants.rsav"))

