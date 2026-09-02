source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(StructuralVariantAnnotation)
library(ggbio)

load(pff("all_patients.rsav"))
MIN_DRIVER_FREQ = 0.01
MIN_PATIENTS_FOR_GLOBAL = MIN_DRIVER_FREQ * length(all_patients)
MIN_PATIENTS_FOR_GENEBASED = 1

# genomic ranges for chromosomes to decorate the circos plots
hg19_seq_info = keepStandardChromosomes(seqinfo(Hsapiens))
gr_chromosomes = GRanges(seqnames(hg19_seq_info), 
		IRanges(rep(1, length(hg19_seq_info)), seqlengths(hg19_seq_info)), 
		seqinfo = hg19_seq_info)
gr_chromosomes$chr_label = gsub("chr", "", seqnames(gr_chromosomes))
	

# SV drivers
load(file = pff("variants_SV.rsav"))
load(file = pff("prepared_variants_SV.rsav"))
load(file = pff("variants_to_elements.rsav"))
load(file = pff("results_signf_merged_annot.rsav"))
load(file = pff("gr_prepared_elements.rsav"))
load(file = pff("ets_based_drivers.rsav"))

gr_CDSgene = gr_prepared_elements[grep("CDSgene", gr_prepared_elements$id)]
start(gr_CDSgene) = start(gr_CDSgene) - 3000
end(gr_CDSgene) = end(gr_CDSgene) + 3000

variants_to_elements = variants_to_elements[['SV']]
results_SV = results_signf_merged_annot[results_signf_merged_annot$mut_type == "SV",]
prepared_variants = prepared_variants_SV

rm(results_signf_merged_annot)
rm(prepared_variants_SV)
gc()


get_sv_2d_events = function(gene, results_SV, variants_to_elements, variants_SV, min_patients, ets_based_drivers) {
	
	elements = results_SV[results_SV$annots_MAIN == gene, "element_ids"]
	
	# these come from a different dfr because for most of the ms we speak of ETS as a single driver
	ets_elements = ets_based_drivers[ets_based_drivers$annots_MAIN == "ETS", "element_ids"]
	ets_elements = unique(unlist(strsplit(ets_elements, split = ",")))		
	
	# all elements if gene unspecified
	if (is.null(gene)[[1]]) {
		elements = unique(c(results_SV$element_ids, ets_elements))
	}
	elements = unique(unlist(strsplit(elements, split = ",")))
	
	if (!is.null(gene)[[1]] && gene == "ETS") {
		elements = ets_elements
	}
	
	TRA_here = variants_to_elements[elements]
	TRA_here = lapply(TRA_here, function(x) x[x$top_mut_signt == "TRA",])

	n_patients_per_TRA = sapply(TRA_here, function(x) length(unique(x$patient)))
	TRA_here = TRA_here[n_patients_per_TRA >= min_patients]
	TRA_here = do.call(rbind, TRA_here)

	tag_TRA_here = unique(TRA_here$SV_tag)

	prep_SVs = prepared_variants[prepared_variants$SV_tag %in% tag_TRA_here, ]
	prep_SVs = prep_SVs[!duplicated(prep_SVs$SV_tag),]
	
	TRA_here_2d = variants_SV[variants_SV$SV_tag %in% tag_TRA_here,]
	TRA_here_2d
}


get_TRA_genes = function(TRA_here_2d, gr_CDSgene) {

	gr_bp1 = GRanges(TRA_here_2d$chr1, IRanges(TRA_here_2d$start1, TRA_here_2d$end1), seqinfo = hg19_seq_info)
	gr_bp2 = GRanges(TRA_here_2d$chr2, IRanges(TRA_here_2d$start2, TRA_here_2d$end2), seqinfo = hg19_seq_info)
	gr_all_bp = c(gr_bp1, gr_bp2)
	
	gr_tra_genes = gr_CDSgene[unique(queryHits(findOverlaps(gr_CDSgene, gr_all_bp)))]
	
	mean_coord = (start(gr_tra_genes) + end(gr_tra_genes)) / 2
	gr_tra_genes_mean = GRanges(seqnames(gr_tra_genes), IRanges(mean_coord, mean_coord), seqinfo = hg19_seq_info, symbol = gr_tra_genes$symbol)
	gr_tra_genes_mean
}



plot_SVs = function(SVs_here, gene, gr_genes_here, MIN_PATIENTS) {

	gr_bp1 = GRanges(SVs_here$chr1, IRanges(SVs_here$start1, SVs_here$end1), seqinfo = hg19_seq_info)
	gr_bp2 = GRanges(SVs_here$chr2, IRanges(SVs_here$start2, SVs_here$end2), seqinfo = hg19_seq_info)	
	mcols(gr_bp1)$to.gr = gr_bp2
	
	n_patients = length(unique(SVs_here$patient))
	n_SVs = nrow(SVs_here)
	plot_title = gene
	plot_subtitle = paste0("n_patients=", n_patients, "; n_SVs=", n_SVs, "; min_patn=", MIN_PATIENTS)
	
	ggplt = ggbio() + 
			circle(gr_bp1, geom = "link", linked.to = "to.gr", alpha = 0.1, size = 0.5) +
			circle(gr_chromosomes, geom = 'text', aes(label = chr_label), vjust = 0, size = 1) + 
			circle(gr_chromosomes, geom = 'scale', size = 2) +
			theme(plot.subtitle = element_text(size = 6)) +
			theme(axis.text.x = element_text(angle = 90)) +
			ggtitle(plot_title, plot_subtitle)

	# some plots have no genes annotated
	if (length(gr_genes_here) > 0) { 
			ggplt = ggplt + 
				circle(gr_genes_here, geom = 'text', aes(label = symbol, color = gene_type), vjust = 0, size = 1.5) +
				circle(gr_genes_here, geom = "point", aes(fill = gene_type), shape = 21, size = 2, color = "black")
	}

			
	return(ggplt)
}


all_genes = unique(results_SV$annots_MAIN)
# red labels are for the genes associated with the SV-drivers; add ETS labels from the other data frame
driver_genes_to_label = c(unique(ets_based_drivers$annots_flank_CGC), all_genes)


# all SVs combined comes first; 
# two plots; at least 10 samples each; or all 
combined_SVs_here_min10 = get_sv_2d_events(NULL, results_SV, variants_to_elements, variants_SV, MIN_PATIENTS_FOR_GLOBAL, ets_based_drivers)
gr_combined_genes_here_min10 = get_TRA_genes(combined_SVs_here_min10, gr_CDSgene)
gr_combined_genes_here_min10$gene_type = "other"
gr_combined_genes_here_min10$gene_type [gr_combined_genes_here_min10$symbol %in% driver_genes_to_label] = "driver"
# for the combined plot, remove the non-driver genes
gr_combined_genes_here_min10 = gr_combined_genes_here_min10[gr_combined_genes_here_min10$gene_type == "driver"]

# 1 samples each
combined_SVs_here_min1 = get_sv_2d_events(NULL, results_SV, variants_to_elements, variants_SV, MIN_PATIENTS_FOR_GENEBASED, ets_based_drivers)
gr_combined_genes_here_min1 = get_TRA_genes(combined_SVs_here_min1, gr_CDSgene)
gr_combined_genes_here_min1$gene_type = "other"
gr_combined_genes_here_min1$gene_type [gr_combined_genes_here_min1$symbol %in% driver_genes_to_label] = "driver"
# for the combined plot, remove the non-driver genes
gr_combined_genes_here_min1 = gr_combined_genes_here_min1[gr_combined_genes_here_min1$gene_type == "driver"]


fname = pff("figures/Circos_plot_translocations.pdf")
pdf(fname)

print(plot_SVs(combined_SVs_here_min10, paste0("all combined_min", MIN_PATIENTS_FOR_GLOBAL), gr_combined_genes_here_min10, MIN_PATIENTS_FOR_GLOBAL))
print(plot_SVs(combined_SVs_here_min1, paste0("all combined", MIN_PATIENTS_FOR_GENEBASED), gr_combined_genes_here_min1, MIN_PATIENTS_FOR_GENEBASED))

# genes in main annots one by one
for (gene in all_genes) {

	# all SVs
	cat(gene, " ")
	SVs_here = get_sv_2d_events(gene, results_SV, variants_to_elements, variants_SV, MIN_PATIENTS_FOR_GENEBASED, ets_based_drivers)
	if (nrow(SVs_here) < 1) {
		cat(" skipping..")
		next
	}
	
	gr_genes_here = get_TRA_genes(SVs_here, gr_CDSgene)
	if (length(gr_genes_here) > 0 ) { 
		gr_genes_here$gene_type = "other"
		gr_genes_here$gene_type [gr_genes_here$symbol %in% driver_genes_to_label] = "driver"
	}
	
	print(plot_SVs(SVs_here, gene, gr_genes_here, MIN_PATIENTS_FOR_GENEBASED))
	cat("\n")
}
dev.off()
file_open_call2(fname)
