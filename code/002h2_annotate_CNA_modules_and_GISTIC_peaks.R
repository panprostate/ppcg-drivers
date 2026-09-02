source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
load(pff("CNA_modules.rsav"))

load(file = pff("gr_prepared_elements.rsav"))
gr_genes = gr_prepared_elements[grep("^CDSgene", gr_prepared_elements$id)]

load(pff("cgc2024.rsav")))
load(pff("prostate_cancer_genes.rsav"))
cancer_genes_combined = unique(c(cgc2024, prostate_cancer_genes))



#add coordinates
add_chr_location_to_module = function(i, CNA_modules, gr_genes) {
	
	genes_here = strsplit(CNA_modules[i, "all_genes"], split = ",")[[1]]
	gr_genes_here = gr_genes[gr_genes$symbol %in% genes_here]
	
	chr = as.character(seqnames(gr_genes_here)[1])
	start = min(start(gr_genes_here)) - 1e6
	end = max(end(gr_genes_here)) + 1e6
	location_str = paste(chr, ":", start, "-", end, sep = "")
	id = paste0("module", "__", location_str, "__", CNA_modules[i, "annot"])
	cat(id, "\n")
	
	cbind(id, chr, start, end, CNA_modules[i,, drop = FALSE], stringsAsFactors = FALSE)
}

CNA_modules_annotated = do.call(rbind, lapply(1:nrow(CNA_modules), add_chr_location_to_module, CNA_modules, gr_genes))
save(CNA_modules_annotated, file = pff("CNA_modules_annotated.rsav"))

fname = pff("figures/CNA_modules_annotated.csv")
write.csv(CNA_modules_annotated, file = fname)
file_open_call2(fname)




# add genes to gistic peaks using consistent formatting with CNA modules
add_genes_to_gistic_peak = function(i, significant_CNAs_GISTIC_with_segs, gr_genes, gene_CNA_enriched, cgc2024, prostate_cancer_genes) {
	
	peak_here = significant_CNAs_GISTIC_with_segs[i,]

	gr_peak = GRanges(peak_here$chr, IRanges(peak_here$start, peak_here$end))
	ov = findOverlaps(gr_peak, gr_genes)
	gr_genes_here = gr_genes[subjectHits(ov)]
	
	genes_here = gr_genes_here$symbol
	
	# select only the genes enriched in matching CNA types - gains or losses
	annot_here = peak_here$CNA_type
	gene_annots_to_select = new_annot = NA
	if (annot_here == "Amp") {
		gene_annots_to_select = c("all_gain", "high_gain")
		new_annot = "all_gain"
	} 
	if (annot_here == "Del") {
		gene_annots_to_select = c("all_loss", "full_loss")
		new_annot = "all_loss"
	}
	if (any(is.na(gene_annots_to_select))) {
		stop("invalid annot")
	}
	
	dfr_genes_enriched_here = gene_CNA_enriched[
			gene_CNA_enriched$gene %in% genes_here & 
			gene_CNA_enriched$fdr < 0.05 & 
			gene_CNA_enriched$annot %in% gene_annots_to_select, ]
	
	genes_here_by_pval = names(sort(by(dfr_genes_enriched_here$pval, dfr_genes_enriched_here$gene, min)))
	all_genes = genes_here_by_pval
	cgc_genes = intersect(genes_here_by_pval, cgc2024)
	prca_genes = intersect(genes_here_by_pval, prostate_cancer_genes)
	TGT_gene = prca_genes
	
	id = paste0("peak", "__", gsub("(.+)_(.+)", "\\1", peak_here$id), "__", peak_here$CNA_type)
	cat(id, "\n")
	
	# pval column is used for the value of FDR_residual
	data.frame(id, chr = peak_here$chr, start = peak_here$start, end = peak_here$end, gene = paste(all_genes, collapse = ","), 
			pval = peak_here$fdr_residual, frac_CNA = NA, fc = NA, n_patients = peak_here$n_patients_segs, patients = peak_here$patients_segs, 
			annot = new_annot, fdr = peak_here$fdr, 
			all_genes = paste(all_genes, collapse = ","), cgc_genes = paste(cgc_genes, collapse = ","), 
			prca_genes = paste(prca_genes, collapse = ","), TGT_gene = paste(prca_genes, collapse = ","), stringsAsFactors = FALSE)
}

load(pff("significant_CNAs_GISTIC_with_segs.rsav"))
load(pff("gene_CNA_enriched.rsav"))

GISTIC_peaks_annotated = do.call(rbind, lapply(1:nrow(significant_CNAs_GISTIC_with_segs), 
		add_genes_to_gistic_peak, significant_CNAs_GISTIC_with_segs, gr_genes, gene_CNA_enriched, cgc2024, prostate_cancer_genes))
save(GISTIC_peaks_annotated, file = pff("GISTIC_peaks_annotated.rsav"))

fname = pff("figures/GISTIC_peaks_annotated.csv")
write.csv(GISTIC_peaks_annotated, file = fname)
file_open_call2(fname)

fname = pff("figures/CNA_modules_annotated_EDIT.csv")
CNA_modules_annotated = read.csv(file = fname, stringsAsFactors = FALSE)[,-1]
CNA_modules_w_tgt_genes = CNA_modules_annotated
save(CNA_modules_w_tgt_genes, file = pff("CNA_modules_w_tgt_genes.rsav"))

fname = pff("figures/GISTIC_peaks_annotated_EDIT.csv")
GISTIC_peaks_annotated = read.csv(file = fname, stringsAsFactors = FALSE)[,-1]
GISTIC_peaks_w_tgt_genes = GISTIC_peaks_annotated
save(GISTIC_peaks_w_tgt_genes, file = pff("GISTIC_peaks_w_tgt_genes.rsav"))
