source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library("patchwork")

# tornado plots of gistic peaks and CNA modules
load(pff("gr_CDS.rsav"))
load(pff("gr_gene_coords.rsav"))
load(pff("gene_CNA_enriched.rsav"))
# make the IDs for CDS and genes the same to make sure the gene_plot function workds
gr_CDS$id = gr_CDS$id2
gr_gene_coords$id = gr_gene_coords$id2

print(load(pff("cgc2024.rsav")))
print(load(pff("prostate_cancer_genes.rsav")))
cancer_genes_combined = unique(c(cgc2024, prostate_cancer_genes))

load(file = pff("all_patients.rsav"))
load(file = pff("prepared_CNAs.rsav"))
gr_CNA_segs = GRanges(prepared_CNAs$chr, IRanges(prepared_CNAs$start, prepared_CNAs$end), 
		patient = prepared_CNAs$patient, annot = prepared_CNAs$annot)
gr_CNA_segs$cn_width = width(gr_CNA_segs)

cna_colors = c("bal" = "white", "gain" = "orange", "high_gain" = "darkred", "loss" = "darkcyan", "full_loss" = "darkblue")


plot_tornado = function(id, CNA_modules_annotated, gr_CNA_segs, all_patients_here, 
		gr_genes, gr_cds, analysis_tag, cna_colors, gene_CNA_enriched, canvas_flank, add_missing_patients) {
	
	cat(id, "\n")		
	this_record = CNA_modules_annotated[CNA_modules_annotated$id == id,]
	CNA_type = this_record$annot
	
	genes_CNAenriched_dfr = gene_CNA_enriched[gene_CNA_enriched$annot == CNA_type & gene_CNA_enriched$fdr < 0.05, ]
	genes_CNAenriched = unique(genes_CNAenriched_dfr$gene)
	all_genes = intersect(strsplit(this_record$all_genes, split = ",")[[1]], genes_CNAenriched)
	cgc_genes = intersect(strsplit(this_record$cgc_genes, split = ",")[[1]], genes_CNAenriched)
	prca_genes = intersect(strsplit(this_record$prca_genes, split = ",")[[1]], genes_CNAenriched)

	# select a subset of genes to be shown as coordinates and counts: those with top signal and known cancer genes
	genes_as_coordinates = genes_as_counts = unique(c(head(all_genes, n = 10), cgc_genes, prca_genes))
	gr_canvas = GRanges(this_record$chr, IRanges(this_record$start - canvas_flank, this_record$end + canvas_flank))
	
	# indicate region with most significant genes
	gr_genes_top = gr_genes[gr_genes$symbol %in% head(all_genes, n = 10)]
	start_min = min(start(gr_genes_top))
	end_min = max(end(gr_genes_top))
	if (length(genes_as_coordinates) == 0) {
		start_min = start(gr_canvas)
		end_min = end(gr_canvas)	
	}
				
	# take the segments overlapping this cytoband, exclude balanced segments
	ov = findOverlaps(gr_canvas, gr_CNA_segs)
	gr_CNA_segs_here = gr_CNA_segs[subjectHits(ov)]
	gr_CNA_segs_here = gr_CNA_segs_here[gr_CNA_segs_here$patient %in% all_patients_here]
	
	# select only gains or only losses based on what Gistic thought
	if (CNA_type == "all_loss") {
		gr_CNA_segs_here = gr_CNA_segs_here[gr_CNA_segs_here$annot %in% c("loss", "full_loss")]
	} else if (CNA_type == "all_gain") {
		gr_CNA_segs_here = gr_CNA_segs_here[gr_CNA_segs_here$annot %in% c("gain", "high_gain")]
	} else {
		gr_CNA_segs_here = gr_CNA_segs_here[gr_CNA_segs_here$annot == CNA_type]
	}
	
	# clip coordinates outside of the cytoband range
	start(gr_CNA_segs_here)[start(gr_CNA_segs_here) < start(gr_canvas)] = start(gr_canvas)
	end(gr_CNA_segs_here)[end(gr_CNA_segs_here) > end(gr_canvas)] = end(gr_canvas)
	
	gr_CNA_segs_here$cn_width_here = end(gr_CNA_segs_here) - start(gr_CNA_segs_here)
	
	# order patients by segment width
	patient_order = names(sort(by(gr_CNA_segs_here$cn_width_here, gr_CNA_segs_here$patient, sum)))
	
	dfr_CNA_segs_here = as.data.frame(gr_CNA_segs_here, stringsAsFactors = FALSE)
	dfr_CNA_segs_here$patient_to_num = as.numeric(factor(dfr_CNA_segs_here$patient, levels = patient_order))
	
	if (add_missing_patients) {
		# patients add patients with balanced CN at this locus 
		missing_patients = setdiff(all_patients_here, patient_order)
		placeholder_CNAs_for_bal_patients = data.frame(
				seqnames = seqnames(gr_canvas),
				start = start(gr_canvas), 
				end = end(gr_canvas), 
				width = width(gr_canvas),
				strand = "*", 
				patient = missing_patients,
				annot = "bal", 
				cn_width = width(gr_canvas),
				patient_to_num = as.numeric(factor(missing_patients)) + max(dfr_CNA_segs_here$patient_to_num), 
				stringsAsFactors = FALSE)
		dfr_CNA_segs_here = rbind(dfr_CNA_segs_here, placeholder_CNAs_for_bal_patients)
	}

	patients_in_record = strsplit(this_record$patients, s = ",")[[1]]
	patients_with_segments = patient_order
	n_pat_rec = length(patients_in_record)
	n_pat_seg = length(patients_with_segments)
	
	plot_title = paste(id, ":", analysis_tag, this_record$TGT_GENE, "\n",
					paste("n_pat_rec=", n_pat_rec, ";\nn_pat_seg=", n_pat_seg, "\n",
					paste("all:", paste(all_genes, collapse = ","), "\n"),
					paste("cgc:", paste(cgc_genes, collapse = ","), "\n"),
					paste("prca:", paste(prca_genes, collapse = ","), "\n")))
						
	plt_segments = ggplot(dfr_CNA_segs_here, 
				aes(xmin = start, xmax = end, fill = annot, color = annot, ymin = patient_to_num, ymax = patient_to_num + 1)) + 
			geom_rect() + 
			geom_vline(xintercept = start_min) +
			geom_vline(xintercept = end_min) +
			scale_fill_manual(values = cna_colors) + 
			scale_color_manual(values = cna_colors) + 
			coord_cartesian(xlim = c(start(gr_canvas), end(gr_canvas))) + 
			plot_theme() + 
			ggtitle(NULL, plot_title) + 
			theme(plot.subtitle = element_text(size = 8), axis.text.x = element_text(size = 5, angle = 90))
			
	gr_genes_here = gr_genes[gr_genes$symbol %in% genes_as_coordinates]
	gr_cds_here = gr_cds[gr_cds$symbol %in% genes_as_coordinates]
	plt_genes = generate_gene_plot(gr_canvas, gr_genes_here, gr_cds_here) + 
			theme(axis.text.x = element_text(size = 5, angle = 90))

	if (length(genes_as_counts) > 0) {
		genes_CNAenriched_to_plot = genes_CNAenriched_dfr[genes_CNAenriched_dfr$gene %in% genes_as_counts,]
		genes_CNAenriched_to_plot$gene = factor(genes_CNAenriched_to_plot$gene, levels =genes_CNAenriched_to_plot$gene)
		
		genes_CNAenriched_to_plot$gene_type = "other"
		genes_CNAenriched_to_plot$gene_type [genes_CNAenriched_to_plot$gene %in% cgc2024] = "cgc"
		genes_CNAenriched_to_plot$gene_type [genes_CNAenriched_to_plot$gene %in% prostate_cancer_genes] = "prca"
		
		plt_count_genes = ggplot(genes_CNAenriched_to_plot, aes(gene, n_patients, fill = gene_type)) +
				geom_bar(stat = "identity") +
				scale_fill_manual( values = c("other" = "grey", "cgc" = "orange", "prca" = "darkred")) +
				plot_theme() +
				coord_flip() +
				theme(axis.text.y = element_text(size = 5))
	} else {
		plt_count_genes = ggplot()	
	}
	
	(plt_segments / plt_genes / plt_count_genes) + plot_layout(heights = unit(c(4, 1, 2), 'null'))
}


load(pff("CNA_modules_annotated.rsav"))
all_plts = lapply(CNA_modules_annotated$id, plot_tornado, 
				CNA_modules_annotated, gr_CNA_segs, all_patients, 
				gr_gene_coords, gr_CDS, "ALL patients", cna_colors, gene_CNA_enriched, canvas_flank = 0, 
				add_missing_patients = FALSE)

fname = pff("figures/Tornado_plots__CNA_modules.pdf")
pdf(fname, height = 8)
lapply(all_plts, function(x) x)
dev.off()
file_open_call2(fname)


load(pff("GISTIC_peaks_annotated.rsav"))
all_plts = lapply(GISTIC_peaks_annotated$id, plot_tornado, 
				GISTIC_peaks_annotated, gr_CNA_segs, all_patients, 
				gr_gene_coords, gr_CDS, "ALL patients", cna_colors, gene_CNA_enriched, canvas_flank = 0, 
				add_missing_patients = FALSE)

fname = pff("figures/Tornado_plots__GISTIC_peaks.pdf")
pdf(fname, height = 8)
lapply(all_plts, function(x) x)
dev.off()
file_open_call2(fname)





# additional genes plotted separately
load(pff("gr_gene_coords.rsav"))
gr_gene_coords$id = gr_gene_coords$id2
load(pff("gene_CNA_enriched.rsav"))

CANVAS_FLANK = 1e7

# gene = "MYC"; annot = "all_gain"
# gene = "SLC30A4"; annot = "all_gain"
format_gene_for_tornado_plots = function(gene, annot, gr_gene_coords, gene_CNA_enriched) {

	gr_this_gene = gr_gene_coords[gr_gene_coords$symbol == gene]
	chr = as.character(seqnames(gr_this_gene))[[1]]
	start = start(gr_this_gene)
	end = end(gr_this_gene)
	id = paste(gene, annot, sep = "__")
	
	cnas_this_gene = gene_CNA_enriched[gene_CNA_enriched$gene == gene & gene_CNA_enriched$annot == annot, ]
	patients = cnas_this_gene$patients

	data.frame(id, chr, start, end, gene, all_genes = gene, cgc_genes = "", prca_genes = "", patients, annot, stringsAsFactors = FALSE)	
}


genes_to_plot = c("FOXA1", "MYC", "TERT", "FAM27C")
driver_gene_cnas = do.call(rbind, lapply(genes_to_plot, format_gene_for_tornado_plots, "all_gain", gr_gene_coords, gene_CNA_enriched))
all_plts = lapply(driver_gene_cnas$id, plot_tornado, 
				driver_gene_cnas, gr_CNA_segs, all_patients, 
				gr_gene_coords, gr_CDS, "ALL patients", cna_colors, gene_CNA_enriched, canvas_flank = CANVAS_FLANK, 
				add_missing_patients = FALSE)
				
fname = pff("figures/Tornado_plots__driver_genes.pdf")
pdf(fname, height = 8)
lapply(all_plts, function(x) x)
dev.off()
file_open_call2(fname)