source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(trackViewer)
library(RColorBrewer)
library(patchwork)


load(pff("cgc2024.rsav"))
load(pff("prostate_cancer_genes.rsav"))
cancer_genes = unique(c(cgc2024, prostate_cancer_genes))

## how large is the canvas around the extreme elements
canvas_site_flank = 10000


create_loops_plot = function(element_id, el2loop, loop2gene, loops, results_signf_merged_annot, 
		gr_genes, gr_CDS, gr_prepared_elements, cancer_genes) {

	cat(element_id, "\n")
	
	this_res = results_signf_merged_annot[results_signf_merged_annot$annots_MAIN == element_id,]
	this_FMRE_IDs = unique(unlist(strsplit(this_res$element_ids, split = ",")))
	this_FMRE_IDs = this_FMRE_IDs[grep("^CDS|^CDSgene|^UTR", this_FMRE_IDs, invert = TRUE)]
	
	if (length(this_FMRE_IDs) == 0 | all(is.na(this_FMRE_IDs)) ) {
		cat("no relevant elements at ", element_id, "\n")
		return(NULL)		
	}
	
	gr_this_FMREs = gr_prepared_elements[gr_prepared_elements$id %in% this_FMRE_IDs]
	
	all_loop_IDs = unique(unlist(el2loop[this_FMRE_IDs]))
	
	dfr_loops_here = loops[loops$tag %in% all_loop_IDs,]
	if (nrow(dfr_loops_here) == 0) {
		cat("no loops at ", element_id, "\n")
		return(NULL)
	}
	
	# take more FMREs from the other ends of the loops	
	other_FMRE_IDs = unique(unlist(strsplit(results_signf_merged_annot$element_ids, split = ",")))
	other_FMRE_IDs = as.matrix(stack(el2loop[other_FMRE_IDs]))
	other_FMRE_IDs = other_FMRE_IDs[other_FMRE_IDs[,1] %in% all_loop_IDs,"ind"]
	gr_other_FMREs = gr_prepared_elements[gr_prepared_elements$id %in% other_FMRE_IDs]
		
	loop_coords = unlist(dfr_loops_here[, c("start1", "end1", "start2", "end2")])
	coord_min = min(loop_coords) - canvas_site_flank
	coord_max = max(loop_coords) + canvas_site_flank
	gr_loop_coord_range = GRanges(dfr_loops_here$chr1[1], IRanges(coord_min, coord_max))
	
	# all genes; not just the genes touched by loop anchors
	#gr_this_genes = gr_genes[subjectHits(findOverlaps(gr_loop_coord_range, gr_genes))]
	# genes associated with loop anchors
	other_gene_symbols_at_loops = unique(gsub("(.+)::(.+)::(.+)::(.+)", "\\3", unique(unlist(loop2gene[all_loop_IDs]))))
	genes_to_show = unique(c(other_gene_symbols_at_loops, element_id))
	gr_this_genes = gr_genes[gr_genes$symbol %in% genes_to_show]
	gr_this_CDS = gr_CDS[gr_CDS$symbol %in% gr_this_genes$symbol]
	
	gr_canvas = gr_loop_coord_range
	dfr_canvas = as.data.frame(gr_canvas)
		
	# plot FMREs 
	dfr_plot_FMREs = unique(rbind(
			as.data.frame(gr_this_FMREs, stringsAsFactors = FALSE),
			as.data.frame(gr_other_FMREs, stringsAsFactors = FALSE)))
	
	dfr_plot_FMREs$site_type = gsub("(.+)::(.+)::(.+)::(.+)", "\\4", dfr_plot_FMREs$id)
	dfr_plot_FMREs$num = as.numeric(factor(dfr_plot_FMREs$id))
	ggplot_sites = ggplot(dfr_plot_FMREs, aes(xmin = start, xmax = end, ymin = num + 0, ymax = num + 0.95, label = site_type)) + 
			coord_cartesian(xlim = c(dfr_canvas$start[1], dfr_canvas$end[1])) + 
			scale_x_continuous(NULL, labels = NULL) + 
			scale_y_continuous("FMREs", breaks = NULL, labels = NULL) + 
			geom_rect(linewidth = 1) + plot_theme() +
			geom_text(aes(x = start, y = num + 0.5), size = 2) + 
			theme(legend.position = "bottom")

	ggplot_genes = generate_gene_plot(gr_canvas, gr_this_genes, gr_this_CDS) + theme(axis.text.x = element_text(angle = 0, size = 8, hjust = 0))

	region_coords = paste0(
			seqnames(gr_canvas), ":", 
			round(start(gr_canvas)/1e6, 1), "-", 
			round(end(gr_canvas)/1e6, 1), " Mbps")
	cancer_genes_here = intersect(gr_this_genes$symbol, cancer_genes)
	title_here = paste(element_id, " | ", region_coords, " | ",
			paste(cancer_genes_here, collapse = ", "))
	subtitle_here = paste(unique(gr_this_genes$symbol), collapse = ", ")
	
	# plot arcs for loops
	ggplot_loops = ggplot(dfr_loops_here) + 
			coord_cartesian(xlim = c(dfr_canvas$start[1], dfr_canvas$end[1]), ylim = c(0, 0.1)) + 
			geom_curve(aes(x = mid1, xend = mid2, y = 0, yend = 0), curvature = -1, alpha = 0.5, linewidth = 1) + 
			geom_segment(aes(x = start1, xend = end1, y = 0, yend = 0), alpha = 0.5, linewidth = 1) + 
			geom_segment(aes(x = start2, xend = end2, y = 0, yend = 0), alpha = 0.5, linewidth = 1) + 
			scale_x_continuous(NULL, labels = NULL) + 
			scale_y_continuous("loops", breaks = NULL, labels = NULL) + 
			plot_theme() + 
			ggtitle(title_here, subtitle_here) +
			theme(plot.subtitle = element_text(size = 8), plot.title = element_text(size = 11)) + 
			theme(legend.position = "bottom",  legend.text = element_text(size = 8))

	ggplot_combined = ggplot_loops / ggplot_sites / (ggplot_genes + scale_y_continuous("genes", breaks = NULL, labels = NULL))
	ggplot_combined = ggplot_combined + plot_layout(heights = unit(c(4, 1, 1), 'null'))	
	ggplot_combined
}


load(pff("gr_prepared_elements.rsav"))
load(pff("gr_CDS.rsav"))
load(pff("gr_gene_coords.rsav"))
# make the IDs for CDS and genes the same to make sure the gene_plot function works
gr_CDS$id = gr_CDS$id2
gr_gene_coords$id = gr_gene_coords$id2


load(file = pff("loops.rsav"))
load(file = pff("results_signf_merged_annot.rsav"))
load(file = pff("dfr_loops_element_to_gene.rsav"))
el2loop = split(dfr_loops_element_to_gene$loop_id, dfr_loops_element_to_gene$nc_id)
el2loop = lapply(el2loop, unique)
loop2gene = split(dfr_loops_element_to_gene$gene_id, dfr_loops_element_to_gene$loop_id)
loop2gene = lapply(loop2gene, unique)


all_el =  unique(results_signf_merged_annot$annots_MAIN)
cat("all_el ", length(all_el), "\n")

fname = pff(c("figures/FMRE_loops.pdf"))
pdf(fname)

lapply(all_el, create_loops_plot, 
		el2loop, loop2gene, loops, results_signf_merged_annot, gr_gene_coords, gr_CDS, gr_prepared_elements, cancer_genes)

dev.off()
file_open_call2(fname)
