source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(gplots)
library(patchwork)
library(ggsignif)
library(gtools)

load(pff("gr_prepared_elements.rsav"))
load(pff("results_signf_merged_annot.rsav"))
load(pff("variants_to_elements.rsav"))
load(pff("gencode_genes.rsav"))
gencode_genes$symbol = gsub("(.+)::(.+)::(.+)", "\\2", gencode_genes$id2)

print(load(file = pff("prepared_variants__SNV_indel.rsav")))
prep_var_SNV = prepared_variants
print(load(file = pff("prepared_variants_SV.rsav")))
prep_var_SV = prepared_variants_SV

# granges objects for quick intersecting
gr_prep_var_SV = GRanges(prep_var_SV$chr, IRanges(prep_var_SV$start, prep_var_SV$end))
gr_prep_var_SNV = GRanges(prep_var_SNV$chr, IRanges(prep_var_SNV$pos1, prep_var_SNV$pos2))
rm(prepared_variants_SV, prepared_variants)


print(load(pff("gr_CDS.rsav")))
print(load(pff("gr_gene_coords.rsav")))
# make the IDs for CDS and genes the same to make sure the gene_plot function workds
gr_CDS$id = gr_CDS$id2
gr_gene_coords$id = gr_gene_coords$id2


# this provides protein-coding (CDS) SNV driver list and CDS SNVs associated with each driver element
# dnds results are the basis of CDS mutations, but we need to focus on genes that were selected by DNDS
results_CDS_SNV_drivers = results_signf_merged_annot[results_signf_merged_annot$mut_type == "SNV_CDS", "annots_MAIN"]
load(pff("dnds_results_annotations_all.rsav"))
dnds_results_annotations_RGA = dnds_results_annotations_all[dnds_results_annotations_all$gene %in% results_CDS_SNV_drivers,] 
rm(dnds_results_annotations_all, results_CDS_SNV_drivers)

mut_colors = c(
		"SV_DEL"     = "royalblue1" ,
		"SV_TRA"     = "darkolivegreen" ,
		"SNV_noncoding" = "darkorange" ,	
		"SNV_CDS_subst"     = "gold" ,
		"SV_INV"     = "darkkhaki" ,
		"SNV_CDS_trunc"  = "azure4" ,
		"SV_DUP"     = "brown2"
)

plot_local_theme = function() {
	theme(
			axis.text.x = element_text(angle = 0, hjust = 0.5),
			axis.title = element_text(size = 12), 
			legend.key.size = unit(0.5, 'cm'), 
			legend.text = element_text(size = 6)
	)	
}

create_needleplot = function(driver, gr_gene_coords, gr_CDS, results_signf_merged_annot, gr_prepared_elements, 
		prep_var_SNV, prep_var_SV, gr_prep_var_SV, gr_prep_var_SNV, dnds_results_annotations_RGA, CANVAS_FLANK, BASEPAIRS_TO_NEEDLES) {
			
	driver_res = results_signf_merged_annot[results_signf_merged_annot$annots_MAIN == driver,]
	
	loc_chr = unique(driver_res$location_chr)
	if (length(loc_chr) != 1) {
		cat("nrow(chr) != 1", driver, "\n")
		stop()
	}
	loc_start = min(driver_res$location_start)
	loc_end = max(driver_res$location_end)
	
	canvas_min = loc_start - CANVAS_FLANK
	canvas_max = loc_end + CANVAS_FLANK
	gr_canvas = GRanges(loc_chr, IRanges(canvas_min, canvas_max))

	plt_genes = generate_gene_plot(gr_canvas, gr_gene_coords, gr_CDS)
		
	SNV_driver_elements = unique(unlist(strsplit(driver_res[driver_res$mut_type %in% c("SNV_CDS", "SNV_NC"), "element_ids"], split = ",")))
	SV_driver_elements = unique(unlist(strsplit(driver_res[driver_res$mut_type == "SV", "element_ids"], split = ",")))
	gr_SNV_driver_element_coords = gr_prepared_elements[gr_prepared_elements$id %in% SNV_driver_elements]
	gr_SV_driver_element_coords = gr_prepared_elements[gr_prepared_elements$id %in% SV_driver_elements]
	
	# driver elements combine the span of the plot
	gr_combined_driver_element_coords = unique(c(gr_SNV_driver_element_coords, gr_SV_driver_element_coords))
	combined_driver_element_coords = as.data.frame(gr_combined_driver_element_coords, stringsAsFactors = FALSE)

	# another plot with non-coding elements - ie exclude lncRNA, CDS, CDSgene
	nc_element_coords = c(gr_SNV_driver_element_coords, gr_SV_driver_element_coords)
	nc_element_coords = nc_element_coords[grep("^CDSgene::|^CDS::|^lncRNA::", nc_element_coords$id, invert = TRUE)]
	nc_element_coords$rownum = as.numeric(factor( nc_element_coords$id))
	nc_element_coords = as.data.frame(nc_element_coords, stringsAsFactors = FALSE)
	nc_element_coords$label = nc_element_coords$id
	nc_element_coords$label[duplicated(nc_element_coords$label)] = NA
	
	# mutations: SNVs and SVs separately
	# first get all mutations in the canvas, then color them based on location in element or flank
	SNVs_here = prep_var_SNV[queryHits(findOverlaps(gr_prep_var_SNV, gr_canvas)),]
	SVs_here = prep_var_SV[queryHits(findOverlaps(gr_prep_var_SV, gr_canvas)),]
	
	if (nrow(SNVs_here) > 0) {
		# label the subset of variants at FMREs	
		gr_SNVs_here = GRanges(SNVs_here$chr[[1]], IRanges(SNVs_here$pos1, SNVs_here$pos2))
		SNVs_here$region_class = "background"
		SNVs_here[queryHits(findOverlaps(gr_SNVs_here, gr_SNV_driver_element_coords)), "region_class"] = "driver"

		# round coordinates to 10 for a better overview
		SNVs_here$mid = round((SNVs_here$pos1 + SNVs_here$pos2) / 2, -log10(BASEPAIRS_TO_NEEDLES))
		SNV_counts = melt(dcast(SNVs_here, mid ~ region_class, fun.aggregate = length, value.var = "patient"),
						id.var = "mid")
		SNV_counts = SNV_counts[SNV_counts$value > 0,]
		
		# focus only on muts in the driver regions; plot their classes	
		SNVs_here_drivers = SNVs_here[SNVs_here$region_class == "driver",]
	} else {
		SNVs_here_drivers = SNV_counts = data.frame()
	}
	
	if (nrow(SVs_here) > 0) {
		# label the subset of variants at FMREs	
		gr_SVs_here = GRanges(SVs_here$chr[[1]], IRanges(SVs_here$start, SVs_here$end))
		SVs_here$region_class = "background"
		SVs_here[queryHits(findOverlaps(gr_SVs_here, gr_SV_driver_element_coords)), "region_class"] = "driver"
		
		# round coordinates to 10 for a better overview
		SVs_here$mid = round((SVs_here$start + SVs_here$end) / 2, -log10(BASEPAIRS_TO_NEEDLES))

		# convert SVs to per-position counts
		SV_counts = melt(dcast(SVs_here, mid ~ region_class, fun.aggregate = length, value.var = "patient"),
						id.var = "mid")
		SV_counts = SV_counts[SV_counts$value > 0,]

		# focus only on muts in the driver regions; plot their classes	
		SVs_here_drivers = SVs_here[SVs_here$region_class == "driver",]
	} else {
		SVs_here_drivers = SV_counts = data.frame()
	}
	
	# curate mutation classes for SNV/indel based on PCG impact
	if (nrow(SNVs_here_drivers) > 0) {
		
		trunc_impacts = c("Nonsense", "indel_delfrshift", "Essential_Splice", "indel_insfrshift", "Stop_loss")
		subst_impacts = c("Missense", "indel_delinframe")
		
		SNVs_here_drivers$mutation_type = "SNV"
			
		# merge current annotations (annovar) with dnds annotations, all.x TRUE to keep NAs for muts that lack annotations
		dnds_results_annotations_RGA$chr_chr = paste0("chr", dnds_results_annotations_RGA$chr)
		SNVs_here_drivers1 = merge(SNVs_here_drivers, dnds_results_annotations_RGA, 
			by.x = c("chr", "pos1", "ref", "alt", "patient"),
			by.y = c("chr_chr", "pos", "ref", "mut", "sampleID"), all.x = TRUE)
		
		# non-syn are excluded. in the CDS analysis dnds does not see these as drivers; non-coding analysis excludes CDS exons
		SNVs_here_drivers1 = SNVs_here_drivers1[!SNVs_here_drivers1$ExonicFunc.refGene == "synonymous SNV",]

		SNVs_here_drivers1$mutation_subtype = ""
		SNVs_here_drivers1[ !is.na(SNVs_here_drivers1$impact) & SNVs_here_drivers1$impact %in% trunc_impacts, "mutation_subtype" ] = "SNV_CDS_trunc"
		SNVs_here_drivers1[ !is.na(SNVs_here_drivers1$impact) & SNVs_here_drivers1$impact %in% subst_impacts, "mutation_subtype" ] = "SNV_CDS_subst"
		SNVs_here_drivers1[ is.na(SNVs_here_drivers1$impact) & SNVs_here_drivers1$mutation_subtype == "", "mutation_subtype" ] = "SNV_noncoding"
		
		SNVs_here_drivers = unique(SNVs_here_drivers1[, c("patient", "mutation_type", "mutation_subtype")])
	} else {
		SNVs_here_drivers = NULL
	}

	if (nrow(SVs_here_drivers) > 0)	{
		# SV mutation subtypes: just collapse inversions	
		SVs_here_drivers$mutation_type = "SV"
		SVs_here_drivers$mutation_subtype = ""
		SVs_here_drivers$mutation_subtype = paste0("SV_", SVs_here_drivers$top_mut_signt)
		SVs_here_drivers$mutation_subtype[ SVs_here_drivers$mutation_subtype %in% c("SV_t2tINV", "SV_h2hINV") ] = "SV_INV"	
		SVs_here_drivers = unique(SVs_here_drivers[, c("patient", "mutation_type", "mutation_subtype")])
	} else {
		SVs_here_drivers = NULL
	}
	
	SV_SNV_drivers = rbind(SVs_here_drivers, SNVs_here_drivers)
	SV_SNV_drivers$mutation_type = factor(SV_SNV_drivers$mutation_type, levels = c("SNV", "SV"))
			
	# count patients based on SNV,SNV, both
	patients_by_mut_type = dcast(SV_SNV_drivers, patient~1, value.var="mutation_type", fun.aggregate = function(x) paste(unique(sort(x)), collapse=","))
	colnames(patients_by_mut_type) = c("patient", "mut_type")
	
	# for plotting, use megabases rather than bp-coordinates
	MBP = 1e6
	
	fdr = signif(min(driver_res$fdr), 2)
	all_fdrs = c(by(driver_res$fdr, driver_res$mut_type, function(x) signif(min(x), 2)))
	fdr_cds = all_fdrs["SNV_CDS"]
	fdr_nc = all_fdrs["SNV_NC"]
	fdr_sv = all_fdrs["SV"]

	fdr_title = paste0("FDR_cds=", fdr_cds, "; FDR_nc=", fdr_nc, "; FDR_sv=", fdr_sv)
	plot_params = paste0("canvas", CANVAS_FLANK, "__NeedleBps", BASEPAIRS_TO_NEEDLES)
	plot_subtitle = paste0(driver, ";    ", loc_chr, "\n", 
			length(unique(patients_by_mut_type$patient)), 
			" patients; flnk=", CANVAS_FLANK, "\n",
			fdr_title, "\n", plot_params)
	cat(plot_subtitle, "\n")
							
	plt_nc_elements = ggplot(nc_element_coords, aes(xmin = start/MBP, xmax = end/MBP, ymin = 0 + rownum, ymax = 1 + rownum, label = label)) + 
			geom_rect() + 
			geom_text(aes(x = (canvas_min+canvas_max)/2/MBP, y = rownum + 0.5), size = 2) + 
			coord_cartesian(xlim = c(canvas_min, canvas_max)/MBP) +
			scale_y_continuous("NC_el", breaks = NULL, labels = NULL) +
			scale_x_continuous(NULL, labels = NULL) +
			plot_theme() + plot_local_theme()

	# if there are no mutations of this type, the object is a null and we need to have an empty plot
	plt_SNV_counts = ggplot()
	if (nrow(SNV_counts) > 0) {
		plt_SNV_counts = ggplot(SNV_counts, aes(x = mid/MBP, y = value, color = variable)) + 
				geom_bar(stat = "identity") +
				geom_point() +
				scale_x_continuous(NULL, labels = NULL) +
				scale_y_continuous(paste0("n. SNVs \n/ ", BASEPAIRS_TO_NEEDLES, " bps")) +
				scale_color_manual(NULL, values = c("driver"= "black", "background" = "grey")) + 
				coord_cartesian(xlim = c(canvas_min, canvas_max)/MBP) +
				plot_theme() + plot_local_theme() + 
				ggtitle(NULL, plot_subtitle)
	}

	plt_SV_counts = ggplot()
	if (nrow(SV_counts) > 0) {
		plt_SV_counts = ggplot(SV_counts, aes(x = mid/MBP, y = value, color = variable)) + 
				geom_bar(stat = "identity") +
				geom_point() +
				scale_y_continuous(paste0("n. SVBPs \n/ ", BASEPAIRS_TO_NEEDLES, " bps")) +
				scale_x_continuous(NULL) +
				scale_color_manual(NULL, values = c("driver"= "black", "background" = "grey")) + 
				coord_cartesian(xlim = c(canvas_min, canvas_max)/MBP) +
				plot_theme() + plot_local_theme()
	}
			
	plt_driver_mut_subtypes = ggplot(SV_SNV_drivers, aes(mutation_type, fill = mutation_subtype)) + 
			geom_bar() + 
			coord_flip() +
			scale_fill_manual(NULL, values = mut_colors) +
			scale_x_discrete(NULL) +
			plot_theme() + plot_local_theme() + theme(legend.position = "bottom")
			
	plt_patient_mut_types = ggplot(patients_by_mut_type, aes(factor(1), fill = mut_type)) + 
			geom_bar() + 
			coord_flip() +
			scale_fill_manual(NULL, values = c("SNV" = "salmon" ,"SV" = "forestgreen", "SNV,SV" = "cyan")) +
			scale_x_discrete(NULL, labels = "patients") +
			scale_y_continuous(NULL) +
			plot_theme() + plot_local_theme() + theme(legend.position = "bottom")
	
	# remove coords from gene plot
	plt_genes = plt_genes + scale_x_continuous(NULL, labels = NULL)
		
	plt_combined = plt_SNV_counts / plt_SV_counts / plt_genes / plt_nc_elements / plt_patient_mut_types / plt_driver_mut_subtypes + 
			plot_layout(heights = unit(c(4, 2, 1, 1, 1, 1), 'null')) 

	plt_combined
}


CANVAS_FLANK = 10000
BASEPAIRS_TO_NEEDLES = 10


# top p-vals first
genes = names(sort(by(results_signf_merged_annot$fdr, results_signf_merged_annot$annots_MAIN, min)))
genes = setdiff(genes, "ETS")

fname = pff(c("figures/needle_plots_new__canvas", CANVAS_FLANK, "__NeedleBps", BASEPAIRS_TO_NEEDLES, ".pdf"))
pdf(fname, height = 8)

lapply(genes, create_needleplot, gr_gene_coords, gr_CDS, results_signf_merged_annot, gr_prepared_elements, 
		prep_var_SNV, prep_var_SV, gr_prep_var_SV, gr_prep_var_SNV, dnds_results_annotations_RGA,
		CANVAS_FLANK = CANVAS_FLANK, BASEPAIRS_TO_NEEDLES = BASEPAIRS_TO_NEEDLES)

dev.off()

file_open_call2(fname)



CANVAS_FLANK = 10000
BASEPAIRS_TO_NEEDLES = 1


# top p-vals first
genes = names(sort(by(results_signf_merged_annot$fdr, results_signf_merged_annot$annots_MAIN, min)))
genes = setdiff(genes, "ETS")

fname = pff(c("figures/needle_plots_new__canvas", CANVAS_FLANK, "__NeedleBps", BASEPAIRS_TO_NEEDLES, ".pdf"))
pdf(fname, height = 8)

lapply(genes, create_needleplot, gr_gene_coords, gr_CDS, results_signf_merged_annot, gr_prepared_elements, 
		prep_var_SNV, prep_var_SV, gr_prep_var_SV, gr_prep_var_SNV, dnds_results_annotations_RGA,
		CANVAS_FLANK = CANVAS_FLANK, BASEPAIRS_TO_NEEDLES = BASEPAIRS_TO_NEEDLES)

dev.off()

file_open_call2(fname)


CANVAS_FLANK = 50000
BASEPAIRS_TO_NEEDLES = 10


# top p-vals first
genes = names(sort(by(results_signf_merged_annot$fdr, results_signf_merged_annot$annots_MAIN, min)))
genes = setdiff(genes, "ETS")

fname = pff(c("figures/needle_plots_new__canvas", CANVAS_FLANK, "__NeedleBps", BASEPAIRS_TO_NEEDLES, ".pdf"))
pdf(fname, height = 8)

lapply(genes, create_needleplot, gr_gene_coords, gr_CDS, results_signf_merged_annot, gr_prepared_elements, 
		prep_var_SNV, prep_var_SV, gr_prep_var_SV, gr_prep_var_SNV, dnds_results_annotations_RGA, 
		CANVAS_FLANK = CANVAS_FLANK, BASEPAIRS_TO_NEEDLES = BASEPAIRS_TO_NEEDLES)

dev.off()

file_open_call2(fname)


CANVAS_FLANK = 100000
BASEPAIRS_TO_NEEDLES = 100

fname = pff(c("figures/needle_plots_new__canvas", CANVAS_FLANK, "__NeedleBps", BASEPAIRS_TO_NEEDLES, ".pdf"))
pdf(fname, height = 8)

lapply(genes, create_needleplot, gr_gene_coords, gr_CDS, results_signf_merged_annot, gr_prepared_elements, 
		prep_var_SNV, prep_var_SV, gr_prep_var_SV, gr_prep_var_SNV, dnds_results_annotations_RGA, 
		CANVAS_FLANK = CANVAS_FLANK, BASEPAIRS_TO_NEEDLES = BASEPAIRS_TO_NEEDLES)

dev.off()

file_open_call2(fname)


CANVAS_FLANK = 1000000
BASEPAIRS_TO_NEEDLES = 1000

fname = pff(c("figures/needle_plots_new__canvas", CANVAS_FLANK, "__NeedleBps", BASEPAIRS_TO_NEEDLES, ".pdf"))
pdf(fname, height = 8)

lapply(genes, create_needleplot, gr_gene_coords, gr_CDS, results_signf_merged_annot, gr_prepared_elements, 
		prep_var_SNV, prep_var_SV, gr_prep_var_SV, gr_prep_var_SNV, dnds_results_annotations_RGA, 
		CANVAS_FLANK = CANVAS_FLANK, BASEPAIRS_TO_NEEDLES = BASEPAIRS_TO_NEEDLES)

dev.off()

file_open_call2(fname)



CANVAS_FLANK = 50000
BASEPAIRS_TO_NEEDLES = 1


# top p-vals first
genes = names(sort(by(results_signf_merged_annot$fdr, results_signf_merged_annot$annots_MAIN, min)))
genes = setdiff(genes, "ETS")

fname = pff(c("figures/needle_plots_new__canvas", CANVAS_FLANK, "__NeedleBps", BASEPAIRS_TO_NEEDLES, ".pdf"))
pdf(fname, height = 8)

lapply(genes, create_needleplot, gr_gene_coords, gr_CDS, results_signf_merged_annot, gr_prepared_elements, 
		prep_var_SNV, prep_var_SV, gr_prep_var_SV, gr_prep_var_SNV, dnds_results_annotations_RGA,
		CANVAS_FLANK = CANVAS_FLANK, BASEPAIRS_TO_NEEDLES = BASEPAIRS_TO_NEEDLES)

dev.off()
file_open_call2(fname)

