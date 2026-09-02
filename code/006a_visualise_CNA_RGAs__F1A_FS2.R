source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(ggrepel)
library(patchwork)
library(gtools)

FDR_CUT = 0.05

load(pff("cgc2024.rsav"))
load(pff("prostate_cancer_genes.rsav"))
cancer_genes = unique(c(cgc2024, prostate_cancer_genes))

load(pff("all_patients.rsav"))
MIN_DRIVER_FREQ = 0.01
MIN_MUT_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)

# plot all gains/amps module genes by absolute CN

get_gene_amps = function(gene, focal_amps, genes_with_high_amps) {
	this_amp = focal_amps[focal_amps$gene == gene,]
	amp_values = genes_with_high_amps[genes_with_high_amps$gene == gene, "rel_cn"]
	dfr = data.frame(gene, amp_values, fdr = this_amp$fdr, n_patients = this_amp$n_patients, stringsAsFactors = FALSE)
	dfr$fdr[-1] = NA
	dfr$n_patients[-1] = NA
	dfr
}

load(pff("gr_gene_coords.rsav"))
load(pff("gene_CNA_enriched.rsav"))
load(pff("gene_cna_annots.rsav"))
load(pff("gene_CNA_module_dfr.rsav"))

gain_modules = gene_CNA_module_dfr[gene_CNA_module_dfr$annot %in% c("all_gain", "high_gain"),]
gain_modules = gain_modules[gain_modules$n_patients >= MIN_MUT_PATIENTS, ]
genes_from_gain_modules = unique(gain_modules$gene)

cancer_gene_colors = c("other" = "cornflowerblue", "PrCa" = "darkred", "CGC" = "darkorange")

# gene FDR-values are those from per-gene estimates; take best of high_gain vs gain
gene_based_FDRs = gene_CNA_enriched[gene_CNA_enriched$annot %in% c("high_gain", "all_gain") & gene_CNA_enriched$gene %in% genes_from_gain_modules,]
gene_based_FDRs = by(gene_based_FDRs$fdr, gene_based_FDRs$gene, min)

# add cytoband labels to genes
cytobands = read.delim("elements/cytoBand.txt", stringsAsFactors = FALSE, header = FALSE)
cytobands$id = paste0(cytobands$V1, cytobands$V4)


# sample counts per gene
gene_based_sample_counts = gene_CNA_enriched[
		gene_CNA_enriched$annot %in% c("high_gain", "all_gain") & gene_CNA_enriched$gene %in% genes_from_gain_modules,]
gene_based_sample_counts = by(gene_based_sample_counts$patients, 
		gene_based_sample_counts$gene, 
		function(x) length(unique(unlist(strsplit(x, s = ",")))))


genes_CN_values = gene_cna_annots[
		gene_cna_annots$annot %in% c("gain", "high_gain") & 
		gene_cna_annots$gene %in% genes_from_gain_modules,]
genes_CN_values = melt(dcast(genes_CN_values, gene~patient, value.var = "total_cn", fun.aggregate = max))
genes_CN_values = genes_CN_values[!genes_CN_values$value %in% c(Inf, -Inf),]
colnames(genes_CN_values) = c("gene", "patient", "total_CN")

# get cytobands 
gr_genes_here = gr_gene_coords[gr_gene_coords$symbol %in% genes_from_gain_modules]
gene2cytoband = get_gene_cytobands(gr_genes_here, cytobands)
gene2cytoband = structure(names = gene2cytoband$gene, gene2cytoband$cytoband)

# assign cytobands
genes_CN_values$cytoband = NA
genes_CN_values$cytoband = gene2cytoband[genes_CN_values$gene]
genes_CN_values$gene_cytoband = paste0(genes_CN_values$gene, " (", genes_CN_values$cytoband, ")")

# get FDR values for genes
genes_CN_values$fdr = NA
genes_CN_values$fdr = gene_based_FDRs[genes_CN_values$gene]
genes_CN_values[duplicated(genes_CN_values$gene), "fdr"] = NA

# get sample counts for genes
genes_CN_values$n_patients = NA
genes_CN_values$n_patients = gene_based_sample_counts[genes_CN_values$gene]
genes_CN_values[duplicated(genes_CN_values$gene), "n_patients"] = NA

# most to least FDR
gene_order =  unique(genes_CN_values$gene_cytoband[order(genes_CN_values$fdr)])
genes_CN_values$gene_cytoband = factor(genes_CN_values$gene_cytoband, levels = gene_order)

# label genes based on evidence
genes_CN_values$gene_type = "other"
genes_CN_values$gene_type [genes_CN_values$gene %in% cgc2024] = "CGC"
genes_CN_values$gene_type [genes_CN_values$gene %in% prostate_cancer_genes] = "PrCa"


n_genes = length(unique(genes_CN_values$gene))
plt_title = paste0("gain/amp genes (",  n_genes, "): FDR<", FDR_CUT)


plt = ggplot(genes_CN_values, aes(gene_cytoband, total_CN, color = gene_type)) +
		scale_color_manual(values = cancer_gene_colors) +
		geom_jitter(alpha = 0.5, width = 0.2, size = 2) + 
		scale_y_continuous("Total CN", trans = "log2", limits = c(2^1, 2^7), breaks = 2^(1:7), labels = 2^(1:7)) +
		geom_text(aes(label = stars.pval(fdr)), y = 7, size = 4, color = "black") + 
		geom_text(aes(label = n_patients), y = 1, size = 2, color = "black") + 
		plot_theme() + 
		ggtitle(NULL, plt_title)
		
fname = pff("figures/focal_amp_values_per_genes.pdf")
ggsave(plt, file = fname, height = 5, width = 8.5)
file_open_call2(fname)


# index genes of modules only, show a bar plot for sample counts
load(pff("gene_CNA_enriched.rsav"))
load(pff("gene_CNA_module_dfr.rsav"))

CNA_modules = gene_CNA_module_dfr

# get cytobands 
gr_genes_here = gr_gene_coords[gr_gene_coords$symbol %in% CNA_modules$gene]
gene2cytoband = get_gene_cytobands(gr_genes_here, cytobands)
gene2cytoband = structure(names = gene2cytoband[,1], gene2cytoband[,2])

CNA_modules$gene_with_cb = paste0(CNA_modules$gene, " (", gene2cytoband[CNA_modules$gene], ")")

gene_order = rev(names(sort(by(CNA_modules$n_patients, CNA_modules$gene_with_cb, sum))))
CNA_modules$gene_with_cb = factor(CNA_modules$gene_with_cb, levels = gene_order)
CNA_modules$gene_type = "other"
CNA_modules$gene_type [CNA_modules$gene %in% cgc2024] = "CGC"
CNA_modules$gene_type [CNA_modules$gene %in% prostate_cancer_genes] = "PrCa"
CNA_modules$annot = factor(CNA_modules$annot, levels = c("all_loss", "all_gain", "high_gain"))

cancer_gene_symbol = c("PrCa" = 1, "CGC" = 4)

plt_CNA_counts = ggplot(CNA_modules, aes(gene_with_cb, n_patients, shape = gene_type, fill = annot, label = stars.pval(fdr))) + 
		geom_bar(stat = "identity") + 
		geom_text(angle = 15) + 
		scale_fill_manual(values = c("all_loss" = "darkcyan", "all_gain" = "gold", "high_gain" = "salmon")) + 
		geom_text(aes(label = n_patients), y = 20, angle = 90, size = 3) + 
		geom_point(y = -10) + 
		scale_shape_manual(values = cancer_gene_symbol) + 
		facet_grid(cols = vars(annot), scale = "free", space = "free") + 
		plot_theme()

fname = pff("figures/CNA_modules_with_main_genes_freq.pdf")
ggsave(plt_CNA_counts, file = fname, width = 8.5)	
file_open_call2(fname)


#
# cna density plot across the chromosomes
#
load(pff("gene_CNA_module_dfr.rsav"))
load(pff("all_patients.rsav"))

load(file = pff("gr_prepared_elements.rsav"))
gr_genes = gr_prepared_elements[grep("^CDSgene", gr_prepared_elements$id)]


cna_colors = c("bal" = "white", 
	"gain" = "orange", "high_gain" = "darkred", "all_gain" = "orange",
	"loss" = "steelblue", "full_loss" = "darkblue", 	"all_loss" = "steelblue"
)


TILE_WIDTH = 1e6

gr_tiles = tileGenome(seqinfo(Hsapiens), tilewidth = TILE_WIDTH, cut.last.tile.in.chrom = TRUE)
gr_tiles = gr_tiles[seqnames(gr_tiles) %in% paste0("chr", c(1:22, "X"))]
dfr_tiles = as.data.frame(gr_tiles)
dfr_tiles$TILES_id = paste(dfr_tiles$seqnames, ":", dfr_tiles$start, "-", dfr_tiles$end, sep = "")

count_patients_by_tiles = function(gr_segs, gr_tiles, annot, gr_genes, gene_CNA_module_dfr, gene_annot, all_patients) {
	
	ov = findOverlaps(gr_segs, gr_tiles)
	dfr_segs_here = as.data.frame(gr_segs[queryHits(ov)])
	colnames(dfr_segs_here) = paste0("SEG_", colnames(dfr_segs_here))
	dfr_tiles_here = as.data.frame(gr_tiles[subjectHits(ov)])
	colnames(dfr_tiles_here) = paste0("TILES_", colnames(dfr_tiles_here))
	dfr_tiles_here$TILES_id = paste(dfr_tiles_here$TILES_seqnames, ":", dfr_tiles_here$TILES_start, "-", dfr_tiles_here$TILES_end, sep = "")
	
	dfr_here = cbind(dfr_segs_here, dfr_tiles_here)
	n_patients_by_tile = c(by(dfr_here$SEG_patient, dfr_here$TILES_id, function(x) length(unique(x))))	
	dfr_patients_by_tile = data.frame(TILES_id = names(n_patients_by_tile), n_patients = n_patients_by_tile, stringsAsFactors = FALSE)
	dfr_tiles_w_patient_counts = merge(dfr_tiles, dfr_patients_by_tile, by = "TILES_id", all = TRUE)
	
	dfr_tiles_w_patient_counts[is.na(dfr_tiles_w_patient_counts$n_patients), "n_patients"] = 0
	dfr_tiles_w_patient_counts$pct_cohort = 100 * dfr_tiles_w_patient_counts$n_patients / length(all_patients)
	dfr_tiles_w_patient_counts$annot = annot

	# add genes here: create tiles gr again and take the first window overlapping with a gene
	dfr_tiles_w_patient_counts$genes_involved = NA
	genes_here = gene_CNA_module_dfr[gene_CNA_module_dfr$annot %in% gene_annot, "gene"]
	gr_genes_here = gr_genes[gr_genes$symbol %in% genes_here]
	gr_tiles_for_genes = GRanges(dfr_tiles_w_patient_counts$seqnames, IRanges(dfr_tiles_w_patient_counts$start, dfr_tiles_w_patient_counts$end),
			id = dfr_tiles_w_patient_counts$TILES_id)
	ov2 = findOverlaps(gr_tiles_for_genes, gr_genes_here)

	if (length(ov2) > 0) {
		tile_gene_map = data.frame(tile_id = gr_tiles_for_genes[queryHits(ov2)]$id, gene = gr_genes_here[subjectHits(ov2)]$symbol)
		tile_gene_map = tile_gene_map[!duplicated(tile_gene_map$gene),]
		tile_gene_map = sapply(split(tile_gene_map$gene, tile_gene_map$tile_id), paste, collapse = ",")
		dfr_tiles_w_patient_counts$genes_involved = tile_gene_map[dfr_tiles_w_patient_counts$TILES_id]
	}
	
	dfr_tiles_w_patient_counts
}



load(file = pff("prepared_CNAs.rsav"))
gr_CNA_segs = GRanges(prepared_CNAs$chr, IRanges(prepared_CNAs$start, prepared_CNAs$end), 
		annot = prepared_CNAs$annot, patient = prepared_CNAs$patient)

gr_gain_segs = gr_CNA_segs[gr_CNA_segs$annot == "gain"]
dfr_tiles_gain = count_patients_by_tiles(gr_gain_segs, gr_tiles, "gain", gr_genes, gene_CNA_module_dfr, "all_gain", all_patients)

gr_loss_segs = gr_CNA_segs[gr_CNA_segs$annot == "loss"]
dfr_tiles_loss = count_patients_by_tiles(gr_loss_segs, gr_tiles, "loss", gr_genes, gene_CNA_module_dfr, "all_loss", all_patients)

gr_hAmp_segs = gr_CNA_segs[gr_CNA_segs$annot == "high_gain"]
dfr_tiles_hAmp = count_patients_by_tiles(gr_hAmp_segs, gr_tiles, "high_gain", gr_genes, gene_CNA_module_dfr, "high_gain", all_patients)

gr_hDel_segs = gr_CNA_segs[gr_CNA_segs$annot == "full_loss"]
dfr_tiles_hDel = count_patients_by_tiles(gr_hDel_segs, gr_tiles, "full_loss", gr_genes, gene_CNA_module_dfr, "full_loss", all_patients)

dfr_tiles_loss$n_patients = -1 * dfr_tiles_loss$n_patients
dfr_tiles_hDel$n_patients = -1 * dfr_tiles_hDel$n_patients

dfr_tiles_loss$pct_cohort = -1 * dfr_tiles_loss$pct_cohort
dfr_tiles_hDel$pct_cohort = -1 * dfr_tiles_hDel$pct_cohort

# one plot for gains, another for losses
dfr_tiles_gain_hAmp = rbind(dfr_tiles_gain, dfr_tiles_hAmp)
dfr_tiles_loss_hDel = rbind(dfr_tiles_loss, dfr_tiles_hDel)

dfr_tiles_loss_hDel$seqnames = gsub("chr", "", dfr_tiles_loss_hDel$seqnames)
dfr_tiles_gain_hAmp$seqnames = gsub("chr", "", dfr_tiles_loss_hDel$seqnames)

seq_levels = c(1:22, "X")
dfr_tiles_gain_hAmp$seqnames = factor(dfr_tiles_gain_hAmp$seqnames, levels = seq_levels)
dfr_tiles_loss_hDel$seqnames = factor(dfr_tiles_loss_hDel$seqnames, levels = seq_levels)


plt_amp = ggplot(dfr_tiles_gain_hAmp, aes(x= start, y = pct_cohort, color = annot, label = genes_involved)) + 
		geom_line() + 
		scale_x_continuous(NULL, breaks = NULL, labels = NULL) +
		facet_grid(~seqnames, scale = "free_x", space = "free_x") + 
		plot_theme() + 
		geom_point(data = dfr_tiles_gain_hAmp[!is.na(dfr_tiles_gain_hAmp$genes_involved),], 
				aes(x = start, y = pct_cohort), shape = 4, color = "black") +
		geom_text_repel(size = 2, angle = 45, max.overlaps = nrow(dfr_tiles_gain_hAmp), min.segment.length = 0) +
		scale_color_manual(values = cna_colors) +
		scale_y_continuous("% cohort", breaks = c(0, 5, 10, 15), labels = c(0, 5, 10, 15), lim = c(0, 18)) + 
		theme(legend.position = "bottom", strip.text = element_text(angle = 90, size = 5))

plt_del = ggplot(dfr_tiles_loss_hDel, aes(x= start, y = pct_cohort, color = annot, label = genes_involved)) + 
		geom_line() + 
		scale_x_continuous(NULL, breaks = NULL, labels = NULL) +
		facet_grid(~seqnames, scale = "free_x", space = "free_x") + 
		plot_theme() + 
		geom_point(data = dfr_tiles_loss_hDel[!is.na(dfr_tiles_loss_hDel$genes_involved),], 
				aes(x = start, y = pct_cohort), shape = 4, color = "black") +
		geom_text_repel(size = 2, angle = 45, max.overlaps = nrow(dfr_tiles_loss_hDel), min.segment.length = 0) +
		scale_color_manual(values = cna_colors) +
		scale_y_continuous("% cohort", breaks = c(0, -30, -60), labels = c(0, 30, 60), lim = c(-65, 0)) + 
		theme(legend.position = "bottom", strip.text = element_text(angle = 90, size = 5))

plt_combined = plt_amp / plt_del	
		
fname = pff("figures/CNA_genomewide_frequency.pdf")
ggsave(plt_combined, file = fname, width = 9, height = 5)
file_open_call2(fname)