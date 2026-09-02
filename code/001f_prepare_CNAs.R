source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
source(paste0("../bin/", date_tag, "/001xx_functions_for_CNA_processing.R"))

# which samples should we select? use only the samples from the tracking sheet: one per patient
load(file = pff("blacklisted_tumor_ids.rsav"))
load(file = pff("tracking_sheet.rsav"))
distinct_tumor_ids = tracking_sheet$WGS_AssayID



#
# collect whole genome duplication info
#
WGD_PLOIDY_THRESHOLD = 2.5

WGD_status = read.csv(
 		paste0("DATA_USED__", this_timestamp, "/SCNA_with_BD_SVs_03_05_2023/Cellularity_Ploidy_Estimates_03_May_2023.csv"),
		stringsAsFactors = FALSE)		
WGD_status$patient = gsub("(.+)._DNA", "\\1", WGD_status$PPCG_Sample_ID)
WGD_status$tumor_id = WGD_status$PPCG_Sample_ID
WGD_status$has_WGD = !is.na(WGD_status$Ploidy) & WGD_status$Ploidy >= WGD_PLOIDY_THRESHOLD

WGD_status = WGD_status[WGD_status$tumor_id %in% distinct_tumor_ids,]
WGD_status = WGD_status[!(WGD_status$tumor_id %in% blacklisted_tumor_ids),]
WGD_status = WGD_status[, c("patient", "tumor_id", "has_WGD")]
WGD_tumors = unique(WGD_status[WGD_status$has_WGD, "tumor_id"])
save(WGD_status, file = pff("WGD_status.rsav"))


#
# read in all CNA segments
#
fnames = list.files(
		paste0("DATA_USED__", this_timestamp, "/SCNA_with_BD_SVs_03_05_2023/Subclonal_SCNA_with_Avg_CN/"), 
		full.names = TRUE)
prepared_CNAs = do.call(rbind, mclapply(fnames, get_CNA_calls, WGD_tumors, CNA_cutoffs, MIN_SEGMENT_SIZE, mc.cores = 16))
prepared_CNAs = prepared_CNAs[!prepared_CNAs$tumor_id %in% blacklisted_tumor_ids,]
# keep one sample per tumor only
prepared_CNAs = prepared_CNAs[prepared_CNAs$tumor_id %in% distinct_tumor_ids,]
prepared_CNAs$patient = gsub("._DNA$", "", prepared_CNAs$tumor_id)
save(prepared_CNAs, file = pff("prepared_CNAs.rsav"))


##
# create TSV file for gistic input
##
load(file = pff("prepared_CNAs.rsav"))
prepared_CNAs_for_gistic = cbind(prepared_CNAs[, c("tumor_id", "chr", "startpos", "endpos")], Num.Markers = NA, 
		Seg.CN = prepared_CNAs[, "foldchange_cn"], 
		stringsAsFactors = FALSE)
prepared_CNAs_for_gistic$Seg.CN = log2(prepared_CNAs_for_gistic$Seg.CN)
prepared_CNAs_for_gistic$Seg.CN[prepared_CNAs_for_gistic$Seg.CN == -Inf] = -10
colnames(prepared_CNAs_for_gistic) = c("Sample", "Chromosome", "Start.bp", "End.bp", "Num.Markers", "Seg.CN")

gistic_folder = pff("GISTIC2_run1")
system(paste("mkdir -pv", gistic_folder))
fname = paste0(gistic_folder, "/CNA_segments_for_gistic.tsv")
write.table(prepared_CNAs_for_gistic, file = fname, quote = F, row.names = F, sep = "\t")
dfr_PGA_estimates = do.call(rbind, lapply(unique(prepared_CNAs$patient), compute_PGA, prepared_CNAs))
save(dfr_PGA_estimates, file = pff("dfr_PGA_estimates.rsav"))



##
# prepared CNA overlaps with genes
##
load(file = pff("gr_CDS.rsav"))

gr_prepared_CNAs = GRanges(prepared_CNAs$chr, IRanges(prepared_CNAs$startpos, prepared_CNAs$endpos))
gr_prepared_CNAs$annot = prepared_CNAs$annot
gr_prepared_CNAs$patient = prepared_CNAs$patient
gr_prepared_CNAs$total_cn = prepared_CNAs$total_cn
gr_prepared_CNAs$rel_cn = prepared_CNAs$rel_cn
gr_prepared_CNAs$CNA_id = prepared_CNAs$CNA_id

ov = findOverlaps(gr_prepared_CNAs, gr_CDS)
gc()

gene_cna_annots = data.frame(
		gene = gr_CDS[subjectHits(ov)]$symbol, 
		annot = gr_prepared_CNAs[queryHits(ov)]$annot,
		patient = gr_prepared_CNAs[queryHits(ov)]$patient,
		total_cn = gr_prepared_CNAs[queryHits(ov)]$total_cn,
		cn_width = width(gr_prepared_CNAs[queryHits(ov)]),
		rel_cn = gr_prepared_CNAs[queryHits(ov)]$rel_cn,
		chr = as.character(seqnames(gr_prepared_CNAs[queryHits(ov)])),
		CNA_id = as.character(gr_prepared_CNAs[queryHits(ov)]$CNA_id),
		stringsAsFactors = FALSE)
gc()

cna_tag = paste(
		gene_cna_annots[, "gene"],
		gene_cna_annots[, "annot"],
		gene_cna_annots[, "patient"],
		gene_cna_annots[, "total_cn"],
		gene_cna_annots[, "cn_width"],
		gene_cna_annots[, "rel_cn"],
		gene_cna_annots[, "chr"],
		sep = ",")
gc()

select_not_duplicated = which(!duplicated(cna_tag))
gene_cna_annots = gene_cna_annots[select_not_duplicated,]
save(gene_cna_annots, file = pff("gene_cna_annots.rsav"))
gc()



##
# percent overlaps of genes and CNAs per patient
##
load(file = pff("prepared_CNAs.rsav"))
load(file = pff("gr_prepared_elements.rsav"))
gr_genes = gr_prepared_elements[grep("^CDSgene", gr_prepared_elements$id)]

genes_with_high_gain = gene_cna_coverage("high_gain", prepared_CNAs, gr_genes, FRACTION_GENE_ALTERED, "high_gain")
gc()
genes_with_full_loss = gene_cna_coverage("full_loss", prepared_CNAs, gr_genes, FRACTION_GENE_ALTERED, "full_loss")
gc()
genes_with_gain = gene_cna_coverage(c("gain", "high_gain"), prepared_CNAs, gr_genes, FRACTION_GENE_ALTERED, "all_gain")
gc()
genes_with_loss = gene_cna_coverage(c("loss", "full_loss"), prepared_CNAs, gr_genes, FRACTION_GENE_ALTERED, "all_loss")
gc()

#
# remove ChrX from full loss genes
chrX_genes = gr_genes[seqnames(gr_genes) == "chrX"]$symbol
genes_with_full_loss = genes_with_full_loss[!genes_with_full_loss$gene %in% chrX_genes,]
genes_with_loss = genes_with_loss[!genes_with_loss$gene %in% chrX_genes,]

genes_with_full_CNAs = rbind(genes_with_high_gain, genes_with_full_loss, genes_with_gain, genes_with_loss)
save(genes_with_full_CNAs, file = pff("genes_with_full_CNAs.rsav"))
