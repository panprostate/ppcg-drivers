source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
source(paste0("../bin/", date_tag, "/001xx_functions_for_CNA_processing.R"))


load(file = pff("blacklisted_tumor_ids.rsav"))
# which samples should we select? use only the samples from the tracking sheet: one per patient
load(file = pff("allmets_tumor_ids.rsav"))


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
WGD_status_mets = WGD_status[WGD_status$tumor_id %in% allmets_tumor_ids,]
WGD_status_mets = WGD_status_mets[!(WGD_status_mets$tumor_id %in% blacklisted_tumor_ids),]
WGD_status_mets = WGD_status_mets[, c("patient", "tumor_id", "has_WGD")]
WGD_tumors_mets = unique(WGD_status_mets[WGD_status_mets$has_WGD, "tumor_id"])
save(WGD_status_mets, file = pff("WGD_status_mets.rsav"))


#
# read in all CNA segments
#
fnames = list.files(
		paste0("DATA_USED__", this_timestamp, "/SCNA_with_BD_SVs_03_05_2023/Subclonal_SCNA_with_Avg_CN/"), 
		full.names = TRUE)
tumor_id_here = gsub(
		paste0("DATA_USED__", this_timestamp, "/SCNA_with_BD_SVs_03_05_2023/Subclonal_SCNA_with_Avg_CN//(.+)_vs_(.+)_Subclonal_SCNA.txt"), 
		"\\1", fnames)
fnames = fnames[tumor_id_here %in% allmets_tumor_ids]
	
prepared_CNAs_mets = do.call(rbind, mclapply(fnames, get_CNA_calls, WGD_tumors_mets, CNA_cutoffs, MIN_SEGMENT_SIZE, mc.cores = 8))
prepared_CNAs_mets = prepared_CNAs_mets[!prepared_CNAs_mets$tumor_id %in% blacklisted_tumor_ids,]
prepared_CNAs_mets$patient = gsub("._DNA$", "", prepared_CNAs_mets$tumor_id)
save(prepared_CNAs_mets, file = pff("prepared_CNAs_mets.rsav"))



## prepared CNA overlaps with genes
load(file = pff("gr_CDS.rsav"))

gr_prepared_CNAs_mets = GRanges(prepared_CNAs_mets$chr, IRanges(prepared_CNAs_mets$startpos, prepared_CNAs_mets$endpos))
gr_prepared_CNAs_mets$annot = prepared_CNAs_mets$annot
gr_prepared_CNAs_mets$patient = prepared_CNAs_mets$patient
gr_prepared_CNAs_mets$total_cn = prepared_CNAs_mets$total_cn
gr_prepared_CNAs_mets$rel_cn = prepared_CNAs_mets$rel_cn
gr_prepared_CNAs_mets$CNA_id = prepared_CNAs_mets$CNA_id

ov = findOverlaps(gr_prepared_CNAs_mets, gr_CDS)
gc()

gene_cna_annots = data.frame(
		gene = gr_CDS[subjectHits(ov)]$symbol, 
		annot = gr_prepared_CNAs_mets[queryHits(ov)]$annot,
		patient = gr_prepared_CNAs_mets[queryHits(ov)]$patient,
		total_cn = gr_prepared_CNAs_mets[queryHits(ov)]$total_cn,
		cn_width = width(gr_prepared_CNAs_mets[queryHits(ov)]),
		rel_cn = gr_prepared_CNAs_mets[queryHits(ov)]$rel_cn,
		chr = as.character(seqnames(gr_prepared_CNAs_mets[queryHits(ov)])),
		CNA_id = as.character(gr_prepared_CNAs_mets[queryHits(ov)]$CNA_id),
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
gene_cna_annots_mets = gene_cna_annots
save(gene_cna_annots_mets, file = pff("gene_cna_annots_mets.rsav"))
gc()





##
# percent overlaps of genes and CNAs per patient
##
load(file = pff("prepared_CNAs_mets.rsav"))
load(file = pff("gr_prepared_elements.rsav"))
gr_genes = gr_prepared_elements[grep("^CDSgene", gr_prepared_elements$id)]

genes_with_high_gain_mets = gene_cna_coverage("high_gain", prepared_CNAs_mets, gr_genes, FRACTION_GENE_ALTERED, "high_gain")
gc()
genes_with_full_loss_mets = gene_cna_coverage("full_loss", prepared_CNAs_mets, gr_genes, FRACTION_GENE_ALTERED, "full_loss")
gc()
genes_with_gain_mets = gene_cna_coverage(c("gain", "high_gain"), prepared_CNAs_mets, gr_genes, FRACTION_GENE_ALTERED, "all_gain")
gc()
genes_with_loss_mets = gene_cna_coverage(c("loss", "full_loss"), prepared_CNAs_mets, gr_genes, FRACTION_GENE_ALTERED, "all_loss")
gc()

# remove ChrX from full loss genes
chrX_genes = gr_genes[seqnames(gr_genes) == "chrX"]$symbol
genes_with_full_loss_mets = genes_with_full_loss_mets[!genes_with_full_loss_mets$gene %in% chrX_genes,]
genes_with_loss_mets = genes_with_loss_mets[!genes_with_loss_mets$gene %in% chrX_genes,]
genes_with_full_CNAs_mets = rbind(genes_with_high_gain_mets, genes_with_full_loss_mets, genes_with_gain_mets, genes_with_loss_mets)
save(genes_with_full_CNAs_mets, file = pff("genes_with_full_CNAs_mets.rsav"))