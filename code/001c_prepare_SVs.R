source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(ActiveDriverWGS)
library(rtracklayer)

# ETS status from a separate curated file
load(file = pff("tracking_sheet.rsav"))
dat_ETS = read.csv(
		paste0("DATA_USED__", this_timestamp, "/main_tracking_sheet_2025-09-30/Main\ Tracking\ Sheet\ ICR\ -\ ETS.csv"),
		stringsAsFactors = FALSE)
dat_ETS = dat_ETS[dat_ETS$WGS_AssayID %in% tracking_sheet$WGS_AssayID & dat_ETS$selected_one_sample_per_donor, ]
save(dat_ETS, file = pff("dat_ETS.rsav"))

# keep ETS-positive patients and top target genes for later curation of RGAs
MIN_N_TO_KEEP_ETS_TGT_GENES = 5
fused_genes =  unlist(dat_ETS[,-1:-9])
fused_genes = fused_genes[fused_genes != "FALSE"]
fused_genes = unlist(strsplit(fused_genes, split = ",|:"))
table_fused_genes = sort(table(fused_genes))
table_fused_genes_top = table_fused_genes[table_fused_genes > MIN_N_TO_KEEP_ETS_TGT_GENES]
ETS_target_genes = names(table_fused_genes_top)

ETS_patients =  dat_ETS[dat_ETS$ETS, "PPCG_Donor_ID"]

save(ETS_target_genes, file = pff("ETS_target_genes.rsav"))
save(ETS_patients, file = pff("ETS_patients.rsav"))

# ETS SV type - deleted or translocated
ETS_SV_status = data.frame(patient = ETS_patients, SV_type = "unknown", stringsAsFactors = FALSE)
ERG_positive_patients = unique(dat_ETS[dat_ETS$ERG, "PPCG_Donor_ID"])
ETV_FLI_positive_patients = unique(dat_ETS[dat_ETS$ETV1 | dat_ETS$ETV4 | dat_ETS$ETV5 | dat_ETS$FLI1, "PPCG_Donor_ID"])
ETS_SV_status[ETS_SV_status$patient %in% ERG_positive_patients, "SV_type"] = "DEL"
ETS_SV_status[ETS_SV_status$patient %in% ETV_FLI_positive_patients, "SV_type"] = "TRA"
save(ETS_SV_status, file = pff("ETS_SV_status.rsav"))

#
# process other SV data
#

fnames = list.files(
		paste0("DATA_USED__", this_timestamp, "/SFTP_SVs_120921/PPCG_SV_Data_Release_15_July_2020__bedpe"), 
		full.names = TRUE)

get_variants = function(fname) {
	cat(".")
	tumor_id = gsub("(.+)/(.+)/(.+)_ppcg_consensus_annotated.somatic.sv.bedpe.gz", "\\3", fname) 
	dat = read.delim(fname, stringsAsFactors = FALSE)	
	
	if (nrow(dat) == 0) {
		return(NULL)
	}
	
	dat$SV_tag = paste(
			paste0("chr", dat$chrom1), dat$start1, dat$end1, 
			paste0("chr", dat$chrom2), dat$start2, dat$end2,
			dat$svclass, dat$samples, sep = "__")
			
	colnames(dat)[colnames(dat) == "chrom1"] = "chr1"
	colnames(dat)[colnames(dat) == "chrom2"] = "chr2"
	colnames(dat)[colnames(dat) == "samples"] = "tumor_id"
	
	dat$chr1 = paste0("chr", dat$chr1)
	dat$chr2 = paste0("chr", dat$chr2)

	standard_chromosomes = paste0("chr", c(1:22, "X", "Y"))
	dat = dat[dat$chr1 %in% standard_chromosomes & dat$chr2 %in% standard_chromosomes,]
	
	dat = data.frame(dat, stringsAsFactors = FALSE)	
}

variants_SV = do.call(rbind, lapply(fnames, get_variants))


# remove blacklisted samples
load(file = pff("blacklisted_tumor_ids.rsav"))
variants_SV = variants_SV[!variants_SV$tumor_id %in% blacklisted_tumor_ids,]

# store these variants for backup before any filtering
variants_SV_prefilter = variants_SV
save(variants_SV_prefilter, file = pff("variants_SV_prefilter.rsav"))


#
# use only the samples from the tracking sheet: one per patient
# 
load(file = pff("tracking_sheet.rsav"))
variants_SV$patient = gsub("._DNA$", "", variants_SV$tumor_id)
variants_SV = variants_SV[variants_SV$tumor_id %in% tracking_sheet$WGS_AssayID, ]
variant_colnames = colnames(variants_SV)
variant_colnames = setdiff(variant_colnames, c("samples", "tumor_id"))
variants_SV = variants_SV[,variant_colnames]


##
# add genes near each breakpoint
##
load(file = pff("gr_prepared_elements.rsav"))
gr_CDSgene = gr_prepared_elements[grep("CDSgene", gr_prepared_elements$id)]
# flanking window to find nearby genes
start(gr_CDSgene) = start(gr_CDSgene) - 3000
end(gr_CDSgene) = end(gr_CDSgene) + 3000


# chr = variants_SV$chr1; start = variants_SV$start1; end = variants_SV$end1; gr_gene_coords = gr_CDSgene
get_genes_per_bp = function(chr, start, end, gr_gene_coords) {
	
	gr_bp = GRanges(chr, IRanges(start, end))
	dfr_bp = as.data.frame(gr_bp)
	dfr_bp$tag = paste(dfr_bp$seqnames, dfr_bp$start, dfr_bp$end, sep = ":")
	dfr_genes = as.data.frame(gr_gene_coords)
	ov = findOverlaps(gr_bp, gr_gene_coords)
	
	res = data.frame(bp = dfr_bp[queryHits(ov), "tag"], gene = dfr_genes[subjectHits(ov), "symbol"], stringsAsFactors= FALSE)
	res = split(res$gene, res$bp)
	res = sapply(res, unique)
	res = sapply(res, paste, collapse = ",")
	res
}


bp1_genes = get_genes_per_bp(variants_SV$chr1, variants_SV$start1, variants_SV$end1, gr_CDSgene)
bp2_genes = get_genes_per_bp(variants_SV$chr2, variants_SV$start2, variants_SV$end2, gr_CDSgene)
variants_SV$bp1_tag = paste(variants_SV$chr1, variants_SV$start1, variants_SV$end1, sep = ":")
variants_SV$bp2_tag = paste(variants_SV$chr2, variants_SV$start2, variants_SV$end2, sep = ":")
variants_SV$genes1 = bp1_genes[variants_SV$bp1_tag]
variants_SV$genes2 = bp2_genes[variants_SV$bp2_tag]
variant_colnames = colnames(variants_SV)
variant_colnames = setdiff(variant_colnames, c("bp1_tag", "bp2_tag"))
variants_SV = variants_SV[,variant_colnames]

save(variants_SV, file = pff("variants_SV.rsav"))





#
# another datasets including only breakpoints
#
SV_bp1 = data.frame(variants_SV[,c("chr1", "start1", "end1", "sv_id", "svclass", "patient", "SV_tag")], anchor = 1, stringsAsFactors = FALSE)
SV_bp2 = data.frame(variants_SV[,c("chr2", "start2", "end2", "sv_id", "svclass", "patient", "SV_tag")], anchor = 2, stringsAsFactors = FALSE)
colnames(SV_bp1) = colnames(SV_bp2) = c("chr", "start", "end", "sv_id", "svclass", "patient", "SV_tag", "anchor")
SV_breakpoints = rbind(SV_bp1, SV_bp2)
prepared_variants_SV = SV_breakpoints
prepared_variants_SV$tag = "indel>X" # for ADWGS
prepared_variants_SV$top_mut_signt = prepared_variants_SV$svclass

gr_prepared_variants = GRanges(prepared_variants_SV$chr,
			                 IRanges(start = prepared_variants_SV$start,
			                         end = prepared_variants_SV$end),
			                 mcols = prepared_variants_SV[,c("patient", "tag", "SV_tag")])

save(gr_prepared_variants, file = pff("gr_prepared_variants__SV.rsav"))
save(prepared_variants_SV, file = pff("prepared_variants_SV.rsav"))