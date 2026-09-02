source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
source(paste0("../bin/", date_tag, "/001xx_functions_for_CNA_processing.R"))
library("dndscv")

# SNVs, indels: mapping using dndscv
load(pff("TCGA_validation/TCGA_SNVs_VarScan2annot.rsav"))
mutations = TCGA_SNVs_VarScan2annot
mutations = mutations[, c("TCGA_sample", "chr", "pos1", "ref", "alt")]
colnames(mutations) = c("sampleID", "chr", "pos", "ref", "alt")
mutations$pos = as.numeric(mutations$pos)
mutations$chr = gsub("chr", "", mutations$chr)
dnds_results_TCGA = dndscv(mutations, refdb = "hg38")
save(dnds_results_TCGA, file = pff("TCGA_validation/dnds_results_TCGA.rsav"))

load(pff("patient_sets_for_drivers.rsav"))
CDS_drivers = grep("CDS", names(patient_sets_for_drivers), value = T)
CDS_drivers = gsub("(.+)__SNV_CDS", "\\1", CDS_drivers)
# use our common function to remove silent mutations and re-label indels
dnds_RGA_annotations_TCGA = dnds_results_TCGA$annotmuts[dnds_results_TCGA$annotmuts$gene %in% CDS_drivers,]
dnds_RGA_annotations_TCGA = dnds_variant_filter_label(dnds_RGA_annotations_TCGA)
save(dnds_RGA_annotations_TCGA, file = pff("TCGA_validation/dnds_RGA_annotations_TCGA.rsav"))



# SVs and non-coding SNVs/indels

# take all the elements of this RGA from ActiveDriverWGS analysis, merge into consecutive regions, perform liftover on the merged regions
get_RGA_regions = function(rga_id, results_signf_merged_annot, chain_hg19_hg38, gr_prepared_elements) {
	
	cat(".")
	mut_type = gsub("(.+)__(.+)", "\\2", rga_id)
	RGA = gsub("(.+)__(.+)", "\\1", rga_id) 
	
	element_ids_altered = results_signf_merged_annot[
			results_signf_merged_annot$annots_MAIN == RGA & 
			results_signf_merged_annot$mut_type == mut_type, "element_ids"]
	element_ids_altered = unique(unlist(strsplit(element_ids_altered, split = ",")))

	gr_elements_altered_hg19 = gr_prepared_elements[gr_prepared_elements$id %in% element_ids_altered]
	gr_elements_altered_hg19 = reduce(gr_elements_altered_hg19)
	
	gr_elements_altered_hg38 = reduce(unlist(liftOver(gr_elements_altered_hg19, chain_hg19_hg38)))
	
	gr_elements_altered_hg38$rga_id = gr_elements_altered_hg19$rga_id = rga_id
	gr_elements_altered_hg38$genome = "hg38"
	gr_elements_altered_hg19$genome = "hg19"
	
	results = c(gr_elements_altered_hg38, gr_elements_altered_hg19)
	
	return(results)
}


# take multi-region merged drivers, their constituent genomic elements
load(pff("results_signf_merged_annot.rsav"))
results_signf_merged_annot = results_signf_merged_annot[
		!is.na(results_signf_merged_annot$element_type) & results_signf_merged_annot$element_type != "CDS",]

# convert hg19 coordinates of consitutent elements to hg38 in which the SV and ncSNV coordinates are presented in TCGA
load(pff("gr_prepared_elements.rsav"))
chain_hg19_hg38 = import.chain(paste0("DATA_USED__", this_timestamp, "/hg19ToHg38.over.chain"))

# SV and ncSNV RGAs follow the same logic; ETS comes from a different table
load(pff("patient_sets_for_drivers.rsav"))
RGA_ids = setdiff(grep("__SV|__SNV_NC", names(patient_sets_for_drivers), value = TRUE), "ETS__SV")

gr_RGA_elements_conv = do.call(c, lapply(RGA_ids, get_RGA_regions, results_signf_merged_annot, chain_hg19_hg38, gr_prepared_elements))
gr_RGA_elements_hg38 = gr_RGA_elements_conv[gr_RGA_elements_conv$genome == "hg38"]
gr_RGA_elements_hg19 = gr_RGA_elements_conv[gr_RGA_elements_conv$genome == "hg19"]

total_width_hg19 = by(gr_RGA_elements_hg19, gr_RGA_elements_hg19$rga_id, function(x) sum(width(x)))
total_width_hg38 = by(gr_RGA_elements_hg38, gr_RGA_elements_hg38$rga_id, function(x) sum(width(x)))

# take ETS SVs : ERG - status from a dedicated table in cbioportal
ETS_status = read.delim(paste0("DATA_USED__", this_timestamp, "/TCGA_PRCA_ERG_cbioportal_2026-06-17/alterations_across_samples.tsv"), 
		stringsAsFactors = FALSE)
ERG_pos_patients = unique(ETS_status[ETS_status$Altered == 1 & !ETS_status$ERG %in% c("AMP", "HOMDEL"), "Patient.ID"])
save(ERG_pos_patients, file = pff("TCGA_validation/ERG_pos_patients.rsav"))


# overlap with variant sets, SVs and SNVs/indels separately
# TCGA SVs mapped to RGA elements, after hg38 mapping of elements
load(pff("TCGA_validation/TCGA_SV_breaks.rsav"))
gr_SV_elements = gr_RGA_elements_hg38[grep("__SV$", gr_RGA_elements_hg38$rga_id)]
gr_SV_breaks = GRanges(TCGA_SV_breaks$chr, IRanges(TCGA_SV_breaks$start, TCGA_SV_breaks$end), 
		sample_id = TCGA_SV_breaks$TCGA_sample, SV_type = TCGA_SV_breaks$SV_type)
ov_SV = findOverlaps(gr_SV_elements, gr_SV_breaks)
dfr_TCGA_SV_RGA_patient = data.frame(
		RGA_id = gr_SV_elements$rga_id[queryHits(ov_SV)], 
		TCGA_sample_id = gr_SV_breaks$sample_id[subjectHits(ov_SV)], 
		SV_type = gr_SV_breaks$SV_type[subjectHits(ov_SV)], 
		stringsAsFactors = FALSE)
save(dfr_TCGA_SV_RGA_patient, file = pff("TCGA_validation/dfr_TCGA_SV_RGA_patient.rsav"))


# TCGA noncoding SNVs/indels mapped to RGA elements, after hg38 mapping of elements
load(pff("TCGA_validation/TCGA_SNVs_VarScan2annot.rsav"))
gr_SNV_elements = gr_RGA_elements_hg38[grep("__SNV_NC$", gr_RGA_elements_hg38$rga_id)]
gr_SNVs = GRanges(TCGA_SNVs_VarScan2annot$chr, IRanges(TCGA_SNVs_VarScan2annot$pos1, TCGA_SNVs_VarScan2annot$pos2), 
		TCGA_sample = TCGA_SNVs_VarScan2annot$TCGA_sample)
ov_SNV = findOverlaps(gr_SNV_elements, gr_SNVs)
dfr_TCGA_ncSNV_RGA_patient = data.frame(
		RGA_id = gr_SNV_elements$rga_id[queryHits(ov_SNV)], 
		TCGA_sample_id = gr_SNVs$TCGA_sample[subjectHits(ov_SNV)], 
		stringsAsFactors = FALSE)
save(dfr_TCGA_ncSNV_RGA_patient, file = pff("TCGA_validation/dfr_TCGA_ncSNV_RGA_patient.rsav"))



# CNAs - annotate segments based on ploidy and the same thresholds as PPCG
chain_hg19_hg38 = import.chain(paste0("DATA_USED__", this_timestamp, "/hg19ToHg38.over.chain"))

load(pff("patient_sets_for_drivers.rsav"))
CNA_drivers = grep("__gain$|__loss$|__hAMP$", names(patient_sets_for_drivers), value = T)
CNA_drivers = gsub("(.+)__gain$|__loss$|__hAMP$", "\\1", CNA_drivers)
CNA_drivers = unique(CNA_drivers)
load(pff("gr_prepared_elements.rsav"))

# take all the elements of this RGA from ActiveDriverWGS analysis, merge into consecutive regions, perform liftover on the merged regions
# convert each gene using liftOver to hg38. we use gene coordinates rather than exons alone.
liftover_gene_coords = function(gene, gr_prepared_elements, chain_hg19_hg38) {
	
	cat(".")
	
	gr_gene_hg19 = gr_prepared_elements[gr_prepared_elements$symbol == gene & grepl("^CDSgene::", gr_prepared_elements$id)]
	gr_gene_hg19 = reduce(gr_gene_hg19)
	gr_gene_hg38 = reduce(unlist(liftOver(gr_gene_hg19, chain_hg19_hg38)))
	
	gr_gene_hg38$gene = gr_gene_hg19$gene = gene
	gr_gene_hg38$genome = "hg38"
	gr_gene_hg19$genome = "hg19"
	
	results = c(gr_gene_hg38, gr_gene_hg19)
	
	return(results)
}

gr_CNA_RGAs_conv = do.call(c, lapply(CNA_drivers, liftover_gene_coords, gr_prepared_elements, chain_hg19_hg38))
gr_CNA_RGAs_hg38 = gr_CNA_RGAs_conv[gr_CNA_RGAs_conv$genome == "hg38"]
gr_CNA_RGAs_hg19 = gr_CNA_RGAs_conv[gr_CNA_RGAs_conv$genome == "hg19"]

total_width_hg19 = by(gr_CNA_RGAs_hg19, gr_CNA_RGAs_hg19$gene, function(x) sum(width(x)))
total_width_hg38 = by(gr_CNA_RGAs_hg38, gr_CNA_RGAs_hg38$gene, function(x) sum(width(x)))

save(gr_CNA_RGAs_hg38, file = pff("TCGA_validation/gr_CNA_RGAs_hg38.rsav"))
save(gr_CNA_RGAs_hg19, file = pff("TCGA_validation/gr_CNA_RGAs_hg19.rsav"))





load(pff("TCGA_validation/TCGA_CNA_segments.rsav"))
# NB! removing small segments supported by 5 or fewer probes
TCGA_CNA_segments = TCGA_CNA_segments[TCGA_CNA_segments$Num_Probes > 5,]


loss_cutoff = -0.2
gain_cutoff = 0.4
high_gain_cutoff = 0.8

TCGA_CNA_segments$annot = "bal"
TCGA_CNA_segments[TCGA_CNA_segments$Segment_Mean <= loss_cutoff, "annot"] = "loss"
TCGA_CNA_segments[TCGA_CNA_segments$Segment_Mean >= gain_cutoff, "annot"] = "gain"
TCGA_CNA_segments[TCGA_CNA_segments$Segment_Mean >= high_gain_cutoff, "annot"] = "high_gain"
TCGA_CNA_segments_annotated = TCGA_CNA_segments
save(TCGA_CNA_segments_annotated, file = pff("TCGA_validation/TCGA_CNA_segments_annotated.rsav"))



# find overlaps of CNA regions and gene coordinates; compute percentages of gene altered by CNA class
# need to relabel columns in data frame to fit existing functions
colnames(TCGA_CNA_segments_annotated) = c("patient1", "patient", "chr", "startpos", "endpos", "num_probes", "rel_cn", "annot")
gr_TCGA_genes = gr_CNA_RGAs_hg38
gr_TCGA_genes$symbol = gr_TCGA_genes$gene



# percent overlaps of genes and CNAs per patient
genes_with_high_gain = gene_cna_coverage("high_gain", TCGA_CNA_segments_annotated, gr_TCGA_genes, FRACTION_GENE_ALTERED, "high_gain")
genes_with_gain = gene_cna_coverage(c("gain", "high_gain"), TCGA_CNA_segments_annotated, gr_TCGA_genes, FRACTION_GENE_ALTERED, "all_gain")
genes_with_loss = gene_cna_coverage(c("loss", "full_loss"), TCGA_CNA_segments_annotated, gr_TCGA_genes, FRACTION_GENE_ALTERED, "all_loss")

# remove ChrX from full loss genes
chrX_genes = gr_TCGA_genes[seqnames(gr_TCGA_genes) == "chrX"]$symbol
genes_with_loss = genes_with_loss[!genes_with_loss$gene %in% chrX_genes,]
TCGA_genes_with_full_CNAs = rbind(genes_with_high_gain, genes_with_gain, genes_with_loss)
save(TCGA_genes_with_full_CNAs, file = pff("TCGA_validation/TCGA_genes_with_full_CNAs.rsav"))


# patient sets by RGAs in TCGA
load(pff("patient_sets_for_drivers.rsav"))
CNA_drivers_by_alt_type = grep("__gain$|__loss$|__hAMP$", names(patient_sets_for_drivers), value = T)
CNA_drivers_by_alt_type = unique(CNA_drivers_by_alt_type)
CNA_drivers_by_alt_type = data.frame(
		do.call(rbind, strsplit(CNA_drivers_by_alt_type, split = "__")), 
		RGA_id = CNA_drivers_by_alt_type, 
		stringsAsFactors = FALSE)
colnames(CNA_drivers_by_alt_type) = c("gene", "alteration", "RGA_id")


which_alterations = list("gain" = "all_gain", "loss" = "all_loss", "hAMP" = "high_gain")

get_patients_per_altered_RGA = function(i, CNA_drivers_by_alt_type, TCGA_genes_with_full_CNAs, which_alterations) {

	gene = CNA_drivers_by_alt_type[i, "gene"]
	altr = CNA_drivers_by_alt_type[i, "alteration"]
	select_CNA_type = which_alterations[[altr]]
	
	CNAs_involved = TCGA_genes_with_full_CNAs[
			TCGA_genes_with_full_CNAs$tag == select_CNA_type & TCGA_genes_with_full_CNAs$gene == gene, ]
	patients_involved = unique(CNAs_involved$patient)
	patients_involved
}

TCGA_CNA_patient_lists = lapply(1:nrow(CNA_drivers_by_alt_type), 
		get_patients_per_altered_RGA, CNA_drivers_by_alt_type, TCGA_genes_with_full_CNAs, which_alterations)
names(TCGA_CNA_patient_lists) = CNA_drivers_by_alt_type$RGA_id
save(TCGA_CNA_patient_lists, file = pff("TCGA_validation/TCGA_CNA_patient_lists.rsav"))
