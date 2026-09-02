source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(ActiveDriverWGS)
library(rtracklayer)

system(paste("mkdir", pff("figures/")))

tracking_sheet = read.csv(
		paste0("DATA_USED__", this_timestamp, "/main_tracking_sheet_2025-09-30/Main\ Tracking\ Sheet\ ICR\ -\ WGS\ v2.4.csv"),
		stringsAsFactors = FALSE)
tracking_sheet = tracking_sheet[
		!is.na(tracking_sheet$selected_one_sample_per_donor) & 
		tracking_sheet$selected_one_sample_per_donor, ]
tracking_sheet = tracking_sheet[tracking_sheet$SampleType == "TUM",]
save(tracking_sheet, file = pff("tracking_sheet.rsav"))


new_blacklist = read.csv(
		paste0("DATA_USED__", this_timestamp, 
				"/main_tracking_sheet_2025-09-30/Main Tracking Sheet\ ICR\ - Exclusion_list v1.0.csv"), 
		stringsAsFactors= FALSE)
new_blacklist = unique(new_blacklist$WGS_AssayID)
blacklisted_tumor_ids = new_blacklist
save(blacklisted_tumor_ids, file = pff("blacklisted_tumor_ids.rsav"))


# which samples should we select? use only the samples from the tracking sheet: one per patient
distinct_tumor_ids = tracking_sheet$WGS_AssayID
save(distinct_tumor_ids, file = pff("distinct_tumor_ids.rsav"))


#
# indel and SNV variant calls
#
fname_snvs = paste0("DATA_USED__", this_timestamp, "/VARIANT_ANNOTATIONS__2023-09-25_release4__SNV+indel/VARIANT_ANNOTATIONS_MERGED_TABLE__2023-09-25_release4__Filtered_SNV_VCFs_20_April_2021.txt.gz")
variants_snv = read.delim(fname_snvs, stringsAsFactors = FALSE)

fname_indels = paste0("DATA_USED__", this_timestamp, "/VARIANT_ANNOTATIONS__2023-09-25_release4__SNV+indel/VARIANT_ANNOTATIONS_MERGED_TABLE__2023-09-25_release4__PPCG_indels_16_Dec_2022_Tier_1.txt.gz")
variants_indel = read.delim(fname_indels, stringsAsFactors = FALSE)

variants = rbind(variants_snv, variants_indel)

# store these variants for backup before any filtering
variants_prefilter = variants
save(variants_prefilter, file = pff("variants_prefilter__SNV_indel.rsav"))

rm(variants_prefilter, variants_snv, variants_indel)
gc()

# remove variants from one pipeline only, and non-standard chromosomes
variants = variants[grep("__", variants$callers),]

standard_chromosomes = paste0("chr", c(1:22, "X", "Y"))
variants = variants[variants$chr %in% standard_chromosomes,]

# keep only patient/donor id here
variants$patient = gsub("._DNA$", "", variants$tumor_id)

#
# keep in variants only the distinct samples selected from the tracking sheet
#
variants = variants[variants$tumor_id %in% tracking_sheet$WGS_AssayID, ]
save(variants, file = pff("variants__SNV_indel.rsav"))


all_patients = unique(tracking_sheet$WGS_AssayID)
all_patients = gsub("._DNA$", "", all_patients)
save(all_patients, file = pff("all_patients.rsav"))



#
# prepare variants object for parallel running of ADWGS
#
	
filter_hyper_MB = 30
this_genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens

prepared_variants = format_muts(variants, 
				filter_hyper_MB = filter_hyper_MB, 
				this_genome = this_genome)

##
# assign mut. signatures 
##
add_mutation_signatures = function(maf, mut_signt, sig_colname, keep_NA_sigs = FALSE) {
	
	by_mut_signt = c("chr", "start", "end", "ref", "alt", "Sample_ID")
	
	maf1 = merge(maf, mut_signt[, c(by_mut_signt, sig_colname)], 
		by.x = c("chr", "pos1", "pos2", "ref", "alt", "tumor_id"),
		by.y = by_mut_signt, 
		all.x = keep_NA_sigs) # keep some mutations where signature is NA
		
	colnames(maf1)[colnames(maf1) == sig_colname] = "top_mutsig"
	cat(paste0("merging with sigs: before=", nrow(maf), 
			"; after=", nrow(maf1), 
			"; equal=", nrow(maf1)==nrow(maf), 
			"; keep_NA_sigs=", keep_NA_sigs, "\n"))
	
	maf1
}
		
SNV_signt = read.delim(
		paste0("DATA_USED__", this_timestamp, "/mut_signature_annotations_Kevin_2024-05-09/SNVs_signature_20240509.tsv.gz"), 
		stringsAsFactors = FALSE)
indel_signt = read.delim(
		paste0("DATA_USED__", this_timestamp, "/mut_signature_annotations_Kevin_2024-05-09/Indel_signature_20240509.tsv.gz"), 
		stringsAsFactors = FALSE)

# assign each variant the top mutational signature
prepared_variants_snv = prepared_variants[prepared_variants$tag != "indel>X",]
prepared_variants_indel = prepared_variants[prepared_variants$tag == "indel>X",]
prepared_variants_snv = add_mutation_signatures(prepared_variants_snv, SNV_signt, "Top_Signature", keep_NA_sigs = TRUE)
prepared_variants_indel = add_mutation_signatures(prepared_variants_indel, indel_signt, "Top_Signature", keep_NA_sigs = TRUE)


#create a genomeRanges object where the trinucleotide tag and patient_ID are indicated
prepared_variants = rbind(prepared_variants_snv, prepared_variants_indel)
save(prepared_variants, file = pff("prepared_variants__SNV_indel.rsav"))

gr_prepared_variants = GRanges(prepared_variants$chr,
		IRanges(
			start = prepared_variants$pos1,
			end = prepared_variants$pos2),
		mcols = prepared_variants[,c("patient", "tag")])

save(gr_prepared_variants, file = pff("gr_prepared_variants__SNV_indel.rsav"))




#
# compile coverage information for WGS data
#
load(file = pff("tracking_sheet.rsav"))

wgs_qa_table = read.delim(
		paste0("DATA_USED__", this_timestamp, "/seq_coverage_2026-06-05/Sanger_QC_Metrics_01_June_2020.txt"),
		stringsAsFactors = FALSE)

wgs_qa_table_denmark = read.delim(
		paste0("DATA_USED__", this_timestamp, "/seq_coverage_2026-06-05/ppcg_sanger_metrics_with_PPCG_ID_DK.tsv"),
		stringsAsFactors = FALSE)

wgs_qa_table = wgs_qa_table[,c("PPCG_Sample_ID", "Tumour.total.depth", "Normal.total.depth")]
wgs_qa_table_denmark = wgs_qa_table_denmark[,c("PPCG_Assay_ID", "Tumour.total.depth", "Normal.total.depth")]

colnames(wgs_qa_table_denmark) = colnames(wgs_qa_table) = c("tumor_id", "tumor_total_depth", "normal_total_depth")

WGS_coverage_tumor_normal = rbind(wgs_qa_table, wgs_qa_table_denmark)
WGS_coverage_tumor_normal = WGS_coverage_tumor_normal[WGS_coverage_tumor_normal$tumor_id %in% tracking_sheet$WGS_AssayID,]

save(WGS_coverage_tumor_normal, file = pff("WGS_coverage_tumor_normal.rsav"))
