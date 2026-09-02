source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(jsonlite)
library(VariantAnnotation)

system(paste("mkdir", pff("TCGA_validation")))

HYPERMUT_COUNT = 3000 * 30


# CNA segments
CNA_file_location = paste0("DATA_USED__", this_timestamp, "/TCGA_PRCA_CNA_data_GATK4__2026-06-12/DATA/")
CNA_files = list.files(CNA_file_location, recursive = TRUE)
CNA_files = CNA_files[grep("_wgs_gdc_realn.cr.igv.reheader.seg.txt$", CNA_files)]

TCGA_CNA_segments = do.call(rbind, lapply(paste0(CNA_file_location, CNA_files), read.delim))
TCGA_CNA_segments = data.frame(
		TCGA_patient =  gsub("(TCGA)-(..)-(....)(.+)", "\\1-\\2-\\3", TCGA_CNA_segments$GDC_Aliquot_ID),
		TCGA_CNA_segments,
		stringsAsFactors = FALSE)
save(TCGA_CNA_segments, file = pff("TCGA_validation/TCGA_CNA_segments.rsav"))


# SV breakpoints
read_sv_files = function(fname) {
	cat(".")
	dat = read.delim(fname)	
	fname_id = strsplit(fname, s = "/")[[1]][[5]]
	dat = data.frame(fname_id = fname_id, dat, stringsAsFactors = FALSE)
	dat
}

get_tcga_id_from_json = function(i, json_records) {
	entity_ids = json_records[i,]$associated_entities[[1]]$entity_submitter_id
	tumor_normal_id = gsub("(TCGA)-(..)-(....)-(..)(.)-(...)-(....)-(..)", "\\4", entity_ids)
	tumor_id = entity_ids[tumor_normal_id == "01"]
	if (nchar(tumor_id) == 0) {
		tumor_id = NA
	}
	tumor_id
}




SV_file_location = paste0("DATA_USED__", this_timestamp, "/TCGA_PRCA_SV_data_Manta__2026-06-12/DATA/")
SV_files = list.files(SV_file_location, recursive = TRUE)
SV_files = SV_files[grep(".somatic.SV.bedpe$", SV_files)]


TCGA_SV_segments = do.call(rbind, lapply(paste0(SV_file_location, SV_files), read_sv_files))

# need to map uuids to TCGA IDs via this metadata file
json_sv_file = paste0("DATA_USED__", this_timestamp, "/TCGA_PRCA_SV_data_Manta__2026-06-12/metadata.repository.2026-06-12(2).json")
json_sv = fromJSON(json_sv_file)


TCGA_SV_map = lapply(1:nrow(json_sv), get_tcga_id_from_json, json_sv)
TCGA_SV_map = sapply(TCGA_SV_map, '[[', 1)
names(TCGA_SV_map) = json_sv$file_name



# add TCGA ID to SV file
TCGA_SV_segments = data.frame(
		TCGA_sample = TCGA_SV_map[TCGA_SV_segments$fname_id], 
		TCGA_patient = gsub("(.+)-(.+)-(.+)-(.+)+", "\\1-\\2-\\3", TCGA_SV_map[TCGA_SV_segments$fname_id]), 
		TCGA_SV_segments, stringsAsFactors = FALSE)
TCGA_SV_segments = TCGA_SV_segments[TCGA_SV_segments$filter == "PASS",]


TCGA_SV_segments$tag = paste(
		TCGA_SV_segments$TCGA_sample, 
		TCGA_SV_segments$TCGA_patient, 
		TCGA_SV_segments$X.chrom1, 
		TCGA_SV_segments$start1,
		TCGA_SV_segments$end1, 
		TCGA_SV_segments$chrom2,
		TCGA_SV_segments$start2, 
		TCGA_SV_segments$end2, 
		TCGA_SV_segments$type, 
		sep = "__"
)

# remove segments into these weird choromosomes
standard_chromosomes = paste0("chr", c(1:22, "X", "Y"))
TCGA_SV_segments = TCGA_SV_segments[TCGA_SV_segments$X.chrom1 %in% standard_chromosomes & TCGA_SV_segments$chrom2 %in% standard_chromosomes,]
save(TCGA_SV_segments, file = pff("TCGA_validation/TCGA_SV_segments.rsav"))

# breakpoints
TCGA_SV_breaks1 = TCGA_SV_segments[, c("TCGA_sample", "TCGA_patient", "X.chrom1", "start1", "end1", "type")]
TCGA_SV_breaks2 = TCGA_SV_segments[, c("TCGA_sample", "TCGA_patient", "chrom2", "start2", "end2", "type")]
colnames(TCGA_SV_breaks1) = colnames(TCGA_SV_breaks2) = c("TCGA_sample", "TCGA_patient", "chr", "start", "end", "SV_type")
TCGA_SV_breaks = rbind(TCGA_SV_breaks1, TCGA_SV_breaks2)
save(TCGA_SV_breaks, file = pff("TCGA_validation/TCGA_SV_breaks.rsav"))


# SNVs, indels; varScan2
read_vcf = function(i, vcf_files, file_loc) {
	
	cat(i, " ")
	
	# read variants from vcf file
	vcf_fname = vcf_files[[i]]
	vcf = readVcf(paste0(file_loc, "/", vcf_fname), "hg38")
	gr_variants = rowRanges(vcf)	
	
	chr = as.character(seqnames(gr_variants))
	pos1 = start(gr_variants)
	pos2 = end(gr_variants)
	ref = as.character(gr_variants$REF)
	alt = sapply(gr_variants$ALT, function(x) paste(as.character(x), collapse = "__"))
	fname_id = gsub("(.+)/(.+)", "\\2", vcf_fname)
	filter = gr_variants$FILTER
	
	data.frame(fname_id, chr, pos1, pos2, ref, alt, filter, stringsAsFactors = FALSE)
}



SNV_Varscan2_file_location = paste0("DATA_USED__", this_timestamp, "/TCGA_PRCA_SNV_data_VarScan2annot_2026-06-12/DATA/")
SNV_Varscan2_files = list.files(SNV_Varscan2_file_location, recursive = TRUE)
SNV_Varscan2_files = SNV_Varscan2_files[grep(".wgs.VarScan2.somatic_annotation.vcf.gz$", SNV_Varscan2_files)]

TCGA_SNVs_VarScan2annot = do.call(rbind, mclapply(1:length(SNV_Varscan2_files), 
		read_vcf, SNV_Varscan2_files, SNV_Varscan2_file_location, mc.cores = 12))
TCGA_SNVs_VarScan2annot = TCGA_SNVs_VarScan2annot[TCGA_SNVs_VarScan2annot$filter == "PASS", ]

# need to map uuids to TCGA IDs via this metadata file
json_snv_Varscan2_file = paste0("DATA_USED__", this_timestamp, "/TCGA_PRCA_SNV_data_VarScan2annot_2026-06-12/metadata.repository.2026-06-12(2).json")
json_snv_Varscan2 = fromJSON(json_snv_Varscan2_file)


TCGA_snv_Varscan2_map = lapply(1:nrow(json_snv_Varscan2), get_tcga_id_from_json, json_snv_Varscan2)
TCGA_snv_Varscan2_map = sapply(TCGA_snv_Varscan2_map, '[[', 1)
names(TCGA_snv_Varscan2_map) = json_snv_Varscan2$file_name

# add TCGA ID to SNVs file
TCGA_SNVs_VarScan2annot = data.frame(
		TCGA_sample = TCGA_snv_Varscan2_map[TCGA_SNVs_VarScan2annot$fname_id], 
		TCGA_patient = gsub("(.+)-(.+)-(.+)-(.+)+", "\\1-\\2-\\3", TCGA_snv_Varscan2_map[TCGA_SNVs_VarScan2annot$fname_id]), 
		TCGA_SNVs_VarScan2annot, 
		stringsAsFactors = FALSE)

# remove hypermutated samples from SNV data
TCGA_hypmut_samples = names(which(sort(table(TCGA_SNVs_VarScan2annot$TCGA_sample)) > HYPERMUT_COUNT))
save(TCGA_hypmut_samples, file = pff("TCGA_validation/TCGA_hypmut_samples.rsav"))
TCGA_SNVs_VarScan2annot_prefilter = TCGA_SNVs_VarScan2annot
save(TCGA_SNVs_VarScan2annot_prefilter, file = pff("TCGA_validation/TCGA_SNVs_VarScan2annot_prefilter.rsav"))		

TCGA_SNVs_VarScan2annot = TCGA_SNVs_VarScan2annot[!TCGA_SNVs_VarScan2annot$TCGA_sample %in% TCGA_hypmut_samples,]
save(TCGA_SNVs_VarScan2annot, file = pff("TCGA_validation/TCGA_SNVs_VarScan2annot.rsav"))

