source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(ActiveDriverWGS)
library(rtracklayer)


#
# variants for mets
# first, collect all mets (top purity by patient)
# second, collect de novo mets from JW list
# third, create variant SNV/indel set for all mets list
#

##
# collect all mets for RGA comparison, keep one met per patient based on highest purity. 
##
purity_data = read.csv(
 		paste0("DATA_USED__", this_timestamp, "/SCNA_with_BD_SVs_03_05_2023/Cellularity_Ploidy_Estimates_03_May_2023.csv"),
		stringsAsFactors = FALSE)

tracking_sheet_mets = read.csv(
		paste0("DATA_USED__", this_timestamp, "/main_tracking_sheet_2025-09-30/Main\ Tracking\ Sheet\ ICR\ -\ WGS\ v2.4.csv"),
		stringsAsFactors = FALSE)
		
load(file = pff("blacklisted_tumor_ids.rsav"))
tracking_sheet_mets = tracking_sheet_mets[tracking_sheet_mets$SampleType == "MET",]
tracking_sheet_mets = tracking_sheet_mets[!tracking_sheet_mets$WGS_AssayID %in% blacklisted_tumor_ids, ]
tracking_sheet_mets = merge(tracking_sheet_mets, purity_data, by.x = "WGS_AssayID", by.y = "PPCG_Sample_ID")

# unique mets by patient id, select highest cellularity sample
tracking_sheet_mets_by_patient = split(tracking_sheet_mets, tracking_sheet_mets$PPCG_Donor_ID)
tracking_sheet_mets = do.call(rbind, lapply(tracking_sheet_mets_by_patient, function(x) x[which.max(x$Cellularity)[[1]],] ))

allmets_tumor_ids = tracking_sheet_mets$WGS_AssayID
save(allmets_tumor_ids, file = pff("allmets_tumor_ids.rsav"))


####
# collect de novo mets for the P-M comparison
####
denovomets = read.delim(paste0("DATA_USED__", this_timestamp, 
				"/DeNovoMet_JW_2026-03-26/ppcg_denovomet.txt"), col.names = FALSE, 
		stringsAsFactors= FALSE)
denovomets_tumor_ids = grep("c_DNA$", denovomets[,1], value = TRUE)
save(denovomets_tumor_ids, file = pff("denovomets_tumor_ids.rsav"))

# take SNV/indel variant file, keep only the variants in selected mets samples
load(file = pff("variants_prefilter__SNV_indel.rsav"))
variants_all_mets = variants_prefilter[!variants_prefilter$tumor_id %in% blacklisted_tumor_ids,]
variants_all_mets = variants_all_mets[variants_all_mets$tumor_id %in% allmets_tumor_ids,]
variants_all_mets = variants_all_mets[grep("__", variants_all_mets$callers),]
standard_chromosomes = paste0("chr", c(1:22, "X", "Y"))
variants_all_mets = variants_all_mets[variants_all_mets$chr %in% standard_chromosomes,]
variants_all_mets$patient = gsub("._DNA$", "", variants_all_mets$tumor_id)
variants_mets__SNV_indel = variants_all_mets
save(variants_mets__SNV_indel, file = pff("variants_mets__SNV_indel.rsav"))
