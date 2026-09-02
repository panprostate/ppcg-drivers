source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")

##
# SVs for all mets
##

load(file = pff("variants_SV_prefilter.rsav"))
load(file = pff("allmets_tumor_ids.rsav"))
variants_SV_prefilter = variants_SV_prefilter[variants_SV_prefilter$tumor_id %in% allmets_tumor_ids,]
variants_SV_prefilter$patient = gsub("._DNA", "", variants_SV_prefilter$tumor_id)

SV_met_bp1 = data.frame(variants_SV_prefilter[,c("chr1", "start1", "end1", "sv_id", "svclass", "patient", "SV_tag")], anchor = 1, 
		stringsAsFactors = FALSE)
SV_met_bp2 = data.frame(variants_SV_prefilter[,c("chr2", "start2", "end2", "sv_id", "svclass", "patient", "SV_tag")], anchor = 2, 
		stringsAsFactors = FALSE)
colnames(SV_met_bp1) = colnames(SV_met_bp2) = c("chr", "start", "end", "sv_id", "svclass", "patient", "SV_tag", "anchor")
SV_met_breakpoints = rbind(SV_met_bp1, SV_met_bp2)

variants_mets__SV = SV_met_breakpoints
variants_mets__SV$tag = "indel>X"
variants_mets__SV$top_mut_signt = variants_mets__SV$svclass

save(variants_mets__SV, file = pff("variants_mets__SV.rsav"))

##
# ETS status for all mets
## 
dat_ETS = read.csv(
		paste0("DATA_USED__", this_timestamp, "/main_tracking_sheet_2025-09-30/Main\ Tracking\ Sheet\ ICR\ -\ ETS.csv"),
		stringsAsFactors = FALSE)

variants_mets__ETS = dat_ETS[dat_ETS$WGS_AssayID %in% allmets_tumor_ids,]
patients_mets_ETS = variants_mets__ETS[variants_mets__ETS$ETS, "PPCG_Donor_ID"]
save(patients_mets_ETS, file = pff("patients_mets_ETS.rsav"))