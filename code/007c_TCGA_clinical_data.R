source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
source(paste0("../bin/", date_tag, "/004xx_functions_for_clinical_analyses.R"))

# select unique samples for each of TCGA tumors
load(pff("TCGA_validation/TCGA_SV_segments.rsav"))
load(pff("TCGA_validation/TCGA_SNVs_VarScan2annot.rsav"))
load(pff("TCGA_validation/TCGA_CNA_segments.rsav"))
load(pff("TCGA_validation/TCGA_hypmut_samples.rsav"))

all_samples = unique(c(TCGA_SV_segments$TCGA_sample, TCGA_SNVs_VarScan2annot$TCGA_sample, TCGA_CNA_segments$GDC_Aliquot_ID))
all_samples = setdiff(all_samples, TCGA_hypmut_samples)

tumor_to_sample = data.frame(
		tumor_id = gsub_TCGA_sample_patient(all_samples), 
		sample_id = all_samples, 
		stringsAsFactors = FALSE)
tumor_to_sample = split(tumor_to_sample$sample_id, tumor_to_sample$tumor_id)

# take the alphabetically first sample
TCGA_tumor_to_sample = sapply(tumor_to_sample, function(x) sort(x)[1])
TCGA_tumor_to_sample = data.frame(
		TCGA_sample = TCGA_tumor_to_sample, 
		TCGA_patient = names(TCGA_tumor_to_sample), 
		stringsAsFactors = FALSE)
rownames(tumor_to_sample) = NULL
save(TCGA_tumor_to_sample, file = pff("TCGA_validation/TCGA_tumor_to_sample.rsav"))


unq_TCGA_samples = unique(TCGA_tumor_to_sample$TCGA_sample)
unq_TCGA_patients = gsub_TCGA_sample_patient(unq_TCGA_samples)
save(unq_TCGA_samples, file = pff("TCGA_validation/unq_TCGA_samples.rsav"))
save(unq_TCGA_patients, file = pff("TCGA_validation/unq_TCGA_patients.rsav"))


# liu 2018 for consolidated tcga clinical stats
TCGA_clinical_data = read.csv(
		paste0("DATA_USED__", this_timestamp, "/TCGA_PRCA_clinical_Liu2018_2026-06-17/mmc1(2).csv"), 
		stringsAsFactors = FALSE)
TCGA_clinical_data = TCGA_clinical_data[TCGA_clinical_data$type == "PRAD",]
TCGA_clinical_data = TCGA_clinical_data[TCGA_clinical_data$bcr_patient_barcode %in% unq_TCGA_patients,]
TCGA_clinical_data$patient = TCGA_clinical_data$bcr_patient_barcode
TCGA_clinical_data$time_PFI = as.numeric(TCGA_clinical_data$PFI.time)
TCGA_clinical_data$status_PFI = as.numeric(TCGA_clinical_data$PFI)
TCGA_clinical_data = TCGA_clinical_data[
		!is.na(TCGA_clinical_data$time_PFI) &
		!is.na(TCGA_clinical_data$status_PFI),]

TCGA_clinical_data = censor_time_to_event(
		time_col = "time_PFI", 
		event_col = "status_PFI", 
		event_value = 1, 
		no_event_value = 0, 
		n_years_censor = 10, 
		new_time_col = "time_PFI_10year", 
		new_event_col = "status_PFI_10year", 
		TCGA_clinical_data)
		
# PFI to show in years instead of days
TCGA_clinical_data$time_PFI_10year = TCGA_clinical_data$time_PFI_10year / 365
TCGA_clinical_data$age = as.numeric(TCGA_clinical_data$age_at_initial_pathologic_diagnosis)
# early-set <= 55
TCGA_clinical_data$early_onset = TCGA_clinical_data$age <= 55
save(TCGA_clinical_data, file = pff("TCGA_validation/TCGA_clinical_data.rsav"))
