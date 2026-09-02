source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
source(paste0("../bin/", date_tag, "/004xx_functions_for_clinical_analyses.R"))

#
# prepare clinical data: remove missing values; GleasonGroup_4plus defined as 4+, EarlyOnset defined as 55+
#
load(pff("all_patients.rsav"))
patient_clinical_data = read.csv(
		paste0("DATA_USED__", this_timestamp, "/clinical_pathology_RedCap__2024-12-05/PPCGCoreDatav282024D_DATA_2024-12-05_1617.csv"), 
		stringsAsFactors = FALSE)
patient_clinical_data = patient_clinical_data[patient_clinical_data$record_id %in% all_patients,]

patient_clinical_data$Gleason_group = patient_clinical_data$gleason_grade_group
patient_clinical_data$GleasonGroup_3plus = patient_clinical_data$Gleason_group >= 3

# age and survival intervals should be numeric values, NAs added otherwise
patient_clinical_data$patient_age_at_tumour_collection = as.numeric(patient_clinical_data$age_at_tumour_collection)
patient_clinical_data$early_onset = patient_clinical_data$age_at_tumour_collection <= 55

patient_clinical_data$donor_interval_of_last_followup = as.numeric(patient_clinical_data$donor_interval_of_last_followup)
patient_clinical_data$donor_relapse_interval = as.numeric(patient_clinical_data$donor_relapse_interval)
patient_clinical_data$donor_first_mets_interval_from_m0 = as.numeric(patient_clinical_data$donor_first_mets_interval_from_m0)

# add T stage, soft-remove T1
patient_clinical_data$T_stage = gsub("(a|b|c)$", "", patient_clinical_data$path_t_stage)
patient_clinical_data[!is.na(patient_clinical_data$T_stage) & !(patient_clinical_data$T_stage %in% c("T2", "T3", "T4")), "T_stage"] = NA

# add PSA
patient_clinical_data$PSA_log2 = log2(patient_clinical_data$psa_at_tumour_collection)

# look at binary variables and deal with missing values
patient_clinical_data$new_metastatic_biology_indicator = 
		as.character(factor(patient_clinical_data$new_metastatic_biology_indicator, levels = c("mets_biol", "no_mets_biol")))

# remove new_pccg_risk_category, damico risk
patient_clinical_data = patient_clinical_data[,!colnames(patient_clinical_data) %in% c("damico_risk", "de_novo_mets")]

patient_clinical_data[,"death_ind"] = as.character(factor(patient_clinical_data[,"death_ind"], levels = c("dead", "alive")))
patient_clinical_data[,"relapse_ind"] = 
		as.character(factor(patient_clinical_data[,"relapse_ind"], levels = c("relapsed", "no relapse")))
patient_clinical_data[,"mets_ind"] = 
		as.character(factor(patient_clinical_data[,"mets_ind"], levels = c("mets", "no mets")))


# adjust zero time to event
adjust_zero_time_to_event = function(time_col, patient_clinical_data, new_nonzero_time_value) {
	which_to_adjust = which(!is.na(patient_clinical_data[, time_col]) & patient_clinical_data[, time_col] == 0)
	patient_clinical_data[which_to_adjust, time_col] = new_nonzero_time_value
	patient_clinical_data
	
}

# soft-remove ignore time =0 patients by setting the new value to NA
new_nonzero_time_value = NA
patient_clinical_data = adjust_zero_time_to_event("donor_interval_of_last_followup", patient_clinical_data, new_nonzero_time_value)
patient_clinical_data = adjust_zero_time_to_event("donor_relapse_interval", patient_clinical_data, new_nonzero_time_value)
patient_clinical_data = adjust_zero_time_to_event("donor_first_mets_interval_from_m0", patient_clinical_data, new_nonzero_time_value)

# censor time to event analyses: 5year overall survival; 10year time to recurrence, 10year time to mets
patient_clinical_data = censor_time_to_event(
		time_col = "donor_interval_of_last_followup", 
		event_col = "death_ind", 
		event_value = "dead", 
		no_event_value = "alive", 
		n_years_censor = 5, 
		new_time_col = "OS_time_5year", 
		new_event_col = "OS_event_5year", 
		patient_clinical_data)

patient_clinical_data = censor_time_to_event(
		time_col = "donor_relapse_interval", 
		event_col = "relapse_ind", 
		event_value = "relapsed", 
		no_event_value = "no relapse", 
		n_years_censor = 10, 
		new_time_col = "relapse_time_10year", 
		new_event_col = "relapse_event_10year", 
		patient_clinical_data)	

patient_clinical_data = censor_time_to_event(
		time_col = "donor_first_mets_interval_from_m0", 
		event_col = "mets_ind", 
		event_value = "mets", 
		no_event_value = "no mets", 
		n_years_censor = 10, 
		new_time_col = "mets_time_10year", 
		new_event_col = "mets_event_10year", 
		patient_clinical_data)


# add percent-genome-altered as a clinical variable
load(file = pff("dfr_PGA_estimates.rsav"))
patient_clinical_data = merge(patient_clinical_data, dfr_PGA_estimates[, c("patient", "PGA")], 
		by.x = "record_id", by.y = "patient", all.x = TRUE)
		

# add SNV count as a clinical variable
load(pff("prepared_variants__SNV_indel.rsav"))
n_variants_tab = table(prepared_variants$patient)
n_variants_dfr = data.frame(record_id = names(n_variants_tab), GTMB_log1p = log1p(as.numeric(n_variants_tab)), 
		stringsAsFactors = FALSE)
patient_clinical_data = merge(patient_clinical_data, n_variants_dfr, all.x = TRUE)
patient_clinical_data[is.na(patient_clinical_data$GTMB_log1p), "GTMB_log1p"] = 0
  
# if ancestry is unknown, assign ancestry as "unknown"
patient_clinical_data$predicted_ancestry [patient_clinical_data$predicted_ancestry == ""] = NA

# rename main id to "patient" to consolidate with other datasets
colnames(patient_clinical_data)[colnames(patient_clinical_data) == "record_id"] = "patient"



# germline_risk_genes
germline_risk_genes = read.delim(
		paste0("DATA_USED__", this_timestamp, "/germline_risk_genes_Daria_2026-06-18/risk_genes_jun2026.tsv"), 
		stringsAsFactors = FALSE)

BRCA1_patients = names(which(unlist(germline_risk_genes[germline_risk_genes$X == "BRCA1", colnames(germline_risk_genes) != "X"]) == 1))
BRCA2_patients = names(which(unlist(germline_risk_genes[germline_risk_genes$X == "BRCA2", colnames(germline_risk_genes) != "X"]) == 1))
non_BRCA_table = germline_risk_genes[, !colnames(germline_risk_genes) %in% c("BRCA1", "BRCA2")]
other_risk_gene_patients = names(which(apply(non_BRCA_table[colnames(non_BRCA_table) != "X"], 2, sum) > 1))

BRCA1_patients = intersect(BRCA1_patients, patient_clinical_data$patient)
BRCA2_patients = intersect(BRCA2_patients, patient_clinical_data$patient)
other_risk_gene_patients = intersect(other_risk_gene_patients, patient_clinical_data$patient)

patient_clinical_data$germline_risk_gene = "none"
patient_clinical_data[patient_clinical_data$patient %in% other_risk_gene_patients, "germline_risk_gene"] = "other_gene"
patient_clinical_data[patient_clinical_data$patient %in% BRCA1_patients, "germline_risk_gene"] = "BRCA1"
patient_clinical_data[patient_clinical_data$patient %in% BRCA2_patients, "germline_risk_gene"] = "BRCA2"

save(patient_clinical_data, file = pff("patient_clinical_data.rsav"))