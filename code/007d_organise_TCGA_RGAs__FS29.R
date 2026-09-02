source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")

load(pff("TCGA_validation/TCGA_CNA_patient_lists.rsav"))
load(pff("TCGA_validation/dfr_TCGA_ncSNV_RGA_patient.rsav"))
load(pff("TCGA_validation/dfr_TCGA_SV_RGA_patient.rsav"))
load(pff("TCGA_validation/dnds_RGA_annotations_TCGA.rsav"))
load(pff("TCGA_validation/ERG_pos_patients.rsav"))
load(pff("TCGA_validation/unq_TCGA_samples.rsav"))
load(pff("TCGA_validation/unq_TCGA_patients.rsav"))

TCGA_CNA_patients = lapply(TCGA_CNA_patient_lists, intersect, unq_TCGA_samples)
TCGA_CNA_patients = lapply(TCGA_CNA_patients, gsub_TCGA_sample_patient)
TCGA_CNA_patients = lapply(TCGA_CNA_patients, unique)

TCGA_ncSNV_patients = dfr_TCGA_ncSNV_RGA_patient[dfr_TCGA_ncSNV_RGA_patient$TCGA_sample_id %in% unq_TCGA_samples,]
TCGA_ncSNV_patients$TCGA_patient = gsub_TCGA_sample_patient(TCGA_ncSNV_patients$TCGA_sample_id)
TCGA_ncSNV_patients = split(TCGA_ncSNV_patients$TCGA_patient, TCGA_ncSNV_patients$RGA_id)

TCGA_SV_patients = dfr_TCGA_SV_RGA_patient[dfr_TCGA_SV_RGA_patient$TCGA_sample_id %in% unq_TCGA_samples,]
TCGA_SV_patients$TCGA_patient = gsub_TCGA_sample_patient(TCGA_SV_patients$TCGA_sample_id)
TCGA_SV_patients = split(TCGA_SV_patients$TCGA_patient, TCGA_SV_patients$RGA_id)

TCGA_cdSNV_patients = dnds_RGA_annotations_TCGA[dnds_RGA_annotations_TCGA$sampleID %in% unq_TCGA_samples,]
TCGA_cdSNV_patients$TCGA_patients = gsub_TCGA_sample_patient(TCGA_cdSNV_patients$sampleID)
TCGA_cdSNV_patients = split(TCGA_cdSNV_patients$TCGA_patient, TCGA_cdSNV_patients$gene) 
names(TCGA_cdSNV_patients) = paste(names(TCGA_cdSNV_patients), "SNV_CDS", sep = "__")

TCGA_ETS_patients = intersect(gsub_TCGA_sample_patient(unq_TCGA_samples), ERG_pos_patients)
TCGA_ETS_patients = data.frame(RGA_id = "ETS__SV", TCGA_patient = TCGA_ETS_patients, stringsAsFactors = FALSE)
TCGA_ETS_patients = split(TCGA_ETS_patients$TCGA_patient, TCGA_ETS_patients$RGA_id)



# collect into sets of patients by drivers, create comparable PPCG and TCGA patient lists and matrices

# biallelic RGAs are excluded from validation
load(pff("patient_sets_for_drivers.rsav"))
PPCG_patient_sets_for_drivers = patient_sets_for_drivers[grep("__BI$", names(patient_sets_for_drivers), invert = TRUE)]

# compile TCGA patient lists; add elements that are not found mutated in TCGA
TCGA_patient_sets_for_drivers = c(TCGA_CNA_patients, TCGA_ncSNV_patients, TCGA_SV_patients, TCGA_ETS_patients, TCGA_cdSNV_patients)

nonmut_RGAs_in_TCGA = setdiff(names(PPCG_patient_sets_for_drivers), names(TCGA_patient_sets_for_drivers))
empty_sets =  lapply(nonmut_RGAs_in_TCGA, function(x) c())
names(empty_sets) = nonmut_RGAs_in_TCGA
TCGA_patient_sets_for_drivers = c(TCGA_patient_sets_for_drivers, empty_sets)
TCGA_patient_sets_for_drivers = TCGA_patient_sets_for_drivers[names(PPCG_patient_sets_for_drivers)]

save(PPCG_patient_sets_for_drivers, file = pff("TCGA_validation/PPCG_patient_sets_for_drivers.rsav"))
save(TCGA_patient_sets_for_drivers, file = pff("TCGA_validation/TCGA_patient_sets_for_drivers.rsav"))


# use this function to include empty sets of RGAs as well
patient_lists_to_matrix = function(patient_sets, unique_patients) {

	RGA_mat = do.call(rbind, lapply(names(patient_sets), function(x) 
			0+(unique_patients %in% patient_sets[[x]]) ))
	RGA_mat = t(RGA_mat)
	rownames(RGA_mat) = unique_patients
	colnames(RGA_mat) = names(patient_sets)
	RGA_mat
}

unq_TCGA_patients = gsub_TCGA_sample_patient(unq_TCGA_samples)
# subtle point: RGA clusters trained on variants with min 1 mut
unq_PPCG_patients = unique(unlist(patient_sets_for_drivers))

TCGA_driver_matrix = patient_lists_to_matrix(TCGA_patient_sets_for_drivers, unq_TCGA_patients)
PPCG_driver_matrix = patient_lists_to_matrix(PPCG_patient_sets_for_drivers, unq_PPCG_patients)
save(TCGA_driver_matrix, file = pff("TCGA_validation/TCGA_driver_matrix.rsav"))
save(PPCG_driver_matrix, file = pff("TCGA_validation/PPCG_driver_matrix.rsav"))


# plot alteration freqs in TCGA and PPCG cohorts
load(file = pff("TCGA_validation/PPCG_patient_sets_for_drivers.rsav"))
load(file = pff("TCGA_validation/TCGA_patient_sets_for_drivers.rsav"))
load(file = pff("TCGA_validation/unq_TCGA_samples.rsav"))


load(pff("all_patients.rsav"))
n_TCGA = length(unq_TCGA_samples)
n_PPCG = length(all_patients)


MIN_DRIVER_FREQ = 0.01
MIN_N_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)
MIN_N_PATIENTS
select_RGAs = names(which(sapply(PPCG_patient_sets_for_drivers, length) > MIN_N_PATIENTS))

TCGA_patient_sets_for_drivers = TCGA_patient_sets_for_drivers[select_RGAs]
PPCG_patient_sets_for_drivers = PPCG_patient_sets_for_drivers[select_RGAs]


get_stats = function(dat, n_cohort) {
	do.call(rbind, lapply(names(dat), function(x) 
		data.frame(RGA = x, 
				percent_cohort = length(dat[[x]])/ n_cohort, 
				n_in_cohort = length(dat[[x]]), 
				stringsAsFactors = FALSE)))
}


TCGA_stats = get_stats(TCGA_patient_sets_for_drivers, n_TCGA)
PPCG_stats = get_stats(PPCG_patient_sets_for_drivers, n_PPCG)

TCGA_stats$dataset = "TCGA"
PPCG_stats$dataset = "PPCG"

combined_stats = rbind(TCGA_stats, PPCG_stats)
combined_stats$alteration_type = gsub("(.+)__(.+)", "\\2", combined_stats$RGA)

combined_stats$RGA1 = gsub("_", " ", combined_stats$RGA)
# order by most to least frequent RGAs
stats_for_ordering = combined_stats[combined_stats$dataset == "PPCG",]
stats_for_ordering = stats_for_ordering[order(stats_for_ordering$percent_cohort),]
combined_stats$RGA1 = factor(combined_stats$RGA1, levels = stats_for_ordering$RGA1)

combined_stats$alteration_type = factor(combined_stats$alteration_type, 
		levels = c("loss", "gain", "hAMP", "SNV_CDS", "SNV_NC", "SV"))

plt = ggplot(combined_stats, aes(RGA1, percent_cohort, fill = dataset)) + 
		geom_bar(stat = "identity", position = "dodge") + 
		facet_grid(alteration_type ~ 1, scales = "free", space = "free") +
		plot_theme() + 
		coord_flip()
		
fname = pff("figures/TCGA_RGA_frequencies.pdf")
ggsave(plt, file = fname, height = 10)
file_open_call2(fname)
