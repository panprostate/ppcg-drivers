source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(patchwork)
library(gtools)

FDR_CUTOFF = 0.05
NULL_ANALYSIS_TAG = "PPCG__No_covar"
load(pff("all_patients.rsav"))
MIN_DRIVER_FREQ = 0.01
MIN_MUT_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)

load(file = pff(c("driver_clinical_stats__", NULL_ANALYSIS_TAG, ".rsav")))
clin_vars = c("early_onset_55")
driver_clinstat = driver_clinical_stats[driver_clinical_stats$clin_var %in% clin_vars,]
driver_clinstat = driver_clinstat[driver_clinstat$n_mut >= MIN_MUT_PATIENTS,]
driver_clinstat$fdr = p.adjust(driver_clinstat$pval, method = "fdr")
driver_clinstat = driver_clinstat[driver_clinstat$fdr < FDR_CUTOFF,]


load(file = pff("patient_clinical_data.rsav"))
load(file = pff("patient_sets_for_drivers.rsav"))
EO_patients = patient_clinical_data[patient_clinical_data$early_onset, "patient"]
EO_patients = setdiff(EO_patients, NA)
LO_patients = patient_clinical_data[!patient_clinical_data$early_onset, "patient"]
LO_patients = setdiff(LO_patients, NA)
global_fraction_EO = length(EO_patients) / (length(EO_patients) + length(LO_patients))

mut_EO_patient_counts = sapply(lapply(patient_sets_for_drivers, intersect, EO_patients), length)
mut_LO_patient_counts = sapply(lapply(patient_sets_for_drivers, intersect, LO_patients), length)
EO_LO_count_labels = paste0(mut_LO_patient_counts[names(mut_EO_patient_counts)], "\n", mut_EO_patient_counts[names(mut_EO_patient_counts)])
names(EO_LO_count_labels) = names(mut_EO_patient_counts)

driver_clinstat$n_mut_LO_patients = mut_LO_patient_counts[driver_clinstat$id]
driver_clinstat$n_mut_EO_patients = mut_EO_patient_counts[driver_clinstat$id]
driver_clinstat$count_label = EO_LO_count_labels[driver_clinstat$id]
driver_clinstat$count_label[duplicated(driver_clinstat$id)] = NA
driver_clinstat$id = gsub("__", " ", driver_clinstat$id)

driver_clinstat = driver_clinstat[,colnames(driver_clinstat) != "patients_mutated_concat"]
driver_clinstat_melt = melt(driver_clinstat, 
		id.var = c("id", "pval", "effect", "effect_lo", "effect_hi", "n_mut", "analysis", "fdr", "clin_var", "count_label"))
driver_clinstat_melt = driver_clinstat_melt[driver_clinstat_melt$variable != "n_mut_patients",]

driver_clinstat_melt$association_type = c("1" = "EO_enriched", "-1" = "LO_enriched")[as.character(sign(log(driver_clinstat_melt$effect)))]
driver_levels = unique(driver_clinstat_melt[order(driver_clinstat_melt$fdr),"id"])
driver_clinstat_melt$id = factor(driver_clinstat_melt$id, levels = driver_levels)
driver_clinstat_melt$variable = factor(driver_clinstat_melt$variable, levels = c("n_mut_LO_patients", "n_mut_EO_patients"))

EO_LO_stats_of_drivers = driver_clinstat_melt
save(EO_LO_stats_of_drivers, file = pff("EO_LO_stats_of_drivers.rsav"))


annot_colors = c("n_mut_LO_patients" = "aquamarine4", "n_mut_EO_patients" = "aquamarine2")

plt = ggplot(driver_clinstat_melt, aes(id, value, fill = variable, label = stars.pval(fdr))) + 
		facet_grid(~association_type, scales = "free", space = "free") +
		geom_bar(stat = "identity", position = "fill") + 
		geom_text(size = 5, y = 1) +
		geom_hline(yintercept = global_fraction_EO, linetype = "dotted") + 
		scale_fill_manual("Age group", values = annot_colors) +
		geom_text(aes(id, label = count_label), y = 0.5) +
		scale_y_continuous("fraction of mutated samples") + 
		plot_theme() + 
		theme(legend.position = "bottom") + 
		ggtitle("EO and LO enriched drivers", paste0("FDR < ", FDR_CUTOFF, "; ", NULL_ANALYSIS_TAG))
		

fname = pff("figures/EO_drivers_barchart.pdf")
ggsave(plt, file = fname, width = 8)
file_open_call2(fname)
