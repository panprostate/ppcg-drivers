source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library("ggrepel")
library(gtools)


FDR_CUTOFF = 0.05
MIN_MUT_FREQ = 0.01

test_RGA_prim_vs_met = function(
		driver_id, patient_sets_for_drivers, patient_sets_for_drivers_mets, prim_patients, mets_patients, prim_met_patients) {
	
	cat(driver_id, " ")
	
	prim_patients_with_driver = patient_sets_for_drivers[[driver_id]]
	mets_patients_with_driver = patient_sets_for_drivers_mets[[driver_id]]

	prob_background = length(prim_patients_with_driver)/length(prim_patients)
	binom_test = binom.test(length(mets_patients_with_driver), length(mets_patients), prob_background)
	pval = binom_test$p.value

	pct_prim = length(prim_patients_with_driver) / length(prim_patients)
	pct_mets = length(mets_patients_with_driver) / length(mets_patients)
	n_prim = length(prim_patients_with_driver)
	n_mets = length(mets_patients_with_driver)
	
	exp_mets_with_driver = prob_background * length(mets_patients)
	obs_mets_with_driver = length(mets_patients_with_driver)
	obs_exp_FC = pct_mets / pct_prim
	
	mut_type = gsub("(.+)__(.+)", "\\2", driver_id)
	driver_id_label = gsub("(.+)__(.+)", "\\1", driver_id) 
	
	dfr = data.frame(driver_id, pval, pct_prim, pct_mets, n_prim, n_mets, obs_exp_FC, mut_type, driver_id_label,
			stringsAsFactors = FALSE)
	rownames(dfr) = NULL
	dfr
}

load(file = pff("all_patients.rsav"))
load(file = pff("allmets_tumor_ids.rsav"))
load(file = pff("patient_sets_for_drivers.rsav"))
load(file = pff("patient_sets_for_drivers_mets.rsav"))

# exclude lower-frequency events
patient_sets_for_drivers = patient_sets_for_drivers[sapply(patient_sets_for_drivers, length) > MIN_MUT_FREQ * length(all_patients)]
alteration_types = gsub("(.+)__(.+)", "\\2", names(patient_sets_for_drivers))

# stat analysis - drivers in prim vs mets
mets_patients = paste0(gsub("._DNA$", "", allmets_tumor_ids), "_M")
prim_patients = paste0(all_patients, "_P")
prim_met_patients = unique(c(prim_patients, mets_patients))
RGAs_to_test = names(patient_sets_for_drivers)

# for each gene run stats
driver_freq_prim_vs_met = do.call(rbind, lapply(RGAs_to_test, test_RGA_prim_vs_met, 
		patient_sets_for_drivers, patient_sets_for_drivers_mets, prim_patients, mets_patients, prim_met_patients))
driver_freq_prim_vs_met$tumor_group = "all_mets"


# repeat the same analysis, but keep only de novo mets
print(load(pff("denovomets_tumor_ids.rsav")))
denovomets_tumor_ids = gsub("._DNA$", "", denovomets_tumor_ids)
patient_sets_for_drivers_mets_denovo = lapply(patient_sets_for_drivers_mets, function(x) intersect(x, denovomets_tumor_ids))
mets_patients_denovo = paste0(denovomets_tumor_ids, "_M")
prim_met_patients_denovo = c(prim_patients, mets_patients_denovo)

# for each gene run stats
driver_freq_prim_vs_met_denovo = do.call(rbind, lapply(RGAs_to_test, test_RGA_prim_vs_met, 
		patient_sets_for_drivers, patient_sets_for_drivers_mets_denovo, prim_patients, mets_patients_denovo, prim_met_patients_denovo))
driver_freq_prim_vs_met_denovo$tumor_group = "denovo_mets"

driver_freq_prim_vs_met_all = rbind(driver_freq_prim_vs_met, driver_freq_prim_vs_met_denovo)
driver_freq_prim_vs_met_all$fdr = p.adjust(driver_freq_prim_vs_met_all$pval, method = "fdr")
save(driver_freq_prim_vs_met_all, file = pff("driver_freq_prim_vs_met_all.rsav"))




load(file = pff("driver_freq_prim_vs_met_all.rsav"))
mets_drivers_signf = driver_freq_prim_vs_met_all[driver_freq_prim_vs_met_all$fdr < 0.05,]

LOG2FC_CAP = 2 
mets_drivers_signf$log2fc_cap = log2(mets_drivers_signf$obs_exp_FC)
mets_drivers_signf$log2fc_cap[mets_drivers_signf$log2fc_cap > LOG2FC_CAP] = LOG2FC_CAP
mets_drivers_signf$log2fc_cap[mets_drivers_signf$log2fc_cap < -LOG2FC_CAP] = -LOG2FC_CAP
mets_drivers_signf$tumor_group = factor(mets_drivers_signf$tumor_group, levels = c("all_mets", "denovo_mets"))
mets_drivers_signf$driver_id = factor(mets_drivers_signf$driver_id, 
		levels = names(sort(by(-log10(mets_drivers_signf$fdr), mets_drivers_signf$driver_id, function(x) sum(x, na.rm = T)))))

mets_drivers_signf$freq_label = paste0(signif(100 * mets_drivers_signf$pct_mets, 2), "% vs ", signif(100 * mets_drivers_signf$pct_prim, 2), "%")

plot_title = paste0("RGAs in mets vs primaries, FDR<", FDR_CUTOFF)
plt3 = ggplot(mets_drivers_signf, aes(tumor_group, driver_id, fill = log2fc_cap, label = stars.pval(fdr))) + 
		geom_tile() +
		geom_text( angle = 0) +
		scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", limits = c(-2, 2)) +
		ggtitle(NULL, plot_title) + 
		geom_text(aes(label = freq_label), color = "darkcyan") + 
		plot_theme()
		

fname = pff("figures/Mets_RGAs_vs_primary_tiles.pdf")		
ggsave(plt3, file = fname, height = 7, width = 5)
file_open_call2(fname)
