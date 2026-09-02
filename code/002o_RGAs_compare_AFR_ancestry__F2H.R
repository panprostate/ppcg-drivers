source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(gtools)

N_AFR_TUMORS = 113
FDR_CUT = 0.05
MIN_COHORT_FREQ = 0.01

load(file = pff("patient_sets_for_drivers.rsav"))
load(file = pff("patient_clinical_data.rsav"))

AFR_patients = unique(patient_clinical_data[
		patient_clinical_data$predicted_ancestry == "AFR" & 
		!is.na(patient_clinical_data$predicted_ancestry),"patient"])

# need to exclude AFR patients from the PPCG sets to do unbiased counting
patient_sets_nonAFR = sapply(patient_sets_for_drivers, function(x) setdiff(x, AFR_patients))
patients_nonAFR = setdiff(patient_clinical_data$patient, AFR_patients)
N_nonAFR_TUMORS = length(patients_nonAFR)

freqs_nonAFR_all = sapply(patient_sets_nonAFR, length)
freqs_nonAFR_all = freqs_nonAFR_all[names(freqs_nonAFR_all) != "AR__hAMP"]

afr_other_muts = read.csv(paste0("DATA_USED__", this_timestamp, 
					"/AFR_RGA_freqs_Vanessa_Veerachai_2026-05-21/RGAs_FOR_VH_Table1_alterations_combined-WJ2.csv"),
		stringsAsFactors = FALSE)
afr_other_muts = afr_other_muts[afr_other_muts$alteration_type %in% c("amplification", "loss", "gain", "SV"),]		
afr_coding_muts = read.csv(paste0("DATA_USED__", this_timestamp, 
					"/AFR_RGA_freqs_Vanessa_Veerachai_2026-05-21/protein-coding_driver_genes_FOR_VH_Table4-WJ.csv"),
		stringsAsFactors = FALSE)

ETS_combined_count = 12 + 2 + 1

afr_other_muts = afr_other_muts[afr_other_muts$RGA != "ETS_combined",]
afr_other_muts = afr_other_muts[afr_other_muts$RGA != "AR",]

# combine MYC,MDM4 gains and amps to avoid thresholding issues
for (rga_to_combine in c("MYC", "MDM4")) {
	# add up counts of gains and amps to gain category in afr
	afr_other_muts[afr_other_muts$RGA == rga_to_combine & afr_other_muts$alteration_type == "gain", "AFR_113"] = 
			sum(as.numeric(afr_other_muts[afr_other_muts$RGA == rga_to_combine & 
					afr_other_muts$alteration_type %in% c("gain", "amplification"), "AFR_113"]))
	afr_other_muts = afr_other_muts[!(afr_other_muts$RGA == rga_to_combine & afr_other_muts$alteration_type == "amplification"),]
	gain_tag = paste0(rga_to_combine, "__gain")
	amp_tag = paste0(rga_to_combine, "__hAMP")
	freqs_nonAFR_all[gain_tag] = freqs_nonAFR_all[gain_tag] + freqs_nonAFR_all[amp_tag]
	freqs_nonAFR_all = freqs_nonAFR_all[names(freqs_nonAFR_all) != amp_tag]
}

# consolidate AFR and nonAFR coding 
freqs_AFR_coding = structure(as.numeric(afr_coding_muts$AFR_113), names = paste0(afr_coding_muts$RGA, "__SNV_CDS"))
freqs_nonAFR_coding = freqs_nonAFR_all[grep("__SNV_CDS$", names(freqs_nonAFR_all))]

# consolidate AFR and nonAFR amps
afr_amplifications = afr_other_muts[afr_other_muts$alteration_type == "amplification",]
freqs_AFR_ampl = structure(as.numeric(afr_amplifications$AFR_113), names = paste0(afr_amplifications$RGA, "__hAMP"))
freqs_nonAFR_ampl = freqs_nonAFR_all[grep("__hAMP$", names(freqs_nonAFR_all))]

# consolidate AFR and nonAFR amps
afr_SVs = afr_other_muts[afr_other_muts$alteration_type == "SV",]
freqs_AFR_SV = structure(as.numeric(afr_SVs$AFR_113), names = paste0(afr_SVs$RGA, "__SV"))
freqs_nonAFR_SV = freqs_nonAFR_all[grep("__SV$", names(freqs_nonAFR_all))]
# manually add ETS freq to AFR 
freqs_AFR_SV[["ETS__SV"]] = ETS_combined_count

afr_gains = afr_other_muts[afr_other_muts$alteration_type == "gain",]
freqs_AFR_gain = structure(as.numeric(afr_gains$AFR_113), names = paste0(afr_gains$RGA, "__gain"))
freqs_nonAFR_gain = freqs_nonAFR_all[grep("__gain$", names(freqs_nonAFR_all))]

afr_losses = afr_other_muts[afr_other_muts$alteration_type == "loss",]
freqs_AFR_loss = structure(as.numeric(afr_losses$AFR_113), names = paste0(afr_losses$RGA, "__loss"))
freqs_nonAFR_loss = freqs_nonAFR_all[grep("__loss$", names(freqs_nonAFR_all))]


test_RGA_afr_vs_rest = function(rga, freqs_AFR, freqs_nonAFR, N_AFR_TUMORS, N_nonAFR_TUMORS, MIN_COHORT_FREQ) {
	
	n_AFR = freqs_AFR[[rga]]
	n_nonAFR = freqs_nonAFR[[rga]]

	freq_AFR = n_AFR / N_AFR_TUMORS
	freq_nonAFR = n_nonAFR / N_nonAFR_TUMORS
	
	if (freq_AFR < MIN_COHORT_FREQ & freq_nonAFR < MIN_COHORT_FREQ) {
		cat(paste("skip lowfreq rga", rga, "\n"))
		return(NULL)
	}
	
	bintest = binom.test(n_AFR, N_AFR_TUMORS, n_nonAFR / N_nonAFR_TUMORS)
	pval = bintest$p.value
	
	pct_afr = n_AFR / N_AFR_TUMORS
	pct_nonafr = n_nonAFR / N_nonAFR_TUMORS
	n_afr = n_AFR
	n_nonafr = n_nonAFR
	
	obs_exp_FC = pct_afr / pct_nonafr
	
	mut_type = gsub("(.+)__(.+)", "\\2", rga)
	driver_id_label = gsub("(.+)__(.+)", "\\1", rga) 
	
	dfr = data.frame(rga, pval, pct_afr, pct_nonafr, n_afr, n_nonafr, obs_exp_FC, mut_type, driver_id_label,
			stringsAsFactors = FALSE)
	rownames(dfr) = NULL
	dfr
}



freqs_AFR_combined = c(freqs_AFR_coding, freqs_AFR_ampl, freqs_AFR_SV, freqs_AFR_gain, freqs_AFR_loss)
freqs_nonAFR_combined = c(freqs_nonAFR_coding, freqs_nonAFR_ampl, freqs_nonAFR_SV, freqs_nonAFR_gain, freqs_nonAFR_loss)
AFR_rga_stats = do.call(rbind, lapply(names(freqs_nonAFR_combined), 
		test_RGA_afr_vs_rest, freqs_AFR_combined, freqs_nonAFR_combined, N_AFR_TUMORS, N_nonAFR_TUMORS, MIN_COHORT_FREQ))

AFR_rga_stats$fdr = p.adjust(AFR_rga_stats$pval, method = "fdr")
AFR_rga_stats = AFR_rga_stats[order(AFR_rga_stats$pval),]
save(AFR_rga_stats, file = pff("AFR_rga_stats.rsav"))


LOG2FC_CAP = 2 

AFR_rga_stats_signf = AFR_rga_stats[AFR_rga_stats$fdr < FDR_CUT,]
AFR_rga_stats_signf$log2fc_cap = log2(AFR_rga_stats_signf$obs_exp_FC)
AFR_rga_stats_signf$log2fc_cap[AFR_rga_stats_signf$log2fc_cap > LOG2FC_CAP] = LOG2FC_CAP
AFR_rga_stats_signf$log2fc_cap[AFR_rga_stats_signf$log2fc_cap < -LOG2FC_CAP] = -LOG2FC_CAP

AFR_rga_stats_signf$fdr_score = -log10(AFR_rga_stats_signf$fdr) * sign(log2(AFR_rga_stats_signf$obs_exp_FC))

AFR_rga_stats_signf$rga = factor(AFR_rga_stats_signf$rga, levels = unique(AFR_rga_stats_signf$rga[order(AFR_rga_stats_signf$fdr_score)]))
AFR_rga_stats_signf$freq_label = paste0(signif(100 * AFR_rga_stats_signf$pct_afr, 2), "% vs ", signif(100 * AFR_rga_stats_signf$pct_nonafr, 2), "%")

plot_title = paste0("RGAs: AFR vs nonAFR, FDR<", FDR_CUT)
plt3 = ggplot(AFR_rga_stats_signf, aes(factor("tumours_afrs"), rga, fill = log2fc_cap, label = stars.pval(fdr))) + 
		geom_tile() +
		geom_text( angle = 0) +
		scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", limits = c(-2, 2)) +
		ggtitle(NULL, plot_title) + 
		geom_text(aes(label = freq_label), color = "darkcyan") + 
		plot_theme()
		
fname = pff("figures/AFR_RGAs_vs_others_tiles.pdf")		
ggsave(plt3, file = fname, height = 7, width = 5)
file_open_call2(fname)
