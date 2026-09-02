source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library("patchwork")

FDR_CUTOFF = 0.05
MIN_N_SIGS = 1000
MIN_N_PATIENTS = 10

add_mut_tag = function(dfr) {
	dfr$mut_tag = paste(dfr$patient, dfr$chr, dfr$pos1, dfr$pos2, dfr$ref, dfr$alt, sep = ":")
	rownames(dfr) = NULL
	dfr
}

get_driver_muts = function(driver_ids, results_signf_merged_annot, variants_to_elements, mut_tag_colname, do_add_mut_tag) {
	this_elements = results_signf_merged_annot[results_signf_merged_annot$annots_MAIN %in% driver_ids, "element_ids"]
	this_elements = unique(unlist(strsplit(this_elements, split = ",")))
	this_muts = do.call(rbind, variants_to_elements[this_elements])
	if (do_add_mut_tag == TRUE) {
		this_muts = add_mut_tag(this_muts)
	}
	this_muts = this_muts[!duplicated(this_muts[, mut_tag_colname]),]
	this_muts
}

test_signt = function(sig, id, this_muts, prepared_variants, mut_tag_colname, mut_sig_colname) {
	is_mutated = factor(prepared_variants[, mut_tag_colname] %in% this_muts[, mut_tag_colname], levels = c("TRUE", "FALSE"))
	is_mutsig = factor(!is.na(prepared_variants[, mut_sig_colname]) & prepared_variants[, mut_sig_colname] == sig, levels = c("TRUE", "FALSE"))
	ft = fisher.test(is_mutated, is_mutsig)
	log2or = log2(ft$estimate)
	pval = ft$p.value
	muts_obs = sum(this_muts[, mut_sig_colname] == sig)
	p_mut = sum(is_mutated == "TRUE")/length(is_mutated)
	p_sig = sum(is_mutsig == "TRUE")/length(is_mutated)
	muts_exp = length(is_mutated) * p_mut * p_sig
	data.frame(id, sig, pval, log2or, muts_obs, muts_exp, stringsAsFactors = FALSE)
}

test_driver_sigs = function(id, results_signf_merged_annot, variants_to_elements, prepared_variants, all_sigs, 
		mut_tag_colname, mut_sig_colname, do_add_mut_tag) {
	cat(id, " ")
	this_muts = get_driver_muts(id, results_signf_merged_annot, variants_to_elements, mut_tag_colname, do_add_mut_tag)
	stats = do.call(rbind, lapply(all_sigs, test_signt, id, this_muts, prepared_variants, mut_tag_colname, mut_sig_colname))
	rownames(stats) = NULL
	stats
}

test_sigs_of_muts = function(id, this_muts, prepared_variants, all_sigs, mut_tag_colname, mut_sig_colname) {
	
	cat(id, " ")
	stats = do.call(rbind, lapply(all_sigs, test_signt, id, this_muts, prepared_variants, mut_tag_colname, mut_sig_colname))
	rownames(stats) = NULL
	stats
}

count_patients_from_results = function(res) {
		c(by(res$patient_ids, res$annots_MAIN, function(x) length(unique(unlist(strsplit(x, s = ","))))))
}


#
# first, SNV and indel signatures
#
load(file = pff("variants__SNV_indel.rsav"))
load(file = pff("prepared_variants__SNV_indel.rsav"))
load(file = pff("variants_to_elements.rsav"))
variants_to_elements_SNV = variants_to_elements[['SNV_indel']]
prepared_variants_SNV = add_mut_tag(prepared_variants)
load(file = pff("results_signf_merged_annot.rsav"))
# drivers correspond to several overlapping elements
results_signf_merged_annot_SNV = results_signf_merged_annot[results_signf_merged_annot$mut_type %in% c("SNV_CDS", "SNV_NC"),]
rm(results_signf_merged_annot, variants_to_elements, prepared_variants)
gc()


SBS_sequencing_artefacts = readLines(paste0("DATA_USED__", this_timestamp, "/SBS_sequencing_artefacts.txt"))
SBS_sequencing_artefacts = c(SBS_sequencing_artefacts, "SBS96D")
SNV_sigs_table = table(prepared_variants_SNV$top_mutsig)
SNV_sigs_table = SNV_sigs_table[SNV_sigs_table > MIN_N_SIGS]
SNV_sigs = names(SNV_sigs_table)
SNV_sigs = setdiff(SNV_sigs, SBS_sequencing_artefacts)

all_SNV_drivers = unique(results_signf_merged_annot_SNV$annots_MAIN)

n_samples_per_SNV_driver = count_patients_from_results(results_signf_merged_annot_SNV)
SNV_drivers_to_test = names(which(n_samples_per_SNV_driver >= MIN_N_PATIENTS))

# for instances testing all driver mutations combined, or subsets of these
SNV_driver_muts = get_driver_muts(all_SNV_drivers, results_signf_merged_annot_SNV, variants_to_elements_SNV, 
				mut_tag_colname = "mut_tag", do_add_mut_tag = TRUE)

# test drivers 1 by 1
SNV_drivers_res = do.call(rbind, mclapply(SNV_drivers_to_test, test_driver_sigs, 
		results_signf_merged_annot_SNV, variants_to_elements_SNV, prepared_variants_SNV, SNV_sigs, "mut_tag", "top_mutsig", 
		do_add_mut_tag = TRUE,
		mc.cores = 4 ))

# protein-truncating mutations, stopgains and frameshifts
trunc_annots = c("stopgain", "frameshift substitution")
trunc_muts = SNV_driver_muts[SNV_driver_muts$ExonicFunc.refGene %in% trunc_annots,]
trunc_res = test_sigs_of_muts("trunc", trunc_muts, prepared_variants_SNV, SNV_sigs, "mut_tag", "top_mutsig")

# substitution muts
subst_annots = c("nonsynonymous SNV", "nonframeshift substitution")
subst_muts = SNV_driver_muts[SNV_driver_muts$ExonicFunc.refGene %in% subst_annots,]
subst_res = test_sigs_of_muts("subst", subst_muts, prepared_variants_SNV, SNV_sigs, "mut_tag", "top_mutsig")

# non-coding mutations
nc_annots = ""
nc_muts = SNV_driver_muts[SNV_driver_muts$ExonicFunc.refGene %in% nc_annots,]
nc_res = test_sigs_of_muts("nc", nc_muts, prepared_variants_SNV, SNV_sigs, "mut_tag", "top_mutsig")
dim(nc_muts)
#[1] 732  29	# 2026-06-26

combined_SNV_res = rbind(SNV_drivers_res, trunc_res, subst_res, nc_res)
combined_SNV_res$mut_type = "SNV_indel"
combined_SNV_res$mut_type2 = "SNV"
combined_SNV_res$mut_type2[grep("^ID", combined_SNV_res$sig)] = "indel"


# test all indels among drivers
assoc_indels_among_drivers = table(prepared_variants_SNV$mut_tag %in% SNV_driver_muts$mut_tag, grepl("^ID", prepared_variants_SNV$top_mutsig))
ft_indels_among_drivers = fisher.test(assoc_indels_among_drivers, alt = "g")

# SV drivers
load(file = pff("variants_SV.rsav"))
load(file = pff("prepared_variants_SV.rsav"))
load(file = pff("variants_to_elements.rsav"))
load(file = pff("results_signf_merged_annot.rsav"))
variants_to_elements_SV = variants_to_elements[['SV']]
results_signf_merged_annot_SV = results_signf_merged_annot[results_signf_merged_annot$mut_type == "SV",]

rm(results_signf_merged_annot, variants_to_elements)
gc()

SV_types = unique(prepared_variants_SV$svclass)
all_SV_drivers = unique(results_signf_merged_annot_SV$annots_MAIN)
n_samples_per_SV_driver = count_patients_from_results(results_signf_merged_annot_SV)
SV_drivers_to_test = names(which(n_samples_per_SV_driver >= MIN_N_PATIENTS))

# testing all driver mutations combined and respective subsets
all_driver_SV_muts = get_driver_muts(all_SV_drivers, results_signf_merged_annot_SV, variants_to_elements_SV, 
				mut_tag_colname = "SV_tag", do_add_mut_tag = FALSE)

# test drivers 1 by 1
SV_drivers_res = do.call(rbind, mclapply(SV_drivers_to_test, test_driver_sigs, 
		results_signf_merged_annot_SV, variants_to_elements_SV, prepared_variants_SV, SV_types, "SV_tag", "top_mut_signt", FALSE,
		mc.cores = 4))
SV_drivers_res$mut_type = "SV"
SV_drivers_res$mut_type2 = "SV"

# test all SVs among drivers
assoc_tra_among_drivers = table(prepared_variants_SV$SV_tag %in% all_driver_SV_muts$SV_tag, prepared_variants_SV$top_mut_signt == "TRA")
ft_tra_among_drivers = fisher.test(assoc_tra_among_drivers, alt = "g")

#
# combine SV and SNV results into one figure
#
combined_res = rbind(SV_drivers_res, combined_SNV_res)
combined_res = combined_res[order(combined_res$pval),]
combined_res$fdr = p.adjust(combined_res$pval, method = "fdr")

combined_res$fdr_capped = combined_res$fdr
combined_res[combined_res$fdr_capped < 1e-16, "fdr_capped"] = 1e-16
combined_res[combined_res$fdr_capped > FDR_CUTOFF, "fdr_capped"] = NA

combined_res$log2or_capped = combined_res$log2or
combined_res[combined_res$log2or_capped < -5, "log2or_capped"] = -5
combined_res[combined_res$log2or_capped > 5, "log2or_capped"] = 5

ids_to_exclude = names(which(by(combined_res$fdr_capped, combined_res$id, function(x) all(is.na(x)))))
sigs_to_exclude = names(which(by(combined_res$fdr_capped, combined_res$sig, function(x) all(is.na(x)))))

combined_res2 = combined_res[!combined_res$id %in% ids_to_exclude,]
combined_res2 = combined_res2[!combined_res2$sig %in% sigs_to_exclude,]
combined_res2$mut_type2 = factor(combined_res2$mut_type2, levels = c("SNV", "indel", "SV"))

combined_res2$id = factor(combined_res2$id, 
	levels = names(sort(by(combined_res2$fdr, combined_res2$id, function(x) -log10(prod(x, na.rm = T))))))
combined_res2$sig = factor(combined_res2$sig, 
	levels = rev(names(sort(by(combined_res2$fdr, combined_res2$sig, function(x) -log10(prod(x, na.rm = T)))))))
	
meta_ids = c("nc", "subst", "trunc")
new_id_levels = c(setdiff(levels(combined_res2$id), meta_ids), meta_ids)
combined_res2$id = factor(combined_res2$id, levels = new_id_levels)
combined_res2$mut_type2 = factor(combined_res2$mut_type2, levels = c("indel", "SV", "SNV"))

plt = ggplot(combined_res2, aes(sig, id, size = -log10(fdr_capped), fill = log2or_capped)) + 
		geom_point(color = "black", shape = 22) + 
		facet_grid(~mut_type2, scale = "free_x", space = "free") + 
		plot_theme() +
		scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
		plot_theme() + 
		ggtitle("mutsigs in drivers", paste0("FDR<", FDR_CUTOFF)) + 
		theme(legend.position = "bottom")

fname = pff("figures/MutSigs_in_SV_drivers_dotplot.pdf")
ggsave(plt, file = fname, width = 6)
file_open_call2(fname)

res_mutsigs_of_drivers = combined_res
save(res_mutsigs_of_drivers, file = pff("res_mutsigs_of_drivers.rsav"))
