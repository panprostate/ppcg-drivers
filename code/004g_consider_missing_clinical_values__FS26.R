source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library("patchwork")
library("gtools")

load(pff("patient_clinical_data.rsav"))
load(pff("all_patients.rsav"))
FDR_CUTOFF = 0.05
MIN_DRIVER_FREQ = 0.01
MIN_MUT_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)

clinical_values = c(MFS_COVARIATES, "mets_time_10year", "mets_event_10year")

# select the samples with at least some missing data
clin_data_missing = patient_clinical_data[!complete.cases(patient_clinical_data[,clinical_values]), c("patient", clinical_values)]
# remove actual values, keep only whether it is an NA or not
clin_data_missing[,setdiff(colnames(clin_data_missing), "patient")] = !is.na(clin_data_missing[,setdiff(colnames(clin_data_missing), "patient")])

# use the oncoprint ordering of rows and columns, convert to numeric matrix for sorting
rownames(clin_data_missing) =  clin_data_missing$patient
clin_data_missing1 = memoSort(0 + t(clin_data_missing[,-1]))
clinvar_order = rownames(clin_data_missing1)
patient_order = colnames(clin_data_missing1)
clin_data_missing = melt(clin_data_missing, id.var = "patient")
clin_data_missing$patient = factor(clin_data_missing$patient, levels = patient_order)
clin_data_missing$variable = factor(clin_data_missing$variable, levels = clinvar_order)


plot_title = paste0("Missing clin.var for ", length(unique(clin_data_missing$patient)), " donors")
plt_matrix = ggplot(clin_data_missing, aes(patient, variable, fill = value)) + 
		geom_tile() + 
		plot_theme() +
		scale_x_discrete("Donors with missing values", labels = NULL) + 
		ggtitle(NULL, plot_title) +
		scale_fill_manual("Variable", values = c("TRUE" = "white", "FALSE" = "darkcyan"), labels = c("TRUE" = "present", "FALSE" = "missing"))
		

# test major drivers for enrichment
patients_w_missing = as.character(unique(clin_data_missing$patient))
all_patients = as.character(unique(patient_clinical_data$patient))

load(pff("patient_sets_for_drivers.rsav"))
all_drivers = names(which(sapply(patient_sets_for_drivers, length) > MIN_MUT_PATIENTS))

# ask if an RGA is enriched/depleted in samples with some missing clinical data
fisher_test = function(dr, all_patients, patient_sets_for_drivers, patients_w_missing) {
	ft = fisher.test(all_patients %in% patient_sets_for_drivers[[dr]], all_patients %in% patients_w_missing)
	pval = ft$p.value
	or = ft$estimate
	n_missing_patients = length(intersect(patients_w_missing, patient_sets_for_drivers[[dr]]))	
	data.frame(RGA = dr, pval, or, n_missing_patients, stringsAsFactors = FALSE)
}


gene_enrichments = do.call(rbind, lapply(all_drivers, fisher_test, all_patients, patient_sets_for_drivers, patients_w_missing))
rownames(gene_enrichments) = NULL
gene_enrichments$fdr = p.adjust(gene_enrichments$pval, method = "fdr")
gene_enrichments = gene_enrichments[order(gene_enrichments$fdr),]

gene_enrichments_signf = gene_enrichments[gene_enrichments$fdr < FDR_CUTOFF,]

# cap log2 OR to +/-4
gene_enrichments_signf$or_log2cap4 = log2(gene_enrichments_signf$or)
gene_enrichments_signf$or_log2cap4[gene_enrichments_signf$or_log2cap4 > 4] = 4
gene_enrichments_signf$or_log2cap4[gene_enrichments_signf$or_log2cap4 < -4] = -4

gene_enrichments_signf = gene_enrichments_signf[order(gene_enrichments_signf$fdr),]
gene_enrichments_signf$RGA = factor(gene_enrichments_signf$RGA, levels = gene_enrichments_signf$RGA)
gene_enrichments_signf$result_label = paste0(stars.pval(gene_enrichments_signf$fdr), "\n\n", gene_enrichments_signf$n_missing_patients)
plot_title = paste0("RGAs enriched in patients with missing values (", nrow(gene_enrichments_signf), ")")

plt_enrichments = ggplot(gene_enrichments_signf, aes(RGA, factor(2), size = -log10(fdr), fill = or_log2cap4, label = result_label)) + 
		geom_point(pch = 21) +
		scale_fill_gradient2(low = "blue", high = "red", mid = "white") + 
		geom_text(size = 3) + 
		plot_theme() + 
		ggtitle(NULL, plot_title)
		
plt_combined = (plt_matrix / plt_enrichments) +   plot_layout(heights = c(2, 0.8))


fname = pff("figures/missing_clinical_data_matrix.pdf")
ggsave(plt_combined, file = fname, width = 7, height = 7)
file_open_call2(fname)


missing_data_gene_enrichments_signf = gene_enrichments_signf
save(missing_data_gene_enrichments_signf, file = pff("missing_data_gene_enrichments_signf.rsav"))
