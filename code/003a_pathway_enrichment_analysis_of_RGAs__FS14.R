source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(ActivePathways)
library(scales)

load(pff("cgc2024.rsav"))
load(pff("prostate_cancer_genes.rsav"))

FDR_CUTOFF = 0.01


# create a matrix for stringent pathway analysis: FDR-selected genes only, combine SNV, SV, BI, hAMP, Gistic

# SNV, SV drivers first
load(file = pff("merged_pvals_counts.rsav"))
combined_results1 = dcast(merged_pvals_counts, gene ~ mut_type, value.var = "fdr")
# !! keep the gene of the most frequent translocation here for pathway analysis purpose
combined_results1[combined_results1$gene == "ETS", "gene"] = "ERG"


# get patient lists for each SV and SNV RGA
expand_patient_list = function(i, merged_pvals_counts, colname_patients) {
	cat(i, " ")
	gene = merged_pvals_counts[i, "gene"]
	mut_type = merged_pvals_counts[i, "mut_type"]
	patient = strsplit(merged_pvals_counts[i, colname_patients], split = ",")[[1]]
	if(length(patient) == 0) {
		return(NULL)
	}
	data.frame(gene, mut_type, patient, stringsAsFactors = FALSE)	
}

SNV_SV_patients = do.call(rbind, lapply(1:nrow(merged_pvals_counts), expand_patient_list, merged_pvals_counts, "concat_patients"))
SNV_SV_patients[SNV_SV_patients$gene == "ETS", "gene"] = "ERG"



# BI drivers
load(file = pff("BI_results.rsav"))
BI_results = BI_results[BI_results$is_cancer_gene == TRUE,]
BI_results_patients = BI_results[, c("gene", "patients")]
BI_results = unique(BI_results[, c("gene", "pval")])

combined_results2 = merge(combined_results1, BI_results, by = "gene", all = TRUE)
colnames(combined_results2) = c("gene", "CDS", "NC", "SV", "BI")

# get patient lists for each BI
BI_results_patients$mut_type = "BI"
BI_patients = do.call(rbind, lapply(1:nrow(BI_results_patients), expand_patient_list, BI_results_patients, "patients"))

# genes with CNAs 
load(file = pff("gene_CNA_module_dfr.rsav"))
hAmp_results = gene_CNA_module_dfr[gene_CNA_module_dfr$annot == "high_gain", c("gene", "fdr")]
combined_results3 = merge(combined_results2, hAmp_results, by = "gene", all = TRUE)
colnames(combined_results3) = c("gene", "CDS", "NC", "SV", "BI", "hAMP")
cna_results = gene_CNA_module_dfr[gene_CNA_module_dfr$annot %in% c("all_gain", "all_loss"), c("gene", "fdr")]
combined_results5 = merge(combined_results3, cna_results, by = "gene", all = TRUE)
colnames(combined_results5) = c("gene", "CDS", "NC", "SV", "BI", "hAMP", "cna")

# patient lists for AMPs
hAmp_results_patients = gene_CNA_module_dfr[gene_CNA_module_dfr$annot == "high_gain", c("gene", "patients")]
hAmp_results_patients$mut_type = "hAMP"
hAmp_patients = do.call(rbind, lapply(1:nrow(hAmp_results_patients), expand_patient_list, hAmp_results_patients, "patients"))

# patient lists for CNA gains
CNA_gain_results_patients = gene_CNA_module_dfr[gene_CNA_module_dfr$annot == "all_gain", c("gene", "patients")]
CNA_gain_results_patients$mut_type = "CNA__gain"
CNA_gain_patients = do.call(rbind, lapply(1:nrow(CNA_gain_results_patients), expand_patient_list, CNA_gain_results_patients, "patients"))

# patient lists for CNA losses
CNA_loss_results_patients = gene_CNA_module_dfr[gene_CNA_module_dfr$annot == "all_loss", c("gene", "patients")]
CNA_loss_results_patients$mut_type = "CNA__loss"
CNA_loss_patients = do.call(rbind, lapply(1:nrow(CNA_loss_results_patients), expand_patient_list, CNA_loss_results_patients, "patients"))


# merge and clean up P-val martrix
combined_results5[,-1][is.na(combined_results5[,-1])] = 1
combined_results5[,-1][combined_results5[,-1] < 1e-300] = 1e-300
driver_pval_matrix_with_CNAs = combined_results5
save(driver_pval_matrix_with_CNAs, file = pff("driver_pval_matrix_with_CNAs.rsav"))

# all combined
driver_patient_sets_for_pathway_analysis = rbind(SNV_SV_patients, BI_patients, hAmp_patients, CNA_gain_patients, CNA_loss_patients)
save(driver_patient_sets_for_pathway_analysis, file = pff("driver_patient_sets_for_pathway_analysis.rsav"))


# pathway enrichment analysis using multiple types of alterations
GMT_fname = pff("gobp_reac.gmt")
gmt_genes = get_gmt(GMT_fname)
load(file = pff("driver_pval_matrix_with_CNAs.rsav"))
pval_mat = as.matrix(driver_pval_matrix_with_CNAs[,-1])
rownames(pval_mat) = driver_pval_matrix_with_CNAs[,1]

mut_type_colors = c(
	"hAMP" = "darkred",
	"cna" = "darkcyan",
	"BI" = "darkblue",
	"CDS" = "gold",
	"NC" = "darkorange",
	"SV" = "darkgreen"
)

APW_results_integrative = ActivePathways(pval_mat, GMT_fname, 
		correction_method = "fdr",
		cutoff = 0.5,
		geneset_filter = c(25, 250),
		significant = FDR_CUTOFF, 
		custom_colors = mut_type_colors,
		cytoscape_file_tag = pff("figures/ActivePathways__integrative__"))

save(APW_results_integrative, file = pff("APW_results_integrative.rsav"))
fname = pff("figures/APW_results_integrative.csv")
export_as_CSV(APW_results_integrative, file = fname)
file_open_call2(fname)
